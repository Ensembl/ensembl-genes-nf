nextflow.enable.dsl=2

include { ENA_COMPUTE_MD5 }   from '../modules/compute_md5.nf'
include { ENA_FTP_UPLOAD }    from '../modules/upload_ftp.nf'
include { ENA_GENERATE_XML }  from '../modules/generate_xml.nf'
include { ENA_SUBMIT_WEBIN as ENA_SUBMIT_PROJECT } from '../modules/submit_webin.nf'
include { ENA_SUBMIT_WEBIN as ENA_SUBMIT_ANALYSIS } from '../modules/submit_webin.nf'
include { ENA_POLL_WEBIN as ENA_POLL_PROJECT } from '../modules/poll_webin.nf'
include { ENA_POLL_WEBIN as ENA_POLL_ANALYSIS } from '../modules/poll_webin.nf'
include { ENA_GENERATE_PROJECT_XML } from '../modules/generate_project_xml.nf'
include { ENA_EXPAND_FILE_MANIFEST } from '../modules/expand_file_manifest.nf'
include { ENA_CONVERT_TO_CRAM } from '../modules/convert_to_cram.nf'

// Parse an annotation-level manifest and dispatch all files for each analysis.
workflow ENA_SUBMIT_WORKFLOW {
    assert params.manifest,       "--manifest is required"
    assert params.webin_user,     "--webin_user is required"
    assert params.webin_password, "--webin_password is required"
    assert params.mode in ['test', 'prod'], "--mode must be 'test' or 'prod'"
    assert params.outdir,        "--outdir is required"

    def ch_webin_user     = Channel.value(params.webin_user)
    def ch_webin_password = Channel.value(params.webin_password)
    // Base endpoint derived from mode unless overridden
    def ch_webin_base = Channel.value(params.webin_base ?: ((params.mode == 'prod') ? 'https://www.ebi.ac.uk/ena/submit/webin-v2' : 'https://wwwdev.ebi.ac.uk/ena/submit/webin-v2'))

    Channel
      .fromPath(params.manifest)
      .splitCsv(header:true, sep:'\t')
      .map { row ->
          if (!row.files_tsv) { throw new RuntimeException("Manifest row missing files_tsv") }
          if (!row.assembly_accession) { throw new RuntimeException("Manifest row missing assembly_accession") }
          if (!row.last_geneset_update && !row.partial_release_label) { throw new RuntimeException("Manifest row missing last_geneset_update/partial_release_label") }
          def release = row.partial_release_label ?: "${row.assembly_accession}-Ensembl-${row.last_geneset_update}"
          def prj_alias = row.project_alias ?: ("prj_${release}").replaceAll('[^A-Za-z0-9._-]', '_')
          // Use derived child project alias as a refname (not accession) unless an existing study is provided.
          if (!row.study) { row.study = prj_alias }
          if (params.umbrella_study && !row.umbrella_study) { row.umbrella_study = params.umbrella_study }
          if (row.umbrella_study && !(row.analysis_attributes ?: '').contains('attr_umbrella_study=')) {
              row.analysis_attributes = [row.analysis_attributes, "attr_umbrella_study=${row.umbrella_study}"].findAll { it }.join('; ')
          }
          def id = row.analysis_alias ?: "rnaseq_alignment_evidence_${row.assembly_accession}_${release}".replaceAll('[^A-Za-z0-9._-]', '_')
          def meta = [ id:id, project_alias: prj_alias, assembly: row.assembly_accession, release: release ]
          tuple(meta, row)
      }
      .set { analyses }

    // Single approach: derive unique Projects from manifest (assembly + release), always
    def ch_proj_rows = analyses
        .map { meta, row ->
            def alias = meta.project_alias
            def title = "Annotation evidence project for ${meta.assembly}, ${meta.release}"
            def description = params.project_description ?: "Annotation evidence project for ${meta.assembly}, ${meta.release}"
            // Key by alias, carry a single meta map
            tuple(alias, [ alias: alias, name: alias, title: title, description: description, hold_until: params.hold_until ?: '' ])
        }
        .groupTuple()
        .map { alias, metas -> metas[0] } // one meta per alias

    ENA_GENERATE_PROJECT_XML(ch_proj_rows)
    ENA_SUBMIT_PROJECT(ENA_GENERATE_PROJECT_XML.out.xml, ch_webin_base, ch_webin_user, ch_webin_password)
    def ch_proj_files = ENA_SUBMIT_PROJECT.out.queued.map { meta, f -> f }

    ENA_POLL_PROJECT(
        ch_proj_files.collect(),
        ch_webin_user,
        ch_webin_password,
        Channel.value(params.poll_interval    ?: 20),
        Channel.value(params.poll_max_attempts ?: 30)
        )

    ENA_EXPAND_FILE_MANIFEST(analyses)

    def file_inputs = ENA_EXPAND_FILE_MANIFEST.out.expanded
        .map { meta, row, expanded_tsv -> expanded_tsv }
        .splitCsv(header:true, sep:'\t')
        .map { file_row ->
            def meta = [
                id: file_row.analysis_id,
                project_alias: file_row.project_alias,
                assembly: file_row.assembly,
                release: file_row.release,
            ]
            def row = [
                study: file_row.study,
                umbrella_study: file_row.umbrella_study,
                analysis_alias: file_row.analysis_alias,
                title: file_row.title,
                description: file_row.description,
                assembly_accession: file_row.assembly_accession,
                last_geneset_update: file_row.last_geneset_update,
                partial_release_label: file_row.partial_release_label,
                species: file_row.species,
                taxon_id: file_row.taxon_id,
                ref_seqs: file_row.ref_seqs,
                analysis_links: file_row.analysis_links,
                analysis_attributes: file_row.analysis_attributes,
                analysis_type: file_row.analysis_type,
                omit_run_refs_in_test: file_row.omit_run_refs_in_test,
            ]
            def file_meta = [
                analysis_id: file_row.analysis_id,
                file_type: file_row.file_type,
                remote_name: file_row.remote_name,
                run_accession: file_row.run_accession,
                sample_accession: file_row.sample_accession,
                experiment_accession: file_row.experiment_accession,
                is_index: false,
            ]
            tuple(meta, row, file_meta, file(file_row.file_path))
        }

    // Optionally convert BAMs to CRAM and retain the CRAI as an upload-only
    // companion file. The CRAI is excluded from the ANALYSIS XML itself.
    def md5_inputs
    def convert_to_cram = params.convert_to_cram?.toString()?.toLowerCase() in ['1', 'true', 'yes']
    if (convert_to_cram) {
        assert params.reference_fasta, "--reference_fasta is required with --convert_to_cram true"
        def reference = file(params.reference_fasta)
        def reference_fai = file("${params.reference_fasta}.fai")
        assert reference.exists(), "Reference FASTA not found: ${params.reference_fasta}"
        assert reference_fai.exists(), "Reference FASTA index not found: ${params.reference_fasta}.fai"

        def bam_inputs = file_inputs.filter { meta, row, file_meta, file ->
            (file_meta.file_type ?: '').toString().toLowerCase() == 'bam'
        }
        def non_bam_inputs = file_inputs.filter { meta, row, file_meta, file ->
            (file_meta.file_type ?: '').toString().toLowerCase() != 'bam'
        }
        ENA_CONVERT_TO_CRAM(bam_inputs, Channel.value(reference), Channel.value(reference_fai))
        def converted_inputs = ENA_CONVERT_TO_CRAM.out.converted.flatMap { meta, row, file_meta, cram, crai ->
            def cram_name = (file_meta.remote_name ?: cram.getName()).replaceFirst(/\.bam$/, '.cram')
            def cram_meta = file_meta + [file_type: 'cram', remote_name: cram_name, is_index: false]
            def crai_meta = file_meta + [file_type: 'crai', remote_name: "${cram_name}.crai", is_index: true]
            [
                tuple(meta, row, cram_meta, cram),
                tuple(meta, row, crai_meta, crai),
            ]
        }
        md5_inputs = converted_inputs.mix(non_bam_inputs)
    } else {
        md5_inputs = file_inputs
    }

    // Compute md5 per file
    md5s = ENA_COMPUTE_MD5(md5_inputs)

    ENA_FTP_UPLOAD(
        md5s,
        params.remote_dir ?: '',
        params.webin_ftp_host ?: 'webin2.ebi.ac.uk',
        ch_webin_user,
        ch_webin_password
        )

    def uploaded_after_projects = ENA_FTP_UPLOAD.out.uploaded
        .combine(ENA_POLL_PROJECT.out.accessions)
        .map { meta, row, file_meta, f, md5, accessions -> tuple(meta.id, meta, row, file_meta, f, md5) }
        .groupTuple(by: 0)
        .map { id, metas, rows, file_metas, files, md5s -> tuple(metas[0], rows[0], file_metas, files, md5s) }

    // Generate one analysis XML after the derived project has been accepted by Webin.
    ENA_GENERATE_XML(
        uploaded_after_projects,
        params.remote_dir ?: '',
        params.hold_until ?: ''
        )

    // Submit to async queue — one POST per annotation analysis, returns immediately with a submission ID.
    ENA_SUBMIT_ANALYSIS(
        ENA_GENERATE_XML.out.xml,
        ch_webin_base,
        ch_webin_user,
        ch_webin_password
        )
    def ch_analysis_files = ENA_SUBMIT_ANALYSIS.out.queued.map { meta, f -> f }

    ENA_POLL_ANALYSIS(
        ch_analysis_files.collect(),
        ch_webin_user,
        ch_webin_password,
        Channel.value(params.poll_interval    ?: 20),
        Channel.value(params.poll_max_attempts ?: 30)
        )

    ENA_POLL_ANALYSIS.out.accessions.view { f -> "Accessions written to: ${f}" }
}
