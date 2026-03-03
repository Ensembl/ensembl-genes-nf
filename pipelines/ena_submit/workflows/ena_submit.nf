nextflow.enable.dsl=2

include { ENA_COMPUTE_MD5 }   from '../modules/compute_md5.nf'
include { ENA_FTP_UPLOAD }    from '../modules/upload_ftp.nf'
include { ENA_GENERATE_XML }  from '../modules/generate_xml.nf'
include { ENA_SUBMIT_WEBIN as ENA_SUBMIT_PROJECT } from '../modules/submit_webin.nf'
include { ENA_SUBMIT_WEBIN as ENA_SUBMIT_ANALYSIS } from '../modules/submit_webin.nf'
include { ENA_POLL_WEBIN } from '../modules/poll_webin.nf'
include { ENA_GENERATE_PROJECT_XML } from '../modules/generate_project_xml.nf'

// Parse manifest and dispatch one analysis per row/file
workflow ENA_SUBMIT_WORKFLOW {
    assert params.manifest,       "--manifest is required"
    assert params.webin_user,     "--webin_user is required"
    assert params.webin_password, "--webin_password is required"
    assert params.mode in ['test', 'prod'], "--mode must be 'test' or 'prod'"
    assert params.outdir,        "--outdir is required"
    assert params.release,       "--release is required (e.g. Ensembl_110)"

    def ch_webin_user     = Channel.value(params.webin_user)
    def ch_webin_password = Channel.value(params.webin_password)
    // Base endpoint derived from mode unless overridden
    def ch_webin_base = Channel.value(params.webin_base ?: ((params.mode == 'prod') ? 'https://www.ebi.ac.uk/ena/submit/webin-v2' : 'https://wwwdev.ebi.ac.uk/ena/submit/webin-v2'))

    Channel
      .fromPath(params.manifest)
      .splitCsv(header:true, sep:'\t')
      .map { row ->
          def f = file(row.file_path)
          if (!f.exists()) { throw new RuntimeException("Missing file: ${row.file_path}") }
          if (!row.assembly_accession) { throw new RuntimeException("Manifest row missing assembly_accession for file: ${row.file_path}") }
          def id = row.analysis_alias ?: f.getBaseName()
          def release = params.release.toString()
          def prj_alias = ("prj_${row.assembly_accession}_${release}").replaceAll('[^A-Za-z0-9._-]', '_')
          // Use derived project alias as a refname (not accession). XML generator will place it as refname, not accession.
          row.study = prj_alias
          def meta = [ id:id, remote_name: (row.remote_name ?: f.getName()), project_alias: prj_alias, assembly: row.assembly_accession, release: release ]
          tuple(meta, row, f)
      }
      .set { inputs }

    // Single approach: derive unique Projects from manifest (assembly + release), always
    def ch_proj_rows = inputs
        .map { meta, row, f ->
            def alias = meta.project_alias
            def title = "Annotation evidence project for ${meta.assembly}, ${meta.release}"
            def description = params.project_description ?: ''
            // Key by alias, carry a single meta map
            tuple(alias, [ alias: alias, name: alias, title: title, description: description, hold_until: params.hold_until ?: '' ])
        }
        .groupTuple()
        .map { alias, metas -> metas[0] } // one meta per alias

    ENA_GENERATE_PROJECT_XML(ch_proj_rows)
    ENA_SUBMIT_PROJECT(ENA_GENERATE_PROJECT_XML.out.xml, ch_webin_base, ch_webin_user, ch_webin_password)
    def ch_proj_files = ENA_SUBMIT_PROJECT.out.queued.map { meta, f -> f }

    // Compute md5 per file
    md5s = ENA_COMPUTE_MD5( inputs )

    ENA_FTP_UPLOAD(
        md5s,
        params.remote_dir ?: '',
        params.webin_ftp_host ?: 'webin2.ebi.ac.uk',
        ch_webin_user,
        ch_webin_password
        )

    // Generate XMLs
    ENA_GENERATE_XML(
        ENA_FTP_UPLOAD.out.uploaded,
        params.remote_dir ?: '',
        params.hold_until ?: ''
        )

    // Submit to async queue — one POST per file, returns immediately with a submission ID
    ENA_SUBMIT_ANALYSIS(
        ENA_GENERATE_XML.out.xml,
        ch_webin_base,
        ch_webin_user,
        ch_webin_password
        )
    def ch_analysis_files = ENA_SUBMIT_ANALYSIS.out.queued.map { meta, f -> f }

    // Collect all queue responses then poll until each submission resolves
    // Merge project and analysis queue files and poll all
    def ch_all_queues = ch_proj_files.mix(ch_analysis_files).collect()

    ENA_POLL_WEBIN(
        ch_all_queues,
        ch_webin_user,
        ch_webin_password,
        Channel.value(params.poll_interval    ?: 20),
        Channel.value(params.poll_max_attempts ?: 30)
        )

    ENA_POLL_WEBIN.out.accessions.view { f -> "Accessions written to: ${f}" }
}
