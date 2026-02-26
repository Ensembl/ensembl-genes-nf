nextflow.enable.dsl=2

include { ENA_COMPUTE_MD5 }   from '../modules/compute_md5.nf'
include { ENA_FTP_UPLOAD }    from '../modules/upload_ftp.nf'
include { ENA_GENERATE_XML }  from '../modules/generate_xml.nf'
include { ENA_SUBMIT_WEBIN }  from '../modules/submit_webin.nf'
include { ENA_POLL_WEBIN }    from '../modules/poll_webin.nf'
include { ENA_GENERATE_PROJECT_XML } from '../modules/generate_project_xml.nf'

// Parse manifest and dispatch one analysis per row/file
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
          def f = file(row.file_path)
          if (!f.exists()) { throw new RuntimeException("Missing file: ${row.file_path}") }
          def id = row.analysis_alias ?: f.getBaseName()
          def meta = [ id:id, remote_name: (row.remote_name ?: f.getName()) ]
          tuple(meta, row, f)
      }
      .set { inputs }

    // Optional project registration: read a projects.tsv if provided
    def ch_proj_queue_list = null
    if (params.projects) {
        Channel
          .fromPath(params.projects)
          .splitCsv(header:true, sep:'\t')
          .map { row ->
              def meta = [
                id:          row.alias,  // satisfy submit module's tag/filename
                alias:       row.alias,
                name:        row.name ?: row.alias,
                title:       row.title ?: row.name ?: row.alias,
                description: row.description ?: '',
                hold_until:  params.hold_until ?: ''
              ]
              tuple(meta)
          }
          .set { project_rows }

        ENA_GENERATE_PROJECT_XML(project_rows)
        ENA_SUBMIT_WEBIN(ENA_GENERATE_PROJECT_XML.out.xml, ch_webin_base, ch_webin_user, ch_webin_password)
        ch_proj_queue_list = ENA_SUBMIT_WEBIN.out.queued.map { meta, f -> f }.collect()
    }

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
    ENA_SUBMIT_WEBIN(
        ENA_GENERATE_XML.out.xml,
        ch_webin_base,
        ch_webin_user,
        ch_webin_password
        )
    def ch_analysis_queue_list = ENA_SUBMIT_WEBIN.out.queued.map { meta, f -> f }.collect()

    // Collect all queue responses then poll until each submission resolves
    def ch_all_queue_list = ch_analysis_queue_list
    if (ch_proj_queue_list) {
        ch_all_queue_list = Channel.combine(ch_proj_queue_list, ch_analysis_queue_list).map { a, b -> a + b }
    }
    ENA_POLL_WEBIN(
        ch_all_queue_list,
        ch_webin_user,
        ch_webin_password,
        Channel.value(params.poll_interval    ?: 20),
        Channel.value(params.poll_max_attempts ?: 30)
        )

    ENA_POLL_WEBIN.out.accessions.view { f -> "Accessions written to: ${f}" }
}
