nextflow.enable.dsl=2

include { ENA_COMPUTE_MD5 }   from '../modules/compute_md5.nf'
include { ENA_ASCP_UPLOAD }   from '../modules/upload_ascp.nf'
include { ENA_FTP_UPLOAD }    from '../modules/upload_ftp.nf'
include { ENA_GENERATE_XML }  from '../modules/generate_xml.nf'
include { ENA_SUBMIT_WEBIN }  from '../modules/submit_webin.nf'

// Parse manifest and dispatch one analysis per row/file
workflow ENA_SUBMIT_WORKFLOW {
    assert params.manifest, "--manifest is required"

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

    // Compute md5 per file
    md5s = ENA_COMPUTE_MD5( inputs )

    // Uploads — broadcast config as value channels
    def ch_remote_dir = Channel.value(params.remote_dir ?: '')
    def ch_ascp_limit = Channel.value(params.ascp_limit ?: '300M')
    def ch_ascp_host  = Channel.value(params.webin_ascp_host ?: 'webin.ebi.ac.uk')
    def ch_ftp_host   = Channel.value(params.webin_ftp_host  ?: 'ftp.webin.ebi.ac.uk')

    def uploads = (params.upload_protocol == 'aspera') ?
        ENA_ASCP_UPLOAD(md5s, ch_remote_dir, ch_ascp_limit, ch_ascp_host) :
        ENA_FTP_UPLOAD(md5s, ch_remote_dir, ch_ftp_host)

    // Generate XMLs
    def ch_hold_until = Channel.value(params.hold_until ?: null)
    xmls = ENA_GENERATE_XML( uploads, ch_hold_until )

    // Submit via Webin — URLs provided by workflow
    def ch_submit_api  = Channel.value(params.submit_api ?: 'v1')
    def ch_v1_base     = Channel.value(params.webin_v1_base ?: ((params.mode == 'prod') ? 'https://www.ebi.ac.uk/ena/submit/drop-box/' : 'https://wwwdev.ebi.ac.uk/ena/submit/drop-box/'))
    def ch_v2_base     = Channel.value(params.webin_v2_base ?: ((params.mode == 'prod') ? 'https://www.ebi.ac.uk/ena/submit/webin-v2/'     : 'https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/'))
    receipts = ENA_SUBMIT_WEBIN( xmls, ch_submit_api, ch_v1_base, ch_v2_base )

    receipts.view { it -> "ENA receipt for ${it[0].id}: ${it[1]}" }
}
