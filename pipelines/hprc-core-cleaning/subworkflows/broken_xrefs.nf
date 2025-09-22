include { CREATE_HGNC_XREFS } from '../modules/create_hgnc_xrefs.nf'
include { UPDATE_GENE_DISPLAY_XREF_IDS } from '../modules/update_display_xrefs.nf'
include { VERIFY_XREF_FIX } from '../modules/check_xref_fix.nf'

workflow FIX_GENE_DISPLAY_XREFS {

    take:
    core_info

    main:
    CREATE_HGNC_XREFS(core_info)
    UPDATE_GENE_DISPLAY_XREF_IDS(CREATE_HGNC_XREFS.out.created_xrefs)
    VERIFY_XREF_FIX(UPDATE_GENE_DISPLAY_XREF_IDS.out.updated_genes)

    emit:
    verification = VERIFY_XREF_FIX.out.verification
    xref_logs = CREATE_HGNC_XREFS.out.created_xrefs
    update_logs = UPDATE_GENE_DISPLAY_XREF_IDS.out.updated_genes

}