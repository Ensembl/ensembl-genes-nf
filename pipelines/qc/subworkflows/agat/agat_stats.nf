
include { AGAT_RUN_STATS } from '../../modules/agat/run_agat_stats.nf'
include { STRIP_GFF_REGIONS } from '../../modules/agat/strip_gff_regions.nf'
include { AGAT_PARSE     } from '../../modules/agat/parse_agat.nf'

workflow AGAT_METRICS {

    take:
        // [meta, gff3]
        gff3_ch
        // feature_levels yaml (may be empty)
        feature_levels_yaml
        // path to ensembl-genes checkout
        ensembl_genes_repo
        // optional explicit parse_agat.py path
        agat_parser

    main:
        stripped_gff = STRIP_GFF_REGIONS(gff3_ch)

        agat_txt = AGAT_RUN_STATS(
            stripped_gff.cleaned_gff,
            feature_levels_yaml
        )

        genebuild = AGAT_PARSE(
            agat_txt.stats_txt,
            ensembl_genes_repo,
            agat_parser,
            feature_levels_yaml
        )

    emit:
        stats_txt     = agat_txt.stats_txt
        genebuild_csv = genebuild.genebuild_csv
        versions      = stripped_gff.versions.mix(agat_txt.versions).mix(genebuild.versions)
}
