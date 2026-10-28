include { AGAT_RUN_STATS } from '../../modules/agat/run_agat_stats.nf'
include { STRIP_GFF_REGIONS } from '../../modules/agat/strip_gff_regions.nf'
include { AGAT_PARSE } from '../../modules/agat/parse_agat.nf'

workflow AGAT_METRICS {
    take:
        // Channel of [meta, gff3] tuples.
        gff3_ch
        // Ensembl-owned AGAT feature definitions used by the statistics step.
        feature_levels_yaml

    main:
        stripped_gff = STRIP_GFF_REGIONS(gff3_ch)
        agat_txt = AGAT_RUN_STATS(
            stripped_gff.cleaned_gff,
            feature_levels_yaml
        )
        genebuild = AGAT_PARSE(
            agat_txt.stats_txt
        )

    emit:
        stats_txt = agat_txt.stats_txt
        genebuild_csv = genebuild.genebuild_csv
        versions = stripped_gff.versions.mix(agat_txt.versions).mix(genebuild.versions)
}
