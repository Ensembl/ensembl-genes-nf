nextflow.enable.dsl = 2

include { COLLATE_RESULTS } from './modules/local/collate_results.nf'
include { GENERATE_SCORES } from './modules/local/generate_scores.nf'


workflow {
    Channel
        .fromPath("${params.outdir}/mirmachine/*/*.heatmap.csv", glob: true)
        .view()
        .unique { it }
        .collect()
        .set { heatmaps }

    COLLATE_RESULTS(heatmaps)
    GENERATE_SCORES(COLLATE_RESULTS.out.heatmap, COLLATE_RESULTS.out.metadata)
}
