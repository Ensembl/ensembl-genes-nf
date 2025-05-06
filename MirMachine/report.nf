nextflow.enable.dsl = 2

include { COLLATE_RESULTS } from './modules/local/collate_results.nf'
include { GENERATE_SINGLE_SCORES } from './modules/local/score_single_heatmap.nf'


workflow {
    Channel
        .fromPath("${params.outdir}/mirmachine/*/*.heatmap.csv", glob: true)
        .set { heatmaps }

    GENERATE_SINGLE_SCORES(heatmaps)

    all_scores = GENERATE_SINGLE_SCORES.out.collect()
    

}
