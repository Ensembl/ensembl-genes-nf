nextflow.enable.dsl=2

include { INPUT_RESOLVER } from './modules/input_resolver'
include { PARSE_ANNOTATION } from './modules/parse_annotation'
include { GENOME_METRICS } from './modules/genome_metrics'
include { FEATURE_METRICS } from './modules/feature_metrics'
include { FASTA_METRICS } from './modules/fasta_metrics'

workflow {

}

