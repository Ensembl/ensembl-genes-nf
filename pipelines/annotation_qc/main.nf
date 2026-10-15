nextflow.enable.dsl=2

//include { INPUT_RESOLVER } from './modules/input_resolver'
include { PARSE_ANNOTATION } from './modules/parse_annotation'
include { GENOME_METRICS } from './modules/genome_metrics'
include { FEATURE_METRICS } from './modules/feature_metrics'
include { FASTA_METRICS } from './modules/fasta_metrics'

workflow {
    // annotation = INPUT_RESOLVER()

    // test:
    annotation = Channel.of(
        tuple(
            'GCA_004011555.1',
            file('/Users/ereboperezsilva/Documents/maintenance_2026.10.13/nf_qc/carp_reference/GCA_004011555.1/genes.gff3.gz'),
            file('/Users/ereboperezsilva/Documents/maintenance_2026.10.13/nf_qc/carp_reference/GCA_004011555.1/softmasked.fa.gz')
        )
    )

    parsed = PARSE_ANNOTATION(annotation)

    fasta_metrics = FASTA_METRICS(parsed)

    feature_input = parsed.map { gca_accession, normalised, _fasta ->
        tuple(gca_accession, normalised)
    }
    feature_metrics = FEATURE_METRICS(feature_input)

    genome_input = parsed
        .map { gca_accession, normalised, _fasta -> tuple(gca_accession, normalised) }
        .join(fasta_metrics, by: 0)
        .map { gca_accession, normalised, fasta_metrics_csv, _masked_overlaps ->
            tuple(gca_accession, normalised, fasta_metrics_csv)
        }
    genome_metrics = GENOME_METRICS(genome_input)

    // Output results
}
