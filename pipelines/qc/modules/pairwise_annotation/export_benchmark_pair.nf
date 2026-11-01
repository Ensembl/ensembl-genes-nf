process EXPORT_BENCHMARK_PAIR {
    tag { meta.pair_id }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_benchmark/${meta.pair_id}" }, mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), path(source_a_models), path(source_b_models), path(intron), path(ends), path(utr), path(structural), path(disagreements), path(evidence_a), path(evidence_b)

    output:
    tuple val(meta), path("${meta.pair_id}.models.tsv"), emit: models
    tuple val(meta), path("${meta.pair_id}.comparisons.tsv"), emit: comparisons
    tuple val(meta), path("${meta.pair_id}.whole_models.tsv"), emit: whole_models
    tuple val(meta), path("${meta.pair_id}.summary.tsv"), emit: summary
    path 'versions.yml', emit: versions

    script:
    def metadata = groovy.json.JsonOutput.toJson(meta)
    def encodedMetadata = metadata.bytes.encodeBase64().toString()
    """
    set -euo pipefail
    printf '%s' '${encodedMetadata}' | base64 -d > benchmark_metadata.json
    build_benchmark_products.py pair \
      --metadata benchmark_metadata.json \
      --source-a '${source_a_models}' --source-b '${source_b_models}' \
      --evidence-a '${evidence_a}' --evidence-b '${evidence_b}' \
      --intron '${intron}' --ends '${ends}' --utr '${utr}' --structural '${structural}' \
      --output-models '${meta.pair_id}.models.tsv' \
      --output-comparisons '${meta.pair_id}.comparisons.tsv' \
      --output-whole-models '${meta.pair_id}.whole_models.tsv' \
      --output-summary '${meta.pair_id}.summary.tsv'
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.1.0
    END_VERSIONS
    """

    stub:
    def stubRows = "${meta.pair_id}\t${meta.assembly}\t${meta.species ?: ''}\t${meta.sample ?: ''}\t${meta.replicate ?: ''}\t${meta.library_protocol ?: ''}\t${meta.source_a}\t${meta.tool_a ?: meta.source_a}\t${meta.tool_version_a ?: ''}\t${meta.technology_a ?: ''}\t${meta.format_a ?: ''}\t${meta.source_b}\t${meta.tool_b ?: meta.source_b}\t${meta.tool_version_b ?: ''}\t${meta.technology_b ?: ''}\t${meta.format_b ?: ''}\tpairwise_concordance\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t0\t\t0\t0\t\t0\t0\t\t0\t0\t\t"
    """
    printf 'pair_id\tassembly\tspecies\tsample\treplicate\tlibrary_protocol\tsource_a\ttool_a\ttool_version_a\ttechnology_a\tformat_a\tsource_b\ttool_b\ttool_version_b\ttechnology_b\tformat_b\tcomparison_scope\ttruth_source\tsource_a_model_count\tsource_b_model_count\tsource_a_locus_count\tsource_b_locus_count\tsource_a_exon_count\tsource_b_exon_count\tn_intron_comparisons\tn_exact_intron_chains\tn_partial_intron_chains\tintron_chain_concordance\tn_end_comparisons\tn_ends_within_threshold\tend_concordance\tn_utr_comparisons\tn_exact_utr_matches\tutr_concordance\n' > '${meta.pair_id}.summary.tsv'
    printf '${stubRows}\n' >> '${meta.pair_id}.summary.tsv'
    printf 'pair_id\tassembly\tspecies\tsample\treplicate\tlibrary_protocol\tsource\ttool\ttool_version\ttechnology\tsource_format\tmodel_role\ttranscript_id\toriginal_transcript_id\toriginal_gene_id\tlocus_id\tgene_assignment_type\tlocus_assignment_method\tmodel_status\tcontig\tstrand\ttranscript_start\ttranscript_end\texon_count\texon_coordinates\tintron_chain\tfive_prime_end\tthree_prime_end\tutr_coordinates\tcds_coordinates\tis_synthetic_locus\tread_count\tjunction_support\tfull_length_support\tevidence_source_tool\tevidence_available\n' > '${meta.pair_id}.models.tsv'
    printf 'pair_id\tassembly\tspecies\tsample\treplicate\tlibrary_protocol\tsource_a\ttool_a\ttool_version_a\ttechnology_a\tsource_b\ttool_b\ttool_version_b\ttechnology_b\tcomparison_scope\tcomparison_type\tlocus_id\tfeature\tsource_a_transcript_id\tsource_b_transcript_id\tclassification\tshared_introns\tjaccard\tdistance_bp\tsource_a_utr_available\tsource_b_utr_available\texact_match\n' > '${meta.pair_id}.comparisons.tsv'
    printf 'pair_id\tassembly\tspecies\tsample\treplicate\tlibrary_protocol\tsource_a\ttool_a\ttool_version_a\ttechnology_a\tsource_b\ttool_b\ttool_version_b\ttechnology_b\tcomparison_scope\tlocus_id\tsource_a_transcript_id\tsource_b_transcript_id\tintron_chain_exact\tfive_prime_within_threshold\tthree_prime_within_threshold\tutr_exact_match\twhole_model_class\n' > '${meta.pair_id}.whole_models.tsv'
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.1.0
    END_VERSIONS
    """
}
