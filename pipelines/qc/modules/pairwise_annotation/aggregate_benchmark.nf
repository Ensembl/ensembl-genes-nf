process AGGREGATE_PAIRWISE_BENCHMARK {
    tag { 'pairwise_benchmark' }
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_benchmark" }, mode: 'copy', pattern: '*.tsv'
    publishDir { "${params.outdir}/qc/pairwise_benchmark" }, mode: 'copy', pattern: 'benchmark_report.md'

    input:
    path model_files
    path comparison_files
    path whole_model_files
    path summary_files
    path status_files

    output:
    path 'models_long.tsv', emit: models
    path 'structural_comparisons_long.tsv', emit: comparisons
    path 'whole_model_comparisons.tsv', emit: whole_models
    path 'benchmark_manifest.tsv', emit: manifest
    path 'tool_technology_summary.tsv', emit: summary
    path 'metric_definitions.tsv', emit: definitions
    path 'metric_status_long.tsv', emit: metric_status
    path 'benchmark_report.md', emit: report
    path 'benchmark_products_complete.txt', emit: complete
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    build_benchmark_products.py aggregate \
      --models ${model_files} \
      --comparisons ${comparison_files} \
      --whole-models ${whole_model_files} \
      --summaries ${summary_files} \
      --status ${status_files}
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.1.0
    END_VERSIONS
    """

    stub:
    """
    printf 'pair_id\tassembly\tspecies\tsample\treplicate\tlibrary_protocol\tsource_a\ttool_a\ttool_version_a\ttechnology_a\tformat_a\tsource_b\ttool_b\ttool_version_b\ttechnology_b\tformat_b\tcomparison_scope\ttruth_source\tsource_a_model_count\tsource_b_model_count\tsource_a_locus_count\tsource_b_locus_count\tsource_a_exon_count\tsource_b_exon_count\tn_intron_comparisons\tn_exact_intron_chains\tn_partial_intron_chains\tintron_chain_concordance\tn_end_comparisons\tn_ends_within_threshold\tend_concordance\tn_utr_comparisons\tn_exact_utr_matches\tutr_concordance\n' > benchmark_manifest.tsv
    printf 'technology_a\ttool_a\ttechnology_b\ttool_b\tn_pairs\n' > tool_technology_summary.tsv
    printf 'metric\tscope\tnumerator\tdenominator\ttruth_required\tdescription\n' > metric_definitions.tsv
    printf 'pair_id\tmetric\tstatus\treason\n' > metric_status_long.tsv
    printf 'pair_id\tassembly\tspecies\tsample\treplicate\tlibrary_protocol\tsource_a\ttool_a\ttool_version_a\ttechnology_a\tsource_b\ttool_b\ttool_version_b\ttechnology_b\tcomparison_scope\tlocus_id\tsource_a_transcript_id\tsource_b_transcript_id\tintron_chain_exact\tfive_prime_within_threshold\tthree_prime_within_threshold\tutr_exact_match\twhole_model_class\n' > whole_model_comparisons.tsv
    printf 'pair_id\tassembly\tspecies\tsample\treplicate\tlibrary_protocol\tsource\ttool\ttool_version\ttechnology\tsource_format\tmodel_role\ttranscript_id\toriginal_transcript_id\toriginal_gene_id\tlocus_id\tgene_assignment_type\tlocus_assignment_method\tmodel_status\tcontig\tstrand\ttranscript_start\ttranscript_end\texon_count\texon_coordinates\tintron_chain\tfive_prime_end\tthree_prime_end\tutr_coordinates\tcds_coordinates\tis_synthetic_locus\n' > models_long.tsv
    printf 'pair_id\tassembly\tspecies\tsample\treplicate\tlibrary_protocol\tsource_a\ttool_a\ttool_version_a\ttechnology_a\tsource_b\ttool_b\ttool_version_b\ttechnology_b\tcomparison_scope\tcomparison_type\tlocus_id\tfeature\tsource_a_transcript_id\tsource_b_transcript_id\tclassification\tshared_introns\tjaccard\tdistance_bp\tsource_a_utr_available\tsource_b_utr_available\texact_match\n' > structural_comparisons_long.tsv
    touch benchmark_products_complete.txt
    printf '# Long-read transcript reconstruction benchmark report\n\nStub run; no benchmark results were calculated.\n' > benchmark_report.md
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.1.0
    END_VERSIONS
    """
}
