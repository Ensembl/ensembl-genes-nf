#!/usr/bin/env nextflow

include { validateParameters } from 'plugin/nf-schema'
include { AGAT_METRICS } from './subworkflows/agat/agat_stats.nf'
include { PAIRWISE_ANNOTATION_COMPARISON } from './subworkflows/pairwise_annotation/pairwise_annotation_comparison.nf'

def is_url(value) {
    value ==~ /(?i)^(https?|ftp):\/\/.+/
}

def is_missing(value) {
    value == null || value.toString().trim() in ['', 'NA', 'null']
}

def validate_params() {
    def errors = []
    def repoRoot = params.ensembl_genes_repo ? file(params.ensembl_genes_repo) : null
    def defaultFeatureLevelsPath = "${repoRoot}/src/python/ensembl/genes/annotation-qc/metrics/config/feature_levels.yaml"

    if (params.run_agat_metrics && !params.gff_csv)
        errors << "  --gff_csv is required"
    else if (params.run_agat_metrics && !file(params.gff_csv).exists())
        errors << "  --gff_csv does not exist: ${params.gff_csv}"

    if (params.run_agat_metrics && !params.ensembl_genes_repo)
        errors << "  --ensembl_genes_repo is required"
    else if (params.run_agat_metrics && !file(params.ensembl_genes_repo).exists())
        errors << "  --ensembl_genes_repo does not exist: ${params.ensembl_genes_repo}"

    if (params.run_agat_metrics && params.feature_levels) {
        if (!file(params.feature_levels).exists())
            errors << "  --feature_levels does not exist: ${params.feature_levels}"
    }
    else if (params.run_agat_metrics && repoRoot) {
        def defaultFeatureLevels = file(defaultFeatureLevelsPath)
        if (!defaultFeatureLevels.exists())
            errors << "  feature_levels.yaml not found under --ensembl_genes_repo; pass --feature_levels explicitly or point --ensembl_genes_repo to a checkout containing it"
    }

    if (params.run_agat_metrics && params.agat_parser) {
        if (!file(params.agat_parser).exists())
            errors << "  --agat_parser does not exist: ${params.agat_parser}"
    }
    else if (params.run_agat_metrics && repoRoot) {
        def defaultParser = file("${repoRoot}/src/python/ensembl/genes/annotation-qc/parsers/parse_agat.py")
        if (!defaultParser.exists())
            errors << "  parse_agat.py not found under --ensembl_genes_repo; pass --agat_parser explicitly or point --ensembl_genes_repo to a checkout containing it"
    }

    if (params.run_pairwise_annotation_comparison && !params.pairwise_csv)
        errors << "  --pairwise_csv is required when --run_pairwise_annotation_comparison is true"
    else if (params.run_pairwise_annotation_comparison && !file(params.pairwise_csv).exists())
        errors << "  --pairwise_csv does not exist: ${params.pairwise_csv}"

    if (params.pairwise_ensg_lookup && !file(params.pairwise_ensg_lookup).exists())
        errors << "  --pairwise_ensg_lookup does not exist: ${params.pairwise_ensg_lookup}"

    if (errors) {
        log.error "Missing or invalid parameters:\n${errors.join('\n')}"
        System.exit(1)
    }
}

workflow {

    validateParameters()
    validate_params()

    /*
     * Run AGAT metrics in parallel for all CSV rows
     */
    if (params.run_agat_metrics) {
        /*
         * Read CSV with columns: sample,gff3
         */
        Channel
            .fromPath(params.gff_csv, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                // row is a map: [sample: 'mouse_1', gff3: '/path/to/file.gff3', ...]
                def gff_path = file(row.gff3, checkIfExists: true)
                def meta = [
                    id     : row.sample,
                    sample : row.sample,
                    gff3   : gff_path.name
                ]
                tuple(meta, gff_path)
            }
            .set { gff_ch }

        /*
         * Singletons
         */
        ensembl_genes_repo = file(params.ensembl_genes_repo)
        agat_parser        = params.agat_parser ? file(params.agat_parser) : ''

        /*
         * Optional: feature levels YAML (single, shared for all samples)
         */
        feature_levels_yaml = params.feature_levels \
            ? file(params.feature_levels) \
            : file("${ensembl_genes_repo}/src/python/ensembl/genes/annotation-qc/metrics/config/feature_levels.yaml")

        agat = AGAT_METRICS(
            gff_ch,
            feature_levels_yaml,
            ensembl_genes_repo,
            agat_parser
        )

        // for now just log the outputs
        agat.genebuild_csv.view { meta, f ->
            "Genebuild metrics for ${meta.id}: ${f}"
        }
    }

    if (params.run_pairwise_annotation_comparison) {
        Channel
            .fromPath(params.pairwise_csv, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def sample = row.sample ?: row.id
                def sourceA = row.source_a ?: 'source_a'
                def sourceB = row.source_b ?: 'source_b'
                def sourceARef = row.gff_a ?: row.source_a_gff ?: row.annotation_a
                def sourceBRef = row.gff_b ?: row.source_b_gff ?: row.annotation_b
                def assemblyReport = row.assembly_report ?: ''

                if (is_missing(sample))
                    error "pairwise_csv row is missing sample/id"
                if (is_missing(sourceARef))
                    error "pairwise_csv row for ${sample} is missing gff_a/source_a_gff/annotation_a"
                if (is_missing(sourceBRef))
                    error "pairwise_csv row for ${sample} is missing gff_b/source_b_gff/annotation_b"

                [sourceARef, sourceBRef, assemblyReport].findAll { !is_missing(it) && !is_url(it) }.each { ref ->
                    if (!file(ref).exists())
                        error "pairwise_csv row for ${sample} references a missing local file: ${ref}"
                }

                def meta = [
                    id       : sample,
                    sample   : sample,
                    source_a : sourceA,
                    source_b : sourceB
                ]
                tuple(meta, sourceARef, sourceBRef, assemblyReport)
            }
            .set { pairwise_samples_ch }

        pairwise_ensg_lookup = params.pairwise_ensg_lookup ?: ''

        pairwise = PAIRWISE_ANNOTATION_COMPARISON(
            pairwise_samples_ch,
            pairwise_ensg_lookup
        )

        pairwise.rbh.view { meta, f ->
            "Pairwise RBH gene pairs for ${meta.id}: ${f}"
        }
    }
}
