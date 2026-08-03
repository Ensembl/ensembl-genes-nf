#!/usr/bin/env nextflow

include { validateParameters } from 'plugin/nf-schema'
include { AGAT_METRICS } from './subworkflows/agat/agat_stats.nf'
include { INTERPRO_SCAN } from './subworkflows/interpro/interproscan.nf'

def validate_params() {
    def errors = []
    def repoRoot = params.ensembl_genes_repo ? file(params.ensembl_genes_repo) : null
    def defaultFeatureLevelsPath = "${repoRoot}/src/python/ensembl/genes/annotation-qc/metrics/config/feature_levels.yaml"

    if (!params.input_csv)
        errors << "  --input_csv is required"
    else if (!file(params.input_csv).exists())
        errors << "  --input_csv does not exist: ${params.input_csv}"

    if (!params.ensembl_genes_repo)
        errors << "  --ensembl_genes_repo is required"
    else if (!file(params.ensembl_genes_repo).exists())
        errors << "  --ensembl_genes_repo does not exist: ${params.ensembl_genes_repo}"

    if (params.run_agat_metrics) {

        if (params.feature_levels) {
            if (!file(params.feature_levels).exists())
                errors << "--feature_levels does not exist: ${params.feature_levels}"
        }
        else if (repoRoot) {
            def defaultFeatureLevels = file(defaultFeatureLevelsPath)
            if (!defaultFeatureLevels.exists())
                errors << "feature_levels.yaml not found under --ensembl_genes_repo; pass --feature_levels explicitly or point --ensembl_genes_repo to a checkout containing it"
        }

        if (params.agat_parser) {
            if (!file(params.agat_parser).exists())
                errors << "--agat_parser does not exist: ${params.agat_parser}"
        }
    }

    def allowed_dbs = [
        'TIGRFAM',
        'SFLD',
        'SUPERFAMILY',
        'PANTHER',
        'Gene3D',
        'Hamap',
        'ProSiteProfiles',
        'Coils',
        'SMART',
        'CDD',
        'PRINTS',
        'PIRSR',
        'ProSitePatterns',
        'AntiFam',
        'Pfam',
        'MobiDBLite',
        'PIRSF'
    ]

    if (params.run_interproscan) {
        if (params.data_file_path) {
            if (!file(params.data_file_path).exists())
                errors << "  --data_file_path does not exist: ${params.data_file_path}"
        }

        if (params.interpro_parser) {
            if (!file(params.interpro_parser).exists())
                errors << "  --interpro_parser does not exist: ${params.interpro_parser}"
        }

        if (!(params.database in allowed_dbs)) {
            errors << "--database must be one of: ${allowed_dbs.join(', ')}"
        }
    }

    if (errors) {
        error "Missing or invalid parameters:\n${errors.join('\n')}"
    }
}

workflow {

    validateParameters()
    println "run_agat_metrics = ${params.run_agat_metrics}"
    println "run_interproscan = ${params.run_interproscan}"
    validate_params()

    /*
     * Read CSV with columns: sample,gff3,protein
     */
    Channel
        .fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->

            def meta = [
                id     : row.sample,
                sample : row.sample
            ]

            def gff_path = row.gff3 ? file(row.gff3, checkIfExists: true) : null
            def protein_path = row.protein ? file(row.protein, checkIfExists: true) : null

            if (gff_path)
                meta.gff3 = gff_path.name

            if (protein_path)
                meta.protein = protein_path.name

            tuple(meta, gff_path, protein_path)
        }
        .set { input_ch }

    // Singletons
    ensembl_genes_repo = file(params.ensembl_genes_repo)
    agat_parser = params.agat_parser ? file(params.agat_parser) : ''

    // Optional: feature levels YAML
    feature_levels_yaml = params.feature_levels \
        ? file(params.feature_levels) \
        : file("${ensembl_genes_repo}/src/python/ensembl/genes/annotation-qc/metrics/config/feature_levels.yaml")

    // Run AGAT metrics
    if (params.run_agat_metrics) {
        agat = AGAT_METRICS(
            input_ch,
            feature_levels_yaml,
            ensembl_genes_repo,
            agat_parser
        )

        agat.genebuild_csv.view { meta, f ->
            "Genebuild metrics for ${meta.id}: ${f}"
        }
    }

    // Run InterProScan
    if (params.run_interproscan) {
        interpro = INTERPRO_SCAN(
            input_ch,
            params.database,
            params.data_file_path,
            params.ensembl_genes_repo
        )
    }
}
