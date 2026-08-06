#!/usr/bin/env nextflow

include { validateParameters; samplesheetToList } from 'plugin/nf-schema'
include { AGAT_METRICS } from './subworkflows/agat/agat_stats.nf'
include { INTERPRO_SCAN } from './subworkflows/interpro/interproscan.nf'

def validate_params() {
    def errors = []
    def repoRoot = params.ensembl_genes_repo ? file(params.ensembl_genes_repo) : null

    if (params.run_agat_metrics && repoRoot) {
        def defaultFeatureLevels = file("${repoRoot}/src/python/ensembl/genes/annotation_qc/config/feature_levels.yaml")
        if (!params.feature_levels && !defaultFeatureLevels.exists())
            errors << "feature_levels.yaml not found under --ensembl_genes_repo; pass --feature_levels explicitly"

        def defaultAgatParser = file("${repoRoot}/src/python/ensembl/genes/annotation_qc/parsers/parse_agat.py")
        if (!params.agat_parser && !defaultAgatParser.exists())
            errors << "parse_agat.py not found under --ensembl_genes_repo; pass --agat_parser explicitly"
    }

    if (params.run_interproscan && repoRoot && !params.interpro_parser) {
        def defaultInterproParser = file("${repoRoot}/src/python/ensembl/genes/annotation_qc/parsers/interpro.py")
        if (!defaultInterproParser.exists())
            errors << "interpro.py not found under --ensembl_genes_repo; pass --interpro_parser explicitly"
    }

    if (errors) {
        error "Missing or invalid parameters:\n${errors.join('\n')}"
    }
}

def samplesheet_schema() {
    if (params.run_agat_metrics && params.run_interproscan)
        return 'assets/samplesheet_agat_interpro.json'
    if (params.run_agat_metrics)
        return 'assets/samplesheet_agat.json'
    return 'assets/samplesheet_interpro.json'
}

workflow {

    validateParameters()
    println "run_agat_metrics = ${params.run_agat_metrics}"
    println "run_interproscan = ${params.run_interproscan}"
    validate_params()

    samples = samplesheetToList(params.input_csv, samplesheet_schema())

    Channel
        .fromList(samples)
        .map { row ->
            def sample = row[0]
            def gff_path = params.run_agat_metrics
                ? file(row[1])
                : null
            def protein_index = params.run_agat_metrics ? 2 : 1
            def protein_path = params.run_interproscan
                ? file(row[protein_index])
                : null

            def meta = [
                id     : sample,
                sample : sample
            ]

            if (gff_path)
                meta.gff3 = gff_path.name

            if (protein_path)
                meta.protein = protein_path.name

            tuple(meta, gff_path, protein_path)
        }
        .set { sample_ch }

    gff_ch = sample_ch
        .filter { meta, gff_path, protein_path -> gff_path }
        .map { meta, gff_path, protein_path -> tuple(meta, gff_path) }

    protein_ch = sample_ch
        .filter { meta, gff_path, protein_path -> protein_path }
        .map { meta, gff_path, protein_path -> tuple(meta, protein_path) }

    // Singletons
    ensembl_genes_repo = file(params.ensembl_genes_repo)
    agat_parser = params.agat_parser ? file(params.agat_parser) : ''
    interpro_parser = params.interpro_parser \
        ? file(params.interpro_parser) \
        : ''

    // Optional: feature levels YAML
    feature_levels_yaml = params.feature_levels \
        ? file(params.feature_levels) \
        : file("${ensembl_genes_repo}/src/python/ensembl/genes/annotation_qc/config/feature_levels.yaml")

    // Run AGAT metrics
    if (params.run_agat_metrics) {
        agat = AGAT_METRICS(
            gff_ch,
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
            protein_ch,
            params.database,
            params.data_file_path,
            params.ensembl_genes_repo,
            interpro_parser
        )
    }
}
