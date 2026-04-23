
include { AGAT_METRICS } from './subworkflows/agat/agat_stats.nf'

def validate_params() {
    def errors = []
    def repoRoot = params.ensembl_genes_repo ? file(params.ensembl_genes_repo) : null
    def defaultFeatureLevelsPath = "${repoRoot}/src/python/ensembl/genes/annotation-qc/metrics/config/feature_levels.yaml"

    if (!params.gff_csv)
        errors << "  --gff_csv is required"
    else if (!file(params.gff_csv).exists())
        errors << "  --gff_csv does not exist: ${params.gff_csv}"

    if (!params.ensembl_genes_repo)
        errors << "  --ensembl_genes_repo is required"
    else if (!file(params.ensembl_genes_repo).exists())
        errors << "  --ensembl_genes_repo does not exist: ${params.ensembl_genes_repo}"

    if (params.feature_levels) {
        if (!file(params.feature_levels).exists())
            errors << "  --feature_levels does not exist: ${params.feature_levels}"
    }
    else if (repoRoot) {
        def defaultFeatureLevels = file(defaultFeatureLevelsPath)
        if (!defaultFeatureLevels.exists())
            errors << "  feature_levels.yaml not found under --ensembl_genes_repo; pass --feature_levels explicitly or point --ensembl_genes_repo to a checkout containing it"
    }

    if (params.agat_parser) {
        if (!file(params.agat_parser).exists())
            errors << "  --agat_parser does not exist: ${params.agat_parser}"
    }
    else if (repoRoot) {
        def defaultParser = file("${repoRoot}/src/python/ensembl/genes/annotation-qc/parsers/parse_agat.py")
        if (!defaultParser.exists())
            errors << "  parse_agat.py not found under --ensembl_genes_repo; pass --agat_parser explicitly or point --ensembl_genes_repo to a checkout containing it"
    }

    if (errors) {
        log.error "Missing or invalid parameters:\n${errors.join('\n')}"
        System.exit(1)
    }
}

workflow {

    validate_params()

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

    /*
     * Run AGAT metrics in parallel for all CSV rows
     */
    if (params.run_agat_metrics) {
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
}
