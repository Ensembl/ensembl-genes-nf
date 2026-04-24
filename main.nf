#!/usr/bin/env nextflow
/*
========================================================================================
    ensembl-genes-nf  —  Full Annotation Master Pipeline
========================================================================================
    Orchestrates all annotation sub-pipelines end-to-end by running each as a child
    `nextflow run` process.  No eHive required for experimental/comparison runs.

    Each stage is a Nextflow process whose script block calls one sub-pipeline.
    Output paths flow between stages as strings derived from the known publishDir
    conventions of each sub-pipeline.  The annotation layers run in parallel;
    their GFF3 outputs are collected before consolidation.

    Usage (minimal — ab initio only):
      nextflow run . \
        --assembly_accession GCA_964261345.1 \
        --assembly_name      mHetGlaV3 \
        --outdir             /hps/scratch/.../hetgla \
        --nf_work_root       /hps/scratch/.../hetgla/nf_work \
        --ab_initio_species  human \
        -profile             cluster

    Usage (full evidence):
      nextflow run . \
        --assembly_accession GCA_964261345.1 \
        --assembly_name      mHetGlaV3 \
        --outdir             /hps/scratch/.../hetgla \
        --nf_work_root       /hps/scratch/.../hetgla/nf_work \
        --sample_sheet       /path/to/rnaseq.csv \
        --uniprot_fasta      /path/to/rodentia_uniprot.fa \
        --assembly_refseq_accession GCF_964261345.1 \
        --ab_initio_species  human \
        --repbase_library    rodentia \
        -profile             cluster

    Required:
      --assembly_accession   GCA accession (e.g. GCA_964261345.1)
      --assembly_name        Assembly name, no spaces (e.g. mHetGlaV3)
      --outdir               Root output directory
      --nf_work_root         Nextflow work-dir root for child sub-pipelines

    RNA-seq evidence (one of):
      --sample_sheet         CSV: id,fastq_1[,fastq_2][,strandedness]
      --rnaseq_bioproject    ENA BioProject accession (reads fetched automatically)
      --rnaseq_run_accessions Comma-separated SRR/ERR accessions

    Long-read evidence (one of):
      --long_read_sample_sheet  TSV with local paths
      --long_read_bioproject    ENA BioProject for long-read data

    Protein homology:
      --uniprot_fasta        Local clade UniProt FASTA
      --uniprot_taxon_id     NCBI taxon ID — fetch clade reviewed proteins automatically

    Other optional evidence:
      --cdna_fasta           cDNA FASTA for best_targeted
      --protein_fasta        Species proteins for best_targeted
      --rfam_cm              Rfam.cm for short_ncrna
      --mirna_fasta          miRNA FASTA for short_ncrna
      --source_fasta         Source genome for projection
      --source_gff3          Source annotation for projection
      --assembly_refseq_accession  GCF accession → refseq_import
      --repbase_library      RepeatMasker library name (default: vertebrates)
      --custom_repeat_library  Path to custom repeat FASTA
      --skip_repeatmodeler   Skip RepeatModeler de-novo library build (bool)
      --selenoprotein_fasta  FASTA for selenoprotein flagging

    Ab initio:
      --ab_initio_species    Augustus species model name (e.g. 'human', 'mus')
                             Use closest trained model for new species.

    Consolidation:
      --layer_priorities     JSON map of stage-name → priority integer
                             Lower = higher quality. Default derived from enabled stages.
                             e.g. '{"rnaseq":2,"ab_initio":5,"refseq_import":4,"genblast_homology":6}'

    Core DB loading (optional — gff3_to_core skipped if not supplied):
      --db_host, --db_port, --db_user, --db_password, --db_name
      --stable_id_prefix     e.g. '' for human, 'ENSHETG' for naked mole-rat
      --stable_id_start      Integer start for new stable IDs
      --species_name         Species scientific name (for core DB loading)
========================================================================================
*/

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Params defaults
// ---------------------------------------------------------------------------
params.assembly_accession        = null
params.assembly_name             = null
params.outdir                    = null
params.nf_work_root              = null

// RNA-seq evidence (local or ENA fetch)
params.sample_sheet              = null
params.rnaseq_bioproject         = null
params.rnaseq_run_accessions     = null

// Long-read evidence
params.long_read_sample_sheet    = null
params.long_read_bioproject      = null

// Protein homology
params.uniprot_fasta             = null
params.uniprot_taxon_id          = null

// Other optional evidence
params.cdna_fasta                = null
params.protein_fasta             = null
params.rfam_cm                   = null
params.mirna_fasta               = null
params.source_fasta              = null
params.source_gff3               = null
params.assembly_refseq_accession = null
params.repbase_library           = 'vertebrates'
params.custom_repeat_library     = null
params.skip_repeatmodeler        = false
params.selenoprotein_fasta       = null

// Ab initio
params.ab_initio_species         = 'human'   // closest Augustus model for new species

// Consolidation priority map — lower integer = higher quality evidence
// Default is set at runtime based on which stages ran; override with explicit JSON.
params.layer_priorities          = null       // JSON string or null → auto-generated

// Core DB loading
params.db_host                   = null
params.db_port                   = 3306
params.db_user                   = null
params.db_password               = ''
params.db_name                   = null
params.stable_id_prefix          = ''
params.stable_id_start           = 1
params.species_id                = 1
params.analysis_logic_name       = 'ensembl'
params.coord_system              = 'chromosome'
params.assembly_version          = null
params.species_name              = null

// Profile propagated to child sub-pipelines
params.nextflow_profile          = 'cluster'

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------
def validate_params() {
    def errors = []
    if (!params.assembly_accession) errors << "  --assembly_accession is required"
    if (!params.assembly_name)      errors << "  --assembly_name is required"
    if (!params.outdir)             errors << "  --outdir is required"
    if (!params.nf_work_root)       errors << "  --nf_work_root is required"
    if (errors) {
        log.error "Missing required parameters:\n${errors.join('\n')}"
        System.exit(1)
    }

    def has_rnaseq     = params.sample_sheet || params.rnaseq_bioproject || params.rnaseq_run_accessions
    def has_longread   = params.long_read_sample_sheet || params.long_read_bioproject
    def has_homology   = params.uniprot_fasta || params.uniprot_taxon_id
    def has_projection = params.source_fasta && params.source_gff3
    def has_targeted   = params.cdna_fasta || params.protein_fasta
    def has_ncrna      = params.rfam_cm
    def has_refseq     = params.assembly_refseq_accession

    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║       ensembl-genes-nf  Full Annotation Pipeline     ║
    ╠══════════════════════════════════════════════════════╣
    ║  Assembly  : ${params.assembly_accession} (${params.assembly_name})
    ║  Outdir    : ${params.outdir}
    ║  Work root : ${params.nf_work_root}
    ║
    ║  Annotation layers:
    ║    Ab initio      : ${params.ab_initio_species} (always)
    ║    RNA-seq        : ${has_rnaseq       ? (params.rnaseq_bioproject ?: params.sample_sheet ?: params.rnaseq_run_accessions) : 'SKIP'}
    ║    Long-read      : ${has_longread     ? (params.long_read_bioproject ?: params.long_read_sample_sheet) : 'SKIP'}
    ║    Protein homol. : ${has_homology     ? (params.uniprot_taxon_id ?: params.uniprot_fasta) : 'SKIP'}
    ║    cDNA targeted  : ${has_targeted     ? 'enabled' : 'SKIP'}
    ║    ncRNA          : ${has_ncrna        ? 'enabled' : 'SKIP'}
    ║    Projection     : ${has_projection   ? 'enabled' : 'SKIP'}
    ║    RefSeq import  : ${has_refseq       ? params.assembly_refseq_accession : 'SKIP'}
    ║
    ║  Core loading    : ${params.db_name ?: 'SKIP'}
    ╚══════════════════════════════════════════════════════╝
    """.stripIndent()
}

// ---------------------------------------------------------------------------
// Helpers — predictable publishDir paths from known sub-pipeline conventions
// ---------------------------------------------------------------------------
def assembly_id() { "${params.assembly_accession}_${params.assembly_name}" }
def stage_outdir(String stage) { "${params.outdir}/${stage}" }
def nf_work(String stage)      { "${params.nf_work_root}/${stage}" }
def pipelines_dir()            { "${projectDir}/pipelines" }

/** Build the base `nextflow run <pipeline>` command for a sub-pipeline. */
def nf_run(String stage, List extra_args = []) {
    ([
        "nextflow run ${pipelines_dir()}/${stage}/main.nf",
        "-profile ${params.nextflow_profile}",
        "-work-dir ${nf_work(stage)}",
        "-resume",
        "--outdir ${stage_outdir(stage)}",
    ] + extra_args).join(" \\\n    ")
}

/**
 * Extract the GFF3 output path from a sub-pipeline's output_manifest.json.
 * Fails loudly with a clear error if the manifest is missing or has no GFF3.
 */
def extract_gff3_from_manifest(String manifest_path) {
    """
    python3 - <<'PYEOF'
import json, sys, os
manifest = '${manifest_path}'
if not os.path.exists(manifest):
    print(f"ERROR: output_manifest.json not found at {manifest}", file=sys.stderr)
    sys.exit(1)
with open(manifest) as fh:
    data = json.load(fh)
gff3_paths = [o['path'] for o in data.get('outputs', []) if o.get('type') == 'gff3']
if not gff3_paths:
    print(f"ERROR: no GFF3 output recorded in {manifest}", file=sys.stderr)
    print(f"Manifest contents: {json.dumps(data, indent=2)}", file=sys.stderr)
    sys.exit(1)
print(gff3_paths[0], end='')
PYEOF
    """
}

// ---------------------------------------------------------------------------
// Default layer priorities — lower = higher quality / confidence
// ---------------------------------------------------------------------------
def default_layer_priorities(Map stages) {
    def prio = [:]
    if (stages.long_read)   prio['long_read']   = 0
    if (stages.targeted)    prio['best_targeted'] = 1
    if (stages.rnaseq)      prio['rnaseq']       = 2
    if (stages.projection)  prio['projection']   = 3
    if (stages.refseq)      prio['refseq_import'] = 4
    prio['ab_initio'] = 5
    if (stages.genblast)    prio['genblast_homology'] = 6
    if (stages.ncrna)       prio['short_ncrna']  = 7
    if (stages.igtr)        prio['igtr']         = 7
    return groovy.json.JsonOutput.toJson(prio)
}

// ---------------------------------------------------------------------------
// Process definitions
// ---------------------------------------------------------------------------

process LOAD_ASSEMBLY {
    label 'process_launcher'
    output:
    val "${stage_outdir('load_assembly')}", emit: outdir

    script:
    """
    ${nf_run('load_assembly', [
        "--assembly_accession ${params.assembly_accession}",
        "--assembly_name      ${params.assembly_name}",
    ])}
    """

    stub:
    def d = stage_outdir('load_assembly')
    """
    mkdir -p ${d}/genome
    touch ${d}/genome/${assembly_id()}_genomic.fna
    touch ${d}/genome/${assembly_id()}_genomic.fna.fai
    touch ${d}/genome/${assembly_id()}.synonyms.tsv
    echo '{"pipeline":"load_assembly","outputs":[{"type":"fasta","path":"${d}/genome/${assembly_id()}_genomic.fna"}]}' > ${d}/output_manifest.json
    """
}

process REPEAT_MASKING {
    label 'process_launcher'
    input:
    val load_outdir

    output:
    val "${stage_outdir('repeat_masking')}", emit: outdir

    script:
    def genome_fasta = "${load_outdir}/genome/${assembly_id()}_genomic.fna"
    def lib_arg = params.custom_repeat_library
        ? "--custom_library     ${params.custom_repeat_library}"
        : "--repbase_library    ${params.repbase_library}"
    def skip_rm_arg = params.skip_repeatmodeler ? "--skip_repeatmodeler true" : ""
    """
    ${nf_run('repeat_masking', [
        "--genome_fasta ${genome_fasta}",
        lib_arg,
        skip_rm_arg,
    ])}
    """

    stub:
    def d = stage_outdir('repeat_masking')
    """
    mkdir -p ${d}/genome ${d}/repeats
    touch ${d}/genome/${assembly_id()}_genomic.softmasked.fa
    touch ${d}/repeats/${assembly_id()}_genomic.repeats.gff3
    echo '{"pipeline":"repeat_masking","outputs":[{"type":"fasta","path":"${d}/genome/${assembly_id()}_genomic.softmasked.fa"}]}' > ${d}/output_manifest.json
    """
}

// ---------------------------------------------------------------------------
// Annotation fan — parallel processes, each emitting its GFF3 path via stdout
// ---------------------------------------------------------------------------

process RNASEQ {
    label 'process_launcher'
    input:
    val repeat_outdir

    output:
    stdout emit: gff3_path

    script:
    def softmasked = "${repeat_outdir}/genome/${assembly_id()}_genomic.softmasked.fa"
    def sheet_arg  = params.sample_sheet
        ? "--sample_sheet             ${params.sample_sheet}"
        : params.rnaseq_bioproject
            ? "--rnaseq_bioproject        ${params.rnaseq_bioproject}"
            : "--rnaseq_run_accessions    ${params.rnaseq_run_accessions}"
    def manifest   = "${stage_outdir('rnaseq')}/output_manifest.json"
    """
    ${nf_run('rnaseq', [
        "--genome_fasta ${softmasked}",
        sheet_arg,
    ])}
    ${extract_gff3_from_manifest(manifest)}
    """

    stub:
    def d = stage_outdir('rnaseq')
    """
    mkdir -p ${d}/rnaseq
    printf '##gff-version 3\\n' > ${d}/rnaseq/merged_rnaseq.gff3
    echo '{"pipeline":"rnaseq","outputs":[{"type":"gff3","path":"${d}/rnaseq/merged_rnaseq.gff3"}]}' > ${d}/output_manifest.json
    python3 -c "print('${d}/rnaseq/merged_rnaseq.gff3', end='')"
    """
}

process LONG_READ {
    label 'process_launcher'
    input:
    val repeat_outdir

    output:
    stdout emit: gff3_path

    script:
    def softmasked = "${repeat_outdir}/genome/${assembly_id()}_genomic.softmasked.fa"
    def sheet_arg  = params.long_read_sample_sheet
        ? "--sample_sheet        ${params.long_read_sample_sheet}"
        : "--long_read_bioproject ${params.long_read_bioproject}"
    def manifest   = "${stage_outdir('long_read')}/output_manifest.json"
    """
    ${nf_run('long_read', [
        "--genome_fasta ${softmasked}",
        sheet_arg,
    ])}
    ${extract_gff3_from_manifest(manifest)}
    """

    stub:
    def d = stage_outdir('long_read')
    """
    mkdir -p ${d}/long_read
    printf '##gff-version 3\\n' > ${d}/long_read/merged_long_read.gff3
    echo '{"pipeline":"long_read","outputs":[{"type":"gff3","path":"${d}/long_read/merged_long_read.gff3"}]}' > ${d}/output_manifest.json
    python3 -c "print('${d}/long_read/merged_long_read.gff3', end='')"
    """
}

process BEST_TARGETED {
    label 'process_launcher'
    input:
    val repeat_outdir

    output:
    stdout emit: gff3_path

    script:
    def softmasked = "${repeat_outdir}/genome/${assembly_id()}_genomic.softmasked.fa"
    def cdna_arg   = params.cdna_fasta    ? "--cdna_fasta    ${params.cdna_fasta}"    : ''
    def prot_arg   = params.protein_fasta ? "--protein_fasta ${params.protein_fasta}" : ''
    def manifest   = "${stage_outdir('best_targeted')}/output_manifest.json"
    """
    ${nf_run('best_targeted', [
        "--genome_fasta ${softmasked}",
        cdna_arg,
        prot_arg,
    ])}
    ${extract_gff3_from_manifest(manifest)}
    """

    stub:
    def d = stage_outdir('best_targeted')
    """
    mkdir -p ${d}/best_targeted
    printf '##gff-version 3\\n' > ${d}/best_targeted/best_targeted.gff3
    echo '{"pipeline":"best_targeted","outputs":[{"type":"gff3","path":"${d}/best_targeted/best_targeted.gff3"}]}' > ${d}/output_manifest.json
    python3 -c "print('${d}/best_targeted/best_targeted.gff3', end='')"
    """
}

process PROJECTION {
    label 'process_launcher'
    input:
    val repeat_outdir

    output:
    stdout emit: gff3_path

    script:
    def softmasked = "${repeat_outdir}/genome/${assembly_id()}_genomic.softmasked.fa"
    def manifest   = "${stage_outdir('projection')}/output_manifest.json"
    """
    ${nf_run('projection', [
        "--query_fasta  ${softmasked}",
        "--source_fasta ${params.source_fasta}",
        "--source_gff3  ${params.source_gff3}",
    ])}
    ${extract_gff3_from_manifest(manifest)}
    """

    stub:
    def d = stage_outdir('projection')
    """
    mkdir -p ${d}/projection
    printf '##gff-version 3\\n' > ${d}/projection/projected.gff3
    echo '{"pipeline":"projection","outputs":[{"type":"gff3","path":"${d}/projection/projected.gff3"}]}' > ${d}/output_manifest.json
    python3 -c "print('${d}/projection/projected.gff3', end='')"
    """
}

process AB_INITIO {
    label 'process_launcher'
    input:
    val repeat_outdir

    output:
    stdout emit: gff3_path

    script:
    def softmasked = "${repeat_outdir}/genome/${assembly_id()}_genomic.softmasked.fa"
    def manifest   = "${stage_outdir('ab_initio')}/output_manifest.json"
    """
    ${nf_run('ab_initio', [
        "--genome_fasta ${softmasked}",
        "--species      ${params.ab_initio_species}",
    ])}
    ${extract_gff3_from_manifest(manifest)}
    """

    stub:
    def d = stage_outdir('ab_initio')
    """
    mkdir -p ${d}/ab_initio
    printf '##gff-version 3\\n' > ${d}/ab_initio/merged_ab_initio.gff3
    echo '{"pipeline":"ab_initio","outputs":[{"type":"gff3","path":"${d}/ab_initio/merged_ab_initio.gff3"}]}' > ${d}/output_manifest.json
    python3 -c "print('${d}/ab_initio/merged_ab_initio.gff3', end='')"
    """
}

process IGTR {
    label 'process_launcher'
    input:
    val repeat_outdir

    output:
    stdout emit: gff3_path

    script:
    def softmasked = "${repeat_outdir}/genome/${assembly_id()}_genomic.softmasked.fa"
    def manifest   = "${stage_outdir('igtr')}/output_manifest.json"
    """
    ${nf_run('igtr', [
        "--genome_fasta  ${softmasked}",
        "--igtr_proteins ${params.igtr_proteins}",
    ])}
    ${extract_gff3_from_manifest(manifest)}
    """

    stub:
    def d = stage_outdir('igtr')
    """
    mkdir -p ${d}/igtr
    printf '##gff-version 3\\n' > ${d}/igtr/igtr.gff3
    echo '{"pipeline":"igtr","outputs":[{"type":"gff3","path":"${d}/igtr/igtr.gff3"}]}' > ${d}/output_manifest.json
    python3 -c "print('${d}/igtr/igtr.gff3', end='')"
    """
}

process GENBLAST_HOMOLOGY {
    label 'process_launcher'
    input:
    val repeat_outdir

    output:
    stdout emit: gff3_path

    script:
    def softmasked = "${repeat_outdir}/genome/${assembly_id()}_genomic.softmasked.fa"
    def prot_arg   = params.uniprot_fasta
        ? "--uniprot_fasta    ${params.uniprot_fasta}"
        : "--uniprot_taxon_id ${params.uniprot_taxon_id}"
    def manifest   = "${stage_outdir('genblast_homology')}/output_manifest.json"
    """
    ${nf_run('genblast_homology', [
        "--genome_fasta  ${softmasked}",
        prot_arg,
    ])}
    ${extract_gff3_from_manifest(manifest)}
    """

    stub:
    def d = stage_outdir('genblast_homology')
    """
    mkdir -p ${d}/genblast_homology
    printf '##gff-version 3\\n' > ${d}/genblast_homology/genblast.merged.gff3
    echo '{"pipeline":"genblast_homology","outputs":[{"type":"gff3","path":"${d}/genblast_homology/genblast.merged.gff3"}]}' > ${d}/output_manifest.json
    python3 -c "print('${d}/genblast_homology/genblast.merged.gff3', end='')"
    """
}

process SHORT_NCRNA {
    label 'process_launcher'
    input:
    val repeat_outdir

    output:
    stdout emit: gff3_path

    script:
    def softmasked = "${repeat_outdir}/genome/${assembly_id()}_genomic.softmasked.fa"
    def mirna_arg  = params.mirna_fasta ? "--mirna_fasta ${params.mirna_fasta}" : ''
    def manifest   = "${stage_outdir('short_ncrna')}/output_manifest.json"
    """
    ${nf_run('short_ncrna', [
        "--genome_fasta ${softmasked}",
        "--rfam_cm      ${params.rfam_cm}",
        mirna_arg,
    ])}
    ${extract_gff3_from_manifest(manifest)}
    """

    stub:
    def d = stage_outdir('short_ncrna')
    """
    mkdir -p ${d}/short_ncrna
    printf '##gff-version 3\\n' > ${d}/short_ncrna/ncrna.gff3
    echo '{"pipeline":"short_ncrna","outputs":[{"type":"gff3","path":"${d}/short_ncrna/ncrna.gff3"}]}' > ${d}/output_manifest.json
    python3 -c "print('${d}/short_ncrna/ncrna.gff3', end='')"
    """
}

process REFSEQ_IMPORT {
    label 'process_launcher'
    input:
    val load_outdir   // needs synonyms TSV from load_assembly

    output:
    stdout emit: gff3_path

    script:
    def synonyms = "${load_outdir}/genome/${assembly_id()}.synonyms.tsv"
    def manifest = "${stage_outdir('refseq_import')}/output_manifest.json"
    """
    ${nf_run('refseq_import', [
        "--assembly_refseq_accession ${params.assembly_refseq_accession}",
        "--assembly_name             ${params.assembly_name}",
        "--synonyms_tsv              ${synonyms}",
    ])}
    ${extract_gff3_from_manifest(manifest)}
    """

    stub:
    def d = stage_outdir('refseq_import')
    """
    mkdir -p ${d}/refseq_import
    printf '##gff-version 3\\n' > ${d}/refseq_import/refseq_parsed.gff3
    echo '{"pipeline":"refseq_import","outputs":[{"type":"gff3","path":"${d}/refseq_import/refseq_parsed.gff3"}]}' > ${d}/output_manifest.json
    python3 -c "print('${d}/refseq_import/refseq_parsed.gff3', end='')"
    """
}

// ---------------------------------------------------------------------------
// Consolidate — waits for all annotation layers
// ---------------------------------------------------------------------------

process CONSOLIDATE {
    label 'process_launcher'
    input:
    val gff3_paths     // List<String> — all annotation GFF3 absolute paths
    val layer_prios    // JSON string of priority map

    output:
    stdout emit: gff3_path

    script:
    // Comma-separated list consumed by consolidate/main.nf's split logic
    def gff3_list = gff3_paths instanceof List
        ? gff3_paths.collect { it.trim() }.join(',')
        : gff3_paths.toString().trim()
    def manifest  = "${stage_outdir('consolidate')}/output_manifest.json"
    """
    ${nf_run('consolidate', [
        "--gff3_files       '${gff3_list}'",
        "--layer_priorities '${layer_prios}'",
    ])}
    ${extract_gff3_from_manifest(manifest)}
    """

    stub:
    def d = stage_outdir('consolidate')
    """
    mkdir -p ${d}/consolidate
    printf '##gff-version 3\\n' > ${d}/consolidate/consolidated.gff3
    echo '{"pipeline":"consolidate","outputs":[{"type":"gff3","path":"${d}/consolidate/consolidated.gff3"}]}' > ${d}/output_manifest.json
    python3 -c "print('${d}/consolidate/consolidated.gff3', end='')"
    """
}

// ---------------------------------------------------------------------------
// Sequential finishing stages
// ---------------------------------------------------------------------------

process UTR_ADDITION {
    label 'process_launcher'
    input:
    val consolidated_gff3
    val repeat_outdir

    output:
    stdout emit: gff3_path

    script:
    // Donor files: list existing stage output paths in priority order.
    // UTR addition pipeline splits on comma, so this is correct.
    // Only include stages that actually ran (paths may not exist if stage was skipped).
    def donor_candidates = [
        params.long_read_sample_sheet || params.long_read_bioproject
            ? "${stage_outdir('long_read')}/long_read/merged_long_read.gff3" : null,
        params.cdna_fasta || params.protein_fasta
            ? "${stage_outdir('best_targeted')}/best_targeted/best_targeted.gff3" : null,
        params.sample_sheet || params.rnaseq_bioproject || params.rnaseq_run_accessions
            ? "${stage_outdir('rnaseq')}/rnaseq/merged_rnaseq.gff3" : null,
    ].findAll { it != null }
    def donor_str  = donor_candidates.join(',')
    def manifest   = "${stage_outdir('utr_addition')}/output_manifest.json"
    """
    ${nf_run('utr_addition', [
        "--consolidated_gff3 ${consolidated_gff3}",
        "--donor_gff3_files  '${donor_str}'",
    ])}
    ${extract_gff3_from_manifest(manifest)}
    """

    stub:
    def d = stage_outdir('utr_addition')
    """
    mkdir -p ${d}/utr_addition
    printf '##gff-version 3\\n' > ${d}/utr_addition/genes.with_utrs.gff3
    echo '{"pipeline":"utr_addition","outputs":[{"type":"gff3","path":"${d}/utr_addition/genes.with_utrs.gff3"}]}' > ${d}/output_manifest.json
    python3 -c "print('${d}/utr_addition/genes.with_utrs.gff3', end='')"
    """
}

process FINALISE_GENESET {
    label 'process_launcher'
    input:
    val utr_gff3
    val repeat_outdir

    output:
    stdout emit: gff3_path

    script:
    def repeat_gff3 = "${repeat_outdir}/repeats/${assembly_id()}_genomic.repeats.gff3"
    def seleno_arg  = params.selenoprotein_fasta
        ? "--selenoprotein_fasta ${params.selenoprotein_fasta}"
        : ''
    def manifest    = "${stage_outdir('finalise_geneset')}/output_manifest.json"
    """
    ${nf_run('finalise_geneset', [
        "--input_gff3  ${utr_gff3}",
        "--repeat_gff3 ${repeat_gff3}",
        seleno_arg,
    ])}
    ${extract_gff3_from_manifest(manifest)}
    """

    stub:
    def d = stage_outdir('finalise_geneset')
    """
    mkdir -p ${d}/finalise_geneset
    printf '##gff-version 3\\n' > ${d}/finalise_geneset/${assembly_id()}.canonical.gff3
    echo '{"pipeline":"finalise_geneset","outputs":[{"type":"gff3","path":"${d}/finalise_geneset/${assembly_id()}.canonical.gff3"}]}' > ${d}/output_manifest.json
    python3 -c "print('${d}/finalise_geneset/${assembly_id()}.canonical.gff3', end='')"
    """
}

process GFF3_TO_CORE {
    label 'process_launcher'
    input:
    val final_gff3
    val load_outdir

    output:
    val "${stage_outdir('gff3_to_core')}/output_manifest.json", emit: manifest

    script:
    def genome_fai = "${load_outdir}/genome/${assembly_id()}_genomic.fna.fai"
    def synonyms   = "${load_outdir}/genome/${assembly_id()}.synonyms.tsv"
    def sp_name    = params.species_name ?: params.assembly_name.toLowerCase()
    def asm_ver    = params.assembly_version ?: params.assembly_name
    """
    ${nf_run('gff3_to_core', [
        "--input_gff3         ${final_gff3}",
        "--db_host            ${params.db_host}",
        "--db_port            ${params.db_port}",
        "--db_user            ${params.db_user}",
        "--db_password        '${params.db_password}'",
        "--db_name            ${params.db_name}",
        "--assembly           ${asm_ver}",
        "--species_name       '${sp_name}'",
        "--stable_id_prefix   '${params.stable_id_prefix}'",
        "--stable_id_start    ${params.stable_id_start}",
        "--genome_fai         ${genome_fai}",
        "--synonyms_tsv       ${synonyms}",
    ])}
    """

    stub:
    def d = stage_outdir('gff3_to_core')
    """
    mkdir -p ${d}
    echo '{"pipeline":"gff3_to_core","outputs":[]}' > ${d}/output_manifest.json
    """
}

// ---------------------------------------------------------------------------
// Main workflow
// ---------------------------------------------------------------------------

workflow {

    validate_params()

    // Determine which stages are enabled (used for layer priority auto-generation)
    def stages = [
        rnaseq:     params.sample_sheet || params.rnaseq_bioproject || params.rnaseq_run_accessions,
        long_read:  params.long_read_sample_sheet || params.long_read_bioproject,
        targeted:   params.cdna_fasta || params.protein_fasta,
        genblast:   params.uniprot_fasta || params.uniprot_taxon_id,
        projection: params.source_fasta && params.source_gff3,
        ncrna:      params.rfam_cm,
        igtr:       params.igtr_proteins,
        refseq:     params.assembly_refseq_accession,
    ]
    def layer_prios = params.layer_priorities ?: default_layer_priorities(stages)

    // ── Stage 1: Download and index genome ────────────────────────────────
    LOAD_ASSEMBLY()

    // ── Stage 2: Repeat masking ───────────────────────────────────────────
    REPEAT_MASKING(LOAD_ASSEMBLY.out.outdir)
    ch_repeat_outdir = REPEAT_MASKING.out.outdir

    // ── Stage 3: Annotation fan (parallel) ───────────────────────────────
    ch_evidence_gff3 = Channel.empty()

    if (stages.rnaseq) {
        RNASEQ(ch_repeat_outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(RNASEQ.out.gff3_path.map { it.trim() })
    }

    if (stages.long_read) {
        LONG_READ(ch_repeat_outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(LONG_READ.out.gff3_path.map { it.trim() })
    }

    if (stages.targeted) {
        BEST_TARGETED(ch_repeat_outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(BEST_TARGETED.out.gff3_path.map { it.trim() })
    }

    if (stages.projection) {
        PROJECTION(ch_repeat_outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(PROJECTION.out.gff3_path.map { it.trim() })
    }

    // Ab initio always runs
    AB_INITIO(ch_repeat_outdir)
    ch_evidence_gff3 = ch_evidence_gff3.mix(AB_INITIO.out.gff3_path.map { it.trim() })

    if (stages.igtr) {
        IGTR(ch_repeat_outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(IGTR.out.gff3_path.map { it.trim() })
    }

    if (stages.genblast) {
        GENBLAST_HOMOLOGY(ch_repeat_outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(GENBLAST_HOMOLOGY.out.gff3_path.map { it.trim() })
    }

    if (stages.ncrna) {
        SHORT_NCRNA(ch_repeat_outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(SHORT_NCRNA.out.gff3_path.map { it.trim() })
    }

    if (stages.refseq) {
        REFSEQ_IMPORT(LOAD_ASSEMBLY.out.outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(REFSEQ_IMPORT.out.gff3_path.map { it.trim() })
    }

    // ── Stage 4: Consolidate ──────────────────────────────────────────────
    CONSOLIDATE(
        ch_evidence_gff3.collect(),
        layer_prios
    )

    // ── Stage 5: UTR addition ─────────────────────────────────────────────
    UTR_ADDITION(
        CONSOLIDATE.out.gff3_path.map { it.trim() },
        ch_repeat_outdir
    )

    // ── Stage 6: Finalise geneset ─────────────────────────────────────────
    FINALISE_GENESET(
        UTR_ADDITION.out.gff3_path.map { it.trim() },
        ch_repeat_outdir
    )

    // ── Stage 7: Load into core DB (optional) ─────────────────────────────
    if (params.db_host && params.db_name && params.db_user) {
        GFF3_TO_CORE(
            FINALISE_GENESET.out.gff3_path.map { it.trim() },
            LOAD_ASSEMBLY.out.outdir
        )
    } else {
        log.info "Skipping GFF3_TO_CORE — set --db_host/--db_name/--db_user to enable"
    }

    // Final GFF3 location for user reference
    FINALISE_GENESET.out.gff3_path.map { it.trim() }.subscribe { path ->
        log.info "\n  ✓ Final annotation: ${path}\n"
    }
}
