#!/usr/bin/env nextflow
/*
========================================================================================
    ensembl-genes-nf  —  Full Annotation Master Pipeline
========================================================================================
    Orchestrates all annotation sub-pipelines end-to-end by running each as a child
    `nextflow run` process.  No eHive required.

    Each stage is a Nextflow process whose script block calls one sub-pipeline.
    Output paths flow between stages as strings derived from the known publishDir
    conventions of each sub-pipeline.  The 9 annotation layers run in parallel;
    their GFF3 outputs are collected before consolidation.

    Usage:
      nextflow run . \
        --assembly_accession GCA_964261345.1 \
        --assembly_name      mHetGlaV3 \
        --outdir             /hps/scratch/.../hetgla_experiment \
        --nf_work_root       /hps/scratch/.../hetgla_experiment/nf_work \
        -profile             cluster

    Required:
      --assembly_accession   GCA accession
      --assembly_name        Assembly name (no spaces)
      --outdir               Root output directory
      --nf_work_root         Nextflow work root for sub-pipelines

    Evidence (optional — stages are skipped if not supplied):
      --sample_sheet         RNA-seq CSV (id,fastq_1,fastq_2,strandedness)
      --long_read_sample_sheet  Long-read TSV
      --cdna_fasta           cDNA FASTA for best_targeted
      --protein_fasta        Species proteins for best_targeted
      --uniprot_fasta        Clade UniProt FASTA for genblast_homology
      --igtr_proteins        IG/TR FASTA for igtr
      --rfam_cm              Rfam.cm for short_ncrna
      --mirna_fasta          miRNA FASTA for short_ncrna
      --source_fasta         Source genome for projection
      --source_gff3          Source annotation for projection
      --assembly_refseq_accession  GCF accession for refseq_import
      --repbase_library      RepeatMasker library name (e.g. rodentia)
      --custom_repeat_library  Path to custom repeat FASTA

    Core DB loading (optional — gff3_to_core stage skipped if not supplied):
      --db_host, --db_port, --db_user, --db_password, --db_name
      --stable_id_prefix     e.g. '' for human, 'ENSTGU' for zebra finch
      --stable_id_start      Integer start for new IDs (from registry)
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

params.sample_sheet              = null
params.long_read_sample_sheet    = null
params.cdna_fasta                = null
params.protein_fasta             = null
params.uniprot_fasta             = null
params.igtr_proteins             = null
params.rfam_cm                   = null
params.mirna_fasta               = null
params.source_fasta              = null
params.source_gff3               = null
params.assembly_refseq_accession = null
params.repbase_library           = 'vertebrates'
params.custom_repeat_library     = null

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
params.assembly_version          = params.assembly_name

params.nextflow_profile          = 'cluster'  // profile passed to sub-pipelines

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
    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║       ensembl-genes-nf  Full Annotation Pipeline     ║
    ╠══════════════════════════════════════════════════════╣
    ║  Assembly  : ${params.assembly_accession} (${params.assembly_name})
    ║  Outdir    : ${params.outdir}
    ║  Work root : ${params.nf_work_root}
    ║
    ║  Evidence layers:
    ║    RNA-seq       : ${params.sample_sheet          ?: 'SKIP'}
    ║    Long-read     : ${params.long_read_sample_sheet ?: 'SKIP'}
    ║    cDNA targeted : ${params.cdna_fasta            ?: 'SKIP'}
    ║    Homology      : ${params.uniprot_fasta         ?: 'SKIP'}
    ║    IGTR          : ${params.igtr_proteins         ?: 'SKIP'}
    ║    ncRNA         : ${params.rfam_cm               ?: 'SKIP'}
    ║    Projection    : ${params.source_gff3           ?: 'SKIP'}
    ║    RefSeq        : ${params.assembly_refseq_accession ?: 'SKIP'}
    ║
    ║  Core loading    : ${params.db_name ?: 'SKIP'}
    ╚══════════════════════════════════════════════════════╝
    """.stripIndent()
}

// ---------------------------------------------------------------------------
// Helpers — predictable publishDir paths derived from known conventions
// ---------------------------------------------------------------------------
def assembly_id() {
    "${params.assembly_accession}_${params.assembly_name}"
}
def stage_outdir(String stage) {
    "${params.outdir}/${stage}"
}
def nf_work(String stage) {
    "${params.nf_work_root}/${stage}"
}
def pipelines_dir() {
    "${projectDir}/pipelines"
}
def nf_run(String stage, List extra_args = []) {
    // Base nextflow run command for any sub-pipeline
    [
        "nextflow run ${pipelines_dir()}/${stage}/main.nf",
        "-profile ${params.nextflow_profile}",
        "-work-dir ${nf_work(stage)}",
        "-resume",
        "--outdir ${stage_outdir(stage)}",
    ].plus(extra_args).join(" \\\n    ")
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
        ? "--custom_library ${params.custom_repeat_library}"
        : "--repbase_library ${params.repbase_library}"
    """
    ${nf_run('repeat_masking', [
        "--genome_fasta ${genome_fasta}",
        lib_arg,
    ])}
    """
}

// ---------------------------------------------------------------------------
// Annotation fan — 9 processes run in parallel
// Each emits the path to its output GFF3 (read from output_manifest.json)
// Skipped stages emit nothing; collect() gathers whatever ran.
// ---------------------------------------------------------------------------

process RNASEQ {
    label 'process_launcher'

    input:
    val repeat_outdir

    output:
    stdout emit: gff3_path

    script:
    def softmasked = "${repeat_outdir}/genome/${assembly_id()}_genomic.softmasked.fa"
    """
    ${nf_run('rnaseq', [
        "--genome_fasta  ${softmasked}",
        "--sample_sheet  ${params.sample_sheet}",
    ])}
    python3 -c "
import json
m = json.load(open('${stage_outdir('rnaseq')}/output_manifest.json'))
print([o['path'] for o in m['outputs'] if o['type']=='gff3'][0], end='')
"
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
    def cdna_arg   = params.cdna_fasta   ? "--cdna_fasta   ${params.cdna_fasta}"   : ''
    def prot_arg   = params.protein_fasta ? "--protein_fasta ${params.protein_fasta}" : ''
    """
    ${nf_run('best_targeted', [
        "--genome_fasta ${softmasked}",
        cdna_arg,
        prot_arg,
    ])}
    python3 -c "
import json
m = json.load(open('${stage_outdir('best_targeted')}/output_manifest.json'))
print([o['path'] for o in m['outputs'] if o['type']=='gff3'][0], end='')
"
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
    """
    ${nf_run('projection', [
        "--query_fasta  ${softmasked}",
        "--source_fasta ${params.source_fasta}",
        "--source_gff3  ${params.source_gff3}",
    ])}
    python3 -c "
import json
m = json.load(open('${stage_outdir('projection')}/output_manifest.json'))
print([o['path'] for o in m['outputs'] if o['type']=='gff3'][0], end='')
"
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
    """
    ${nf_run('ab_initio', [
        "--genome_fasta ${softmasked}",
        "--species      ${params.assembly_name}",
    ])}
    python3 -c "
import json
m = json.load(open('${stage_outdir('ab_initio')}/output_manifest.json'))
print([o['path'] for o in m['outputs'] if o['type']=='gff3'][0], end='')
"
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
    """
    ${nf_run('igtr', [
        "--genome_fasta  ${softmasked}",
        "--igtr_proteins ${params.igtr_proteins}",
    ])}
    python3 -c "
import json
m = json.load(open('${stage_outdir('igtr')}/output_manifest.json'))
print([o['path'] for o in m['outputs'] if o['type']=='gff3'][0], end='')
"
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
    """
    ${nf_run('long_read', [
        "--genome_fasta ${softmasked}",
        "--sample_sheet ${params.long_read_sample_sheet}",
    ])}
    python3 -c "
import json
m = json.load(open('${stage_outdir('long_read')}/output_manifest.json'))
print([o['path'] for o in m['outputs'] if o['type']=='gff3'][0], end='')
"
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
    """
    ${nf_run('genblast_homology', [
        "--genome_fasta  ${softmasked}",
        "--uniprot_fasta ${params.uniprot_fasta}",
    ])}
    python3 -c "
import json
m = json.load(open('${stage_outdir('genblast_homology')}/output_manifest.json'))
print([o['path'] for o in m['outputs'] if o['type']=='gff3'][0], end='')
"
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
    """
    ${nf_run('short_ncrna', [
        "--genome_fasta ${softmasked}",
        "--rfam_cm      ${params.rfam_cm}",
        mirna_arg,
    ])}
    python3 -c "
import json
m = json.load(open('${stage_outdir('short_ncrna')}/output_manifest.json'))
print([o['path'] for o in m['outputs'] if o['type']=='gff3'][0], end='')
"
    """
}

process REFSEQ_IMPORT {
    label 'process_launcher'

    input:
    val load_outdir   // needs synonyms TSV from load_assembly, not softmasked genome

    output:
    stdout emit: gff3_path

    script:
    def synonyms = "${load_outdir}/genome/${assembly_id()}.synonyms.tsv"
    """
    ${nf_run('refseq_import', [
        "--assembly_refseq_accession ${params.assembly_refseq_accession}",
        "--assembly_name             ${params.assembly_name}",
        "--synonyms_tsv              ${synonyms}",
    ])}
    python3 -c "
import json
m = json.load(open('${stage_outdir('refseq_import')}/output_manifest.json'))
print([o['path'] for o in m['outputs'] if o['type']=='gff3'][0], end='')
"
    """
}

// ---------------------------------------------------------------------------
// Consolidate — waits for all annotation layers
// ---------------------------------------------------------------------------

process CONSOLIDATE {
    label 'process_launcher'

    input:
    val gff3_paths   // list of all annotation GFF3 paths

    output:
    stdout emit: gff3_path

    script:
    // Pass as comma-separated list; consolidate main.nf splits on comma
    def gff3_list = gff3_paths instanceof List
        ? gff3_paths.join(',')
        : gff3_paths
    """
    ${nf_run('consolidate', [
        "--gff3_files '${gff3_list}'",
    ])}
    python3 -c "
import json
m = json.load(open('${stage_outdir('consolidate')}/output_manifest.json'))
print([o['path'] for o in m['outputs'] if o['type']=='gff3'][0], end='')
"
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
    // Donor evidence: all long-read, best_targeted, rnaseq GFF3s
    def donor_dirs = [
        stage_outdir('long_read'),
        stage_outdir('best_targeted'),
        stage_outdir('rnaseq'),
    ]
    def donor_globs = donor_dirs.collect { "${it}/**/*.gff3" }.join(',')
    """
    ${nf_run('utr_addition', [
        "--consolidated_gff3 ${consolidated_gff3}",
        "--donor_gff3_files  '${donor_globs}'",
    ])}
    python3 -c "
import json
m = json.load(open('${stage_outdir('utr_addition')}/output_manifest.json'))
print([o['path'] for o in m['outputs'] if o['type']=='gff3'][0], end='')
"
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
    """
    ${nf_run('finalise_geneset', [
        "--input_gff3  ${utr_gff3}",
        "--repeat_gff3 ${repeat_gff3}",
        seleno_arg,
    ])}
    python3 -c "
import json
m = json.load(open('${stage_outdir('finalise_geneset')}/output_manifest.json'))
print([o['path'] for o in m['outputs'] if o['type']=='gff3'][0], end='')
"
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
    """
    ${nf_run('gff3_to_core', [
        "--input_gff3         ${final_gff3}",
        "--db_host            ${params.db_host}",
        "--db_port            ${params.db_port}",
        "--db_user            ${params.db_user}",
        "--db_password        '${params.db_password}'",
        "--db_name            ${params.db_name}",
        "--assembly           ${params.assembly_version ?: params.assembly_name}",
        "--species_name       ${params.species_name ?: 'unknown'}",
        "--stable_id_prefix   '${params.stable_id_prefix}'",
        "--stable_id_start    ${params.stable_id_start}",
        "--genome_fai         ${genome_fai}",
        "--synonyms_tsv       ${synonyms}",
    ])}
    """
}

// ---------------------------------------------------------------------------
// Main workflow
// ---------------------------------------------------------------------------

workflow {

    validate_params()

    // ── Stage 1: Download and index genome ────────────────────────────────
    LOAD_ASSEMBLY()

    // ── Stage 2: Repeat masking ───────────────────────────────────────────
    REPEAT_MASKING(LOAD_ASSEMBLY.out.outdir)

    ch_repeat_outdir = REPEAT_MASKING.out.outdir

    // ── Stage 3: Annotation fan (parallel) ───────────────────────────────
    // Each conditional branch only runs if the required evidence is supplied.
    // Empty channels are silently skipped in the downstream collect().

    ch_evidence_gff3 = Channel.empty()

    if (params.sample_sheet) {
        RNASEQ(ch_repeat_outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(RNASEQ.out.gff3_path.map { it.trim() })
    }

    if (params.cdna_fasta || params.protein_fasta) {
        BEST_TARGETED(ch_repeat_outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(BEST_TARGETED.out.gff3_path.map { it.trim() })
    }

    if (params.source_fasta && params.source_gff3) {
        PROJECTION(ch_repeat_outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(PROJECTION.out.gff3_path.map { it.trim() })
    }

    // Ab initio always runs (uses RepeatMasker-trained species model)
    AB_INITIO(ch_repeat_outdir)
    ch_evidence_gff3 = ch_evidence_gff3.mix(AB_INITIO.out.gff3_path.map { it.trim() })

    if (params.igtr_proteins) {
        IGTR(ch_repeat_outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(IGTR.out.gff3_path.map { it.trim() })
    }

    if (params.long_read_sample_sheet) {
        LONG_READ(ch_repeat_outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(LONG_READ.out.gff3_path.map { it.trim() })
    }

    if (params.uniprot_fasta) {
        GENBLAST_HOMOLOGY(ch_repeat_outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(GENBLAST_HOMOLOGY.out.gff3_path.map { it.trim() })
    }

    if (params.rfam_cm) {
        SHORT_NCRNA(ch_repeat_outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(SHORT_NCRNA.out.gff3_path.map { it.trim() })
    }

    if (params.assembly_refseq_accession) {
        REFSEQ_IMPORT(LOAD_ASSEMBLY.out.outdir)
        ch_evidence_gff3 = ch_evidence_gff3.mix(REFSEQ_IMPORT.out.gff3_path.map { it.trim() })
    }

    // ── Stage 4: Consolidate ──────────────────────────────────────────────
    // .collect() blocks until every evidence channel has emitted, then passes
    // the full list to CONSOLIDATE in a single call.
    CONSOLIDATE(ch_evidence_gff3.collect())

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

    // ── Stage 7: Load into core DB (skipped if no DB params) ─────────────
    if (params.db_host && params.db_name && params.db_user) {
        GFF3_TO_CORE(
            FINALISE_GENESET.out.gff3_path.map { it.trim() },
            LOAD_ASSEMBLY.out.outdir
        )
    } else {
        log.info "Skipping GFF3_TO_CORE (--db_host/--db_name/--db_user not set)"
    }
}
