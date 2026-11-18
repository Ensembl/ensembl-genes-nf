def shell_quote(String value) {
    "'${value.replace("'", "'\\''")}'"
}

process RUN_FUNANNOTATE {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_high'
    container 'docker://nextgenusfs/funannotate'
    containerOptions '--cleanenv'
    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${meta.assembly_accession}/${meta.safe_id}/${filename}" }
    input:
    tuple val(meta), path(staged)
    val cpus

    output:
    tuple val(meta), path('predict_results'), emit: results
    path 'predict_results/*.gff3', emit: gff3
    path 'funannotate.stdout.log', emit: stdout
    path 'funannotate.stderr.log', emit: stderr
    path 'versions.yml', emit: versions

    script:
    def args = [
        '--input', "staged/${meta.genome_name}",
        '--out', 'predict_results',
        '--species', shell_quote(meta.species),
        '--rna_bam staged/all_runs.merged.sorted.bam',
        meta.files.transcripts.withIndex().collect { pair -> "--transcript_evidence staged/transcript_${String.format('%02d', pair[1] + 1)}.fa" },
        meta.files.scallop ? '--transcript_evidence staged/scallop_transcripts.fa' : null,
        meta.files.stringtie ? '--stringtie staged/stringtie.gtf' : null,
        meta.files.protein_alignments ? '--protein_alignments staged/genblast_protein_alignments.gff3' : null,
        meta.files.trnascan ? '--trnascan staged/trnascan.gff' : null,
        '--cpus', cpus
    ].flatten().findAll { value -> value != null && value.toString() }.collect { value -> value.toString() }.join(' ')
    """
    # The host-side launcher in .venv/bin/funannotate starts this same image.
    # Nextflow is already inside the image, so call its executable directly.
    command -v funannotate
    funannotate predict ${args} > funannotate.stdout.log 2> funannotate.stderr.log
    cat > versions.yml <<VERSIONS
"${task.process}":
    funannotate: \$(funannotate --version 2>&1 | head -n1)
VERSIONS
    """

    stub:
    """
    mkdir -p predict_results
    touch predict_results/${meta.safe_id}.gff3 funannotate.stdout.log funannotate.stderr.log versions.yml
    """
}
