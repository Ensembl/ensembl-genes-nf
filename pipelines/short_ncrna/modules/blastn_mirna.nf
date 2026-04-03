// BLASTN_MIRNA — BLASTN search of miRBase sequences against genome.
// Uses the BLAST suite container (same as blast_blastp.nf in long_read).
// Outputs tabular format (-outfmt 6) for parsing by filter_ncrna.py.

process BLASTN_MIRNA {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::blast=2.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.15.0--pl5321h6f7f691_1' :
        'biocontainers/blast:2.15.0--pl5321h6f7f691_1' }"

    input:
    tuple val(meta),  path(mirna_fasta)    // miRBase sequences (query)
    tuple val(meta2), path(genome_db, stageAs: 'blast_db/*')  // genome BLAST DB

    output:
    tuple val(meta), path("*.blastn.tsv"), emit: tsv
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: meta.id
    def args    = task.ext.args   ?: '-num_threads 3 -perc_identity 80 -evalue 0.01'
    def db_name = genome_db[0].baseName
    """
    blastn \\
        -query   ${mirna_fasta} \\
        -db      blast_db/${db_name} \\
        -out     ${prefix}.blastn.tsv \\
        -outfmt  '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | head -1 | sed 's/blastn: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf 'hsa-mir-21\tchr1\t95.0\t88\t4\t0\t1\t88\t1000\t1088\t1.0e-20\t150.0\n' \\
        > ${prefix}.blastn.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: 2.15.0
    END_VERSIONS
    """
}
