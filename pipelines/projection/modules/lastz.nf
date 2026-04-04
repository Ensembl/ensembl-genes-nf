// LASTZ
// Whole-genome pairwise alignment: target (reference/source) vs query (new genome).
// Produces a UCSC chain file via lastz → axtChain → chainSort → chainMergeSort.
//
// Expects:
//   target_fasta  — softmasked or hard-masked source genome (split by chrom for parallelism)
//   query_fasta   — target (new) genome FASTA
//
// Output: <prefix>.chain  (one chain per target sequence)

process LASTZ {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::lastz=1.04.22 bioconda::ucsc-axtchain bioconda::ucsc-chainantirepeat bioconda::ucsc-chainsort"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d9785a536a3c6f5c6f55ca96ccd51f3d83b2ba77:latest' :
        'quay.io/biocontainers/mulled-v2-d9785a536a3c6f5c6f55ca96ccd51f3d83b2ba77:latest' }"

    input:
    tuple val(meta), path(target_fasta)
    path  query_fasta

    output:
    tuple val(meta), path("*.chain"), emit: chain
    path  "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: ''
    // lastz options: softmask lower-case bases, output as axt for chaining
    def lastz_opts = 'T=2 C=2 H=2000 Y=3400 L=6000 K=2200 --format=axt'
    """
    lastz \\
        ${target_fasta}[multiple] \\
        ${query_fasta}[multiple] \\
        ${lastz_opts} \\
        ${args} \\
        > ${prefix}.axt

    axtChain \\
        -linearGap=medium \\
        ${prefix}.axt \\
        ${target_fasta} \\
        ${query_fasta} \\
        stdout \\
    | chainAntiRepeat \\
        ${target_fasta} \\
        ${query_fasta} \\
        stdin \\
        stdout \\
    | chainSort stdin ${prefix}.chain

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lastz: \$(lastz --version 2>&1 | grep -oP 'version \\K[\\d.]+' || echo 'unknown')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    # Minimal valid chain header
    printf 'chain 1000 chr1 10000 + 0 10000 chr1 10000 + 0 10000 1\\n' > ${prefix}.chain
    printf '10000\\n\\n' >> ${prefix}.chain

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lastz: 1.04.22
    END_VERSIONS
    """
}
