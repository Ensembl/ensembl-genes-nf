process GEDI_PRICE {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/native_outputs", mode: 'copy', saveAs: { filename -> "price/${meta.id}/${filename}" }

    conda "${moduleDir}/environment.yml"
    container "${ params.container_gedi_price ?: (workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cd/cd008e5721759d5909909254c77ec449778e0fc7c669b7c926b68f0c9059f510/data' :
        'community.wave.seqera.io/library/gedi_price:2392624d5f803049') }"

    input:
    tuple val(meta), path(bams, stageAs: 'bams/*'), path(bais, stageAs: 'bams/*')
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("${prefix}.orfs.tsv")                                                 , emit: orfs_tsv
    tuple val(meta), path("${prefix}.orfs.cit")                                                 , emit: orfs_cit, optional: true
    tuple val(meta), path("${prefix}.orfs.cit.metadata.json")                                   , emit: orfs_metadata, optional: true
    tuple val(meta), path("${prefix}.codons.cit")                                               , emit: codons_cit, optional: true
    tuple val(meta), path("${prefix}.model")                                                    , emit: model, optional: true
    tuple val(meta), path("${prefix}.signal.tsv")                                               , emit: signal, optional: true
    tuple val(meta), path("${prefix}.param")                                                    , emit: param, optional: true
    tuple val("${task.process}"), val('gedi'), eval("gedi -e Version 2>&1 | sed -n 's/.*Gedi version \\([^ ]*\\).*/\\1/p' | head -n 1"), topic: versions, emit: versions_gedi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def price_prefix = prefix
    def index_dir = index.toString()
    def reference_id = meta2.id ?: 'reference'
    def threads = task.cpus ?: 1
    """
    mkdir -p prepared_bams
    reference_fasta=\$(ls ${index_dir}/*.fa ${index_dir}/*.fasta ${index_dir}/*.fna 2>/dev/null | head -n 1)
    test -n "\$reference_fasta" || { echo 'GEDI PRICE index does not contain its reference FASTA' >&2; exit 1; }
    for bam in bams/*.bam; do
        samtools calmd -b -@ ${threads} "\$bam" "\$reference_fasta" > "prepared_bams/\$(basename "\$bam")"
        samtools index -@ ${threads} "prepared_bams/\$(basename "\$bam")"
    done
    ls -1 prepared_bams/*.bam > price_input.bamlist
    # Ribo-seq BAMs from the upstream pipeline may not carry MD tags.  The
    # registry's bamlist2cit wrapper enables variation reconstruction and
    # fails in that case; PRICE does not need variation annotations, so use
    # The prepared BAMs now carry MD/NM tags, so retain mismatch information
    # for PRICE's model-estimation stage.
    gedi -e Bam2CIT -p price_input.bamlist.cit prepared_bams/*.bam

    # .oml member paths are absolute; repoint them at the staged index
    sed "s|file=\\"[^\\"]*/|file=\\"\$PWD/${index_dir}/|g" ${index_dir}/${reference_id}.oml > ${price_prefix}.genomic.oml

    gedi -e Price \\
        -reads price_input.bamlist.cit \\
        -genomic ${price_prefix}.genomic.oml \\
        -prefix ${price_prefix} \\
        -nthreads ${threads} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.orfs.tsv
    touch ${prefix}.orfs.cit
    touch ${prefix}.orfs.cit.metadata.json
    touch ${prefix}.codons.cit
    touch ${prefix}.model
    touch ${prefix}.signal.tsv
    touch ${prefix}.param
    """
}
