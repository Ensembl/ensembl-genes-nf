process FLAIR_JUNCTIONS {
    tag "${meta.id}:flair-junctions"
    label 'process_medium'
    container 'docker://community.wave.seqera.io/library/flair_coreutils_python:9744ce1f6049d5bd'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('flair_junctions.bed'), emit: bed
    path 'versions.yml', emit: versions

    script:
    """
    test -s '${bam}' || { echo 'FLAIR junction extraction received an empty BAM' >&2; exit 1; }
    test -s '${bam}.bai' || { echo 'FLAIR junction extraction received no BAM index' >&2; exit 1; }
    bam_to_junction_bed.py '${bam}' flair_junctions.bed
    test -s flair_junctions.bed || { echo 'No splice junctions were found in the BAM' >&2; exit 1; }
    printf '"%s":\n    flair_junctions: \$(flair --version 2>&1 | tail -n1 || true)\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'chrStub\t1\t3\t.\t1\t+\n' > flair_junctions.bed
    printf '"%s":\n    flair_junctions: stub\n' '${task.process}' > versions.yml
    """
}

process FLAIR_TRANSCRIPTOME {
    tag "${meta.id}:flair-transcriptome"
    label 'process_high_memory'
    container 'docker://community.wave.seqera.io/library/flair_coreutils_python:9744ce1f6049d5bd'

    input:
    tuple val(meta), path(bam), path(bai), path(junctions)
    path reference

    output:
    tuple val(meta), path('flair_products.tar.gz'), emit: products
    tuple val(meta), path('flair_products/isoforms.gtf'), emit: gtf
    path 'versions.yml', emit: versions

    script:
    """
    test -s '${bam}' || { echo 'FLAIR received an empty BAM' >&2; exit 1; }
    test -s '${bam}.bai' || { echo 'FLAIR received no BAM index' >&2; exit 1; }
    flair transcriptome --genomealignedbam '${bam}' --genome '${reference}' \
        --junction_bed '${junctions}' --junction_support 1 --noaligntoannot \
        --threads ${task.cpus} --output flair_result ${params.flair_args ?: ''}
    normalise_flair_outputs.py . flair_products flair_result
    printf 'backend\\tflair\\nmethod\\ttranscriptome\\nscope\\taccession\\n' > flair_products/backend_manifest.tsv
    tar -czf flair_products.tar.gz flair_products
    printf '"%s":\n    flair: \$(flair --version 2>&1 | tail -n1 || true)\n' '${task.process}' > versions.yml
    """

    stub:
    """
    mkdir -p flair_products
    printf 'chrStub\\t0\\t4\\t${meta.id}.flair.1\\t0\\t+\\t0\\t4\\t0\\t1\\t4,\\t0,\\n' > flair_products/isoforms.bed
    printf 'chrStub\\tFLAIR\\texon\\t1\\t4\\t.\\t+\\t.\\tgene_id "${meta.id}.g1"; transcript_id "${meta.id}.flair.1";\\n' > flair_products/isoforms.gtf
    printf '>flair.1\\nACGT\\n' > flair_products/isoforms.fa
    printf 'read_id\\tisoform_id\\n${meta.id}.read1\\t${meta.id}.flair.1\\n' > flair_products/read.map.txt
    printf 'stable_name\\tsource_name\\nisoforms.bed\\tstub\\nisoforms.gtf\\tstub\\nisoforms.fa\\tstub\\nread.map.txt\\tstub\\n' > flair_products/product_manifest.tsv
    printf 'backend\\tflair\\nmethod\\ttranscriptome\\nscope\\taccession\\n' > flair_products/backend_manifest.tsv
    tar -czf flair_products.tar.gz flair_products
    printf '"%s":\\n    flair: stub\\n' '${task.process}' > versions.yml
    """
}

process FLAIR_COMBINE {
    tag "${meta.id}:flair-combine"
    label 'process_high_memory'
    container 'docker://community.wave.seqera.io/library/flair_coreutils_python:9744ce1f6049d5bd'

    input:
    tuple val(meta), path(archives, stageAs: 'input??.tar.gz')

    output:
    tuple val(meta), path('flair_products.tar.gz'), emit: products
    tuple val(meta), path('flair_products/isoforms.gtf'), emit: gtf
    path 'versions.yml', emit: versions

    script:
    """
    : > flair_manifest.tsv
    for archive in input??.tar.gz; do
        product="\${archive%.tar.gz}"
        mkdir -p "\$product"
        tar -xzf "\$archive" --strip-components=1 -C "\$product"
        test -s "\$product/isoforms.bed" || { echo "Missing FLAIR BED in \$product" >&2; exit 1; }
        printf '%s\\tisoform\\t%s/isoforms.bed\\t%s/isoforms.fa\\t%s/read.map.txt\\n' \\
            "\$(basename "\$product")" "\$product" "\$product" "\$product" >> flair_manifest.tsv
    done
    test -s flair_manifest.tsv || { echo 'No FLAIR accession products were provided' >&2; exit 1; }
    flair combine --manifest flair_manifest.tsv --output_prefix flair_cohort \\
        --minpercentusage 0 --filter none --include_se --convert_gtf
    normalise_flair_outputs.py . flair_products flair_cohort
    printf 'backend\\tflair\\nmethod\\tcombine\\nscope\\tcohort\\n' > flair_products/backend_manifest.tsv
    tar -czf flair_products.tar.gz flair_products
    printf '"%s":\\n    flair: \$(flair --version 2>&1 | tail -n1 || true)\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    mkdir -p flair_products
    printf 'chrStub\\t0\\t4\\tflair.cohort.1\\t0\\t+\\t0\\t4\\t0\\t1\\t4,\\t0,\\n' > flair_products/isoforms.bed
    printf 'chrStub\\tFLAIR\\texon\\t1\\t4\\t.\\t+\\t.\\tgene_id "flair.g1"; transcript_id "flair.cohort.1";\\n' > flair_products/isoforms.gtf
    printf '>flair.cohort.1\\nACGT\\n' > flair_products/isoforms.fa
    printf 'read_id\\tisoform_id\\nread1\\tflair.cohort.1\\n' > flair_products/read.map.txt
    printf 'stable_name\\tsource_name\\nisoforms.bed\\tstub\\nisoforms.gtf\\tstub\\nisoforms.fa\\tstub\\nread.map.txt\\tstub\\n' > flair_products/product_manifest.tsv
    printf 'backend\\tflair\\nmethod\\tcombine\\nscope\\tcohort\\n' > flair_products/backend_manifest.tsv
    tar -czf flair_products.tar.gz flair_products
    printf '"%s":\\n    flair: stub\\n' '${task.process}' > versions.yml
    """
}
