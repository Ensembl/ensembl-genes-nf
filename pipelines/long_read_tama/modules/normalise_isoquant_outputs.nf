process NORMALISE_ISOQUANT_OUTPUTS {
    tag "${meta.id}:${scope}:isoquant-normalise"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), val(scope), path(native_dir), path(input_manifest)

    output:
    tuple val(meta), val(scope), path('isoquant_products'), emit: products
    tuple val(meta), val(scope), path('isoquant_products/transcript_models.gtf'), emit: gtf
    path 'versions.yml', emit: versions

    script:
    """
    normalise_isoquant_outputs.py ${native_dir} isoquant_products '${scope}'
    cp ${input_manifest} isoquant_products/isoquant_input_manifest.tsv
    test -s isoquant_products/transcript_models.gtf || { echo 'IsoQuant normalisation produced no transcript models' >&2; exit 1; }
    printf '"%s":\n    isoquant_output_normalisation: python\n' '${task.process}' > versions.yml
    """

    stub:
    """
    mkdir -p isoquant_products
    printf '# stub\nchrStub\tIsoQuant\texon\t1\t4\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' > isoquant_products/transcript_models.gtf
    printf 'read\ttranscript\nread1\tt1\n' | gzip -c > isoquant_products/transcript_model_reads.tsv.gz
    printf 'read\tassignment\nread1\tunique\n' | gzip -c > isoquant_products/read_info.tsv.gz
    printf 'scope\tstable_name\tsource_name\tsource_type\n${scope}\ttranscript_models.gtf\tstub\tnative\n${scope}\tread_info.tsv.gz\tstub\tnative\n' > isoquant_products/product_manifest.tsv
    printf '"%s":\n    isoquant_output_normalisation: stub\n' '${task.process}' > versions.yml
    """
}
