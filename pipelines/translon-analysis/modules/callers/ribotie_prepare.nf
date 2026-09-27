process PREPARE_RIBOTIE_DATA {
    tag "${meta.id}"
    label 'process_high_memory'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-ribotie:1.0.0'

    input:
    tuple val(meta), path(bam), path(bai), path(gtf), path(fasta)

    output:
    tuple val(meta), path('raw/ribotie.yml'), path("raw/${meta.id}.h5"), path(bam), path(gtf), path(fasta), emit: prepared
    path 'versions.yml', emit: versions, topic: versions

    script:
    """
    set -euo pipefail
    mkdir -p raw
    printf '%s\\n' \\
        'gtf_path: ${gtf}' \\
        'fa_path: ${fasta}' \\
        'ribo_paths:' \\
        '  ${meta.id}: ${bam}' \\
        'h5_path: raw/${meta.id}.h5' > raw/ribotie.yml
    # Parse reference features and mapped reads before the GPU process. This
    # creates the HDF5 data store consumed by RiboTIE's fine-tuning/inference
    # stage and keeps GPU time focused on model execution.
    ribotie raw/ribotie.yml --data
    test -s raw/${meta.id}.h5 || { echo 'RiboTIE data preparation did not create the HDF5 store' >&2; exit 1; }
    printf '"%s":\n    RiboTIE: source-pinned\n    stage: data-preparation\n' '${task.process}' > versions.yml
    """

    stub:
    """
    mkdir -p raw
    printf 'gtf_path: stub\nfa_path: stub\nh5_path: raw/${meta.id}.h5\n' > raw/ribotie.yml
    touch raw/${meta.id}.h5
    printf '"stub":\n    RiboTIE: stub\n    stage: data-preparation\n' > versions.yml
    """
}
