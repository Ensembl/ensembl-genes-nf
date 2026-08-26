process RUN_RIBOTIE {
    tag "${meta.id}:${meta.shard_id ?: 'all'}"
    label 'process_high'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container params.ribotie_gpu.toString().toBoolean() ? 'ghcr.io/jackcurragh/translon-ribotie-cuda:1.0.0' : 'ghcr.io/jackcurragh/translon-ribotie:1.0.0'
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    script:
    def args = task.ext.args ?: params.args_ribotie ?: ''
    """
    mkdir -p raw
    python - <<'PY'
    import torch
    if not torch.cuda.is_available():
        raise SystemExit('RiboTIE requires CUDA, but torch.cuda.is_available() is false')
    print('RiboTIE CUDA device:', torch.cuda.get_device_name(0))
    PY
    cat > raw/ribotie.yml <<-END_CONFIG
    gtf_path: ${gtf}
    fa_path: ${fasta}
    ribo_paths:
      ${meta.id}: ${bam}
    h5_path: raw/${meta.id}.h5
    END_CONFIG
    ribotie raw/ribotie.yml ${args}
    test -n "\$(find raw -type f \( -name '*.csv' -o -name '*.gtf' \) | head -1)" || { echo 'RiboTIE produced no native result table' >&2; exit 1; }
    printf '"%s":\n    RiboTIE: source-pinned\n' '${task.process}' > versions.yml
    """
    stub:
    """
    mkdir -p raw
    printf 'id,chrom,start,end,strand,start_codon,score\nRT1,chr1,100,200,+,ATG,0.9\n' > raw/ribotie.csv
    printf '"stub":\n    RiboTIE: stub\n' > versions.yml
    """
}
