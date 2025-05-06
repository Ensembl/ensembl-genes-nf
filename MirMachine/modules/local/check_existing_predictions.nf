

process CHECK_EXISTING_PREDICTIONS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(fasta), val(node), val(model)

    output:
    tuple val(meta), path("*.heatmap.csv"), optional: true, emit: heatmap_csv
    tuple val(meta), path("not_found.txt"), optional: true, emit: not_found

    script:
    """
    predictions_dir="${params.outdir}/mirmachine/${meta.id}/results/predictions/"

    if [ -d "\$predictions_dir" ] && find "\$predictions_dir/heatmap" -name "*.heatmap.csv" | grep -q .; then
        ln -s "\$predictions_dir"/heatmap/*.heatmap.csv ./${meta.id}_${meta.species}.heatmap.csv
    else
        touch not_found.txt
    fi
    """
}
