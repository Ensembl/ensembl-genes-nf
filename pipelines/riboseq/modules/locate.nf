process LOCATE {
    tag "${meta.id}"
    label 'process_light'

    input:
    tuple val(meta), val(run)

    output:
    tuple val(meta), path("${run}.collapsed.fa.gz"), emit: collapsed_reads, optional: true
    tuple val(meta), val(run), path("${run}_needs_processing"), emit: needs_processing, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Accept both the historical *_rpfs name and the current getRPF gated
    # *_rpf_20_40 name. The latter is the normal input for rerunning analysis
    # on an existing extracted subset with --fetch false.
    candidates=(
        "${params.collapsed_read_path}/${run}_rpf_20_40.collapsed.fa"
        "${params.collapsed_read_path}/${run}_rpf_20_40.collapsed.fa.gz"
        "${params.collapsed_read_path}/${run}_rpfs.collapsed.fa"
        "${params.collapsed_read_path}/${run}_rpfs.collapsed.fa.gz"
        "${params.collapsed_read_path}/${run}_1_rpfs.collapsed.fa.gz"
    )

    found=""
    for candidate in "\${candidates[@]}"; do
        if [ -f "\$candidate" ]; then
            found="\$candidate"
            break
        fi
    done

    if [ -n "\$found" ]; then
        case "\$found" in
            *.gz) ln -s "\$found" "${run}.collapsed.fa.gz" ;;
            *) gzip -c "\$found" > "${run}.collapsed.fa.gz" ;;
        esac
        echo "Found collapsed file for $run: \$found"
    else
        echo "Collapsed file not found for $run. Needs processing."
        touch "${run}_needs_processing"
    fi
    """

    stub:
    """
    touch "${run}.collapsed.fa.gz"
    touch "${run}_needs_processing"
    """
}
