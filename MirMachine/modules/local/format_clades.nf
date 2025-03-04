

process FORMAT_CLADES {
    input:
    path clades

    output:
    path "formatted_clades.txt", emit: formatted_clades

    script:
    """
    awk 'NR>1 {for (i=1; i<=NF; i++) print \$i}' ${clades} > formatted_clades.txt
    """
}