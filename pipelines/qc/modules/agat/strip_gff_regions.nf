process STRIP_GFF_REGIONS {

    tag { meta.id }

    input:
        tuple val(meta), path(gff3)

    output:
        tuple val(meta), path("${meta.sample ?: meta.id ?: gff3.simpleName}.no_region.gff3"), emit: cleaned_gff
        path "versions.yml", emit: versions

    script:
        def stem = meta.sample ?: meta.id ?: gff3.simpleName

        """
        grep -v \$'\\tregion\\t' ${gff3} > ${stem}.no_region.gff3

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            grep: \$(grep --version | head -n 1 | sed 's/^grep (GNU grep) //')
        END_VERSIONS
        """

    stub:
        def stem = meta.sample ?: meta.id ?: gff3.simpleName
        """
        touch ${stem}.no_region.gff3

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            grep: 3.8
        END_VERSIONS
        """
}
