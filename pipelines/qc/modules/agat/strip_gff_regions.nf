process STRIP_GFF_REGIONS {

    tag { meta.id }

    input:
        tuple val(meta), path(gff3)

    output:
        tuple val(meta), path("${meta.sample ?: meta.id ?: gff3.simpleName}.no_region.gff3")

    script:
        def stem = meta.sample ?: meta.id ?: gff3.simpleName

        """
        grep -v \$'\\tregion\\t' ${gff3} > ${stem}.no_region.gff3
        """
}
