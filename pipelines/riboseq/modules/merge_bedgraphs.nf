process MERGE_BEDGRAPHS {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    publishDir "${params.outdir}/bedgraphs", mode: 'copy'

    input:
    tuple val(meta), path(bedgraphs)

    output:
    tuple val(meta), path("*.merged.sorted.bedgraph"), emit: bedgraph
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Use glob patterns to collect all bedgraph files by strand
    # Check if we have stranded data
    if ls *.forward.sorted.bedgraph 1> /dev/null 2>&1; then
        # Stranded data - merge forward and reverse separately
        bedtools unionbedg -i *.forward.sorted.bedgraph | \\
            awk 'BEGIN {OFS="\\t"} {sum=0; for(i=4; i<=NF; i++) sum+=\$i; print \$1, \$2, \$3, sum}' | \\
            sort -k1,1 -k2,2n > ${prefix}.merged.forward.sorted.bedgraph

        bedtools unionbedg -i *.reverse.sorted.bedgraph | \\
            awk 'BEGIN {OFS="\\t"} {sum=0; for(i=4; i<=NF; i++) sum+=\$i; print \$1, \$2, \$3, sum}' | \\
            sort -k1,1 -k2,2n > ${prefix}.merged.reverse.sorted.bedgraph
    else
        # Unstranded data - merge all together
        bedtools unionbedg -i *.sorted.bedgraph | \\
            awk 'BEGIN {OFS="\\t"} {sum=0; for(i=4; i<=NF; i++) sum+=\$i; print \$1, \$2, \$3, sum}' | \\
            sort -k1,1 -k2,2n > ${prefix}.merged.sorted.bedgraph
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def has_stranded = bedgraphs.any { it.name.contains('.forward.') || it.name.contains('.reverse.') }

    if (has_stranded) {
        """
        touch ${prefix}.merged.forward.sorted.bedgraph
        touch ${prefix}.merged.reverse.sorted.bedgraph

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: 2.31.1
        END_VERSIONS
        """
    } else {
        """
        touch ${prefix}.merged.sorted.bedgraph

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: 2.31.1
        END_VERSIONS
        """
    }
}
