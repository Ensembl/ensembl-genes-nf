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
    // Convert bedgraphs to list if it's a single file
    def bg_list = bedgraphs instanceof List ? bedgraphs : [bedgraphs]

    // Check if we have stranded bedgraphs (forward/reverse pairs)
    def has_stranded = bg_list.any { it.name.contains('.forward.') || it.name.contains('.reverse.') }

    if (has_stranded) {
        // Separate forward and reverse bedgraphs
        def forward_bgs = bg_list.findAll { it.name.contains('.forward.') }
        def reverse_bgs = bg_list.findAll { it.name.contains('.reverse.') }

        // Create space-separated file lists
        def forward_files = forward_bgs.collect { it.toString() }.join(' ')
        def reverse_files = reverse_bgs.collect { it.toString() }.join(' ')

        """
        # Merge forward strand bedgraphs
        bedtools unionbedg -i ${forward_files} | \\
            awk 'BEGIN {OFS="\\t"} {sum=0; for(i=4; i<=NF; i++) sum+=\$i; print \$1, \$2, \$3, sum}' | \\
            sort -k1,1 -k2,2n > ${prefix}.merged.forward.sorted.bedgraph

        # Merge reverse strand bedgraphs
        bedtools unionbedg -i ${reverse_files} | \\
            awk 'BEGIN {OFS="\\t"} {sum=0; for(i=4; i<=NF; i++) sum+=\$i; print \$1, \$2, \$3, sum}' | \\
            sort -k1,1 -k2,2n > ${prefix}.merged.reverse.sorted.bedgraph

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed 's/bedtools v//')
        END_VERSIONS
        """
    } else {
        // Unstranded bedgraphs
        def all_files = bg_list.collect { it.toString() }.join(' ')
        """
        bedtools unionbedg -i ${all_files} | \\
            awk 'BEGIN {OFS="\\t"} {sum=0; for(i=4; i<=NF; i++) sum+=\$i; print \$1, \$2, \$3, sum}' | \\
            sort -k1,1 -k2,2n > ${prefix}.merged.sorted.bedgraph

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed 's/bedtools v//')
        END_VERSIONS
        """
    }

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
