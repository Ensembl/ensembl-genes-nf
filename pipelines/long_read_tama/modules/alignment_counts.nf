process ALIGNMENT_COUNTS {
    tag "${meta.id}:alignment-counts"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1'
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path('alignment_counts.tsv'), emit: stats
    path 'flagstat.txt', emit: flagstat
    path 'idxstats.txt', emit: idxstats
    path 'versions.yml', emit: versions
    script:
    """
    samtools quickcheck -v "${bam}"
    samtools flagstat "${bam}" > flagstat.txt
    samtools idxstats "${bam}" > idxstats.txt
    awk -v sample='${meta.id}' '
        / in total/ { total = \$1 }
        / mapped \(/ && !/primary mapped/ { mapped = \$1 }
        \$NF == "primary" { primary = \$1 }
        \$NF == "secondary" { secondary = \$1 }
        \$NF == "supplementary" { supplementary = \$1 }
        / unmapped \(/ { unmapped = \$1; saw_unmapped = 1 }
        END {
            if (!saw_unmapped) unmapped = total - mapped
            printf "sample_id\\ttotal_alignments\\tmapped_alignments\\tprimary_alignments\\tsecondary_alignments\\tsupplementary_alignments\\tunmapped_alignments\\tmapping_rate_percent\\tprimary_mapping_rate_percent\\n"
            printf "%s\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\t%.6f\\t%.6f\\n", sample, total, mapped, primary, secondary, supplementary, unmapped, (total ? 100 * mapped / total : 0), (total ? 100 * primary / total : 0)
        }
    ' flagstat.txt > alignment_counts.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1)
    END_VERSIONS
    """
    stub:
    """
    printf 'sample_id\\ttotal_alignments\\tmapped_alignments\\tprimary_alignments\\tsecondary_alignments\\tsupplementary_alignments\\tunmapped_alignments\\tmapping_rate_percent\\tprimary_mapping_rate_percent\\n${meta.id}\\t10\\t9\\t9\\t0\\t0\\t1\\t90.000000\\t90.000000\\n' > alignment_counts.tsv
    touch flagstat.txt idxstats.txt
    printf '"%s":\\n    samtools: 1.20\\n' '${task.process}' > versions.yml
    """
}
