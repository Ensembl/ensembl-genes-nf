process STAGE_ANNOTATION_PAIR {
    tag { meta.id }
    label 'process_low'

    container "https://depot.galaxyproject.org/singularity/curl:8.8.0--he654da7_1"

    publishDir "${params.outdir}/qc/pairwise_annotation",
        mode: 'copy',
        pattern: "*"

    input:
        tuple val(meta), val(source_a_ref), val(source_b_ref), val(assembly_report_ref)
        val cache_dir

    output:
        tuple val(meta), path("source_a.gff3"), path("source_b.gff3"), path("assembly_report.txt"), emit: staged
        path "versions.yml", emit: versions

    script:
        def cache = cache_dir ?: "${params.outdir}/cache/pairwise_annotation"
        """
        set -euo pipefail

        mkdir -p "${cache}"

        fetch_or_link() {
            src="\$1"
            dest="\$2"

            if [[ "\$src" =~ ^https?://|^ftp:// ]]; then
                base=\$(basename "\$src")
                cached="${cache}/\$base"
                if [ ! -s "\$cached" ]; then
                    curl -L --fail --retry 3 --retry-delay 5 -o "\$cached" "\$src"
                fi
                ln -s "\$cached" "\$dest"
            elif [ "\$src" = "" ] || [ "\$src" = "NA" ] || [ "\$src" = "null" ]; then
                : > "\$dest"
            else
                ln -s "\$(readlink -f "\$src")" "\$dest"
            fi
        }

        fetch_or_link "${source_a_ref}" source_a.gff3
        fetch_or_link "${source_b_ref}" source_b.gff3
        fetch_or_link "${assembly_report_ref}" assembly_report.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            curl: \$(curl --version | head -n1 | sed 's/curl //g' | cut -d' ' -f1)
        END_VERSIONS
        """

    stub:
        """
        touch source_a.gff3
        touch source_b.gff3
        touch assembly_report.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            curl: 8.8.0
        END_VERSIONS
        """
}
