process RESOLVE_BAM_REFERENCE {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1'
    publishDir "${params.outdir}/reference", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(bam), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.gca_accession}.${meta.safe_id}.chrom.sizes"), emit: resolved
    path "${meta.gca_accession}.${meta.safe_id}.reference_resolution.tsv", emit: report
    path "${meta.gca_accession}.${meta.safe_id}.reference_resolution.versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    samtools view -H ${bam} | awk -F '\\t' '
        \$1 == "@SQ" {
            name = ""
            seq_length = ""
            for (i = 2; i <= NF; i++) {
                if (\$i ~ /^SN:/) name = substr(\$i, 4)
                if (\$i ~ /^LN:/) seq_length = substr(\$i, 4)
            }
            if (name == "" || seq_length !~ /^[0-9]+\$/ || seq_length < 1) {
                print "Invalid BAM sequence dictionary row: " \$0 > "/dev/stderr"
                bad = 1
            }
            if (seen[name]++) {
                print "Duplicate BAM sequence name: " name > "/dev/stderr"
                bad = 1
            }
            print name "\\t" seq_length
            rows++
        }
        END { if (bad || rows == 0) exit 1 }
    ' > ${meta.gca_accession}.${meta.safe_id}.bam.chrom.sizes

    awk '
        NF != 2 || \$1 == "" || \$2 !~ /^[0-9]+\$/ || \$2 < 1 {
            print "Invalid chrom.sizes row: " \$0 > "/dev/stderr"
            bad = 1
        }
        seen[\$1]++ { print "Duplicate chromosome in chrom.sizes: " \$1 > "/dev/stderr"; bad = 1 }
        END { if (bad || NR == 0) exit 1 }
    ' ${chrom_sizes}

    LC_ALL=C sort -k1,1 -k2,2n ${chrom_sizes} > ${meta.gca_accession}.${meta.safe_id}.reference.chrom.sizes
    LC_ALL=C sort -k1,1 -k2,2n ${meta.gca_accession}.${meta.safe_id}.bam.chrom.sizes > ${meta.gca_accession}.${meta.safe_id}.bam.sorted.chrom.sizes

    if cmp -s ${meta.gca_accession}.${meta.safe_id}.reference.chrom.sizes ${meta.gca_accession}.${meta.safe_id}.bam.sorted.chrom.sizes; then
        cp ${meta.gca_accession}.${meta.safe_id}.reference.chrom.sizes ${meta.gca_accession}.${meta.safe_id}.chrom.sizes
        resolution='input_names_match'
    elif cmp -s \\
        <(cut -f1 ${meta.gca_accession}.${meta.safe_id}.reference.chrom.sizes) \\
        <(cut -f1 ${meta.gca_accession}.${meta.safe_id}.bam.sorted.chrom.sizes); then
        echo 'BAM and chrom.sizes names match but sequence lengths differ.' >&2
        exit 1
    elif diff -q \\
        <(cut -f2 ${meta.gca_accession}.${meta.safe_id}.reference.chrom.sizes | LC_ALL=C sort -n) \
        <(cut -f2 ${meta.gca_accession}.${meta.safe_id}.bam.sorted.chrom.sizes | LC_ALL=C sort -n) >/dev/null; then
        cp ${meta.gca_accession}.${meta.safe_id}.bam.sorted.chrom.sizes ${meta.gca_accession}.${meta.safe_id}.chrom.sizes
        resolution='bam_names_used_lengths_match'
        echo 'BAM sequence names differ from chrom.sizes; using BAM names because all sequence lengths match.' >&2
    else
        echo 'BAM sequence lengths do not match chrom.sizes; refusing to create tracks against a different reference.' >&2
        exit 1
    fi

    printf 'gca_accession\\tassembly_release\\tentity_id\\tresolution\\tchromosome_count\\n' > ${meta.gca_accession}.${meta.safe_id}.reference_resolution.tsv
    printf '%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        '${meta.gca_accession}' '${meta.assembly_release}' '${meta.id}' \\
        "\$resolution" "\$(wc -l < ${meta.gca_accession}.${meta.safe_id}.chrom.sizes | tr -d ' ')" \\
        >> ${meta.gca_accession}.${meta.safe_id}.reference_resolution.tsv

    printf '"%s":\\n    samtools: %s\\n    reference_resolution: %s\\n' \\
        '${task.process}' "\$(samtools --version | head -n1 | awk '{print \$2}')" \\
        "\$resolution" > ${meta.gca_accession}.${meta.safe_id}.reference_resolution.versions.yml
    """

    stub:
    """
    cp ${chrom_sizes} ${meta.gca_accession}.${meta.safe_id}.chrom.sizes
    printf 'gca_accession\\tassembly_release\\tentity_id\\tresolution\\tchromosome_count\\n${meta.gca_accession}\\t${meta.assembly_release}\\t${meta.id}\\tstub\\t0\\n' > ${meta.gca_accession}.${meta.safe_id}.reference_resolution.tsv
    printf '"%s":\\n    samtools: stub\\n    reference_resolution: stub\\n' '${task.process}' > ${meta.gca_accession}.${meta.safe_id}.reference_resolution.versions.yml
    """
}
