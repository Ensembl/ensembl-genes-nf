process INDEX_AND_VALIDATE_BAMS {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_medium'
    input:
    tuple val(meta), path(staged)
    val stage_dir

    output:
    tuple val(meta), path('staged'), emit: indexed
    path 'bam_validation.tsv', emit: validation
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    genome=\$(find staged -maxdepth 1 -type f -name 'genome.fa*' | head -n1)
    test -n "\$genome" && test -s "\$genome"
    samtools faidx "\$genome"
    printf 'bam\\tstatus\\n' > bam_validation.tsv
    for bam in staged/bam_*.bam; do
        test -s "\$bam"
        samtools quickcheck -v "\$bam"
        samtools view -H "\$bam" | grep -q 'SO:coordinate' || { echo "BAM is not coordinate sorted: \$bam" >&2; exit 1; }
        bai="\${bam}.bai"
        if [ ! -s "\$bai" ]; then samtools index -@ ${task.cpus} "\$bam" "\$bai"; fi
        samtools idxstats "\$bam" >/dev/null
        printf '%s\\tVALID\\n' "\$(basename "\$bam")" >> bam_validation.tsv
    done
    python3 - "\$genome" <<'PY'
import gzip, sys
path = sys.argv[1]
opener = gzip.open if path.endswith('.gz') else open
lower = 0
with opener(path, 'rt') as handle:
    for line in handle:
        if not line.startswith('>'):
            lower += sum(char.islower() for char in line)
if lower == 0:
    raise SystemExit('Genome appears not to be soft-masked: no lowercase bases found')
PY
    cat > versions.yml <<VERSIONS
"\${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """

    stub:
    """
    mkdir -p staged
    touch staged/bam_01.bam.bai staged/bam_02.bam.bai
    printf 'bam\\tstatus\\n bam_01.bam\\tVALID\\n bam_02.bam\\tVALID\\n' > bam_validation.tsv
    touch versions.yml
    """
}
