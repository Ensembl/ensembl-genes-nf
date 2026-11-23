process FILTER_PROFILE_ELIGIBLE_BAM {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.getrpf_container ?: 'ghcr.io/jackcurragh/get-rpf:sha-e9dce4a'}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.profile_eligible.bam"), path("*.profile_eligible.bam.bai"), emit: bam_and_bai
    tuple val(meta), path("*.profile_eligibility.json"), emit: stats

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - ${bam} ${prefix}.profile_eligible.bam ${prefix}.profile_eligibility.json <<'PY'
import json
import sys
import pysam

source, target, report = sys.argv[1:]
stats = {"mapped_reads": 0, "eligible_reads": 0, "excluded_5p_clip": 0,
         "excluded_length": 0, "excluded_secondary": 0}
with pysam.AlignmentFile(source, "rb") as inp, pysam.AlignmentFile(target, "wb", template=inp) as out:
    for read in inp.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        count = int(read.query_name.rsplit("_x", 1)[1]) if "_x" in read.query_name else 1
        stats["mapped_reads"] += count
        if read.is_secondary or read.is_supplementary or read.cigartuples is None:
            stats["excluded_secondary"] += count
            continue
        cigar = list(read.cigartuples)
        aligned_length = sum(n for op, n in cigar if op in (0, 1, 7, 8))
        oriented = list(reversed(cigar)) if read.is_reverse else cigar
        clip5 = 0
        for op, n in oriented:
            if op in (5, 6):
                continue
            clip5 = n if op == 4 else 0
            break
        if clip5:
            stats["excluded_5p_clip"] += count
            continue
        if not 20 <= aligned_length <= 40:
            stats["excluded_length"] += count
            continue
        # Replace soft clips with hard clips and retain only aligned query bases.
        # This keeps the provider's read length equal to M/=X/I consumption.
        aligned_sequence = read.query_alignment_sequence
        aligned_qualities = read.query_alignment_qualities
        if aligned_sequence is None or len(aligned_sequence) != aligned_length:
            stats["excluded_length"] += count
            continue
        read.cigartuples = [(5 if op == 4 else op, n) for op, n in cigar]
        read.query_sequence = aligned_sequence
        if aligned_qualities is not None:
            read.query_qualities = aligned_qualities
        out.write(read)
        stats["eligible_reads"] += count

pysam.index(target)
stats["eligible_fraction"] = (stats["eligible_reads"] / stats["mapped_reads"]
                              if stats["mapped_reads"] else 0.0)
with open(report, "w") as handle:
    json.dump(stats, handle, indent=2)
    handle.write("\\n")
PY
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.profile_eligible.bam
    touch ${prefix}.profile_eligible.bam.bai
    cat > ${prefix}.profile_eligibility.json <<'PROFILE_STATS'
    {"mapped_reads":0,"eligible_reads":0,"eligible_fraction":0.0}
    PROFILE_STATS
    """
}
