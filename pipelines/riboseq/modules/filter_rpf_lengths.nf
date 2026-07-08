process FILTER_RPF_LENGTHS {
    tag "${meta.id}"
    label 'process_light'

    conda "conda-forge::python=3.10"
    container "ghcr.io/jackcurragh/get-rpf:0.2.2"

    publishDir "${params.outdir}/getRPF/gated", mode: 'copy', pattern: "*.{collapsed.fa,summary.tsv}"

    input:
    tuple val(meta), path(trimmed_collapsed_fasta)

    output:
    tuple val(meta), path("*_rpf_20_40.collapsed.fa"), emit: collapsed_fasta
    tuple val(meta), path("*_rpf_20_40.summary.tsv"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - <<'PY'
from pathlib import Path

input_path = Path("${trimmed_collapsed_fasta}")
output_path = Path("${prefix}_rpf_20_40.collapsed.fa")
summary_path = Path("${prefix}_rpf_20_40.summary.tsv")

def parse_count(header: str) -> int:
    if "_x" in header:
        try:
            return int(header.rsplit("_x", 1)[1].split()[0])
        except ValueError:
            return 1
    return 1

total_unique = total_reads = kept_unique = kept_reads = 0
current_header = None
current_seq = []

with input_path.open() as handle, output_path.open("w") as out:
    def flush():
        global total_unique, total_reads, kept_unique, kept_reads
        if current_header is None:
            return
        seq = "".join(current_seq).strip().upper()
        if not seq:
            return
        count = parse_count(current_header)
        total_unique += 1
        total_reads += count
        if 20 <= len(seq) <= 40:
            kept_unique += 1
            kept_reads += count
            out.write(">{}\\n{}\\n".format(current_header, seq))

    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            current_header = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    flush()

with summary_path.open("w") as summary:
    summary.write("sample\\ttotal_unique\\ttotal_reads\\tkept_unique_20_40\\tkept_reads_20_40\\tkept_read_fraction\\n")
    fraction = kept_reads / total_reads if total_reads else 0
    summary.write("{}\\t{}\\t{}\\t{}\\t{}\\t{:.6f}\\n".format(
        "${prefix}", total_unique, total_reads, kept_unique, kept_reads, fraction
    ))
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rpf_20_40.collapsed.fa
    cat > ${prefix}_rpf_20_40.summary.tsv <<'EOF'
sample	total_unique	total_reads	kept_unique_20_40	kept_reads_20_40	kept_read_fraction
${prefix}	0	0	0	0	0.000000
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10.0
    END_VERSIONS
    """
}
