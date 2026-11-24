/*
 * EXPAND_COLLAPSED_FASTQ
 * Recreate a FASTQ stream from getRPF's already-trimmed collapsed FASTA.
 * This is used only for optional contaminant filtering, before the reads are
 * collapsed again after filtering.
 */

process EXPAND_COLLAPSED_FASTQ {
    tag "${meta.id}"
    label 'process_high'

    conda "conda-forge::python=3.10"

    input:
    tuple val(meta), path(collapsed_fasta)

    output:
    tuple val(meta), path("*_trimmed.fastq"), emit: trimmed_fastq
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - <<'PY'
from pathlib import Path

input_path = Path("${collapsed_fasta}")
output_path = Path("${prefix}_trimmed.fastq")

def parse_count(header: str) -> int:
    if "_x" in header:
        try:
            return int(header.rsplit("_x", 1)[1].split()[0])
        except ValueError:
            pass
    return 1

header = None
sequence = []

with input_path.open() as source, output_path.open("w") as output:
    def flush():
        if header is None:
            return
        seq = "".join(sequence).strip().upper()
        if not seq:
            return
        for copy in range(parse_count(header)):
            read_id = f"{header}_r{copy + 1}"
            output.write(f"@{read_id}\n{seq}\n+\n{'I' * len(seq)}\n")

    for raw_line in source:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header = line[1:]
            sequence = []
        else:
            sequence.append(line)
    flush()
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_trimmed.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10.0
    END_VERSIONS
    """
}
