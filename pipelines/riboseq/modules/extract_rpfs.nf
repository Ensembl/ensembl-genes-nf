process EXTRACT_RPFS {
    tag "${meta.id}"
    label 'process_high'

    conda "conda-forge::python=3.10 conda-forge::biopython"
    container "ghcr.io/jackcurragh/get-rpf:sha-e9dce4a"

    input:
    tuple val(meta), path(input_file)
    path star_index

    output:
    // No reads when getRPF withholds the transform (the read structure could
    // not be resolved): the sample stops here, but its reports are published.
    tuple val(meta), path("*.collapsed.fa"), emit: trimmed_collapsed, optional: true
    tuple val(meta), path("*.seqspec.yaml"), emit: seqspec, optional: true
    tuple val(meta), path("*.extraction_report.json"), emit: report
    tuple val(meta), path("*.structure.{json,txt}"), emit: structure, optional: true
    tuple val(meta), path("*_getrpf_hold.fastq.gz"), emit: getrpf_hold_reads, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // 'structure': infer the read structure from the first getrpf_infer_reads
    // reads, then apply it to every read (--star-index adds an alignment check).
    // 'legacy': the pre-0.4 architecture matching and boundary estimation.
    def structure = params.getrpf_method == 'legacy' ? '' : "--infer-structure --infer-reads ${params.getrpf_infer_reads}"
    def hold_export = params.align_getrpf_holds ? """
    python3 - <<'PY'
import gzip
import json
import shutil
from pathlib import Path

structure_path = Path("${prefix}_trimmed.structure.json")
source_path = Path("${input_file}")
target_path = Path("${prefix}_getrpf_hold.fastq.gz")

if structure_path.exists():
    with structure_path.open() as handle:
        structure = json.load(handle)
else:
    structure = {}

if structure and not (structure.get("transform") or {}).get("emit", False):
    source = gzip.open(source_path, "rt") if source_path.name.endswith(".gz") else source_path.open("rt")
    try:
        with gzip.open(target_path, "wt") as target:
            shutil.copyfileobj(source, target)
    finally:
        source.close()
PY
    """ : ''
    """
    getRPF extract \\
        ${input_file} \\
        ${prefix}_trimmed.fastq \\
        -f fastq \\
        --generate-seqspec \\
        --output-format json \\
        --star-index ${star_index} \\
        --star-threads ${task.cpus} \\
        --collapsed-only \\
        ${structure} \\
        $args

    $hold_export

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getRPF: \$(getRPF --version 2>&1 | head -n1 | sed 's/^.*version //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_trimmed.collapsed.fa
    touch ${prefix}_trimmed.seqspec.yaml
    touch ${prefix}_trimmed.extraction_report.json
    touch ${prefix}_trimmed.structure.json
    touch ${prefix}_trimmed.structure.txt
    ${params.align_getrpf_holds ? "touch ${prefix}_getrpf_hold.fastq.gz" : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getRPF: 0.5.0
    END_VERSIONS
    """
}
