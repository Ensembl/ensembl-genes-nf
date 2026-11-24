process EXTRACT_RPFS {
    tag "${meta.id}"
    label 'process_high'

    conda "conda-forge::python=3.10 conda-forge::biopython"
    container "${params.getrpf_container ?: 'ghcr.io/jackcurragh/get-rpf:sha-e9dce4a'}"

    input:
    tuple val(meta), path(input_file), val(trim_recipe_json)
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
    def prepare = params.getrpf_prepare_for_alignment ? '--prepare-for-alignment' : ''
    def structure = (prepare || params.getrpf_method == 'legacy') ? '' : "--infer-structure --infer-reads ${params.getrpf_infer_reads}"
    // Preparation skips its diagnostic STAR pass; the production BAM is the
    // sole alignment and supplies clipping evidence downstream.
    def star_check = prepare ? '' : "--star-index ${star_index} --star-threads ${task.cpus}"
    def recipe_encoded = trim_recipe_json && trim_recipe_json != 'null'
        ? trim_recipe_json.getBytes('UTF-8').encodeBase64().toString()
        : ''
    def recipe_write = recipe_encoded ? "python3 -c \"import base64; from pathlib import Path; Path('${prefix}.trim_recipe.json').write_bytes(base64.b64decode('${recipe_encoded}'))\"" : ''
    def recipe_arg = trim_recipe_json && trim_recipe_json != 'null' ? "--trim-recipe ${prefix}.trim_recipe.json" : ''
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
    $recipe_write

    # Use the recommended 'extract' command (alignment-based extraction)
    # Keep the default collapsed-only output; optional contaminant filtering
    # expands this already-trimmed output in a separate process.
    getRPF extract \\
        ${input_file} \\
        ${prefix}_trimmed.fastq \\
        -f fastq \\
        --generate-seqspec \\
        --output-format json \\
        $star_check \\
        --collapsed-only \\
        $prepare \\
        $recipe_arg \\
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
    cat > ${prefix}_trimmed.extraction_report.json <<'GETRPF_STUB_REPORT'
    {"schema":"getrpf.structure.extraction/1","status":"prepared","preparation":{"requires_local_alignment":false}}
    GETRPF_STUB_REPORT
    touch ${prefix}_trimmed.structure.json
    touch ${prefix}_trimmed.structure.txt
    ${params.align_getrpf_holds ? "touch ${prefix}_getrpf_hold.fastq.gz" : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getRPF: 0.6.0
    END_VERSIONS
    """
}
