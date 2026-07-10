process RIBOMETRIC {
    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.10 conda-forge::biopython bioconda::pysam"
    container "ghcr.io/jackcurragh/ribometric:1.4.2"

    publishDir "${params.outdir}/RiboMetric", mode: 'copy', pattern: "*RiboMetric.{html,json,csv}"
    publishDir "${params.outdir}/RiboMetric/offsets", mode: 'copy', pattern: "*.{offsets.tsv,best_offset.txt}"

    input:
    tuple val(meta), path(transcriptome_bam), path(transcriptome_bam_index)
    path ribometric_annotation
    path offset_file  // Optional: external offset file (e.g., calculated offsets)

    output:
    tuple val(meta), path("*RiboMetric.html"), emit: html
    tuple val(meta), path("*RiboMetric.json"), emit: json
    tuple val(meta), path("*RiboMetric.csv"), emit: csv
    tuple val(meta), path("*.offsets.tsv"), optional: true, emit: offsets_audit
    tuple val(meta), path("*.best_offset.txt"), emit: offsets
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_size = params.ribometric_sample_size ?: 10000000

    // Determine offset configuration
    // Priority: 1) external offset file, 2) global offset, 3) calculation method
    def offset_args = ''
    if (offset_file && offset_file.name != 'NO_OFFSET_FILE') {
        offset_args = "--offset-read-length ${offset_file}"
    } else if (params.ribometric_offset_global) {
        offset_args = "--offset-global ${params.ribometric_offset_global}"
    } else {
        def offset_method = params.ribometric_offset_method ?: 'tripsviz'
        offset_args = "--offset-calculation-method ${offset_method}"
    }
    // Only output offsets when calculating internally (not when using an external file)
    def output_offsets_arg = (offset_file && offset_file.name != 'NO_OFFSET_FILE') ? '' : "--output-offsets ${prefix}.offsets.tsv"
    def offset_target = params.ribometric_offset_target ?: 'a_site'
    """
    RiboMetric run \\
        --bam ${transcriptome_bam} \\
        --annotation ${ribometric_annotation} \\
        --threads $task.cpus \\
        --html \\
        --json \\
        --csv \\
        --offset-target ${offset_target} \\
        ${offset_args} \\
        ${output_offsets_arg} \\
        -S ${sample_size} \\
        $args

    if [ "${offset_file.name}" != "NO_OFFSET_FILE" ]; then
        cp ${offset_file} ${prefix}.best_offset.txt
    else
        python3 - <<'PY'
    import csv
    from pathlib import Path

    audit_path = Path("${prefix}.offsets.tsv")
    out_path = Path("${prefix}.best_offset.txt")

    def clean_int(value):
        if value is None:
            return None
        value = str(value).strip()
        if not value:
            return None
        try:
            return int(float(value))
        except ValueError:
            return None

    with audit_path.open(newline="") as handle, out_path.open("w") as out:
        reader = csv.DictReader(handle, delimiter="\\t")
        fields = set(reader.fieldnames or [])
        out.write("length\\toffset\\n")

        if {"read_len", "offset"}.issubset(fields):
            for row in reader:
                read_length = clean_int(row.get("read_len"))
                offset = clean_int(row.get("offset"))
                if read_length is not None and offset is not None:
                    out.write(f"{read_length}\\t{offset}\\n")
        else:
            for row in reader:
                read_length = clean_int(row.get("read_length"))
                if read_length is None:
                    continue

                offset = clean_int(row.get("new_offset"))
                if offset is None:
                    offset = clean_int(row.get("computed_offset"))
                if offset is None:
                    applied = [clean_int(v) for v in str(row.get("applied_offsets", "")).split("|")]
                    applied = [v for v in applied if v is not None]
                    if len(set(applied)) == 1:
                        offset = applied[0]
                if offset is None:
                    offset = clean_int(row.get("min_offset"))

                if offset is not None:
                    out.write(f"{read_length}\\t{offset}\\n")
    PY
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribometric: \$(RiboMetric --version 2>&1 | sed 's/RiboMetric version //g' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_RiboMetric.html
    touch ${prefix}_RiboMetric.json
    touch ${prefix}_RiboMetric.csv
    cat <<-END_OFFSETS > ${prefix}.offsets.tsv
    sample	offset_source	offset_target	offset_calculation_method	read_length	n_reads	n_unique_offsets	applied_offsets	min_offset	max_offset	computed_offset	global_offset	frame_adjusted	old_offset	new_offset	dominant_frame	dominant_fraction	frame_adjustment_reads
    ${prefix}	calculated	a_site	tripsviz	28	1	1	15	15	15	15		False
    END_OFFSETS
    cat <<-END_BEST_OFFSET > ${prefix}.best_offset.txt
    length	offset
    28	15
    END_BEST_OFFSET

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribometric: 1.4.2
    END_VERSIONS
    """
}
