
process STANDARDISE_BED12 {
    label 'process_high'

    container "oras://community.wave.seqera.io/library/pybedtools_pysam_pip_biopython:1dbd8151223e518c"
    
    tag "${meta.id}"

    publishDir "${params.outdir}/standardised_bed12s/${tool}", mode: 'copy'

    input:
    tuple val(meta), val(tool), path(bed_file)
    path genome_fasta
    path genome_fasta_fai

    output:
    tuple val(meta), val(tool), path("*.valid.bed12")

    script:    
    """
    mkdir -p out
    # Ensure converter is available on PATH
    if ! command -v translonscorer_to_bed12.py >/dev/null 2>&1; then
        cp ${projectDir}/pipelines/translon-consensus/bin/translonscorer_to_bed12.py /usr/local/bin/
        chmod +x /usr/local/bin/translonscorer_to_bed12.py || true
    fi

    # If input is a TranslonScorer CSV, convert first
    if echo "${bed_file}" | grep -qi '\\.csv$'; then
        translonscorer_to_bed12.py ${bed_file} out/
        inbed=$(ls out/*.bed12 | head -n1)
    else
        inbed=${bed_file}
    fi

    standardise_bed12.py \\
        -i ${inbed} \\
        -f ${genome_fasta} \\
        --output_prefix ${meta.id} --verbose
    """
    
    stub:
    """
    mkdir -p standardised_bed12s
    touch standardised_bed12s/${meta.id}_0based.bed12
    """
}
