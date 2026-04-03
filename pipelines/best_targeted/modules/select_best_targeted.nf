// SELECT_BEST_TARGETED
// Merge cDNA and protein exonerate GFF3 results (after per-batch filtering),
// cluster overlapping models, and pick the best non-redundant set.

process SELECT_BEST_TARGETED {
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir "${params.outdir}/best_targeted", mode: 'copy', pattern: "*.best_targeted.gff3"

    input:
    path cdna_gff3s    // merged list of per-batch cdna filtered GFF3 files (or [])
    path protein_gff3s // merged list of per-batch protein filtered GFF3 files (or [])

    output:
    path "best_targeted.gff3", emit: gff3
    path "versions.yml",       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def min_cov    = params.best_targeted_min_coverage ?: 50
    def min_pid    = params.best_targeted_min_pid      ?: 50
    def cdna_arg   = cdna_gff3s    ? "--cdna_gff3 ${cdna_gff3s instanceof List ? cdna_gff3s.join(' ') : cdna_gff3s}" : ''
    def prot_arg   = protein_gff3s ? "--protein_gff3 ${protein_gff3s instanceof List ? protein_gff3s.join(' ') : protein_gff3s}" : ''
    """
    # Concatenate per-batch files if needed
    if [ -n "${cdna_gff3s ? cdna_gff3s : ''}" ]; then
        cat ${cdna_gff3s} > merged_cdna.gff3
    fi
    if [ -n "${protein_gff3s ? protein_gff3s : ''}" ]; then
        cat ${protein_gff3s} > merged_protein.gff3
    fi

    select_best_targeted.py \\
        ${cdna_gff3s    ? '--cdna_gff3 merged_cdna.gff3'       : ''} \\
        ${protein_gff3s ? '--protein_gff3 merged_protein.gff3'  : ''} \\
        --out         best_targeted.gff3 \\
        --min_coverage ${min_cov} \\
        --min_pid      ${min_pid} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf '##gff-version 3\\n' > best_targeted.gff3
    printf 'chr1\\texonerate\\tgene\\t1000\\t5000\\t500\\t+\\t.\\tID=stub_bt_gene_00000001;Name=NM_001;biotype=best_targeted;coverage=95.0;pid=97.0\\n' \\
        >> best_targeted.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
