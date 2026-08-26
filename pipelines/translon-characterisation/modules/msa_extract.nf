process MSA_EXTRACT {
    label 'process_medium'
    // Requires mafExtract (Kent tools) plus Python; production supplies an immutable image.
    container params.msa_extract_container
    tag "${meta.id}"
    publishDir "${params.outdir}/06_msa", mode: 'copy', pattern: '*.msa.fasta'

    input:
    tuple val(meta), path(instances)
    path maf_dir

    output:
    tuple val(meta), path('msas'), path('msa_identity.jsonl'), emit: alignments
    tuple val(meta), path('msa_extract.axis.jsonl'), emit: axis
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    extract_msas.py --instances ${instances} --maf-dir ${maf_dir} --outdir msas --axis-output msa_extract.axis.jsonl --identity-output msa_identity.jsonl --reference ${params.alignment_reference} --species-map ${params.phylocsf_species_map}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafExtract: \$(mafExtract 2>&1 | head -n1 || echo unknown)
        stitch_maf_alignment: 1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    mkdir msas
    printf '>Human\nATGTAA\n' > msas/${prefix}.msa.fasta
    printf '{"instance_id":"stub","alignment":"${prefix}.msa.fasta"}\n' > msa_identity.jsonl
    touch msa_extract.axis.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafExtract: stub
        stitch_maf_alignment: 1
    END_VERSIONS
    """
}
