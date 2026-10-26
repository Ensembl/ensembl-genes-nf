process MINIMAP2_ALIGN {
    publishDir path: { "${params.outdir}/per_run/${meta.id}/alignment" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id}:alignment"
    label 'process_high'
    // Keep minimap2 and samtools in one image so SAM streams directly into
    // samtools sort instead of being materialized as a large work-dir file.
    container 'docker://niemasd/minimap2_samtools:2.28_1.20'

    input:
    tuple val(meta), path(reads)
    path minimap_index

    output:
    tuple val(meta), path('*.sorted.bam'), path('*.sorted.bam.bai'), emit: bam
    tuple val(meta), path("${meta.id}_alignment_config.tsv"), emit: config
    path 'versions.yml', emit: versions

    script:
    def secondary = task.ext.secondary ?: (params.secondary_mode == 'yes' ? '--secondary=yes' : '--secondary=no')
    def preset = meta.minimap2_preset ?: params.minimap2_preset
    def preset_source = meta.minimap2_preset ? 'manifest' : 'params_fallback'
    def minimap2_args = meta.minimap2_args ?: ''
    def args_source = meta.containsKey('minimap2_args') ? 'manifest' : 'empty_fallback'
    def technology = meta.technology ?: (meta.classification?.toString()?.startsWith('ONT_') ? 'ont' : 'pacbio')
    def representation = meta.read_representation ?: (technology == 'ont' ? 'ont_cdna' : 'ccs')
    def preprocessing = meta.preprocessing_state ?: 'processed'
    def isoquant_type = meta.containsKey('isoquant_data_type') ? (meta.isoquant_data_type ?: 'not_eligible') : (technology == 'ont' ? 'nanopore' : 'pacbio_ccs')
    def observation_unit = meta.observation_unit ?: (technology == 'ont' ? 'ont_read' : 'consensus_molecule')
    """
    minimap2 -t ${task.cpus} -ax ${preset} ${minimap2_args} ${secondary} --cs=long "${minimap_index}" "${reads}" \\
        | samtools sort -@ ${task.cpus} -m ${params.samtools_sort_memory} -O BAM -o "${meta.id}.sorted.bam" -
    test -s "${meta.id}.sorted.bam" || { echo "alignment produced no BAM for ${meta.id}" >&2; exit 1; }
    samtools index "${meta.id}.sorted.bam"
    test -s "${meta.id}.sorted.bam.bai" || { echo "alignment produced no BAI for ${meta.id}" >&2; exit 1; }
    printf 'run_accession\\ttechnology\\tread_representation\\tpreprocessing_state\\tminimap2_preset\\tminimap2_preset_source\\tminimap2_args\\tminimap2_args_source\\tsecondary_alignment\\tcs_tag\\tisoquant_data_type\\tobservation_unit\\n${meta.id}\\t${technology}\\t${representation}\\t${preprocessing}\\t${preset}\\t${preset_source}\\t${minimap2_args}\\t${args_source}\\t${secondary}\\t--cs=long\\t${isoquant_type}\\t${observation_unit}\\n' > ${meta.id}_alignment_config.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -n1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai
    printf 'run_accession\\ttechnology\\tread_representation\\tpreprocessing_state\\tminimap2_preset\\tminimap2_preset_source\\tminimap2_args\\tminimap2_args_source\\tsecondary_alignment\\tcs_tag\\tisoquant_data_type\\tobservation_unit\\n${meta.id}\\t${meta.technology ?: 'pacbio'}\\t${meta.read_representation ?: 'ccs'}\\t${meta.preprocessing_state ?: 'processed'}\\t${meta.minimap2_preset ?: params.minimap2_preset}\\t${meta.minimap2_preset ? 'manifest' : 'params_fallback'}\\t${meta.minimap2_args ?: ''}\\t${meta.containsKey('minimap2_args') ? 'manifest' : 'empty_fallback'}\\t${params.secondary_mode == 'yes' ? '--secondary=yes' : '--secondary=no'}\\t--cs=long\\t${meta.isoquant_data_type ?: (meta.read_representation == 'raw_subread' ? 'not_eligible' : (meta.technology == 'ont' ? 'nanopore' : 'pacbio_ccs'))}\\t${meta.observation_unit ?: 'consensus_molecule'}\\n' > ${meta.id}_alignment_config.tsv
    printf '"%s":\n    minimap2: 2.28\n    samtools: 1.20\n' '${task.process}' > versions.yml
    """
}
