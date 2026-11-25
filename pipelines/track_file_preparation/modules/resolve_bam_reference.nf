process RESOLVE_BAM_REFERENCE {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1'
    publishDir "${params.outdir}/reference", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(bam), path(chrom_sizes)
    path resolver_script

    output:
    tuple val(meta), path("${meta.gca_accession}.${meta.safe_id}.chrom.sizes"), emit: resolved
    path "${meta.gca_accession}.${meta.safe_id}.reference_resolution.tsv", emit: report
    path "${meta.gca_accession}.${meta.safe_id}.reference_resolution.versions.yml", emit: versions

    script:
    """
    python3 ${resolver_script} \\
        --bam ${bam} \\
        --chrom-sizes ${chrom_sizes} \\
        --gca '${meta.gca_accession}' \\
        --assembly-release '${meta.assembly_release}' \\
        --entity-id '${meta.id}' \\
        --output-prefix '${meta.gca_accession}.${meta.safe_id}'
    """

    stub:
    """
    cp ${chrom_sizes} ${meta.gca_accession}.${meta.safe_id}.chrom.sizes
    printf 'gca_accession\\tassembly_release\\tentity_id\\tresolution\\tchromosome_count\\n${meta.gca_accession}\\t${meta.assembly_release}\\t${meta.id}\\tstub\\t0\\n' > ${meta.gca_accession}.${meta.safe_id}.reference_resolution.tsv
    printf '"%s":\\n    samtools: stub\\n    reference_resolution: stub\\n' '${task.process}' > ${meta.gca_accession}.${meta.safe_id}.reference_resolution.versions.yml
    """
}
