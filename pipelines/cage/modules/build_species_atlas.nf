process BUILD_SPECIES_ATLAS {
    tag "${meta.species}:${meta.assembly}"
    label 'cage_r'
    publishDir { "${params.outdir}/atlas/${meta.species}/${meta.assembly}" }, mode: 'copy'

    input:
    tuple val(meta), val(sample_records), path(ctss_bed), path(ctss_bedgraph), path(gtf)
    path helper

    output:
    tuple val(meta), path("atlas_${meta.assembly}"), emit: species_atlas

    script:
    def beds = ctss_bed instanceof List ? ctss_bed : [ctss_bed]
    def bedgraphs = ctss_bedgraph instanceof List ? ctss_bedgraph : [ctss_bedgraph]
    def sample_rows = sample_records.collect { record -> [record.sample_id, record.species, record.assembly, record.assay_class, record.protocol, record.layout, record.tissue, record.condition, record.biological_replicate].join('\t') }.join('\n')
    def bed_args = beds.join(' ')
    def bedgraph_copy = bedgraphs.collect { bedgraph -> "cp ${bedgraph} atlas_${meta.assembly}/${bedgraph.baseName}" }.join('\n')
    """
    mkdir atlas_${meta.assembly}
    cat > sample_metadata.tsv <<'EOF_SAMPLE_METADATA'
    sample_id\tspecies\tassembly\tassay_class\tprotocol\tlayout\ttissue\tcondition\tbiological_replicate
    ${sample_rows}
    EOF_SAMPLE_METADATA
    ${bedgraph_copy}
    Rscript ${helper} --assembly ${meta.assembly} --species ${meta.species} --metadata sample_metadata.tsv --beds ${bed_args} --gtf ${gtf} --output atlas_${meta.assembly}
    """

    stub:
    """
    mkdir atlas_${meta.assembly}
    touch atlas_${meta.assembly}/sample_metadata.tsv atlas_${meta.assembly}/tss_clusters.bed atlas_${meta.assembly}/enhancers.bed atlas_${meta.assembly}/tss_annotations.tsv atlas_${meta.assembly}/ctss_counts.tsv atlas_${meta.assembly}/tss_counts.tsv atlas_${meta.assembly}/promoter_counts.tsv atlas_${meta.assembly}/atlas.qc.tsv
    """
}
