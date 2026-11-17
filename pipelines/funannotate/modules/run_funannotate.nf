def shell_quote(String value) {
    "'${value.replace("'", "'\\''")}'"
}

process RUN_FUNANNOTATE {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_high'
    container params.funannotate_container ?: '/nfs/production/flicek/ensembl/shared_data/singularity_images/cache/funannotate_latest.sif'
    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${meta.assembly_accession}/${meta.safe_id}/${filename}" }
    input:
    tuple val(meta), path(staged)
    val funannotate_db
    val cpus
    val augustus_species
    val soft_mask
    val min_training_models
    val busco_db
    val parameters
    val extra_args

    output:
    tuple val(meta), path('predict_results'), emit: results
    path 'predict_results/*.gff3', emit: gff3
    path 'funannotate.stdout.log', emit: stdout
    path 'funannotate.stderr.log', emit: stderr
    path 'versions.yml', emit: versions

    script:
    def args = [
        '--input', "staged/${meta.genome_name}",
        '--out', 'predict_results',
        '--species', shell_quote(meta.species),
        meta.isolate ? "--isolate ${shell_quote(meta.isolate)}" : null,
        '--rna_bam staged/all_runs.merged.sorted.bam',
        meta.files.transcripts.withIndex().collect { pair -> "--transcript_evidence staged/transcript_${String.format('%02d', pair[1] + 1)}.fa" },
        meta.files.stringtie ? '--transcript_evidence staged/stringtie_transcripts.fa' : null,
        meta.files.proteins.withIndex().collect { pair -> "--protein_evidence staged/protein_${String.format('%02d', pair[1] + 1)}.fa" },
        meta.files.pasa_gff ? '--pasa_gff staged/pasa.gff3' : null,
        meta.files.parameters_json && !parameters ? '--parameters staged/parameters.json' : null,
        parameters ? "--parameters ${shell_quote(parameters)}" : null,
        meta.files.stringtie ? '--stringtie staged/stringtie.gtf' : null,
        meta.files.protein_alignments ? '--protein_alignments staged/genblast_protein_alignments.gff3' : null,
        meta.files.trnascan ? '--trnascan staged/trnascan.gff' : null,
        augustus_species ? "--augustus_species ${shell_quote(augustus_species)}" : null,
        '--soft_mask', soft_mask,
        '--min_training_models', min_training_models,
        busco_db ? "--busco_db ${shell_quote(busco_db)}" : null,
        '--cpus', cpus,
        funannotate_db ? "--database ${shell_quote(funannotate_db)}" : null,
        extra_args ?: null
    ].flatten().findAll { it != null && it.toString() }.collect { it.toString() }.join(' ')
    """
    funannotate predict ${args} > funannotate.stdout.log 2> funannotate.stderr.log
    cat > versions.yml <<VERSIONS
"\${task.process}":
    funannotate: \$(funannotate --version 2>&1 | head -n1)
VERSIONS
    """

    stub:
    """
    mkdir -p predict_results
    touch predict_results/${meta.safe_id}.gff3 funannotate.stdout.log funannotate.stderr.log versions.yml
    """
}
