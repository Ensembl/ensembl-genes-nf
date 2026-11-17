process RUN_FUNANNOTATE {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_high'
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
    path 'funannotate.command.txt', emit: command
    path 'funannotate.stdout.log', emit: stdout
    path 'funannotate.stderr.log', emit: stderr
    path 'funannotate.status.tsv', emit: status
    path 'versions.yml', emit: versions

    script:
    def bam_args = (1..meta.bam_count).collect { index -> "--rna_bam staged/bam_${String.format('%02d', index)}.bam" }.join(' ')
    def transcript_args = meta.files.transcripts.withIndex().collect { pair -> "--transcript_evidence staged/transcript_${String.format('%02d', pair[1] + 1)}.fa" }.join(' ')
    def protein_args = meta.files.proteins.withIndex().collect { pair -> "--protein_evidence staged/protein_${String.format('%02d', pair[1] + 1)}.fa" }.join(' ')
    def pasa_arg = meta.files.pasa_gff ? '--pasa_gff staged/pasa.gff3' : ''
    def manifest_parameters = meta.files.parameters_json ? '--parameters staged/parameters.json' : ''
    def parameter_arg = parameters ? "--parameters ${parameters}" : manifest_parameters
    def isolate_arg = meta.isolate ? "--isolate '${meta.isolate.replace("'", "'\\''")}'" : ''
    def augustus_arg = augustus_species ? "--augustus_species '${augustus_species.replace("'", "'\\''")}'" : ''
    def db_arg = funannotate_db ? "--database '${funannotate_db}'" : ''
    def busco_arg = busco_db ? "--busco_db '${busco_db}'" : ''
    def command = "funannotate predict --input staged/${meta.genome_name} --out predict_results --species '${meta.species.replace("'", "'\\''")}' ${isolate_arg} ${bam_args} ${transcript_args} ${protein_args} ${pasa_arg} ${parameter_arg} ${augustus_arg} --soft_mask ${soft_mask} --min_training_models ${min_training_models} ${busco_arg} --cpus ${cpus} ${db_arg} ${extra_args}".replaceAll(/\s+/, ' ').trim()
    """
    set -euo pipefail
    ${command} > funannotate.stdout.log 2> funannotate.stderr.log
    test -s predict_results/*.gff3
    printf '%s\\tCOMPLETE\\n' '${meta.id}' > funannotate.status.tsv
    printf '%s\\n' '${command.replace("'", "'\\''")}' > funannotate.command.txt
    cat > versions.yml <<VERSIONS
"\${task.process}":
    funannotate: \$(funannotate --version 2>&1 | head -n1)
VERSIONS
    """

    stub:
    """
    mkdir -p predict_results
    printf '##gff-version 3\\n' > predict_results/${meta.safe_id}.gff3
    printf '%s\\tSTUB\\n' '${meta.id}' > funannotate.status.tsv
    printf '%s\\n' 'funannotate predict (stub)' > funannotate.command.txt
    printf 'stub\\n' > funannotate.stdout.log
    printf 'stub\\n' > funannotate.stderr.log
    printf '"%s":\\n    funannotate: stub\\n' '${task.process}' > versions.yml
    """
}
