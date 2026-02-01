/*
 * CREATE_ANNOTATION_SQLITE
 * Create SQLite database from GTF annotation for fast lookups
 */

process CREATE_ANNOTATION_SQLITE {
    tag "${organism}_${version}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffutils:0.12--pyh7cba7a3_0' :
        'quay.io/biocontainers/gffutils:0.12--pyh7cba7a3_0' }"

    publishDir "${params.outdir}/organism_setup/${organism}/${version}", mode: 'copy'

    input:
    path(gtf)
    path(transcriptome_fasta)
    val(organism)
    val(version)

    output:
    path "*.sqlite",       emit: sqlite_db
    path "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${organism}_${version}"
    """
    #!/usr/bin/env python3

    import gffutils
    import os

    # Create database from GTF
    db_path = "${prefix}_annotation.sqlite"

    # Remove existing database if present
    if os.path.exists(db_path):
        os.remove(db_path)

    # Create the database with reasonable defaults
    db = gffutils.create_db(
        "${gtf}",
        db_path,
        force=True,
        keep_order=True,
        merge_strategy="create_unique",
        sort_attribute_values=True,
        disable_infer_genes=False,
        disable_infer_transcripts=False
    )

    print(f"Created SQLite database: {db_path}")
    print(f"Number of features: {db.count_features_of_type()}")

    # Write versions
    with open("versions.yml", "w") as f:
        f.write(f'"${task.process}":\\n')
        f.write(f'    gffutils: {gffutils.__version__}\\n')
        f.write(f'    python: 3.10\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${organism}_${version}"
    """
    touch ${prefix}_annotation.sqlite

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffutils: 0.12
        python: 3.10
    END_VERSIONS
    """
}
