process AUDIT_FASTQ {
    tag "${meta.id}:${meta.classification}"
    label 'process_light'

    input:
    tuple val(meta), path(probe), path(validation)
    path auditor

    output:
    tuple val(meta), path('molecule_audit.tsv'), emit: molecule_audit
    path 'versions.yml', emit: versions

    script:
    """
    ./${auditor} ${probe} ${validation} '${meta.id}' '${meta.classification ?: 'UNCLASSIFIED'}' molecule_audit.tsv
    printf '"%s":\n    audit: generated\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'run_accession\tclassification\theader_representation\trecords_sampled\tdistinct_ids\tdistinct_molecules\tmalformed\tstatus\n${meta.id}\t${meta.classification ?: 'UNCLASSIFIED'}\tUNKNOWN\t1\t1\t1\t0\tstub\n' > molecule_audit.tsv
    printf '"%s":\n    audit: generated\n' '${task.process}' > versions.yml
    """
}
