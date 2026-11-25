process BUILD_SPLICE_MATRIX {
    tag { "${meta.gca_accession}:splice_junction" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'
    publishDir "${params.outdir}/splice_junctions", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(beds), val(run_ids)
    path builder

    output:
    path "${meta.safe_id}.splice_junctions.tsv", emit: matrix
    path "${meta.safe_id}.splice_junctions.metadata.json", emit: metadata
    path "${meta.safe_id}.splice_junctions.descriptor.tsv", emit: descriptor
    path "${meta.safe_id}.splice_junctions.versions.yml", emit: versions

    script:
    def input_args = []
    beds.eachWithIndex { bed, index -> input_args << "--input '${run_ids[index]}=${bed}'" }
    """
    python3 ${builder} ${input_args.join(' ')} \
        --output ${meta.safe_id}.splice_junctions.tsv \
        --metadata-output ${meta.safe_id}.splice_junctions.metadata.json
    printf 'output_kind\toutput_file\trun_ids\nsplice_junction\t${meta.safe_id}.splice_junctions.tsv\t%s\n' '${run_ids.join(',')}' > ${meta.safe_id}.splice_junctions.descriptor.tsv
    printf '"%s":\n    build_splice_matrix: python\n' '${task.process}' > ${meta.safe_id}.splice_junctions.versions.yml
    """

    stub:
    """
    printf 'chrom\tchromStart\tchromEnd\tname\tstrand\tsplice_motif\t%s.unique_reads\t%s.multi_reads\nchr1\t0\t100\tchr1:1-100:+\t+\tGT-AG\t10\t0\n' '${run_ids[0]}' '${run_ids[0]}' > ${meta.safe_id}.splice_junctions.tsv
    printf '{"run_ids":["%s"],"row_count":1}\n' '${run_ids.join('","')}' > ${meta.safe_id}.splice_junctions.metadata.json
    printf 'output_kind\toutput_file\trun_ids\nsplice_junction\t${meta.safe_id}.splice_junctions.tsv\t%s\n' '${run_ids.join(',')}' > ${meta.safe_id}.splice_junctions.descriptor.tsv
    printf '"%s":\n    build_splice_matrix: stub\n' '${task.process}' > ${meta.safe_id}.splice_junctions.versions.yml
    """
}
