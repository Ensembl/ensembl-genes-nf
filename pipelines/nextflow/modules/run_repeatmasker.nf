process RUN_REPEATMASKER {
    label "python"
    tag "$gca:genome"

    // Multiple publishDir statements as needed
    publishDir "${params.outDir}/repeatmasker/", pattern: "*.fa", mode: "move"
    publishDir "${params.outDir}/repeatmasker/", pattern: "*.fa.cat", mode: "move"
    publishDir "${params.outDir}/repeatmasker/", pattern: "*.fa.masked", mode: "move"
    publishDir "${params.outDir}/repeatmasker/", pattern: "*.fa.ori.out", mode: "move"
    publishDir "${params.outDir}/repeatmasker/", pattern: "*.fa.out", mode: "move"
    publishDir "${params.outDir}/repeatmasker/", pattern: "*.fa.tbl", mode: "move"
    publishDir "${params.outDir}/repeatmasker/", pattern: "*.fa.rm.gtf", mode: "move"
    publishDir "${params.outDir}/repeatmasker/", pattern: "*.gtf", mode: "move"

    input:
    tuple val(gca), path "${params.outDir}/${gca}/ncbi_dataset/${gca}*.fna"  // Input genome file path

    output:
    tuple val(gca),
    path "*.fa", emit: slice_fasta,
    path "*.fa.cat", emit: slice_fasta_cat,
    path "*.fa.masked", emit: slice_fasta_masked,
    path "*.fa.ori.out", emit: slice_fasta_ori_out,
    path "*.fa.out", emit: slice_fa_out,
    path "*.fa.tbl", emit: slice_fa_tbl,
    path "*.fa.rm.gtf", emit: slice_fa_rm_gtf,
    path "*.gtf", emit: slice_gtf

    script:
    """
    chmod +x ${projectDir}/bin/repeatmasker.py
    chmod +x ${projectDir}/bin/_utils.py

    ${projectDir}/bin/repeatmasker.py --genome_file ${params.outDir}/${gca}/ncbi_dataset/${gca}*.fna \
                                      --output_dir ${params.outDir}/repeatmasker \
                                      --repeatmasker_bin ${params.repeatmasker_path} \
                                      --library ${params.outDir}/${gca}/rm_library/${gca}.repeatmodeler.fa \
                                      --repeatmasker_engine ${params.engine_repeatmasker}
    """
}

