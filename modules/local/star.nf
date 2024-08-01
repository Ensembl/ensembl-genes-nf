process buildSTARIndex {

    publishDir genome_dir, mode: 'copy'

    input:
        path genome
        path genome_dir

    output:
        path "star_index"

    script:
        """
        mkdir star_index
        STAR --runMode genomeGenerate \
             --genomeDir ${genome_dir}/star_index \
             --genomeFastaFiles $genome
        """
}

process alignReads {

    publishDir output_dir, mode: 'copy'

    input:
        tuple val(sample_id), path(reads)
        path star_index
        path output_dir

    output:
        tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")

    script:
        """
        STAR --genomeDir $star_index \
             --readFilesIn ${reads[0]} ${reads[1]} \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix ${sample_id}_
        """
}
