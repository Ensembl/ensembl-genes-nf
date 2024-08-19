process buildSTARIndex {

    //publishDir params.genome_file.getParent(), mode: 'copy'
    label 'heavy'

    input:
        file genome

    output:
        path "star_index"

    script:
        """
        STAR --runMode genomeGenerate \
             --genomeDir star_index \
             --genomeFastaFiles $genome
        """
}

process alignReads {

    //publishDir output_dir, mode: 'copy'
    label 'heavy'

    input:
        tuple val(sample_id), path(reads)
        path star_index
        path output_dir

    output:
        tuple val(sample_id), path("${output_dir}/${sample_id}_Aligned.sortedByCoord.out.bam")

    script:
        """
        STAR --genomeDir $star_index \
             --readFilesIn ${reads[0]} ${reads[1]} \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix ${output_dir}/${sample_id}_
        """
}
