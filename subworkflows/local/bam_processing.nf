include { buildSTARIndex; alignReads } from "${projectDir}/modules/local/star.nf"
include { indexBAM; generateBigWig } from "${projectDir}/modules/local/coverage.nf"

workflow bam_processing {

    take:
        species
        assembly
        parent_dir

    main:

        // Validate inputs
        if (!params.species || !params.assembly || !params.parent_dir) {
            error "Must provide species, assembly, and parent_dir parameters."
        }

        // Define file paths
        genome_file = file("${params.parent_dir}/${params.species}/${params.assembly}/genome_dumps/${params.species}_softmasked_toplevel.fa")
        reads_dir = "${params.parent_dir}/${params.species}/${params.assembly}/rnaseq/input/"
        genome_dir = genome_file.getParent()
        output_dir = "${params.parent_dir}/${params.species}/${params.assembly}/rnaseq/merged/"

        // Validate paths
        if (!genome_file.exists()) error "Genome file not found: $genome_file"
        if (!file(reads_dir).exists()) error "Reads directory not found: $reads_dir"

        // Create a channel for read pairs
        read_pairs = Channel
            .fromFilePairs("${reads_dir}/*_{1,2}.fastq.gz", checkIfExists: true)
            .ifEmpty { error "No read pairs found in ${reads_dir}" }

        // Build STAR index
        star_index = buildSTARIndex(genome_file, genome_dir)

        // Align reads
        aligned_bams = alignReads(read_pairs, star_index, output_dir)

        // Index BAM files
        indexed_bams = indexBAM(aligned_bams, output_dir)

        // Generate BigWig files
        bigwig_files = generateBigWig(indexed_bams, output_dir)
}
