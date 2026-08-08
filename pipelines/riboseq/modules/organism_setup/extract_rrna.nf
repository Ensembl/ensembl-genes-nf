/*
 * EXTRACT_RRNA
 * Extract rRNA sequences from GTF annotation
 */

process EXTRACT_RRNA {
    tag "${organism}_${version}"
    label 'process_low'

    container "quay.io/biocontainers/gffread:0.9.12--0"

    input:
    path(gtf)
    path(genome_fasta)
    val(organism)
    val(version)

    output:
    path "rRNA.fa",        emit: rrna_fasta, optional: true
    path "versions.yml",   emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Extract rRNA annotations from GTF
    # Common rRNA types: rRNA, Mt_rRNA, rRNA_pseudogene, 5S_rRNA, etc.
    grep -E 'rRNA|rRNA_pseudogene|Mt_rRNA|5S_rRNA|5_8S_rRNA|12S_rRNA|16S_rRNA|28S_rRNA|18S_rRNA' ${gtf} > rRNA.gtf || true

    # Also try biotype-based extraction
    grep -E 'gene_biotype "rRNA"|gene_type "rRNA"|transcript_biotype "rRNA"' ${gtf} >> rRNA.gtf || true
    grep -E 'gene_biotype "Mt_rRNA"|gene_type "Mt_rRNA"' ${gtf} >> rRNA.gtf || true
    grep -E 'gene_biotype "rRNA_pseudogene"|gene_type "rRNA_pseudogene"' ${gtf} >> rRNA.gtf || true

    # Remove duplicates
    sort -u rRNA.gtf > rRNA_unique.gtf
    mv rRNA_unique.gtf rRNA.gtf

    # Check if any rRNA annotations were found
    if [ ! -s rRNA.gtf ]; then
        echo "WARNING: No rRNA annotations found in GTF."
        touch rRNA.fa  # Create empty file for optional output
    else
        echo "Found \$(wc -l < rRNA.gtf) rRNA annotation lines"

        # Extract rRNA sequences using gffread
        gffread -w rRNA.fa -g ${genome_fasta} rRNA.gtf ${args}

        # Check if sequences were extracted
        if [ ! -s rRNA.fa ]; then
            echo "WARNING: No rRNA sequences could be extracted from GTF."
        else
            echo "Successfully extracted \$(grep -c '^>' rRNA.fa) rRNA sequences"
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch rRNA.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: 0.12.7
    END_VERSIONS
    """
}
