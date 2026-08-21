process PREP_AA_FASTA {
    tag { meta.id }
    label 'process_medium'
    container 'community.wave.seqera.io/library/gffread:0.12.7--33b95f1cfcc0e572'

    input:
        tuple val(meta), path(genome_fasta), path(gff3)

    output:
        tuple val(meta), path("*.clean.faa"), emit: aa_fasta
        path 'versions.yml', emit: versions

    script:
        """
        gffread ${gff3} -g ${genome_fasta} -x ${meta.id}.cds.fa -y ${meta.id}.raw.faa

        awk '
            /^>/ { print; next }
            {
                sequence = toupper(\$0)
                gsub(/[^ACDEFGHIKLMNPQRSTVWY]/, "X", sequence)
                print sequence
            }
        ' ${meta.id}.raw.faa > ${meta.id}.clean.faa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gffread: \$(gffread --version 2>&1 | sed 's/^gffread v//')
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}.faa
        touch ${meta.id}.cds.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gffread: 0.12.7
        END_VERSIONS
        """
}
