process BAMBU_DISCOVERY {
    tag "${meta.id}:${scope}:bambu"
    label 'process_high_memory'
    container "${params.bambu_container ?: 'docker://quay.io/biocontainers/bioconductor-bambu:3.12.1--r45ha27e39d_1'}"

    input:
    tuple val(meta), path(bams), path(bais), val(scope), val(data_type)
    path reference

    output:
    tuple val(meta), val(scope), path('bambu_models.gtf'), path('bambu_models.rds'), path('bambu_read_assignments.tsv'), path('bambu_input_manifest.tsv'), path('bambu_output'), emit: products
    tuple val(meta), val(scope), path('bambu_models.gtf'), emit: gtf
    path 'versions.yml', emit: versions

    script:
    def bam_files = bams instanceof List ? bams.sort { bam_file -> bam_file.name } : [bams]
    def read_vector = bam_files.collect { bam_file -> '"' + bam_file.name + '"' }.join(', ')
    def r_lib = params.bambu_r_lib ?: ''
    """
    mkdir -p bambu_home
    export HOME="\$PWD/bambu_home"
    if [ -n '${r_lib}' ]; then export R_LIBS_USER='${r_lib}'; else export R_LIBS_USER="\$PWD/bambu_home/Rlib"; fi
    mkdir -p "\$R_LIBS_USER"
    if ${params.bambu_install_biocmanager ? 'true' : 'false'} && ! Rscript -e 'quit(status=if (requireNamespace("BiocManager", quietly=TRUE)) 0 else 1)'; then
        Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org", lib=Sys.getenv("R_LIBS_USER"))'
    fi
    for bam in ${bam_files.collect { bam_file -> "'${bam_file.name}'" }.join(' ')}; do
        test -s "\$bam" || { echo "Missing or empty Bambu BAM: \$bam" >&2; exit 1; }
        test -s "\${bam}.bai" || { echo "Missing or empty Bambu BAI for: \$bam" >&2; exit 1; }
    done
    Rscript -e 'suppressPackageStartupMessages(library(bambu)); reads <- c(${read_vector}); result <- bambu(reads=reads, annotations=NULL, genome="${reference}", NDR=1, quant=${params.bambu_quantify ? 'TRUE' : 'FALSE'}, trackReads=${params.bambu_track_reads ? 'TRUE' : 'FALSE'}, ncore=${task.cpus}); dir.create("bambu_output", showWarnings=FALSE); if (${params.bambu_quantify ? 'TRUE' : 'FALSE'}) { writeBambuOutput(result, path="bambu_output"); file.copy("bambu_output/extended_annotations.gtf", "bambu_models.gtf", overwrite=TRUE) } else { writeToGTF(result, "bambu_models.gtf") }; saveRDS(result, "bambu_models.rds"); maps <- metadata(result)\$readToTranscriptMaps; if (is.null(maps) || !length(maps)) { writeLines("sample\\tread_id\\tequal_matches\\tcompatible_matches", "bambu_read_assignments.tsv") } else { rows <- do.call(rbind, lapply(seq_along(maps), function(i) { x <- as.data.frame(maps[[i]]); if (!nrow(x)) return(NULL); data.frame(sample=i, read_id=as.character(x\$readId), equal_matches=vapply(x\$equalMatches, paste, collapse=",", FUN.VALUE=character(1)), compatible_matches=vapply(x\$compatibleMatches, paste, collapse=",", FUN.VALUE=character(1)), stringsAsFactors=FALSE) })); if (is.null(rows)) writeLines("sample\\tread_id\\tequal_matches\\tcompatible_matches", "bambu_read_assignments.tsv") else write.table(rows, "bambu_read_assignments.tsv", sep="\\t", quote=FALSE, row.names=FALSE) }'
    test -s bambu_models.gtf || { echo 'Bambu did not produce candidate models' >&2; exit 1; }
    test -s bambu_read_assignments.tsv || { echo 'Bambu did not produce read assignments' >&2; exit 1; }
    {
        printf '%s\n' 'key\tvalue' \\
            'backend\tbambu' \\
            'scope\t${scope}' \\
            'accession_or_cohort\t${meta.id}' \\
            "bambu_version\t\$(Rscript -e 'cat(as.character(packageVersion("bambu")))')" \\
            'data_type\t${data_type}' \\
            'mode\tannotation_free' \\
            'quantify\t${params.bambu_quantify}' \\
            'track_reads\t${params.bambu_track_reads}' \\
            "bambu_read_assignment_count\t\$(tail -n +2 bambu_read_assignments.tsv | wc -l | tr -d ' ')" \\
            "reference_fasta_sha256\t\$(sha256sum '${reference}' | cut -d ' ' -f1)" \\
            'input_bam_names\t${bam_files.collect { bam_file -> bam_file.name }.join(',')}' \\
            "input_bam_checksums\t\$(sha256sum ${bam_files.collect { bam_file -> "'${bam_file.name}'" }.join(' ')} | tr '\\n' ';')"
    } > bambu_input_manifest.tsv
    bambu_version=\$(Rscript -e 'cat(as.character(packageVersion("bambu")))')
    printf '"%s":\\n    bambu: %s\\n' '${task.process}' "\$bambu_version" > versions.yml
    """

    stub:
    """
    printf '# bambu stub\nchrStub\tbambu\texon\t1\t4\t.\t+\t.\tgene_id "bambu_g1"; transcript_id "bambu_t1";\n' > bambu_models.gtf
    printf 'stub\n' > bambu_models.rds
    mkdir bambu_output
    printf 'stub\n' > bambu_output/extended_annotations.gtf
    printf 'sample\tread_id\tequal_matches\tcompatible_matches\n1\tread1\tbambu_t1\tbambu_t1\n' > bambu_read_assignments.tsv
    printf 'key\tvalue\nbackend\tbambu\nscope\t${scope}\naccession_or_cohort\t${meta.id}\nbambu_version\t3.12.1\ndata_type\t${data_type}\nmode\tannotation_free\nquantify\t${params.bambu_quantify}\ntrack_reads\t${params.bambu_track_reads}\nbambu_read_assignment_count\t1\n' > bambu_input_manifest.tsv
    printf '"%s":\n    bambu: 3.12.1\n' '${task.process}' > versions.yml
    """
}
