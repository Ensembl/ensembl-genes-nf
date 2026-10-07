process BAMCOVERAGE {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/deeptools:3.5.5--pyhdfd78af_0'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        "${meta.gca_accession}/${meta.safe_id}/coverage/${filename}"
    }
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(meta), path(bam), path(bai), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.safe_id}.bw"), path("${meta.safe_id}.coverage.track_result.tsv"), path("${meta.safe_id}.coverage.provenance.tsv"), path("${meta.safe_id}.coverage.exceptions.tsv"), emit: results
    path "${meta.safe_id}.coverage.versions.yml", emit: versions

    script:
    def bin_size = task.ext.bin_size ?: params.coverage_bin_size
    def normalization = task.ext.normalization ?: params.coverage_normalization
    def destination = "${params.outdir}/${meta.gca_accession}/${meta.safe_id}/coverage/${meta.safe_id}.bw"
    """
    echo 'Starting BAM validation for ${meta.gca_accession}:${meta.id}' >&2
    samtools quickcheck -v ${bam}
    samtools view -H ${bam} | grep -q 'SO:coordinate' || { echo 'BAM is not coordinate sorted' >&2; exit 1; }
    bamCoverage --bam ${bam} --outFileName ${meta.safe_id}.bw --binSize ${bin_size} --normalizeUsing ${normalization} --numberOfProcessors ${task.cpus} --effectiveGenomeSize 0
    test -s ${meta.safe_id}.bw || { echo 'Coverage BigWig is empty' >&2; exit 1; }
    sha=\$(sha256sum ${meta.safe_id}.bw | awk '{print \$1}')
    printf 'gca_accession\tassembly_release\tentity_id\tentity_type\ttrack_type\ttrack_path\tsource_path\tstatus\tsha256\ttool_versions\tnormalization_parameters\n' > ${meta.safe_id}.coverage.track_result.tsv
    printf '%s\t%s\t%s\t%s\tcoverage\t%s\t%s\tcomplete\t%s\tdeepTools=%s\tbin_size=%s;normalization=%s\n' '${meta.gca_accession}' '${meta.assembly_release}' '${meta.id}' '${meta.entity_type}' '${destination}' '${bam}' "\$sha" "\$(bamCoverage --version 2>&1 | head -n1)" '${bin_size}' '${normalization}' >> ${meta.safe_id}.coverage.track_result.tsv
    cp ${meta.safe_id}.coverage.track_result.tsv ${meta.safe_id}.coverage.provenance.tsv
    : > ${meta.safe_id}.coverage.exceptions.tsv
    cat <<-END_VERSIONS > ${meta.safe_id}.coverage.versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1)
        deepTools: \$(bamCoverage --version 2>&1 | head -n1)
    END_VERSIONS
    """

    stub:
    def stub_bin_size = task.ext.bin_size ?: params.coverage_bin_size
    def stub_normalization = task.ext.normalization ?: params.coverage_normalization
    """
    printf 'stub\n' > ${meta.safe_id}.bw
    printf 'gca_accession\tassembly_release\tentity_id\tentity_type\ttrack_type\ttrack_path\tsource_path\tstatus\tsha256\ttool_versions\tnormalization_parameters\n${meta.gca_accession}\t${meta.assembly_release}\t${meta.id}\t${meta.entity_type}\tcoverage\t${params.outdir}/${meta.gca_accession}/${meta.safe_id}/coverage/${meta.safe_id}.bw\t${bam}\tcomplete\tstub\tdeepTools=stub\tbin_size=${stub_bin_size};normalization=${stub_normalization}\n' > ${meta.safe_id}.coverage.track_result.tsv
    cp ${meta.safe_id}.coverage.track_result.tsv ${meta.safe_id}.coverage.provenance.tsv
    : > ${meta.safe_id}.coverage.exceptions.tsv
    printf '"%s":\n    samtools: stub\n    deepTools: stub\n' '${task.process}' > ${meta.safe_id}.coverage.versions.yml
    """
}
