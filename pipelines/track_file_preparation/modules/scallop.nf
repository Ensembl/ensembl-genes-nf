process SCALLOP {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_high_memory'
    container 'https://depot.galaxyproject.org/singularity/scallop:0.10.5--hea69786_9'

    input:
    tuple val(meta), path(bam), path(bai), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.safe_id}.scallop.gtf"), path(chrom_sizes), emit: gtf
    path "${meta.safe_id}.scallop.versions.yml", emit: versions

    script:
    def library_type = params.scallop_library_type ?: 'auto'
    def library_arg = library_type == 'auto' ? '' : "--library_type ${library_type}"
    def extra_args = params.scallop_args ?: ''
    """
    set -euo pipefail
    test -s ${bam} || { echo 'Scallop input BAM is empty' >&2; exit 1; }
    scallop -i ${bam} -o ${meta.safe_id}.scallop.gtf ${library_arg} ${extra_args}
    test -s ${meta.safe_id}.scallop.gtf || { echo 'Scallop produced an empty GTF' >&2; exit 1; }
    printf '"%s":\\n    scallop: %s\\n' '${task.process}' "\$(scallop --version 2>&1 | head -n1)" > ${meta.safe_id}.scallop.versions.yml
    """

    stub:
    """
    cat > ${meta.safe_id}.scallop.gtf <<'EOF'
    chr1\tscallop\ttranscript\t1\t10\t.\t+\t.\tgene_id "stub_gene"; transcript_id "stub_tx";
    chr1\tscallop\texon\t1\t10\t.\t+\t.\tgene_id "stub_gene"; transcript_id "stub_tx"; exon_number "1";
    EOF
    printf '"%s":\\n    scallop: stub\\n' '${task.process}' > ${meta.safe_id}.scallop.versions.yml
    """
}
