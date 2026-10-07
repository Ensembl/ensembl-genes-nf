process GTF_TO_BIGBED {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_medium'
    // The production UCSC tool image can be supplied through the pipeline
    // config once the site-approved image is selected; Python is retained as
    // the portable pinned development image for stub and contract testing.
    container 'https://depot.galaxyproject.org/singularity/python:3.11'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        "${meta.gca_accession}/${meta.safe_id}/gene_model/${filename}"
    }
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(gtf), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.safe_id}.bb"), path("${meta.safe_id}.gene_model.track_result.tsv"), path("${meta.safe_id}.gene_model.provenance.tsv"), path("${meta.safe_id}.gene_model.exceptions.tsv"), emit: results
    path "${meta.safe_id}.gene_model.versions.yml", emit: versions

    script:
    def destination = "${params.outdir}/${meta.gca_accession}/${meta.safe_id}/gene_model/${meta.safe_id}.bb"
    """
    echo 'Starting GTF validation for ${meta.gca_accession}:${meta.id}' >&2
    if ! (
        awk 'NF != 9 { print "Malformed GTF row: " \$0 > "/dev/stderr"; bad=1 } END { if (bad || NR == 0) exit 1 }' ${gtf}
        gtfToGenePred -genePredExt ${gtf} models.genePred
        genePredToBed models.genePred models.bed12
        LC_ALL=C sort -k1,1 -k2,2n -k3,3n models.bed12 > sorted.bed12
        bedToBigBed -type=bed12 -tab sorted.bed12 ${chrom_sizes} ${meta.safe_id}.bb
        test -s ${meta.safe_id}.bb
    ); then
        echo 'Gene-model track failed; recording an entity-level exception' >&2
        : > ${meta.safe_id}.bb
        printf 'gca_accession\tassembly_release\tentity_id\tentity_type\ttrack_type\ttrack_path\tsource_path\tstatus\tsha256\ttool_versions\tnormalization_parameters\n%s\t%s\t%s\t%s\tgene_model\t\t%s\tfailed\t\tUCSC=runtime\tvalidation_or_conversion_failed\n' '${meta.gca_accession}' '${meta.assembly_release}' '${meta.id}' '${meta.entity_type}' '${gtf}' > ${meta.safe_id}.gene_model.track_result.tsv
        cp ${meta.safe_id}.gene_model.track_result.tsv ${meta.safe_id}.gene_model.provenance.tsv
        printf 'gca_accession\tentity_id\ttrack_type\tstage\tprocess_name\tfailure_class\n%s\t%s\tgene_model\tGTF_TO_BIGBED\t%s\tvalidation_or_conversion\n' '${meta.gca_accession}' '${meta.id}' '${task.process}' > ${meta.safe_id}.gene_model.exceptions.tsv
        printf '"%s":\n    ucsc: runtime\n' '${task.process}' > ${meta.safe_id}.gene_model.versions.yml
        exit 0
    fi
    sha=\$(sha256sum ${meta.safe_id}.bb | awk '{print \$1}')
    printf 'gca_accession\tassembly_release\tentity_id\tentity_type\ttrack_type\ttrack_path\tsource_path\tstatus\tsha256\ttool_versions\tnormalization_parameters\n' > ${meta.safe_id}.gene_model.track_result.tsv
    printf '%s\t%s\t%s\t%s\tgene_model\t%s\t%s\tcomplete\t%s\tUCSC=runtime\tsort=LC_ALL=C\n' '${meta.gca_accession}' '${meta.assembly_release}' '${meta.id}' '${meta.entity_type}' '${destination}' '${gtf}' "\$sha" >> ${meta.safe_id}.gene_model.track_result.tsv
    cp ${meta.safe_id}.gene_model.track_result.tsv ${meta.safe_id}.gene_model.provenance.tsv
    : > ${meta.safe_id}.gene_model.exceptions.tsv
    printf '"%s":\n    ucsc: runtime\n' '${task.process}' > ${meta.safe_id}.gene_model.versions.yml
    """

    stub:
    """
    printf 'stub\n' > ${meta.safe_id}.bb
    printf 'gca_accession\tassembly_release\tentity_id\tentity_type\ttrack_type\ttrack_path\tsource_path\tstatus\tsha256\ttool_versions\tnormalization_parameters\n${meta.gca_accession}\t${meta.assembly_release}\t${meta.id}\t${meta.entity_type}\tgene_model\t${params.outdir}/${meta.gca_accession}/${meta.safe_id}/gene_model/${meta.safe_id}.bb\t${gtf}\tcomplete\tstub\tUCSC=stub\tsort=LC_ALL=C\n' > ${meta.safe_id}.gene_model.track_result.tsv
    cp ${meta.safe_id}.gene_model.track_result.tsv ${meta.safe_id}.gene_model.provenance.tsv
    : > ${meta.safe_id}.gene_model.exceptions.tsv
    printf '"%s":\n    ucsc: stub\n' '${task.process}' > ${meta.safe_id}.gene_model.versions.yml
    """
}
