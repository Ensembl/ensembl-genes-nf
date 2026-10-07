process PREPARE_SPLICE_JUNCTIONS {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_1'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        "${meta.gca_accession}/${meta.safe_id}/splice_junction/${filename}"
    }
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(sj_out_tab, stageAs: 'input.sj'), path(bam, stageAs: 'input.bam'), path(chrom_sizes, stageAs: 'chrom.sizes')

    output:
    tuple val(meta), path("${meta.safe_id}.bb"), path("${meta.safe_id}.splice_junction.track_result.tsv"), path("${meta.safe_id}.splice_junction.provenance.tsv"), path("${meta.safe_id}.splice_junction.exceptions.tsv"), emit: results
    path "${meta.safe_id}.splice_junction.versions.yml", emit: versions

    script:
    def derive = meta.derive_junctions
    def source = derive ? bam : sj_out_tab
    def destination = "${params.outdir}/${meta.gca_accession}/${meta.safe_id}/splice_junction/${meta.safe_id}.bb"
    """
    echo 'Starting splice-junction validation for ${meta.gca_accession}:${meta.id}' >&2
    if ! (
        ${derive ? "regtools junctions extract -o derived.sj.out.tab ${bam}" : "awk 'NF != 9 || \\$2 < 1 || \\$3 < \\$2 { print \"Malformed SJ row: \" \\$0 > \"/dev/stderr\"; bad=1 } END { if (bad || NR == 0) exit 1 }' ${sj_out_tab}"}
    ); then
        echo 'Splice-junction input failed validation or derivation; recording an exception' >&2
        : > ${meta.safe_id}.bb
        printf 'gca_accession\tassembly_release\tentity_id\tentity_type\ttrack_type\ttrack_path\tsource_path\tstatus\tsha256\ttool_versions\tnormalization_parameters\n%s\t%s\t%s\t%s\tsplice_junction\t\t%s\tfailed\t\tUCSC=runtime\tinput_validation_or_derivation_failed\n' '${meta.gca_accession}' '${meta.assembly_release}' '${meta.id}' '${meta.entity_type}' '${source}' > ${meta.safe_id}.splice_junction.track_result.tsv
        cp ${meta.safe_id}.splice_junction.track_result.tsv ${meta.safe_id}.splice_junction.provenance.tsv
        printf 'gca_accession\tentity_id\ttrack_type\tstage\tprocess_name\tfailure_class\n%s\t%s\tsplice_junction\tPREPARE_SPLICE_JUNCTIONS\t%s\tinput_validation_or_derivation\n' '${meta.gca_accession}' '${meta.id}' '${task.process}' > ${meta.safe_id}.splice_junction.exceptions.tsv
        printf '"%s":\n    ucsc: runtime\n' '${task.process}' > ${meta.safe_id}.splice_junction.versions.yml
        exit 0
    fi
    input_sj=${derive ? 'derived.sj.out.tab' : sj_out_tab}
    awk 'BEGIN { OFS="\\t" } { strand=(\$4 == 1 ? "+" : (\$4 == 2 ? "-" : ".")); print \$1, \$2-1, \$3, \$1":"\$2"-"\$3, \$7, strand }' \"\$input_sj\" | LC_ALL=C sort -k1,1 -k2,2n -k3,3n > junctions.bed
    bedToBigBed -type=bed6+1 -tab junctions.bed ${chrom_sizes} ${meta.safe_id}.bb
    sha=\$(sha256sum ${meta.safe_id}.bb | awk '{print \$1}')
    printf 'gca_accession\tassembly_release\tentity_id\tentity_type\ttrack_type\ttrack_path\tsource_path\tstatus\tsha256\ttool_versions\tnormalization_parameters\n' > ${meta.safe_id}.splice_junction.track_result.tsv
    printf '%s\t%s\t%s\t%s\tsplice_junction\t%s\t%s\tcomplete\t%s\tUCSC=runtime\tSTAR_start_to_BED_start=start-1;end=end\n' '${meta.gca_accession}' '${meta.assembly_release}' '${meta.id}' '${meta.entity_type}' '${destination}' '${source}' "\$sha" >> ${meta.safe_id}.splice_junction.track_result.tsv
    cp ${meta.safe_id}.splice_junction.track_result.tsv ${meta.safe_id}.splice_junction.provenance.tsv
    : > ${meta.safe_id}.splice_junction.exceptions.tsv
    printf '"%s":\n    ucsc: runtime\n' '${task.process}' > ${meta.safe_id}.splice_junction.versions.yml
    """

    stub:
    def stub_source = meta.derive_junctions ? bam : sj_out_tab
    """
    printf 'stub\n' > ${meta.safe_id}.bb
    printf 'gca_accession\tassembly_release\tentity_id\tentity_type\ttrack_type\ttrack_path\tsource_path\tstatus\tsha256\ttool_versions\tnormalization_parameters\n${meta.gca_accession}\t${meta.assembly_release}\t${meta.id}\t${meta.entity_type}\tsplice_junction\t${params.outdir}/${meta.gca_accession}/${meta.safe_id}/splice_junction/${meta.safe_id}.bb\t${stub_source}\tcomplete\tstub\tUCSC=stub\tSTAR_start_to_BED_start=start-1;end=end\n' > ${meta.safe_id}.splice_junction.track_result.tsv
    cp ${meta.safe_id}.splice_junction.track_result.tsv ${meta.safe_id}.splice_junction.provenance.tsv
    : > ${meta.safe_id}.splice_junction.exceptions.tsv
    printf '"%s":\n    ucsc: stub\n' '${task.process}' > ${meta.safe_id}.splice_junction.versions.yml
    """
}
