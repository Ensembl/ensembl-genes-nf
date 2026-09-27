process IRIBO_GENERATE_TRANSLATOME {
    tag "${meta.id}:${meta.shard_id ?: 'all'}"
    label 'process_ultra_high'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-iribo:1.0.0'
    input:
    tuple val(meta), path(profile), path(candidates)
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    script:
    """
    mkdir -p raw
    set +e
    Rscript \$IRIBO_HOME/GenerateTranslatome.R --TranslationCalls=${profile}/translation_calls --NullDistribution=${profile}/null_distribution --CandidateORFs=${candidates}/candidate_orfs --Output=raw/translatome --Threads=${task.cpus ?: 2} --Scrambles=1 > raw/translatome.log 2>&1
    rc=\$?
    set -e
    if [ \$rc -ne 0 ]; then
        if grep -Eqi 'No reads detected|empty scrambled|empty statistic' raw/translatome.log; then
            # Preserve a genuine native zero-call result; do not invent ORFs.
            printf 'id,chrom,start,end,strand,start_codon,score\n' > raw/translatome/translated_orfs.csv
            printf 'iRibo completed with zero calls: upstream statistics had no valid reads\n' > raw/translatome/zero_call.txt
        else
            cat raw/translatome.log >&2
            exit \$rc
        fi
    fi
    test -s raw/translatome/translated_orfs.csv || { echo 'iRibo produced no translatome output' >&2; exit 1; }
    printf '"%s":\n    iRibo: source-pinned\n    stage: translatome\n' '${task.process}' > versions.yml
    """
    stub:
    """
    mkdir -p raw/translatome
    printf 'id,chrom,start,end,strand,start_codon,score\nIR1,chr1,100,200,+,ATG,0.9\n' > raw/translatome/translated_orfs.csv
    printf '"stub":\n    iRibo: stub\n' > versions.yml
    """
}
