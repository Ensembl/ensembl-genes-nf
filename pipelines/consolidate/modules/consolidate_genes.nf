// CONSOLIDATE_GENES
// Merge gene models from multiple annotation sources using layer priority.
// Lower priority integer = higher evidence quality.
//
// Inputs are collected GFF3 paths; priorities come from params.layer_priorities
// which maps filenames (or source labels) to integer priorities.

process CONSOLIDATE_GENES {
    label 'process_medium'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir path: "${params.outdir}/consolidate", mode: 'copy', overwrite: true,
               pattern: 'consolidated.gff3'

    input:
    // List of tuples: [path, priority]
    path gff3_files

    output:
    path "consolidated.gff3", emit: gff3
    path "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Build the --inputs args from the staged files and the priority map
    // Priority map passed as a JSON string via params.layer_priorities
    // Format: {"filename_pattern": priority}
    // Fallback: files not in map get priority 99
    def args = task.ext.args ?: ''
    """
    python3 -c "
import json, os, sys

priority_map = json.loads('''${params.layer_priorities}''')
files = '${gff3_files}'.split()
specs = []
for f in files:
    prio = 99
    for pattern, p in priority_map.items():
        if pattern in f:
            prio = int(p)
            break
    specs.append(f'{f}:{prio}')

cmd = 'consolidate_genes.py --inputs ' + ' '.join(specs) + ' --out consolidated.gff3'
os.system(cmd)
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf '##gff-version 3\\n' > consolidated.gff3
    printf 'chr1\\tconsolidation\\tgene\\t1000\\t5000\\t.\\t+\\t.\\tID=consolidated_gene_00000001;biotype=protein_coding\\n' \\
        >> consolidated.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
