process EXAMPLE_MODULE {
    // Resource label for appropriate compute allocation
    label 'process_medium'
    
    // Software dependencies - currently only support singularity as a requirement d
    container "https://depot.galaxyproject.org/singularity/mulled-v2-example:latest"
    
    // Dynamic tagging for better process identification
    tag "${meta.id}"

    errorStrategy 'ignore'
    
    // Publish results to organized output directory
    publishDir "${params.outdir}/example_module", 
        pattern: "*.{txt,log,html}",
        saveAs: { filename -> 
            if (filename.endsWith('.log')) "logs/${meta.id}_${filename}"
            else if (filename.endsWith('.html')) "reports/${meta.id}_${filename}" 
            else filename 
        }
    
    input:
    tuple val(meta), path(reads)            // Meta map with sample info + input files
    path reference                          // Reference file (shared across samples)
    val mode                                // Required parameter for tool execution
    
    output:
    tuple val(meta), path("${prefix}.txt"),             emit: results
    tuple val(meta), path("${prefix}.log"),             emit: logs
    tuple val(meta), path("${prefix}_report.html"),     emit: reports, optional: true
    tuple val(meta), path("${prefix}_stats.json"),      emit: stats,   optional: true
    path "versions.yml",                                emit: versions
    
    when:
    // Use ext.when instead of hardcoded when conditions
    // this enables the module to be conditionally executed based on external parameters
    // and these parameters come from config rather than $params.<param> meaning they are 
    // more flexible for reuse
    task.ext.when == null || task.ext.when
    
    script:
    // Define variables for flexibility and reusability
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_paired = meta.single_end ? false : true
    def input_files = is_paired ? "${reads[0]} ${reads[1]}" : "${reads}"
    
    // Check for input/output name conflicts
    if ("${reads}" == "${prefix}.txt") {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }
    
    """
    # Main tool execution with proper argument handling
    example_tool \\
        --input ${input_files} \\
        --reference ${reference} \\
        --mode ${mode} \\
        --threads ${task.cpus} \\
        --output ${prefix}.txt \\
        ${args}
    
    # Secondary tool for additional processing (if needed)
    if [[ -s ${prefix}.txt ]]; then
        secondary_tool \\
            --input ${prefix}.txt \\
            --output ${prefix}_stats.json \\
            ${args2} \\
            || echo "Secondary tool failed, continuing..." >&2
    fi
    
    # Generate report if conditions are met
    if [[ "${mode}" == "detailed" ]]; then
        generate_report \\
            --input ${prefix}.txt \\
            --output ${prefix}_report.html \\
            --sample-id ${meta.id}
    fi
    
    # Capture versions for reproducibility
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        example_tool: \\$(example_tool --version | sed 's/example_tool v//g')
        secondary_tool: \\$(secondary_tool --version | grep -oP '\\d+\\.\\d+\\.\\d+')
        samtools: \\$(samtools --version | sed '1!d ; s/samtools //')
    END_VERSIONS
    """
    
    stub:
    // Stub block for testing - must create all expected outputs
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Create main output files
    touch ${prefix}.txt
    touch ${prefix}.log
    
    # Create optional outputs based on conditions
    if [[ "${mode}" == "detailed" ]]; then
        touch ${prefix}_report.html
    fi
    
    # Create compressed stub files properly
    echo "stub stats data" | gzip > ${prefix}_stats.json.gz || touch ${prefix}_stats.json
    
    # Generate versions file (same as main script)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        example_tool: 0.11.9
        secondary_tool: 1.2.3
        samtools: 1.17
    END_VERSIONS
    """
}