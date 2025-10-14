// TODO: Rename this process to match your tool/module name (e.g., MY_TOOL)
process EXAMPLE_MODULE {
    // TODO: Update the process label based on resource requirements (process_low, process_medium, process_high)
    label 'process_medium'

    // TODO: Update container URL to point to your tool's Singularity/Docker container
    container "https://depot.galaxyproject.org/singularity/mulled-v2-example:latest"
    
    tag "${meta.id}"

    errorStrategy 'ignore'

    // TODO: Update publishDir path to reflect your module's output directory name
    publishDir "${params.outdir}/example_module",
        saveAs: 'example_module/${meta.id}_${filename}'

    // TODO: Update input channels to match your tool's requirements
    input:
    tuple val(meta), path(input_file)

    // TODO: Update output files and emit names to match what your tool produces
    output:
    tuple val(meta), path("${prefix}.txt"),             emit: results 
    path "versions.yml",                                emit: versions 
    
    when:
    task.ext.when == null || task.ext.when

    // TODO: Implement your tool's actual command-line execution here
    script:
    """
    # Main tool execution with proper argument handling
    # TODO: Replace 'example_tool' with your actual tool name and update all arguments
    example_tool \\
        --input ${input_files} \\
        --reference ${reference} \\
        --mode ${mode} \\
        --threads ${task.cpus} \\
        --output ${prefix}.txt \\
        ${args}

    # Capture versions for reproducibility
    # TODO: Update version extraction command to match your tool's --version output format
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        example_tool: \\$(example_tool --version | sed 's/example_tool v//g')
    END_VERSIONS
    """
    
    // TODO: Update stub section to create mock output files matching your tool's actual outputs
    stub:
    """
    touch ${prefix}.txt
    touch ${prefix}.log

    # TODO: Update tool name and version number in stub
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        example_tool: 0.11.9
    END_VERSIONS
    """
}