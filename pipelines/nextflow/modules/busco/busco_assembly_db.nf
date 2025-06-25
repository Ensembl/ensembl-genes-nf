#!/usr/bin/env nextflow
/*
See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

process BUSCO_ASSEMBLY_DB {

    label 'python'
    tag "$gca"

    input:
    tuple val(gca), val(dbname),path(summary_file)
    
    script:
    """
    # Check if Python dependencies are installed
    # Read each line in the requirements file
    while read -r package; do \\
    if ! pip show -q "\$package" &>/dev/null; then 
        echo "\$package is not installed" 
        pip install "\$package"
    else
        echo "\$package is already installed"
    fi
    done < ${projectDir}/bin/requirements.txt

    chmod +x $projectDir/bin/busco_metakeys_to_asmdb.py
    busco_metakeys_to_asmdb.py -gca ${gca} -file ${summary_file} -asm_metadata ${params.asm_metadata}
    """  
}



