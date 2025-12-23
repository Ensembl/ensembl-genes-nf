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

// run Busco in genome mode 

process BUSCO_GENOME_LINEAGE {
    label "busco"
    tag "$meta.gca:busco_genome"
    storeDir "${params.cacheDir}/$meta.gca/busco_genome/" 
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 10

    input:
    //val(busco_dataset)
    tuple val(meta), path(genome_file)
    //tuple val(gca), val(dbname), path(genome_file), val(busco_dataset), val(species_id)

    output:
    //tuple val(gca), val(dbname), path("genome_output/*.txt"), val(species_id)
    tuple val(meta), path("genome_output/*.txt"), val("${params.cacheDir}/$meta.gca/busco_genome/"), emit: busco_genome_lineage_output
    path "versions.yml", emit: versions_file

    script:
    //def buscoDataset = params.busco_dataset ? params.busco_dataset.trim() : meta.busco_dataset.trim() 

    log.info("Selected BUSCO dataset: $meta.busco_dataset")

    """
    busco -f \
    -i ${genome_file} \
    --mode genome \
    -l ${meta.busco_dataset} \
    -c ${task.cpus} \
    --out genome_output \
    --offline \
    --download_path ${params.download_path}


    cat <<EOF > versions.yml
    BUSCO_GENOME_LINEAGE:
    busco: \$(busco --version 2>&1 | head -n1 | awk '{print \$2}')
    EOF
    """
}
