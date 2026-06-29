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


process FETCH_GENOME {
  tag "$gca:genome"
  label 'fetch_file'
  storeDir "${params.cacheDir}/$gca/ncbi_dataset/"
  afterScript "sleep $params.files_latency"  // Needed because of file system latency
  maxForks 10
  input:
  tuple val(gca), val(dbname), val(species_id), val(busco_dataset)

  output:
  tuple val(gca), val(dbname), path("*.fna"), val(busco_dataset), val(species_id)
  
script:
"""
set -euo pipefail

# 1. Try NCBI Datasets API first
if curl -fsSL \
  -H "Accept: application/zip" \
  "${params.ncbiBaseUrl}/${gca}/download?include_annotation_type=GENOME_FASTA&hydrated=FULLY_HYDRATED" \
  --output genome_file.zip \
  && unzip -j genome_file.zip '*_genomic.fna'
then
  echo "Downloaded genome FASTA from NCBI Datasets" >&2
else
  echo "Datasets API failed or zip had no genomic FASTA for ${gca}; falling back to FTP" >&2

  gca_acc="${gca}"

  prefix=\${gca_acc%%_*}
  acc_nover=\${gca_acc#*_}
  acc_digits=\${acc_nover%%.*}

  d1=\${acc_digits:0:3}
  d2=\${acc_digits:3:3}
  d3=\${acc_digits:6:3}

  ftp_parent="https://ftp.ncbi.nlm.nih.gov/genomes/all/\${prefix}/\${d1}/\${d2}/\${d3}"

  asm_dir=\$(curl -fsSL "\${ftp_parent}/" \
    | grep -oE "${gca}_[^/\\\"]+" \
    | head -n 1 || true)

  if [[ -z "\${asm_dir}" ]]; then
    echo "Could not find FTP assembly directory for ${gca} under \${ftp_parent}" >&2
    exit 1
  fi

  fasta_url="\${ftp_parent}/\${asm_dir}/\${asm_dir}_genomic.fna.gz"

  echo "Downloading \${fasta_url}" >&2
  curl -fL "\${fasta_url}" -o "\${asm_dir}_genomic.fna.gz"
  gunzip "\${asm_dir}_genomic.fna.gz"
fi

ls -lh *.fna
"""

  //ncbi_dataset/data/GCA_963576655.1/GCA_963576655.1_icGasPoly1.1_genomic.fna 

}
