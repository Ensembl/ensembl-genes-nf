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
/*This process uploads the RepeatModeler library files (FASTA, STK, and log files)
to the Ensembl FTP server for a given species. It creates a directory 
for the species on the FTP server and copies the files into that directory. 
The process also generates a versions.yml file containing the version 
of the upload process used.*/
process UPLOAD_REPEATS_INTO_FTP {
    tag "$meta.gca:upload_repeats_into_ftp"
        label 'ensembl_ftp'

input:
tuple val(meta),path(fasta_file),path(stk_file),path(log_file)

output:
tuple val(meta),path(fasta_file), emit: library_out                                                            
path "versions.yml", emit: versions_file

script:
"""
sudo -u genebuild rsync -ahvW #output_path#/ftp_release/ /nfs/ftp/public/databases/ensembl/pre-release
sudo -u genebuild mkdir -p /nfs/ftp/public/databases/ensembl/pre-release/repeats/species/${meta.species_name}
sudo -u genebuild cp ${fasta_file} /nfs/ftp/public/databases/ensembl/pre-release/repeats/species/${meta.species_name}
sudo -u genebuild cp ${stk_file} /nfs/ftp/public/databases/ensembl/pre-release/repeats/species/${meta.species_name}
sudo -u genebuild cp ${log_file} /nfs/ftp/public/databases/ensembl/pre-release/repeats/species/${meta.species_name}
cat <<'EOF' > versions.yml
    "${task.process}":
        upload_into_ftp: "1.0"
            EOF
"""

}
