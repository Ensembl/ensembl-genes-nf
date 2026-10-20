#!/usr/bin/env nextflow
/*
This process uploads the RepeatModeler library files (FASTA, STK, and log files)
to the Ensembl FTP server for a given species. It creates a directory 
for the species on the FTP server and copies the files into that directory. 
The process also generates a versions.yml file containing the version 
of the upload process used.
*/
process UPLOAD_REPEAT_LIBRARY_INTO_FTP {
    tag "$meta.gca:upload_repeat_library_into_ftp"
        label 'ensembl_ftp'

input:
tuple val(meta),path(fasta_file),path(stk_file),path(log_file)
val ftp_base

output:
tuple val(meta),path(fasta_file), emit: library_out                                                            
path "versions.yml", emit: versions_file

script:
"""
sudo -u genebuild mkdir -p ${ftp_base}/species/${meta.species_name}
sudo -u genebuild cp ${fasta_file} ${ftp_base}/species/${meta.species_name}
sudo -u genebuild cp ${stk_file} ${ftp_base}/species/${meta.species_name}
sudo -u genebuild cp ${log_file} ${ftp_base}/species/${meta.species_name}
cat <<'EOF' > versions.yml
    "${task.process}":
    upload_repeat_library_into_ftp: "1.0"
EOF
"""

}
