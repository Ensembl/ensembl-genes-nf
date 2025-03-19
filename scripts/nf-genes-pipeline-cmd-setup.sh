#!/usr/bin/env bash
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#### A helper script to help formulate and dump the full nextflow run command 'pipelines/nextflow/workflows/main.nf'

### Works in combination with a pre-installed ensembl-mod-env environment setup.
### Requires env conf setup: https://github.com/Ensembl/ensembl-mod-env/blob/main/confs/ensembl/busco.yml
## First setup the environment: (See: https://github.com/Ensembl/ensembl-mod-env/wiki)
# `modenv_create ensembl/busco nf-busco-pipeline`
# `module load ensembl/nf-busco-pipeline`

## Users are required to add the various options for the specific pipeline they want to setup.
## Here the default is: 'BUSCO', 'GENOME ASM STATS' and 'APPLY BETA METAKEYS' are set to run.

HOSTFLAG=$1
HOST=$2
PORTFLAG=$3
PORT=$4
USERFLAG=$5
USER=$6
CSV_INPUT=$7
OUTDIR_NAME=$8

if [[ -z $HOSTFLAG ]] || [[ -z $HOST ]] || [[ -z $PORTFLAG ]] || [[ -z $PORT ]] || [[ -z $USERFLAG ]] || [[ -z $USER ]] || [[ -z $CSV_INPUT ]] || [[ -z $OUTDIR_NAME ]]; then
echo "Usage: nf-genes-cmd-setup.sh (host[READ ONLY] details script) <INPUT CSV> <OUTDIR NAME>"
exit 0
fi

WORKDIR=`readlink -f $PWD`
mkdir -p ${WORKDIR}/$OUTDIR

# Check for modenv VAR:
if [[ ! -d $ENSEMBL_ROOT_DIR ]]; then
    echo "Unable to detect mod env variable 'ENSEMBL_ROOT_DIR'. Are you using an 'ensembl-mod-env' configurated env setup ? Exiting..."
    exit
fi

## Check mod env setup
MISSING_REPO="False"
for REPO in ensembl ensembl-analysis ensembl-genes ensembl-genes-nf ensembl-io
do
    if [[ ! -d "${ENSEMBL_ROOT_DIR}/$REPO" ]]; then
        echo "Unable to locate dependency repo: '$REPO'."
        MISSING_REPO="True"
    fi
done
if [[ "$MISSING_REPO" == "True" ]]; then
    echo "exiting...."
fi

# Check appropriate Bioperl available:
if [[ ! -d "${ENSEMBL_ROOT_DIR}/bioperl-1.6.924" ]]; then
    echo "Unable to locate BioPerl dependency 'bioperl-1.6.924'"
    ls -l $ENSEMBL_ROOT_DIR
    echo "Attempting to download and extract now..."
    wget https://github.com/bioperl/bioperl-live/archive/release-1-6-924.zip -P $ENSEMBL_ROOT_DIR
    unzip -q $ENSEMBL_ROOT_DIR/release-1-6-924.zip -d $ENSEMBL_ROOT_DIR
    mv $ENSEMBL_ROOT_DIR/bioperl-live-release-1-6-924 $ENSEMBL_ROOT_DIR/bioperl-1.6.924
    readlink -f $ENSEMBL_ROOT_DIR/bioperl-1.6.924
    BIOPERLLIB="$ENSEMBL_ROOT_DIR/bioperl-1.6.924"
else
    BIOPERLLIB="$ENSEMBL_ROOT_DIR/bioperl-1.6.924"
fi

# Some important ENV variables
ENSCODE=$ENSEMBL_ROOT_DIR
GENES_NF_REPO="${ENSEMBL_ROOT_DIR}/ensembl-genes-nf"
SERVER_SET=""
PASSWORD=""
TEAM=$TEAM_NAME

if [[ $PASSWORD == "" ]]; then
    echo "Critical... Writable host server password is not defined! Initialise setup script variable: 'PASSWORD' and rerun."
    exit 0
fi
if [[ $TEAM_NAME == "" ]]; then
    echo "Warning... Team responsible name is not defined! Suggest to define setup script variable 'TEAM' and rerun."
    exit 0
fi
if [[ $SERVER_SET == "" ]]; then
    echo "Critical... 'Server set' is not defined in this setup! Hint: Server set usually == writable host user !! Initialise setup script variable: 'SERVER_SET' and rerun."
    exit 0
fi

echo -e -n "Setup configured properly!\n\n"

echo -e -n "## Nextflow run command: ##\nnextflow run ${GENES_NF_REPO}/pipelines/nextflow/workflows/main.nf \
--bioperl $BIOPERLLIB \
--enscode $ENSCODE \
--outDir ${WORKDIR}/${OUTDIR_NAME} \
--csvFile $WORKDIR/$CSV_INPUT \
--host $HOST \
--port $PORT \
--user_r $USER \
--user_w $SERVER_SET \
--password $PASSWORD \
--server_set $SERVER_SET \
--run_busco_core true \
--apply_busco_metakeys true \
--run_ensembl_stats true \
--apply_ensembl_stats true \
--run_ensembl_beta_metakeys true \
--apply_ensembl_beta_metakeys true \
--team $TEAM \
-profile slurm\n"


# ## Pipeline parameters reference:
# Input/output parameters [REQUIRED]:
#   --enscode                     [string] env ENSCODE path.
#   --bioperl                     [string] Path to the directory containing the BioPerl library. [default: ${params.enscode}/bioperl-1.6.924]
#   --outDir                      [string] Output directory used to store finalized pipeline data.
#   --csvFile                     [string] Path for the input csv file containing the db name(s).

# General workflow parameters:
#   --run_busco_core              [boolean] Run BUSCO given a mysql db. [default: false]
#   --run_busco_ncbi              [boolean] Run BUSCO given a assembly_accession and taxonomy id in genome mode only. [default: false]
#   --run_omark                   [boolean] Run OMARK given a mysql db, default false [default: false]
#   --run_ensembl_stats           [boolean] Run Ensembl statistics given a mysql db [default: false]
#   --apply_busco_metakeys        [boolean] Create JSON file with BUSCO metakeys and load it into the db. [default: false]
#   --apply_ensembl_stats         [boolean] Insert Ensembl statistics into a mysql db [default: false]
#   --run_ensembl_beta_metakeys   [boolean] Run Ensembl beta metakeys given a mysql db. [default: false]
#   --apply_ensembl_beta_metakeys [boolean] Insert Ensembl beta metakeys into a mysql db. [default: false]
#   --host                        [string]  Full mysql host database name.
#   --port                        [integer] Four digit port number for specified mysql DB host.
#   --password                    [string]  Password for MYSQL write-enabled host connection.
#   --server_set                  [string]  Specific MYSQL write user.  (accepted: ensadmin, ensrw)
#   --team                        [string]  Required CoreDB meta_key for data team responsible.

# Optional parameters:
#   --user_w                      [string]  DB write user name.  (accepted: ensadmin, ensrw) [default: ensadmin]
#   --user_r                      [string]  DB read_only user. [default: ensro]
#   --project                     [string]  Project for the formatting of the output.  (accepted: ensembl, brc) [default: ensembl]
#   --copyToFtp                   [boolean] Copy output in Ensembl ftp. [default: false]
#   --cacheDir                    [string]  The path where downloaded files will be cached.

# Advanced parameters:
#   --canonical_only              [string] Enable 'canonical only' sequence dumps. [default: --canonical_only 1]
#   --genomio_version             [string] Ensembl-genomio Python library container. [default: latest]

# BUSCO specific config setup:
#   --busco_mode                  [array] The mode in which to run BUSCO. [default: [protein, genome]]