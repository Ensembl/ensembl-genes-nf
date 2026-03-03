#!/bin/bash -ue
set -euo pipefail
hold_until_arg=""
if [ -n "" ]; then hold_until_arg="--hold-until "; fi
python3 $(command -v generate_project_xml.py)       --alias prj_GCA_000001405.28_2026_02       --name "prj_GCA_000001405.28_2026_02"       --title "Annotation evidence project for GCA_000001405.28, 2026_02"       --description ""       ${hold_until_arg}       --outdir .
