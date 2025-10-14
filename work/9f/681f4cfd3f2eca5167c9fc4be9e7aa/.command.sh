#!/bin/bash -euo pipefail
cat sample2_A.txt sample2_B.txt > sample2_combined.txt

cat <<-END_VERSIONS > versions.yml
"SUBWORKFLOW_EXAMPLE:SIMPLE_SEQUENTIAL_EXAMPLE:COMBINE_OUTPUTS":
    cat: $(cat --version | head -n1 | sed 's/cat (GNU coreutils) //g')
END_VERSIONS
