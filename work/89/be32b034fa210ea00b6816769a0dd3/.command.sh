#!/bin/bash -euo pipefail
wc -l sample1_combined.txt > sample1_line_count.txt

cat <<-END_VERSIONS > versions.yml
"SUBWORKFLOW_EXAMPLE:SIMPLE_SEQUENTIAL_EXAMPLE:COUNT_LINES":
    wc: $(wc --version | head -n1 | sed 's/wc (GNU coreutils) //g')
END_VERSIONS
