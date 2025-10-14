#!/bin/bash -euo pipefail
echo "A output" > sample2_A.txt

cat <<-END_VERSIONS > versions.yml
"SUBWORKFLOW_EXAMPLE:MINIMAL_SUBWORKFLOW_EXAMPLE:TOOL_A":
    tool_a: 1.0.0
END_VERSIONS
