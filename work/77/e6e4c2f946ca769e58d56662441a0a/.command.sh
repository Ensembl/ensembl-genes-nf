#!/bin/bash -euo pipefail
echo "B output" > sample1_B.txt

cat <<-END_VERSIONS > versions.yml
"SUBWORKFLOW_EXAMPLE:MINIMAL_SUBWORKFLOW_EXAMPLE:TOOL_B":
    tool_b: 1.0.0
END_VERSIONS
