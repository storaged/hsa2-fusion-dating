#!/bin/bash
# Copy all workflow scripts into the flat $HSA2_ROOT/scripts directory that the pipeline calls.
set -euo pipefail
ROOT=${HSA2_ROOT:?set HSA2_ROOT to the analysis directory (data/, results/, logs/ are created there)}
HERE=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$ROOT"/{scripts,data,results,logs}
find "$HERE/workflow" "$HERE/utils" -type f \( -name '*.py' -o -name '*.sh' -o -name '*.R' -o -name '*.slim' \) -exec cp {} "$ROOT/scripts/" \;
cp "$HERE"/config/*.json "$ROOT/scripts/"
echo "installed $(ls "$ROOT/scripts" | wc -l) scripts into $ROOT/scripts"
