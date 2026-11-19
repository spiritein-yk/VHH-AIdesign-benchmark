#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
[[ $# -eq 2 ]] || { echo 'usage: bash score/score_af2.sh INPUT.silent OUTPUT_DIRECTORY' >&2; exit 2; }
INPUT=$(realpath "$1")
mkdir -p "$2"; OUT=$(cd "$2" && pwd)
[[ -s "$INPUT" ]] || { echo 'missing input silent' >&2; exit 2; }
[[ ! -e "$OUT/out.sc" && ! -e "$OUT/out.silent" && ! -e "$OUT/check.point" ]] || { echo 'output directory already contains AF2 results' >&2; exit 2; }
export PYTHONNOUSERSITE=1
export PYTHONPATH="${AF2_REPO:?Configure the initial-guess fork}:${AF2_SILENT_TOOLS:?Configure silent_tools}${PYTHONPATH:+:$PYTHONPATH}"
export AF2_DATA_DIR="${AF2_DATA_DIR:?Configure AF2 parameters}"
cd "$OUT"
"${SCORE_PYTHON:?Configure the AF2 interpreter}" -u "$ROOT/score/af2_predictor.py" -silent "$INPUT" -recycle 3 > model.log 2>&1
