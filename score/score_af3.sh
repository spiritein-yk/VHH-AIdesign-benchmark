#!/usr/bin/env bash
set -Eeuo pipefail
[[ $# -eq 2 ]] || { echo 'usage: bash score/score_af3.sh INPUT.json OUTPUT_DIRECTORY' >&2; exit 2; }
INPUT=$(realpath "$1")
mkdir -p "$2"; OUT=$(cd "$2" && pwd)
"${SCORE_PYTHON:?Configure the AF3 interpreter}" - "$INPUT" "$OUT" <<'PY'
import json,sys,re
from pathlib import Path
x=json.load(open(sys.argv[1]))
assert x['dialect']=='alphafold3' and x['version']==1
assert x['modelSeeds']==[1], 'This protocol uses seed1'
assert re.fullmatch(r'[a-z0-9_]+',x['name']), 'Use canonical AF3 task name'
assert len(x['sequences'])==2, 'Expected VHH and target proteins'
assert all('protein' in p for p in x['sequences'])
if (Path(sys.argv[2]) / x['name']).exists():
    raise SystemExit('AF3 prediction directory already exists')
PY
export PYTHONNOUSERSITE=1
# Activate the installed AF3 3.0.0 environment before invoking this script.
"$SCORE_PYTHON" -u "${AF3_REPO:?Configure AF3}/run_alphafold.py" \
 --json_path="$INPUT" --output_dir="$OUT" \
 --model_dir="${AF3_WEIGHTS:?Configure AF3 parameters}" --db_dir="${AF3_DATABASES:?Configure the AF3 databases}" \
 --pdb_database_path="${AF3_PDB_DATABASE:?Configure the AF3 PDB database}" \
 --run_data_pipeline=true --run_inference=true

# This archived AF3 3.0.0 entrypoint has no diffusion-sample CLI option.
# Verify its native five-sample result rather than passing an unsupported flag.
"${SCORE_PYTHON:?Configure the AF3 interpreter}" - "$INPUT" "$OUT" <<'CHECK'
import csv,json,sys
from pathlib import Path
name=json.load(open(sys.argv[1]))['name']
p=Path(sys.argv[2])/name/'ranking_scores.csv'
with p.open() as f:rows=list(csv.DictReader(f))
assert len(rows)==5 and all(int(row['seed'])==1 for row in rows), 'Expected seed1 and five diffusion samples'
CHECK
