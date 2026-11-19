#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
[[ $# -ge 2 && $# -le 3 ]] || { echo 'usage: bash score/score_rosetta.sh INPUT.silent OUTPUT.sc [flow]' >&2; exit 2; }
INPUT=$(realpath "$1"); mkdir -p "$(dirname "$2")"
OUT=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
[[ -s "$INPUT" && ! -e "$OUT" ]] || { echo 'missing input or output exists' >&2; exit 2; }
# The optional legacy flow argument is accepted but does not select a protocol.
# Input must retain the normalized chain order: VHH first (A), target second (B).
XML="$ROOT/score/rosetta/score.xml"
ROSETTA="${ROSETTA_EXE:?Configure the Rosetta executable}"
DALPHABALL="${DALPHABALL_EXE:?Configure DAlphaBall}"
PSIPRED="${PSIPRED_EXE:?Configure PSIPRED}"
[[ -x "$ROSETTA" && -x "$DALPHABALL" && -x "$PSIPRED" ]] || { echo 'missing installed Rosetta/DAlphaBall/PSIPRED links' >&2; exit 2; }
WORK=$(mktemp -d "$(dirname "$OUT")/.rosetta_work.XXXXXX")
export VHH_SCORE_WORK="$WORK"
"${MPI_RUN:?Configure mpirun}" -np "${NCPU:-8}" bash "$ROOT/score/rosetta_rank.sh" \
 "$ROSETTA" @"$ROOT/score/rosetta/scoring.flags" \
 -script_vars "runpsipred_single=$PSIPRED" -dalphaball "$DALPHABALL" \
 -parser:protocol "$XML" -in:file:silent "$INPUT" -out:file:score_only "$OUT" \
 > "$OUT.log" 2>&1
