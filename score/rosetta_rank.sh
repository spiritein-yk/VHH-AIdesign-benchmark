#!/usr/bin/env bash
set -Eeuo pipefail
RANK=${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-${SLURM_PROCID:-0}}}
RANK_DIR="${VHH_SCORE_WORK:?}/rank_$RANK"
mkdir -p "$RANK_DIR/tmp"
export TMPDIR="$RANK_DIR/tmp" TMP="$RANK_DIR/tmp" TEMP="$RANK_DIR/tmp"
cd "$RANK_DIR"
exec "$@"
