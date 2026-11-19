#!/usr/bin/env python3
"""Run one configured benchmark scoring stage."""
import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from runtime import ROOT, invoke, location, main

DEPENDENCIES = {
    "rosetta": {
        "exe": "ROSETTA_EXE", "mpirun": "MPI_RUN",
        "dalphaball": "DALPHABALL_EXE", "psipred": "PSIPRED_EXE",
    },
    "af2": {
        "python": "SCORE_PYTHON", "repo": "AF2_REPO",
        "silent_tools": "AF2_SILENT_TOOLS", "weights": "AF2_DATA_DIR",
    },
    "af3": {
        "python": "SCORE_PYTHON", "repo": "AF3_REPO",
        "weights": "AF3_WEIGHTS", "databases": "AF3_DATABASES",
        "pdb_database": "AF3_PDB_DATABASE",
    },
}


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=DEPENDENCIES)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    env = dict(os.environ)
    for key, variable in DEPENDENCIES[args.stage].items():
        env[variable] = str(location(args.stage, key, check=not args.dry_run))
    if args.dry_run:
        for variable in DEPENDENCIES[args.stage].values():
            print(f"{variable}={env[variable]}")
    script = ROOT / "score" / f"score_{args.stage}.sh"
    invoke(["bash", script, args.input.resolve(), args.output.resolve()],
           dry_run=args.dry_run, env=env)


if __name__ == "__main__":
    main(run)
