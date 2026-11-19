#!/usr/bin/env python3
"""Run the IgGM backbone atom repair with the configured PDBFixer interpreter."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from runtime import ROOT, invoke, location, main


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    python = location("pdbfixer", "python", check=not args.dry_run)
    invoke([python, ROOT / "mpnn/repair_iggm_backbone.py",
            "--input-dir", args.input_dir.resolve(), "--output-dir", args.output_dir.resolve()],
           dry_run=args.dry_run)


if __name__ == "__main__":
    main(run)
