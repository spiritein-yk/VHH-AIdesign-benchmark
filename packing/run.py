#!/usr/bin/env python3
"""Run the released side-chain packer with its configured interpreter."""
import argparse
import fcntl
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from runtime import FLOWS, MODELS, ROOT, TARGETS, invoke, location, main

CHECKPOINT = "model_params/ligandmpnn_sc_v_32_002_16.pt"
CHECKPOINT_SHA256 = "799a9de6c1c72bb0c7ceb37998391a6fd0e3e21cb42928ee31a28313a3b1b46a"


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", choices=TARGETS, required=True)
    parser.add_argument("--flow", choices=FLOWS, required=True)
    parser.add_argument("--scaffold", required=True)
    parser.add_argument("--model", choices=MODELS, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    check = not args.dry_run
    python = location("packing", "python", check=check)
    repo = location("packing", "repo", check=check)
    output_root = args.output_root.resolve()
    command = [
        python, ROOT / "packing/pack_only_batch.py",
        "--target", args.target, "--flow", args.flow, "--scaffold", args.scaffold,
        "--sequence-model", args.model, "--input-root", args.input_root.resolve(),
        "--output-root", output_root, "--ligandmpnn-repo", repo, "--validate-existing",
    ]
    if args.dry_run:
        invoke(command, dry_run=True)
        return
    with (repo / CHECKPOINT).open("rb") as weight:
        if hashlib.file_digest(weight, "sha256").hexdigest() != CHECKPOINT_SHA256:
            raise ValueError(f"Unexpected checkpoint checksum: {repo / CHECKPOINT}")
    group = output_root / args.flow / args.target / args.scaffold / args.model
    group.mkdir(parents=True, exist_ok=True)
    with (group / ".packing.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        invoke(command)


if __name__ == "__main__":
    main(run)
