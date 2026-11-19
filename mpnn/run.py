#!/usr/bin/env python3
"""Generate sequences with one of the three benchmark MPNNs."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from runtime import MODELS, invoke, location, main


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=MODELS, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    model = args.model
    check = not args.dry_run
    python = location(model, "python", check=check)
    repo = location(model, "repo", check=check)
    inputs = args.inputs.resolve()
    out = args.output_dir.resolve()
    if check:
        for filename in ("parsed.jsonl", "chains.jsonl", "fixed_positions.jsonl", "manifest.json"):
            if not (inputs / filename).is_file():
                raise ValueError(f"Missing {inputs / filename}. Run mpnn/prepare_inputs.py first")
        if out.exists():
            raise ValueError(f"Output directory exists: {out}")
    anti = model == "antibmpnn"
    cwd = repo / "example" if anti else repo
    entry = repo / ("Running_AntiBMPNN_run.py" if anti else "protein_mpnn_run.py")
    command = [
        python, entry, "--jsonl_path", inputs / "parsed.jsonl",
        "--chain_id_jsonl", inputs / "chains.jsonl",
        "--fixed_positions_jsonl", inputs / "fixed_positions.jsonl",
        "--num_seq_per_target", 4, "--sampling_temp", 0.1, "--batch_size", 4,
        "--backbone_noise", 0.0, "--omit_AAs", "CX",
        "--model_name", "antibmpnn_000" if anti else "v_48_020", "--out_folder", out,
    ]
    if not anti:
        weights = location(model, "weights", check=check)
        command += ["--path_to_model_weights", weights, "--seed", 43 if model == "solublempnn" else 42]
        if model == "solublempnn":
            command.append("--use_soluble_model")
    invoke(command, cwd=cwd, dry_run=args.dry_run)
    if check:
        manifest = json.loads((inputs / "manifest.json").read_text())
        if len(list((out / "seqs").glob("*.fa"))) != len(manifest):
            raise ValueError(f"MPNN output count differs from {inputs / 'manifest.json'}")


if __name__ == "__main__":
    main(run)
