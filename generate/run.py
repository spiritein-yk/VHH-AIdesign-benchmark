#!/usr/bin/env python3
"""Run one benchmark generator for one of the released scaffolds."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from runtime import FLOWS, ROOT, TARGETS, invoke, location, main


def generate(args, release, metadata):
    flow = args.flow
    out = args.output_dir
    check = not args.dry_run

    def external(key):
        return location(flow, key, check=check)

    scaffold = args.scaffold
    antigen = release / "common" / f"{args.target}.pdb"
    hotspots = metadata["hotspot_positions_1based"]
    chain_hotspots = ",".join(f"B{position}" for position in hotspots)
    loops = metadata["scaffolds"][scaffold]["variable_loop_length_ranges"]
    ranges = [f"{low}-{high}" for low, high in loops]
    if flow == "boltzgen":
        command = [
            external("exe"), "run", release / flow / f"{scaffold}.yaml",
            "--output", out, "--protocol", "nanobody-anything",
            "--num_designs", args.count, "--diffusion_batch_size", 1,
            "--inverse_fold_num_sequences", 4, "--config", "design", "sampling_steps=500",
            "--steps", "design", "inverse_folding", "--cache", external("cache"),
        ]
        return command, ROOT
    repo = external("repo")
    if flow == "ppiflow":
        command = [
            external("python"), repo / "tools/PPIFlow/sample_antibody_nanobody.py",
            "--antigen_pdb", antigen, "--antigen_chain", "B",
            "--framework_pdb", release / flow / f"{scaffold}_framework.pdb", "--heavy_chain", "A",
            "--specified_hotspots", chain_hotspots, "--cdr_length",
            ",".join(f"CDRH{i},{item}" for i, item in enumerate(ranges, 1)),
            "--samples_per_target", args.count, "--model_weights", external("weights"),
            "--output_dir", out, "--name", f"{args.target}_{scaffold}",
            "--config", release / flow / "inference_nanobody.yaml",
        ]
    elif flow == "rfantibody":
        command = [
            external("exe"), "--target", antigen,
            "--framework", release / flow / f"{scaffold}_framework_H.pdb",
            "--output", out / "design", "--num-designs", args.count,
            "--design-loops", ",".join(f"H{i}:{item}" for i, item in enumerate(ranges, 1)),
            "--hotspots", chain_hotspots, "--diffuser-t", 50, "--final-step", 1,
            "--no-trajectory", "-e", "inference.design_startnum=0", "-e", "inference.cautious=True",
        ]
    else:
        command = [
            external("python"), repo / "design.py",
            "--fasta", release / flow / f"{scaffold}.fasta",
            "--antigen", release / flow / f"{args.target}.pdb", "--output", out,
            "--epitope", *hotspots, "--steps", 10, "--chunk_size", 64,
            "--temperature", 1.0, "--num_samples", args.count, "--max_antigen_size", 2000,
            "--run_task", "design",
        ]
    return command, repo


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", choices=TARGETS, required=True)
    parser.add_argument("--scaffold", required=True)
    parser.add_argument("--flow", choices=FLOWS, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--count", type=int, default=5000)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    args.output_dir = args.output_dir.resolve()
    if args.count < 1:
        parser.error("--count must be positive")
    if not args.dry_run and args.output_dir.exists():
        parser.error(f"Output directory exists: {args.output_dir}")
    release = ROOT / "configs" / args.target
    metadata = json.loads((release / "validation.json").read_text())
    if args.scaffold not in metadata["scaffolds"]:
        parser.error(f"Unknown scaffold for {args.target}: {args.scaffold}")
    command, cwd = generate(args, release, metadata)
    if not args.dry_run:
        args.output_dir.mkdir(parents=True)
    invoke(command, cwd=cwd, dry_run=args.dry_run)


if __name__ == "__main__":
    main(run)
