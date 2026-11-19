#!/usr/bin/env python3
"""Run the optional native folding/analysis or cross-template filtering."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from runtime import ROOT, TARGETS, invoke, location, main

FILTERING = ["--config", "filtering", "modality=antibody", "top_budget=500",
             "filter_biased=true", "filter_cysteine=true", "filter_designfolding=false",
             "filter_bindingsite=false", "use_affinity=false"]


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("analyse", "select"))
    parser.add_argument("--target", choices=TARGETS, required=True)
    parser.add_argument("--scaffold", required=True, help="A released BoltzGen specification")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--source", type=Path, action="append", required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    exe = location("boltzgen", "exe", check=not args.dry_run)
    cache = location("boltzgen", "cache", check=not args.dry_run)
    spec = ROOT / "configs" / args.target / "boltzgen" / f"{args.scaffold}.yaml"
    if not spec.is_file() or spec.parent != (ROOT / "configs" / args.target / "boltzgen"):
        parser.error(f"Unknown released specification: {args.scaffold}")
    sources = [path.resolve() for path in args.source]
    out = args.output_dir.resolve()
    if not args.dry_run and out.exists():
        raise FileExistsError(out)
    if args.action == "analyse":
        if len(sources) != 1:
            parser.error("analyse expects one generated inverse-fold output directory")
        raw = sources[0] / "intermediate_designs_inverse_folded"
        paired = {p.stem for p in raw.glob("*.cif")}
        if not args.dry_run:
            if not paired or paired != {p.stem for p in raw.glob("*.npz")}:
                raise ValueError(f"Expected paired CIF and NPZ in {raw}")
            (out / raw.name).mkdir(parents=True)
            for path in raw.iterdir():
                if path.stem in paired and path.suffix in (".cif", ".npz"):
                    (out / raw.name / path.name).symlink_to(path)
        steps = ["folding", "analysis"]
        extra = [
            "--num_workers", 4, "--config", "folding", "diffusion_samples=5", "recycling_steps=3",
            "sampling_steps=200", "--config", "analysis", "num_processes=8", "liability_modality=antibody",
        ]
    else:
        if len(sources) != 4 or len(set(sources)) != 4:
            parser.error("select expects four distinct analysed scaffold directories")
        invoke([exe, "merge", *sources, "--output", out], dry_run=args.dry_run)
        steps = ["filtering"]
        extra = []
    invoke([exe, "run", spec, "--output", out, "--protocol", "nanobody-anything",
            "--steps", *steps, "--budget", 20, "--alpha", 0.001, "--reuse",
            "--cache", cache, *extra, *FILTERING], dry_run=args.dry_run)


if __name__ == "__main__":
    main(run)
