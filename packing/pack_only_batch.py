#!/usr/bin/env python3
"""Pack side chains for many PDBs with one persistent LigandMPNN packer.

The input amino-acid identities and backbone coordinates are retained.  All
protein side chains, on every protein chain, are rebuilt by the side-chain
packer.  Different PDBs with the same residue count are combined into a real
GPU batch; ``run.py --batch_size`` is intentionally not used because that
option replicates one structure rather than batching different structures.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch


SEQUENCE_MODELS = ("antibmpnn", "proteinmpnn", "solublempnn")


@dataclass
class Sample:
    input_path: Path
    output_path: Path
    feature: dict[str, torch.Tensor]
    chain_letters: np.ndarray
    icodes: np.ndarray
    other_atoms: object


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--flow",
        help="optional leading flow directory; mirrors it in the output layout",
    )
    parser.add_argument("--sequence-model", choices=SEQUENCE_MODELS, required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--scaffold", required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--ligandmpnn-repo", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--max-files-per-model", type=int)
    parser.add_argument("--denoising-steps", type=int, default=3)
    parser.add_argument("--num-samples", type=int, default=16)
    parser.add_argument("--seed", type=int, default=20260913)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--validate-existing", action="store_true")
    parser.add_argument("--parse-zero-occupancy", action="store_true")
    parser.add_argument("--force-hetatm", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    if args.batch_size < 1:
        parser.error("--batch-size must be positive")
    if args.max_files_per_model is not None and args.max_files_per_model < 1:
        parser.error("--max-files-per-model must be positive")
    return args


def load_ligandmpnn(repo: Path):
    repo = repo.resolve()
    if not (repo / "sc_utils.py").is_file():
        raise FileNotFoundError(f"invalid LigandMPNN repository: {repo}")
    sys.path.insert(0, str(repo))
    from data_utils import featurize, parse_PDB, write_full_PDB
    from sc_utils import Packer, pack_side_chains

    return featurize, parse_PDB, write_full_PDB, Packer, pack_side_chains


def output_is_complete(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    size = path.stat().st_size
    with path.open("rb") as handle:
        head = handle.read(4096)
        handle.seek(max(0, size - 4096))
        tail = handle.read()
    return b"ATOM  " in head and any(line.startswith(b"END") for line in tail.splitlines())


def discover(args: argparse.Namespace) -> list[tuple[Path, Path]]:
    pairs: list[tuple[Path, Path]] = []
    for sequence_model in (args.sequence_model,):
        input_dir = args.input_root
        output_dir = args.output_root
        if args.flow:
            input_dir /= args.flow
            output_dir /= args.flow
        input_dir = input_dir / args.target / args.scaffold / sequence_model
        if not input_dir.is_dir():
            raise FileNotFoundError(f"input directory does not exist: {input_dir}")
        paths = sorted(input_dir.glob("*.pdb"))
        if args.max_files_per_model is not None:
            paths = paths[: args.max_files_per_model]
        output_dir = output_dir / args.target / args.scaffold / sequence_model
        for input_path in paths:
            output_path = output_dir / f"{input_path.stem}_packed.pdb"
            if output_path.exists() and not args.overwrite:
                if not args.validate_existing or output_is_complete(output_path):
                    continue
            pairs.append((input_path, output_path))
    return pairs


def make_model(Packer, checkpoint_path: Path, device: torch.device):
    model = Packer(
        node_features=128,
        edge_features=128,
        num_positional_embeddings=16,
        num_chain_embeddings=16,
        num_rbf=16,
        hidden_dim=128,
        num_encoder_layers=3,
        num_decoder_layers=3,
        atom_context_num=16,
        lower_bound=0.0,
        upper_bound=20.0,
        top_k=32,
        dropout=0.0,
        augment_eps=0.0,
        atom37_order=False,
        device=device,
        num_mix=3,
    )
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.to(device).eval()
    return model


def parse_sample(pair, args, parse_PDB, featurize) -> Sample:
    input_path, output_path = pair
    protein, _backbone, other_atoms, icodes, _water = parse_PDB(
        str(input_path),
        device="cpu",
        parse_all_atoms=False,
        parse_atoms_with_zero_occupancy=args.parse_zero_occupancy,
    )
    if protein["S"].numel() == 0:
        raise ValueError("no protein residues")
    if torch.any(protein["S"] == 20):
        raise ValueError("non-standard/unknown protein residue found")
    if not torch.all(protein["mask"] == 1):
        missing = int(torch.sum(protein["mask"] == 0))
        raise ValueError(f"{missing} residues lack one or more N/CA/C/O atoms")
    protein["chain_mask"] = torch.ones_like(protein["mask"])
    # OpenFold's atom14 conversion calls one_hot, which requires int64 indices.
    protein["S"] = protein["S"].long()
    feature = featurize(
        protein,
        cutoff_for_score=8.0,
        use_atom_context=True,
        number_of_ligand_atoms=16,
        model_type="ligand_mpnn",
    )
    return Sample(
        input_path=input_path,
        output_path=output_path,
        feature=feature,
        chain_letters=np.asarray(protein["chain_letters"]),
        icodes=np.asarray(icodes),
        other_atoms=other_atoms,
    )


def collate(samples: list[Sample], device: torch.device) -> dict[str, torch.Tensor]:
    keys = samples[0].feature.keys()
    if any(sample.feature.keys() != samples[0].feature.keys() for sample in samples):
        raise ValueError("feature keys differ within batch")
    batch: dict[str, torch.Tensor] = {}
    for key in keys:
        values = [sample.feature[key] for sample in samples]
        if not all(torch.is_tensor(value) for value in values):
            continue
        batch[key] = torch.cat(values, dim=0).to(device, non_blocking=True)
    return batch


def stable_batch_seed(base_seed: int, samples: list[Sample]) -> int:
    names = "\n".join(str(sample.input_path) for sample in samples).encode()
    digest = int.from_bytes(hashlib.sha256(names).digest()[:8], "little")
    return (base_seed + digest) % (2**63 - 1)


def append_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    new_file = not path.exists()
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        if new_file:
            writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())


def write_batch(samples, packed, write_full_PDB, args, completed_path):
    fields = ["input_pdb", "output_pdb", "residues", "batch_size", "seconds"]
    started = time.monotonic()
    completed_rows = []
    for index, sample in enumerate(samples):
        sample.output_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = sample.output_path.with_name(sample.output_path.stem + ".tmp.pdb")
        write_full_PDB(
            str(tmp_path),
            packed["X"][index].detach().cpu().numpy(),
            packed["X_m"][index].detach().cpu().numpy(),
            packed["b_factors"][index].detach().cpu().numpy(),
            sample.feature["R_idx_original"][0].cpu().numpy(),
            sample.chain_letters,
            sample.feature["S"][0].cpu().numpy(),
            other_atoms=sample.other_atoms,
            icodes=sample.icodes,
            force_hetatm=args.force_hetatm,
        )
        if not output_is_complete(tmp_path):
            raise RuntimeError(f"incomplete temporary output: {tmp_path}")
        os.replace(tmp_path, sample.output_path)
        completed_rows.append({
            "input_pdb": sample.input_path,
            "output_pdb": sample.output_path,
            "residues": sample.feature["S"].shape[1],
            "batch_size": len(samples),
            "seconds": f"{time.monotonic() - started:.3f}",
        })
    append_tsv(completed_path, fields, completed_rows)


def main() -> None:
    args = parse_args()
    featurize, parse_PDB, write_full_PDB, Packer, pack_side_chains = load_ligandmpnn(
        args.ligandmpnn_repo
    )
    if not torch.cuda.is_available():
        raise SystemExit("CUDA GPU is required")
    device = torch.device("cuda")
    checkpoint = args.ligandmpnn_repo / "model_params/ligandmpnn_sc_v_32_002_16.pt"
    if not checkpoint.is_file():
        raise FileNotFoundError(checkpoint)

    pairs = discover(args)
    run_root = args.output_root
    if args.flow:
        run_root /= args.flow
    run_root = run_root / args.target / args.scaffold / args.sequence_model
    completed_path = run_root / "completed.tsv"
    failed_path = run_root / "failed.tsv"
    print(
        f"flow={args.flow or '(none)'} target={args.target} scaffold={args.scaffold} "
        f"pending={len(pairs)} "
        f"batch_size={args.batch_size} device={torch.cuda.get_device_name(0)}",
        flush=True,
    )
    if not pairs:
        return

    model = make_model(Packer, checkpoint, device)
    pending: dict[int, list[Sample]] = defaultdict(list)
    completed = 0
    failed = 0
    start = time.monotonic()

    def process(samples: list[Sample]) -> None:
        nonlocal completed, failed
        if not samples:
            return
        batch_started = time.monotonic()
        batch = None
        packed = None
        input_backbone = None
        try:
            torch.manual_seed(stable_batch_seed(args.seed, samples))
            batch = collate(samples, device)
            input_backbone = batch["X"].clone()
            with torch.inference_mode():
                packed = pack_side_chains(
                    batch,
                    model,
                    args.denoising_steps,
                    args.num_samples,
                    repack_everything=True,
                )
                # OpenFold rebuilds atom14 coordinates from rigid frames and
                # can perturb N/C/O by tenths of an Angstrom.  This workflow
                # is side-chain-only, so restore all four input backbone atoms.
                packed["X"][:, :, :4] = input_backbone
            if not torch.isfinite(packed["X"]).all():
                raise RuntimeError("packer returned non-finite coordinates")
            write_batch(samples, packed, write_full_PDB, args, completed_path)
            completed += len(samples)
            elapsed = time.monotonic() - batch_started
            peak_gib = torch.cuda.max_memory_allocated() / 1024**3
            print(
                f"completed={completed}/{len(pairs)} batch={len(samples)} "
                f"L={samples[0].feature['S'].shape[1]} seconds={elapsed:.2f} "
                f"peak_cuda_GiB={peak_gib:.2f}",
                flush=True,
            )
            del batch, packed, input_backbone
            torch.cuda.reset_peak_memory_stats()
        except torch.cuda.OutOfMemoryError:
            del batch, packed, input_backbone
            torch.cuda.empty_cache()
            if len(samples) == 1:
                raise
            midpoint = len(samples) // 2
            print(f"CUDA OOM: retrying {len(samples)} as {midpoint}+{len(samples)-midpoint}", flush=True)
            process(samples[:midpoint])
            process(samples[midpoint:])
        except Exception as exc:
            failed += len(samples)
            for sample in samples:
                append_tsv(failed_path, ["input_pdb", "error"], [{
                    "input_pdb": sample.input_path,
                    "error": f"{type(exc).__name__}: {exc}",
                }])
            print(f"batch failed ({len(samples)} files): {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)

    for pair in pairs:
        try:
            sample = parse_sample(pair, args, parse_PDB, featurize)
            length = sample.feature["S"].shape[1]
            pending[length].append(sample)
            if len(pending[length]) >= args.batch_size:
                batch_samples = pending.pop(length)
                process(batch_samples)
        except Exception as exc:
            failed += 1
            append_tsv(failed_path, ["input_pdb", "error"], [{
                "input_pdb": pair[0],
                "error": f"{type(exc).__name__}: {exc}",
            }])
            print(f"parse failed: {pair[0]}: {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)

    for length in sorted(pending):
        process(pending[length])

    elapsed = time.monotonic() - start
    print(
        f"DONE completed={completed} failed={failed} elapsed_seconds={elapsed:.2f} "
        f"structures_per_second={completed / elapsed if elapsed else 0:.4f}",
        flush=True,
    )
    if failed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
