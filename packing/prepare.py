#!/usr/bin/env python3
"""Validate sampled FASTA records, backfill the normalized PDBs and index them."""
import argparse
import csv
import json
import re
import sys
import tempfile
from pathlib import Path

from backfill import backfill_pdb, read_template

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from runtime import MODELS, main


def fasta(path):
    records = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            records.append((line[1:], []))
        elif line.strip():
            if not records:
                raise ValueError(f"Sequence before FASTA header: {path}")
            records[-1][1].append(line.strip())
    return [(header, "".join(seq)) for header, seq in records]


def prepare():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True, help="Prepared MPNN input directory")
    parser.add_argument("--sequences", type=Path, required=True, help="MPNN output directory")
    parser.add_argument("--model", choices=MODELS, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    inp = args.inputs.resolve()
    sources = json.loads((inp / "manifest.json").read_text())
    fixed = json.loads((inp / "fixed_positions.jsonl").read_text())
    parsed = {row["name"]: row for row in map(json.loads, (inp / "parsed.jsonl").read_text().splitlines())}
    files = {}
    for path in (args.sequences / "seqs").glob("*.fa"):
        name = path.stem.split("|")[0]
        if name in files:
            raise ValueError(f"Repeated FASTA for {name}")
        files[name] = path
    if not sources or set(files) != set(parsed) or set(parsed) != set(fixed):
        raise ValueError("FASTA, parsed backbones and fixed positions contain different ids")
    if len(sources) != len(parsed) or {row["name"] for row in sources} != set(parsed):
        raise ValueError("The source manifest must contain exactly one row per backbone")
    groups = {(row["flow"], row["target"], row["scaffold"]) for row in sources}
    if len(groups) != 1:
        raise ValueError("Use one flow, target and scaffold per packing input directory")
    flow, target, scaffold = next(iter(groups))
    dest = args.output_root.resolve() / flow / target / scaffold / args.model
    if dest.exists():
        raise FileExistsError(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sequence_id", "flow", "target", "scaffold", "model", "sequence", "target_sequence",
              "backbone_id", "sample", "source_pdb", "backfilled_pdb"]
    with tempfile.TemporaryDirectory(prefix=".backfill-", dir=dest.parent) as tmp:
        work = Path(tmp) / "output"
        work.mkdir()
        with (work / "manifest.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t")
            writer.writeheader()
            for row in sources:
                name = row["name"]
                native = parsed[name]["seq_chain_A"]
                records = fasta(files[name])
                if len(records) != 5 or records[0][1] != native:
                    raise ValueError(f"Expected native plus four samples for {name}")
                fixed_positions = {position - 1 for position in fixed[name]["A"]}
                template = read_template(inp / "pdbs" / f"{name}.pdb", "A")
                for index, (header, seq) in enumerate(records[1:], 1):
                    if (not re.search(r"\bsample=", header) or len(seq) != len(native)
                            or set(seq) - set("ACDEFGHIKLMNPQRSTVWY")
                            or any(seq[i] != native[i] for i in fixed_positions)
                            or any(aa == "C" for i, aa in enumerate(seq) if i not in fixed_positions)):
                        raise ValueError(f"Invalid sampled sequence {index} for {name}")
                    sample = f"sample_{index:04d}"
                    file = f"{name}__{sample}.pdb"
                    identifier = "__".join((flow, target, scaffold, args.model, name, sample))
                    backfill_pdb(template, seq, "A", work / file)
                    writer.writerow({
                        "sequence_id": identifier, "flow": flow, "target": target,
                        "scaffold": scaffold, "model": args.model, "sequence": seq,
                        "target_sequence": parsed[name]["seq_chain_B"], "backbone_id": name,
                        "sample": index, "source_pdb": row["source"], "backfilled_pdb": file,
                    })
        work.rename(dest)
    print(f"Prepared {4 * len(sources)} packed-input structures in {dest}")


if __name__ == "__main__":
    main(prepare)
