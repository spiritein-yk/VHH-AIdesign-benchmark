#!/usr/bin/env python3
"""Prepare the Rosetta-to-AF2 and AF2-to-AF3 inputs and preserve source IDs."""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from table_io import read_table, value, write_table

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from runtime import TARGETS, location, main


def af2_input(args):
    target = args.output.resolve()
    if target.exists():
        raise FileExistsError(target)
    ids = [value(row, "source_id") for row in read_table(args.passed)]
    if not ids:
        raise ValueError("No Rosetta-passing candidates to send to AF2")
    tools = location("silent", "tools")
    python = location("silent", "python")
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join((str(python.parent), str(tools), env.get("PATH", "")))
    target.parent.mkdir(parents=True, exist_ok=True)
    src = args.silent.resolve(strict=True)
    with tempfile.TemporaryDirectory(dir=target.parent) as temp:
        work = Path(temp)
        (work / "input.silent").symlink_to(src)
        with (work / "af2.silent").open("wb") as out:
            subprocess.run([tools / "silentslice", "input.silent"], check=True, cwd=work, env=env,
                           input="\n".join(ids) + "\n", text=True, stdout=out)
        result = subprocess.run([tools / "silentls", "af2.silent"], check=True, cwd=work,
                                env=env, text=True, capture_output=True)
        selected = result.stdout.splitlines()
        if len(selected) != len(ids) or set(selected) != set(ids):
            raise ValueError("The extracted archive differs from the Rosetta selection")
        os.link(work / "af2.silent", target)
    print(f"Prepared {len(ids)} initial-guess inputs")


def indexed(path, key):
    rows = {}
    for row in read_table(path):
        identifier = key(row)
        if identifier in rows:
            raise ValueError(f"Duplicate ID {identifier} in {path}")
        rows[identifier] = row
    return rows


AF3_FIELDS = ["source_id", "flow", "target", "scaffold", "mpnn", "pae_interaction", "ddg",
              "vhh_sequence", "target_sequence", "original_design_pdb"]


def af3_scores(args):
    rows = []
    seen = set()
    for tag_map, rosetta, af2 in args.group:
        root = Path(tag_map).resolve().parent
        tags = indexed(tag_map, lambda row: row["sequence_id"])
        scores = indexed(rosetta, lambda row: value(row, "source_id"))
        for af2_row in read_table(af2):
            tag = value(af2_row, "source_id")
            if not tag.endswith("_af2pred"):
                raise ValueError(f"Unrecognized AF2 prediction tag: {tag}")
            source = tag.removesuffix("_af2pred")
            if source in seen or source not in tags or source not in scores:
                raise ValueError(f"Duplicate or unmapped AF2 source: {source}")
            seen.add(source)
            row = tags[source]
            if row["target"] != args.target:
                raise ValueError(f"Unexpected target in {tag_map}: {row['target']}")
            rows.append({
                "source_id": source, "flow": row["flow"], "target": row["target"],
                "scaffold": row["scaffold"], "mpnn": row["model"],
                "pae_interaction": value(af2_row, "pae_interaction"),
                "ddg": value(scores[source], "ddg"), "vhh_sequence": row["sequence"],
                "target_sequence": row["target_sequence"],
                "original_design_pdb": os.path.relpath(root / row["packed_pdb"], args.output.resolve().parent),
            })
    write_table(args.output, rows, AF3_FIELDS)
    print(f"Joined {len(rows)} AF2 scores to original packed structures")


def af3_input(args):
    dest = args.output_dir.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    rows = read_table(args.candidates)
    if not rows:
        raise ValueError("No candidates were admitted to AF3")
    msa = args.target_msa.read_text() if args.target_msa else None
    targets = {row["target_sequence"] for row in rows}
    if len(targets) != 1 or {row["target"] for row in rows} != {args.target}:
        raise ValueError("An AF3 input batch must use one released target and target sequence")
    target_sequence = next(iter(targets))
    if msa is not None:
        blocks = msa.split(">")
        query = "".join(blocks[1].splitlines()[1:]).strip() if len(blocks) > 1 else ""
        if query != target_sequence:
            raise ValueError("Target MSA query does not match the target sequence")
    dest.parent.mkdir(parents=True, exist_ok=True)
    sources = []
    tasks = {}
    with tempfile.TemporaryDirectory(prefix=".af3-", dir=dest.parent) as tmp:
        work = Path(tmp) / "output"
        (work / "jobs").mkdir(parents=True)
        for row in rows:
            vhh = row["vhh_sequence"]
            digest = hashlib.sha256(f"{vhh}\n{target_sequence}".encode()).hexdigest()[:16]
            name = f"{args.target.lower()}_{digest}"
            if name in tasks and tasks[name] != vhh:
                raise ValueError(f"AF3 task name collision: {name}")
            if name not in tasks:
                tasks[name] = vhh
                antigen = {"id": "B", "sequence": target_sequence}
                if msa is not None:
                    antigen.update(unpairedMsa=msa, pairedMsa="", templates=[])
                payload = {
                    "name": name, "modelSeeds": [1], "sequences": [
                        {"protein": {"id": "A", "sequence": vhh}}, {"protein": antigen},
                    ], "dialect": "alphafold3", "version": 1,
                }
                (work / "jobs" / f"{name}.json").write_text(json.dumps(payload, indent=2) + "\n")
            source = {field: row[field] for field in AF3_FIELDS}
            packed = args.candidates.resolve().parent / source["original_design_pdb"]
            source["original_design_pdb"] = os.path.relpath(packed, dest)
            source["af3_result_dir"] = f"predictions/{name}"
            sources.append(source)
        write_table(work / "sources.tsv", sources, [*AF3_FIELDS, "af3_result_dir"])
        work.rename(dest)
    print(f"Prepared {len(tasks)} AF3 tasks for {len(sources)} source designs")


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="stage", required=True)
    extract = sub.add_parser("af2-input")
    extract.add_argument("--passed", type=Path, required=True)
    extract.add_argument("--silent", type=Path, required=True)
    extract.add_argument("--output", type=Path, required=True)
    extract.set_defaults(action=af2_input)
    join = sub.add_parser("af3-scores")
    join.add_argument("--target", choices=TARGETS, required=True)
    join.add_argument("--group", nargs=3, action="append", required=True,
                      metavar=("TAG_MAP", "ROSETTA_SCORES", "AF2_SCORES"))
    join.add_argument("--output", type=Path, required=True)
    join.set_defaults(action=af3_scores)
    jobs = sub.add_parser("af3-input")
    jobs.add_argument("--target", choices=TARGETS, required=True)
    jobs.add_argument("--candidates", type=Path, required=True)
    jobs.add_argument("--target-msa", type=Path)
    jobs.add_argument("--output-dir", type=Path, required=True)
    jobs.set_defaults(action=af3_input)
    args = parser.parse_args()
    args.action(args)


if __name__ == "__main__":
    main(run)
