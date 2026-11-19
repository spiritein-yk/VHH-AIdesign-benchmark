#!/usr/bin/env python3
"""Write a binary silent archive and its canonical design ID mapping."""
import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mpnn"))
from input_adapter import atoms_from, chains_and_sequence
from runtime import location, main


def backbone(atoms):
    return {
        (atom.chain, atom.residue_id, atom.atom): (atom.x, atom.y, atom.z)
        for atom in atoms
        if atom.atom in {"N", "CA", "C", "O"}
    }


def archive():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--packed-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    dest = args.output_dir.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    tools = location("silent", "tools")
    python = location("silent", "python")
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join((str(python.parent), str(tools), env.get("PATH", "")))
    with args.manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    names = {Path(row["backfilled_pdb"]).stem + "_packed": row for row in rows}
    if not rows or len(names) != len(rows) or len({row["sequence_id"] for row in rows}) != len(rows):
        raise ValueError("Expected a nonempty manifest with unique design ids and structure names")
    packed = {path.stem: path.resolve() for path in args.packed_dir.glob("*_packed.pdb")}
    if set(packed) != set(names):
        raise ValueError("Packed structure names do not match the manifest")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".silent-", dir=dest.parent) as tmp:
        work = Path(tmp)
        release = work / "archive"
        release.mkdir()
        staged = work / "pdbs"
        staged.mkdir()
        for tag, path in packed.items():
            row = names[tag]
            atoms = atoms_from(path)
            _, sequence = chains_and_sequence(atoms)
            if sequence != {"A": row["sequence"], "B": row["target_sequence"]}:
                raise ValueError(f"Packed sequence or chain identity differs: {path}")
            expected_coords = backbone(atoms_from(args.manifest.parent / row["backfilled_pdb"]))
            packed_coords = backbone(atoms)
            if (
                not expected_coords
                or expected_coords.keys() != packed_coords.keys()
                or any(
                    abs(a - b) > 0.002
                    for key, xyz in expected_coords.items()
                    for a, b in zip(xyz, packed_coords[key])
                )
            ):
                raise ValueError(f"Packed backbone changed beyond coordinate tolerance: {path}")
            (staged / path.name).symlink_to(path)

        def tool(name, *args, stdin=None, stdout=subprocess.PIPE):
            return subprocess.run([tools / name, *args], cwd=work, env=env, check=True,
                                  text=True, input=stdin, stdout=stdout).stdout

        (work / "flags").write_text("-mute all\n")
        with (work / "raw.silent").open("wb") as handle:
            tool("silentfrompdbs", "@flags", stdin="".join(f"pdbs/{name}.pdb\n" for name in names),
                 stdout=handle)
        original = tool("silentls", "raw.silent").splitlines()
        if len(original) != len(rows) or set(original) != set(names):
            raise ValueError("The converted archive has missing or repeated tags")
        expected = [names[name]["sequence_id"] for name in original]
        with (release / "packed.silent").open("wb") as handle:
            tool("silentrename", "raw.silent", stdin="\n".join(expected) + "\n", stdout=handle)
        if tool("silentls", "archive/packed.silent").splitlines() != expected:
            raise ValueError("The final silent archive has unexpected tags")
        fields = [*rows[0], "packed_pdb"]
        with (release / "tag_map.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t")
            writer.writeheader()
            for name in original:
                writer.writerow({**names[name], "packed_pdb": os.path.relpath(packed[name], dest)})
        (release / "validation.json").write_text(json.dumps({"count": len(expected), "status": "PASS"}) + "\n")
        release.rename(dest)


if __name__ == "__main__":
    main(archive)
