"""Insert sampled VHH sequences into the normalized benchmark backbone."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


AA1_TO_AA3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}
BACKBONE_ATOMS = {"N", "CA", "C", "O"}
PDB_RECORDS = {"ATOM  ", "HETATM"}


@dataclass
class Template:
    lines: list[str]
    vhh_residues: list[tuple[str, str, str]]


def residue_key(line: str) -> tuple[str, str, str]:
    """Return chain, residue number, insertion code using PDB columns."""
    return line[21:22], line[22:26], line[26:27]


def read_template(path: Path, vhh_chain: str) -> Template:
    lines = path.read_text().splitlines(keepends=True)
    residues: list[tuple[str, str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for line in lines:
        if not line.startswith("ATOM  "):
            continue
        if line[21:22] != vhh_chain or line[12:16].strip() != "CA":
            continue
        key = residue_key(line)
        if key not in seen:
            seen.add(key)
            residues.append(key)
    if not residues:
        raise ValueError(f"no VHH CA residues found in {path} on chain {vhh_chain!r}")
    return Template(lines=lines, vhh_residues=residues)


def backfill_pdb(template: Template, sequence: str, vhh_chain: str, output: Path) -> int:
    if any(aa not in AA1_TO_AA3 for aa in sequence):
        bad = sorted(set(sequence) - set(AA1_TO_AA3))
        raise ValueError(f"unsupported amino-acid code(s) {bad} for {output}")
    if len(sequence) != len(template.vhh_residues):
        raise ValueError(
            f"sequence/PDB VHH length mismatch for {output}: "
            f"sequence={len(sequence)}, PDB={len(template.vhh_residues)}"
        )

    residue_to_aa3 = dict(zip(template.vhh_residues, map(AA1_TO_AA3.get, sequence)))
    output_lines: list[str] = []
    for line in template.lines:
        if not line.startswith(tuple(PDB_RECORDS)):
            output_lines.append(line)
            continue
        if line[21:22] != vhh_chain:
            output_lines.append(line)
            continue
        if line.startswith("ATOM  "):
            key = residue_key(line)
            aa3 = residue_to_aa3.get(key)
            if aa3 is None:
                raise ValueError(f"VHH atom has no CA residue mapping in {output}: {line!r}")
            # Keep only the backbone of the designed VHH.  LigandMPNN will
            # reconstruct/pack its side chains for the backfilled sequence.
            if line[12:16].strip() not in BACKBONE_ATOMS:
                continue
            output_lines.append(line[:17] + aa3 + line[20:])
        else:
            # Do not discard an unusual HETATM on the VHH chain.
            output_lines.append(line)

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(output_lines))
    return len(sequence)


