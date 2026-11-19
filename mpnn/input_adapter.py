#!/usr/bin/env python3
"""Convert first-stage VHH--target structure outputs into MPNN inputs.

The four supported structure generators use different file names and chain IDs.
This adapter writes one deliberately small, common contract:

  output/
    pdbs/<design>.pdb        VHH=A; target=B
    manifest.json            actual, per-backbone CDR sequence positions
    parsed.jsonl
    chains.jsonl
    fixed_positions.jsonl

``chains.jsonl`` and ``fixed_positions.jsonl`` can be supplied unchanged to
ProteinMPNN, SolubleMPNN, or AntiBMPNN.  The latter three differ only in the
runner/model weights, not in this structural-input contract.

CDR positions are read from each generator's own output metadata, never
inferred by sequence alignment to a template.
"""
from __future__ import annotations

import argparse
import ast
import csv
import json
import re
import shlex
import struct
import subprocess
import sys
import zipfile
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
STRUCTURE_MODELS = ("boltzgen", "iggm", "ppiflow", "rfantibody")
CHAIN_IDS = "BCDEFGIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz"


@dataclass(frozen=True)
class Atom:
    group: str
    atom: str
    alt: str
    residue: str
    chain: str
    residue_id: str
    x: float
    y: float
    z: float
    occupancy: float = 1.0
    bfactor: float = 0.0
    element: str = ""


def fail(message: str) -> None:
    raise ValueError(message)


def read_pdb(path: Path) -> list[Atom]:
    atoms: list[Atom] = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        alt = line[16:17]
        # Avoid duplicate alternate conformers.  Blank/A is the conventional
        # primary conformer and is what ProteinMPNN's parser expects.
        if alt not in (" ", "A", "1"):
            continue
        try:
            atoms.append(Atom(
                line[:6].strip(), line[12:16].strip(), alt, line[17:20].strip(),
                line[21:22], line[22:27].strip(), float(line[30:38]),
                float(line[38:46]), float(line[46:54]),
                float(line[54:60] or 1.0), float(line[60:66] or 0.0),
                line[76:78].strip(),
            ))
        except (ValueError, IndexError) as exc:
            fail(f"{path}: cannot parse PDB ATOM record: {line!r} ({exc})")
    if not atoms:
        fail(f"{path}: no ATOM/HETATM records")
    return atoms


def read_cif(path: Path) -> list[Atom]:
    """Read the simple atom_site loop written by BoltzGen without extra deps."""
    lines = path.read_text().splitlines()
    start = next((i for i, line in enumerate(lines) if line.strip() == "loop_" and
                  i + 1 < len(lines) and lines[i + 1].startswith("_atom_site.")), None)
    if start is None:
        fail(f"{path}: mmCIF _atom_site loop not found")
    columns, cursor = [], start + 1
    while cursor < len(lines) and lines[cursor].startswith("_atom_site."):
        columns.append(lines[cursor].strip().split(".", 1)[1])
        cursor += 1
    needed = {"group_PDB", "label_atom_id", "label_alt_id", "label_comp_id",
              "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv",
              "auth_seq_id", "auth_asym_id"}
    missing = needed - set(columns)
    if missing:
        fail(f"{path}: atom_site columns absent: {sorted(missing)}")
    ix = {name: i for i, name in enumerate(columns)}
    atoms: list[Atom] = []
    for line in lines[cursor:]:
        if not line or line.startswith("#") or line == "loop_" or line.startswith("_"):
            break
        values = line.split() if not any(c in line for c in ("\"", "'", "\\")) else shlex.split(line)
        if len(values) != len(columns):
            fail(f"{path}: unsupported multiline mmCIF atom row: {line!r}")
        if values[ix["group_PDB"]] not in ("ATOM", "HETATM"):
            continue
        alt = values[ix["label_alt_id"]]
        if alt not in (".", "?", "A", "1"):
            continue
        chain = values[ix["auth_asym_id"]]
        if len(chain) != 1:
            fail(f"{path}: chain ID {chain!r} cannot be represented in PDB")
        atoms.append(Atom(
            values[ix["group_PDB"]], values[ix["label_atom_id"]], alt,
            values[ix["label_comp_id"]], chain, values[ix["auth_seq_id"]],
            float(values[ix["Cartn_x"]]), float(values[ix["Cartn_y"]]),
            float(values[ix["Cartn_z"]]), float(values[ix["occupancy"]]),
            float(values[ix["B_iso_or_equiv"]]),
            values[ix.get("type_symbol", ix["label_atom_id"])],
        ))
    if not atoms:
        fail(f"{path}: no supported atom_site rows")
    return atoms


def atoms_from(path: Path) -> list[Atom]:
    if path.suffix.lower() == ".pdb":
        return read_pdb(path)
    if path.suffix.lower() in (".cif", ".mmcif"):
        return read_cif(path)
    fail(f"unsupported structure extension: {path}")


def chains_and_sequence(atoms: list[Atom]) -> tuple[list[str], dict[str, str]]:
    chains, residues, seen = [], {}, set()
    for atom in atoms:
        if atom.chain not in chains:
            chains.append(atom.chain)
        if atom.atom != "CA":
            continue
        key = (atom.chain, atom.residue_id)
        if key not in seen:
            seen.add(key)
            residues.setdefault(atom.chain, []).append(AA3_TO_1.get(atom.residue.upper(), "X"))
    return chains, {chain: "".join(seq) for chain, seq in residues.items()}


def residue_order(atoms: list[Atom]) -> list[tuple[str, str]]:
    """Residues in output order, one entry per CA atom."""
    order, seen = [], set()
    for atom in atoms:
        key = (atom.chain, atom.residue_id)
        if atom.atom == "CA" and key not in seen:
            seen.add(key)
            order.append(key)
    return order


def local_residue_positions(atoms: list[Atom], vhh: str) -> dict[int, int]:
    """Map source PDB residue numbers to continuous 1-based VHH positions."""
    result = {}
    for chain, residue_id in residue_order(atoms):
        if chain != vhh:
            continue
        match = re.match(r"-?\d+", residue_id)
        if not match:
            fail(f"VHH residue ID {residue_id!r} is not numeric")
        number = int(match.group())
        if number in result:
            fail(f"duplicate VHH residue number {number}; cannot map metadata safely")
        result[number] = len(result) + 1
    return result


def split_three_runs(positions: list[int], source: Path, label: str) -> list[list[int]]:
    if not positions:
        fail(f"{source}: {label} contains no VHH CDR positions")
    positions = sorted(set(positions))
    groups = [[positions[0]]]
    for position in positions[1:]:
        if position == groups[-1][-1] + 1:
            groups[-1].append(position)
        else:
            groups.append([position])
    if len(groups) != 3:
        fail(f"{source}: expected three CDR runs in {label}, found {groups}")
    return groups


def load_npz_vector(path: Path, member: str) -> list[float]:
    """Read the primitive one-dimensional NPZ metadata vector without NumPy."""
    try:
        with zipfile.ZipFile(path) as archive:
            data = archive.read(member + ".npy")
    except (FileNotFoundError, KeyError, zipfile.BadZipFile) as exc:
        fail(f"{path}: cannot read {member!r} metadata ({exc})")
    if data[:6] != b"\x93NUMPY":
        fail(f"{path}: invalid NPY magic for {member}")
    major = data[6]
    if major == 1:
        header_len, offset = struct.unpack("<H", data[8:10])[0], 10
    elif major in (2, 3):
        header_len, offset = struct.unpack("<I", data[8:12])[0], 12
    else:
        fail(f"{path}: unsupported NPY version {major}")
    header = ast.literal_eval(data[offset:offset + header_len].decode("latin1"))
    if header.get("fortran_order") or len(header.get("shape", ())) != 1:
        fail(f"{path}: {member} must be a C-order 1-D array")
    formats = {"|b1": "?", "|u1": "B", "|i1": "b", "<u1": "B", "<i1": "b",
               "<f4": "f", "<f8": "d", "<i4": "i", "<i8": "q", "<u4": "I", "<u8": "Q"}
    dtype, count = header["descr"], header["shape"][0]
    if dtype not in formats:
        fail(f"{path}: unsupported {member} dtype {dtype!r}")
    fmt, payload = formats[dtype], data[offset + header_len:]
    expected = struct.calcsize("<" + fmt) * count
    if len(payload) != expected:
        fail(f"{path}: malformed {member} payload")
    return list(struct.unpack("<" + fmt * count, payload))


def boltzgen_cdrs(path: Path, atoms: list[Atom], vhh: str) -> list[list[int]]:
    """Use paired BoltzGen NPZ design_mask aligned to mmCIF token order."""
    mask = load_npz_vector(path.with_suffix(".npz"), "design_mask")
    order = residue_order(atoms)
    if len(mask) != len(order):
        fail(f"{path}: BoltzGen design_mask length {len(mask)} != residue count {len(order)}")
    local, positions = 0, []
    for designed, (chain, _) in zip(mask, order):
        if chain == vhh:
            local += 1
            if designed:
                positions.append(local)
    return split_three_runs(positions, path, "BoltzGen design_mask")


@lru_cache(maxsize=4)
def ppiflow_inputs(input_dir: Path) -> dict[str, dict[str, str]]:
    inputs = {}
    for csv_path in sorted((input_dir / "input").glob("*_input.csv")):
        with csv_path.open() as handle:
            for row in csv.DictReader(handle):
                name = row.get("pdb_name")
                if not name or name in inputs:
                    fail(f"{csv_path}: missing or repeated PPIFlow pdb_name {name!r}")
                inputs[name] = row
    return inputs


def ppiflow_input_row(path: Path, input_dir: Path) -> dict[str, str]:
    row = ppiflow_inputs(input_dir).get(path.stem)
    if row is None:
        fail(f"{path}: PPIFlow metadata missing in {input_dir / 'input'}")
    return row


def ppiflow_cdrs(path: Path, input_dir: Path, atoms: list[Atom], vhh: str) -> list[list[int]]:
    """Read cdr_H1/H2/H3 from PPIFlow's sample_metrics.csv."""
    metrics = input_dir / "sample_metrics.csv"
    if not metrics.is_file():
        fail(f"{path}: required PPIFlow metadata is missing: {metrics}")
    with metrics.open() as handle:
        row = next((row for row in csv.DictReader(handle) if row.get("sample") == path.name), None)
    if row is None:
        fail(f"{path}: no sample_metrics.csv row for this PDB")
    source_to_local = local_residue_positions(atoms, vhh)
    cdrs = []
    for key in ("cdr_H1", "cdr_H2", "cdr_H3"):
        try:
            source_positions = json.loads(row[key])
            cdrs.append([source_to_local[int(position)] for position in source_positions])
        except (KeyError, json.JSONDecodeError) as exc:
            fail(f"{metrics}: invalid {key} for {path.name}: {exc}")
    return cdrs


def rfantibody_cdrs(path: Path, atoms: list[Atom], vhh: str) -> list[list[int]]:
    """Read RFantibody's explicit H1/H2/H3 PDBinfo-LABEL records."""
    labels = {"H1": [], "H2": [], "H3": []}
    expression = re.compile(r"^REMARK PDBinfo-LABEL:\s*(\d+)\s+(H[123])\s*$")
    for line in path.read_text().splitlines():
        match = expression.match(line)
        if match:
            labels[match.group(2)].append(int(match.group(1)))
    if not all(labels.values()):
        fail(f"{path}: RFantibody H1/H2/H3 REMARK labels are required but missing")
    length = len(local_residue_positions(atoms, vhh))
    cdrs = [sorted(set(labels[name])) for name in ("H1", "H2", "H3")]
    if any(min(cdr) < 1 or max(cdr) > length for cdr in cdrs):
        fail(f"{path}: RFantibody HLT CDR labels do not map to its VHH chain")
    return cdrs


def detect_vhh_chain(source: str, path: Path, input_dir: Path, atoms: list[Atom],
                     iggm_design_fasta: Path | None) -> str:
    """Determine the VHH from source-owned metadata, not a fixed chain map."""
    source_chains, _ = chains_and_sequence(atoms)
    if source == "boltzgen":
        mask = load_npz_vector(path.with_suffix(".npz"), "design_mask")
        order = residue_order(atoms)
        if len(mask) != len(order):
            fail(f"{path}: BoltzGen design_mask length does not match its structure")
        designed_chains = {chain for designed, (chain, _) in zip(mask, order) if designed}
        if len(designed_chains) != 1:
            fail(f"{path}: BoltzGen design_mask identifies {designed_chains}, not one VHH chain")
        return designed_chains.pop()
    if source == "ppiflow":
        chain = ppiflow_input_row(path, input_dir).get("binder_id")
        if not chain or chain not in source_chains:
            fail(f"{path}: PPIFlow binder_id {chain!r} is absent from output chains {source_chains}")
        return chain
    if source == "rfantibody":
        # RFantibody's H1/H2/H3 PDBinfo labels define its heavy/VHH chain.
        if "H" not in source_chains:
            fail(f"{path}: RFantibody CDR labels require H chain; found {source_chains}")
        return "H"
    if iggm_design_fasta is None:
        fail(f"{path}: IgGM output needs --iggm-design-fasta to identify its VHH X-mask")
    if "H" not in read_fasta_sequences(iggm_design_fasta) or "H" not in source_chains:
        fail(f"{path}: IgGM VHH H chain is absent from its FASTA or output PDB")
    return "H"


def read_fasta_sequences(path: Path) -> dict[str, str]:
    sequences, name, chunks = {}, None, []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(chunks).upper()
            name, chunks = line[1:].split()[0], []
        else:
            if name is None:
                fail(f"{path}: sequence occurs before a FASTA header")
            chunks.append(line)
    if name is not None:
        sequences[name] = "".join(chunks).upper()
    return sequences


def iggm_cdrs(path: Path, atoms: list[Atom], vhh: str, design_fasta: Path | None) -> list[list[int]]:
    """Use the X mask from the exact FASTA supplied to IgGM design.py."""
    if design_fasta is None:
        fail(f"{path}: IgGM output does not retain its input X-mask; pass --iggm-design-fasta")
    if not design_fasta.is_file():
        fail(f"IgGM design FASTA does not exist: {design_fasta}")
    sequence = read_fasta_sequences(design_fasta).get(vhh)
    if sequence is None:
        fail(f"{design_fasta}: expected IgGM VHH FASTA header {vhh!r}")
    length = len(local_residue_positions(atoms, vhh))
    if len(sequence) != length:
        fail(f"{design_fasta}: H sequence length {len(sequence)} != generated VHH length {length}")
    return split_three_runs([i + 1 for i, aa in enumerate(sequence) if aa == "X"], path, "IgGM FASTA X-mask")


def source_cdrs(source: str, path: Path, input_dir: Path, atoms: list[Atom], vhh: str,
                iggm_design_fasta: Path | None) -> list[list[int]]:
    if source == "boltzgen":
        return boltzgen_cdrs(path, atoms, vhh)
    if source == "ppiflow":
        return ppiflow_cdrs(path, input_dir, atoms, vhh)
    if source == "rfantibody":
        return rfantibody_cdrs(path, atoms, vhh)
    return iggm_cdrs(path, atoms, vhh, iggm_design_fasta)


def cdr_metadata_reference(source: str, path: Path, input_dir: Path,
                           iggm_design_fasta: Path | None) -> tuple[str, str]:
    """Human-readable provenance recorded alongside every normalized design."""
    if source == "boltzgen":
        return "boltzgen_npz_design_mask", str(path.with_suffix(".npz"))
    if source == "ppiflow":
        return "ppiflow_sample_metrics_cdr_H1_H2_H3", str(input_dir / "sample_metrics.csv")
    if source == "rfantibody":
        return "rfantibody_pdbinfo_label_remarks", str(path) + "#H1/H2/H3"
    assert iggm_design_fasta is not None
    return "iggm_design_fasta_X_mask", str(iggm_design_fasta)


def write_pdb(path: Path, atoms: list[Atom], chain_map: dict[str, str]) -> None:
    serial, current_chain, output = 1, None, []
    for atom in atoms:
        chain = chain_map[atom.chain]
        if current_chain is not None and chain != current_chain:
            output.append(f"TER   {serial:5d}      {atom.residue:>3} {current_chain:1}\n")
            serial += 1
        current_chain = chain
        residue_number = int(re.match(r"-?\d+", atom.residue_id).group())
        insertion = atom.residue_id[len(str(residue_number)):][:1] or " "
        atom_name = atom.atom.rjust(4)
        element = (atom.element or atom.atom[:1]).strip().upper()[:2]
        output.append(
            f"{atom.group:<6}{serial:5d} {atom_name}{' ':1}{atom.residue:>3} {chain:1}"
            f"{residue_number:4d}{insertion:1}   {atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
            f"{atom.occupancy:6.2f}{atom.bfactor:6.2f}          {element:>2}\n"
        )
        serial += 1
    if current_chain is not None:
        output.append(f"TER   {serial:5d}      {atoms[-1].residue:>3} {current_chain:1}\n")
    output.append("END\n")
    path.write_text("".join(output))

