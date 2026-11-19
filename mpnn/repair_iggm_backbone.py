#!/usr/bin/env python3
"""Add missing N/CA/C/O only; preserve every existing atom and coordinate."""
import argparse
from pathlib import Path
import input_adapter as a
from pdbfixer import PDBFixer
from openmm import unit

def main():
 p=argparse.ArgumentParser();p.add_argument('--input-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);args=p.parse_args()
 files=sorted(args.input_dir.glob('*_CDR_All_*.pdb'))
 if not files:raise ValueError('No IgGM structures')
 if args.output_dir.exists() and any(args.output_dir.iterdir()):raise ValueError('Refusing to overwrite existing files')
 args.output_dir.mkdir(parents=True,exist_ok=True)
 for path in files:
  atoms=a.atoms_from(path)
  fixer=PDBFixer(filename=str(path));fixer.findMissingResidues();fixer.missingResidues={};fixer.findMissingAtoms()
  fixer.missingAtoms={r:[x for x in ats if x.name in {'N','CA','C','O'}] for r,ats in fixer.missingAtoms.items()}
  fixer.missingAtoms={r:ats for r,ats in fixer.missingAtoms.items() if ats};fixer.missingTerminals={}
  wanted={(r.chain.id,str(r.id)+(r.insertionCode or ''),x.name) for r,ats in fixer.missingAtoms.items() for x in ats}
  if wanted:
   fixer.addMissingAtoms(seed=42)
   coords=fixer.positions.value_in_unit(unit.angstrom)
   for atom in fixer.topology.atoms():
    rid=str(atom.residue.id)+(atom.residue.insertionCode or '')
    if (atom.residue.chain.id,rid,atom.name) in wanted:
     xyz=coords[atom.index]
     atoms.append(a.Atom('ATOM',atom.name,' ',atom.residue.name,atom.residue.chain.id,rid,float(xyz[0]),float(xyz[1]),float(xyz[2]),element=atom.element.symbol))
  order=a.residue_order(atoms);grouped={key:[] for key in order}
  for atom in atoms:grouped.setdefault((atom.chain,atom.residue_id),[]).append(atom)
  if any(not {'N','CA','C','O'}<={x.atom for x in ats} for ats in grouped.values()):raise ValueError('Incomplete backbone after repair')
  # Write from original records and append only newly created atom records,
  # so existing numeric coordinates and metadata remain exactly as supplied.
  text=[s for s in path.read_text().splitlines() if s not in ('END','ENDMDL')]
  additions=[atom for atom in atoms if (atom.chain,atom.residue_id,atom.atom) in wanted]
  tmp=args.output_dir/(path.name+'.newatoms')
  if additions:
   a.write_pdb(tmp,additions,{x.chain:x.chain for x in additions})
   serial=max((int(s[6:11]) for s in text if s.startswith(('ATOM  ','HETATM','TER   ')) and s[6:11].strip().isdigit()),default=0)
   for line in tmp.read_text().splitlines():
    if line.startswith('ATOM  '):
     serial+=1;text.append(line[:6]+f'{serial:5d}'+line[11:])
   tmp.unlink()
  (args.output_dir/path.name).write_text('\n'.join(text+['END'])+'\n')
 print('Prepared',len(files),'IgGM structures; existing coordinates unchanged')
if __name__=='__main__':main()
