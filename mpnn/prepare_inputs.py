#!/usr/bin/env python3
"""Prepare common MPNN JSONL from actual per-design generator metadata."""
import argparse,dataclasses,hashlib,json,math,re
from pathlib import Path
import input_adapter as a

def prepare(args):
 root=args.config_dir.resolve()
 info=json.loads((root/'validation.json').read_text())
 if args.scaffold not in info['scaffolds']:raise ValueError('Unknown scaffold')
 ti=info['scaffolds'][args.scaffold]
 target_atoms=a.atoms_from(root/'common'/f'{args.target}.pdb')
 target_seq=a.chains_and_sequence(target_atoms)[1]['B']
 scaffold_seq=a.chains_and_sequence(a.atoms_from(root/f'common/{args.scaffold}_A.pdb'))[1]['A']
 template_design={p for lo,hi in ti['cdr_ranges_1based'] for p in range(lo,hi+1)}
 template_framework=''.join(x for i,x in enumerate(scaffold_seq,1) if i not in template_design)
 default=(Path(__file__).resolve().parents[1]/'outputs'/args.target)/('iggm_fixed' if args.flow=='iggm' else args.flow)/args.scaffold
 if args.flow=='boltzgen':default=default/'intermediate_designs'
 sources=args.input_dir or [default]
 pattern={'boltzgen':'*.cif','ppiflow':'*.pdb','rfantibody':'design_*.pdb','iggm':'*_CDR_All_*.pdb'}[args.flow]
 paths=[]
 for src in sources:
  if not src.is_dir():raise ValueError('Generator output directory is missing: '+str(src))
  batch=sorted(src.glob(pattern),key=lambda p:[int(s) if s.isdigit() else s for s in re.split(r'(\d+)',p.name)])
  if not batch:raise ValueError('No generator structures matched: '+str(src/pattern))
  paths.extend(batch)
 if len({p.resolve() for p in paths})!=len(paths):raise ValueError('Repeated source structures')
 if args.max_structures:
  if args.max_structures<1:raise ValueError('max-structures must be positive')
  paths=paths[:args.max_structures]
 dest=args.output_dir or (Path(__file__).resolve().parents[1]/'outputs'/args.target/'mpnn/input')/args.flow/args.scaffold
 if dest.exists() and any(dest.iterdir()):raise ValueError('Refusing to overwrite existing inputs: '+str(dest))
 fasta=(root/'iggm')/(args.scaffold+'.fasta')
 parsed=[];chains={};fixed={};manifest=[];normalized=[]
 for k,path in enumerate(paths,1):
  atoms=a.atoms_from(path);source_chains,seqs=a.chains_and_sequence(atoms)
  vhh=a.detect_vhh_chain(args.flow,path,path.parent,atoms,fasta)
  targets=[c for c in source_chains if c!=vhh]
  if len(targets)!=1 or seqs.get(targets[0])!=target_seq:raise ValueError('Target sequence mismatch: '+str(path))
  target=targets[0];seq=seqs[vhh]
  cdrs=a.source_cdrs(args.flow,path,path.parent,atoms,vhh,fasta)
  design={p for ps in cdrs for p in ps}
  if len(design)!=sum(map(len,cdrs)) or any(p<1 or p>len(seq) for p in design):raise ValueError('Invalid CDR positions')
  if a.split_three_runs(sorted(design),path,'CDR')!=cdrs:raise ValueError('CDR order/contiguity mismatch')
  if args.flow=='iggm':
   lengths=ti['iggm_cdr_lengths']
   if list(map(len,cdrs))!=lengths:raise ValueError('IgGM fixed CDR length mismatch')
  elif any(not lo<=len(ps)<=hi for ps,(lo,hi) in zip(cdrs,ti['variable_loop_length_ranges'])):raise ValueError('CDR length outside configured range')
  if ''.join(x for i,x in enumerate(seq,1) if i not in design)!=template_framework:raise ValueError('Fixed framework sequence changed: '+str(path))
  order=a.residue_order(atoms);idx={};counts={vhh:0,target:0}
  byres={}
  for chain,rid in order:
   counts[chain]+=1;idx[(chain,rid)]=counts[chain]
  for atom in atoms:
   key=(atom.chain,atom.residue_id)
   if key not in idx:raise ValueError('Residue without CA')
   if atom.atom in byres.setdefault(key,{}):raise ValueError('Duplicate atom')
   if not all(math.isfinite(v) for v in (atom.x,atom.y,atom.z)):raise ValueError('Non-finite coordinates')
   byres[key][atom.atom]=atom
  if any(not {'N','CA','C','O'}<=set(ats) for ats in byres.values()):raise ValueError('Incomplete N/CA/C/O; repair IgGM backbone first')
  name=f'backbone_{k:04d}'
  rec={'name':name,'num_of_chains':2,'seq':seq+target_seq}
  for old,new in [(vhh,'A'),(target,'B')]:
   chainorder=[key for key in order if key[0]==old]
   rec['seq_chain_'+new]=seqs[old]
   rec['coords_chain_'+new]={at+'_chain_'+new:[[getattr(byres[key][at],q) for q in ('x','y','z')] for key in chainorder] for at in ['N','CA','C','O']}
  parsed.append(rec);chains[name]=[['A'],['B']]
  fixed[name]={'A':[i for i in range(1,len(seq)+1) if i not in design],'B':list(range(1,len(target_seq)+1))}
  norm=[dataclasses.replace(x,residue_id=str(idx[(x.chain,x.residue_id)])) for ch in (vhh,target) for x in atoms if x.chain==ch]
  normalized.append((name,norm,{vhh:'A',target:'B'}))
  manifest.append({'name':name,'source':str(path.resolve()),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'flow':args.flow,'target':args.target,'scaffold':args.scaffold,'cdr_positions_1based':cdrs})
 dest.mkdir(parents=True);(dest/'pdbs').mkdir()
 for name,norm,mapping in normalized:a.write_pdb(dest/'pdbs'/f'{name}.pdb',norm,mapping)
 (dest/'parsed.jsonl').write_text(''.join(json.dumps(row,allow_nan=False)+'\n' for row in parsed))
 for file,data in [('chains.jsonl',chains),('fixed_positions.jsonl',fixed),('manifest.json',manifest)]:
  (dest/file).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n' if file.endswith('.json') else json.dumps(data,allow_nan=False)+'\n')
 (dest/'validation.json').write_text(json.dumps({'status':'PASS','backbones':len(parsed),'chain_A':'VHH','chain_B':args.target,'per_backbone_CDR_metadata':True,'fixed_framework_sequence_verified':True,'target_sequence_verified':True,'N_CA_C_O_complete':True},indent=2))
 print('Prepared',len(parsed),'backbones in',dest)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--flow',choices=['boltzgen','ppiflow','rfantibody','iggm'],required=True)
 p.add_argument('--scaffold',required=True);p.add_argument('--target',required=True)
 p.add_argument('--config-dir',type=Path,required=True,help='Target configuration directory containing validation.json, common/ and iggm/')
 p.add_argument('--input-dir',type=Path,action='append',help='Repeat in shard order for independently generated shards');p.add_argument('--output-dir',type=Path);p.add_argument('--max-structures',type=int)
 args=p.parse_args()
 prepare(args)
