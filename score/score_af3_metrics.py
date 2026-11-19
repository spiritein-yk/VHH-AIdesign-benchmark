#!/usr/bin/env python3
"""Score native AF3 top versus original packed design, after target fitting."""
import argparse,json,math
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser,MMCIFParser
from table_io import read_table,value,number,write_table
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
def ca_chains(path):
 parser=MMCIFParser(QUIET=True) if path.suffix.lower() in ['.cif','.mmcif'] else PDBParser(QUIET=True)
 model=next(parser.get_structure('s',str(path)).get_models());out=[]
 for chain in model:
  seq=[];coords=[]
  for r in chain:
   if r.id[0]!=' ':continue
   if r.resname not in AA or 'CA' not in r:raise ValueError('Unknown or incomplete residue: '+str(path))
   atom=r['CA']
   if atom.is_disordered():raise ValueError('Ambiguous alternate CA atom')
   seq.append(AA[r.resname]);coords.append(atom.coord.astype(float))
  if seq:out.append((''.join(seq),np.array(coords)))
 if len(out)!=2:raise ValueError('Expected exactly two protein chains')
 return out

def roles(chains,vhh,target):
 v=[xyz for seq,xyz in chains if seq==vhh];t=[xyz for seq,xyz in chains if seq==target]
 if len(v)!=1 or len(t)!=1 or vhh==target:raise ValueError('Chain sequence identity/coverage mismatch')
 if not np.isfinite(v[0]).all() or not np.isfinite(t[0]).all():raise ValueError('Non-finite CA coordinates')
 return v[0],t[0]

def fitted_rmsd(v,t,v3,t3):
 a=t.mean(0);b=t3.mean(0);u,s,vh=np.linalg.svd((t-a).T@(t3-b))
 if np.linalg.det(u@vh)<0:u[:,-1]*=-1
 rotation=u@vh
 rmsd=lambda x,y:float(np.sqrt(np.mean(np.sum((x-y)**2,axis=1))))
 return rmsd((v-a)@rotation+b,v3),rmsd((t-a)@rotation+b,t3)

def main():
 p=argparse.ArgumentParser();p.add_argument('--sources',required=True);p.add_argument('--output',required=True);args=p.parse_args()
 rows=read_table(args.sources);scored=[]
 new=['af3_iptm','af3_ptm','af3_ranking_score','rmsd_vhh_after_target_fit_A','rmsd_target_fit_A']
 for row in rows:
  if set(row)&set(new):raise ValueError('AF3 metrics already present')
  base=Path(args.sources).resolve().parent
  result=base/value(row,'af3_result_dir');original=base/value(row,'original_design_pdb')
  cifs=list(result.glob('*_model.cif'));summaries=list(result.glob('*_summary_confidences.json'))
  if len(cifs)!=1 or len(summaries)!=1:raise ValueError('Native top model/summary is missing or ambiguous')
  score=json.loads(summaries[0].read_text());iptm=float(score['iptm']);ranking=float(score['ranking_score']);ptm=float(score['ptm'])
  if not all(math.isfinite(x) for x in [iptm,ranking,ptm,number(row,'ddg')]):raise ValueError('Non-finite required score')
  vhh=value(row,'vhh_sequence');target=value(row,'target_sequence')
  v,t=roles(ca_chains(original),vhh,target);v3,t3=roles(ca_chains(cifs[0]),vhh,target)
  vr,tr=fitted_rmsd(v,t,v3,t3)
  scored.append(dict(row,af3_iptm=iptm,af3_ptm=ptm,af3_ranking_score=ranking,rmsd_vhh_after_target_fit_A=vr,rmsd_target_fit_A=tr))
 fields=list(rows[0])+new if rows else ['source_id']+new
 write_table(args.output,scored,fields)
 print('Scored',len(scored),'source designs against AF3 native top models')
if __name__=='__main__':main()
