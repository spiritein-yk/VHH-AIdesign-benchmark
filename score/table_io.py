"""Read strict score tables; preserve source identities and all columns."""
import csv,gzip,json,math
from pathlib import Path
ALIASES={'source_id':['source_id','source_tag','source_design_tag','sequence_id','design_id','id','description','tag'], 'flow':['flow','structure_model','structure_method'], 'mpnn':['mpnn','mpnn_method','sequence_model','sequence_method'], 'pae_interaction':['pae_interaction','af2_pae_interaction'], 'ddg':['ddg','rosetta_ddg'], 'iptm':['iptm','af3_iptm'], 'rmsd':['rmsd_vhh_after_target_fit_A','vhh_rmsd_after_target_fit_A','rmsd_vhh_after_target_fit','rmsd_vhh_after_target_fit_angstrom']}
ALIASES.update({
 'ddg_per_sasa':['rosetta_ddg_per_sasa','ddg_per_sasa'],
 'buns_heavy_ball_1.1D':['rosetta_buns_heavy_ball_1.1D','buns_heavy_ball_1.1D'],
 'contact_molecular_surface':['rosetta_contact_molecular_surface','contact_molecular_surface'],
 'rmsd':['af3_rmsd','rmsd_vhh_after_target_fit_A','vhh_rmsd_after_target_fit_A','rmsd_vhh_after_target_fit','rmsd_vhh_after_target_fit_angstrom','rmsd']})
def read_table(path):
 path=Path(path)
 with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:lines=[s for s in f if s.strip()]
 if not lines:raise ValueError('Empty score table')
 if any(s.startswith('SCORE:') for s in lines):
  header=None;rows=[]
  for s in lines:
   if not s.startswith('SCORE:'):continue
   parts=s.split()[1:]
   if parts and parts[-1]=='description':
    if header and header!=parts:raise ValueError('Inconsistent score headers')
    header=parts;continue
   if header is None or len(parts)!=len(header):raise ValueError('Invalid score row')
   rows.append(dict(zip(header,parts)))
 else:
  delimiter='\t' if '\t' in lines[0] else ','
  reader=csv.DictReader(lines,delimiter=delimiter)
  if len(reader.fieldnames)!=len(set(reader.fieldnames)):raise ValueError('Duplicate columns')
  rows=list(reader)
  if any(None in r or None in r.values() for r in rows):raise ValueError('Malformed table row')
 seen=set()
 for r in rows:
  key=value(r,'source_id')
  if key in seen:raise ValueError('Duplicate source ID: '+key)
  seen.add(key)
 return rows

def value(row,key):
 
 if key=='source_id':
  for canonical in ('design_id','source_id'):
   if row.get(canonical):return row[canonical]
 names=ALIASES.get(key,[key]);found=[row[n] for n in names if n in row and row[n]!='']
 if not found:raise ValueError('Missing field: '+key)
 if len(set(found))>1:raise ValueError('Conflicting aliases: '+key)
 return found[0]

def number(row,key):return float(value(row,key))
def finite(row,key):return math.isfinite(number(row,key))
def group(row):
 flow=value(row,'flow').lower().replace('-','').replace('_','')
 mpnn=value(row,'mpnn').lower().replace('-','').replace('_','')
 flow={'iggmfixed':'iggm','rfantibody':'rfantibody','boltzgen':'boltzgen','ppiflow':'ppiflow','iggm':'iggm'}.get(flow)
 mpnn={'proteinmpnn':'proteinmpnn','solublempnn':'solublempnn','antibmpnn':'antibmpnn'}.get(mpnn)
 if flow is None or mpnn is None:raise ValueError('Unknown flow/MPNN combination')
 return flow,mpnn

def write_table(path,rows,fields):
 path=Path(path)
 if path.exists():raise ValueError('Refusing to overwrite '+str(path))
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)

def write_selection(args,rows,flags,stage,rules):
 fields=list(rows[0]) if rows else ['source_id']
 new=['selection_stage','selection_pass','selection_reason']
 if set(fields)&set(new):raise ValueError('Selection fields already exist')
 annotated=[dict(r,selection_stage=stage,selection_pass=("" if ok is None else int(ok)),selection_reason=reason) for r,(ok,reason) in zip(rows,flags)]
 outputs=[Path(args.output),Path(args.annotated),Path(args.output+'.summary.json')]
 if any(p.exists() for p in outputs):raise ValueError('Output already exists')
 write_table(args.output,[r for r in annotated if r['selection_pass']],fields+new)
 write_table(args.annotated,annotated,fields+new)
 outputs[2].write_text(json.dumps({'stage':stage,'input_count':len(rows),'passed_count':sum(ok is True for ok,_ in flags),'rules':rules},indent=2)+'\n')
 print(stage,'passed',sum(ok is True for ok,_ in flags),'/',len(rows))
