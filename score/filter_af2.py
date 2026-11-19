#!/usr/bin/env python3
import argparse,math
from table_io import read_table,number,write_selection
p=argparse.ArgumentParser();p.add_argument('--scores',required=True);p.add_argument('--output',required=True);p.add_argument('--annotated',required=True);args=p.parse_args()
rows=read_table(args.scores);flags=[]
for row in rows:
 if row.get('af2_score_available','').lower()=='false' or ('af2_pae_interaction' in row and not row['af2_pae_interaction']):
  flags.append((None,'AF2 score unavailable; not evaluated'));continue
 pae=number(row,'pae_interaction')
 if not math.isfinite(pae):raise ValueError('Non-finite AF2 PAE: input integrity error')
 for key,val in row.items():
  # Interface pLDDT NaN is allowed. Other supplied confidence values must be finite.
  if key.startswith(('plddt','pae_','af2_plddt','af2_pae_')):
   num=float(val)
   if math.isnan(num) and key in ['plddt_interface','af2_plddt_interface']:continue
   if not math.isfinite(num):raise ValueError('Non-finite AF2 score: '+key)
 flags.append((pae<10,'pass' if pae<10 else 'pae_interaction >= 10'))
write_selection(args,rows,flags,'AF2',{'pae_interaction':'< 10','plddt_interface_NaN':'allowed; does not count as prediction failure'})
