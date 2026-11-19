#!/usr/bin/env python3
import argparse,math
from table_io import read_table,number,value,group,write_selection
p=argparse.ArgumentParser();p.add_argument('--scores',required=True);p.add_argument('--output',required=True);p.add_argument('--annotated',required=True);p.add_argument('--target',choices=['IL7RA','PDGFRB'],required=True);args=p.parse_args()
rows=read_table(args.scores);eligible=[]
for row in rows:
 x=number(row,'pae_interaction');flow,mpnn=group(row)
 if not math.isfinite(x):raise ValueError('Non-finite AF2 PAE')
 eligible.append(x<10 if args.target=='PDGFRB' else x<8 or (x<10 and (flow=='iggm' or mpnn=='antibmpnn' or (flow=='ppiflow' and mpnn=='proteinmpnn'))))
chosen={i for i,ok in enumerate(eligible) if ok}
if args.target=='PDGFRB':
 groups={}
 for i in chosen:groups.setdefault(group(rows[i]),[]).append(i)
 chosen={i for indices in groups.values() for i in sorted(indices,key=lambda j:(number(rows[j],'pae_interaction'),value(rows[j],'source_id')))[:500]}
flags=[(i in chosen,'pass' if i in chosen else 'outside AF3 eligibility or group top500') for i in range(len(rows))]
rules={'IL7RA':'PAE < 8 OR [PAE < 10 AND (IgGM OR AntiB-MPNN OR PPIFlow × ProteinMPNN)]; no new topK truncation','PDGFRB':'PAE < 10 then ascending PAE top500 per flow × MPNN across templates; ties by source ID; all if fewer than500'}
write_selection(args,rows,flags,'AF3 admission',{'target':args.target,'rule':rules[args.target],'sequence_deduplication':False,'tie_break':'source ID for this reusable version'})
