#!/usr/bin/env python3
import argparse,math
from table_io import read_table,number,write_selection
p=argparse.ArgumentParser();p.add_argument('--scores',required=True);p.add_argument('--output',required=True);p.add_argument('--annotated',required=True);args=p.parse_args()
rows=read_table(args.scores);flags=[]
for row in rows:
 vals=[number(row,k) for k in ['ddg','ddg_per_sasa','buns_heavy_ball_1.1D','contact_molecular_surface']]
 ok=all(math.isfinite(x) for x in vals) and vals[0]<-20 and vals[1]<=-0.02 and vals[2]<=2 and vals[3]>=200
 flags.append((ok,'pass' if ok else 'Rosetta criteria not met or required score undefined'))
write_selection(args,rows,flags,'Rosetta',{'ddg':'< -20','ddg_per_sasa':'<= -0.02','buns_heavy_ball_1.1D':'<= 2','contact_molecular_surface':'>= 200','logic':'AND'})
