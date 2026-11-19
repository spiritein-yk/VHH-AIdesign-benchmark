#!/usr/bin/env python3
import argparse,math
from table_io import read_table,number,write_selection
p=argparse.ArgumentParser();p.add_argument('--scores',required=True);p.add_argument('--output',required=True);p.add_argument('--annotated',required=True);args=p.parse_args()
rows=read_table(args.scores);flags=[]
for row in rows:
 if 'af3_iptm' in row and (not row['af3_iptm'] or not row.get('af3_rmsd',row.get('rmsd_vhh_after_target_fit_A',''))):
  flags.append((None,'AF3 score unavailable; outcome unknown'));continue
 ddg,rmsd,iptm=[number(row,k) for k in ['ddg','rmsd','iptm']]
 if not all(math.isfinite(x) for x in [ddg,rmsd,iptm]):raise ValueError('Non-finite required AF3/source score')
 if rmsd<0 or not 0<=iptm<=1:raise ValueError('Invalid RMSD or ipTM range')
 ok=ddg<-30 and rmsd<5 and iptm>0.8
 flags.append((ok,'pass' if ok else 'AF3 joint criteria not met'))
write_selection(args,rows,flags,'AF3 final',{'original_design_ddg':'< -30','VHH_CA_RMSD_after_target_fit_A':'< 5','native_top_iptm':'> 0.8','logic':'AND','additional_AF2_PAE_lt8':False})
