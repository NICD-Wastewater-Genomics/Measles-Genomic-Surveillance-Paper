import argparse
import os
import time
import pandas as pd
from datetime import datetime
## script to extract clades from metadata where possible and correct entries that are wrong on genbank
meta_fn = 'measles_bg_metadata.csv'
meta = pd.read_csv(meta_fn,index_col=0)
meta = meta[['organism','isolate','note','strain','host','geo_loc_name','collection_date_iso','collection_date']]
if len(meta[meta['host']!='Homo sapiens'])>0:
    print('Some non-human host specimens',meta[meta['host']!='Homo sapiens'])
meta = meta[meta['host']=='Homo sapiens']
genotypes = []
for iso in meta.index:
    if 'genotype' in meta.loc[iso,'organism']:
        genotypes.append(meta.loc[iso,'organism'].split(' ')[-1])
    elif '[' in meta.loc[iso,'organism']:
        genotypes.append(meta.loc[iso,'organism'].split('[')[-1].split(']')[0])
    elif 'genotype:' in meta.loc[iso,'note']:
        genotypes.append(meta.loc[iso,'note'].split('genotype:')[-1].strip())
    else:
        genotypes.append(None)
meta['clade'] = genotypes
corrections = [['MZ712082.1','B3'],['MZ031229.1','B3'],['MZ031227.1','D8']]
for c in corrections:
    meta.loc[c[0],'clade'] = c[1]
print('Missing clade information for:', meta[meta['clade'].isnull()])
meta['country'] = meta['geo_loc_name'].apply(lambda x: x.split(':')[0]  if (isinstance(x,str)) else None)
meta['locality'] = meta['geo_loc_name'].apply(lambda x: None if (not isinstance(x,str)) or (':' not in x) else x.split(':')[1].strip())
meta.to_csv('measles_bg_metadata_cleaned.csv')
print(meta['geo_loc_name'].unique())
meta = meta.dropna(subset=['clade'])
meta = meta[['clade']].reset_index()
meta = meta.iloc[:,::-1]
meta.to_csv('background/measles_clades_cleaned.tsv',sep='\t',index=False)