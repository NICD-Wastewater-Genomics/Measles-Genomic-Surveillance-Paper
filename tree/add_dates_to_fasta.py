import pandas as pd
import sys
from Bio import SeqIO
import numpy as np
## script for adding date information to fasta files
dates_fn = sys.argv[1]
print("Getting date info from:",dates_fn)
date_col = sys.argv[4] #name of date column
if len(sys.argv)>=6:
    with open(sys.argv[5], 'r') as file:
        drops = file.readlines()
    drops = [d.strip() for d in drops]
    print('Dropping Files: ',drops)
else:
    drops = []
# filter_genotype = sys.argv[5] #name of date column
if dates_fn.endswith('tsv'):
    meta = pd.read_csv(dates_fn,index_col=0,sep='\t')
elif dates_fn.endswith('csv'):
    meta = pd.read_csv(dates_fn,index_col=0)
else:
    print('Unknown file format for metadata')
    sys.exit()
meta['country'] = meta['country'].fillna('UNKNOWN')
original_file = sys.argv[2]
corrected_file = sys.argv[3]
j=0
k=0
l=0
filter_only = False
record_info = []
with open(original_file) as original, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        k+=1
        if record.id in drops:
            l+=1
            # skip sequences with likely errors
            continue
        if record.id in meta.index and (isinstance(meta.loc[record.id,date_col], str)):
            if not filter_only:
                orig_id = record.id
                record.id = record.id +'|' + meta.loc[record.id,'country'].replace(' ','_')+'|' + meta.loc[record.id,date_col]
            record_info.append([record.id,meta.loc[orig_id,date_col]])
        else:
            #skip if not in current grouping 
            print('skipping', record.id)
            continue
        record.description = ''
        # record.id = record.id.replace('_threshold_0.8_quality_20','').replace('Consensus_','')
        j+=1
        SeqIO.write(record, corrected, 'fasta')
print(f'Including {j} of {k-l} records.')

new_meta = pd.DataFrame(data = record_info,columns=['strain','date'])
new_meta.to_csv(corrected_file.replace('.fasta','_meta.tsv'),sep='\t')

