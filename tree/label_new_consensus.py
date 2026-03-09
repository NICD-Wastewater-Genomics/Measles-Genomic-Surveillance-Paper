import pandas as pd
import sys
from Bio import SeqIO
import numpy as np

# add date information to new sequences from South Africa
dates_fn = sys.argv[1]
print("Getting date info from:",dates_fn)
index_col = sys.argv[4] #name of index column
date_col = sys.argv[5] #name of date column
if len(sys.argv)>=7:
    with open(sys.argv[6], 'r') as file:
        drops = file.readlines()
    drops = [d.split('\t')[0].split('|')[0].strip() for d in drops[1:]]
    print('seqs to drop: ',drops)
else:
    drops = []
if dates_fn.endswith('tsv'):
    meta = pd.read_csv(dates_fn,index_col=index_col,sep='\t')
elif dates_fn.endswith('csv'):
    meta = pd.read_csv(dates_fn,index_col=index_col)
else:
    print('Unknown file format for metadata')
    sys.exit()

meta.index = [mi.replace('-','_').split('_S')[0] if isinstance(mi,str) else mi for mi in meta.index]
if "OtherID (e.g. Sipho's trial IDs)" in meta.columns: #special case for samples sequenced multiple times with different concentration methods
    meta.index = [ind if other.startswith('Not') else other.replace('-','_') for ind,other in zip(meta.index,meta["OtherID (e.g. Sipho's trial IDs)"])]
# asf
#all RSA for this paper. 
meta['Country'] = 'SouthAfrica' 
original_file = sys.argv[2]
corrected_file = sys.argv[3]
if corrected_file == original_file:
    print('Need to specify different name for date-containing file')
    exit()
j=0
k=0
filter_only = False
record_info = []
print(meta)
with open(original_file) as original, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        k+=1
        if record.id.startswith('Consensus'):
            record.id = record.id.split('Consensus_')[1]
        if '_S' in record.id:
            record.id = record.id.split('_S')[0]
        if '_threshold' in record.id:
            record.id = record.id.split('_threshold')[0]
            
        if record.id in drops:
            print("dropping:", record.id)
            # skip sequences with likely errors
            continue
        if '-' in record.id:
            record.id = record.id.replace('-','_')

        og = record.id
        record.id = '_'.join(record.id.split('_')[0:5])  if record.id.endswith('_C') or record.id.endswith('_D') else '_'.join(record.id.split('_')[0:4])

        if record.id in meta.index and (isinstance(meta.loc[record.id,date_col], str)):
            if not filter_only:
                orig_id = record.id
                if "X" not in meta.loc[record.id,date_col]:
                    meta.loc[record.id,date_col] = pd.to_datetime(meta.loc[record.id,date_col]).strftime('%Y-%m-%d')
                record.id = record.id +'|' + meta.loc[record.id,'Country'].replace(' ','_')+'|' + meta.loc[record.id,date_col]
            record_info.append([record.id,meta.loc[orig_id,date_col]])
        else:
            #skip if not in current grouping
            print('skipping', record.id, og)
            # print(record.id,meta.loc[record.id,'Country'], meta.loc[record.id,date_col])
            continue
        record.description = ''
        # record.id = record.id.replace('_threshold_0.8_quality_20','').replace('Consensus_','')
        j+=1
        SeqIO.write(record, corrected, 'fasta')
print(f'Including {j} of {k} records.')

new_meta = pd.DataFrame(data = record_info,columns=['strain','date'])
new_meta.to_csv(corrected_file.replace('.fasta','_meta.tsv'),sep='\t')

