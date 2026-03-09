import pandas as pd
import pickle
import json
import requests
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.dates as mdates
from scipy import signal
from datetime import date,timedelta
import yaml
import copy
import numpy as np 
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point, Polygon
from matplotlib_scalebar.scalebar import ScaleBar
from adjustText import adjust_text
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['savefig.bbox'] = 'tight'

meta_ww = pd.read_csv('../metadata/Measles_seqdata_02022026.csv')
meta_ww.iloc[:,2] = meta_ww.iloc[:,2].astype(str)
meta_ww.index = [meta_ww['SampleID'].iloc[j] if meta_ww.iloc[j,2].startswith('Not') else meta_ww.iloc[j,2] for j in range(meta_ww.shape[0])]
meta_ww = meta_ww[[mc for mc in meta_ww.columns if 'Unnamed' not in mc]]
meta_ww.index = [mi.replace('_','-').split('-S')[0] for mi in meta_ww.index]
meta_ww = meta_ww[~meta_ww['SampleID'].str.startswith('AIR')]
### get number of unique positives
ww_pos = meta_ww[(meta_ww['MeaslesResult']=='Positive')]
ww_pos['date']= pd.to_datetime(ww_pos['SampleCollectionDate'])
ww_pos = ww_pos[ww_pos['date']>='2024-10-01']
ww_pos = ww_pos[(~ww_pos['SampleID'].duplicated())]
print('Number of unique positive WW samples: ',  ww_pos.shape[0])
ww_pos['Sentforseq (Y/N)'] = ww_pos['Sentforseq (Y/N)'].str.strip()
ww_seqd = ww_pos[ww_pos['Sentforseq (Y/N)']=='Yes']
print('number of samples sent for sequencing',ww_seqd.shape[0])
agg_df = pd.read_csv('../agg_demixed.tsv', skipinitialspace=True, sep='\t',index_col=0)
agg_df = agg_df[~agg_df.index.str.contains('PC_')]
from freyja.utils import prepLineageDict, prepSummaryDict
agg_df = prepSummaryDict(agg_df)
agg_df['abundances'] = agg_df['abundances'].astype(str)
agg_df_likes = prepLineageDict(agg_df, thresh=0.0000000001,  mergeLikes=False)
agg_df_likes.index =  [dsi.split('.tsv')[0].split('_')[-1] for dsi in agg_df_likes.index]
agg_df = prepLineageDict(agg_df, thresh=0.0000000001,  mergeLikes=True)

agg_df['region'] = agg_df['linDict'].apply(lambda x:list(x.keys())[0].split('-')[0].split('MEASLES')[1])
agg_df2 = agg_df.copy()
agg_df2 = agg_df2[(agg_df2.index.str.contains('ENV') | agg_df2.index.str.contains('NAT') | agg_df2.index.str.contains('CST') | agg_df2.index.str.contains('ART_MEV') )]
agg_df2.index = [agi.split('_S')[0] for agi in agg_df2.index]
agg_df2['ID'] = ['-'.join(agi.split('_')[0:5])  if agi.endswith('_C') or agi.endswith('_D') else '-'.join(agi.split('_')[0:4]) for agi in agg_df2.index]
agg_df2['num_lineages'] = [len([x for x in agg_df2['abundances'].iloc[j] if x>0.001]) for j in range(agg_df2.shape[0])]
agg_df2['stype'] = ['&'.join([x.split('MEASLESgenome-')[-1] for x in agg_df2['lineages'].iloc[j]]) for j in range(agg_df2.shape[0])]
agg_df2['ID_min'] = ['-'.join(agi.split('-')[0:4])  if agi.endswith('-C') or agi.endswith('-D') else agi for agi in agg_df2['ID']]

idx = agg_df2.groupby('ID_min')['coverage'].idxmax()

agg_df2 = agg_df2.loc[idx]

colors = plt.cm.Set2

color_list = [colors(0),colors(1),colors(0)]
hatches = ['','','///']

# generate WW sample type (pure/mixed) counts, Fig 2b
fig,ax = plt.subplots(figsize = (2,4))
sns.countplot(x='stype',data=agg_df2,ax=ax,palette=color_list,edgecolor=colors(1),linewidth=0)

for i, bar in enumerate(ax.patches):
    if len(hatches[i])>0:
        bar.set_hatch(hatches[i])

ax.set_xlabel('Genotype present')
ax.set_ylabel('Count')
ax.spines[['right', 'top']].set_visible(False)
plt.savefig('../figures/mixture_genotypes_counts.pdf',bbox_inches='tight')
plt.close('all')

fig,ax = plt.subplots(figsize = (2,4))
sns.countplot(x='num_lineages',data=agg_df2,ax=ax)
ax.set_xlabel('Genotype present')
ax.set_ylabel('Count')
ax.spines[['right', 'top']].set_visible(False)
plt.savefig('../figures/number_of_genotypes_present.pdf',bbox_inches='tight')
plt.close('all')

# check to make sure nothing is missing from the metadata
print('Missing from metadata!: ',agg_df2[~agg_df2['ID'].isin(meta_ww.index)]['ID'].values)

agg_df2 = agg_df2[agg_df2['ID'].isin(meta_ww.index)]#drop samples from outside of Gauteng for now. 
agg_df2['site'] = meta_ww.loc[agg_df2['ID'].values,'SiteName'].values
agg_df2['Province'] = meta_ww.loc[agg_df2['ID'].values,'SiteProvince'].values
agg_df2['date'] = meta_ww.loc[agg_df2['ID'].values,'SampleCollectionDate'].values

# now pull in clinical data
clin_df = pd.read_csv('../metadata/mev clinical metadata_updated 040226 ks_JL_090226.csv',sep=",")

clin_df = clin_df[~(clin_df['Sample ID'].str.startswith('5') | clin_df['Sample ID'].str.startswith('UNK'))]
clin_df['Onset date'] = clin_df['Onset date'].fillna("")
clin_df['Onset date'] = clin_df['Onset date'].apply(lambda x: x if (len(x.split('/')[0])>2) else "/".join(x.split('/')[::-1]))
clin_df['Onset date'] = clin_df['Onset date'].apply(lambda x: x if len(x)>0 else np.nan)
clin_df['Onset date'] = pd.to_datetime(clin_df['Onset date'])

clin_df = clin_df[['Sample ID','Seq ID','Onset date','Province','City/district','WHO naming','WGS','Genotype']]

clin_df = clin_df.set_index('Onset date')
clin_df = clin_df.groupby(pd.Grouper(freq='MS'))['Genotype'].value_counts()
clin_df = clin_df.reset_index()
clin_df = clin_df.pivot(index='Onset date',columns='Genotype',values='count')
clin_df = clin_df[clin_df.index>=pd.to_datetime('2024-10-01')]

clin_gt_totals = clin_df.groupby(pd.Grouper(freq="QE")).sum().fillna(0)
print('2025 sample counts: ',clin_gt_totals.iloc[1:].sum())
clin_totals = clin_df.sum(axis=1)
clin_df = clin_df.div(clin_df.sum(axis=1),axis=0)

types = ['genome']
for t0 in types:
    #convert to abundances matrix (need to make into freyja function)
    queryTypes = ['linDict']
    for queryType in queryTypes:
        df_abundances = pd.DataFrame()
        for i, sampLabel in enumerate(agg_df2.index):
            dat = agg_df2.loc[sampLabel, queryType]
            if isinstance(dat, list):
                df_abundances = pd.concat([
                    df_abundances,
                    pd.Series(
                        agg_df2.loc[sampLabel, queryType][0],
                        name=meta_ww.loc[agg_df2.loc[sampLabel,'ID'],'SampleCollectionDate'])
                ], axis=1)
            else:
                df_abundances = pd.concat([
                    df_abundances,
                    pd.Series(
                        agg_df2.loc[sampLabel, queryType],
                        name=meta_ww.loc[agg_df2.loc[sampLabel,'ID'],'SampleCollectionDate'])
                ], axis=1)
        df_abundances = df_abundances.T.fillna(0.)

df_abundances.index = pd.to_datetime(df_abundances.index)
df_ = df_abundances.groupby(pd.Grouper(freq='MS')).mean()

ww_gt_totals = df_abundances.round().groupby(pd.Grouper(freq="QE")).sum().fillna(0)
ww_totals = df_abundances.groupby(pd.Grouper(freq='MS'))['MEASLESgenome-B3'].count()
print('Unique WW samples successfully sequenced: ',ww_totals.sum())

all_inds = ww_totals.index.union(clin_totals.index)
ww_totals = ww_totals.reindex(all_inds).asfreq("MS").fillna(0)
clin_df = clin_df.reindex(all_inds).asfreq("MS") 
df_ = df_.reindex(all_inds).asfreq("MS") 
clin_totals = clin_totals.reindex(all_inds).asfreq("MS").fillna(0) 
colors = plt.cm.Set2

## let's plot by quarter (Fig 2c)
clin_df = clin_df.groupby(pd.Grouper(freq="QE")).mean().fillna(0)
clin_totals = clin_totals.groupby(pd.Grouper(freq="QE")).sum().fillna(0)
df_ = df_.groupby(pd.Grouper(freq="QE")).mean().fillna(0)
ww_totals = ww_totals.groupby(pd.Grouper(freq="QE")).sum().fillna(0)

fig,ax = plt.subplots(figsize=(5,5))
from matplotlib.colors import colorConverter

width=30
for i in range(0, df_.shape[1]):
    bar0 = ax.bar(df_.index + timedelta(days=width/2+0.5) , df_.iloc[:, i],
           width=width, bottom=df_.iloc[:, 0:i].sum(axis=1),
           label=df_.columns[i], color=colors(i), hatch='..')
    #make hatch lines white
    for bc in bar0:
        bc._hatch_color = colorConverter.to_rgba('white')
        bc.stale = True

for i in range(0, clin_df.shape[1]): 
    ax.bar(clin_df.index-timedelta(days=width/2+0.5), clin_df.iloc[:, i],
        width=width, bottom=clin_df.iloc[:, 0:i].sum(axis=1),
        label=clin_df.columns[i], color=colors(i))

ax.set_ylim(0,1)
locator = mdates.MonthLocator(bymonthday=1)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
ax.set_ylabel('Genotype Prevalence',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
ax.set_xticks(clin_totals.index)
labels = [f"Q{d.quarter}\n{d.year}" for d in clin_totals.index]
ax.set_xticklabels([labels[j0] + f"\n({int(clin_totals.iloc[j0])},{int(ww_totals.iloc[j0])})" for j0,label in enumerate(ax.get_xticklabels())])

fig.tight_layout()
plt.savefig('../figures/MeV_stackplot_Quarterly.pdf')

## run CMH test comparing WW and Clincal genotype prevalence estimates over time. 
from statsmodels.stats.contingency_tables import StratifiedTable
tables = [ [ww_gt_totals.loc[ind,['MEASLESgenome-B3','MEASLESgenome-D8']],clin_gt_totals.loc[ind,['B3','D8']]] for ind in ww_gt_totals.index]
cmh = StratifiedTable(tables)
print(cmh.test_null_odds())