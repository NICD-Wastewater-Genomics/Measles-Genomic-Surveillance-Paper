import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['savefig.bbox'] = 'tight'

# Figures produced here:
# SFig 5a,b

dirD = '../depths/'
dirV = '../variants/'
depcut=20
covcut=40
nVs = []
nVsDominant = []
samps = []
regions = []
coverages = []
for f in os.listdir(dirD):
    ## only consider wastewater samples, which can have multiple different naming conventions 
    if ('ENV' not in f) and ('CST' not in f) and ('NAT' not in f) and ('ART_MEV' not in f):
        continue
    df_depth = pd.read_csv(dirD +f, sep='\t', header=None, index_col=1)[[3]]
    cov = 100*sum(df_depth.loc[:, 3] > depcut)/float(df_depth.shape[0])
    if cov<covcut: # require minimum overall genome coverage for diversity calculation
        continue
    df_vars = pd.read_csv(dirV + f,sep='\t')
    df_vars = df_vars[df_vars['TOTAL_DP']>depcut]
    bottom = 0.10
    top = 0.90
    numVars = sum((df_vars['ALT_FREQ']>bottom) &(df_vars['ALT_FREQ']<top))
    numVars_2= sum((df_vars['ALT_FREQ']>top))
    df_vars = df_vars[(df_vars['POS']>=1233) & (df_vars['POS']<=1682)] #n450 coords
    df_depth_N450 = df_depth[(df_depth.index>=1233) & (df_depth.index<=1682)] #n450 coords
    cov_450 = 100*sum(df_depth_N450.loc[:, 3] > covcut)/float(df_depth_N450.shape[0])
    numVars_N450 = sum((df_vars['ALT_FREQ']>bottom) &(df_vars['ALT_FREQ']<top))
    numVars_N450_2= sum((df_vars['ALT_FREQ']>top))
    if numVars_N450>0:
        print('N450 isnvs', df_vars[(df_vars['ALT_FREQ']>bottom) &(df_vars['ALT_FREQ']<top)])
    nVs.extend([numVars,numVars_N450])
    nVsDominant.extend([numVars_2,numVars_N450_2])
    regions.extend(['Whole Genome','N450'])
    samps.extend([f,f])
    coverages.extend([cov,cov_450])

agg_df2 = pd.DataFrame({'sname':samps,'numVars':nVs,'region':regions,'numVarsDominant':nVsDominant,'coverage':coverages})
agg_df2 = agg_df2[~agg_df2['sname'].str.startswith('CST-WGS-25-0402')] #filter out the sample with a known genotype mix

meta_ww = pd.read_csv('../metadata/Measles_seqdata_02022026.csv')
meta_ww.iloc[:,2] = meta_ww.iloc[:,2].astype(str)
meta_ww.index = [meta_ww['SampleID'].iloc[j] if meta_ww.iloc[j,2].startswith('Not') else meta_ww.iloc[j,2] for j in range(meta_ww.shape[0])]
meta_ww = meta_ww[[mc for mc in meta_ww.columns if 'Unnamed' not in mc]]
meta_ww.index = [mi.replace('_','-') for mi in meta_ww.index]
# clean up namings from comparison study of  concentration methods
agg_df2['ID'] = ['-'.join(agi.split('_')[0:5]) if '_C_' in agi or '_D_' in agi else '-'.join(agi.split('_')[0:4]) for agi in agg_df2.sname]

meta_ww['Concentration method'] = meta_ww['ConcMethod'].astype(str).apply(lambda x:x.replace(' ','')).apply(lambda x:x.replace('CeresNanotrapmicrobiomeParticles&Dynabeads','2-bead'))
meta_ww['Concentration method'] = [m if m!='nan' else 'Dynabeads' for m in meta_ww['Concentration method']]
meta_ww.to_csv('cleaned_ww_metadata.csv')
meta_ww['Measles Concentration (copies/uL)'] = meta_ww['MeaslesConc']
agg_df2['Measles Concentration (copies/uL)'] = [meta_ww.loc[id,'Measles Concentration (copies/uL)'] if id in meta_ww.index else None for id in agg_df2['ID']]
agg_df2['log concentration'] = np.log10(agg_df2['Measles Concentration (copies/uL)']+1)
agg_df2['coverage'] = agg_df2['coverage']/100.
agg_df2['Concentration method'] = [meta_ww.loc[id,'Concentration method'] if id in meta_ww.index else None for j,id in enumerate(agg_df2['ID'])]


# build SFig 5b
fig,ax = plt.subplots(figsize = (2,4))
sns.swarmplot(x='region',y='numVars',data=agg_df2,ax=ax,color='black',clip_on=False,size=3.)
sns.boxplot(x='region',y='numVars',data=agg_df2,ax=ax,color='lightgrey',fliersize=0,medianprops=dict(linewidth=2),zorder=-10)
ax.set_ylim([0,agg_df2.numVars.max()+1])
ax.set_xticklabels(['Whole\ngenome','N450'])
ax.set_xlabel('')
ax.set_ylabel('Number of non-dominant SNVs')
ax.spines[['right', 'top']].set_visible(False)
plt.savefig('../figures/diversity_by_region.pdf',bbox_inches='tight')
plt.close('all')

medians = agg_df2.groupby('region')['numVars'].median()
print("non-dominant medians: ",medians)

# build SFig 5a
fig,ax = plt.subplots(figsize = (2,4))
sns.swarmplot(x='region',y='numVarsDominant',data=agg_df2,ax=ax,color='black',clip_on=False,size=3.)
sns.boxplot(x='region',y='numVarsDominant',data=agg_df2,ax=ax,color='lightgrey',fliersize=0,medianprops=dict(linewidth=2),zorder=-10)

ax.set_xlabel('')
ax.set_ylabel('Number of strongly dominant SNVs')
ax.set_ylim([0,agg_df2.numVarsDominant.max()+1])

ax.spines[['right', 'top']].set_visible(False)
plt.savefig('../figures/dominant_diversity_by_region.pdf',bbox_inches='tight')
plt.close('all')

medians = agg_df2.groupby('region')['numVarsDominant'].median()
print("non-dominant medians: ",medians)