import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import matplotlib
import statsmodels.api as sm
from statsmodels.formula.api import ols

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['savefig.bbox'] = 'tight'

dir0 = '../depths/'
covcut=10
## Figures produced here: 
# SFig 1b and SFig 3b

# read in depth files for coverage analyses. 
covs = []
regions = []
samps = []
for f in os.listdir(dir0):
    print(f)
    df_depth = pd.read_csv(dir0 +f, sep='\t', header=None, index_col=1)[[3]]
    cov = sum(df_depth.loc[:, 3] > covcut)/float(df_depth.shape[0])
    region = 'Whole Genome'
    covs.append(cov*100.)
    regions.append(region)
    samps.append(f.split('_S')[0])

agg_df2 = pd.DataFrame({'sname':samps,'region':regions,'coverage':covs})
agg_df2 = agg_df2[agg_df2['region']!='N450']
agg_df2.to_csv('all_coverages.csv')
agg_df2 = agg_df2[agg_df2['sname'].str.contains('ENV') | agg_df2['sname'].str.contains('NAT') | agg_df2['sname'].str.contains('CST') | agg_df2['sname'].str.contains('ART_MEV') ]


meta_ww = pd.read_csv('../metadata/Measles_seqdata_02022026.csv')
meta_ww.iloc[:,2] = meta_ww.iloc[:,2].astype(str)
meta_ww.index = [meta_ww['SampleID'].iloc[j] if meta_ww.iloc[j,2].startswith('Not') else meta_ww.iloc[j,2] for j in range(meta_ww.shape[0])] #,index_col='SampleID')#pd.read_csv('../metadata/MeV Wastewater Sequences Final Metadata_05092025.csv',index_col='Sample Number')

meta_ww = meta_ww[[mc for mc in meta_ww.columns if 'Unnamed' not in mc]]
meta_ww.index = [mi.replace('_','-').split('-S')[0] for mi in meta_ww.index]

agg_df2['ID'] = ['-'.join(agi.split('_')[0:5])  if agi.endswith('_C') or agi.endswith('_D') else '-'.join(agi.split('_')[0:4]) for agi in agg_df2.sname]
meta_ww['Concentration method'] = meta_ww['ConcMethod'].astype(str).apply(lambda x:x.replace(' ','')).apply(lambda x:x.replace('CeresNanotrapmicrobiomeParticles&Dynabeads','2-bead'))
meta_ww['Concentration method'] = [m if m!='nan' else 'Dynabeads' for m in meta_ww['Concentration method']]
meta_ww.to_csv('cleaned_ww_metadata.csv')
meta_ww['Measles Concentration (copies/uL)'] = meta_ww['MeaslesConc']
agg_df2['Measles Concentration (copies/uL)'] = [meta_ww.loc[id,'Measles Concentration (copies/uL)'] if id in meta_ww.index else None for id in agg_df2['ID']]
agg_df2['log concentration'] = np.log10(agg_df2['Measles Concentration (copies/uL)'])
agg_df2['coverage'] = agg_df2['coverage']/100.
agg_df2['Concentration method'] = [meta_ww.loc[id,'Concentration method'] if id in meta_ww.index else None for j,id in enumerate(agg_df2['ID'])]

agg_df_comparison = agg_df2.copy()

## Generate SFig 3B
agg_df2 = agg_df2[(agg_df2['Concentration method']=='Dynabeads')]# | (agg_df2['Concentration method']=='Ceres+Dynabeads')]
agg_df2 = agg_df2.dropna()
for cut in [0.05,0.10,0.25,0.5,0.75]:
    passing = agg_df2[agg_df2['coverage']>cut]
    print(f'over threshold{cut}:', passing.shape[0], ' of ', agg_df2.shape[0])


xlim=[-0.5,3.1]
fig, ax = plt.subplots(figsize=(3, 4))
ax.set_xlim(xlim)
ax.set_ylim([0, 1])

# Define bins in log10 space
bin_edges = [-np.inf, 0, 1, np.inf]
log_conc = agg_df2['log concentration'].values
coverage = agg_df2['coverage'].values

# Add alternating bin shading first (so it sits behind points)
for i, (lo, hi) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):
    x_lo = max(lo, xlim[0]) if lo != -np.inf else xlim[0]
    x_hi = min(hi, xlim[1]) if hi != np.inf else xlim[1]
    if i % 2 == 0:
        ax.axvspan(x_lo, x_hi, color='gray', alpha=0.08, zorder=0)

ax.scatter(agg_df2['log concentration'], agg_df2["coverage"], clip_on=False, color=plt.cm.Accent(1), zorder=2)

for i, (lo, hi) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):
    mask = (log_conc > lo) & (log_conc <= hi)
    if mask.sum() == 0:
        continue
    median_val = np.median(coverage[mask])
    mad_val = np.median(np.abs(coverage[mask] - median_val))

    x_lo = max(lo, xlim[0]) if lo != -np.inf else xlim[0]
    x_hi = min(hi, xlim[1]) if hi != np.inf else xlim[1]
    x_mid = (x_lo + x_hi) / 2

    ax.plot([x_lo, x_hi], [median_val, median_val], color='black', linewidth=1.5, zorder=3)
    ax.errorbar(x_mid, median_val, yerr=mad_val, fmt='none', color='black', capsize=3, linewidth=1.2, zorder=3)

ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_xlabel('Log10(Measles Concentration (copies/uL))')
ax.set_ylabel('Fraction genome coverage')
log_labels = [f'$10^{label.get_text()}$' for label in ax.get_xticklabels()]
ax.set_xticklabels(log_labels)
ax.spines[['right', 'top']].set_visible(False)
plt.savefig('../figures/coverage_vs_concentration.pdf', bbox_inches='tight', transparent=True)
plt.close()

# build SFig 1B
fig,ax = plt.subplots(figsize=(3.5,4))
sns.swarmplot( x='Concentration method', y='log concentration', data=agg_df_comparison, ax=ax,color='black',clip_on=False)
sns.boxplot(data=agg_df_comparison, x='Concentration method', y="log concentration",ax=ax,showfliers=False,color='lightgrey',widths=0.6)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_ylabel('Log10(Measles Concentration (copies/uL))')
ax.set_ylim([-0.75,3.25])
ax.spines[['right', 'top']].set_visible(False)
plt.savefig('../figures/conc_method_vs_conc.pdf',bbox_inches='tight',transparent=True)
plt.close()
