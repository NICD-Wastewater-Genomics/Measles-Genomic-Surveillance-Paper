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
fig,ax = plt.subplots(figsize=(3,4))
ax.set_xlim(xlim)
ax.scatter(agg_df2['log concentration'], agg_df2["coverage"],clip_on=False,color=plt.cm.Accent(1))
import statsmodels.api as sm
# build design matrix 
X = sm.add_constant(agg_df2['log concentration'])
y = agg_df2["coverage"].values

# perform initial fit
model = sm.GLM(y, X, family=sm.families.Binomial())
res = model.fit()

# iterative robust reweighting 
weights = np.ones_like(y, dtype=float)

for _ in range(10):
    res_w = sm.GLM(
        y, X,
        family=sm.families.Binomial(),
        var_weights=weights
    ).fit()

    # use pearson residuals for reweighting
    p_hat = res_w.predict(X)
    resid = (y - p_hat) / np.sqrt(p_hat * (1 - p_hat) + 1e-8)
    new_weights = 1.0 / (1.0 + (resid ) ** 2)
    new_weights *= len(new_weights) / new_weights.sum()

    # check for convergence
    if np.max(np.abs(new_weights - weights)) < 1e-3:
        weights = new_weights
        break

    weights = new_weights

x_grid = np.linspace(
    xlim[0],
    xlim[1],
    200
)
Xg = sm.add_constant(x_grid)
p = res_w.predict(Xg)

ax.plot(x_grid, p, linewidth=2, color=plt.cm.Accent(1))
B = 1000          
n_iter = 10
rng = np.random.default_rng(0)

X_mat = X.values if hasattr(X, "values") else np.asarray(X) 
y_vec = y if isinstance(y, np.ndarray) else np.asarray(y)

Xg_mat = Xg.values if hasattr(Xg, "values") else np.asarray(Xg)

p_boot = np.empty((B, len(x_grid)), dtype=float)

n = len(y_vec)
for b in range(B):
    idx = rng.integers(0, n, size=n)  # bootstrap indices
    Xb = X_mat[idx, :]
    yb = y_vec[idx]

    # perform initial fit
    weights_b = np.ones_like(yb, dtype=float)
    res_w_b = sm.GLM(yb, Xb, family=sm.families.Binomial()).fit()

    # robust reweighting
    for _ in range(n_iter):
        res_w_b = sm.GLM(
            yb, Xb,
            family=sm.families.Binomial(),
            var_weights=weights_b
        ).fit()

        p_hat_b = res_w_b.predict(Xb)
        resid_b = (yb - p_hat_b) / np.sqrt(p_hat_b * (1 - p_hat_b) + 1e-8)

        new_w_b = 1.0 / (1.0 + (resid_b ) ** 2)
        new_w_b *= len(new_w_b) / new_w_b.sum()

        if np.max(np.abs(new_w_b - weights_b)) < 1e-3:
            weights_b = new_w_b
            break
        weights_b = new_w_b

    p_boot[b, :] = res_w_b.predict(Xg_mat)

p_lo = np.percentile(p_boot, 2.5, axis=0)
p_hi = np.percentile(p_boot, 97.5, axis=0)

ax.fill_between(x_grid, p_lo, p_hi, alpha=0.2, linewidth=0,color=plt.cm.Accent(1))
ax.set_ylim([0,1])
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_xlabel('Log10(Measles Concentration (copies/uL))')
ax.set_ylabel('Fraction genome coverage')
print('labels',ax.get_xticklabels())
log_labels = [f'$10^{label.get_text()}$' for label in ax.get_xticklabels()]
ax.set_xticklabels(log_labels)
ax.spines[['right', 'top']].set_visible(False)
plt.savefig('../figures/coverage_vs_concentration.pdf',bbox_inches='tight',transparent=True)
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
