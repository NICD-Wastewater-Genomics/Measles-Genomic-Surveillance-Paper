import matplotlib.pyplot as plt
from BCBio import GFF
# generate ED Fig 9
# read in info from the ref gff file
gff_file = "../assets/measles_ref.gff"
genes = []
with open(gff_file) as handle:
    for rec in GFF.parse(handle):
        for feature in rec.features:
            if feature.type == "gene":
                name = feature.qualifiers.get("Name", ["unknown"])[0]
                start = int(feature.location.start)
                end = int(feature.location.end)
                genes.append((name, start, end))


import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['font.family'] = 'Arial'

df_barcodes = pd.read_csv('../assets/MEASLESwholegenome_barcodes.csv',index_col=0)
from freyja.convert_paths2barcodes import sortFun
dictMuts = {}
for lin in df_barcodes.index:
    muts = sorted([df_barcodes.columns[m0]
                    for m0, v in enumerate(df_barcodes.loc[lin])
                    if v > 0], key=sortFun)
    dictMuts[lin] = muts

b3muts = dictMuts['MEASLESgenome-B3']
d8muts = dictMuts['MEASLESgenome-D8']

dirD = '../depths/'
dirV = '../variants/'
covcut=10

nVs = []
nVsDominant = []
samps = []
regions = []
coverages = []
minDepth = covcut
df_vars = []
for f in os.listdir(dirD):
    if ('ENV' not in f) and ('CST' not in f) and ('NAT' not in f) and ('ART_MEV' not in f):
        continue
    if 'n450' in f:
        continue
    if '_C_' in f:
        print(f)
        continue
    df_depth = pd.read_csv(dirD +f, sep='\t', header=None, index_col=1)[[3]]
    cov = 100*sum(df_depth.loc[:, 3] > covcut)/float(df_depth.shape[0])
    if cov<40: #require minimum overall genome coverage for diversity calculation
        continue
    df_var = pd.read_csv(dirV + f,sep='\t')
    df_var = df_var[df_var['TOTAL_DP']>=minDepth]
    df_var = df_var[['POS','REF','ALT','ALT_FREQ','TOTAL_DP']]
    df_var['sample'] = f
    df_vars.append(df_var)

df_vars = pd.concat(df_vars)
# remove all frame shift indels
df_vars = df_vars[df_vars['ALT'].apply(lambda x: False if (('+' in x) and (len(x) % 3)!=1) else True)]
df_vars = df_vars[df_vars['ALT'].apply(lambda x: False if (('-' in x) and (len(x) % 3)!=1) else True)]

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
agg_df2 = agg_df2[agg_df2.index.str.contains('ENV') | agg_df2.index.str.contains('NAT') | agg_df2.index.str.contains('CST') | agg_df2.index.str.contains('ART_MEV') ]
agg_df2 = agg_df2[~agg_df2.index.str.contains('_C_')]
agg_df2['ID'] = [agi.split('_S')[0] for agi in agg_df2.index]
agg_df2['num_lineages'] = [len([x for x in agg_df2['abundances'].iloc[j] if x>0.001]) for j in range(agg_df2.shape[0])]
agg_df2['stype'] = ['&'.join([x.split('MEASLESgenome-')[-1] for x in agg_df2['lineages'].iloc[j]]) for j in range(agg_df2.shape[0])]

df_vars['stype'] = agg_df2.loc[df_vars['sample'],'stype'].values
df_vars['num_lineages'] = agg_df2.loc[df_vars['sample'],'num_lineages'].values
df_vars = df_vars[df_vars['num_lineages']==1]
df_vars_ = df_vars.copy()
colors = plt.cm.Set2
colors_list = [colors(0),colors(1)]
def iqr(x):
    q25, q75 = np.percentile(x, [25, 75])
    return q25, q75

for st in ['B3','D8']:
    df_vars = df_vars_.copy()
    df_vars = df_vars[df_vars['stype']==st]
    df_vars = df_vars[(df_vars['ALT_FREQ']>0.1) & (df_vars['ALT_FREQ']<0.9)]
    df_vars = df_vars.groupby(['POS','REF','ALT'])['ALT_FREQ'].agg(median='median',iqr=iqr,count='count').reset_index()
    df_vars['q25'] = df_vars['iqr'].apply(lambda x:x[0])
    df_vars['q75'] = df_vars['iqr'].apply(lambda x:x[1])
    df_vars = df_vars[df_vars['count']>1]
    # df_vars = df_vars[(df_vars['median']>0.1) & (df_vars['median']<0.9)]
    df_vars['mutName'] = df_vars['REF'] + df_vars['POS'].astype(str) + df_vars['ALT']
    # check which lineage each belongs to!
    df_vars['coreB3'] = [True if mut in b3muts else False for mut in df_vars['mutName']]
    df_vars['coreD8'] = [True if mut in d8muts else False for mut in df_vars['mutName']]
    df_vars['coreBoth'] =  df_vars['coreB3'] & df_vars['coreD8']
    df_vars['core_type'] = ['both' if df_vars['coreBoth'].iloc[j] else 'B3' if df_vars['coreB3'].iloc[j] else 'D8' if df_vars['coreD8'].iloc[j] else None for j in range(df_vars.shape[0]) ]
    df_vars['colors'] = [colors_list[0] if df_vars['core_type'].iloc[j]=='B3' else colors_list[1] if df_vars['core_type'].iloc[j]=='D8' else 'brown' if df_vars['core_type'].iloc[j]=='both' else 'black' for j in range(df_vars.shape[0]) ]
    df_vars = df_vars[(df_vars['core_type']!=st)] 
    df_vars = df_vars[(df_vars['core_type']!='both')]

    df_loc_muts = pd.read_csv(f"../tree/{st}/{st}_mutation_counts.csv")
    df_loc_muts = df_loc_muts[df_loc_muts['mutation'].isin(df_vars['mutName'])]
    df_vars['colors'] = ['cornflowerblue' if n in df_loc_muts['mutation'].values else c for n,c in zip(df_vars['mutName'],df_vars['colors'])]
    # --- Set up figure ---
    fig, ax = plt.subplots(figsize=(7, 3))
    # --- Parameters ---
    box_height = 0.1   # height of gene boxes
    gene_track_top = box_height   # top of boxes = baseline for bars

    # --- Plot gene boxes (flush with x-axis) ---
    for name, start, end in genes:
        ax.add_patch(plt.Rectangle(
            (start, 0),               # bottom-left corner
            end - start,              # width
            box_height,               # height
            facecolor="gray",
            edgecolor="white",
            alpha=0.8
        ))
        ax.text(
            (start + end) / 2,               # center of box
            box_height / 2,                  # vertical midpoint
            name,
            ha="center", va="center",
            fontsize=12, color="white", fontweight="bold", rotation=0,
            clip_on=True                     # ensure text stays within axes
        )
    
    # --- Plot markers (centered at ALT_FREQ values) ---
    ax.scatter(
        df_vars["POS"],
        gene_track_top + df_vars["median"],  # vertically position at ALT_FREQ above the gene track
        c=df_vars["colors"],
        s=df_vars['count']*7,                # marker size
        marker="o",          # change to "^", "s", etc. for other shapes
        zorder=3,
        clip_on=False,
        alpha=0.75
    )

    for _, row in df_vars.iterrows():
        # Vertical bar for IQR
        ax.vlines(
            row['POS'],
            gene_track_top + row['q25'],
            gene_track_top + row['q75'],
            color=row['colors'],
            lw=1,
            alpha=0.75
        )

    for c in [1, 5, 10]:
        ax.scatter([], [], s=c * 9, color='grey', edgecolor='grey', label=f'{c}')

    legend = ax.legend(
        title='Count',
        scatterpoints=1,
        frameon=False,
        labelspacing=0.2,
        borderpad=0.0,
        handletextpad=0.1,
        loc='upper left',
    )
    # legend.get_title().set_color('grey')
    legend.get_title().set_fontweight('bold')

    new_lookup = {colors_list[0]:'Known B3', colors_list[1]:'Known D8','cornflowerblue':'SA-specific'}
    l0 = []
    for c in new_lookup.keys():
        l0.append(ax.scatter([], [], s=18, color=c, edgecolor=c, label=f'{new_lookup[c]}'))

    legend2 = plt.legend(handles=l0,
        scatterpoints=1,
        frameon=False,
        labelspacing=0.2,
        borderpad=0.,
        handletextpad=0.0,
        loc='upper right',
        title='SNV Type'
    )
    legend2.get_title().set_fontweight('bold')
    ax.add_artist(legend)


    # --- Format axis ---
    ax.set_xlim(1, 15894)
    ax.set_ylim(0, gene_track_top + 1.0)  # extend y-limit for bars
    ax.set_xlabel("Genome position (nt)")
    ax.set_ylabel("SNV frequency")
    ax.set_yticks([gene_track_top, 0.5 + gene_track_top, 1 + gene_track_top])
    ax.set_yticklabels([0, 0.5, 1])
    ax.axhline(gene_track_top,zorder=1000,color='black',linewidth=1)
    # Hide top/right spines
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    # Keep bottom spine at the bottom of the gene boxes
    ax.spines["bottom"].set_position(("data", 0))
    ax.spines["left"].set_bounds(box_height, ax.get_ylim()[1])

    plt.tight_layout()
    plt.savefig(f'../figures/ww_genome_variation_{st}.pdf',transparent=True,bbox_inches='tight')
    plt.close('all')
