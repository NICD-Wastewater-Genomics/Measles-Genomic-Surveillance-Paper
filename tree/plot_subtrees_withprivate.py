import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from Bio import SeqIO
import requests
from io import StringIO as sio
from matplotlib.transforms import Affine2D
import baltic as bt
import pandas as pd
import numpy as np
import dendropy
from dendropy import Tree

# Script for generating genotype specific tree plots and
# subtrees containing clusters with sequences from south africa 

df_clin = pd.read_csv('../metadata/mev clinical metadata_updated 040226 ks_JL_090226.csv')
df_clin = df_clin.dropna(subset=['Seq ID']).set_index('Seq ID')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.sans-serif'] = ['Arial']


meta_ww = pd.read_csv('../metadata/Measles_seqdata_02022026.csv')
meta_ww.iloc[:,2] = meta_ww.iloc[:,2].astype(str)
meta_ww.index = [meta_ww['SampleID'].iloc[j] if meta_ww.iloc[j,2].startswith('Not') else meta_ww.iloc[j,2] for j in range(meta_ww.shape[0])] #,index_col='SampleID')#pd.read_csv('../metadata/MeV Wastewater Sequences Final Metadata_05092025.csv',index_col='Sample Number')
meta_ww = meta_ww[[mc for mc in meta_ww.columns if 'Unnamed' not in mc]]
meta_ww.index = [mi.replace('_','-').split('-S')[0] for mi in meta_ww.index]
meta_ww = meta_ww[~meta_ww['SampleID'].str.startswith('AIR')]

def LegendVertical(Ax, Rotation=90, XPad=0, YPad=0, **LegendArgs):
    if Rotation not in (90,270):
        raise NotImplementedError('Rotation must be 90 or 270.')

    # Extra spacing between labels is needed to fit the rotated labels;
    # and since the frame will not adjust to the rotated labels, it is
    # disabled by default
    DefaultLoc = 'center left' if Rotation==90 else 'center right'
    ArgsDefaults = dict(loc=DefaultLoc, labelspacing=4, frameon=False)
    Args = {**ArgsDefaults, **LegendArgs}

    Handles, Labels = Ax.get_legend_handles_labels()
    if Rotation==90:
        # Reverse entries
        Handles, Labels = (reversed(_) for _ in (Handles, Labels))
    AxLeg = Ax.legend(Handles, Labels, **Args)

    LegTexts = AxLeg.get_texts()
    LegHandles = AxLeg.legend_handles

    for L,Leg in enumerate(LegHandles):
        if type(Leg) == mpl.patches.Rectangle:
            BBounds = np.ravel(Leg.get_bbox())
            BBounds[2:] = BBounds[2:][::-1]
            Leg.set_bounds(BBounds)

            LegPos = (
                # Ideally,
                #    `(BBounds[0]+(BBounds[2]/2)) - AxLeg.handletextpad`
                # should be at the horizontal center of the legend patch,
                # but for some reason it is not. Therefore the user will
                # need to specify some padding.
                (BBounds[0]+(BBounds[2]/2)) - AxLeg.handletextpad + XPad,

                # Similarly, `(BBounds[1]+BBounds[3])` should be at the vertical
                # top of the legend patch, but it is not.
                (BBounds[1]+BBounds[3])+YPad
            )

        elif type(Leg) == mpl.lines.Line2D:
            LegXY = Leg.get_xydata()[:,::-1]
            Leg.set_data(*(LegXY[:,_] for _ in (0,1)))

            LegPos = (
                LegXY[0,0] - AxLeg.handletextpad + XPad,
                max(LegXY[:,1]) + YPad
            )

        elif type(Leg) == mpl.collections.PathCollection:
            LegPos = (
                Leg.get_offsets()[0][0] + XPad,
                Leg.get_offsets()[0][1] + YPad,
            )
        else:
            raise NotImplementedError('Legends should contain Rectangle, Line2D or PathCollection.')

        PText = LegTexts[L]
        PText.set_verticalalignment('bottom')
        PText.set_rotation(Rotation)
        PText.set_x(LegPos[0])
        PText.set_y(LegPos[1])

    return(None)

# load nexus trees and make plots
color_dict = {}
cInd = 0
countries_set = set()
subtrees = ['B3','D8']
for st in subtrees:
    tree_path = f"{st}/{st}_timetree_withprivate/timetree.nexus"
    datefile_path = f"{st}/{st}_timetree/dates.tsv"
    tree = bt.loadNexus(tree_path,treestring_regex='tree1=',absoluteTime=True,verbose=True,dateFile=datefile_path)

    fig,ax = plt.subplots(figsize=(5,12),facecolor='w')
    all_ys=tree.getParameter('y')
    all_xs=[tO.absoluteTime for tO in tree.Objects]
    all_names=[tO.name if (tO.is_leaf() and "|" in tO.name) else f"{tO.name}|South Africa|" if (tO.is_leaf() and "|" not in tO.name)  else '' for tO in tree.Objects]
    all_countries=[tO.name.split('|')[1] if (tO.is_leaf() and '|' in tO.name) else 'SouthAfrica' if (tO.is_leaf() and '|' not in tO.name) else '' for tO in tree.Objects]
    node_xys = [[tO.absoluteTime,tO.y] for tO in tree.Objects if tO.is_node()]
    node_ids=[j for j,tO in enumerate(tree.Objects) if tO.is_node()]
    df = pd.DataFrame([all_names,all_countries,all_xs,all_ys],index=['name','country','x','y']).T
    df0 = df[df['name'].str.contains("\|")]
    x_attr=lambda k: k.absoluteTime
    tree.plotTree(ax,x_attr=x_attr,colour='black',zorder=20000) ## tree
    treePlot= ax.plot() ## need to call plot when only drawing the tree to force drawing of line collections
    countries = df0.country.unique()
    for j,country in enumerate(countries):
        if country is not None:
            df0_ = df0[df0['country']==country]
            if country == 'SouthAfrica':
                df0_['sample_type'] = ['Clinical' if ('CVI' in name or name.startswith('Consensus_') or name.startswith('5'))else 'Wastewater' for name in df0_['name']]
                stypes = df0_['sample_type'].unique()
                c0 = {'Clinical':'red','Wastewater':'royalblue','Wastewater(other)':'darkturquoise'}
                for st0 in stypes:
                    df0__ =df0_[df0_['sample_type']==st0]
                    if st0 =='Clinical':
                        inUnk,inGau,inOther = False,False,False
                        df0__['ID'] = df0__['name'].apply(lambda x: x.split('|')[0])

                        df_clin_ = df_clin.loc[df0__['ID']]
                        for k in range(df0__.shape[0]):
                            if isinstance(df_clin_.iloc[k]['Province'], float):
                                if np.isnan(df_clin_.iloc[k]['Province']):
                                    if inUnk:
                                        ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='green',alpha=0.8)
                                    else:
                                        inUnk=True
                                        ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='green',alpha=0.8,label='South Africa (C, unk. province)')
                            elif df_clin_.iloc[k]['Province'].lower() =='gauteng':
                                if inGau:
                                    ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='red',alpha=0.8)         
                                else:
                                    inGau=True
                                    ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='red',alpha=0.8,label='South Africa (C, Gauteng)')
                            elif df_clin_.iloc[k]['Province'].lower() =='unknown':
                                if inUnk:
                                    ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='green',alpha=0.8)
                                else:
                                    inUnk=True
                                    ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='green',alpha=0.8,label='South Africa (C, unk. province)')
                            else:
                                if inOther:
                                    ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='purple',alpha=0.8)
                                else:
                                    inOther=True
                                    ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='purple',alpha=0.8,label='South Africa (C, outside Gauteng)')
                    else:
                        print('checking',df0__)
                        seqNames = [dfi0.split('|')[0].replace('_','-') for dfi0 in df0__['name']]
                        seq_locs = meta_ww.loc[seqNames,['SiteProvince', 'DistrictName']]
                        seq_locs['type'] = ['Wastewater' if sp=='Gauteng' else 'Wastewater(other)' for sp in seq_locs['SiteProvince']]
                        seq_locs['color'] = [c0[sl] for sl in seq_locs['type']]
                        ax.scatter(df0__['x'],df0__['y'],s=20,zorder=20002,label=f"{country.replace('_',' ')} ({st0})",color=seq_locs['color'],alpha=0.8)
            else:
                ax.scatter(df0_['x'],df0_['y'],s=20,zorder=20002,color='lightgrey',alpha=0.8)

    LegendVertical(ax, 90, XPad=-32, YPad=8)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.setp(ax.get_xticklabels(), rotation=90, ha='center',va='top')
    ax.set_yticks([])
    ax.set_yticklabels([])
    [ax.spines[loc].set_visible(False) for loc in ax.spines if loc not in ['bottom']]
    if 0:
        for n in range(0,len(node_ids)):
            ax.annotate(node_ids[n],node_xys[n],color='black')
    plt.savefig(f'{st}/{st}_time_resolved_tree_withprivate.pdf',bbox_inches='tight',transparent=True)
    plt.close('all')

    # now get more focused subtrees by pulling from specific nodes in the tree
    import copy
    nodeNames = {'B3':['NODE_0000163'],'D8':['NODE_0000253','NODE_0000021','NODE_0000317']}
    # nodeNames = {'B3':[],'D8':[]}
    tree0 = copy.deepcopy(tree)
    for j0 in range(len(nodeNames[st])):
        ## get the internal nodes we want. 
        nodes =[tO for tO in tree0.Objects if tO.is_node()]
        n = [n0 for n0 in nodes if n0.traits['label'] in nodeNames[st][j0]]
        # newRoot = tree.Objects[nodes[st]]
        tree = tree.subtree(starting_node=n[0],stem=False)
        ## reset basic variables for plotting
        all_ys=tree.getParameter('y')
        all_xs=[tO.absoluteTime for tO in tree.Objects]
        all_names=[tO.name if (tO.is_leaf() and "|" in tO.name) else f"{tO.name}|South Africa|" if (tO.is_leaf() and "|" not in tO.name)  else '' for tO in tree.Objects]
        all_countries=[tO.name.split('|')[1] if (tO.is_leaf() and '|' in tO.name) else 'SouthAfrica' if (tO.is_leaf() and '|' not in tO.name) else '' for tO in tree.Objects]
        for country in all_countries:
            if country not in color_dict.keys():
                if cInd <8:
                    color_dict[country] = plt.cm.Dark2(cInd)
                else:
                    color_dict[country] = plt.cm.Accent(cInd-8)
                cInd+=1
        print(color_dict)
        df = pd.DataFrame([all_names,all_countries,all_xs,all_ys],index=['name','country','x','y']).T
        df0 = df[df['name'].str.contains("\|")]


        fig,ax = plt.subplots(figsize=(5,len(all_ys)//22+1))
        x_attr=lambda k: k.absoluteTime
        tree.plotTree(ax,x_attr=x_attr,colour='black',zorder=20000) ## tree
        treePlot= ax.plot() ## need to call plot when only drawing the tree to force drawing of line collections
        countries = df0.country.unique()
        for j,country in enumerate(countries):
            if country is not None:
                df0_ = df0[df0['country']==country]
                if country == 'SouthAfrica':
                    df0_['sample_type'] = ['Clinical' if ('CVI' in name or name.startswith('Consensus_') or name.startswith('5')) else 'Wastewater' for name in df0_['name']]
                    stypes = df0_['sample_type'].unique()
                    c0 = {'Clinical':'red','Wastewater':'royalblue','Wastewater(other)':'darkturquoise'}
                    for st0 in stypes:
                        df0__ =df0_[df0_['sample_type']==st0]
                        if st0 =='Clinical':
                            inUnk,inGau,inOther = False,False,False
                            df0__['ID'] = df0__['name'].apply(lambda x: x.split('|')[0])
                            df_clin_ = df_clin.loc[df0__['ID']]
                            for k in range(df0__.shape[0]):
                                if isinstance(df_clin_.iloc[k]['Province'], float):
                                    if np.isnan(df_clin_.iloc[k]['Province']):
                                        if inUnk:
                                            ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='green',alpha=0.8)
                                        else:
                                            inUnk=True
                                            ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='green',alpha=0.8,label='South Africa (Clinical, unk. province)')
                                elif df_clin_.iloc[k]['Province'].lower() =='gauteng':
                                    if inGau:
                                        ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='red',alpha=0.8)         
                                    else:
                                        inGau=True
                                        ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='red',alpha=0.8,label='South Africa (Clinical, Gauteng)')
                                elif df_clin_.iloc[k]['Province'].lower() =='unknown':
                                    if inUnk:
                                        ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='green',alpha=0.8)
                                    else:
                                        inUnk=True
                                        ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='green',alpha=0.8,label='South Africa (Clinical, unk. province)')
                                else:
                                    if inOther:
                                        ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='purple',alpha=0.8)
                                    else:
                                        inOther=True
                                        ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='purple',alpha=0.8,label='South Africa (Clinical, outside Gauteng)')
                        else:
                            seqNames = [dfi0.split('|')[0].replace('_','-') for dfi0 in df0__['name']]
                            seq_locs = meta_ww.loc[seqNames,['SiteProvince', 'DistrictName']]
                            seq_locs['type'] = ['Wastewater' if sp=='Gauteng' else 'Wastewater(other)' for sp in seq_locs['SiteProvince']]
                            seq_locs['color'] = [c0[sl] for sl in seq_locs['type']]
                            ax.scatter(df0__['x'],df0__['y'],s=20,zorder=20002,label=f"{country.replace('_',' ')} ({st0})",color=seq_locs['color'],alpha=0.8)
              elif country !='UNKNOWN':
                    ax.scatter(df0_['x'],df0_['y'],s=20,zorder=20002,label=country.replace('_',' '),color=color_dict[country],alpha=0.8)
                else:
                    ax.scatter(df0_['x'],df0_['y'],s=20,zorder=20002,color='lightgrey',alpha=0.8)

        plt.gcf().canvas.draw()
        leg = ax.legend(ncol=1,loc='best',fontsize=5)

        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.setp(ax.get_xticklabels(), rotation=90, ha='center',va='top')
        ax.set_yticks([])
        ax.set_yticklabels([])
        [ax.spines[loc].set_visible(False) for loc in ax.spines if loc not in ['bottom']]
        if 0:
            for n in range(0,len(node_ids)):
                ax.annotate(node_ids[n],node_xys[n],color='black')
        plt.savefig(f'{st}/{st}_focused_tree{str(j0)}_withprivate.pdf',bbox_inches='tight',transparent=True)
        plt.close('all')

        countries_set = countries_set | set(countries)
    if 1:# generate unique colorset
        k0 = 0
        country_colors = {}
        for c in list(countries_set):
            if c != 'South_Africa':
                country_colors[c] = plt.cm.tab20(k0)
                k0+=1

