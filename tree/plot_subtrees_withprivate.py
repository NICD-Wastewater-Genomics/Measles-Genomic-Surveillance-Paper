import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from skbio.stats.distance import mantel
import re
import seaborn as sns
from Bio import SeqIO

import requests
from io import StringIO as sio
# from adjustText import adjust_text
from matplotlib.transforms import Affine2D
from matplotlib.patches import ConnectionPatch
from matplotlib_scalebar.scalebar import ScaleBar


import baltic as bt
import pandas as pd
import numpy as np
import dendropy
from dendropy import Tree
import geopandas as gpd

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

from itertools import combinations
def get_children(tree, node):
    return [k for k in tree.Objects if getattr(k, "parent", None) == node]


def get_descendant_tip_names(obj):
    if obj.is_leaf():
        return [obj.name]
    else:
        return list(obj.leaves)


def get_pairs_across_child_branches(tree,max_clade_size=4):
    pairs = []

    for node in tree.Objects:
        if node.is_leaf():
            continue

        children = get_children(tree, node)

        # keep child branches that contain at least one tip
        child_tip_sets = [
            get_descendant_tip_names(child)
            for child in children
        ]

        child_tip_sets = [
            tips for tips in child_tip_sets
            if len(tips) > 0
        ]

        all_desc_tips = get_descendant_tip_names(node)
        if len(all_desc_tips) > max_clade_size:
            continue

        # pair tips from different child branches
        for tips1, tips2 in combinations(child_tip_sets, 2):
            for a in tips1:
                for b in tips2:
                    pairs.append(tuple(sorted([a.split('|')[0].replace('-','_'), b.split('|')[0].replace('-','_')])))

    return sorted(set(pairs))

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
# read in metadata
# meta = pd.read_csv('../clinical_seqs/sequences_measles_meta.tsv',sep='\t',index_col='Accession')
# meta = meta[~meta.Collection_Date.isna()]
color_dict = {}
cInd = 0
countries_set = set()
subtrees = ['B3','D8']
gdf = gpd.read_file("../assets/map_files/zaf_admbnda_adm2_sadb_ocha_20201109.shp").set_index("ADM2_EN")
gdf= gdf[gdf['ADM1_EN']=='Gauteng']
gdf = gdf[~gdf.index.isin(['Sedibeng','West Rand'])]
gdf_r = gdf.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0)).simplify(20)
for st in subtrees:
    tree_path = f"{st}/{st}_timetree_withprivate/timetree.nexus"
    datefile_path = f"{st}/{st}_timetree/dates.tsv"
    tree = bt.loadNexus(tree_path,treestring_regex='tree1=',absoluteTime=True,verbose=False,dateFile=datefile_path)

    dist_mat = pd.read_csv(f"../src/{st.lower()}_ww_aligned_distance_matrix.csv",index_col=0)
    dist_mat.index = [dmi.replace('_WGS_25_','-WGS-25-').replace('_WGS_24_', '-WGS-24-').replace('_COV_25_','-COV-25-').replace('_COV_24_', '-COV-24-').replace('_MEV_25_','-MEC-25-').replace('_MEV_24_', '-MEV-24-').split('_')[1] for dmi in dist_mat.index]
    dist_mat.columns = [dmi.replace('_WGS_25_','-WGS-25-').replace('_WGS_24_', '-WGS-24-').replace('_COV_25_','-COV-25-').replace('_COV_24_', '-COV-24-').replace('_MEV_25_','-MEC-25-').replace('_MEV_24_', '-MEV-24-').split('_')[1] for dmi in dist_mat.columns]


    fig,ax = plt.subplots(figsize=(5,12),facecolor='w')

    all_ys=tree.getParameter('y')
    all_xs=[tO.absoluteTime for tO in tree.Objects]
    all_names=[tO.name if (tO.is_leaf() and "|" in tO.name) else f"{tO.name}|South Africa|" if (tO.is_leaf() and "|" not in tO.name)  else '' for tO in tree.Objects]
    all_countries=[tO.name.split('|')[1] if (tO.is_leaf() and '|' in tO.name) else 'SouthAfrica' if (tO.is_leaf() and '|' not in tO.name) else '' for tO in tree.Objects]
    node_xys = [[tO.absoluteTime,tO.y] for tO in tree.Objects if tO.is_node()]
    node_ids=[j for j,tO in enumerate(tree.Objects) if tO.is_node()]
    df = pd.DataFrame([all_names,all_countries,all_xs,all_ys],index=['name','country','x','y']).T
    df0 = df[df['name'].str.contains("\|")]
    # fig,ax = plt.subplots(figsize=(5,15))
    x_attr=lambda k: k.absoluteTime
    tree.plotTree(ax,x_attr=x_attr,colour='black',zorder=20000) ## tree
    # if color_option == 'lineage':
    # tree.plotPoints(ax,x_attr=x_attr,size=20,zorder=20002) ## plot coloured points at tips
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
                        # private_samps = [dfi for dfi in df0__['ID'] if dfi.startswith('5')]
                        # df_private = pd.DataFrame(columns=df_clin.columns, index=private_samps)
                        # df_private['Genotype'] = st
                        # df_private['Province'] = 'unknown'
                        # df_clin = pd.concat((df_clin,df_private),axis=0)
                        df_clin_ = df_clin.loc[df0__['ID']]
                        for k in range(df0__.shape[0]):
                            if isinstance(df_clin_.iloc[k]['Province'], float):
                                if np.isnan(df_clin_.iloc[k]['Province']):
                                    if inUnk:
                                        ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='limegreen',alpha=0.8)
                                    else:
                                        inUnk=True
                                        ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='limegreen',alpha=0.8,label='South Africa (C, unk. province)')
                            elif df_clin_.iloc[k]['Province'].lower() =='gauteng':
                                if inGau:
                                    ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='red',alpha=0.8)         
                                else:
                                    inGau=True
                                    ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='red',alpha=0.8,label='South Africa (C, Gauteng)')
                            elif df_clin_.iloc[k]['Province'].lower() =='unknown':
                                if inUnk:
                                    ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='limegreen',alpha=0.8)
                                else:
                                    inUnk=True
                                    ax.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='limegreen',alpha=0.8,label='South Africa (C, unk. province)')
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
    # leg = ax.legend(loc='lower left')
    # # leg = ax.legend(loc='best',ncol=2,fontsize=4.5,columnspacing=0.3)
    # # angle = 90  # degrees
    # leg.set_bbox_to_anchor((0.5, 0.5))  # position
    # leg._legend_box.set_transform(Affine2D().rotate_deg(90) + ax.transAxes)
    LegendVertical(ax, 90, XPad=-32, YPad=8)

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.setp(ax.get_xticklabels(), rotation=90, ha='center',va='top')
    # ax.grid(axis='x',ls='--',color='grey')
    ax.set_yticks([])
    ax.set_yticklabels([])
    [ax.spines[loc].set_visible(False) for loc in ax.spines if loc not in ['bottom']]
    if 0:
        for n in range(0,len(node_ids)):
            ax.annotate(node_ids[n],node_xys[n],color='black')
    plt.savefig(f'{st}/{st}_time_resolved_tree_withprivate.pdf',bbox_inches='tight',transparent=True)
    plt.close('all')

    #get more focused SA specific subtrees
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
        sibling_pairs = get_pairs_across_child_branches(tree)
        # sibling_pairs = get_sibling_pairs_including_polytomies(tree)
        
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

        df_test = df0[df0['country']=='SouthAfrica']
        df_test['sample_type'] = ['Clinical' if ('CVI' in name or name.startswith('Consensus_') or name.startswith('5')) else 'Wastewater' for name in df_test['name']]
        if 'Wastewater' not in df_test['sample_type'].to_list():
            fig,ax_tree = plt.subplots(figsize=(5,len(all_ys)//22+1))
            map=False
        else:
            if st=='B3':
                fig, (ax_tree, ax_map) = plt.subplots(1, 2, figsize=(11, len(all_ys)//22+1),
                                                gridspec_kw={'width_ratios': [2, 1.5], 'wspace': 0.05})
            else:
                 fig, (ax_tree, ax_map) = plt.subplots(1, 2, figsize=(11, len(all_ys)//22+1),
                                                gridspec_kw={'width_ratios': [2, 1.5], 'wspace': 0.05})               
            map=True
        x_attr=lambda k: k.absoluteTime
        tree.plotTree(ax_tree,x_attr=x_attr,colour='black',zorder=20000) ## tree
        treePlot= ax_tree.plot() ## need to call plot when only drawing the tree to force drawing of line collections
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
                                            ax_tree.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='limegreen',alpha=0.8)
                                        else:
                                            inUnk=True
                                            ax_tree.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='limegreen',alpha=0.8,label='South Africa (Clinical, unk. province)')
                                elif df_clin_.iloc[k]['Province'].lower() =='gauteng':
                                    if inGau:
                                        ax_tree.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='red',alpha=0.8)
                                    else:
                                        inGau=True
                                        ax_tree.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='red',alpha=0.8,label='South Africa (Clinical, Gauteng)')
                                elif df_clin_.iloc[k]['Province'].lower() =='unknown':
                                    if inUnk:
                                        ax_tree.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='limegreen',alpha=0.8)
                                    else:
                                        inUnk=True
                                        ax_tree.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='limegreen',alpha=0.8,label='South Africa (Clinical, unk. province)')
                                else:
                                    if inOther:
                                        ax_tree.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='purple',alpha=0.8)
                                    else:
                                        inOther=True
                                        ax_tree.scatter(df0__.iloc[k]['x'],df0__.iloc[k]['y'],s=20,zorder=30000,color='purple',alpha=0.8,label='South Africa (Clinical, outside Gauteng)')
                        else:
                            seqNames = [dfi0.split('|')[0].replace('_','-') for dfi0 in df0__['name']]
                            seq_locs = meta_ww.loc[seqNames,['SiteProvince', 'DistrictName','Latitude', 'Longitude']]
                            seq_locs['type'] = ['Wastewater' if sp=='Gauteng' else 'Wastewater(other)' for sp in seq_locs['SiteProvince']]
                            seq_locs['color'] = [c0[sl] for sl in seq_locs['type']]
                            if map:
                                seq_locs_for_map = seq_locs[seq_locs['SiteProvince']=='Gauteng']
                                seq_names_for_map = seq_locs_for_map.index.to_list()
                                pts = gpd.GeoDataFrame(seq_locs_for_map,
                                    geometry=gpd.points_from_xy(seq_locs_for_map['Longitude'], seq_locs_for_map['Latitude']),
                                    crs='EPSG:4326').to_crs('EPSG:2053')
                                
                                from scipy.spatial.distance import pdist, squareform
                                coords = pts.geometry.apply(lambda p: (p.x, p.y)).to_list()
                                # pairwise Euclidean distances in meters
                                D_m = squareform(pdist(coords, metric="euclidean"))
                                clean_inds = [re.sub(r"-(D|C)$", "",pi) for pi in  pts.index.to_list()]
                                dist_mat_geo = pd.DataFrame(D_m, index=clean_inds, columns=clean_inds)/1000
                                dist_mat_ = dist_mat.loc[dist_mat_geo.index,:]
                                dist_mat_ = dist_mat_.loc[:,dist_mat_geo.columns]
                                map_x = -pts.geometry.x.values
                                map_y = -pts.geometry.y.values
                                r, p, n = mantel(dist_mat_, dist_mat_geo, method="spearman", permutations=9999)
                                print("mantel test",r,p,n)
                                
                                sister_dist = []
                                for a, b in combinations(seq_names_for_map, 2):
                                    clean_a_name = re.sub(r"-(D|C)$", "",a)
                                    clean_b_name = re.sub(r"-(D|C)$", "",b)
                                    geo_dist = dist_mat_geo.loc[clean_a_name,clean_b_name]
                                    sameSite = True if geo_dist==0. else False
                                    a_mod = a.replace('-','_')
                                    b_mod = b.replace('-','_')
                                    areSiblings = ((a_mod, b_mod) in sibling_pairs) or ((b_mod, a_mod) in sibling_pairs)
                                    sister_dist.append([a,b,areSiblings,geo_dist,sameSite])
                                sister_df = pd.DataFrame(sister_dist,columns=['Seq1','Seq2','are_sisters','geo_dist','same_site'])
                                figSis, axSis = plt.subplots()
                                sns.boxplot(data=sister_df,x='are_sisters',y='geo_dist',ax=axSis,color='lightgrey')
                                sns.swarmplot(data=sister_df,x='are_sisters',y='geo_dist',ax=axSis,color='black',alpha=0.5)
                                axSis.set_xticklabels(['Non-Sibling','Sibling'])
                                axSis.set_xlabel('')
                                axSis.set_ylabel("Geographic distance (km)")
                                figSis.tight_layout()
                                figSis.savefig(f'{st}_{j0}_sister_dists.pdf')
                                sis_rows = sister_df[sister_df['are_sisters']]
                                obs_stat = sis_rows['same_site'].mean()
                                n_sis = len(sis_rows)
                                ## sample background of non-siblings and whether or not they're from the same site
                                background = sister_df[~sister_df['are_sisters']]['same_site'].values
                                rng = np.random.default_rng(1)
                                n_perm = 100000
                                perm_means = rng.choice(background, size=(n_perm, n_sis), replace=True).mean(axis=1)
                                # calculate finite-sampling adjusted p value for the permutation test. 
                                perm_p = (np.sum(perm_means >= obs_stat) + 1) / (n_perm + 1)
                                print(f"Permutation test: obs mean={obs_stat:.2f} km, p={perm_p:.4f} (n_perm={n_perm}, n_sis={n_sis}, n_all_pairs={len(background)})")
                                gdf_r.plot(ax=ax_map, color='#e8e8e8', edgecolor='black', lw=0.3)
                                ax_map.scatter(map_x, map_y, s=20, c=seq_locs_for_map['color'],edgecolor='black', zorder=5, alpha=0.8)
                                b = gdf_r.total_bounds
                                ax_map.set_xlim(b[0], b[2])
                                ax_map.set_ylim(b[1], b[3])
                                ax_map.set_axis_off()
                                ax_map.add_artist(ScaleBar(1,box_alpha=0,location='lower right'))

                            ax_tree.scatter(df0__['x'],df0__['y'],s=20,zorder=20002,label=f"{country.replace('_',' ')} ({st0})",color=seq_locs['color'],alpha=0.8)
                            if map:
                                df0__['shortname'] = [dfi0.split('|')[0].replace('_','-') for dfi0 in df0__['name']]
                                df0__ = df0__.set_index('shortname')
                                for i in range(len(seq_names_for_map)):
                                    con = ConnectionPatch(
                                        xyA=(df0__.loc[seq_names_for_map[i]]['x'], df0__.loc[seq_names_for_map[i]]['y']),
                                        xyB=(map_x[i], map_y[i]),
                                        coordsA='data', coordsB='data',
                                        axesA=ax_tree, axesB=ax_map,
                                        color=seq_locs_for_map['color'].iloc[i],
                                        lw=1, alpha=0.8, clip_on=False,
                                    )
                                    fig.add_artist(con)

                elif country !='UNKNOWN':
                    ax_tree.scatter(df0_['x'],df0_['y'],s=20,zorder=20002,label=country.replace('_',' '),color=color_dict[country],alpha=0.8)
                else:
                    ax_tree.scatter(df0_['x'],df0_['y'],s=20,zorder=20002,color='lightgrey',alpha=0.8)

        plt.gcf().canvas.draw()
        leg = ax_tree.legend(ncol=1,loc='lower left',fontsize=5)

        ax_tree.xaxis.set_major_locator(MaxNLocator(integer=True))
        # plt.setp(ax_tree.get_xticklabels(), rotation=90, ha='center',va='top')
        ax_tree.set_yticks([])
        ax_tree.set_yticklabels([])
        [ax_tree.spines[loc].set_visible(False) for loc in ax_tree.spines if loc not in ['bottom']]
        if 0:
            for n in range(0,len(node_ids)):
                ax_tree.annotate(node_ids[n],node_xys[n],color='black')
        fig.savefig(f'{st}/{st}_focused_tree{str(j0)}_withprivate.pdf',bbox_inches='tight',transparent=True)
        plt.close('all')

        countries_set = countries_set | set(countries)
    if 1:# generate unique colorset
        k0 = 0
        country_colors = {}
        for c in list(countries_set):
            if c != 'South_Africa':
                country_colors[c] = plt.cm.tab20(k0)
                k0+=1

