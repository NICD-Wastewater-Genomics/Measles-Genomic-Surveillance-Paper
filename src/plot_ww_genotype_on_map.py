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
import geopandas as gpd
from shapely.geometry import Point, Polygon
from matplotlib_scalebar.scalebar import ScaleBar
from adjustText import adjust_text
import matplotlib

# This script generates Fig 2d,e and Fig 3b,c

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['font.family'] = 'Arial'

meta_ww = pd.read_csv('../metadata/Measles_seqdata_02022026.csv')
meta_ww.iloc[:,2] = meta_ww.iloc[:,2].astype(str)
meta_ww.index = [meta_ww['SampleID'].iloc[j] if meta_ww.iloc[j,2].startswith('Not') else meta_ww.iloc[j,2] for j in range(meta_ww.shape[0])]
meta_ww.index = [mi.replace('_','-').split('-S')[0] for mi in meta_ww.index]

agg_df = pd.read_csv('../agg_demixed.tsv', skipinitialspace=True, sep='\t',index_col=0)
agg_df = agg_df[~agg_df.index.str.contains('PC_')]
from freyja.utils import prepLineageDict, prepSummaryDict
agg_df = prepSummaryDict(agg_df)
agg_df['abundances'] = agg_df['abundances'].astype(str)
agg_df_likes = prepLineageDict(agg_df, thresh=0.0000000001,  mergeLikes=False)
agg_df_likes.index =  [(dsi.split('.tsv')[0].split('_')[-1]) for dsi in agg_df_likes.index]
agg_df = prepLineageDict(agg_df, thresh=0.0000000001,  mergeLikes=True)

agg_df['region'] = agg_df['linDict'].apply(lambda x:list(x.keys())[0].split('-')[0].split('MEASLES')[1])
agg_df2 = agg_df.copy()
agg_df2 = agg_df2[agg_df2.index.str.contains('ENV') | agg_df2.index.str.contains('NAT') | agg_df2.index.str.contains('CST') | agg_df2.index.str.contains('ART_MEV') ]

agg_df2.index = [agi.split('_S')[0] for agi in agg_df2.index]
agg_df2['ID'] = ['-'.join(agi.split('_')[0:5])  if agi.endswith('_C') or agi.endswith('_D') else '-'.join(agi.split('_')[0:4]) for agi in agg_df2.index]
agg_df2['ID_min'] = ['-'.join(agi.split('-')[0:4])  if agi.endswith('-C') or agi.endswith('-D') else agi for agi in agg_df2['ID']]

idx = agg_df2.groupby('ID_min')['coverage'].idxmax()
agg_df2 = agg_df2.loc[idx]

## check to make sure nothing is missing from metadata
print('Missing from metadata!: ',agg_df2[~agg_df2['ID'].isin(meta_ww.index)]['ID'].values)
agg_df2 = agg_df2[agg_df2['ID'].isin(meta_ww.index)] 
agg_df2['site'] = meta_ww.loc[agg_df2['ID'].values,'SiteName'].values
agg_df2['Province'] = meta_ww.loc[agg_df2['ID'].values,'SiteProvince'].values
agg_df2['date'] = pd.to_datetime(meta_ww.loc[agg_df2['ID'].values,'SampleCollectionDate'].values)

#read in location of each wwtp
df = pd.read_csv('../assets/All sampling sites_v5_300502024.csv',skipinitialspace=True)
df['Latitude'] = df['GPS latitude']
df['Longitude'] = df['GPS longitude']

## restrict to Gauteng sites, since data is too limited elsewhere
df = df[df['Site'].isin(meta_ww['SiteName'])]
df = df[df.Province=='Gauteng']

sites = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude), crs="EPSG:4326")

seq_sites = sites
district_lookup= {"City of Johannesburg":'Johannesburg MM',
"City of Tshwane":"Tshwane MM",
"Ekurhuleni":"Ekurhuleni MM",
}
### then add country map
gdf = gpd.read_file("../assets/map_files/zaf_admbnda_adm2_sadb_ocha_20201109.shp").set_index("ADM2_EN")
covered_gdf = gdf[gdf['ADM1_EN']=='Gauteng']
covered_gdf = covered_gdf[covered_gdf.index.isin(district_lookup.keys())]

#read in location of each wwtp
dfCities = pd.read_csv('../assets/cities.tsv',sep='\t')
dfCities['Latitude'] = dfCities['Coords'].apply(lambda x:x.split(',')[0])
dfCities['Longitude'] = dfCities['Coords'].apply(lambda x:x.split(',')[1])
dfCities = dfCities.iloc[0:2]

cities = gpd.GeoDataFrame(dfCities, geometry=gpd.points_from_xy(dfCities.Longitude, dfCities.Latitude), crs="EPSG:4326")
cities_r = cities.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0))

### 
# convert map to projected WGS 84. 
covered_gdf_r = covered_gdf.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0)).simplify(20)

## make pie chart markers
def plot_pie(locX,locY,ratios,ax0,scalar=1,edgecolor='black',c0 = None,zorder=100):
    # first define the ratios
    r = np.cumsum(ratios)
    # define some sizes of the scatter marker
    xy,s1 = [],[]
    for j,r0 in enumerate(r):
        if r0==0:
            continue
        if j>0:
            if r[j]==r[j-1]:
                continue
        if j==0:
            x = [0] + np.cos(np.linspace(0, 2 * np.pi * r0, 100)).tolist()
            y = [0] + np.sin(np.linspace(0, 2 * np.pi * r0, 100)).tolist()
        else:
            x = [0] + np.cos(np.linspace(2 * np.pi * r[j-1], 2 * np.pi * r0, 100)).tolist()
            y = [0] + np.sin(np.linspace(2 * np.pi * r[j-1], 2 * np.pi * r0, 100)).tolist()
        xy= np.column_stack([x, y])
        s1 = np.abs(xy).max()
        if c0 is None:
            ax0.scatter(locX, locY, marker=xy,s=s1 ** 2 * 30*scalar, facecolor=plt.cm.Set2(j),alpha=0.7,edgecolor='none',zorder=zorder)
            ax0.scatter(locX, locY, marker='o',s=s1 ** 2 * 30*scalar,facecolor='none',edgecolor=edgecolor,linewidth=0.2,zorder=zorder+0.1)
        else: 
            ax0.scatter(locX, locY, marker='o',s=s1 ** 2 * 30*scalar,facecolor=c0,edgecolor=edgecolor,linewidth=0.2,zorder=zorder+0.1)

# convert to projected WGS 84. 
# Hartebeesthoek94 has a -1 scale factor built in, so we reverse it. 
seq_sites_r = seq_sites.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0))
seq_sites['geometry']= seq_sites_r
seq_sites_r = seq_sites

## Build Figures 2d and 2e
df = df.set_index('Site')
seq_sites = seq_sites.set_index('Site')
seq_sites = seq_sites[seq_sites.index.isin(df.index)]
for type0 in seq_sites.Type.unique():
    fig,ax = plt.subplots()
    ax = covered_gdf_r.plot(facecolor='papayawhip',edgecolor='silver')
    cities_r.plot(ax=ax,edgecolor='forestgreen',facecolor='forestgreen',markersize=7)
    site_grps = agg_df2.groupby('site').agg({'lineages':list,'abundances':list})
    avgs = []
    for l,y in zip(site_grps['lineages'],site_grps['abundances']):
        avg = {}
        for ls,ys in zip(l,y):
            for l0,y0 in zip(ls,ys):
                if l0 in avg.keys():
                    avg[l0]+=y0
                else:
                    avg[l0]=y0
        for l0 in avg.keys():
            avg[l0]/=len(l)
        avgs.append(avg)
    site_grps['avgs'] = avgs
    seq_sites0 = seq_sites[seq_sites['Type']==type0]
    allKeys = [list(avg.keys()) for avg in avgs]
    allKeys = np.unique([a for ak in allKeys for a in ak])
    jit = 2500
    seq_sites0['num_detects'] = [ len(site_grps.loc[ind,'lineages']) if ind in site_grps.index else 0 for ind in seq_sites0.index ]
    seq_sites0 = seq_sites0.sort_values(by='num_detects',ascending=False)
    for jj,ind in enumerate(seq_sites0.index):
        if ind in site_grps.index:
            r = site_grps.loc[ind,'avgs']
            r0 = []
            for key in allKeys:
                if key in r:
                    r0.append(r[key])
                else:
                    r0.append(0)
            if seq_sites.loc[ind,'Type']=='National':     
                plot_pie(seq_sites.loc[ind,'geometry'].x,seq_sites.loc[ind,'geometry'].y,r0,ax,scalar=len(site_grps.loc[ind,'lineages']))
            else:
                plot_pie(seq_sites.loc[ind,'geometry'].x + np.random.uniform(-jit, jit),seq_sites.loc[ind,'geometry'].y+ np.random.uniform(-jit, jit),r0,ax,scalar=len(site_grps.loc[ind,'lineages']),edgecolor='white',zorder=100+jj)
 
 
    df_key = pd.DataFrame([[1,-26.15, 28.739551]],columns=['count','Latitude','Longitude'])
    df_key_ = gpd.GeoDataFrame(df_key, geometry=gpd.points_from_xy(df_key.Longitude, df_key.Latitude), crs="EPSG:4326")
    df_key_ = df_key_.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0))
    df_key['geometry'] = df_key_
    # add in 5, 10

    from shapely.geometry import Point
    for j,v in enumerate([5,10]):
        newRow = df_key.iloc[j]
        newRow['geometry'] = Point(newRow['geometry'].x,newRow['geometry'].y-1E4)
        newRow['count'] = v
        df_key = pd.concat((df_key,pd.DataFrame([newRow])),axis=0)

    df_key = df_key.reset_index()
    for j in df_key.index:
        plot_pie(df_key.loc[j,'geometry'].x,df_key.loc[j,'geometry'].y,1,ax,edgecolor='black',c0='white',scalar=df_key.loc[j,'count'])
        ax.text(df_key.loc[j,'geometry'].x+4000,df_key.loc[j,'geometry'].y,df_key.loc[j,'count'],color='black',horizontalalignment='left',verticalalignment='center')
    ax.text(df_key.loc[0,'geometry'].x,df_key.loc[0,'geometry'].y+6e3,'Samples',color='black',horizontalalignment='center',verticalalignment='center')

    labels = []
    for x, y, label in zip(cities_r.geometry.x, cities_r.geometry.y, dfCities.City):
        labels.append(plt.text(x, y,label.split(' ')[0],fontsize=8,color='forestgreen'))

    adjust_text(labels,arrowprops=dict(arrowstyle="-", color='k', lw=0.5),explode_radius=200)

    ax.add_artist(ScaleBar(1,box_alpha=0,location='lower right'))
    plt.axis('off')
    plt.savefig(f'../figures/map_lineages_{type0}.pdf')
    plt.close('all')

## Build Figures 3b and 3c
agg_df2['quarter'] = agg_df2['date'].dt.to_period('Q')
fig,ax1 = plt.subplots(2,4)
all_quarters = agg_df2['quarter'].unique()
for jj, q in enumerate(all_quarters):
    #### same, but now broken down by QE
    for ii, type0 in enumerate(seq_sites.Type.unique()):
        ax1[ii,jj] = covered_gdf_r.plot(facecolor='papayawhip',edgecolor='silver',ax=ax1[ii,jj])#,label='Sequencing Data')
        agg_df2_q = agg_df2[agg_df2['quarter']==q]
        site_grps = agg_df2_q.groupby('site').agg({'lineages':list,'abundances':list})
        avgs = []
        for l,y in zip(site_grps['lineages'],site_grps['abundances']):
            avg = {}
            for ls,ys in zip(l,y):
                for l0,y0 in zip(ls,ys):
                    if l0 in avg.keys():
                        avg[l0]+=y0
                    else:
                        avg[l0]=y0
            for l0 in avg.keys():
                avg[l0]/=len(l)
            avgs.append(avg)
        site_grps['avgs'] = avgs
        seq_sites0 = seq_sites[seq_sites['Type']==type0]
        allKeys = [list(avg.keys()) for avg in avgs]
        allKeys = np.unique([a for ak in allKeys for a in ak])
        jit = 2500
        seq_sites0['num_detects'] = [ len(site_grps.loc[ind,'lineages']) if ind in site_grps.index else 0 for ind in seq_sites0.index ]
        seq_sites0 = seq_sites0.sort_values(by='num_detects',ascending=False)
        for jj_,ind in enumerate(seq_sites0.index):
            if ind in site_grps.index:
                r = site_grps.loc[ind,'avgs']
                r0 = []
                for key in allKeys:
                    if key in r:
                        r0.append(r[key])
                    else:
                        r0.append(0)
                if seq_sites.loc[ind,'Type']=='National':     
                    plot_pie(seq_sites.loc[ind,'geometry'].x,seq_sites.loc[ind,'geometry'].y,r0,ax1[ii,jj],scalar=len(site_grps.loc[ind,'lineages']))
                else:
                    plot_pie(seq_sites.loc[ind,'geometry'].x + np.random.uniform(-jit, jit),seq_sites.loc[ind,'geometry'].y+ np.random.uniform(-jit, jit),r0,ax1[ii,jj],scalar=len(site_grps.loc[ind,'lineages']),edgecolor='white',zorder=100+jj_)
    
        if ii==0 and jj==0:
            df_key = pd.DataFrame([[1,-26.15, 28.739551]],columns=['count','Latitude','Longitude'])
            df_key_ = gpd.GeoDataFrame(df_key, geometry=gpd.points_from_xy(df_key.Longitude, df_key.Latitude), crs="EPSG:4326")
            df_key_ = df_key_.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0))
            df_key['geometry'] = df_key_
            # add in 5, 10

            from shapely.geometry import Point
            for j,v in enumerate([2,5]):
                newRow = df_key.iloc[j]
                newRow['geometry'] = Point(newRow['geometry'].x,newRow['geometry'].y-1.5E4)
                newRow['count'] = v
                df_key = pd.concat((df_key,pd.DataFrame([newRow])),axis=0)

            df_key = df_key.reset_index()
            for j in df_key.index:
                plot_pie(df_key.loc[j,'geometry'].x,df_key.loc[j,'geometry'].y,1,ax1[ii,jj],edgecolor='black',c0='white',scalar=df_key.loc[j,'count'])
                ax1[ii,jj].text(df_key.loc[j,'geometry'].x+10000,df_key.loc[j,'geometry'].y,df_key.loc[j,'count'],color='black',horizontalalignment='left',verticalalignment='center')
            ax.text(df_key.loc[0,'geometry'].x,df_key.loc[0,'geometry'].y+6e3,'Samples',color='black',horizontalalignment='center',verticalalignment='center')

        labels = []
        if ii==1 and jj==3:
            ax1[ii,jj].add_artist(ScaleBar(1,box_alpha=0,location='lower right'))
            ax1[ii,jj].legend(loc='upper left')
        ax1[ii,jj].set_axis_off()
        if ii==0:
            ax1[ii,jj].set_title(f"Q{q.quarter}\n{q.year}")
        if jj==0:
            ax1[ii,jj].annotate(
            type0,
            xy=(-0.02, 0.5),
            xycoords='axes fraction',
            ha='right',
            va='center',
            fontsize=13,
            fontweight='bold',
            rotation=90
    )

fig.subplots_adjust(
    left=0.0,
    right=1.0,
    top=1.,
    bottom=0.0,
    wspace=-0.08,
    hspace=-0.42
)      
fig.savefig(f'../figures/map_lineages_all_quarters.pdf')