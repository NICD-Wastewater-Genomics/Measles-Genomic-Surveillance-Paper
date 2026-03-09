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
matplotlib.rcParams['font.family'] = 'Arial'

clin_df = pd.read_csv('../metadata/mev clinical metadata_updated 040226 ks_JL_090226.csv',sep=";")

clin_df['Onset date'] = pd.to_datetime(clin_df['Onset date'])
clin_df = clin_df[['Sample ID','Seq ID','Onset date','Province','City/district','WHO naming','WGS','Genotype']]
clin_df = clin_df[clin_df['Onset date']>='2024-10-01']
clin_df = clin_df[clin_df['Onset date']<'2025-10-01']
asdf
district_lookup= {"City of Johannesburg":"city of johannesburg",
"City of Tshwane":"city of tshwane",
"Ekurhuleni":"ekurhuleni"
}

district_lookup_r = {v: k for k, v in district_lookup.items()}

### then add country map
gdf = gpd.read_file("../assets/map_files/zaf_admbnda_adm2_sadb_ocha_20201109.shp").set_index("ADM2_EN")
covered_gdf = gdf[gdf['ADM1_EN']=='Gauteng']
### excude metros 
covered_gdf = covered_gdf[covered_gdf.index.isin(district_lookup.keys())]
# convert map to projected WGS 84. 
covered_gdf_r = covered_gdf.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0)).simplify(20)


## define function to make pie chart markers
def plot_pie(locX,locY,ratios,ax0,scalar=1,edgecolor='black',c0 = None,zorder=100):
    # first define the ratios
    r = np.cumsum(ratios)
    # define some sizes of the scatter marker
    xy,s1 = [],[]
    for j0,r0 in enumerate(r):
        if r0==0:
            continue
        if j0>0:
            if r[j0]==r[j0-1]:
                continue
        if j0==0:
            x = [0] + np.cos(np.linspace(0, 2 * np.pi * r0, 100)).tolist()
            y = [0] + np.sin(np.linspace(0, 2 * np.pi * r0, 100)).tolist()
        else:
            x = [0] + np.cos(np.linspace(2 * np.pi * r[j0-1], 2 * np.pi * r0, 100)).tolist()
            y = [0] + np.sin(np.linspace(2 * np.pi * r[j0-1], 2 * np.pi * r0, 100)).tolist()
        xy= np.column_stack([x, y])
        s1 = np.abs(xy).max()
        if c0 is None:
            ax0.scatter(locX, locY, marker=xy,s=s1 ** 2 * 30*scalar, facecolor=plt.cm.Set2(j0),alpha=0.7,edgecolor='none',zorder=zorder)
            ax0.scatter(locX, locY, marker='o',s=s1 ** 2 * 30*scalar,facecolor='none',edgecolor=edgecolor,linewidth=0.2,zorder=zorder+0.1)
        else: 
            ax0.scatter(locX, locY, marker='o',s=s1 ** 2 * 30*scalar,facecolor=c0,edgecolor=edgecolor,linewidth=0.2,zorder=zorder+0.1)


## build Fig 3a
fig,ax = plt.subplots()
ax = covered_gdf_r.plot(facecolor='papayawhip',edgecolor='silver')
seq_sites = covered_gdf_r.centroid
clin_df['City/district'] = clin_df['City/district'].str.lower()
clin_df['City/district'] = clin_df['City/district'].str.split(' metro').str[0]
clin_df['City/district'] = clin_df['City/district'].str.split(' mm').str[0]
clin_df = clin_df[clin_df['City/district'].isin(district_lookup.values())]

overall = clin_df.groupby('City/district')['Genotype'].value_counts()
overall_fracs = overall/overall.groupby(level=0).sum()
dists = overall_fracs.reset_index(level=1).index
gts = overall_fracs.reset_index(level=0).index.unique().sort_values()
all_gts = gts
for jj,ind in enumerate(dists):
    r0 = []
    for gt in gts:
        r0.append(overall_fracs.loc[(ind,gt)])
    plot_pie(seq_sites.loc[district_lookup_r[ind]].x,seq_sites.loc[district_lookup_r[ind]].y,r0,ax,scalar=overall.loc[ind].sum()/7)

df_key = pd.DataFrame([[10,-26.15, 28.739551]],columns=['count','Latitude','Longitude'])
df_key_ = gpd.GeoDataFrame(df_key, geometry=gpd.points_from_xy(df_key.Longitude, df_key.Latitude), crs="EPSG:4326")
df_key_ = df_key_.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0))
df_key['geometry'] = df_key_

from shapely.geometry import Point
for j,v in enumerate([20,40]):
    newRow = df_key.iloc[j]
    newRow['geometry'] = Point(newRow['geometry'].x,newRow['geometry'].y-1E4)
    newRow['count'] = v
    df_key = pd.concat((df_key,pd.DataFrame([newRow])),axis=0)

df_key = df_key.reset_index()
for j in df_key.index:
    plot_pie(df_key.loc[j,'geometry'].x,df_key.loc[j,'geometry'].y,1,ax,edgecolor='black',c0='white',scalar=df_key.loc[j,'count']/7)
    ax.text(df_key.loc[j,'geometry'].x+4000,df_key.loc[j,'geometry'].y,df_key.loc[j,'count'],color='black',horizontalalignment='left',verticalalignment='center')
ax.text(df_key.loc[0,'geometry'].x,df_key.loc[0,'geometry'].y+6e3,'Clinical\nSequences',color='black',horizontalalignment='center',verticalalignment='center')
ax.add_artist(ScaleBar(1,box_alpha=0,location='lower right'))
plt.axis('off')
plt.savefig(f'../figures/map_lineages_clinical.pdf')
plt.close('all')


print('OVERALL',overall_fracs)
# now break down by quarter
clin_df['quarter'] = clin_df['Onset date'].dt.to_period('Q')
fig,ax1 = plt.subplots(1,4)
all_quarters = clin_df['quarter'].unique()
all_quarters = clin_df['quarter'].unique()[clin_df['quarter'].unique().argsort()]

for jj, q in enumerate(all_quarters):
    #### same, but now broken down by QE
    ax1[jj] = covered_gdf_r.plot(facecolor='papayawhip',edgecolor='silver',ax=ax1[jj])
    clin_df_q = clin_df[clin_df['quarter']==q]
    overall = clin_df_q.groupby('City/district')['Genotype'].value_counts()
    overall_fracs = overall/overall.groupby(level=0).sum()
    print(q,overall_fracs,overall)
    dists = overall_fracs.reset_index(level=1).index
    for jj_,ind in enumerate(dists):
        r0 = []
        for gt in all_gts:
            if (ind,gt) in overall_fracs.index:
                r0.append(overall_fracs.loc[(ind,gt)])
            else:
                r0.append(0)
        plot_pie(seq_sites.loc[district_lookup_r[ind]].x,seq_sites.loc[district_lookup_r[ind]].y,r0,ax1[jj],scalar=overall.loc[ind].sum()/2)

    if jj==0:
        # adding in the key
        df_key = pd.DataFrame([[1,-26.15, 28.739551]],columns=['count','Latitude','Longitude'])
        df_key_ = gpd.GeoDataFrame(df_key, geometry=gpd.points_from_xy(df_key.Longitude, df_key.Latitude), crs="EPSG:4326")
        df_key_ = df_key_.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0))
        df_key['geometry'] = df_key_
        from shapely.geometry import Point
        for j,v in enumerate([5,25]):
            newRow = df_key.iloc[j]
            newRow['geometry'] = Point(newRow['geometry'].x,newRow['geometry'].y-1.5E4)
            newRow['count'] = v
            df_key = pd.concat((df_key,pd.DataFrame([newRow])),axis=0)

        df_key = df_key.reset_index()
        for j in df_key.index:
            plot_pie(df_key.loc[j,'geometry'].x,df_key.loc[j,'geometry'].y,1,ax1[jj],edgecolor='black',c0='white',scalar=df_key.loc[j,'count']/2)
            ax1[jj].text(df_key.loc[j,'geometry'].x+10000,df_key.loc[j,'geometry'].y,df_key.loc[j,'count'],color='black',horizontalalignment='left',verticalalignment='center')
    labels = []
    if jj==0:
            ax1[jj].annotate(
            'Clinical',
            xy=(-0.02, 0.5),          # position left of subplot
            xycoords='axes fraction',
            ha='right',
            va='center',
            fontsize=13,
            fontweight='bold',
            rotation=90
            )
    if jj==3:
        ax1[jj].add_artist(ScaleBar(1,box_alpha=0,location='lower right'))
        # ax1[jj].legend(loc='upper left')
    ax1[jj].set_axis_off()
    ax1[jj].set_title(f"Q{q.quarter}\n{q.year}")

fig.subplots_adjust(
    left=0.0,
    right=1.0,
    top=1.,
    bottom=0.0,
    wspace=-0.08,   # ↓ column gap
    hspace=-0.42    # ↓ row gap
)      
fig.savefig(f'../figures/map_lineages_all_quarters_clinical.pdf')


