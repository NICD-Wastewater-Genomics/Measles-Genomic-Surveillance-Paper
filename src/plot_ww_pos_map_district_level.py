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
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['font.sans-serif'] = ['Arial']

# cases = pd.read_csv('../assets/province_case_counts_2025.csv',index_col='Province')
#counts by Province
vl_df = pd.read_excel('../metadata/Wastewater metadata_20250218_SG.xlsx',sheet_name="Sample Testing Results")
# vl_df = pd.read_csv('../metadata/wastewater_positivity.csv')
vl_df.columns = vl_df.columns.str.strip()
vl_df = vl_df[vl_df['Measles Result']!='Failed']

vl_province = pd.DataFrame(vl_df.groupby('Site Province')['Measles Result'].value_counts())
vl_province = vl_province.pivot_table(index=['Site Province'], columns='Measles Result',values='count',aggfunc='sum',fill_value=0)
vl_province['Total'] = vl_province['Negative'] + vl_province['Positive'] 
vl_province['Positivity'] = vl_province['Positive'] /vl_province['Total'] 

# vl_df = vl_df[vl_df['Site Province']=='Gauteng']
vl_district = pd.DataFrame(vl_df.groupby(['Site Province','District Name'])['Measles Result'].value_counts())
vl_district = vl_district.pivot_table(index=['District Name'], columns='Measles Result',values='count',aggfunc='sum',fill_value=0)
vl_district['Total'] = vl_district['Negative'] + vl_district['Positive'] 
vl_district['Positivity'] = vl_district['Positive'] /vl_district['Total'] 

print('TOTAL SAMPLE COUNT:', vl_district.Total.sum())
print('TOTAL POSITIVE COUNT:', vl_district.Positive.sum())
print('AVERAGE TPR:', vl_district.Positive.sum()/vl_district.Total.sum())

local_district_map ={'Rustenburg Local Municipality':'Bojanala Platinum DM',
                     'Sol Plaatjie Local Municipality':'Frances Baard DM',
                     'Musina Local Municipality':'Vhembe DM'}
for l in local_district_map.keys():
    if l in vl_district.index:
        vl_district.loc[l] = vl_district.loc[l] + vl_district.loc[l]
        vl_district = vl_district.drop(index=l)

district_lookup= {"City of Johannesburg":'Johannesburg MM',
                  "City of Tshwane":"Tshwane MM",
                  "Ekurhuleni":"Ekurhuleni MM",
                  "Bojanala": "Bojanala Platinum DM",
                  "Buffalo City":"Buffalo City MM",
                  "City of Cape Town":"Cape Town MM",
                  "Ehlanzeni":"Ehlanzeni DM",
                  "eThekwini":"Ethekwini MM",
                  "Frances Baard":"Frances Baard DM",
                  "Mangaung":"Mangaung MM",
                  "Nelson Mandela Bay":"Nelson Mandela Bay MM",
                  "Ngaka Modiri Molema":"Ngaka Modiri Molema DM",
                  "Umkhanyakude":"Umkhanyakude DM",
                  "Vhembe":"Vhembe DM",
                  "o r tambo":"o.r.tambo",
                  "thabo mofutsanyana":"thabo mofutsanyane"}
print('missing districts: ',sum(~vl_district.index.isin(district_lookup.values())))
vl_district = vl_district.loc[vl_district.index.isin(district_lookup.values())]

### then add country map
gdf = gpd.read_file("../assets/map_files/zaf_admbnda_adm1_sadb_ocha_20201109.shp").set_index("ADM1_EN")
gdf_Gauteng = gdf[gdf.index=='Gauteng']
gdf_others = gdf[gdf.index!='Gauteng']
# convert map to projected WGS 84. 
gdfG_r = gdf_Gauteng.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0)).simplify(20)
gdfO_r = gdf_others = gdf[gdf.index!='Gauteng'].to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0)).simplify(20)
#fix northern cape spelling
gdfO_r.loc['Northern Cape'] = gdfO_r.loc['Nothern Cape'] 
gdfO_r = gdfO_r[gdfO_r.index !='Nothern Cape'] 


# norm = matplotlib.colors.Normalize(vmin=0, vmax=vl_province['Positivity'].max()) 
# create a scalarmappable from the colormap

base = plt.cm.Blues
colors = base(np.linspace(0.1, 1, 256))  # skip pale colors
cmap = mcolors.LinearSegmentedColormap.from_list("Blues_trunc", colors)
norm = matplotlib.colors.Normalize(vmin=0, vmax=vl_district['Positivity'].max())


fig,ax = plt.subplots()
#Gauteng first, with red outline
gdfG_r.plot(facecolor=cmap(norm(vl_province.loc['Gauteng','Positivity'])),edgecolor='red',ax=ax,zorder=100000)#,label='Sequencing Data')
for province in gdfO_r.index:
    gdfO_r.loc[[province]].plot(facecolor=cmap(norm(vl_province.loc[province,'Positivity'])),edgecolor='silver',ax=ax)#,label='Sequencing Data')

# ax.add_artist(ScaleBar(1,box_alpha=0,location='lower right'))
# ax.legend(loc='upper left')
plt.axis('off')
plt.savefig('../figures/ww_positivity_map_provinces.pdf')
plt.close('all')

vl_province.to_csv('ww_positivity_provinces.csv')
### now district map, just for Gauteng. 

# cases = pd.read_csv('../assets/gauteng_district_case_counts_2025.csv',index_col='District')

gdf = gpd.read_file("../assets/map_files/zaf_admbnda_adm2_sadb_ocha_20201109.shp").set_index("ADM2_EN")

## wastewater data, and city locations

meta_ww = pd.read_csv('../metadata/MeV Wastewater Sequences Final Metadata_05092025.csv',index_col='Sample Number')

#read in location of each wwtp
df = pd.read_csv('../assets/All sampling sites_v5_300502024.csv',skipinitialspace=True)
df['Latitude'] = df['GPS latitude']
df['Longitude'] = df['GPS longitude']
meta_ww[~meta_ww['Site Name'].isin(df['Site'])]
hh = df[~df['Site'].isin(meta_ww['Site Name'])]
meta_ww['Site Name'] = meta_ww['Site Name'].str.lower()
vl_df['Site Name'] = vl_df['Site Name'].str.lower()
df['Site'] = df['Site'].str.lower()
df = df[df['Site'].isin(vl_df['Site Name'])]
df = df[df.Type.str.contains('National')]


sites = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude), crs="EPSG:4326")
seq_sites = sites#[sites['Seq']=="Yes"]
seq_sites_r = seq_sites.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0))
#################

# convert map to projected WGS 84. 
gdf_r = gdf.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0)).simplify(20)

#counts by Province
cases_2024 = pd.read_csv('../metadata/MEV_clinical_numbers_2024.csv',index_col='Districts')
cases_2024 = cases_2024.drop(index=['Grand Total'],columns=['Grand Total']).fillna(0)
cases_2024['Cases'] = cases_2024.sum(axis=1)
cases_2025 = pd.read_csv('../metadata/MEV_clinical_numbers_2025.csv',index_col='Districts')
cases_2025 = cases_2025.drop(index=['Grand Total'],columns=['Grand Total']).fillna(0)
cases_2025['Cases'] = cases_2025.sum(axis=1)

cases = pd.concat((cases_2024['Cases'],cases_2025['Cases']),axis=1).fillna(0)
cases['All_Cases'] = cases.sum(axis=1)
cases.index = [ci.lower().replace(' metro','').replace(' platinum','') for ci in cases.index]


district_lookup_= {"City of Johannesburg":'Johannesburg MM',
                  "City of Tshwane":"Tshwane MM",
                  "Ekurhuleni":"Ekurhuleni MM",
                  "Bojanala": "Bojanala Platinum DM",
                  "Buffalo City":"Buffalo City MM",
                  "City of Cape Town":"Cape Town MM",
                  "Ehlanzeni":"Ehlanzeni DM",
                  "eThekwini":"Ethekwini MM",
                  "Frances Baard":"Frances Baard DM",
                  "Mangaung":"Mangaung MM",
                  "Nelson Mandela Bay":"Nelson Mandela Bay MM",
                  "Ngaka Modiri Molema":"Ngaka Modiri Molema DM",
                  "Umkhanyakude":"Umkhanyakude DM",
                  "Vhembe":"Vhembe DM",
                  "o.r.tambo":"o r tambo",
                  "thabo mofutsanyane":"thabo mofutsanyana",
                  "eden":"garden route",
                  "sisonke":"harry gwala",
                  "uthungulu":"king cetshwayo",
                  "cacadu":"sarah baartman",
                  "z f mgcawu":"zf mgcawu"}
district_lookup_ = {k.lower(): v for k, v in district_lookup_.items()}
district_lookup_flip_ = {value: key for key, value in district_lookup_.items()}

vl_district_ = vl_district.copy()
vl_district_.index = [district_lookup_flip_[ci] for ci in vl_district_.index]


all_dat = pd.merge(cases,vl_district_,left_index=True,right_index=True)
import seaborn as sns

# generate Fig 1c
fig, ax  = plt.subplots(figsize=(3.,4))
ax.set_xlim([1,1000])
ax.set_ylim([0,0.15])
all_dat = all_dat[all_dat['Positivity']>0]
sns.regplot(data=all_dat, x='All_Cases', y='Positivity',ax=ax,logx=True, scatter_kws={'clip_on':False,"color": "black"},truncate=False, line_kws={"color": "grey"})
# ax.scatter(all_dat['All_Cases'],all_dat['Positivity'])
ax.set_xlabel('Cases per district')
ax.set_ylabel('Wastewater TPR per district')
ax.spines[['right', 'top']].set_visible(False)
ax.set_xscale('log')
fig.tight_layout()
plt.savefig('../figures/case_count_vs_ww_positivity.pdf')
plt.close('all')

import scipy.stats as stats
r, p = stats.spearmanr(np.log10(all_dat['All_Cases']),all_dat['Positivity'])
r2 = r**2
print('spearman p: ', r)
print('R2: ',r2)
print("P value:", p)



## generate Fig 1b, now with sampling sites shown
fig,ax = plt.subplots()

gdfG_r.plot(facecolor='none',edgecolor='red',ax=ax,zorder=100000)#,label='Sequencing Data')
for province in gdfO_r.index:
    gdfO_r.loc[[province]].plot(facecolor='none',edgecolor='black',ax=ax,zorder=9999)#,label='Sequencing Data')

for district in gdf_r.index:
    if district in vl_district.index:
        gdf_r.loc[[district]].plot(facecolor=plt.cm.Blues(norm(vl_district.loc[district,'Positivity'])),edgecolor='grey',ax=ax,linewidth=0.8,zorder=100)#,label='Sequencing Data')
        print(district,vl_district.loc[district,'Positive'],vl_district.loc[district,'Total'],vl_district.loc[district,'Positivity'])
        vl_district = vl_district.drop(index=district)
    elif district in district_lookup.keys():
        district0 = district_lookup[district]
        gdf_r.loc[[district]].plot(facecolor=plt.cm.Blues(norm(vl_district.loc[district0,'Positivity'])),edgecolor='grey',ax=ax,linewidth=0.8,zorder=100)#,label='Sequencing Data')
        print(district,vl_district.loc[district0,'Positive'],vl_district.loc[district0,'Total'],vl_district.loc[district0,'Positivity'])
        vl_district = vl_district.drop(index=district0)
    else:
        gdf_r.loc[[district]].plot(facecolor='none',edgecolor='silver',ax=ax,linewidth=0.1)#,label='Sequencing Data')
        print(district)
seq_sites_r.plot(ax=ax,edgecolor='black',facecolor='white',markersize=5,zorder=10000000)

ax.add_artist(ScaleBar(1,box_alpha=0,location='lower right'))
# ax.legend(loc='upper left')
plt.axis('off')
plt.savefig('../figures/ww_positivity_map_district_level.pdf')
plt.close('all')

###

gdf= gdf[gdf['ADM1_EN']=='Gauteng']
gdf_r = gdf.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0)).simplify(20)

df = df[df.Province=='Gauteng']
df = df[df.Type.str.contains('National')]

sites = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude), crs="EPSG:4326")
seq_sites = sites#[sites['Seq']=="Yes"]
seq_sites_r = seq_sites.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0))

## generate ED Fig 3a

fig,ax = plt.subplots()
gdfG_r.plot(facecolor='none',edgecolor='red',ax=ax,zorder=100000)#,label='Sequencing Data')
for district in ['City of Johannesburg','City of Tshwane','Ekurhuleni','Sedibeng','West Rand']:
    # district0 = district_lookup[district]
    if district.lower() in vl_district_.index:
        gdf_r.loc[[district]].plot(facecolor=plt.cm.Blues(norm(vl_district_.loc[district.lower(),'Positivity'])),edgecolor='silver',ax=ax)#,label='Sequencing Data')
    else:
        gdf_r.loc[[district]].plot(facecolor='none',edgecolor='silver',ax=ax)#,label='Sequencing Data')
seq_sites_r.plot(ax=ax,edgecolor='black',facecolor='white',markersize=20)

ax.add_artist(ScaleBar(1,box_alpha=0,location='lower right'))
# ax.legend(loc='upper left')
plt.axis('off')
plt.savefig('../figures/ww_positivity_map_gauteng.pdf')
plt.close('all')

## generate ED Figs 3b-d
fig, ax = plt.subplots(figsize=(0.5, 3), layout='constrained')

# Colorbar
cbar = fig.colorbar(
    plt.cm.ScalarMappable(norm=norm, cmap=cmap),
    cax=ax,
    orientation='vertical',
    label='WW Positivity'
)
cbar.ax.tick_params(labelsize=14)
cbar.ax.yaxis.label.set_size(18)
# --- No data patch ---
nodata_patch = mpatches.Patch(
    facecolor="white",
    edgecolor="black",
    label="No data"
)

# Place legend below colorbar
fig.legend(
    handles=[nodata_patch],
    loc="lower left",
    bbox_to_anchor=(0.05, -0.05),
    frameon=False
)
cbar.ax.yaxis.set_label_position('left')
fig.savefig('../figures/ww_colorbar_map_districts.pdf',transparent=True)

vl_df['District Name'] = [district_lookup_flip_[ci] for ci in vl_df['District Name']]

igm = pd.read_csv('../metadata/IgM stats_080626.csv')
igm.columns = igm.columns.str.strip()
igm['Date'] = pd.to_datetime({'year': igm['Year'], 'month': igm['Month'], 'day': 1})

igm['District'] = [id0.split(' METRO')[0] for id0 in igm['District']]

igm['District'] = igm['District'].str.lower()
vl_df = vl_df[vl_df['Site Province']=='Gauteng']
districts = vl_df['District Name'].unique()

district_data = []
for district in districts:
    vl_loc = vl_df[vl_df['District Name']==district]
    vl_loc = vl_loc.set_index('Sample Collection Date')
    vc = vl_loc.groupby(pd.Grouper(freq='MS'))['Measles Result'].value_counts().unstack(fill_value=0)
    pos_fraction = vc.get('Positive', 0) / vc.sum(axis=1)

    igm_loc = igm[igm['District']==district]
    merged = pd.merge(
        pos_fraction.rename('WW Positivity').reset_index().rename(columns={'Sample Collection Date': 'Date'}),
        igm_loc[['Date', 'IgM Positives', 'Total Tests']],
        on='Date', how='outer'
    ).sort_values('Date')

    merged['Clin Positivity'] = merged['IgM Positives'] / merged['Total Tests']

    paired = merged[['WW Positivity', 'IgM Positives']].dropna()
    if len(paired) >= 2:
        r, p = stats.spearmanr(paired['WW Positivity'], paired['IgM Positives'])
    else:
        r, p = float('nan'), float('nan')

    district_data.append((district, merged, r, p))

fig, axes = plt.subplots(len(districts), 1, figsize=(8, 3 * len(districts)), sharex=True)
if len(districts) == 1:
    axes = [axes]

for i, (district, merged, r, p) in enumerate(district_data[::-1]):
    ax1 = axes[i]
    ax2 = ax1.twinx()

    ax1.plot(merged['Date'], merged['WW Positivity'], '-o',color='steelblue', label='Wastewater TPR',clip_on=False)
    ax2.plot(merged['Date'], merged['IgM Positives'].fillna(0), '-o',color='firebrick', label='Clinical IgM Positives',clip_on=False)

    ax1.set_ylabel('Wastewater TPR', color='steelblue')
    ax2.set_ylabel('Clinical IgM Positives', color='firebrick')
    ax1.tick_params(axis='y', labelcolor='steelblue')
    ax2.tick_params(axis='y', labelcolor='firebrick')
    locator = mdates.MonthLocator(bymonthday=1)
    ax1.xaxis.set_major_locator(locator)
    ax1.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

    ax1.set_title(f"{district}  (ρ = {r:.2f}, p = {p:.3f})")
    ax1.tick_params(labelbottom=True)
    ax1.set_xlim(pd.to_datetime(['2023-12-15','2025-08-15']))
    ax1.set_ylim([0, 0.4])
    ax2.set_ylim([0, merged['IgM Positives'].max() + 2])

plt.tight_layout()
plt.savefig('../figures/ww_clin_positivity_districts.pdf')
plt.close('all')
