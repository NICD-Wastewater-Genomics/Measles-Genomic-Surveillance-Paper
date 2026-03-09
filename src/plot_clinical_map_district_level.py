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
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

# Figures produced here:
# Fig 1a

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['savefig.bbox'] = 'tight'

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

gauteng_cases = cases.loc[['city of johannesburg','city of tshwane','ekurhuleni','sedibeng','west rand'],'All_Cases'].sum()
fraction_gauteng = gauteng_cases/cases['All_Cases'].sum()
print("Fraction of cases in Gauteng: ",fraction_gauteng)
gdf = gpd.read_file("../assets/map_files/zaf_admbnda_adm1_sadb_ocha_20201109.shp").set_index("ADM1_EN")
gdf_Gauteng = gdf[gdf.index=='Gauteng']
gdf_others = gdf[gdf.index!='Gauteng']
# convert map to projected WGS 84. 
gdfG_r = gdf_Gauteng.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0)).simplify(20)
gdfO_r = gdf_others = gdf[gdf.index!='Gauteng'].to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0)).simplify(20)
#fix northern cape spelling
gdfO_r.loc['Northern Cape'] = gdfO_r.loc['Nothern Cape'] 
gdfO_r = gdfO_r[gdfO_r.index !='Nothern Cape'] 


cmap = plt.cm.Reds
norm = matplotlib.colors.Normalize(vmin=0, vmax=cases['All_Cases'].max()) 
## some differences in spelling across metadata an map. 
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
                  "o.r.tambo":"o r tambo",
                  "thabo mofutsanyane":"thabo mofutsanyana",
                  "eden":"garden route",
                  "sisonke":"harry gwala",
                  "uthungulu":"king cetshwayo",
                  "cacadu":"sarah baartman",
                  "z f mgcawu":"zf mgcawu"}
district_lookup = {k.lower(): v for k, v in district_lookup.items()}

gdf = gpd.read_file("../assets/map_files/zaf_admbnda_adm2_sadb_ocha_20201109.shp").set_index("ADM2_EN")
gdf_r = gdf.to_crs("EPSG:2053").scale(xfact=-1,yfact=-1,origin=(0,0)).simplify(20)
gdf_r.index = [gri.lower() for gri in gdf_r.index]


fig,ax = plt.subplots()
gdfG_r.plot(facecolor='none',edgecolor='red',ax=ax,zorder=100000)
for province in gdfO_r.index:
    gdfO_r.loc[[province]].plot(facecolor='none',edgecolor='black',ax=ax,zorder=9999)


for district in gdf_r.index:
    if district in cases.index:
        gdf_r.loc[[district]].plot(facecolor=cmap(norm(cases.loc[district,'All_Cases'])),edgecolor='black',ax=ax)
        cases = cases.drop(index=district)
    elif district in district_lookup.keys():
        district0 = district_lookup[district]
        gdf_r.loc[[district]].plot(facecolor=cmap(norm(cases.loc[district0,'All_Cases'])),edgecolor='black',ax=ax)
        cases = cases.drop(index=district0)
    else:
        gdf_r.loc[[district]].plot(facecolor='none',edgecolor='silver',ax=ax,linewidth=0.1)
        # check to make sure we haven't left out any districts that should be in there. 
        print(district)

ax.add_artist(ScaleBar(1,box_alpha=0,location='lower right'))
plt.axis('off')
plt.savefig('../figures/map_all_district_cases.pdf')
plt.close('all')


fig, ax = plt.subplots(figsize=(0.5, 3), layout='constrained')
cbar = fig.colorbar(
    plt.cm.ScalarMappable(norm=norm, cmap=cmap),
    cax=ax,
    orientation='vertical',
    label='Cases'
)

cbar.ax.yaxis.set_label_position('left')
cbar.ax.tick_params(labelsize=14)
cbar.ax.yaxis.label.set_size(18)
fig.savefig('../figures/colorbar_map_districts_all.pdf',transparent=True)