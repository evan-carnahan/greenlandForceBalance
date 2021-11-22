#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 14:36:13 2021

@author: evan
"""
#%% import dependencies
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from matplotlib_scalebar.scalebar import ScaleBar
mpl.rcParams.update({'font.size':10})

#%%
cwLand = gpd.read_file('glac3LandMask.shp')

rnkT = pd.read_csv('rnkTerm.csv',names=['yr','T'])
umiT = pd.read_csv('umiTerm.csv',names=['yr','T'])
ingT = pd.read_csv('ingTerm.csv',names=['yr','T'])

rnkT = rnkT[rnkT.yr > 1985]
rnkT['T'] = rnkT['T'] - rnkT['T'].iloc[0]
umiT = umiT[umiT.yr > 1985]
umiT['T'] = umiT['T'] - umiT['T'].iloc[0]
ingT = ingT[ingT.yr > 1985]
ingT['T'] = ingT['T'] - ingT['T'].iloc[0]


glac = 'umi'
yr1 = 1985
vdb85 = gpd.read_file('vdbShp/'+glac+'Vdb'+str(yr1)+'.shp')
yr1 = 2002
vdb02 = gpd.read_file('vdbShp/'+glac+'Vdb'+str(yr1)+'.shp')
yr1 = 2007
vdb07 = gpd.read_file('vdbShp/'+glac+'Vdb'+str(yr1)+'.shp')
yr1 = 2015
vdb15 = gpd.read_file('vdbShp/'+glac+'Vdb'+str(yr1)+'.shp')
umi8502 = gpd.sjoin(vdb85,vdb02)
umi0207 = gpd.sjoin(vdb02,vdb07)
umi0715 = gpd.sjoin(vdb07,vdb15)

# 02 - 85
umi8502['dz'] = (umi8502['z_right'] - umi8502['z_left'])/(2002-1985)
umi0207['dz'] = (umi0207['z_right'] - umi0207['z_left'])/(2007-2002)
umi0715['dz'] = (umi0715['z_right'] - umi0715['z_left'])/(2013.5-2007)
umi8502['dzDyn'] = (umi8502['z_right'] - umi8502['z_left'] - umi8502['dzSmb'])/(2002-1985)


umi85 = vdb85
umi85['Hab'] = (vdb85.z - vdb85.zBed) + 1028/917*vdb85.zBed



glac = 'ing'
yr1 = 1985
vdb85 = gpd.read_file('vdbShp/'+glac+'Vdb'+str(yr1)+'.shp')
yr1 = 2003
vdb02 = gpd.read_file('vdbShp/'+glac+'Vdb'+str(yr1)+'.shp')
yr1 = 2007
vdb07 = gpd.read_file('vdbShp/'+glac+'Vdb'+str(yr1)+'.shp')
yr1 = 2015
vdb15 = gpd.read_file('vdbShp/'+glac+'Vdb'+str(yr1)+'.shp')
ing8502 = gpd.sjoin(vdb85,vdb02)
ing0207 = gpd.sjoin(vdb02,vdb07)
ing0715 = gpd.sjoin(vdb07,vdb15)

# 02 - 85
ing8502['dz'] = (ing8502['z_right'] - ing8502['z_left'])/(2003-1985)
ing0207['dz'] = (ing0207['z_right'] - ing0207['z_left'])/(2007-2003)
ing0715['dz'] = (ing0715['z_right'] - ing0715['z_left'])/(2013.5-2007)
ing8502['dzDyn'] = (ing8502['z_right'] - ing8502['z_left'] - ing8502['dzSmb'])/(2003-1985)

ing85 = vdb85
ing85['Hab'] = (vdb85.z - vdb85.zBed) + 1028/917*vdb85.zBed


glac = 'rnk'
yr1 = 1985
vdb85 = gpd.read_file('vdbShp/'+glac+'Vdb'+str(yr1)+'.shp')
yr1 = 2001
vdb02 = gpd.read_file('vdbShp/'+glac+'Vdb'+str(yr1)+'.shp')
yr1 = 2007
vdb07 = gpd.read_file('vdbShp/'+glac+'Vdb'+str(yr1)+'.shp')
yr1 = 2015
vdb15 = gpd.read_file('vdbShp/'+glac+'Vdb'+str(yr1)+'.shp')

rnk8502 = gpd.sjoin(vdb85,vdb02)
rnk0207 = gpd.sjoin(vdb02,vdb07)
rnk0715 = gpd.sjoin(vdb07,vdb15)
# 02 - 85
rnk8502['dz'] = (rnk8502['z_right'] - rnk8502['z_left'])/(2001-1985)
rnk0207['dz'] = (rnk0207['z_right'] - rnk0207['z_left'])/(2007-2001)
rnk0715['dz'] = (rnk0715['z_right'] - rnk0715['z_left'])/(2013.5-2007)
rnk8502['dzDyn'] = (rnk8502['z_right'] - rnk8502['z_left'] - rnk8502['dzSmb'])/(2001-1985)

rnk85 = vdb85
rnk85['Hab'] = (vdb85.z - vdb85.zBed) + 1028/917*vdb85.zBed


habMin = min(rnk85['Hab'].min(),ing85['Hab'].min(),umi85['Hab'].min())
habMax = min(rnk85['Hab'].max(),ing85['Hab'].max(),umi85['Hab'].max())
 
vAbs = 15 
habMax = 300
#%%
fig, ax = plt.subplots(2,3,figsize = (7.5,6.27))

wid = 0.3
xBr = 0.05
ht = 0.5

i = 0
j = 0
ax[j,i].set_position([0.005+0.33*i,0.5,0.33,0.49])
scAx = cwLand.plot(ax = ax[j,i],color = 'g')
ing8502.plot('dz',ax = ax[j,i],markersize = 1,cmap = 'coolwarm_r',vmin = -vAbs,vmax = vAbs)
umi8502.plot('dz',ax = ax[j,i],markersize = 1,cmap = 'coolwarm_r',vmin = -vAbs,vmax = vAbs)
rnk8502.plot('dz',ax = ax[j,i],markersize = 1,cmap = 'coolwarm_r',vmin = -vAbs,vmax = vAbs)
ax[j,i].set_xticks(ticks = [])
ax[j,i].set_yticks(ticks = [])
scAx.add_artist(ScaleBar(1))

x, y, arrow_length = 0.87, 0.86, 0.2
scAx.annotate(' ',xy=(x+0.07, y), xytext=(x, y-arrow_length),
            arrowprops=dict(facecolor='black', width=2, headwidth=10),
            ha='center', va='center', fontsize=14, fontweight = 'bold',
            xycoords=scAx.transAxes) 
scAx.annotate('N', xy=(x, y - arrow_length), rotation = -13,
              ha='center', va='center', fontsize=14, fontweight = 'bold',
            xycoords=scAx.transAxes)


scAx.add_artist(ScaleBar(1,border_pad = 0.1))

i = 1
cwLand.plot(ax = ax[j,i],color = 'g')
ing0207.plot('dz',ax = ax[j,i],markersize = 1,cmap = 'coolwarm_r',vmin = -vAbs,vmax = vAbs)
umi0207.plot('dz',ax = ax[j,i],markersize = 1,cmap = 'coolwarm_r',vmin = -vAbs,vmax = vAbs)
rnk0207.plot('dz',ax = ax[j,i],markersize = 1,cmap = 'coolwarm_r',vmin = -vAbs,vmax = vAbs)
ax[j,i].set_xticks(ticks = [])
ax[j,i].set_yticks(ticks = [])
ax[j,i].set_position([0.005+0.33*i,0.5,0.33,0.49])

i = 2
cwLand.plot(ax = ax[j,i],color = 'g')
ing0715.plot('dz',ax = ax[j,i],markersize = 1,cmap = 'coolwarm_r',vmin = -vAbs,vmax = vAbs)
umi0715.plot('dz',ax = ax[j,i],markersize = 1,cmap = 'coolwarm_r',vmin = -vAbs,vmax = vAbs)
rnk0715.plot('dz',ax = ax[j,i],markersize = 1,cmap = 'coolwarm_r',vmin = -vAbs,vmax = vAbs)
ax[j,i].set_xticks(ticks = [])
ax[j,i].set_yticks(ticks = [])
ax[j,i].set_position([0.005+0.33*i,0.5-0.49*j,0.33,0.49])


j = 1
i = 0
cwLand.plot(ax = ax[j,i],color = 'g')
ing8502.plot('dzDyn',ax = ax[j,i],markersize = 1,cmap = 'coolwarm_r',vmin = -vAbs,vmax = vAbs)
umi8502.plot('dzDyn',ax = ax[j,i],markersize = 1,cmap = 'coolwarm_r',vmin = -vAbs,vmax = vAbs)
rnk8502.plot('dzDyn',ax = ax[j,i],markersize = 1,cmap = 'coolwarm_r',vmin = -vAbs,vmax = vAbs)
ax[j,i].set_xticks(ticks = [])
ax[j,i].set_yticks(ticks = [])
ax[j,i].set_position([0.005+0.33*i,0.01,0.33,0.49])

i = 1
cwLand.plot(ax = ax[j,i],color = 'g')
ing85.plot('Hab',ax = ax[j,i],markersize = 1,cmap = 'cividis_r',vmin = 0,vmax = habMax)
umi85.plot('Hab',ax = ax[j,i],markersize = 1,cmap = 'cividis_r',vmin=0,vmax = habMax)
rnk85.plot('Hab',ax = ax[j,i],markersize = 1,cmap = 'cividis_r',vmin = 0,vmax = habMax)
ax[j,i].set_xticks(ticks = [])
ax[j,i].set_yticks(ticks = [])
ax[j,i].set_position([0.005+0.33*i,0.01,0.33,0.49])

i = 2
ax[j,i].plot(ingT['yr'],ingT['T']/1e3,'r',label = 'Ing')
ax[j,i].plot(umiT['yr'],umiT['T']/1e3,'g',label = 'Umi')
ax[j,i].plot(rnkT['yr'],rnkT['T']/1e3,'b',label = 'Rnk')
ax[j,i].set_xlabel('Year')
ax[j,i].set_ylabel('Retreat, km')
ax[j,i].fill([2007,2007,2008,2008],[max(ingT['T'])/1e3 + 0.3,min(ingT['T'])/1e3,min(ingT['T'])/1e3,max(ingT['T'])/1e3 + 0.3],'k',alpha = 0.25, 
             label = 'DEM years')
ax[j,i].fill([2001,2001,2004,2004],[max(ingT['T'])/1e3 + 0.3,min(ingT['T'])/1e3,min(ingT['T'])/1e3,max(ingT['T'])/1e3 + 0.3],'k',alpha = 0.25)
ax[j,i].fill([2014,2014,2015,2015],[max(ingT['T'])/1e3 + 0.3,min(ingT['T'])/1e3,min(ingT['T'])/1e3,max(ingT['T'])/1e3 + 0.3],'k',alpha = 0.25)
ax[j,i].fill([1985,1985,1986,1986],[max(ingT['T'])/1e3 + 0.3,min(ingT['T'])/1e3,min(ingT['T'])/1e3,max(ingT['T'])/1e3 + 0.3],'k',alpha = 0.25)
ax[j,i].set_xlim([1985,2015])
ax[j,i].set_ylim([-7,max(ingT['T'])/1e3 + 0.3])
ax[j,i].legend(loc = 'lower left',frameon= False)

ax[j,i].set_position([0.67+0.055,0.08,0.265,0.4])
ax[j,i].set_yticklabels(('7','6','5','4','3','2','1','0'))

cbar_ax = fig.add_axes([0.9, 0.665, 0.01, 0.3])
sm = plt.cm.ScalarMappable(cmap='coolwarm_r', norm=plt.Normalize(vmin=-vAbs, vmax=vAbs))
fig.colorbar(sm, cax=cbar_ax,label = 'rate of elevation change, m/yr')


cbar_ax2 = fig.add_axes([0.57, 0.17, 0.01, 0.3])
sm2 = plt.cm.ScalarMappable(cmap='cividis_r', norm=plt.Normalize(vmin = 0,vmax = habMax))
fig.colorbar(sm2, cax=cbar_ax2,label = 'height above bouyancy, m')

ax[j,i].annotate('Ing', xy = (0.18,0.95), xycoords = 'figure fraction',fontsize = 14,fontweight = 'bold')
ax[j,i].annotate('Umi', xy = (0.22,0.77), xycoords = 'figure fraction',fontsize = 14,fontweight = 'bold')
ax[j,i].annotate('Rnk', xy = (0.25,0.65), xycoords = 'figure fraction',fontsize = 14,fontweight = 'bold')


ax[j,i].annotate('a', xy = (0.02,0.96), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
ax[j,i].annotate('b', xy = (0.018+0.33,0.96), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
ax[j,i].annotate('c', xy = (0.015+0.66,0.96), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
ax[j,i].annotate('d', xy = (0.02,0.47), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
ax[j,i].annotate('e', xy = (0.018+0.33,0.47), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
ax[j,i].annotate('f', xy = (0.015+0.66,0.47), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')

