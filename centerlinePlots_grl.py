#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 14:26:19 2021

@author: evan
"""

#%% dependencies 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def add_interval(ax, xdata, ydata, c, caps="  ",lSty = "-"):
    line = ax.add_line(mpl.lines.Line2D(xdata, ydata, color = c,linewidth = 1.5,linestyle = lSty))
    anno_args = {
        'ha': 'center',
        'va': 'center',
        'size': 8,
        'color': line.get_color()
    }
    a0 = ax.annotate(caps[0], xy=(xdata[0], ydata[0]), **anno_args)
    a1 = ax.annotate(caps[1], xy=(xdata[1], ydata[1]), **anno_args)
    return (line,(a0,a1))

fs = 10
glac = 'rnk'#'ing'#'umi'#

df = pd.read_csv(glac + 'ClAllInfo.csv')


if glac == 'ing' or glac == 'umi':
    xBd = [-0.5,35]
    indScl = 0
elif glac == 'rnk':
    xBd = [-0.5,25]
    indScl = 75

#%% define stress coupling length
scl = np.mean((df.z85 - df.bed)[indScl])*4/1e3

termInd85 = np.argwhere(~np.isnan(df.Td85.values))[0][0]
scl85Log = np.logical_and(df.dis > df.dis[termInd85],df.dis < df.dis[termInd85] + scl)
abScl85Log = df.dis > df.dis[termInd85] + scl

termInd02 = np.argwhere(~np.isnan(df.Td02.values))[0][0]
scl02Log = np.logical_and(df.dis > df.dis[termInd02],df.dis < df.dis[termInd02] + scl)
abScl02Log = df.dis > df.dis[termInd02] + scl

termInd07 = np.argwhere(~np.isnan(df.Td07.values))[0][0]
scl07Log = np.logical_and(df.dis > df.dis[termInd07],df.dis < df.dis[termInd07] + scl)
abScl07Log = df.dis > df.dis[termInd07] + scl

termInd15 = np.argwhere(~np.isnan(df.Td15.values))[0][0]
scl15Log = np.logical_and(df.dis > df.dis[termInd15],df.dis < df.dis[termInd15] + scl)
abScl15Log = df.dis > df.dis[termInd15] + scl

df['scl85'] = scl85Log
df['scl02'] = scl02Log
df['scl07'] = scl07Log
df['scl15'] = scl15Log

#%% ingia plot

if glac == 'ing':
    yrs = ['85','02','07','15']
    clr = ['r','k','g','c']
    labels = ['1985','2003','2007','2015']
    
    fig, ax1 = plt.subplots(5,1,figsize = (8,10))
    ax2 = ax1
    ax3 = ax1
    i = 0
    ax1[i].plot(df.dis,df.bed,'b')
    for j,iYr in enumerate(yrs):
        ax1[i].plot(df.dis,df['z'+iYr],clr[j])
    ax1[i].set_ylabel('z, m', color='k')
    ax1[i].tick_params('y', colors='k')
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    ax1[i].set_xticklabels([])
    
    ax1[i].plot([df.dis[df.t85],df.dis[df.t85]],[df.bed[df.t85],df.z85[df.t85]],'r')
    ax1[i].plot([df.dis[df.t02],df.dis[df.t02]],[df.bed[df.t02],df.z02[df.t02]],'k')
    ax1[i].plot([df.dis[df.t07],df.dis[df.t07]],[df.bed[df.t07],df.z07[df.t07]],'g')
    ax1[i].plot([df.dis[df.t15],df.dis[df.t15]],[df.bed[df.t15],df.z15[df.t15]],'c')
    
    
    ax2[i] = ax1[i].twinx()
    for j,iYr in enumerate(yrs):
        ax2[i].plot(df.dis,df['v' + iYr],clr[j],label = labels[j])
    ax2[i].set_ylabel('vel, km/yr', color='k')
    ax2[i].tick_params('y', colors='k')
    ax2[i].legend(fontsize = fs)
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    ax2[i].set_ylim([0,1.75])
    ax1[i].set_xlim(xBd)
    
    
    i = 1
    for j,iYr in enumerate(yrs):
        ax1[i].plot(df.dis,df['Td' + iYr],clr[j])
        add_interval(ax1[i], (df.dis[df['scl' + iYr]].iloc[0],df.dis[df['scl' + iYr]].iloc[-1]),
                  (np.mean(df['Td' + iYr][df['scl' + iYr]]),np.mean(df['Td' + iYr][df['scl' + iYr]])), clr[j] ,"[]")
    ax1[i].set_ylim([0,300])
    ax1[i].set_ylabel('driving stress, KPa', color='k')
    ax1[i].tick_params('y', colors='k')
    ax1[i].set_xlim(xBd)
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    ax1[i].set_xticklabels([])
    
    
    i = 2
    for j,iYr in enumerate(yrs):
        ax1[i].plot(df.dis,df['Tb' + iYr],clr[j])
        add_interval(ax1[i], (df.dis[df['scl' + iYr]].iloc[1],df.dis[df['scl' + iYr]].iloc[-1]),
                 (np.mean(df['Tb' + iYr][df['scl' + iYr]]),np.mean(df['Tb' + iYr][df['scl' + iYr]])), clr[j] ,"[]")
    ax1[i].set_ylim([-65,225])
    ax1[i].set_ylabel('basal drag, KPa', color='k')
    ax1[i].tick_params('y', colors='k')
    ax1[i].set_xlim(xBd)
    ax1[i].set_xticklabels([])
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    
    
    i = 3
    for j,iYr in enumerate(yrs):
        ax1[i].plot(df.dis,df['TLat' + iYr],clr[j])
        add_interval(ax1[i], (df.dis[df['scl' + iYr]].iloc[1],df.dis[df['scl' + iYr]].iloc[-1]),
                 (np.mean(df['TLat' + iYr][df['scl' + iYr]]),np.mean(df['TLat' + iYr][df['scl' + iYr]])), clr[j] ,"[]")
    ax1[i].set_ylim([-20,110])
    ax1[i].set_ylabel('Lateral drag, KPa', color='k')
    ax1[i].tick_params('y', colors='k')
    ax1[i].set_xlim(xBd)
    ax1[i].set_xticklabels([])
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    
    
    i = 4
    for j,iYr in enumerate(yrs):
        ax1[i].plot(df.dis,df['TLon' + iYr],clr[j])
        add_interval(ax1[i], (df.dis[df['scl' + iYr]].iloc[1],df.dis[df['scl' + iYr]].iloc[-1]),
                  (np.mean(df['TLon' + iYr][df['scl' + iYr]]),np.mean(df['TLon' + iYr][df['scl' + iYr]])), clr[j] ,"[]")
    ax1[i].set_ylabel('long. coupling, KPa', color='k')
    ax1[i].set_xlim(xBd)
    ax1[i].set_ylim([-30,40])
    ax1[i].set_xlabel('Flowline distance, km')
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    
    
    
    ax3[i] = ax1[i].twinx()
    ax3[i].plot(df.dis,df.fTLon85,'r--',alpha = 0.5)
    ax3[i].plot(df.dis,df.fTLon02,'k--',label = 'percent resistance',alpha = 0.5)
    ax3[i].plot(df.dis,df.fTLon07,'g--',alpha = 0.5)
    ax3[i].plot(df.dis,df.fTLon15,'c--',alpha = 0.5)
    for j,iYr in enumerate(yrs):
        add_interval(ax3[i], (df.dis[df['scl' + iYr]].iloc[1],df.dis[df['scl' + iYr]].iloc[-1]),
              (np.mean(df['fTLon' + iYr][df['scl' + iYr]]),np.mean(df['fTLon' + iYr][df['scl' + iYr]])), clr[j] ,"[]",lSty = "--")
    ax3[i].set_ylabel('driving stress resistance, %', color='k')
    ax3[i].tick_params('y', colors='k')
    ax3[i].legend(fontsize = fs,loc = 'lower right')
    fig.tight_layout()
    
    
    ax1[i].annotate('a', xy = (0.01,0.97), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
    ax1[i].annotate('b', xy = (0.01,0.78), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
    ax1[i].annotate('c', xy = (0.01,0.58), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
    ax1[i].annotate('d', xy = (0.01,0.39), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
    ax1[i].annotate('e', xy = (0.01,0.2), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
    
    ax1[i].annotate('Retreat onset: 2002', xy = (0.81,0.832), xycoords = 'figure fraction',fontsize = 10,fontweight = 'bold',ha = 'center')
    


#%% umiamako plot
elif glac == 'umi':
    yrs = ['85','02','07','15']
    clr = ['r','k','g','c']
    labels = ['1985','2002','2007','2015']
    
    fig, ax1 = plt.subplots(5,1,figsize = (8,10))
    ax2 = ax1
    ax3 = ax1
    i = 0
    ax1[i].plot(df.dis,df.bed,'b')
    for j,iYr in enumerate(yrs):
        ax1[i].plot(df.dis,df['z'+iYr],clr[j])
    ax1[i].set_ylabel('z, m', color='k')
    ax1[i].tick_params('y', colors='k')
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    ax1[i].set_xticklabels([])
    
    ax1[i].plot([df.dis[df.t85],df.dis[df.t85]],[df.bed[df.t85],df.z85[df.t85]],'r')
    ax1[i].plot([df.dis[df.t02],df.dis[df.t02]],[df.bed[df.t02],df.z02[df.t02]],'k')
    ax1[i].plot([df.dis[df.t07],df.dis[df.t07]],[df.bed[df.t07],df.z07[df.t07]],'g')
    ax1[i].plot([df.dis[df.t15],df.dis[df.t15]],[df.bed[df.t15],df.z15[df.t15]],'c')
    
    
    ax2[i] = ax1[i].twinx()
    for j,iYr in enumerate(yrs):
        ax2[i].plot(df.dis,df['v' + iYr],clr[j],label = labels[j])
    ax2[i].set_ylabel('vel, km/yr', color='k')
    ax2[i].tick_params('y', colors='k')
    ax2[i].legend(fontsize = fs)
    ax2[i].set_ylim([0,1.9])
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    ax1[i].set_xlim(xBd)
    
    i = 1
    for j,iYr in enumerate(yrs):
        ax1[i].plot(df.dis,df['Td' + iYr],clr[j])
        
        add_interval(ax1[i], (df.dis[df['scl' + iYr]].iloc[0],df.dis[df['scl' + iYr]].iloc[-1]),
                 (np.mean(df['Td' + iYr][df['scl' + iYr]]),np.mean(df['Td' + iYr][df['scl' + iYr]])), clr[j] ,"[]")
        
    ax1[i].set_ylabel('driving stress, KPa', color='k')
    ax1[i].tick_params('y', colors='k')
    ax1[i].set_xlim(xBd)
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    ax1[i].set_xticklabels([])
    
    
    i = 2
    for j,iYr in enumerate(yrs):
        ax1[i].plot(df.dis,df['Tb' + iYr],clr[j])
        add_interval(ax1[i], (df.dis[df['scl' + iYr]].iloc[0],df.dis[df['scl' + iYr]].iloc[-1]),
                 (np.mean(df['Tb' + iYr][df['scl' + iYr]]),np.mean(df['Tb' + iYr][df['scl' + iYr]])), clr[j] ,"[]")
    ax1[i].set_ylabel('basal drag, KPa', color='k')
    ax1[i].tick_params('y', colors='k')
    ax1[i].set_xlim(xBd)
    ax1[i].set_xticklabels([])
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    
    
    i = 3
    for j,iYr in enumerate(yrs):
        ax1[i].plot(df.dis,df['TLat' + iYr],clr[j])
        add_interval(ax1[i], (df.dis[df['scl' + iYr]].iloc[0],df.dis[df['scl' + iYr]].iloc[-1]),
                 (np.mean(df['TLat' + iYr][df['scl' + iYr]]),np.mean(df['TLat' + iYr][df['scl' + iYr]])), clr[j] ,"[]")
    ax1[i].set_ylabel('Lateral drag, KPa', color='k')
    ax1[i].tick_params('y', colors='k')
    ax1[i].set_xlim(xBd)
    ax1[i].set_xticklabels([])
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    
    
    i = 4
    for j,iYr in enumerate(yrs):
        ax1[i].plot(df.dis,df['TLon' + iYr],clr[j])
        add_interval(ax1[i], (df.dis[df['scl' + iYr]].iloc[1],df.dis[df['scl' + iYr]].iloc[-1]),
                  (np.mean(df['TLon' + iYr][df['scl' + iYr]]),np.mean(df['TLon' + iYr][df['scl' + iYr]])), clr[j] ,"[]")
    ax1[i].set_ylabel('long. coupling, KPa', color='k')
    ax1[i].set_xlim(xBd)
    ax1[i].set_xlabel('Flowline distance, km')
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    
    
    
    ax3[i] = ax1[i].twinx()
    ax3[i].plot(df.dis,df.fTLon85,'r--',alpha = 0.5)
    ax3[i].plot(df.dis,df.fTLon02,'k--',label = 'percent resistance',alpha = 0.5)
    ax3[i].plot(df.dis,df.fTLon07,'g--',alpha = 0.5)
    ax3[i].plot(df.dis,df.fTLon15,'c--',alpha = 0.5)
    for j,iYr in enumerate(yrs):
        add_interval(ax3[i], (df.dis[df['scl' + iYr]].iloc[0],df.dis[df['scl' + iYr]].iloc[-1]),
              (np.mean(df['fTLon' + iYr][df['scl' + iYr]]),np.mean(df['fTLon' + iYr][df['scl' + iYr]])), clr[j] ,"[]",lSty = "--")
    ax3[i].set_ylim([-29,121])
    ax3[i].set_ylabel('driving stress resistance, %', color='k')
    ax3[i].tick_params('y', colors='k')
    ax3[i].legend(fontsize = fs,loc = 'lower right')
    fig.tight_layout()
    
    
    ax1[i].annotate('a', xy = (0.01,0.97), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
    ax1[i].annotate('b', xy = (0.01,0.78), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
    ax1[i].annotate('c', xy = (0.01,0.58), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
    ax1[i].annotate('d', xy = (0.01,0.39), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
    ax1[i].annotate('e', xy = (0.01,0.2), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
    
    ax1[i].annotate('Retreat: 2001 - 2009', xy = (0.805,0.825), xycoords = 'figure fraction',fontsize = 10,fontweight = 'bold',ha = 'center')


#%% rnk plot
elif glac == 'rnk':
    yrs = ['85','02','07','15']
    clr = ['r','k','g','c']
    labels = ['1985','2001','2007','2015']
    
    fig, ax1 = plt.subplots(5,1,figsize = (8,10))
    ax2 = ax1
    ax3 = ax1
    i = 0
    ax1[i].plot(df.dis,df.bed,'b')
    for j,iYr in enumerate(yrs):
        ax1[i].plot(df.dis,df['z'+iYr],clr[j])
    ax1[i].set_ylabel('z, m', color='k')
    ax1[i].tick_params('y', colors='k')
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    ax1[i].set_xticklabels([])
    
    ax1[i].plot([df.dis[df.t85],df.dis[df.t85]],[df.bed[indScl],df.z85[df.t85]],'r')
    ax1[i].plot([df.dis[df.t02],df.dis[df.t02]],[df.bed[df.t02],df.z02[df.t02]],'k')
    ax1[i].plot([df.dis[df.t07],df.dis[df.t07]],[df.bed[df.t07],df.z07[df.t07]],'g')
    ax1[i].plot([df.dis[df.t15],df.dis[df.t15]],[df.bed[df.t15],df.z15[df.t15]],'c')
    
    
    ax2[i] = ax1[i].twinx()
    for j,iYr in enumerate(yrs):
        ax2[i].plot(df.dis,df['v' + iYr],clr[j],label = labels[j])
    ax2[i].set_ylabel('vel, km/yr', color='k')
    ax2[i].tick_params('y', colors='k')
    ax2[i].legend(fontsize = fs)
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    ax1[i].set_xlim(xBd)
    
    i = 1
    for j,iYr in enumerate(yrs):
        ax1[i].plot(df.dis,df['Td' + iYr],clr[j])
        add_interval(ax1[i], (df.dis[df['scl' + iYr]].iloc[0],df.dis[df['scl' + iYr]].iloc[-1]),
                 (np.mean(df['Td' + iYr][df['scl' + iYr]]),np.mean(df['Td' + iYr][df['scl' + iYr]])), clr[j] ,"[]")
    ax1[i].set_ylabel('driving stress, KPa', color='k')
    ax1[i].tick_params('y', colors='k')
    ax1[i].set_xlim(xBd)
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    ax1[i].set_xticklabels([])
    
    
    i = 2
    for j,iYr in enumerate(yrs):
        ax1[i].plot(df.dis,df['Tb' + iYr],clr[j])
        add_interval(ax1[i], (df.dis[df['scl' + iYr]].iloc[0],df.dis[df['scl' + iYr]].iloc[-1]),
                 (np.mean(df['Tb' + iYr][df['scl' + iYr]]),np.mean(df['Tb' + iYr][df['scl' + iYr]])), clr[j] ,"[]")
    ax1[i].set_ylabel('basal drag, KPa', color='k')
    ax1[i].tick_params('y', colors='k')
    ax1[i].set_xlim(xBd)
    ax1[i].set_xticklabels([])
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    
    
    i = 3
    for j,iYr in enumerate(yrs):
        ax1[i].plot(df.dis,df['TLat' + iYr],clr[j])
        add_interval(ax1[i], (df.dis[df['scl' + iYr]].iloc[0],df.dis[df['scl' + iYr]].iloc[-1]),
                 (np.mean(df['TLat' + iYr][df['scl' + iYr]]),np.mean(df['TLat' + iYr][df['scl' + iYr]])), clr[j] ,"[]")
    ax1[i].set_ylabel('Lateral drag, KPa', color='k')
    ax1[i].tick_params('y', colors='k')
    ax1[i].set_xlim(xBd)
    ax1[i].set_xticklabels([])
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    
    
    i = 4
    for j,iYr in enumerate(yrs):
        ax1[i].plot(df.dis,df['TLon' + iYr],clr[j])
        add_interval(ax1[i], (df.dis[df['scl' + iYr]].iloc[1],df.dis[df['scl' + iYr]].iloc[-1]),
                  (np.mean(df['TLon' + iYr][df['scl' + iYr]]),np.mean(df['TLon' + iYr][df['scl' + iYr]])), clr[j] ,"[]")
    ax1[i].set_ylabel('long. coupling, KPa', color='k')
    ax1[i].set_xlim(xBd)
    ax1[i].set_xlabel('Flowline distance, km')
    ax1[i].tick_params(bottom = True,top = True,direction = 'in')
    
    
    
    ax3[i] = ax1[i].twinx()
    ax3[i].plot(df.dis,df.fTLon85,'r--',alpha = 0.5)
    ax3[i].plot(df.dis,df.fTLon02,'k--',label = 'percent resistance',alpha = 0.5)
    ax3[i].plot(df.dis,df.fTLon07,'g--',alpha = 0.5)
    ax3[i].plot(df.dis,df.fTLon15,'c--',alpha = 0.5)
    for j,iYr in enumerate(yrs):
        add_interval(ax3[i], (df.dis[df['scl' + iYr]].iloc[0],df.dis[df['scl' + iYr]].iloc[-1]),
              (np.mean(df['fTLon' + iYr][df['scl' + iYr]]),np.mean(df['fTLon' + iYr][df['scl' + iYr]])), clr[j] ,"[]")
    ax3[i].set_ylim([-20,80])
    ax3[i].set_ylabel('driving stress resistance, %', color='k')
    ax3[i].tick_params('y', colors='k')
    ax3[i].legend(fontsize = fs,loc = 'lower right')
    fig.tight_layout()
    
    
    ax1[i].annotate('a', xy = (0.01,0.97), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
    ax1[i].annotate('b', xy = (0.01,0.78), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
    ax1[i].annotate('c', xy = (0.01,0.58), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
    ax1[i].annotate('d', xy = (0.01,0.39), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
    ax1[i].annotate('e', xy = (0.01,0.2), xycoords = 'figure fraction',fontsize = 12,fontweight = 'bold')
