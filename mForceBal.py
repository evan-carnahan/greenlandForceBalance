#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 11:27:46 2019

@author: evan
"""
#module for claculating force balances for glaciers
import pandas as pd
import numpy as np
import math

#glac= 'umi'
#gridSp = 500
#yr = '1985'
#fbDf = pd.read_csv(glac + 'VelDemBed_grid'+str(gridSp)+'_yr' +yr+'.csv')

def strainCalc(glac,gridSp,yr,method,fol=''):
    
    dx = gridSp
    dy = gridSp
    #calculate strain rates
    #dx = 150
    #dy = 150
    n = 3 # glens flow law exponent
    B = 300*1e3 #Pa a^(1/3) for temperate ice ~ -5 ^oC
    if fol == '':
        fbDf = pd.read_csv(glac + 'VelDemBed_grid'+str(gridSp)+method+'_yr' +yr+'.csv')
    else:
        fbDf = pd.read_csv(fol+'/'+glac + 'VelDemBed_grid'+str(gridSp)+method+'_yr' +yr+'.csv')
    nameLis = ['dux_dx','duy_dy','dux_dy','duy_dx','eDot_xx','eDot_yy','eDot_xy','eDot_e','R_xx','R_yy','R_xy']
    for na in nameLis:
        fbDf[na] = np.nan
    
    for i in range(fbDf.shape[0]):
        x_i = fbDf.x[i]
        y_j = fbDf.y[i]
        rXP150 = fbDf.query('@x_i+@dx == x and @y_j == y')
        rXN150 = fbDf.query('@x_i-@dx == x and @y_j == y')
        rYP150 = fbDf.query('@x_i == x and @y_j+@dy == y')
        rYN150 = fbDf.query('@x_i == x and @y_j-@dy == y')
        # ignore boundary values
        if rYP150.empty or rYN150.empty or rXP150.empty or rXN150.empty:
            continue
        # normal strain
        dux_dx = (rXP150.iloc[0,fbDf.columns.get_loc('vX')] - rXN150.iloc[0,fbDf.columns.get_loc('vX')])/(2*dx)
        duy_dy = (rYP150.iloc[0,fbDf.columns.get_loc('vY')] - rYN150.iloc[0,fbDf.columns.get_loc('vY')])/(2*dy)
        # cross strain
        dux_dy = (rYP150.iloc[0,fbDf.columns.get_loc('vX')] - rYN150.iloc[0,fbDf.columns.get_loc('vX')])/(2*dy)
        duy_dx = (rXP150.iloc[0,fbDf.columns.get_loc('vY')] - rXN150.iloc[0,fbDf.columns.get_loc('vY')])/(2*dx)
        
        eDot_xx = dux_dx
        eDot_yy = duy_dy
        eDot_xy = 1/2*(dux_dy+duy_dx)
        
        # from O'Neel, 2005
#        eDot_zz = - eDot_xx - eDot_yy
#        eDot_e = (eDot_xx ** 2 + eDot_yy ** 2 + eDot_zz ** 2 + 2 * eDot_xy ** 2) ** (1/2)
        
        # vanDerVeen, 2013, should be (2 * (.)) ** (1/2)
        eDot_e = (eDot_xx ** 2 + eDot_yy **2 + eDot_xx * eDot_yy + eDot_xy ** 2) ** (1/2)
        
        R_xx = B*eDot_e ** (1/n-1) * (2*eDot_xx+eDot_yy)
        R_yy = B*eDot_e ** (1/n-1) * (eDot_xx+2*eDot_yy)
        R_xy = B*eDot_e ** (1/n-1) * eDot_xy
        varLis = [dux_dx,duy_dy,dux_dy,duy_dx,eDot_xx,eDot_yy,eDot_xy,eDot_e,R_xx,R_yy,R_xy]
        
        fbDf.loc[i,nameLis] = varLis
    fbDf.dropna(0,inplace=True)
    fbDf.reset_index(drop=True,inplace=True)
        
    fbDf.to_csv(glac + 'Strain_grid'+str(gridSp)+method+'_yr' +yr+'.csv',index=False)
    return fbDf


#%%
#dx = 150
#dy = 150
        
def stressCalc(glac,gridSp,yr,method):
    dx = dy = gridSp
    rho = 917 #kg/m^3
    g = 9.81 #m/s^2
    fbDf= pd.read_csv(glac + 'Strain_grid'+str(gridSp)+method+'_yr' +yr+'.csv')
    
    nameLis = ['t_dx','t_dy','t_latx','t_lonx','t_laty','t_lony','t_bx','t_by']
    for na in nameLis:
        fbDf[na] = np.nan
    
    for i in range(fbDf.shape[0]):
        x_i = fbDf.x[i]
        y_j = fbDf.y[i]
        rXP150 = fbDf.query('@x_i+@dx == x and @y_j == y')
        rXN150 = fbDf.query('@x_i-@dx == x and @y_j == y')
        rYP150 = fbDf.query('@x_i == x and @y_j+@dy == y')
        rYN150 = fbDf.query('@x_i == x and @y_j-@dy == y')
        row = fbDf.loc[i,:]
        if rYP150.empty or rYN150.empty or rXP150.empty or rXN150.empty:
            continue
        # x-direction forces
        t_dx = (- rho * g * row.H * (rXP150.h.values-rXN150.h.values)/(2*dx))[0]                
        t_lonx = ((rXP150.H.values * rXP150.R_xx.values - rXN150.H.values * rXN150.R_xx.values)/(2*dx))[0]
        t_latx = ((rYP150.H.values * rYP150.R_xy.values - rXN150.H.values * rYN150.R_xy.values)/(2*dy))[0]
                
        # y-direction forces
        t_dy = (- rho * g * row.H * (rYP150.h.values-rYN150.h.values)/(2*dy))[0]        
        t_lony = ((rYP150.H.values * rYP150.R_yy.values - rYN150.H.values * rYN150.R_yy.values)/(2*dy))[0]
        t_laty = ((rXP150.H.values * rXP150.R_xy.values - rXN150.H.values * rXN150.R_xy.values)/(2*dx))[0]
        
        # directional force balance 
        t_bx = t_dx+t_latx+t_lonx # negative lat and lon resist flow, positive acts in direction of flow
        t_by = t_dy+t_laty+t_lony
        varLis = [t_dx,t_dy,t_latx,t_lonx,t_laty,t_lony,t_bx,t_by]
        fbDf.loc[i,nameLis] = varLis
    
    fbDf.dropna(axis=0,inplace=True)
    fbDf.reset_index(drop=True,inplace=True)
    fbDf.to_csv(glac + 'Stress_gridXY'+str(gridSp)+method+'_yr' +yr+'.csv')
    return fbDf
    
    
#%% claculate on s,n coordinates
    
def stressOrient(glac,gridSp,yr,method):
    fbDf = pd.read_csv(glac + 'Stress_gridXY'+str(gridSp)+method+'_yr' +yr+'.csv')
    fbSn = fbDf[['x','y','vX','vY']]

    nameLis = ['t_ds','t_dn','t_lats','t_lons','t_latn','t_lonn','t_bs','t_bn']
    for na in nameLis:
        fbSn[na] = np.nan
    
    # role reversal in coordinates
    rads = np.arctan2(fbSn.vY.values,fbSn.vX.values)
    
    for i in range(fbSn.shape[0]):
        rad = rads[i]
        # make rotation matrix
        rotMat = np.array([[math.cos(rad),math.sin(rad)],[-math.sin(rad),math.cos(rad)]])
        
        t_dxy = np.array([[fbDf.t_dx[i]],[fbDf.t_dy[i]]])
        t_lonxy = np.array([[fbDf.t_latx[i]],[fbDf.t_laty[i]]])
        t_latxy = np.array([[fbDf.t_lonx[i]],[fbDf.t_lony[i]]])
        t_bxy = np.array([[fbDf.t_bx[i]],[fbDf.t_by[i]]])
        
        t_dsn = np.matmul(rotMat,t_dxy)
        t_bsn = np.matmul(rotMat,t_bxy)
        t_lonsn = np.matmul(rotMat,t_lonxy)
        t_latsn = np.matmul(rotMat,t_latxy)
        
        t_ds = t_dsn[0][0]
        t_dn = t_dsn[1][0]
        t_bs = t_bsn[0][0]
        t_bn = t_bsn[1][0]    
        t_lons = t_lonsn[0][0]
        t_lonn = t_lonsn[1][0]    
        t_lats = t_latsn[0][0]
        t_latn = t_latsn[1][0]
        
        varLis = [t_ds,t_dn,t_lats,t_lons,t_latn,t_lonn,t_bs,t_bn]
        fbSn.loc[i,nameLis] = varLis
    
    fbSn.to_csv(glac + 'Stress_gridSN'+str(gridSp)+method+'_yr' +yr+'.csv',index=False)
    return fbSn

#%% orient strain
def strainOrient(glac,gridSp,yr,method):
    fbDf = pd.read_csv(glac + 'Strain_grid'+str(gridSp)+method+'_yr' +yr+'.csv')
    fbSn = fbDf[['x','y','vX','vY','eDot_xx','eDot_yy','eDot_xy','eDot_e']]

    nameLis = ['eDot_ss','eDot_nn','eDot_sn']
    for na in nameLis:
        fbSn[na] = np.nan
    
    # role reversal in coordinates
    rads = np.arctan2(fbSn.vY.values,fbSn.vX.values)
    
    for i in range(fbSn.shape[0]):
        rad = rads[i]
        
        eDot_ss = fbDf['eDot_xx'][i]*math.cos(rad)**2 + fbDf['eDot_yy'][i]*math.sin(rad)**2 + 2*fbDf['eDot_xy'][i]*math.cos(rad)*math.sin(rad)
        eDot_nn = fbDf['eDot_xx'][i]*math.sin(rad)**2 + fbDf['eDot_yy'][i]*math.cos(rad)**2 - 2*fbDf['eDot_xy'][i]*math.cos(rad)*math.sin(rad)
        eDot_sn = (fbDf['eDot_xx'][i] - fbDf['eDot_yy'][i])*math.sin(rad)*math.cos(rad) + fbDf['eDot_xy'][i]*(math.cos(rad)**2 - math.sin(rad)**2)
        
        varLis = [eDot_ss,eDot_nn,eDot_sn]
        fbSn.loc[i,nameLis] = varLis
    
    fbSn.to_csv(glac + 'Strain_gridSN'+str(gridSp)+method+'_yr' +yr+'.csv',index=False)
    return fbSn
        
        

#%% full force balance pipeline 
def forceBalFull(glac,gridSp,yr,method,fol):
    dfStrain = strainCalc(glac,gridSp,yr,method,fol)
    dfStress = stressCalc(glac,gridSp,yr,method)
    dfStrainSn = strainOrient(glac,gridSp,yr,method)
    dfSn = stressOrient(glac,gridSp,yr,method)
    return dfStrain,dfStress,dfSn
    
    
    
if __name__ == "__main__":    

    glacs = ['ing','umi','rnk']
    for glac in glacs:
        yrs = [1985,2003,2007,2015] 
        if glac == 'umi':
          yrs[1] = 2002
        elif glac == 'rnk':
          yrs[1] = 2001
        for yr in yrs:
            forceBalFull(glac,250,str(yr),method='',fol='./vdbCsv')
#    strainOrient(glac,250,str(yr),method = '')
            
            















    
    