#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import xarray as xr
from xmitgcm import open_mdsdataset
from MITgcmutils import jmd95 as jmd95
from scipy.optimize import linear_sum_assignment
import scipy as scipy
from prompt_toolkit import prompt

def ocape(par_amnt,theta_file,depth_file,sal_file):
    
    interpdepth=np.linspace(0,np.max(depth_file),par_amnt)
   
    interptheta=np.interp(interpdepth, depth_file, theta_file)
    
    interpsalt=np.interp(interpdepth, depth_file, sal_file)
    
    def jmd95buo(S,T,z):
        rhoConst=1000

        buo = 9.8 * ((rhoConst * (1 / jmd95.densjmd95(S, T, z) )) - 1)
        return buo

    def jmd95enth(S,T,z):
        enthalpy=np.nansum(jmd95buo(S, T, z)*(interpdepth[2]-interpdepth[1]))
        return enthalpy

    enthalpy_array=np.zeros([par_amnt,par_amnt])
    enthalpy_current=np.zeros([par_amnt])

    for k in np.arange(0,len(interpdepth)):
        enthalpy_current[k]=jmd95enth(interpsalt[k], interptheta[k], interpdepth[:k])
    
    for i in np.arange(0,len(interpdepth)):
        for j in np.arange(0,len(interpdepth)):
            enthalpy_array[i,j]= jmd95enth(interpsalt[i], interptheta[i], interpdepth[:j])
    
    row_ind, col_ind = linear_sum_assignment(enthalpy_array)
    
    return (enthalpy_current-enthalpy_array[row_ind, col_ind]).sum()/(par_amnt)

