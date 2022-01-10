#!/usr/bin/env python
# coding: utf-8

# In[1]:


#loading libraries/modules and setting up the visualization details, not all of these are needed for using this ocape function


import numpy as np
import xarray as xr
from xmitgcm import open_mdsdataset
from MITgcmutils import jmd95 as jmd95
from scipy.optimize import linear_sum_assignment
import scipy as scipy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from prompt_toolkit import prompt

plt.rc('font', family='serif')
plt.rc('mathtext',**{'default':'regular'})
plt.rc('axes', linewidth=1.5)
get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('config', "InlineBackend.figure_format='retina'")


# In[2]:


def calrho(S,T,z): #just a function for calling the density based on Jackett and McDougall (1995)
    rho=jmd95.densjmd95(S,T,z)
    return rho


# In[3]:


#we calculate the salinity change at the interface using Newton's method
steps = 20000
high_sal = 35
low_sal = 34
average_sal = (high_sal + low_sal) / 2
analytical_dens = jmd95.densjmd95(34.47, -1.6, 300) + 3e-3 
compute_dens = jmd95.densjmd95(average_sal, 0.9, 300)
for k in range (0,steps):
        #calculate density from the linear step given in the Su et. al paper
        #and then calculate with average of high and low salinity until the average is close to the density change given in paper
    average_sal = (high_sal + low_sal) / 2
    compute_dens = jmd95.densjmd95(average_sal, 0.9, 300)
    if compute_dens < analytical_dens:
        low_sal = average_sal
    else:
        high_sal = average_sal
    if k == (steps-1):
        print(str(average_sal))


# In[4]:


#now, we calculate the salinity profiles in Case 2 from Table 1 of Su et. al

epsilon = .00001 #psu
sal_depth_arr = []

for i in range (0,701):
    if i == 0:
        sal_depth_arr.append(34.650137077196746)
    else:
        change_sal_depth = (1e-7) / (9.8*((1/calrho(sal_depth_arr[i-1],-1.6,300+(i)))* (calrho(sal_depth_arr[i-1]+epsilon,-1.6,300+i)-calrho(sal_depth_arr[i-1]-epsilon,-1.6,300+i))/(2*epsilon)))
        sal_depth_arr.append(sal_depth_arr[i-1]+(change_sal_depth))


# In[13]:


#finally, we make arrays of the depths, thetas, and salinities
sal_case_2 = np.empty([1000])
theta_case_2 = np.empty([1000])
depth_case_2 = np.empty([1000])


for i in range(0,301):
    sal_case_2[i] = 34.47
    depth_case_2[i] = i
    theta_case_2[i] = -1.6

for b in range(301,1001):
    sal_case_2[b-1] = sal_depth_arr[b-300]
    depth_case_2[b-1] = b
    theta_case_2[b-1] = 0.9


# In[6]:


#now, all that is left is to call the function
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

#run the scripts for the southern ocean


# In[15]:


ocape(100,theta_case_2, depth_case_2, sal_case_2)


# In[ ]:




