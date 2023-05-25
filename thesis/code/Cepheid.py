#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 21:10:14 2023

@author: adri
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, join
from astropy.utils.data import get_pkg_data_filename
plt.rcParams.update({'font.size': 20})


##### epoch photometry ######

hdulist = fits.open('EPOCH_PHOTOMETRY-Gaia DR3 4301612233202961024.fits')
hdulist.info()
dat = hdulist[1].data

b=dat.field('band')
mag = dat.field("mag")
time=dat.field("time")


p1 = 0.21859 #  Period corresponding to the first overtone pulsation mode in the G band time series 


G=[]
BP=[]
RP=[]
t_g =[]
t_bp=[]
t_rp=[]
f_g=[]
f_bp=[]
f_rp=[]

i=0
for i in range(len(b)):
    if (b[i] =="G"):
        G.append(mag[i])
        t_g.append(time[i])
       
    if (b[i] =="BP"):
        BP.append(mag[i])
        t_bp.append(time[i])
        
    if (b[i] =="RP"):
        RP.append(mag[i])
        t_rp.append(time[i])
      





# unfolded magnitude

plt.figure(figsize=(15, 10)) 
plt.scatter(t_rp,RP,color='red',marker = 'o',  label = 'RP band')
plt.scatter(t_g,G,color='green',marker = 'o',  label = 'G band')
plt.scatter(t_bp,BP,color='blue',marker = 'o',  label = 'BP band')
plt.text(x = 2200, y = 9.6, s = "Cepheid Gaia DR3 4301612233202961024", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
plt.legend(fontsize = 15,loc = "best")
plt.xlabel("Time (BJD)")
plt.ylabel('Magnitude (mag)')
plt.ylim(11,9.4)
#plt.xlim(2200,2400)
plt.title('Epoch photometry')
plt.show()




# folded magnitude 


plt.figure(figsize=(15, 10)) 

t_rp=np.array(t_rp)      
phase =  ((t_rp - 1725.73221 - p1) / p1/2) % 1 - 0.5
plt.scatter(phase,RP,s= 80,color='red',marker = 'o',  label = 'RP band') 

t_g=np.array(t_g)        
phase =  ((t_g - 1725.73182 - p1) / p1/2) % 1  - 0.5
plt.scatter(phase,G,s= 80,color='green',marker = 'o',  label = 'G band') 

t_bp=np.array(t_bp)    
phase =  ((t_bp - 1725.73212 - p1) / p1/2) % 1 - 0.5
plt.scatter(phase,BP,s= 80,color='blue',marker = 'o',  label = 'BP band') 

plt.text(x = 0, y = 9.6, s = "Cepheid Gaia DR3 4301612233202961024", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
#plt.xlim(-0.4,0.4)
plt.ylim(11,9.4) # BP band
#plt.ylim(11.85,11.75) # G band
#plt.ylim(11.5,11.2) # RP band
#plt.axhline(11.37, c ='C2', alpha=1)
plt.legend(fontsize = 15, loc = "lower right")
plt.xlabel("Phase")
plt.ylabel('Magnitude (mag)')
plt.title('Epoch photometry')
plt.show()










