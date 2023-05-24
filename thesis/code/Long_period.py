#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 16:35:22 2023

@author: adri
"""


from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle

plt.rcParams.update({'font.size': 20})

hdulist = fits.open('EPOCH_PHOTOMETRY-Gaia DR3 3005222633152619904.fits')
hdulist.info()
dat = hdulist[1].data

b=dat.field('band')
mag = dat.field("mag")
time=dat.field("time")
flux=dat.field("flux")

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
        f_g.append(flux[i])
        
    if (b[i] =="BP"):
        BP.append(mag[i])
        t_bp.append(time[i])
        f_bp.append(flux[i])
        
    if (b[i] =="RP"):
        RP.append(mag[i])
        t_rp.append(time[i])
        f_rp.append(flux[i])
        







# unfolded flux

# erg.s-1.cm-2.µm-1
# e-/s electron charge per second

plt.figure(figsize=(15, 10)) 
plt.scatter(t_rp,f_rp,color='red',marker = 'o',  label = 'RP band')
plt.scatter(t_g,f_g,color='green',marker = 'o',  label = 'G band')
plt.scatter(t_bp,f_bp,color='blue',marker = 'o',  label = 'BP band')
#for ii in range(0,10):
#    plt.axvline(t_ref + p1*ii, c='C3', alpha=0.25)
#plt.text(x = 2200, y = 11.2, s = "Long period variable Gaia DR3 3005222633152619904", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
plt.legend(fontsize = 15, loc = "lower right")
plt.xlabel("Time (BJD)")
plt.ylabel('Flux (e-/s)')
#plt.ylim(13,11)
#plt.xlim(2200,2400)
plt.title('Flux time series')
plt.show()



# unfolded magnitude

plt.figure(figsize=(15, 10)) 
plt.scatter(t_rp,RP,color='red',marker = 'o',  label = 'RP band')
plt.scatter(t_g,G,color='green',marker = 'o',  label = 'G band')
plt.scatter(t_bp,BP,color='blue',marker = 'o',  label = 'BP band')
#for ii in range(0,10):
#    plt.axvline(t_ref + p1*ii, c='C3', alpha=0.25)
plt.text(x = 2200, y = 11.1, s = "Long period variable Gaia DR3 3005222633152619904", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
plt.legend(fontsize = 15, loc = "lower right")
plt.xlabel("Time (BJD)")
plt.ylabel('Magnitude (mag)')
plt.ylim(13,11)
#plt.xlim(2200,2400)
plt.title('Epoch photometry')
plt.show()




# folded magnitude 

i=0
P=[100]    # guess 

for i in range(len(P)):
    p1 = P[i]
    
    plt.figure(figsize=(15, 10)) 

    t_rp=np.array(t_rp)      
    phase =  ((t_rp - 1720.93511 - p1/2) / p1) % 1 #- 0.5
    plt.scatter(phase,RP,s= 80,color='red',marker = 'o',  label = 'RP band') 

    t_g=np.array(t_g)        
    phase =  ((t_g - 1720.93468 - p1/2) / p1) % 1  #- 0.5
    plt.scatter(phase,G,s= 80,color='green',marker = 'o',  label = 'G band') 

    t_bp=np.array(t_bp)    
    phase =  ((t_bp - 1720.93503 - p1/2) / p1) % 1  #- 0.5
    plt.scatter(phase,BP,s= 80,color='blue',marker = 'o',  label = 'BP band') 

    plt.text(x = 0.5, y = 11, s = "Long period variable Gaia DR3 3005222633152619904", horizontalalignment='center', color = 'gray')
    plt.gca().invert_yaxis()
#plt.xlim(-0.4,0.4)
    plt.ylim(13.2,10.8)
#plt.ylim(11.85,11.75) # G band
#plt.ylim(11.5,11.2) # RP band
#plt.axhline(11.37, c ='C2', alpha=1)
    plt.legend(fontsize = 15, loc = "lower right")
    plt.xlabel("Phase")
    plt.ylabel('Magnitude (mag)')
    plt.title('Epoch photometry')
    plt.show()


# vizier epoch photometry: cdsarc.cds.unistra.fr//vizier/vizgraph.gml?-s=I/355&-i=.graph_sql_epphot&Pos=092.35626437760-09.64698864270&Source=3005222633152619904
