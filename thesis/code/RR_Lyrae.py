#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 14:06:49 2023

@author: adri
"""

### RR Lyrae stars

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle
import matplotlib.cm as cm
import matplotlib
import math

plt.rcParams.update({'font.size': 20})

### All-sky Aitoff projection in Galactic coordinates

hdulist5 = fits.open('rrlyrae1.fits')
hdulist5.info()
dat5 = hdulist5[1].data
l5=dat5.field('l')
b5=dat5.field('b')
l5 = Angle(l5 * u.deg)
l5.wrap_at('180d', inplace=True)
print(l5)
b5 = Angle(b5 * u.deg)
b5.wrap_at('90d', inplace=True)
print(b5)
l5 = np.deg2rad(np.array(l5))
b5 = np.deg2rad(np.array(b5))


hdulist6 = fits.open('rrlyrae2.fits')
hdulist6.info()
dat6 = hdulist6[1].data
l6=dat6.field('l')
b6=dat6.field('b')
l6 = Angle(l6 * u.deg)
l6.wrap_at('180d', inplace=True)
print(l6)
b6 = Angle(b6 * u.deg)
b6.wrap_at('90d', inplace=True)
print(b6)
l6 = np.deg2rad(np.array(l6))
b6 = np.deg2rad(np.array(b6))

fig = plt.figure(figsize=(15,12))
ax = fig.add_subplot(1,1,1, projection='aitoff')
ax.scatter(l5, b5, s=50, color='red', alpha=1, marker ='*', label = 'RR Lyrae in LOPN1')   
ax.scatter(l6, b6, s=50, color='blue', alpha=1,marker ='*', label = 'RR Lyrae in LOPS2')   
ax.legend(fontsize='15')
plt.title("Aitoff projection in Galactic coordinates\n")
plt.tight_layout()
plt.grid()


#####################

# Absolute Color-Magnitude Diagram



#rrlyrae
hdulist2 = fits.open('rrlyrae1.fits')
hdulist2.info()
dat2 = hdulist2[1].data

g2=dat2.field('median_mag_g_fov')
bp2=dat2.field('median_mag_bp')
rp2=dat2.field('median_mag_rp')
p2=dat2.field("parallax")
ag2 = dat2.field("ag_gspphot") 

print(g2)
print(bp2)
print(rp2)
print(p2)
print(ag2)

N = len(g2)

color2=np.zeros(N,float)
G2 = np.zeros(N,float)
i=0
for i in range(N):
        color2[i]=bp2[i]-rp2[i]
        G2[i] = g2[i]+5 *( 1+np.log10 (p2[i]/1000) ) - ag2[i]


print(G2)



plt.figure(figsize=(15, 10))

plt.scatter(color2,G2,color='red',marker = 'o',  label = 'RR Lyrae in LOPN1')


hdulist2 = fits.open('rrlyrae2.fits')
hdulist2.info()
dat2 = hdulist2[1].data

g2=dat2.field('median_mag_g_fov')
bp2=dat2.field('median_mag_bp')
rp2=dat2.field('median_mag_rp')
p2=dat2.field("parallax")
ag2=dat2.field("ag_gspphot")

print(g2)
print(bp2)
print(rp2)
print(p2)
print(ag2)

N = len(g2)

color2=np.zeros(N,float)
G2 = np.zeros(N,float)
i=0
for i in range(N):
        color2[i]=bp2[i]-rp2[i]
        G2[i] = g2[i]+5 *( 1+np.log10 (p2[i]/1000) ) - ag2[i]


print(G2)




plt.scatter(color2,G2,color='blue',marker = 'o',  label = 'RR Lyrae in LOPS2')
plt.legend(fontsize="15")
plt.gca().invert_yaxis() 
plt.xlabel('BP-RP')
plt.ylabel('G')
plt.title('aCMD')
#p.axis([0,1.4,26,18])
plt.show()






###############################

# period - amplitude diagram for rr lyrae stars


hdulist = fits.open('rrlyrae1.fits')

hdulist.info()
dat = hdulist[1].data

amp=dat.field('peak_to_peak_g')
pf=dat.field('pf')
p1=dat.field('p1_o')
print(amp,pf,p1)
p = pf
z=dat.field('metallicity')


i = 0
for i in range(len(amp)):
    if math.isnan(pf[i]):
        pf[i] = 0
        
    if math.isnan(p1[i]):
        p1[i] = 0        
   
i = 0    
for i in range(len(amp)):
    if (pf[i]==0):
        p[i] = p1[i]
    if (p1[i]==0): 
        p[i] = pf[i]
    else:
        p[i] = (pf[i] + p1[i])/2


print(p) 


plt.figure(figsize=(15, 10))
#plt.scatter(p,amp,color='red',marker = 'o',  label = 'RR Lyrae in LOPN1')
plt.scatter(p,amp,c=z,s = 50, cmap = 'viridis')





hdulist = fits.open('rrlyrae2.fits')

hdulist.info()
dat = hdulist[1].data

amp=dat.field('peak_to_peak_g')
pf=dat.field('pf')
p1=dat.field('p1_o')
print(amp,pf,p1)
p = pf
z=dat.field('metallicity')


i = 0
for i in range(len(amp)):
    if math.isnan(pf[i]):
        pf[i] = 0
        
    if math.isnan(p1[i]):
        p1[i] = 0        
   
i = 0    
for i in range(len(amp)):
    if (pf[i]==0):
        p[i] = p1[i]
    if (p1[i]==0): 
        p[i] = pf[i]
    else:
        p[i] = (pf[i] + p1[i])/2


print(p)    


#plt.scatter(p,amp,color='blue',marker = 'o',  label = 'RR Lyrae in LOPS2')
plt.scatter(p,amp,c=z,s = 50, cmap = 'viridis')
#plt.legend()
plt.xlabel('Period (d)')
plt.ylabel('Amp(G) (mag)')
plt.title('P-A')
clb=plt.colorbar()
#clb.ax.set_title('Metallicity (dex)')
clb.ax.set_ylabel('[Fe/H] (dex)')
plt.show()



#It results that the stars with longer periods (RRab) have lower metallicity 



#####################

'''
Metrics targeting non-periodic variations, such as
skewness_mag_g_fov vs abbe_mag_g_fov;
metrics targeting periodic variations of pulsating stars,
such as log 10 (std_dev_mag_bp/std_dev_mag_rp) versus
median_mag_g_fov.
'''

hdulist = fits.open('rrlyrae1.fits')
hdulist.info()
dat = hdulist[1].data

skew=dat.field('skewness_mag_g_fov')
abbe = dat.field("abbe_mag_g_fov")

plt.figure(figsize=(15, 10)) 
#plt.scatter(abbe,skew,s= 80, color='blue',marker = 'o') 
plt.scatter(abbe,skew,s= 80, c= abbe, cmap = 'plasma') 
plt.gca().invert_yaxis()
#plt.legend(fontsize = 15, loc = "best")
plt.xlabel("abbe_mag_g_fov")
plt.xlim(0.3,1.5)
#plt.gca().invert_yaxis()
plt.ylabel('skewness_mag_g_fov')
#plt.ylim(-2,0.75) 
plt.title('Metrics targeting non-periodic variations')
#plt.show()



hdulist = fits.open('rrlyrae2.fits')
hdulist.info()
dat = hdulist[1].data

skew=dat.field('skewness_mag_g_fov')
abbe = dat.field("abbe_mag_g_fov")

#plt.figure(figsize=(15, 10)) 
#plt.scatter(abbe,skew,s= 80, color='blue',marker = 'o') 
plt.scatter(abbe,skew,s= 80, c= abbe, cmap = 'plasma') 
plt.gca().invert_yaxis()
#plt.legend(fontsize = 15, loc = "best")
plt.xlabel("abbe_mag_g_fov")
plt.xlim(0.3,1.5)
#plt.gca().invert_yaxis()
plt.ylabel('skewness_mag_g_fov')
plt.ylim(0.75,-2) 
plt.title('Metrics targeting non-periodic variations')
plt.show()



hdulist = fits.open('rrlyrae1.fits')
hdulist.info()
dat = hdulist[1].data

g=dat.field('median_mag_g_fov')
s_bp=dat.field("std_dev_mag_bp")
s_rp=dat.field("std_dev_mag_rp")
y=np.log10(s_bp/s_rp)  ### ratio between the standard deviations in magnitude of FoV transit observations in the G BP and G RP bands

plt.figure(figsize=(15, 10)) 
#plt.scatter(abbe,skew,s= 80, color='blue',marker = 'o') 
plt.scatter(g,y,s= 80, c= g, cmap = 'plasma_r') 
#plt.legend(fontsize = 15, loc = "best")
plt.xlabel("median_mag_g_fov")
plt.gca().invert_xaxis()
plt.xlim(8.5,14)
plt.ylabel('std_dev_mag_bp/std_dev_mag_rp')
plt.ylim(0.05,0.35) 
plt.title('Metrics targeting periodic variations')
#plt.show()



hdulist = fits.open('rrlyrae2.fits')
hdulist.info()
dat = hdulist[1].data

g=dat.field('median_mag_g_fov')
s_bp=dat.field("std_dev_mag_bp")
s_rp=dat.field("std_dev_mag_rp")
y=np.log10(s_bp/s_rp)  ### ratio between the standard deviations in magnitude of FoV transit observations in the G BP and G RP bands

#plt.figure(figsize=(15, 10)) 
#plt.scatter(abbe,skew,s= 80, color='blue',marker = 'o') 
plt.scatter(g,y,s= 80, c= g, cmap = 'plasma_r') 
#plt.legend(fontsize = 15, loc = "best")
plt.xlabel("median_mag_g_fov")
plt.gca().invert_xaxis()
plt.xlim(8.5,14)
plt.ylabel('std_dev_mag_bp/std_dev_mag_rp')
plt.ylim(0.05,0.35) 
plt.title('Metrics targeting periodic variations')
plt.show()

######################

# EPOCH PHOTOMETRY FOR RR LYRAE IN THE G, G_BP, G_RP BAND

# source_id 1985008784705347456
# this is a RRc star pulsating with fundamental period p1 = 0.38031




hdulist = fits.open('EPOCH_PHOTOMETRY-Gaia DR3 1985008784705347456.fits')

hdulist.info()
dat = hdulist[1].data

p1 = 0.38031
b=dat.field('band')
mag = dat.field("mag")
time=dat.field("time")
G=[]
BP=[]
RP=[]
t_g =[]
t_bp=[]
t_rp=[]

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
        

# folded

        

plt.figure(figsize=(15, 10)) 

t_rp=np.array(t_rp)        # time of first transit in RP band 1689.14785
phase =  ((t_rp - 1689.14785 - p1/2) / p1/2) % 1 - 0.5
plt.scatter(phase,RP,s= 80,color='red',marker = 'o',  label = 'RP band') 
#plt.scatter(phase,RP,s= 80, c= phase, cmap = 'autumn',  label = 'RP band') 


t_g=np.array(t_g)        # time of first transit in G band 1689.14744
phase =  ((t_g - 1689.14744 - p1/2) / p1/2) % 1 - 0.5
plt.scatter(phase,G,s= 80,color='green',marker = 'o',  label = 'G band') 
#plt.scatter(phase,G,s= 80, c= phase, cmap = 'summer',  label = 'G band') 

t_bp=np.array(t_bp)        # time of first transit in BP band 1689.14776
phase =  ((t_bp - 1689.14776 - p1/2) / p1/2) % 1 - 0.5
plt.scatter(phase,BP,s= 80, color='blue',marker = 'o',  label = 'BP band') 
#plt.scatter(phase,BP,s= 80, c= phase, cmap = 'winter',  label = 'BP band') 

plt.text(x = 0, y = 11.9, s = "RRc 1985008784705347456 P=0.380 (d)", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
plt.legend(fontsize = 15, loc = "best")
plt.xlabel("Phase")
plt.ylabel('Magnitude (mag)')
plt.title('Epoch photometry')
plt.show()







# source_id 863983115982280448
# this is a RRab star pulsating with fundamental period pf = 0.57183



plt.rcParams['font.size'] = 20

hdulist = fits.open('EPOCH_PHOTOMETRY-Gaia DR3 863983115982280448.fits')

hdulist.info()
dat = hdulist[1].data

pf = 0.57183
b=dat.field('band')
mag = dat.field("mag")
time=dat.field("time")
G=[]
BP=[]
RP=[]
t_g =[]
t_bp=[]
t_rp=[]

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
        
 

# folded


plt.figure(figsize=(15, 10)) 

t_rp=np.array(t_rp)        # time of first transit in RP band 
phase =  ((t_rp - 1706.11532 - pf/2) / pf/2) % 1 - 0.5
plt.scatter(phase,RP,s= 80,color='red',marker = 'o',  label = 'RP band') 
#plt.scatter(phase,RP,s= 80, c= phase, cmap = 'autumn',  label = 'RP band') 

t_g=np.array(t_g)        # time of first transit in G band 
phase =  ((t_g - 1706.11492 - pf/2) / pf/2) % 1 - 0.5
plt.scatter(phase,G,s= 80,color='green',marker = 'o',  label = 'G band') 
#plt.scatter(phase,G,s= 80, c= phase, cmap = 'summer',  label = 'G band') 

t_bp=np.array(t_bp)        # time of first transit in BP band 
phase =  ((t_bp - 1706.11523 - pf/2) / pf/2) % 1 - 0.5
plt.scatter(phase,BP,s= 80, color='blue',marker = 'o',  label = 'BP band') 
#plt.scatter(phase,BP,s= 80, c= phase, cmap = 'winter',  label = 'BP band') 

plt.text(x = 0, y = 11.7, s = "RRab 863983115982280448 P=0.572 (d)", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
plt.legend(fontsize = 15, loc = "best")
plt.xlabel("Phase")
plt.ylim(13.5, 11.5)
plt.ylabel('Magnitude (mag)')
plt.title('Epoch photometry')
plt.show()





############################

# ### EPOCH RADIAL VELOCITY FOR RR LYRAE STARS



hdu_list = fits.open('rr1_rv.fits')
hdu_list.info()
dat = hdu_list[1].data

s=dat.field("source_id")
rad = dat.field("radial_velocity")
rad_err=dat.field("radial_velocity_error")
time=dat.field("rv_obs_time")

print(s)
print(rad)
print(rad_err)
print(time)



### RRc 1956531880222667904  folded


RV=[]
RV_err=[]
t=[]


i=0
for i in range(len(rad)):
    if (s[i] == 1956531880222667904):    ### RRc
        RV.append(rad[i])
        RV_err.append(rad_err[i])
        t.append(time[i])
        
RV=np.array(RV)
RV_err=np.array(RV_err)
t=np.array(t)
N=len(RV)
print(N)
print(RV)
print(RV_err)
print(t)


plt.figure(figsize=(15, 10))  

p1=0.2537 # 1956531880222667904
phase =  ((t - 2689.29186 - p1/2) / p1/2) % 1 - 0.5 # 1956531880222667904  

#convert time to a color tuple using the colormap used for scatter

norm = matplotlib.colors.Normalize(vmin=min(RV), vmax=max(RV), clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap='inferno')
phase_color = np.array([(mapper.to_rgba(v)) for v in RV])

plt.figure(figsize=(15, 10))    
plt.text(x = 0, y = 20, s = "RRc 1956531880222667904 P=0.254 (d)", horizontalalignment='center', color = 'gray')
#loop over each data point to plot
for x, y, e, color in zip(phase, RV, RV_err, phase_color):  
   plt.scatter(x,y, c=color, s = 150)
   plt.errorbar(x, y, e, lw=1, capsize=4, color=color)
   plt.legend(fontsize = 15, loc = "best")
   plt.xlabel("Phase")
   plt.ylabel('Radial velocity (km/s)')
   plt.ylim(-30,25) 
   plt.title('Epoch RV')
#plt.show() 

# this source shows the sinusoidal shape






hdu_list = fits.open('rr2_rv.fits')
hdu_list.info()
dat = hdu_list[1].data

s=dat.field("source_id")
rad = dat.field("radial_velocity")
rad_err=dat.field("radial_velocity_error")
time=dat.field("rv_obs_time")

print(s)
print(rad)
print(rad_err)
print(time)

### RRab 2926381228470699392 folded

RV=[]
RV_err=[]
t=[]

i=0
for i in range(len(rad)):
    if (s[i] == 2926381228470699392): #  2926381228470699392   ### Rab
        RV.append(rad[i])
        RV_err.append(rad_err[i])
        t.append(time[i])
        
RV=np.array(RV)
t=np.array(t)
N=len(RV)
print(N)
print(RV)
print(RV_err)
print(t)


plt.figure(figsize=(15, 10))  

pf=0.53739 # 2926381228470699392
phase =  ((t - 1803.98457 - pf/2) / pf/2) % 1 - 0.5 # 2926381228470699392

#convert time to a color tuple using the colormap used for scatter

norm = matplotlib.colors.Normalize(vmin=min(RV), vmax=max(RV), clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap='inferno')
phase_color = np.array([(mapper.to_rgba(v)) for v in RV])

plt.figure(figsize=(15, 10))    
plt.text(x = 0, y = 70, s = "RRab 2926381228470699392 P=0.537 (d)", horizontalalignment='center', color = 'gray')
#loop over each data point to plot
for x, y, e, color in zip(phase, RV, RV_err, phase_color):  
   plt.scatter(x,y, c=color, s = 150)
   plt.errorbar(x, y, e, lw=1, capsize=4, color=color)
   plt.legend(fontsize = 15, loc = "best")
   plt.xlabel("Phase")
   plt.ylabel('Radial velocity (km/s)')
   plt.ylim(-30,80) 
   plt.title('Epoch RV')
#plt.show() 


# this source shows the saw-tooth shape

























