#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 19:15:13 2023

@author: adri
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle

plt.rcParams.update({'font.size': 20})

### All-sky Aitoff projection in Galactic coordinates

hdulist = fits.open('planetary_transits1.fits')
hdulist.info()
dat = hdulist[1].data
l=dat.field('l')
b=dat.field('b')
l = Angle(l * u.deg)
l.wrap_at('180d', inplace=True)
print(l)
b = Angle(b * u.deg)
b.wrap_at('90d', inplace=True)
print(b)
l = np.deg2rad(np.array(l))
b = np.deg2rad(np.array(b))
print(np.deg2rad(np.array(l)))
print(np.deg2rad(np.array(b)))


fig = plt.figure(figsize=(15,12))
ax = fig.add_subplot(1,1,1, projection='aitoff')
ax.scatter(l, b, s=50, marker='*', color='red', alpha=1, label = 'planetary transits in LOPN1')  


hdulist2 = fits.open('planetary_transits2.fits')
hdulist2.info()
dat2 = hdulist2[1].data
l2=dat2.field('l')
b2=dat2.field('b')
l2 = Angle(l2 * u.deg)
l2.wrap_at('180d', inplace=True)
print(l2)
b2 = Angle(b2 * u.deg)
b2.wrap_at('90d', inplace=True)
print(b2)
l2 = np.deg2rad(np.array(l2))
b2 = np.deg2rad(np.array(b2))
 

ax.scatter(l2, b2, s=50, marker='*', color='blue', alpha=1,  label = 'planetary transits in LOPN1')
ax.grid()
plt.legend(fontsize="15", loc ="lower right")
plt.title('All-sky map in Galactic coordinates\n')


#######################à



# Absolute Color-Magnitude Diagram


hdulist1 = fits.open('planetary_transits1.fits')
hdulist1.info()
dat1 = hdulist1[1].data

g1=dat1.field('median_mag_g_fov')
bp1=dat1.field('median_mag_bp')
rp1=dat1.field('median_mag_rp')
p1=dat1.field("parallax")
ag1 = dat1.field("ag_gspphot")  

print(g1)
print(bp1)
print(rp1)
print(p1)
print(ag1)

N = len(g1)

color1=np.zeros(N,float)
G1 = np.zeros(N,float)
i=0
for i in range(N):
        color1[i]=bp1[i]-rp1[i]
        G1[i] = g1[i]+5 *( 1+np.log10 (p1[i]/1000) ) - ag1[i]


print(G1)


plt.figure(figsize=(15, 10))
plt.scatter(color1,G1,color='red',marker = 'o',  label = 'planetary transits in LOPN1')


hdulist1 = fits.open('planetary_transits2.fits')
hdulist1.info()
dat1 = hdulist1[1].data

g1=dat1.field('median_mag_g_fov')
bp1=dat1.field('median_mag_bp')
rp1=dat1.field('median_mag_rp')
p1=dat1.field("parallax")
ag1=dat1.field("ag_gspphot")

print(g1)
print(bp1)
print(rp1)
print(p1)
print(ag1)

N = len(g1)

color1=np.zeros(N,float)
G1 = np.zeros(N,float)
i=0
for i in range(N):
        color1[i]=bp1[i]-rp1[i]
        G1[i] = g1[i]+5 *( 1+np.log10 (p1[i]/1000) ) - ag1[i]


print(G1)


plt.scatter(color1,G1,color='blue',marker = 'o',  label = 'planetary transits in LOPS2')
plt.legend(fontsize="15")
plt.gca().invert_yaxis() 
plt.xlabel('BP-RP')
plt.ylabel('G')
plt.title('aCMD')
#p.axis([0,1.4,26,18])
plt.show()




####### non periodic variations #########



hdulist = fits.open('planetary_transits1.fits')
hdulist.info()
dat = hdulist[1].data

skew=dat.field('skewness_mag_g_fov')
abbe = dat.field("abbe_mag_g_fov")

plt.figure(figsize=(15, 10)) 
plt.scatter(abbe,skew,s= 100, color='blue',marker = '*') 


hdulist = fits.open('planetary_transits2.fits')
hdulist.info()
dat = hdulist[1].data

skew=dat.field('skewness_mag_g_fov')
abbe = dat.field("abbe_mag_g_fov")

#plt.figure(figsize=(15, 10)) 
plt.scatter(abbe,skew,s= 100, color='blue',marker = '*') 
#plt.scatter(abbe,skew,s= 80, c= abbe, cmap = 'plasma') 
#plt.legend(fontsize = 15, loc = "best")
plt.xlabel("abbe_mag_g_fov")
#plt.xlim(0.3,1.5)
#plt.gca().invert_yaxis()
plt.ylabel('skewness_mag_g_fov')
#plt.ylim(0.75,-2) 
plt.title('Metrics targeting non-periodic variations')
plt.show()





####### periodic variations #######


hdulist = fits.open('planetary_transits1.fits')
hdulist.info()
dat = hdulist[1].data

g=dat.field('median_mag_g_fov')
s_bp=dat.field("std_dev_mag_bp")
s_rp=dat.field("std_dev_mag_rp")
y=np.log10(s_bp/s_rp)  ### ratio between the standard deviations in magnitude of FoV transit observations in the G BP and G RP bands

plt.figure(figsize=(15, 10)) 
#plt.scatter(g,y,s= 80, c= g, cmap = 'plasma_r') 
plt.scatter(g,y,s= 100,color='purple',marker = '*' ) 


hdulist = fits.open('planetary_transits2.fits')
hdulist.info()
dat = hdulist[1].data

g=dat.field('median_mag_g_fov')
s_bp=dat.field("std_dev_mag_bp")
s_rp=dat.field("std_dev_mag_rp")
y=np.log10(s_bp/s_rp)  ### ratio between the standard deviations in magnitude of FoV transit observations in the G BP and G RP bands

#plt.figure(figsize=(15, 10)) 
plt.scatter(g,y,s= 100, color='purple',marker = '*') 
#plt.scatter(g,y,s= 80, c= g, cmap = 'plasma_r') 
#plt.legend(fontsize = 15, loc = "best")
plt.xlabel("median_mag_g_fov")
#plt.gca().invert_xaxis()
#plt.xlim(8.5,14)
plt.ylabel('std_dev_mag_bp/std_dev_mag_rp')
#plt.ylim(0.05,0.35) 
plt.title('Metrics targeting periodic variations')
plt.show()





############  time series standard deviation as a function of magnitude  ############



hdulist = fits.open('planetary_transits1.fits')
hdulist.info()
dat = hdulist[1].data

g=dat.field('median_mag_g_fov')
s_g=dat.field("std_dev_mag_g_fov")


plt.figure(figsize=(15, 10)) 
#plt.scatter(g,y,s= 80, c= g, cmap = 'plasma_r') 
plt.scatter(g,s_g ,s= 100,color='green',marker = '*' ) 


hdulist = fits.open('planetary_transits2.fits')
hdulist.info()
dat = hdulist[1].data

g=dat.field('median_mag_g_fov')
s_g=dat.field("std_dev_mag_g_fov")


plt.scatter(g,s_g, s= 100, color='green',marker = '*') 
#plt.scatter(g,y,s= 80, c= g, cmap = 'plasma_r') 
#plt.legend(fontsize = 15, loc = "best")
plt.xlabel("median_mag_g_fov")
#plt.gca().invert_xaxis()
#plt.xlim(8.5,14)
plt.ylabel('std_dev_mag_g_fov')
#plt.ylim(0.05,0.35) 
plt.title('Time series variations')
plt.show()





####### Epoch photometry #######

# Gaia DR3 1864885215233116032 
# TOI 1282.01
# Star KELT 16 Yellow-White star (spectral class F7V)
# https://www.stellarcatalog.com/stars/kelt-16

hdulist = fits.open('EPOCH_PHOTOMETRY-Gaia DR3 1864885215233116032.fits')

hdulist.info()
dat = hdulist[1].data

b=dat.field('band')
mag = dat.field("mag")
time=dat.field("time")
flux=dat.field("flux")

p1=0.9689828586932296 # transit period
t_ref= 2203.8143098946257    # planetary transit reference time
# transit depth: 13.955157 n transits = 35
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





# unfolded magnitude

plt.figure(figsize=(15, 10)) 
plt.scatter(t_rp,RP,color='red',marker = 'o',  label = 'RP band')
plt.scatter(t_g,G,color='green',marker = 'o',  label = 'G band')
plt.scatter(t_bp,BP,color='blue',marker = 'o',  label = 'BP band')
#for ii in range(0,10):
#    plt.axvline(t_ref + p1*ii, c='C3', alpha=0.25)
plt.text(x = 2200, y = 11.5, s = "Exoplanet host Gaia DR3 1864885215233116032", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
plt.legend(loc = "best")
plt.xlabel("Time (BJD)")
plt.ylabel('Magnitude (mag)')
#plt.ylim(10.16,10)
#plt.xlim(2200,2400)
plt.title('Epoch photometry')
plt.show()




# folded magnitude 


plt.figure(figsize=(15, 10)) 

t_rp=np.array(t_rp)      
phase =  ((t_rp - t_ref - p1/2) / p1) % 1 - 0.5
plt.scatter(phase,RP,s= 80,color='red',marker = 'o',  label = 'RP band') 

t_g=np.array(t_g)        
phase =  ((t_g - t_ref - p1/2) / p1) % 1 - 0.5
plt.scatter(phase,G,s= 80,color='green',marker = 'o',  label = 'G band') 

t_bp=np.array(t_bp)    
phase =  ((t_bp - t_ref - p1/2) / p1) % 1 - 0.5
plt.scatter(phase,BP,s= 80,color='blue',marker = 'o',  label = 'BP band') 

plt.text(x = 0, y = 11.5, s = "Exoplanet host Gaia DR3 1864885215233116032", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
#plt.xlim(-0.4,0.1)
#plt.ylim(12.08,12.02) # BP band
#plt.ylim(11.85,11.75) # G band
#plt.ylim(11.4,11.3) # RP band
#plt.axhline(11.37, c ='C2', alpha=1)
plt.legend(fontsize = 15, loc = "center right")
plt.xlabel("Phase")
plt.ylabel('Magnitude (mag)')
plt.title('Epoch photometry')
plt.show()





# folded magnitude with zoom


plt.figure(figsize=(15, 10)) 

t_rp=np.array(t_rp)      
phase =  ((t_rp - t_ref - p1/2) / p1) % 1 - 0.5
RP=np.array(RP)
RP=RP+0.6
plt.scatter(phase,RP,s= 80,color='red',marker = 'o',  label = 'RP + 0.6 mag') 
#plt.text(x = 0.4, y = 11.96, s = "RP + 0.6 mag", horizontalalignment='center', color = 'red')

t_g=np.array(t_g)        
phase =  ((t_g - t_ref - p1/2) / p1) % 1 - 0.5
G=np.array(G)
G=G+0.22
plt.scatter(phase,G,s= 80,color='green',marker = 'o',  label = 'G + 0.2 mag') 
#plt.text(x = 0.4, y = 12, s = "G + 0.2 mag", horizontalalignment='center', color = 'green')

t_bp=np.array(t_bp)    
phase =  ((t_bp - t_ref - p1/2) / p1) % 1 - 0.5
plt.scatter(phase,BP,s= 80,color='blue',marker = 'o',  label = 'BP band') 

plt.text(x = 0, y = 11.945, s = "Exoplanet host Gaia DR3 1864885215233116032", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
#plt.xlim(-0.4,0.1)
#plt.ylim(12.08,12.02) # BP band
#plt.ylim(11.85,11.75) # G band
#plt.ylim(11.4,11.3) # RP band
#plt.axhline(11.37, c ='C2', alpha=1)
plt.legend(fontsize = 15, loc = "upper right")
plt.xlabel("Phase")
plt.ylabel('Magnitude (mag)')
plt.title('Epoch photometry')
plt.show()









# epoch photometry for the exoplanet candidate #

hdulist = fits.open('EPOCH_PHOTOMETRY-Gaia DR3 2177849036734191744.fits')

hdulist.info()
dat = hdulist[1].data

b=dat.field('band')
mag = dat.field("mag")
time=dat.field("time")
flux=dat.field("flux")

p1=3.80952 # transit period
t_ref= 2223.76861    # planetary transit reference time
# transit depth: 13.955157 n transits = 35
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





# unfolded magnitude

plt.figure(figsize=(15, 10)) 
plt.scatter(t_rp,RP,color='red',marker = 'o',  label = 'RP band')
plt.scatter(t_g,G,color='green',marker = 'o',  label = 'G band')
plt.scatter(t_bp,BP,color='blue',marker = 'o',  label = 'BP band')
#for ii in range(0,10):
#    plt.axvline(t_ref + p1*ii, c='C3', alpha=0.25)
plt.text(x = 2200, y = 10.8, s = "Candidate exoplanet host Gaia DR3 2177849036734191744", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
plt.legend(loc = "best")
plt.xlabel("Time (BJD)")
plt.ylabel('Magnitude (mag)')
#plt.ylim(10.16,10)
#plt.xlim(2200,2400)
plt.title('Epoch photometry')
plt.show()




# folded magnitude 


plt.figure(figsize=(15, 10)) 

t_rp=np.array(t_rp)      
phase =  ((t_rp - t_ref - p1/2) / p1) % 1 - 0.5
plt.scatter(phase,RP,s= 80,color='red',marker = 'o',  label = 'RP band') 

t_g=np.array(t_g)        
phase =  ((t_g - t_ref - p1/2) / p1) % 1 - 0.5
plt.scatter(phase,G,s= 80,color='green',marker = 'o',  label = 'G band') 

t_bp=np.array(t_bp)    
phase =  ((t_bp - t_ref - p1/2) / p1) % 1 - 0.5
plt.scatter(phase,BP,s= 80,color='blue',marker = 'o',  label = 'BP band') 

plt.text(x = 0, y = 10.8, s = "Candidate exoplanet host Gaia DR3 2177849036734191744", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
#plt.xlim(-0.4,0.4)
#plt.ylim(12.08,12.02) # BP band
#plt.ylim(11.85,11.75) # G band
#plt.ylim(11.4,11.3) # RP band
#plt.axhline(11.37, c ='C2', alpha=1)
plt.legend(fontsize = 15, loc = "center right")
plt.xlabel("Phase")
plt.ylabel('Magnitude (mag)')
plt.title('Epoch photometry')
plt.show()




# folded magnitude with zoom


plt.figure(figsize=(15, 10)) 

t_rp=np.array(t_rp)      
phase =  ((t_rp - t_ref - p1/2) / p1) % 1 - 0.5
RP=np.array(RP)
RP=RP+0.6
plt.scatter(phase,RP,s= 80,color='red',marker = 'o',  label = 'RP + 0.6 mag') 
#plt.text(x = 0.4, y = 11.96, s = "RP + 0.6 mag", horizontalalignment='center', color = 'red')

t_g=np.array(t_g)        
phase =  ((t_g - t_ref - p1/2) / p1) % 1 - 0.5
G=np.array(G)
G=G+0.22
plt.scatter(phase,G,s= 80,color='green',marker = 'o',  label = 'G + 0.2 band') 
#plt.text(x = 0.4, y = 12, s = "G + 0.2 mag", horizontalalignment='center', color = 'green')

t_bp=np.array(t_bp)    
phase =  ((t_bp - t_ref - p1/2) / p1) % 1 - 0.5
plt.scatter(phase,BP,s= 80,color='blue',marker = 'o',  label = 'BP band') 

plt.text(x = 0, y = 11.22, s = "Candidate exoplanet Gaia host DR3 2177849036734191744", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
plt.xlim(-0.4,0.4)
plt.ylim(11.36,11.20) # BP band
#plt.ylim(11.85,11.75) # G band
#plt.ylim(11.4,11.3) # RP band
#plt.axhline(11.37, c ='C2', alpha=1)
plt.legend(fontsize = 15, loc = "lower right")
plt.xlabel("Phase")
plt.ylabel('Magnitude (mag)')
plt.title('Epoch photometry')
plt.show()







############################

# XP SPECTRA

# Gaia DR3 1864885215233116032

hdulist = fits.open('LANETARY/XP_SAMPLED-Gaia DR3 1864885215233116032.fits')

hdulist.info()
dat = hdulist[1].data


wave=dat.field('wavelength')
flux = dat.field("flux")
flux_err=dat.field("flux_error")

plt.figure(figsize=(15, 10)) 
plt.scatter(wave,flux, c='blue', s = 80)
plt.errorbar(wave, flux, flux_err, lw=1, capsize=4)
plt.text(x = 800, y = 7.5*1e-16, s = "Exoplanet host Gaia DR3 1864885215233116032", horizontalalignment='center', color = 'gray')
plt.legend(fontsize = 15, loc = "best")
plt.xlabel("$\lambda$ (nm)")
plt.ylabel(r'$ \mathrm{ \nu F(\nu) (erg \cdot m^{-2} s^{-1} nm^{-1}) }$')
#plt.ylim(-30,25) 
plt.title('XP mean spectrum')
plt.show()




############################

# RVS SPECTRA

# Gaia DR3 2102117871259036672    Kepler-7

hdulist = fits.open('RVS-Gaia DR3 2102117871259036672.fits')

hdulist.info()
dat = hdulist[1].data


wave=dat.field('wavelength')
flux = dat.field("flux")
flux_err=dat.field("flux_error")

plt.figure(figsize=(15, 10)) 
plt.scatter(wave,flux, c='blue', s = 40)
plt.errorbar(wave, flux, flux_err, lw=0.5, capsize=3)
plt.text(x = 858, y = 1.3, s = "Exoplanet host Gaia DR3 2102117871259036672", horizontalalignment='center', color = 'gray')
plt.legend(fontsize = 15, loc = "best")
plt.xlabel("$\lambda$ (nm)")
plt.ylabel("normalized flux")
#plt.ylabel(r'$ \mathrm{ \nu F(\nu) (erg \cdot m^{-2} s^{-1} nm^{-1}) }$')
plt.ylim(0.2,1.4) 
plt.title('RVS mean spectrum')
plt.show()





