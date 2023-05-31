#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 13:07:03 2023

@author: adri
"""

###### short timescale variables #######

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle

plt.rcParams.update({'font.size': 20})


### All-sky Aitoff projection in Galactic coordinates


hdulist7 = fits.open('short_timescale1.fits')
hdulist7.info()
dat7 = hdulist7[1].data
l7=dat7.field('l')
b7=dat7.field('b')
l7 = Angle(l7 * u.deg)
b7 = Angle(b7 * u.deg)
print(l7)
print(b7)
l7.wrap_at('180d', inplace=True)
b7.wrap_at('90d', inplace=True)
l7 = np.deg2rad(np.array(l7))
b7 = np.deg2rad(np.array(b7))


fig = plt.figure(figsize=(15,10))  #
ax = fig.add_subplot(1,1,1, projection='aitoff')
ax.scatter(l7, b7, s=25, color='red', marker = '.', alpha=1, label = 'Short-timescale variables in LOPN1')   


hdulist8 = fits.open('short_timescale2.fits')
hdulist8.info()
dat8 = hdulist8[1].data
l8=dat8.field('l')
b8=dat8.field('b')
l8 = Angle(l8 * u.deg)
l8.wrap_at('180d', inplace=True)
print(l8)
b8 = Angle(b8 * u.deg)
b8.wrap_at('90d', inplace=True)
print(b8)
l8 = np.deg2rad(np.array(l8))
b8 = np.deg2rad(np.array(b8))



ax.scatter(l8, b8, s=25, color='blue',marker = '.', alpha=1,label = 'Short-timescale variables in LOPS2')   
ax.legend(fontsize='15')
plt.title("Aitoff projection in Galactic coordinates\n")
plt.tight_layout()
plt.grid()





########### Absolute Color-Magnitude Diagram ##########




hdulist5 = fits.open('short_timescale1.fits')
hdulist5.info()
dat5 = hdulist5[1].data

g5=dat5.field('median_mag_g_fov')
bp5=dat5.field('median_mag_bp')
rp5=dat5.field('median_mag_rp')
p5=dat5.field("parallax")
ag5 = dat5.field("ag_gspphot") 

print(g5)
print(bp5)
print(rp5)
print(p5)
print(ag5)

N = len(g5)

color5=np.zeros(N,float)
G5 = np.zeros(N,float)
i=0
for i in range(N):
        color5[i]=bp5[i]-rp5[i]
        G5[i] = g5[i] +5 *( 1+np.log10 (p5[i]/1000) ) - ag5[i]



plt.figure(figsize=(15, 10))
plt.scatter(color5,G5,color='red',marker = '.',  label = 'Short-timescale variables in LOPN1')


hdulist5 = fits.open('short_timescale2.fits')
hdulist5.info()
dat5 = hdulist5[1].data

g5=dat5.field('median_mag_g_fov')
bp5=dat5.field('median_mag_bp')
rp5=dat5.field('median_mag_rp')
p5=dat5.field("parallax")
ag5=dat5.field("ag_gspphot")

print(g5)
print(bp5)
print(rp5)
print(p5)
print(ag5)

N = len(g5)

color5=np.zeros(N,float)
G5 = np.zeros(N,float)
i=0
for i in range(N):
        color5[i]=bp5[i]-rp5[i]
        G5[i] = g5[i] +5 *( 1+np.log10 (p5[i]/1000) ) - ag5[i]



plt.scatter(color5,G5,color='blue',marker = '.',  label = 'Short-timescale variables in LOPS2')

plt.legend(fontsize="15")
plt.gca().invert_yaxis() 
plt.xlabel('BP-RP')
plt.ylabel('G')
plt.title('aCMD')
#p.axis([0,1.4,26,18])
plt.show()





########### periods distribution histo  ##########


    
hdulist = fits.open('short_timescale2.fits')
hdulist.info()
dat = hdulist[1].data

f2=dat.field('frequency')
p2=1/f2
print(len(p2))

fig, ax = plt.subplots(figsize=(15, 10)) 

ax.hist(p2, bins=50, color = "blue", label ="LOPS2")

hdulist = fits.open('short_timescale1.fits')
hdulist.info()
dat = hdulist[1].data

f1=dat.field('frequency')
p1=1/f1
print(len(p1))


ax.hist(p1, bins=50,color = 'red',label ="LOPN1")
ax.set_xlabel("Period (d)")
#ax.set_xlim(0,10)
#ax.set_ylim(0,750)
ax.set_ylabel("Counts")
ax.set_title("Period histogram")
ax.legend(loc = "best", fontsize = "15")


axins1 = ax.inset_axes([0.5, 0.1, 0.35, 0.35])
axins1.hist(p2, bins=50, color = "blue", label ="LOPS2" )
axins1.hist(p1, bins=50,color = 'red',label ="LOPN1")
# subregion of the original image
x1, x2, y1, y2 = 0.1, 1, 0, 250
axins1.set_xlim(x1, x2)
axins1.set_ylim(y1, y2)
axins1.set_xlabel("Period (d)")
axins1.set_ylabel("Counts")


axins2 = ax.inset_axes([0.5, 0.6, 0.35, 0.35])  # 5.5, 6.6, 2.4, 2.4
axins2.hist(p2, bins=500, color = "blue", label ="LOPS2" )
axins2.hist(p1, bins=500,color = 'red',label ="LOPN1")
# subregion of the original image
x1, x2, y1, y2 = 0, 0.1, 0, 450
axins2.set_xlim(x1, x2)
axins2.set_ylim(y1, y2)
axins2.set_xlabel("Period (d)")
axins2.set_ylabel("Counts")

plt.show()
   
    

i=0
k=0
for i in range(len(p1)):
    if (p1[i] < 0.1):
        k+=1
    
print(k)   # 1448 out of 2996 are longer than 0.1 d, 1548 are shorter than 0.1 d.


i=0
k=0
for i in range(len(p2)):
    if (p2[i] < 0.1):
        k+=1

print(k)  # 1818 out of 4081 are longer than 0.1 d, 2263 are shorter than 0.1 d





######### freq-amplitude #######


hdulist = fits.open('short_timescale1.fits')
hdulist.info()
dat = hdulist[1].data

f1=dat.field('frequency')
a1=dat.field('amplitude_estimate')


fig=plt.figure(figsize=(15, 10)) 

plt.scatter(f1, a1,s=50, marker = '*', color = "red", label ="LOPN1")

hdulist = fits.open('short_timescale2.fits')
hdulist.info()
dat = hdulist[1].data

f2=dat.field('frequency')
a2=dat.field("amplitude_estimate")

plt.scatter(f2, a2 , s = 50, marker = '*', color = "blue", label ="LOPS2")
plt.legend(loc="best")
plt.xlabel("Frequency $(d^{-1})$")
plt.ylabel('Amplitude (mag)')
plt.title('Frequency-amplitude diagram')
plt.show()







######### epoch photometry ##########



hdulist = fits.open('EPOCH_PHOTOMETRY-Gaia DR3 4594961454434214144.fits')
hdulist.info()
dat = hdulist[1].data

b=dat.field('band')
mag = dat.field("mag")
time=dat.field("time")

f1 = 5.54246  #  frequency
p1 = 1/f1 #  period
print(p1)
#t_ref = 2140  # reference time 2140.13085 (sudden decrease in mag)

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
#for ii in range(0,2):
#    plt.axvline(t_ref + p1*ii, c='C2', alpha=0.5)  # what does it mean?
#plt.axvline(t_ref - p1, c ='C2', alpha =0.5)
plt.text(x = 2200, y = 10.2, s = "Short-timescale variable Gaia DR3  4594961454434214144", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
plt.legend(fontsize = 15, loc = "lower right")
plt.xlabel("Time (BJD)")
plt.ylabel('Magnitude (mag)')
plt.ylim(11.8,10)
#plt.xlim(2200,2400)
plt.title('Epoch photometry')
plt.show()




# folded magnitude 

plt.figure(figsize=(15, 10)) 

t_rp=np.array(t_rp)       # first transit in RP band  
phase =  ((t_rp - 1722.25021 - p1/2) / p1/2) % 1 -0.5
plt.scatter(phase,RP,s= 80,color='red',marker = 'o',  label = 'RP band') 

t_g=np.array(t_g)        # first transit in G band   
phase =  ((t_g - 1722.24983 - p1/2) / p1/2) % 1 - 0.5
plt.scatter(phase,G,s= 80,color='green',marker = 'o',  label = 'G band') 

t_bp=np.array(t_bp)    # first transit in BP band   
phase =  ((t_bp - 1722.25012 - p1/2) / p1/2) % 1 - 0.5
plt.scatter(phase,BP,s= 80,color='blue',marker = 'o',  label = 'BP band') 

plt.text(x = 0, y = 10.2, s = "Short-timescale variable Gaia DR3  4594961454434214144", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
#plt.xlim(-0.2,0.2)
plt.ylim(11.8,10) 
#plt.ylim(11.85,11.75) # G band
#plt.ylim(11.5,11.2) # RP band
#plt.axhline(11.37, c ='C2', alpha=1)
plt.legend(fontsize = 15, loc = "lower right")
plt.xlabel("Phase")
plt.ylabel('Magnitude (mag)')
plt.title('Epoch photometry')
plt.show()





########################################à

# xp spectra

# Gaia DR3 4594961454434214144


hdulist = fits.open('XP_SAMPLED-Gaia DR3 4594961454434214144.fits')

hdulist.info()
dat = hdulist[1].data


wave=dat.field('wavelength')
flux = dat.field("flux")
flux_err=dat.field("flux_error")

plt.figure(figsize=(15, 10)) 
plt.scatter(wave,flux, c='blue', s = 80)
plt.errorbar(wave, flux, flux_err, lw=1, capsize=4)
plt.text(x = 800, y = 1.5*1e-15, s = "Short timescale Gaia DR3 4594961454434214144", horizontalalignment='center', color = 'gray')
plt.legend(fontsize = 15, loc = "best")
plt.xlabel("$\lambda$ (nm)")
plt.ylabel(r'$ \mathrm{ \nu F(\nu) (erg \cdot m^{-2} s^{-1} nm^{-1}) }$')
#plt.ylim(-30,25) 
plt.title('XP mean spectrum')
plt.show()



########################################

# rvs spectra


# Gaia DR3 6344493361434378624

hdulist = fits.open('RVS-Gaia DR3 6344493361434378624.fits')

hdulist.info()
dat = hdulist[1].data


wave=dat.field('wavelength')
flux = dat.field("flux")
flux_err=dat.field("flux_error")

plt.figure(figsize=(15, 10)) 
plt.scatter(wave,flux, c='blue', s = 40)
plt.errorbar(wave, flux, flux_err, lw=0.5, capsize=3)
plt.text(x = 858, y = 1.15, s = "Short timescale Gaia DR3 6344493361434378624", horizontalalignment='center', color = 'gray')
plt.legend(fontsize = 15, loc = "best")
plt.xlabel("$\lambda$ (nm)")
plt.ylabel("normalized flux")
#plt.ylabel(r'$ \mathrm{ \nu F(\nu) (erg \cdot m^{-2} s^{-1} nm^{-1}) }$')
plt.ylim(0.75,1.2) 
plt.title('RVS mean spectrum')
plt.show()


