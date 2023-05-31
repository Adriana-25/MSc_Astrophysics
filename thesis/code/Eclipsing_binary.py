#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 22:51:36 2023

@author: adri
"""

######## eclipsing binaries ##########

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle

plt.rcParams.update({'font.size': 20})


### All-sky Aitoff projection in Galactic coordinates


hdulist7 = fits.open('eclipsing_binaries1.fits')
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
ax.scatter(l7, b7, s=25, color='red', marker = '.', alpha=1, label = 'Eclipsing binaries in LOPN1')   


hdulist8 = fits.open('eclipsing_binaries2.fits')
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



ax.scatter(l8, b8, s=25, color='blue',marker = '.', alpha=1,label = 'Eclipsing binaries in LOPS2')   
ax.legend(fontsize='15')
plt.title("Aitoff projection in Galactic coordinates\n")
plt.tight_layout()
plt.grid()




########### Absolute Color-Magnitude Diagram ##########




hdulist5 = fits.open('eclipsing_binaries1.fits')
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
plt.scatter(color5,G5,color='red',marker = '.',  label = 'eclipsing binaries in LOPN1')


hdulist5 = fits.open('eclipsing_binaries2.fits')
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



plt.scatter(color5,G5,color='blue',marker = '.',  label = 'eclipsing binaries in LOPS2')

plt.legend(fontsize="15")
plt.gca().invert_yaxis() 
plt.xlabel('BP-RP')
plt.ylabel('G')
plt.title('aCMD')
#p.axis([0,1.4,26,18])
plt.show()




##### period histogram #######



hdulist = fits.open('eclipsing_binaries1.fits')
hdulist.info()
dat = hdulist[1].data

f1=dat.field('frequency')
p1=1/f1
print(len(p1))

fig, ax = plt.subplots(figsize=(15, 10)) 

ax.hist(p1, bins=2977,color = 'red',label ="LOPN1")

hdulist = fits.open('eclipsing_binaries2.fits')
hdulist.info()
dat = hdulist[1].data

f2=dat.field('frequency')
p2=1/f2
print(len(p2))


ax.hist(p2, bins=3466, color = "blue", label ="LOPS2")
ax.set_xlabel("Period (d)")
ax.set_xlim(0,10)
ax.set_ylim(0,750)
ax.set_ylabel("Counts")
ax.set_title("Period histogram")



axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])
axins.hist(p1,bins=92,color = 'red',label ="LOPN1" )
axins.hist(p2, bins=131, color = "blue", label ="LOPS2")
# subregion of the original image
x1, x2, y1, y2 = 10, 150, 0, 20
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xlabel("Period (d)")
axins.set_ylabel("Counts")

ax.legend(loc = "best", fontsize = "15")

plt.show()





i=0
k=0
for i in range(len(p1)):
    if (p1[i] > 50):
        k+=1
 
    
print(k)   # 92 out of 2977 are longer than 10 d, 41 are longer than 20d
           # 15 are longer than 50 d
i=0
k=0
for i in range(len(p2)):
    if (p2[i] > 50):
        k+=1

print(k)  # 131 out of 3466 are longer than 10 d, 64 are longer than 20 d
          # 19 are longer than 50 d
          
          
      
          
          

##### epoch photometry #######

hdulist = fits.open('EPOCH_PHOTOMETRY-Gaia DR3 393454948901343360.fits')
hdulist.info()
dat = hdulist[1].data

b=dat.field('band')
mag = dat.field("mag")
time=dat.field("time")

f1 = 0.5829491038768037  # eclipse frequency
p1 = 1/f1 # eclipse period
print(p1)
t_ref = 2100.833606832954  # reference time from vari table

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
#for ii in range(0,10):
plt.axvline(t_ref + p1, c='C3', alpha=0.25)
plt.text(x = 2200, y = 11.7, s = "Eclipsing binary Gaia DR3 393454948901343360", horizontalalignment='center', color = 'gray')
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

#plt.text(x = 0, y = 11.65, s = "Eclipsing binary Gaia DR3 393454948901343360", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
#plt.xlim(-0.5,0.1)
#plt.ylim(12.08,12.02) # BP band
#plt.ylim(11.85,11.75) # G band
#plt.ylim(11.5,11.2) # RP band
#plt.axhline(11.37, c ='C2', alpha=1)
plt.legend(fontsize = 15, loc = "lower right")
plt.xlabel("Phase")
plt.ylabel('Magnitude (mag)')
plt.title('Epoch photometry')
plt.show()



########################

# XP SPECTRA

# Gaia DR3 393454948901343360

hdulist = fits.open('XP_SAMPLED-Gaia DR3 393454948901343360.fits')

hdulist.info()
dat = hdulist[1].data


wave=dat.field('wavelength')
flux = dat.field("flux")
flux_err=dat.field("flux_error")

plt.figure(figsize=(15, 10)) 
plt.scatter(wave,flux, c='blue', s = 80)
plt.errorbar(wave, flux, flux_err, lw=1, capsize=4)
plt.text(x = 800, y = 7*1e-16, s = "Eclipsing binary Gaia DR3 393454948901343360", horizontalalignment='center', color = 'gray')
plt.legend(fontsize = 15, loc = "best")
plt.xlabel("$\lambda$ (nm)")
plt.ylabel(r'$ \mathrm{ \nu F(\nu) (erg \cdot m^{-2} s^{-1} nm^{-1}) }$')
#plt.ylim(-30,25) 
plt.title('XP mean spectrum')
plt.show()


############################À

# RVS SPECTRA


# Gaia DR3 2887834618539920256

hdulist = fits.open('RVS-Gaia DR3 2887834618539920256.fits')

hdulist.info()
dat = hdulist[1].data


wave=dat.field('wavelength')
flux = dat.field("flux")
flux_err=dat.field("flux_error")

plt.figure(figsize=(15, 10)) 
plt.scatter(wave,flux, c='blue', s = 40)
plt.errorbar(wave, flux, flux_err, lw=0.5, capsize=3)
plt.text(x = 858, y = 1.15, s = "Eclipsing binary Gaia DR3 2887834618539920256", horizontalalignment='center', color = 'gray')
plt.legend(fontsize = 15, loc = "best")
plt.xlabel("$\lambda$ (nm)")
plt.ylabel("normalized flux")
#plt.ylabel(r'$ \mathrm{ \nu F(\nu) (erg \cdot m^{-2} s^{-1} nm^{-1}) }$')
plt.ylim(0.6,1.2) 
plt.title('RVS mean spectrum')
plt.show()




