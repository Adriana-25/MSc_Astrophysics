#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 22:33:45 2023

@author: adri
"""

### Main Sequence Oscillators

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle

plt.rcParams.update({'font.size': 20})





### All-sky Aitoff projection in Galactic coordinates


hdulist7 = fits.open('ms_oscillator1.fits')
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


hdulist8 = fits.open('ms_oscillator2.fits')
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


fig = plt.figure(figsize=(15,10))  #
ax = fig.add_subplot(1,1,1, projection='aitoff')
ax.scatter(l7, b7, s=50, color='red', marker = '*', alpha=1, label = 'MS oscillators in LOPN1')   
ax.scatter(l8, b8, s=50, color='blue',marker = '*', alpha=1,label = 'MS oscillators in LOPS2')   
ax.legend(fontsize='15')
plt.title("Aitoff projection in Galactic coordinates\n")
plt.tight_layout()
plt.grid()




###############

# Absolute Color-Magnitude Diagram


hdulist3 = fits.open('ms_oscillator1.fits')
hdulist3.info()
dat3 = hdulist3[1].data

g3=dat3.field('median_mag_g_fov')
bp3=dat3.field('median_mag_bp')
rp3=dat3.field('median_mag_rp')
p3=dat3.field("parallax")
ag3 = dat3.field("ag_gspphot") 

print(g3)
print(bp3)
print(rp3)
print(p3)
print(ag3)

N = len(g3)

color3=np.zeros(N,float)
G3 = np.zeros(N,float)
i=0
for i in range(N):
        color3[i]= bp3[i]-rp3[i]
        G3[i] = g3[i] +5 *( 1+np.log10 (p3[i]/1000) ) - ag3[i]


print(G3)



plt.figure(figsize=(15, 10))
plt.scatter(color3,G3,color='red',marker = 'o',  label = 'MS oscillators in LOPN1')


hdulist3 = fits.open('ms_oscillator2.fits')
hdulist3.info()
dat3 = hdulist3[1].data

g3=dat3.field('median_mag_g_fov')
bp3=dat3.field('median_mag_bp')
rp3=dat3.field('median_mag_rp')
p3=dat3.field("parallax")
ag3=dat3.field("ag_gspphot")

print(g3)
print(bp3)
print(rp3)
print(p3)
print(ag3)

N = len(g3)

color3=np.zeros(N,float)
G3 = np.zeros(N,float)
i=0
for i in range(N):
        color3[i]= bp3[i]-rp3[i]
        G3[i] = g3[i] +5 *( 1+np.log10 (p3[i]/1000) ) - ag3[i]


print(G3)



plt.scatter(color3,G3,color='blue',marker = 'o',  label = 'MS oscillators in LOPS2')

plt.legend(fontsize="15")
plt.gca().invert_yaxis() 
plt.xlabel('BP-RP')
plt.ylabel('G')
plt.title('aCMD')
#p.axis([0,1.4,26,18])
plt.show()



#################

# Amplitude-frequency distribution

hdulist1 = fits.open('ms_oscillator1.fits')
hdulist1.info()
dat1 = hdulist1[1].data

f1=dat1.field('frequency1')
a1=dat1.field("amplitude_g_freq1")


plt.figure(figsize=(15, 10)) 
#plt.scatter(f1,a1,s= 10, c=f1, cmap= "viridis", marker = 'o') 
plt.scatter(f1,a1,s= 100, c='red', marker = '*',label ="LOPN1") 
#plt.xlabel("Main frequency f1 $(d^{-1})$")
#plt.ylabel("Amplitude G (mag)")
#plt.tight_layout()
#plt.show()



hdulist1 = fits.open('ms_oscillator2.fits')
hdulist1.info()
dat1 = hdulist1[1].data

f1=dat1.field('frequency1')
a1=dat1.field("amplitude_g_freq1")


#plt.figure(figsize=(15, 10)) 
#plt.scatter(f1,a1,s= 10, c=f1, cmap= "viridis", marker = 'o') 
plt.scatter(f1,a1,s= 100, c='blue', marker = '*', label ="LOPS2") 

plt.xlabel("Main frequency f1 $(d^{-1})$")
plt.ylabel("Amplitude G (mag)")
plt.title("Amplitude-frequency diagram")
plt.legend(loc = "best")
plt.tight_layout()
plt.show()


#################################

# Primary frequency histogram

hdulist1 = fits.open('ms_oscillator1.fits')
hdulist1.info()
dat1 = hdulist1[1].data

f1=dat1.field('frequency1')
a1=dat1.field("amplitude_g_freq1")


plt.figure(figsize=(15, 10)) 
plt.hist(f1, bins=153, color = 'red',label ="LOPN1")
#plt.xlabel("Main frequency f1 $(d^{-1})$")
#plt.ylabel("Counts")
#plt.tight_layout()
#plt.show()



hdulist1 = fits.open('ms_oscillator2.fits')
hdulist1.info()
dat1 = hdulist1[1].data

f1=dat1.field('frequency1')
a1=dat1.field("amplitude_g_freq1")


#plt.figure(figsize=(15, 10)) 
plt.hist(f1, bins=438, color = "blue", label ="LOPS2")
plt.xlabel("Main frequency f1 $(d^{-1})$")
plt.ylabel("Counts")
plt.title("Primary frequency histogram")
plt.tight_layout()
plt.legend(loc = "best")
plt.show()

#########################

# Photometric time-seris

hdulist = fits.open('EPOCH_PHOTOMETRY-Gaia DR3 393326065521087616.fits')

hdulist.info()
dat = hdulist[1].data

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


# unfolded

plt.figure(figsize=(15, 10)) 
plt.scatter(t_rp,RP,color='red',marker = 'o',  label = 'RP band')
plt.scatter(t_g,G,color='green',marker = 'o',  label = 'G band')
plt.scatter(t_bp,BP,color='blue',marker = 'o',  label = 'BP band')

plt.text(x = 2200, y = 9.9, s = "MS oscillator Gaia DR3 393326065521087616", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
plt.legend(loc = "best")
plt.xlabel("Time (BJD)")
plt.ylabel('Magnitude (mag)')
plt.title('Epoch photometry')
plt.show()




# folded

phase1 = 2.6748378   # phase_g_freq1 for this source
f1=1.13638        # main frequency for this source
p1=1/f1        #  inverse of main frequency
plt.figure(figsize=(15, 10)) 

t_rp=np.array(t_rp)        # time of first transit in RP band 
phase =  ((t_rp - 1696.90631 - p1/2) / p1/2) % 1 - 0.5
plt.scatter(phase,RP,s= 80,color='red',marker = 'o',  label = 'RP band') 

t_g=np.array(t_g)        # time of first transit in G band
phase =  ((t_g - 1696.90593 - p1/2) / p1/2) % 1 - 0.5
plt.scatter(phase,G,s= 80,color='green',marker = 'o',  label = 'G band') 

t_bp=np.array(t_bp)        # time of first transit in BP band
phase =  ((t_bp - 1696.90622 - p1/2) / p1/2) % 1 - 0.5
plt.scatter(phase,BP,s= 80,color='blue',marker = 'o',  label = 'BP band') 

plt.text(x = 0, y = 9.9, s = "MS oscillator Gaia DR3 393326065521087616 f$_{1}$ = 1.136 $(d^{-1})$", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
#plt.ylim(10.16,10)
plt.legend(fontsize = 15, loc = "center right")
plt.xlabel("Phase")
plt.ylabel('Magnitude (mag)')
plt.title('Epoch photometry')
plt.show()



#########################à

# Period–Wesenheit relation: W index over logP
# W = M G − 1.9 (G BP − G RP )    mg from tables, G BP and G RP from median values in tables
# P = 1/f1

hdulist1 = fits.open('ms_oscillator1.fits')
hdulist1.info()
dat1 = hdulist1[1].data

f1=dat1.field('frequency1') # 153
p1=1/f1 # 153
P=np.log(p1)
mg=dat1.field("mg_gspphot")
g_bp=dat1.field("median_mag_bp") # 153
g_rp=dat1.field("median_mag_rp") # 153
W=mg -1.9*(g_bp - g_rp)

print(len(p1))
print(len(W))

plt.figure(figsize=(15, 10)) 
plt.scatter(P,W,s= 100,color='red',marker = '*',  label = 'LOPN1') 



hdulist1 = fits.open('ms_oscillator2.fits')
hdulist1.info()
dat1 = hdulist1[1].data

f1=dat1.field('frequency1') # 438
p1=1/f1 # 438
P=np.log(p1)
mg=dat1.field("mg_gspphot")
g_bp=dat1.field("median_mag_bp") # 438
g_rp=dat1.field("median_mag_rp") # 438
W=mg -1.9*(g_bp - g_rp)

print(len(p1))
print(len(W))


plt.scatter(P,W,s= 100,color='blue',marker = '*',  label = 'LOPS2') 


plt.gca().invert_yaxis()
plt.legend(loc = "best")
plt.xlabel("logP (d)")
plt.ylabel('Wesenheit G (mag)')
plt.title('Period–Wesenheit relation')
plt.show()

# what information from this plot? 
# It would probably be more significative if the sample were interely composed of delta scuti stars





#################

# Amplitude- rotation relation for delta scuti stars
# vsini_esphs from astrophysical parameters table 
# not a good plot: too few data points

hdulist1 = fits.open('ms_oscillator1.fits')
hdulist1.info()
dat1 = hdulist1[1].data

v1=dat1.field('vsini_esphs')
a1=dat1.field("amplitude_g_freq1")
a1=1000*a1

plt.figure(figsize=(15, 10)) 
plt.scatter(v1,a1,s= 300, c='red', marker = '*',label ="LOPN1") 
#plt.xlabel("Main frequency f1 $(d^{-1})$")
#plt.ylabel("Amplitude G (mag)")
#plt.tight_layout()
#plt.show()


hdulist1 = fits.open('ms_oscillator2.fits')
hdulist1.info()
dat1 = hdulist1[1].data

v1=dat1.field('vsini_esphs')
a1=dat1.field("amplitude_g_freq1")
a1=1000*a1

#plt.figure(figsize=(15, 10)) 
plt.scatter(v1,a1,s= 300, c='blue', marker = '*',label ="LOPS2") 
plt.xlabel("vsini $(km s^{-1})$")
plt.ylabel("Amplitude G (mmag)")
plt.title("Amplitude-velocity diagram")
plt.legend(loc = "best")
plt.tight_layout()
plt.show()


# From ADQL query interface in Gaia Archive, all MS pulsators result DSCT|GDOR|SXPHE
