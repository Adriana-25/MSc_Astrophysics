#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 23:40:06 2023

@author: adri
"""


from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np


### absolute Color-Magnitude Diagram for LOPN1 variable objects

 
#planetary transits
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


#MS oscillator
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



#rotation modulation
hdulist4 = fits.open('rotation_modulation1.fits')
hdulist4.info()
dat4 = hdulist4[1].data

g4=dat4.field('median_mag_g_fov')
bp4=dat4.field('median_mag_bp')
rp4=dat4.field('median_mag_rp')
p4=dat4.field("parallax")
ag4 = dat4.field("ag_gspphot") 

print(g4)
print(bp4)
print(rp4)
print(p4)
print(ag4)

N = len(g4)

color4=np.zeros(N,float)
G4 = np.zeros(N,float)
i=0
for i in range(N):
        color4[i]= bp4[i]-rp4[i]
        G4[i] = g4[i] +5 *( 1+np.log10 (p4[i]/1000) ) - ag4[i]
  
print(G4)



#eclipsing binaries
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



#short timescale
hdulist6 = fits.open('short_timescale1.fits')
hdulist6.info()
dat6 = hdulist6[1].data

g6=dat6.field('median_mag_g_fov')
bp6=dat6.field('median_mag_bp')
rp6=dat6.field('median_mag_rp')
p6=dat6.field("parallax")
ag6 = dat6.field("ag_gspphot") 

print(g6)
print(bp6)
print(rp6)
print(p6)
print(ag6)

N6 = len(g6)

color6=np.zeros(N6,float)
G6=np.zeros(N6,float)

i=0
for i in range(N6):
        color6[i]=bp6[i]-rp6[i]
        G6[i] = g6[i] +5 *( 1+np.log10 (p6[i]/1000) ) - ag6[i]


#cepheids
hdulist7 = fits.open('cepheid1.fits')
hdulist7.info()
dat7 = hdulist7[1].data

g7=dat7.field('median_mag_g_fov')
bp7=dat7.field('median_mag_bp')
rp7=dat7.field('median_mag_rp')
p7=dat7.field("parallax")
ag7 = dat7.field("ag_gspphot") 

print(g7)
print(bp7)
print(rp7)
print(p7)
print(ag7)

N7 = len(g7)

color7=np.zeros(N7,float)
G7=np.zeros(N7,float)

i=0
for i in range(N7):
        color7[i]=bp7[i]-rp7[i]
        G7[i] = g7[i] +5 *( 1+np.log10 (p7[i]/1000) ) - ag7[i]




### plots
plt.figure(figsize=(15, 10))
plt.scatter(color5,G5,color='navy',marker = '.',  label = 'eclipsing binaries')
plt.scatter(color2,G2,color='purple',marker = '.',  label = 'RR Lyrae')
plt.scatter(color6,G6,color='red',marker = '.',  label = ' short time variables')
plt.scatter(color3,G3,color='forestgreen',marker = '.',  label = 'MS oscillators')
plt.scatter(color4,G4,color='orange',marker = '.',  label = 'solar-like variables')
plt.scatter(color1,G1,color='aqua',marker = '.',  label = 'planetary transits')
plt.scatter(color7,G7,color='black',marker = '.',  label = 'cepheids')


plt.legend(fontsize="15")
plt.gca().invert_yaxis() 
plt.xlabel('BP-RP')
plt.ylabel('G')
plt.title('aCMD LOPN1')
#p.axis([0,1.4,26,18])
plt.show()


###############################################

### absolute Color-Magnitude Diagram for LOPS2 variable objects



#planetary transits
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





#rrlyrae
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


#MS oscillator
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



#rotation modulation
hdulist4 = fits.open('rotation_modulation2.fits')
hdulist4.info()
dat4 = hdulist4[1].data

g4=dat4.field('median_mag_g_fov')
bp4=dat4.field('median_mag_bp')
rp4=dat4.field('median_mag_rp')
p4=dat4.field("parallax")
ag4=dat4.field("ag_gspphot")

print(g4)
print(bp4)
print(rp4)
print(p4)
print(ag4)

N = len(g4)

color4=np.zeros(N,float)
G4 = np.zeros(N,float)
i=0
for i in range(N):
        color4[i]= bp4[i]-rp4[i]
        G4[i] = g4[i] +5 *( 1+np.log10 (p4[i]/1000) ) - ag4[i]
  
print(G4)



#eclipsing binaries
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



#short timescale
hdulist6 = fits.open('short_timescale2.fits')
hdulist6.info()
dat6 = hdulist6[1].data

g6=dat6.field('median_mag_g_fov')
bp6=dat6.field('median_mag_bp')
rp6=dat6.field('median_mag_rp')
p6=dat6.field("parallax")
ag6=dat6.field("ag_gspphot")

print(g6)
print(bp6)
print(rp6)
print(p6)
print(ag6)

N6 = len(g6)

color6=np.zeros(N6,float)
G6=np.zeros(N6,float)

i=0
for i in range(N6):
        color6[i]=bp6[i]-rp6[i]
        G6[i] = g6[i] +5 *( 1+np.log10 (p6[i]/1000) ) - ag6[i]




#long period variables
hdulist7 = fits.open('long_period2.fits')
hdulist7.info()
dat7 = hdulist7[1].data

g7=dat7.field('median_mag_g_fov')
bp7=dat7.field('median_mag_bp')
rp7=dat7.field('median_mag_rp')
p7 =dat7.field('parallax')
ag7=dat7.field("ag_gspphot")

print(g7)
print(bp7)
print(rp7)
print(p7)
print(ag7)

N = len(g7)

color7=np.zeros(N,float)
G7 = np.zeros(N,float)
i=0
for i in range(N):
        color7[i]=bp7[i]-rp7[i]
        G7[i] = g7[i] +5 *( 1+np.log10 (p7[i]/1000) ) - ag7[i]


### plots
plt.figure(figsize=(15, 10))
plt.scatter(color5,G5,color='red',marker = '.',  label = 'eclipsing binaries')
plt.scatter(color2,G2,color='blue',marker = '.',  label = 'RR Lyrae')
plt.scatter(color6,G6,color='yellow',marker = '.',  label = ' short time variables')
plt.scatter(color3,G3,color='green',marker = '.',  label = 'MS oscillators')
plt.scatter(color4,G4,color='orange',marker = '.',  label = 'solar-like variables')
plt.scatter(color1,G1,color='aqua',marker = '.',  label = 'planetary transits')
plt.scatter(color7,g7,color='fuchsia',marker = '.',  label = 'long period variables')

plt.legend(fontsize="15")
plt.gca().invert_yaxis() 
plt.xlabel('BP-RP')
plt.ylabel('G')
plt.title('aCMD LOPS2')
#p.axis([0,1.4,26,18])
plt.show()