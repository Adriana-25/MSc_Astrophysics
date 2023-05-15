#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 15:31:52 2023

@author: adri
"""



from astropy.io import fits
from astropy.table import Table, join
from astropy.utils.data import get_pkg_data_filename
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import scipy
import matplotlib.colors as colors
from astropy import units as u
from astropy.coordinates import Angle
plt.rcParams['font.size'] = '20'



### All-sky Aitoff projection in Galactic coordinates


hdulist11 = fits.open('rotation_modulation1.fits')
hdulist11.info()
dat11 = hdulist11[1].data

l11=dat11.field('l')
b11=dat11.field('b')


l11 = Angle(l11 * u.deg)
l11.wrap_at('180d', inplace=True)

print(l11)


b11 = Angle(b11 * u.deg)
b11.wrap_at('90d', inplace=True)

print(b11)

l11 = np.deg2rad(np.array(l11))
b11 = np.deg2rad(np.array(b11))


hdulist12 = fits.open('rotation_modulation2.fits')
hdulist12.info()
dat12 = hdulist12[1].data

l12=dat12.field('l')
b12=dat12.field('b')


l12 = Angle(l12 * u.deg)
l12.wrap_at('180d', inplace=True)

print(l12)


b12 = Angle(b12 * u.deg)
b12.wrap_at('90d', inplace=True)

print(b12)

l12 = np.deg2rad(np.array(l12))
b12 = np.deg2rad(np.array(b12))



fig = plt.figure(figsize=(15,12))
ax = fig.add_subplot(1,1,1, projection='aitoff')
 
ax.scatter(l11, b11, s=5, marker='o',color='red', alpha=1, label = 'solar-like variables in LOPN1')   
ax.scatter(l12, b12, s=5, marker='o',color='blue', alpha=1, label = 'solar-like variables in LOPS2')   

ax.grid()
plt.legend(fontsize="15", loc ="lower right")
plt.title('All-sky map in Galactic coordinates\n')



########################

### absolute Color-Magnitude Diagram for solar-like variables


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




plt.figure(figsize=(15, 10))
plt.scatter(color4,G4,color='red',marker = '.',  label = 'solar-like variables in LOPN1')



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




plt.scatter(color4,G4,color='blue',marker = '.',  label = 'solar-like variables in LOPS2')

plt.legend(fontsize="15")
plt.gca().invert_yaxis() 
plt.xlabel('BP-RP')
plt.ylabel('G')
plt.title('aCMD')
#p.axis([0,1.4,26,18])
plt.show()





######## period - amplitude (activity) diagram for solar-like stars displaying rotational modulation 

# Two families: shorter and longer period: 
# High/Low amplitude Fast Rotators and Low Amplitude Slow Rotators



hdulist1 = fits.open('rotation_modulation1.fits')
hdulist1.info()
dat1 = hdulist1[1].data

p1=dat1.field('best_rotation_period')
a1=dat1.field("max_activity_index_g")


hdulist2 = fits.open('rotation_modulation2.fits')
hdulist2.info()
dat2 = hdulist2[1].data

p2=dat2.field('best_rotation_period')
a2=dat2.field("max_activity_index_g")



P=[]
A=[]
i=0
j=0

for i in range(len(p1)):
    P.append(p1[i])
    A.append(a1[i])

for j in range(len(p2)):
    P.append(p2[j])
    A.append(a2[j])


N = len(P)
P =np.array(P)
A =np.array(A)

    
min=P[0]    
for i in range (1,len(P)):
    if(P[i]<min):
        min=P[i]
   
        
print("The minimum of P is: ", min)


max=P[0]
for j in range(1,len(P)):
   if(P[j]>max):
        max=P[j]
        
print("The maximum of P is: ", max)
        
#nbin = (max - min)/N
#print(nbin)           
            
min2=A[0]    
for i in range (1,len(P)):
    if(A[i]<min2):
        min2=A[i]
   
        
print("The minimum of A is: ", min2)


max2=A[0]
for j in range(1,len(P)):
   if(A[j]>max2):
        max2=A[j]
        
print("The maximum of A is: ", max2)


plt.figure(figsize=(15, 10)) 

plt.hist2d(P, A, bins=60, norm=colors.LogNorm(), cmap="inferno")
plt.colorbar()
plt.xlim(-0.5,45)
plt.ylim(-0.0001,0.15)
plt.xlabel("P (d)")
plt.ylabel("A (mag)")
plt.title("P-A diagram")
plt.tight_layout()
plt.show()












######## epoch photometry for solar-like stars with rotational modulation 

hdulist = fits.open('EPOCH_PHOTOMETRY-Gaia DR3 393450310336818176.fits')
                                                                                              ##   Gaia DR3 
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




### unfolded

p = 3.59365  # best rotation period from tables for 393450310336818176 solar like variable 

plt.figure(figsize=(15, 10)) 
plt.scatter(t_rp,RP,color='red',marker = 'o',  label = 'RP band')
plt.scatter(t_g,G,color='green',marker = 'o',  label = 'G band')
plt.scatter(t_bp,BP,color='blue',marker = 'o',  label = 'BP band')

plt.text(x = 2200, y = 10, s = "solar-like 393450310336818176 P=3.594 (d)", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
plt.ylim(11,9.8)
plt.legend(fontsize = 15, loc = "best")
plt.xlabel("Time (BJD)")
plt.ylabel('Magnitude (mag)')
plt.title('Epoch photometry')
plt.show()


plt.figure(figsize=(15, 10))

### folded full time series

 
t_rp=np.array(t_rp)        # time of first transit in RP band 
phase =  ((t_rp - 1697.08228 - p/2) / p) % 1 - 0.5
plt.scatter(phase,RP,s= 80,color='red',marker = 'o',  label = 'RP band') 
#plt.scatter(phase,RP,s= 80, c= phase, cmap = 'autumn',  label = 'RP band') 



t_g=np.array(t_g)        # time of first transit in G band 
phase =  ((t_g - 1697.08192 - p/2) / p) % 1 - 0.5
plt.scatter(phase,G,s= 80,color='green',marker = 'o',  label = 'G band') 
#plt.scatter(phase,G,s= 80, c= phase, cmap = 'summer',  label = 'G band') 



t_bp=np.array(t_bp)        # time of first transit in BP band 
phase =  ((t_bp - 1697.08219 - p/2) / p) % 1 -0.5
plt.scatter(phase,BP,s= 80, color='blue',marker = 'o',  label = 'BP band') 
#plt.scatter(phase,BP,s= 80, c= phase, cmap = 'winter',  label = 'BP band') 



plt.text(x = 0, y = 10, s = "solar-like 393450310336818176 P=3.594 (d)", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
plt.ylim(11,9.8)
plt.legend(fontsize = 15, loc = "best")
plt.xlabel("Phase")
plt.ylabel('Magnitude (mag)')
plt.title('Epoch photometry')
plt.show()



### The rotational modulation signal looses coherence across the full time-series because of the intrinsic
### evolution of MARs and can be detected only in the shorter sub-series.
###  The segmentation of time-series is needed because the spots due to the stellar magnetic activity have a life-time 
### shorter than the whole Gaia time-series. The rotational modulation induced by spots can therefore be detected only 
### in segments whose duration is comparable with the spots life-time. 


#segments_rotation_period   
#segments_start_time
#segments_end_time



### unfolded full time series with zoom on G and vertical lines delimiting the segments

plt.figure(figsize=(15, 10)) 
plt.scatter(t_g,G,color='green',marker = 'o',  label = 'G band')
#plt.text(x = 2200, y = 10, s = "solar-like 393450310336818176 P=3.59365 (d)", horizontalalignment='center', color = 'gray')
plt.gca().invert_yaxis()
plt.ylim(10.68,10.64)
'''
plt.axvline(x = 2054.52916475, color = 'b', label = 'segment 5')    # p = 3.59615166
plt.axvline(x = 2151.73314403, color = 'b')

plt.axvline(x = 1826.3264414366593, color = 'black', label = 'segment 1')  # p = 0.58436493
plt.axvline(x = 1944.876554855351, color = 'black')

plt.axvline(x =1880.55698048 , color = 'r', label = 'segment 2')  # p = 0.59543512 
plt.axvline(x =1998.04344347, color = 'r')
'''
plt.axvline(x = 1998.04344347, color = 'purple', label = 'segment 4')  # p = 3.53940262
plt.axvline(x =2103.72930378, color = 'purple')
plt.text(x = 2200, y = 10.645, s = "solar-like 393450310336818176 P=3.539 (d)", horizontalalignment='center', color = 'gray')

plt.legend(fontsize = 15, loc = "best")
plt.xlabel("Time (BJD)")
plt.ylabel('Magnitude (mag)')
plt.title('Epoch photometry')
plt.show()





### folded fourth segment plus sinusoidal fit
# start 1998.04344347
# end 2103.72930378
# period 3.53940262


t_g_seg=[]
g_seg=[]

for i in range(24):
    t_g_seg.append(t_g[i+43])
    g_seg.append(G[i+43])


print(t_g_seg)
print(g_seg)


plt.figure(figsize=(15, 10)) 
p = 3.53940262
t_g_seg=np.array(t_g_seg)        # time of first transit in this segment in G band
phase =  ((t_g_seg - 1998.04344347 - p/2) / p) % 1
print(phase)

plt.scatter(phase,g_seg,color='purple',marker = '*', s=200, label = 'G band') 
plt.gca().invert_yaxis()
plt.ylim(10.68,10.64)
plt.text(x = 0.5, y = 10.645, s = "solar-like 393450310336818176 P=3.539 (d)", horizontalalignment='center', color = 'gray')
plt.xlim(0,1)
plt.xlabel("Phase")
plt.ylabel('Magnitude (mag)')
plt.title('Epoch photometry 4th segment')

time = np.arange(-3, 3, 0.01);
amplitude = (-np.sin(0.25*time+1.44)+11.659)
plt.plot(time, amplitude, label = 'fit ')
plt.legend(fontsize = 15, loc = "best")

plt.show()
