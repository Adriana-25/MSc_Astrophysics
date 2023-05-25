#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 18:28:32 2023

@author: adri
"""

### All-sky Aitoff projection in Galactic coordinates



from astropy.io import fits    
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import Angle
import numpy as np

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


hdulist3 = fits.open('eclipsing_binaries1.fits')
hdulist3.info()
dat3 = hdulist3[1].data

l3=dat3.field('l')
b3=dat3.field('b')


l3 = Angle(l3 * u.deg)
l3.wrap_at('180d', inplace=True)

print(l3)


b3 = Angle(b3 * u.deg)
b3.wrap_at('90d', inplace=True)

print(b3)

l3 = np.deg2rad(np.array(l3))
b3 = np.deg2rad(np.array(b3))


hdulist4 = fits.open('eclipsing_binaries2.fits')
hdulist4.info()
dat4 = hdulist4[1].data

l4=dat4.field('l')
b4=dat4.field('b')


l4 = Angle(l4 * u.deg)
l4.wrap_at('180d', inplace=True)

print(l4)


b4 = Angle(b4 * u.deg)
b4.wrap_at('90d', inplace=True)

print(b4)

l4 = np.deg2rad(np.array(l4))
b4 = np.deg2rad(np.array(b4))


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
b8 = np.deg2rad(np.array(b8
                         ))


hdulist9 = fits.open('short_timescale1.fits')
hdulist9.info()
dat9 = hdulist9[1].data

l9=dat9.field('l')
b9=dat9.field('b')


l9 = Angle(l9 * u.deg)
l9.wrap_at('180d', inplace=True)

print(l9)


b9 = Angle(b9 * u.deg)
b9.wrap_at('90d', inplace=True)

print(b9)

l9 = np.deg2rad(np.array(l9))
b9 = np.deg2rad(np.array(b9))


hdulist10 = fits.open('short_timescale2.fits')
hdulist10.info()
dat10 = hdulist10[1].data

l10=dat10.field('l')
b10=dat10.field('b')


l10 = Angle(l10 * u.deg)
l10.wrap_at('180d', inplace=True)

print(l10)


b10 = Angle(b10 * u.deg)
b10.wrap_at('90d', inplace=True)

print(b10)

l10 = np.deg2rad(np.array(l10))
b10 = np.deg2rad(np.array(b10))

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


hdulist13 = fits.open('long_period2.fits')
hdulist13.info()
dat13 = hdulist13[1].data

l13=dat13.field('l')
b13=dat13.field('b')


l13 = Angle(l13 * u.deg)
l13.wrap_at('180deg', inplace=True)

print(l13)


b13 = Angle(b13 * u.deg)
b13.wrap_at('90d', inplace=True)

print(b13)

l13 = np.deg2rad(np.array(l13))
b13 = np.deg2rad(np.array(b13))


hdulist14 = fits.open('cepheid1.fits')
hdulist14.info()
dat14 = hdulist14[1].data

l14=dat14.field('l')
b14=dat14.field('b')


l14 = Angle(l14 * u.deg)
l14.wrap_at('180deg', inplace=True)

print(l14)


b14 = Angle(b14 * u.deg)
b14.wrap_at('90d', inplace=True)

print(b14)

l14 = np.deg2rad(np.array(l14))
b14 = np.deg2rad(np.array(b14))


fig = plt.figure(figsize=(15,12))
ax = fig.add_subplot(1,1,1, projection='aitoff')

ax.scatter(l, b, s=3, color='red', alpha=1, label = 'planetary transits')   
ax.scatter(l2, b2, s=3, color='red', alpha=1)
ax.scatter(l3, b3, s=2, color='navy', alpha=0.5, label = 'eclipsing binaries')   
ax.scatter(l4, b4, s=2, color='navy', alpha=0.5)
ax.scatter(l5, b5, s=3, color='green', alpha=1, label = 'RR Lyrae')   
ax.scatter(l6, b6, s=3, color='green', alpha=1)   
ax.scatter(l7, b7, s=3, color='magenta', alpha=1, label = 'MS oscillators')   
ax.scatter(l8, b8, s=3, color='magenta', alpha=1)   
ax.scatter(l9, b9, s=2, color='aqua', alpha=1, label = 'short-timescale variables')   
ax.scatter(l10, b10, s=2, color='aqua', alpha=1)   
ax.scatter(l11, b11, s=2, color='yellow', alpha=1, label = 'solar-like variables')   
ax.scatter(l12, b12, s=2, color='yellow', alpha=1)   
ax.scatter(l13, b13, s=3, color='black', alpha=1, label = 'long-period variables')   
ax.scatter(l14, b14, s=3, color='orange', alpha=1, label = 'Cepheids')   

ax.grid()
plt.legend(fontsize="10", loc ="lower right")
plt.title('All-sky map in Galactic coordinates\n')
