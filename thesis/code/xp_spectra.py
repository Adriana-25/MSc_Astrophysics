#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 27 13:39:40 2023

@author: adri
"""

from astropy.io import fits
from astropy.table import Table
from astropy.utils.data import get_pkg_data_filename
import matplotlib.pyplot as plt
import numpy as np
import math

plt.rcParams['font.size'] = '20'


# Sampled BP/RP (XP) mean spectra:  Spectral (energy) distribution to investigate multifrequential variability
#	W m**-2 nm**-1 as a function of nm  ( \nu F (\nu) over \lambda )
    

### RR Lyrae
# source_id 1985008784705347456
# this is a RRc star pulsating with fundamental period p1 = 0.38031 


hdulist = fits.open('/media/adri/ANNA/Master_Thesis/xp_spectra/RR/XP_SAMPLED-Gaia DR3 1985008784705347456.fits')

hdulist.info()
dat = hdulist[1].data


wave=dat.field('wavelength')
flux = dat.field("flux")
flux_err=dat.field("flux_error")

plt.figure(figsize=(15, 10)) 
plt.scatter(wave,flux, c='blue', s = 80)
plt.errorbar(wave, flux, flux_err, lw=1, capsize=4)
plt.text(x = 800, y = 3.7*1e-16, s = "RRc Lyrae Gaia DR3 1985008784705347456", horizontalalignment='center', color = 'gray')

plt.legend(fontsize = 15, loc = "best")
plt.xlabel("$\lambda$ (nm)")
plt.ylabel(r'$ \mathrm{ \nu F(\nu) (erg \cdot m^{-2} s^{-1}) }$')
#plt.ylim(-30,25) 
plt.title('XP mean spectrum')
plt.show()




# source_id 863983115982280448
# this is a RRab star pulsating with fundamental period pf = 0.57183 

hdulist = fits.open('/media/adri/ANNA/Master_Thesis/xp_spectra/RR/XP_SAMPLED-Gaia DR3 863983115982280448.fits')

hdulist.info()
dat = hdulist[1].data


wave=dat.field('wavelength')
flux = dat.field("flux")
flux_err=dat.field("flux_error")

plt.figure(figsize=(15, 10)) 
plt.scatter(wave,flux, c='blue', s = 80)
plt.errorbar(wave, flux, flux_err, lw=1, capsize=4)
plt.text(x = 800, y = 3.7*1e-16, s = "RRab Lyrae Gaia DR3 863983115982280448", horizontalalignment='center', color = 'gray')
plt.legend(fontsize = 15, loc = "best")
plt.xlabel("$\lambda$ (nm)")
plt.ylabel(r'$ \mathrm{ \nu F(\nu) (erg \cdot m^{-2} s^{-1}) }$')
#plt.ylim(-30,25) 
plt.title('XP mean spectrum')
plt.show()





### SOLAR-LIKE ROTATION MODULATION
# Gaia DR3 393450310336818176

hdulist = fits.open('/media/adri/ANNA/Master_Thesis/xp_spectra/SOLAR-LIKE/XP_SAMPLED-Gaia DR3 393450310336818176.fits')

hdulist.info()
dat = hdulist[1].data


wave=dat.field('wavelength')
flux = dat.field("flux")
flux_err=dat.field("flux_error")

plt.figure(figsize=(15, 10)) 
plt.scatter(wave,flux, c='blue', s = 80)
plt.errorbar(wave, flux, flux_err, lw=1, capsize=4)
plt.text(x = 800, y = 2.15*1e-15, s = "Solar-like Gaia DR3 393450310336818176", horizontalalignment='center', color = 'gray')
plt.legend(fontsize = 15, loc = "best")
plt.xlabel("$\lambda$ (nm)")
plt.ylabel(r'$ \mathrm{ \nu F(\nu) (erg \cdot m^{-2} s^{-1}) }$')
#plt.ylim(-30,25) 
plt.title('XP mean spectrum')
plt.show()



### MAIN SEQUENCE OSCILLATORS
# Gaia DR3 393326065521087616

hdulist = fits.open('/media/adri/ANNA/Master_Thesis/xp_spectra/MS/XP_SAMPLED-Gaia DR3 393326065521087616.fits')

hdulist.info()
dat = hdulist[1].data


wave=dat.field('wavelength')
flux = dat.field("flux")
flux_err=dat.field("flux_error")

plt.figure(figsize=(15, 10)) 
plt.scatter(wave,flux, c='blue', s = 80)
plt.errorbar(wave, flux, flux_err, lw=1, capsize=4)
plt.text(x = 800, y = 4.25*1e-15, s = "MS oscillator Gaia DR3 393326065521087616", horizontalalignment='center', color = 'gray')
plt.legend(fontsize = 15, loc = "best")
plt.xlabel("$\lambda$ (nm)")
plt.ylabel(r'$ \mathrm{ \nu F(\nu) (erg \cdot m^{-2} s^{-1}) }$')
#plt.ylim(-30,25) 
plt.title('XP mean spectrum')
plt.show()




#### PLANETARY TRANSITS
# Gaia DR3 1864885215233116032

hdulist = fits.open('/media/adri/ANNA/Master_Thesis/xp_spectra/PLANETARY/XP_SAMPLED-Gaia DR3 1864885215233116032.fits')

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
plt.ylabel(r'$ \mathrm{ \nu F(\nu) (erg \cdot m^{-2} s^{-1}) }$')
#plt.ylim(-30,25) 
plt.title('XP mean spectrum')
plt.show()






#### ECLIPSING BINARIES
# Gaia DR3 393454948901343360

hdulist = fits.open('/media/adri/ANNA/Master_Thesis/xp_spectra/EB/XP_SAMPLED-Gaia DR3 393454948901343360.fits')

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
plt.ylabel(r'$ \mathrm{ \nu F(\nu) (erg \cdot m^{-2} s^{-1}) }$')
#plt.ylim(-30,25) 
plt.title('XP mean spectrum')
plt.show()




#### SHORT TIMESCALE
# Gaia DR3 4594961454434214144


hdulist = fits.open('/media/adri/ANNA/Master_Thesis/xp_spectra/SHORT/XP_SAMPLED-Gaia DR3 4594961454434214144.fits')

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
plt.ylabel(r'$ \mathrm{ \nu F(\nu) (erg \cdot m^{-2} s^{-1}) }$')
#plt.ylim(-30,25) 
plt.title('XP mean spectrum')
plt.show()





#### CEPHEIDS
# Gaia DR3 4301612233202961024

hdulist = fits.open('/media/adri/ANNA/Master_Thesis/xp_spectra/CEPHEID/XP_SAMPLED-Gaia DR3 4301612233202961024.fits')

hdulist.info()
dat = hdulist[1].data


wave=dat.field('wavelength')
flux = dat.field("flux")
flux_err=dat.field("flux_error")

plt.figure(figsize=(15, 10)) 
plt.scatter(wave,flux, c='blue', s = 80)
plt.errorbar(wave, flux, flux_err, lw=1, capsize=4)
plt.text(x = 800, y = 2.75*1e-15, s = "Cepheid Gaia DR3 4301612233202961024", horizontalalignment='center', color = 'gray')
plt.legend(fontsize = 15, loc = "best")
plt.xlabel("$\lambda$ (nm)")
plt.ylabel(r'$ \mathrm{ \nu F(\nu) (erg \cdot m^{-2} s^{-1}) }$')
#plt.ylim(-30,25) 
plt.title('XP mean spectrum')
plt.show()






#### LONG PERIOD : NO XP SPECTRUM FOUND in Gaia Archive, only from CDS Portal
# Gaia DR3 3005222633152619904










