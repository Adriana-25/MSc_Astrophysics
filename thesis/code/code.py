#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 13:39:48 2022

@author: adri
"""
from astropy.io import fits
from astropy.table import Table, join
from astropy.utils.data import get_pkg_data_filename
import pandas as pd

# PIC from DR3 in csv format 
# vari-tables from Gaia Archive in fits format 

dir = 'Master_Thesis/'

path1 = "PICtargetDR3/"
path2 = "GDR3_Variability/"
path3 = "GDR3_PIC_output/LOPN1/"
path4 = "GDR3_PIC_output/LOPS2/"


# PIC LOPN1

df1 = pd.read_csv(dir + path1 + 'PICtargetNSRinputLOPN1.csv')
pic_1 = Table.from_pandas(df1)
pic_1.write(dir + path1 +'plato1.fits', format='fits', overwrite = True)
pic_data1 = get_pkg_data_filename(dir + path1 +'plato1.fits')
fits.info(pic_data1)
table1a = Table.read(pic_data1, hdu = 1)
table1a.rename_column('sourceId','source_id')
print(table1a)

# PIC LOPS2

df2 = pd.read_csv(dir + path1 + 'PICtargetNSRinputLOPS2.csv')
pic_2 = Table.from_pandas(df2)
pic_2.write(dir + path1 + 'plato2.fits', format='fits', overwrite = True)
pic_data2 = get_pkg_data_filename(dir + path1 +'plato2.fits')
fits.info(pic_data2)
table1b = Table.read(pic_data2, hdu = 1)
table1b.rename_column('sourceId','source_id')
print(table1b)


# Variability tables from Gaia Archive 
# http://cdn.gea.esac.esa.int/Gaia/gdr3/Variability/
# https://gea.esac.esa.int/archive/  (advanced ADQL search)
#    SELECT  ALL *
#    FROM gaiadr3.vari_classname
#    FROM gaiadr3.nss_two_body_orbit

a = [str("agn-result.fits"), str("planetary_transits-result.fits"), str("microlensing-result.fits"), str("compact_companion-result.fits"), 
     str("cepheid-result.fits"), str("ms_oscillator-result.fits"), str("rrlyrae-result.fits"), str("short_timescale-result.fits"), 
     str("rotation_modulation-result.fits"), str("long_period_variable-result.fits"), str("eclipsing_binary-result.fits"), str("nss_two_body_orbit-results.fits")]


b =[str("agn.fits"), str("planetary_transits.fits"), str("microlensing.fits"), str("compact_companion.fits"), 
     str("cepheid.fits"), str("ms_oscillator.fits"), str("rrlyrae.fits"), str("short_timescale-result.fits"), 
     str("rotation_modulation.fits"), str("long_period_variable.fits"), str("eclipsing_binary.fits"), str("nss_two_body_orbit.fits")]


c = [str("agn.ecsv"), str("planetary_transits.ecsv"), str("microlensing.ecsv"), str("compact_companion.ecsv"), 
     str("cepheid.ecsv"), str("ms_oscillator.ecsv"), str("rrlyrae.ecsv"), str("short_timescale-result.ecsv"), 
     str("rotation_modulation.ecsv"), str("long_period_variable.ecsv"), str("eclipsing_binary.ecsv"), str("nss_two_body_orbit.ecsv")]


for i in range(12):
    path = dir + path2 + a[i]
    gaia_data = get_pkg_data_filename(path)
    fits.info(gaia_data)
    table2 = Table.read(gaia_data, hdu = 1)
    print("\n%s\n" % a[i])
    print(table2)
    
    # find common rows by key
    # my aim is to check which source_id values from table2 (GDR3 Variability types) can be found in table1a and table1b 
    # for LOPN1 and LOPS2 (PLATO Input Catalogue) 
    
    #join
    table1a_table2 = join(table2, table1a)
    print(table1a_table2)
    table1b_table2 = join(table2, table1b)
    print(table1b_table2)
   

    # save the output in an ecsv file    

    table1a_table2.write(dir + path3 + c[i], format = "ascii.ecsv", overwrite = True)  
    table1b_table2.write(dir + path4 + c[i], format = "ascii.ecsv", overwrite = True)  

    # save the output in a fits file 

    table1a_table2.write(dir + path3 + b[i], format = "fits", overwrite = True) 
    table_filename_a = get_pkg_data_filename(dir + path3 + b[i])
    fits.info(table_filename_a)
    table3a = Table.read(table_filename_a, hdu = 1)
    print(table3a)

    table1b_table2.write(dir + path4 + b[i], format = "fits", overwrite = True)
    table_filename_b = get_pkg_data_filename(dir + path4 + b[i])
    fits.info(table_filename_b)
    table3b = Table.read(table_filename_b, hdu = 1)
    print(table3b)



#Results

print("Cross-matching between GDR3 Variability tables and PLATO Input Catalogue\n")
print("%-30s %-30s %-30s" % ("SOS type", "LOPN1 ", "LOPS2 "))
print("%-30s %-30s %-30s" % ("planetary transits", 21, 31))
print("%-30s %-30s %-30s" % ("microlensing", 0, 0))
print("%-30s %-30s %-30s" % ("compact companion", 0, 0))
print("%-30s %-30s %-30s" % ("cepheid", 0, 0))
print("%-30s %-30s %-30s" % ("ms oscillator", 153, 438))
print("%-30s %-30s %-30s" % ("rrlyrae", 84, 78))
print("%-30s %-30s %-30s" % ("short timescale", 2996, 4081))
print("%-30s %-30s %-30s" % ("rotation modulation", 1794, 1849))
print("%-30s %-30s %-30s" % ("agn", 0, 0))
print("%-30s %-30s %-30s" % ("long period variable", 0,1))
print("%-30s %-30s %-30s" % ("eclipsing binary", 2977, 3467))  
print("%-30s %-30s %-30s" % ("nss two body orbit", 28885, 28960 ))    
