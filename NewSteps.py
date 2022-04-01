# =============================================================================
# Basics Modules 
# =============================================================================
import os
import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# Astropy
# =============================================================================
from astropy.wcs import WCS
from astropy.io import fits


# =============================================================================
# Originals Clean steps
# =============================================================================
"""
1. Reading the data and establishing properties of the Python object,
2. Creating a mask with the valid data (masking the edges of the CCD). 
3. Cleaning the CCD defects NaN y Negativos
4. Fixing small wavelength shifts W
5. Applying throughput 2D T
6. Correcting for extinction X
7. Telluric correction (only for red data) U

8. Sky subtraction S

9. Correcting negative values N Sube algunas cuentas (negativos) del cielo
10. Emission lines identification E 
11. Correcting small CCD / sky residuals R  (Cosmic)
12. Saving the processed / cleaned RSS file.


List and description of RSS 

modular_koala/input_data/sample_RSS/385R/27feb20025red.fits HD60753 C 3W WIDE
modular_koala/input_data/sample_RSS/385R/27feb20026red.fits HD60753 B 3S WIDE
modular_koala/input_data/sample_RSS/385R/27feb20027red.fits HD60753 A WIDE

modular_koala/input_data/sample_RSS/385R/27feb20028red.fits HILT600 A
modular_koala/input_data/sample_RSS/385R/27feb20029red.fits HILT600 B 3S
modular_koala/input_data/sample_RSS/385R/27feb20030red.fits HILT600 C 3W

modular_koala/input_data/sample_RSS/385R/27feb20031red.fits Tol30 A PA0

modular_koala/input_data/sample_RSS/385R/27feb20032red.fits He2-10 B 1.5S
modular_koala/input_data/sample_RSS/385R/27feb20033red.fits He2-10 C 1.5S 3E
modular_koala/input_data/sample_RSS/385R/27feb20034red.fits He2-10 D 18S 3E
modular_koala/input_data/sample_RSS/385R/27feb20035red.fits He2-10 E 1.5E
modular_koala/input_data/sample_RSS/385R/27feb20036red.fits He2-10 F 1.5E 1.5S

modular_koala/input_data/sample_RSS/385R/combined_skyflat_red.fits SKYFLAT
"""






# =============================================================================
# 1. Reading the data
# =============================================================================
# Main class for reading RSS files into RSS class
from modular_koala.rss import read_rss 

# Ancilliary function to match the current structure of a PyKoala header
from modular_koala.ancillary import py_koala_header 

# Ancilliary function to match the current structure of a PyKoala spaxels table
from modular_koala.ancillary import py_koala_spaxels_table

# COMMENT: The code repits any time we read a RSS file. But at the end of the day
# for reading KOALA RSS we are only telling the task a file path. In that way
# we keep RSS class simple and useful for any RSS data (only file_path and WCS
# are required arguments). At the same time is easy to write a function to read
# KOALA RSS, incluiding headers and spaxels table passing only file_path argument



# Standard star: HD60753 ------------------------------------------------------
file_path = 'modular_koala/input_data/sample_RSS/385R/27feb20025red.fits'

# Constructing PyKoala header from 2dfdr (header[0] and header[2]) 
header = fits.getheader(file_path,0) + fits.getheader(file_path,2)
koala_header = py_koala_header(header)

# WCS
koala_wcs = WCS(header)

# Constructing Pykoala Spaxels table from 2dfdr spaxels table (data[2])
spaxels_table = fits.getdata(file_path,2)
koala_spaxels_table = py_koala_spaxels_table(spaxels_table)

# List of bad spaxels from 2dfdr spaxels table
bad_spaxels_list = (spaxels_table['SPEC_ID'][spaxels_table['SELECTED']==0]-1).tolist() 
# -1 to start in 0 rather than in 1


# Read RSS file into PyKoala structure
std_rss = read_rss(file_path, 
                   wcs= koala_wcs, 
                   bad_spaxels_list = bad_spaxels_list,
                   header=koala_header,
                   spaxels_table=koala_spaxels_table,
                   )


# combined_skyflat_red --------------------------------------------------------
file_path = 'modular_koala/input_data/sample_RSS/385R/combined_skyflat_red.fits'

# Constructing PyKoala header from 2dfdr (header[0] and header[2]) 
header = fits.getheader(file_path,0) + fits.getheader(file_path,2)
koala_header = py_koala_header(header)

# WCS
koala_wcs = WCS(header)

# Constructing Pykoala Spaxels table from 2dfdr spaxels table (data[2])
spaxels_table = fits.getdata(file_path,2)
koala_spaxels_table = py_koala_spaxels_table(spaxels_table)

# List of bad spaxels from 2dfdr spaxels table
bad_spaxels_list = (spaxels_table['SPEC_ID'][spaxels_table['SELECTED']==0]-1).tolist()
# -1 to start in 0 rather than in 1


# Read RSS file into PyKoala structure
skyflat_rss = read_rss(file_path, 
                   wcs= koala_wcs, 
                   bad_spaxels_list = bad_spaxels_list,
                   header=koala_header,
                   spaxels_table=koala_spaxels_table,
                   )


# Object: Tol30 ---------------------------------------------------------------
file_path = 'modular_koala/input_data/sample_RSS/385R/27feb20031red.fits'

# Constructing PyKoala header from 2dfdr (header[0] and header[2]) 
header = fits.getheader(file_path,0) + fits.getheader(file_path,2)
koala_header = py_koala_header(header)

# WCS
koala_wcs = WCS(header)

# Constructing Pykoala Spaxels table from 2dfdr spaxels table (data[2])
spaxels_table = fits.getdata(file_path,2)
koala_spaxels_table = py_koala_spaxels_table(spaxels_table)

# List of bad spaxels from 2dfdr spaxels table
bad_spaxels_list = (spaxels_table['SPEC_ID'][spaxels_table['SELECTED']==0]-1).tolist() 
# -1 to start in 0 rather than in 1


# Read RSS file into PyKoala structure
tol30_rss = read_rss(file_path, 
                   wcs= koala_wcs, 
                   bad_spaxels_list = bad_spaxels_list,
                   header=koala_header,
                   spaxels_table=koala_spaxels_table,
                   )

# Object: He2_100 -------------------------------------------------------------
file_path = 'modular_koala/input_data/sample_RSS/385R/27feb20032red.fits'

# Constructing PyKoala header from 2dfdr (header[0] and header[2]) 
header = fits.getheader(file_path,0) + fits.getheader(file_path,2)
koala_header = py_koala_header(header)

# WCS
koala_wcs = WCS(header)

# Constructing Pykoala Spaxels table from 2dfdr spaxels table (data[2])
spaxels_table = fits.getdata(file_path,2)
koala_spaxels_table = py_koala_spaxels_table(spaxels_table)

# List of bad spaxels from 2dfdr spaxels table
bad_spaxels_list = (spaxels_table['SPEC_ID'][spaxels_table['SELECTED']==0]-1).tolist() 
# -1 to start in 0 rather than in 1


# Read RSS file into PyKoala structure
he2_100_rss = read_rss(file_path, 
                   wcs= koala_wcs, 
                   bad_spaxels_list = bad_spaxels_list,
                   header=koala_header,
                   spaxels_table=koala_spaxels_table,
                   )

#------------------------------------------------------------------------------

# List of RSS corresponding to objects to iterate (as an example) further.
# More elegant to use dict here, buy we keep simple for geting the main idea.
 
sci_rss = [tol30_rss,he2_100_rss]


# =============================================================================
# 2. Creating a mask with the valid data (masking the edges of the CCD).  
# =============================================================================
# Edges are automatically detected. Only NaNs pixels are corrected
from modular_koala.clean_residuals import fix_edges

# skyflat
skyflat_edges = fix_edges(skyflat_rss,verbose=True)
skyflat_edges.show()


# standard star
std_edges = fix_edges(std_rss,verbose=True)
std_edges.show()


# Iterating over all the science objects
sci_edge = []

for rss in sci_rss:
    sci_edge.append( fix_edges(rss,verbose=True) )



# =============================================================================
# 3. Cleaning the CCD defects 
# =============================================================================
# NaNs
from modular_koala.clean_residuals import clean_nans


# skyflat
skyflat_clean = clean_nans(skyflat_edges,verbose=True)
skyflat_clean.show()


# standard star
std_clean = clean_nans(std_edges,verbose=True)
std_clean.show()


# Iterating over all the science objects
sci_clean = []

for rss in sci_edge:
    sci_clean.append( clean_nans(rss,verbose=True) )



# =============================================================================
# 4. Fixing small wavelength shifts 
# =============================================================================
# We compute and apply the wavelenght drift solution on the science RSS and then 
# we apply one of this solution to the skyflat and the std   
from modular_koala.fix_wavelengths import fix_wavelengths


# Iterating over all the science objects
sci_drift = []

for rss in sci_clean:
    sci_drift.append( fix_wavelengths(rss, grating='385R',verbose=True) )


# Compare the solution in the two science RSS
from modular_koala.fix_wavelengths import compare_fix_wavelengths
compare_fix_wavelengths(sci_drift[0],sci_drift[1])


# We apply the wavelength correction from Tol30 rss (i.e. first rss in sci_drift )
# to skyflat and std
sol = sci_drift[0].log['wavelenght fix']['sol']

# skyflat
skyflat_drift = fix_wavelengths(skyflat_clean, grating='385R',sol = sol) 

#std
std_drift = fix_wavelengths(std_clean, grating='385R',sol = sol) 





# =============================================================================
# 5. Applying throughput
# =============================================================================
# Constructing the 2D throughput from sky_flat
from modular_koala.throughput import get_from_sky_flat
throughput_2D = get_from_sky_flat(skyflat_drift)














# We apply the throughput to the std and the science rss
from modular_koala.throughput import apply_throughput

# From 2D array
rss_throughput = apply_throughput(fixed_wl,throughput = throughput,verbose=True)
rss_throughput.tofits('test_throughput.fits',overwrite=True)

# From fits
throughput_2D_path = "test_throughput.fits"
throughput = fits.getdata(throughput_2D_path)
rss_throughput_from_file = apply_throughput(fixed_wl,throughput = throughput)





# =============================================================================
# 6. Correcting for extinction X
# =============================================================================
from modular_koala.extinction import extinction
from modular_koala.ancillary import airmass_from_header

airmass = airmass_from_header(header)

rss_extintion = extinction(rss_throughput,
               airmass,
               extinction_file=None,
               verbose=True)
               

plt.close('all')

plt.figure(1)
rss_extintion.show()


# =============================================================================
# 7. Telluric correction (only for red data) U
# =============================================================================
# Read standard star
# Reading sky_flat

# 
file_path = "/Users/felipe/DataCentral/KOALA/sample_RSS/385R/27feb20027red.fits"



# Parameters
header = fits.getheader(file_path)
koala_wcs = WCS(header)

spaxels_table = fits.getdata(file_path,2)

bad_spaxels_list = np.arange(len(spaxels_table['SELECTED']))[~np.bool_(spaxels_table['SELECTED'])].tolist()       



# Reading data
 #print("\n> Reading a COMBINED skyflat / domeflat to get the 2D throughput...")
standard = read_rss(file_path, 
                    wcs= koala_wcs,
                    bad_spaxels_list = bad_spaxels_list)


# Applying wavelength correction
sol = fixed_wl.log['wavelenght fix']['sol']
standard_wl = fix_wavelengths(standard, grating='385R',sol = sol) 

# Correnting edges and NaN
standard_edges = fix_edges(standard_wl,verbose=True)
standard_clean = clean_nans(standard_edges)




from modular_koala.tellurics import Tellurics

tellurric_correction = Tellurics(rss_extintion, exclude_wlm= [[6245,6390],[6450,6850],[6840,7000],[7140,7400],[7550,7720],[8050,9000]])



rss_telluric = tellurric_correction.apply(rss_extintion)

from modular_koala.tellurics import tellurics_from_file

telluric_file = '/Users/felipe/DataCentral/KOALA/sample_RSS/385R/telluric_correction_20180227_385R.dat'
rss_telluric_from_file = tellurics_from_file(rss_clean_nans, telluric_file=telluric_file)


w = np.genfromtxt('/Users/felipe/DataCentral/KOALA/sample_RSS/385R/telluric_correction_20180227_385R.dat')[:,0]

t = np.genfromtxt('/Users/felipe/DataCentral/KOALA/sample_RSS/385R/telluric_correction_20180227_385R.dat')[:,1]

plt.close('all')

plt.figure(1)
plt.plot(w,t)
plt.plot(w,tellurric_correction.telluric_correction,label='new')

# plt.close('all')

# plt.figure(1)
# rss_extintion.show()    

# plt.figure(2)
# rss_telluric.show()


# plt.figure(3)
# rss_telluric_from_file.show()




# =============================================================================
# Checing header
# =============================================================================
"""

# Valores agregados por PyKoala
                
        header["SOL0"] = sol[0]
        header["SOL1"] = sol[1]
        header["SOL2"] = sol[2]
     
        header['DESCRIP'] = description

    for item in rss.history_RSS:
        if item == "- Created fits file (this file) :":
            header['HISTORY'] = "- Created fits file :"
        else:    
            header['HISTORY'] = item        
    header['FILE_IN'] = rss.filename     
    header['HISTORY'] = '-- RSS processing using PyKOALA '+ version
    #header['HISTORY'] = 'Developed by Angel Lopez-Sanchez, Yago Ascasibar, Lluis Galbany et al.'
    #header['HISTORY'] =  version #'Version 0.10 - 12th February 2019'    
    now=datetime.datetime.now()
    
    header['HISTORY'] = now.strftime("File created on %d %b %Y, %H:%M:%S using input file:")
    header['DATE'] = now.strftime("%Y-%m-%dT%H:%M:%S") #'2002-09-16T18:52:44'   # /Date of FITS file creation
    #header['HISTORY'] = 'using input file:'
    header['HISTORY'] = rss.filename

    for item in rss.history:
        header['HISTORY'] = item

    header['HISTORY'] = "- Created fits file (this file) :"
    header['HISTORY'] = " "+fits_file   
    header['FILE_OUT'] = fits_file





# TODO: Why add again this information
py_koala_spaxels_table.header['CENRA']  =  rss.RA_centre_deg  / ( 180/np.pi )   # Must be in radians 
py_koala_spaxels_table.header['CENDEC']  =  rss.DEC_centre_deg / ( 180/np.pi )

fits_image_hdu = fits.PrimaryHDU(data)
# TO BE DONE    
errors = [0]  ### TO BE DONE                
error_hdu = fits.ImageHDU(errors)


hdu_list = fits.HDUList([fits_image_hdu,error_hdu, header2_hdu]) 

hdu_list.writeto(fits_file, overwrite=True) 




"""




# =============================================================================
# 11. Correcting small CCD / sky residuals R  (Cosmic)
# =============================================================================

#    
# # Cosmics
# from modular_koala.clean_residuals import kill_cosmics

# rss_cosmic = kill_cosmics(rss_clean_nans,sigclip=15)


# plt.close('all')

# plt.figure(1)
# vmin,vmax = np.nanpercentile(rss_cosmic.intensity,(5,95))
# plt.imshow(rss_cosmic.intensity,vmin=vmin,vmax=vmax,cmap='viridis')
    

# plt.figure(2)
# vmin,vmax = np.nanpercentile(rss_cosmic.intensity_corrected,(5,95))
# rss_cosmic.show(vmin=vmin,vmax=vmax,cmap='viridis')
    


# plt.figure(3)
# cosmic_mask = rss_cosmic.mask_layer(4)
# plt.imshow(cosmic_mask)
# 





