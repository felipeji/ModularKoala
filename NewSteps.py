# =============================================================================
# Basics Modules 
# =============================================================================
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors

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

NOT IMPLEMENTED YET: 8. Sky subtraction S

9. Correcting negative values N Sube algunas cuentas (negativos) del cielo

NOT IMPLEMENTED YET: 10. Emission lines identification E 


11. Correcting small CCD / sky residuals R  (Cosmic)
12. Saving the processed / cleaned RSS file.


List and description of the sample RSS in the red arm

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



# Save figures to PDF
pdf = PdfPages('test_output_data/figures.pdf')



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
# Original rss
# =============================================================================
plt.close('all')
plt.figure(1,figsize=(14,8))
plt.suptitle('Original RSS',size=25)
ax1 = plt.subplot(2,2,1)
ax1.set_title('Combined sky flat')
skyflat_rss.show()

ax2 = plt.subplot(2,2,2,sharex=ax1,sharey=ax1)
ax2.set_title('Standard star (HD60753)')
std_rss.show()

ax3 = plt.subplot(2,2,3,sharex=ax1,sharey=ax1)
sci_rss[0].show()    
ax3.set_title(sci_rss[0].header['OBJECT'])

ax4 = plt.subplot(2,2,4,sharex=ax1,sharey=ax1)
ax4.set_title(sci_rss[1].header['OBJECT'])
sci_rss[1].show()    
plt.tight_layout()

pdf.savefig(1)


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
# NaNs and edges corrected
# =============================================================================
plt.close('all')
plt.figure(1,figsize=(14,8))
plt.suptitle('NaNs and edges corrected',size=25)
ax1 = plt.subplot(2,2,1)
ax1.set_title('Combined sky flat')
skyflat_clean.show()

ax2 = plt.subplot(2,2,2,sharex=ax1,sharey=ax1)
ax2.set_title('Standard star (HD60753)')
std_clean.show()

ax3 = plt.subplot(2,2,3,sharex=ax1,sharey=ax1)
sci_clean[0].show()    
ax3.set_title(sci_rss[0].header['OBJECT'])

ax4 = plt.subplot(2,2,4,sharex=ax1,sharey=ax1)
ax4.set_title(sci_rss[1].header['OBJECT'])
sci_clean[1].show()    
plt.tight_layout()

pdf.savefig(1)







# =============================================================================
# 4. Fixing small wavelength shifts 
# =============================================================================
# We compute and apply the wavelenght drift solution on the science RSS and then 
# we apply one of this solution to the skyflat and the standadr star.   
from modular_koala.fix_wavelengths import fix_wavelengths


# Loop over all the science objects
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

# std
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
std_throughput = apply_throughput(std_drift,throughput = throughput_2D,verbose=True)

# Loop over all the science objects
sci_throughput = []
for rss in sci_drift:
    sci_throughput.append( apply_throughput(rss,throughput = throughput_2D,verbose=True) )



# =============================================================================
# 6. Correcting for extinction 
# =============================================================================
from modular_koala.extinction import extinction

# This function compute the extinction from the KOALA headers
from modular_koala.ancillary import airmass_from_header


# std
airmass = airmass_from_header(std_throughput.header)
std_extinction = extinction(std_throughput,airmass=airmass,verbose=True,plot=True)

# Loop over all the science objects
sci_extinction = []
for rss in sci_throughput:
    airmass = airmass_from_header(rss.header)
    sci_extinction.append( extinction(rss,airmass = airmass,verbose=True,plot=True) )


     

# =============================================================================
# Throughput and extinction corrected
# =============================================================================
plt.close('all')
plt.figure(1,figsize=(14,8))
plt.suptitle('Throughput and extinction corrected',size=25)

ax1 = plt.subplot(2,2,1)
ax1.set_title('Combined sky flat')
ax1.text(0.5,0.5,'NO CORRECTED',ha='center',va='center',size = 20)
ax1.axis('off')

ax2 = plt.subplot(2,2,2)
ax2.set_title('Standard star (HD60753)')
std_extinction.show()

ax3 = plt.subplot(2,2,3)
sci_extinction[0].show()    
ax3.set_title(sci_rss[0].header['OBJECT'])

ax4 = plt.subplot(2,2,4)
ax4.set_title(sci_rss[1].header['OBJECT'])
sci_extinction[1].show()    
plt.tight_layout()

pdf.savefig(1)






          

# =============================================================================
# 7. Telluric correction (only for red data) U
# =============================================================================
from modular_koala.tellurics import Tellurics

# This class compute the telluric correction from a std star. Once computed we 
# can apply this solution to RSS or save the correction to a file.
 
telluric_correction = Tellurics(std_extinction, 
                                 exclude_wlm= [[6245,6390],[6450,6750],[6840,7000],
                                               [7140,7400],[7550,7720],[8050,8450]],
                                 weight_fit_median = 1)

# This is the step for applying the computed telluric correction to the RSS
# rss_telluric = tellurric_correction.apply(XXX_extinction)


# Comparison between the two solutions. 
# Computed
w = std_extinction.wavelength
corr = telluric_correction.telluric_correction

# From file
w_txt = np.genfromtxt('/Users/felipe/DataCentral/KOALA/sample_RSS/385R/telluric_correction_20180227_385R.dat')[:,0]
corr_txt = np.genfromtxt('/Users/felipe/DataCentral/KOALA/sample_RSS/385R/telluric_correction_20180227_385R.dat')[:,1]

plt.close('all')
plt.figure(1)
plt.plot(w,corr,label='Computed correction')
plt.plot(w_txt,corr_txt,label='From telluric_correction_20180227_385R.dat')
plt.legend()

# Some diferecnces between the existing and the computed solution.
# Are the parameters considered in telluric_correction correct?
# We apply the solution from file to the science RSS.

# We read an existing telluric solution
from modular_koala.tellurics import tellurics_from_file
telluric_file = '/Users/felipe/DataCentral/KOALA/sample_RSS/385R/telluric_correction_20180227_385R.dat'




# Loop over all the science objects
sci_tellurics = []
for rss in sci_extinction:
    sci_tellurics.append( tellurics_from_file(rss, telluric_file=telluric_file) )



# =============================================================================
# Tellurics corrected
# =============================================================================
plt.close('all')
plt.figure(1,figsize=(14,8))
plt.suptitle('Tellurics correction',size=25)

ax1 = plt.subplot(2,2,1)
ax1.set_title('Combined sky flat')
ax1.text(0.5,0.5,'NO CORRECTED',ha='center',va='center',size = 20)
ax1.axis('off')

ax2 = plt.subplot(2,2,2)
ax2.set_title('Standard star (HD60753)')
ax2.text(0.5,0.5,'NO CORRECTED',ha='center',va='center',size = 20)
ax2.axis('off')

ax3 = plt.subplot(2,2,3)
sci_tellurics[0].show()    
ax3.set_title(sci_rss[0].header['OBJECT'])

ax4 = plt.subplot(2,2,4)
ax4.set_title(sci_rss[1].header['OBJECT'])
sci_tellurics[1].show()    
plt.tight_layout()

pdf.savefig(1)

















# =============================================================================
# 9. Correcting negative values 
# =============================================================================
from modular_koala.clean_residuals import extreme_negatives

# Loop over all the science objects
sci_extreme = []
for rss in sci_tellurics:
    sci_extreme.append( extreme_negatives(rss, 
                                          fibre_list='all',
                                          percentile_min=0.5,
                                          plot=True,
                                          verbose=True)
                        )


# =============================================================================
# 11. Correcting small CCD / sky residuals R  (Cosmic)
# =============================================================================
# This is a new (test) implementation of the task using L.A.Cosmic
# Cosmics

# # 2D version
# from modular_koala.clean_residuals import kill_cosmics_2D
# for rss in sci_extreme:
#     sci_cosmic.append( kill_cosmics_2D(rss, cleantype='idw',sigclip=15, objlim=7,verbose=True) )


from modular_koala.clean_residuals import kill_cosmics

Halpha = 6564.6 # Wavelengh of the brightest line 

sci_cosmic = []
for rss in sci_extreme:
    sci_cosmic.append( kill_cosmics(rss, Halpha) )



cosmic_mask = []
for rss in sci_cosmic:
    cosmic_mask.append(rss.mask_layer(4))




# =============================================================================
# Cosmics ray corrected
# =============================================================================
plt.close('all')
plt.figure(1,figsize=(14,8))
plt.suptitle('Extreme negatives and cosmic rays correction',size=25)

ax1 = plt.subplot(2,2,1)
ax1.set_title('Combined sky flat')
ax1.text(0.5,0.5,'NO CORRECTED',ha='center',va='center',size = 20)
ax1.axis('off')

ax2 = plt.subplot(2,2,2)
ax2.set_title('Standard star (HD60753)')
ax2.text(0.5,0.5,'NO CORRECTED',ha='center',va='center',size = 20)
ax2.axis('off')

ax3 = plt.subplot(2,2,3)
sci_cosmic[0].show()    
ax3.set_title(sci_rss[0].header['OBJECT'])

ax4 = plt.subplot(2,2,4)
ax4.set_title(sci_rss[1].header['OBJECT'])
sci_cosmic[1].show()    
plt.tight_layout()

pdf.savefig(1)




# =============================================================================
# Mask
# =============================================================================
plt.close('all')
plt.figure(1,figsize=(14,8))
plt.suptitle('Masks',size=25)

ax1 = plt.subplot(2,2,1)
ax1.set_title('Combined sky flat')
skyflat_clean.show_mask()


ax2 = plt.subplot(2,2,2)
ax2.set_title('Standard star (HD60753)')
std_throughput.show_mask()


ax3 = plt.subplot(2,2,3)
ax3.set_title(sci_rss[0].header['OBJECT'])
sci_cosmic[0].show_mask()


ax4 = plt.subplot(2,2,4)
ax4.set_title(sci_rss[1].header['OBJECT'])
sci_cosmic[1].show_mask()

plt.tight_layout()

pdf.savefig(1)









# =============================================================================
# Just saving the colorbar
# =============================================================================
plt.close('all')
plt.figure(1,figsize=(14,8))
cmap = colors.ListedColormap(['darkgreen','blue','red','darkviolet','black','orange'])
figure = plt.imshow(sci_clean[0].mask,vmin=0,vmax=6,  cmap=cmap, interpolation='none')
figure.remove()
plt.axis('off')
cb = plt.colorbar(figure,location="top",aspect=20)
ticks_position = np.arange(6)+0.5
ticks_name = ['From file','Blue edge','Red edge','NaN\'s','Cosmics','Extreme \n negatives']
cb.set_ticks(ticks_position)
cb.set_ticklabels(ticks_name)
pdf.savefig(1)
pdf.close()





# =============================================================================
# 12. Saving the processed / cleaned RSS file.
# =============================================================================
for rss in sci_cosmic:
    # We save the fist with the name of the object + CLEAN
    path = 'test_output_data'
    name = rss.header['OBJECT'].replace(' ','_')+'_CLEAN.fits'
    file_path = os.path.join(path,name)
    rss.tofits(file_path)    






"""

# # =============================================================================
# # Checking header
# # =============================================================================

 
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

"""



#Checking the values of skylines in





# # TODO: Why add again this information
# py_koala_spaxels_table.header['CENRA']  =  rss.RA_centre_deg  / ( 180/np.pi )   # Must be in radians 
# py_koala_spaxels_table.header['CENDEC']  =  rss.DEC_centre_deg / ( 180/np.pi )

# fits_image_hdu = fits.PrimaryHDU(data)
# # TO BE DONE    
# errors = [0]  ### TO BE DONE                
# error_hdu = fits.ImageHDU(errors)


# hdu_list = fits.HDUList([fits_image_hdu,error_hdu, header2_hdu]) 

# hdu_list.writeto(fits_file, overwrite=True) 




# """









