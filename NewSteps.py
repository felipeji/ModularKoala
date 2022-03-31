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
"""









# =============================================================================
# 1. Reading the data
# =============================================================================
from modular_koala.rss import read_rss

# Data path
path = "/Users/felipe/DataCentral/KOALA/sample_RSS/385R/"
file_path = os.path.join(path, "27feb20033red.fits")


# Basics parameter for reading RSS
# Parameters
header = fits.getheader(file_path)
koala_wcs = WCS(header)

spaxels_table = fits.getdata(file_path,2)

bad_spaxels_list = np.arange(len(spaxels_table['SELECTED']))[~np.bool_(spaxels_table['SELECTED'])].tolist()       


# Reading data
rss_raw = read_rss(file_path, 
                   wcs= koala_wcs, 
                   bad_spaxels_list = bad_spaxels_list,
                   header=header
                   )

rss_raw.show()


# =============================================================================
# 2. Masking the edges of the CCD and replacing NaN with the median. 
# =============================================================================

from modular_koala.clean_residuals import fix_edges

rss_edges = fix_edges(rss_raw,verbose=True)
rss_edges.show()



# =============================================================================
# 3. Cleaning the CCD defects 
# =============================================================================
# NaNs
from modular_koala.clean_residuals import clean_nans
rss_clean_nans = clean_nans(rss_edges)


"""
plt.figure(1)
rss_edges.show()

plt.figure(2)
rss_clean_nans.show()
"""

# =============================================================================
# 4. Fixing small wavelength shifts 
# =============================================================================
# Edges
from modular_koala.fix_wavelengths import fix_wavelengths_edges
fixed_wl_edge = fix_wavelengths_edges(rss_clean_nans)

# Fix wavelengths 
from modular_koala.fix_wavelengths import fix_wavelengths
fixed_wl = fix_wavelengths(rss_clean_nans, grating='385R') 

# Comparison
from modular_koala.fix_wavelengths import compare_fix_wavelengths
compare_fix_wavelengths(fixed_wl_edge,fixed_wl)


# =============================================================================
# 5. Applying throughput
# =============================================================================
# Relative
from modular_koala.throughput import relative_throughput

rss_relative_throughput = relative_throughput(fixed_wl) 


# Constructing the 2D throughput from sky_flat
# Reading sky_flat
file_path = "/Users/felipe/DataCentral/KOALA/sample_RSS/385R/throughput_2D_20180227_385R.fits"

# Parameters
header = fits.getheader(file_path)
koala_wcs = WCS(header)


# Reading data
 #print("\n> Reading a COMBINED skyflat / domeflat to get the 2D throughput...")

skyflat = read_rss(file_path, 
                   wcs= koala_wcs)
# Applying wavelength correction
sol = fixed_wl.log['wavelenght fix']['sol']
skyflat_wl = fix_wavelengths(skyflat, grating='385R',sol = sol) 

# Correnting edges and NaN
skyflat_edges = fix_edges(skyflat_wl,verbose=True)
skyflat_clean = clean_nans(skyflat_edges)


from modular_koala.throughput import get_from_sky_flat
throughput = get_from_sky_flat(skyflat_clean)




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

header = fits.getheader('/Users/felipe/DataCentral/KOALA/sample_RSS/385R/27feb20031red.fits')

# Valores ya existentes en el header
    header['BITPIX']  =  16  #TODO: En el header pone -32, se fuerza el cambio por algo?
    header["ORIGIN"]  = 'AAO'    #    / Originating Institution                        
    header["TELESCOP"]= 'Anglo-Australian Telescope'    # / Telescope Name  
    header["ALT_OBS"] =                 1164 # / Altitude of observatory in metres              
    header["LAT_OBS"] =            -31.27704 # / Observatory latitude in degrees                
    header["LONG_OBS"]=             149.0661 # / Observatory longitude in degrees 
    header["INSTRUME"] = "AAOMEGA-KOALA"             # / Instrument in use  
    header["GRATID"]  = '385R'      # / Disperser ID 
    header["SPECTID"] = 'RD'                        # / Spectrograph ID                                
    header["DICHROIC"]= 'X5700'                        # / Dichroic name   ---> CHANGE if using X6700!!    
    header['OBJECT'] = 'HILT600 A'
    header["EXPOSED"] = 120
    header["ZDSTART"]= 35.5060963215581 
    header["ZDEND"]= 35.3457657340589
    header['NAXIS']   =   2                              # / number of array dimensions                       
    header['NAXIS1']  =   2048                 
    header['NAXIS2']  =   1000                 
    header['TEL_PA'] = 89.9177910093928
    header["CTYPE2"] = 'Fibre number'          # / Label for axis 2  
    header["CUNIT2"] = ' '           # / Units for axis 2     
    header["CTYPE1"] = 'Wavelength'          # / Label for axis 2  
    header["CUNIT1"] = 'Angstroms'           # / Units for axis 2     

    header["CRVAL1"] = 7692.370367828 #  / Co-ordinate value of axis 2
    header["CDELT1"] = 1.57518231234 # 
    header["CRPIX1"] = 1024 # 1024. / Reference pixel along axis 2
    header["CRVAL2"] = 5.000000000000E-01 # / Co-ordinate value of axis 2  
    header["CDELT2"] = 1.000000000000E+00 # / Co-ordinate increment along axis 2
    header["CRPIX2"] = 1.000000000000E+00 # / Reference pixel along axis 2 
 
    
    # Aca se cambia el valor de la cabecera solo para RA y DEC (originalmente CENRA, CENDEC) 
    # Sería mejor respetar esto y dentro de la rutina llamar segun este valor o el que sea que defina RA y Dec en la cabecera.
    header['RAcen'] = rss.RA_centre_deg #TODO: not in header  ???
    header['DECcen'] = rss.DEC_centre_deg # TODO: not in header ???
 


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







# Spaxels Table 

header2_all_fibres = fits.getdata('/Users/felipe/DataCentral/KOALA/sample_RSS/385R/27feb20031red.fits',2)



    # Header 2 with the RA and DEC info!    
    header2_all_fibres = rss.header2_data  
    header2_good_fibre = []
    header2_original_fibre = []
    header2_new_fibre = []
    header2_delta_RA=[]
    header2_delta_DEC=[]
    header2_2048 =[]
    header2_0 =[]
    
    fibre = 1
    for i in range (len(header2_all_fibres)):
        if header2_all_fibres[i][1]  == 1:
            header2_original_fibre.append(i+1)
            header2_new_fibre.append(fibre)
            header2_good_fibre.append(1)
            header2_delta_RA.append(header2_all_fibres[i][5])
            header2_delta_DEC.append(header2_all_fibres[i][6])
            header2_2048.append(2048)
            header2_0.append(0)
            fibre = fibre + 1
     
#    header2_=[header2_new_fibre, header2_good_fibre, header2_good_fibre, header2_2048, header2_0,  header2_delta_RA,  header2_delta_DEC,  header2_original_fibre]
#    header2 = np.array(header2_).T.tolist()   
#    header2_hdu = fits.ImageHDU(header2)
            
    col1 = fits.Column(name='Fibre', format='I', array=np.array(header2_new_fibre))
    col2 = fits.Column(name='Status', format='I', array=np.array(header2_good_fibre))
    col3 = fits.Column(name='Ones', format='I', array=np.array(header2_good_fibre))
    col4 = fits.Column(name='Wavelengths', format='I', array=np.array(header2_2048))
    col5 = fits.Column(name='Zeros', format='I', array=np.array(header2_0))
    col6 = fits.Column(name='Delta_RA', format='D', array=np.array(header2_delta_RA))
    col7 = fits.Column(name='Delta_Dec', format='D', array=np.array(header2_delta_DEC))
    col8 = fits.Column(name='Fibre_OLD', format='I', array=np.array(header2_original_fibre))
    
    cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8])
    header2_hdu = fits.BinTableHDU.from_columns(cols)
    
    header2_hdu.header['CENRA']  =  rss.RA_centre_deg  / ( 180/np.pi )   # Must be in radians
    header2_hdu.header['CENDEC']  =  rss.DEC_centre_deg / ( 180/np.pi )
    
    hdu_list = fits.HDUList([fits_image_hdu,error_hdu, header2_hdu]) #  hdu_list = fits.HDUList([fits_image_hdu, wavelengths_hdu, flux_correction_hdu])

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





