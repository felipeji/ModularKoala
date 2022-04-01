# =============================================================================
# Basics packages
# =============================================================================
import numpy as np
from scipy import interpolate

# =============================================================================
# Astropy and associated packages
# =============================================================================
from astropy.io import fits
# =============================================================================



def vprint(*arg):
    """
    Prints the arguments only if vprint.verbose is previously setted True.
    """
    vprint.verbose
    if vprint.verbose:
        for i in arg:
            print(i,end=' ')


def airmass_from_header(header):
    """
    Compute the airmass extracting the parameters from KOALAS's header'
    """
    # Get ZD, airmass
    ZDSTART = header['ZDSTART']
    ZDEND = header['ZDEND']
    ZD = (ZDSTART + ZDEND) / 2
    airmass = 1 / np.cos(np.radians(ZD))
    return airmass        




def interpolate_nan(spectrum):
    """
    Replace any NaN value in a spectrum by interpolating non-NaNs sideways.
    """
    
    inds = np.arange(spectrum.shape[0])
    good = np.where(np.isfinite(spectrum))
    f = interpolate.interp1d(inds[good], spectrum[good],bounds_error=False)
    filled_spectrum = np.where(np.isfinite(spectrum),spectrum,f(inds))
    return filled_spectrum


def nearest(array, value):
    """
    Returns the index of the element in <array> nearest to <value>.
    
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def py_koala_header(header):
    """
    Copy 2dfdr headers values from extensions 0 and 2 needed for the initial 
    header for PyKoala. (based in the header constructed in  save_rss_fits in 
    koala.io)
    """
    
    # To fit actual PyKoala header format
    header.rename_keyword('CENRA','RACEN')
    header.rename_keyword('CENDEC','DECCEN')
    
    cards = [header.cards['BITPIX'],
            header.cards["ORIGIN"],                        
            header.cards["TELESCOP"],  
            header.cards["ALT_OBS"],              
            header.cards["LAT_OBS"],                
            header.cards["LONG_OBS"], 
            header.cards["INSTRUME"],      
            header.cards["GRATID"],         
            header.cards["SPECTID"],                                
            header.cards["DICHROIC"],
            header.cards['OBJECT'],    
            header.cards["EXPOSED"],
            header.cards["ZDSTART"], 
            header.cards["ZDEND"],                                           
            header.cards['NAXIS'],                       
            header.cards['NAXIS1'],                 
            header.cards['NAXIS2'],                  
            header.cards['RACEN'],
            header.cards['DECCEN'], 
            header.cards['TEL_PA'],
            header.cards["CTYPE2"],  
            header.cards["CUNIT2"],     
            header.cards["CTYPE1"],  
            header.cards["CUNIT1"],    
            header.cards["CRVAL1"],
            header.cards["CDELT1"],
            header.cards["CRPIX1"],
            header.cards["CRVAL2"],  
            header.cards["CDELT2"],
            header.cards["CRPIX2"],
            ]
    py_koala_header = fits.header.Header(cards=cards, copy=False)

    return py_koala_header


def py_koala_spaxels_table(spaxels_table):
    """
    Generates the spaxels tables needed for PyKoala from the 2dfdr spaxels table.
    """
    # Filtering only selected (in use) fibres 
    spaxels_table = spaxels_table[spaxels_table['SELECTED']==1]
    
    # Defining new arrays
    arr1 = np.arange(len(spaxels_table)) + 1 # +  for starting in 1
    arr2 = np.ones(len(spaxels_table))
    arr3 = np.ones(len(spaxels_table))
    arr4 = np.ones(len(spaxels_table)) * 2048
    arr5 = np.zeros(len(spaxels_table))
    arr6 = spaxels_table['XPOS']
    arr7 = spaxels_table['YPOS']
    arr8 = spaxels_table['SPEC_ID']
         
    # Defining new columns
    col1 = fits.Column(name='Fibre', format='I', array=arr1)
    col2 = fits.Column(name='Status', format='I', array=arr2)
    col3 = fits.Column(name='Ones', format='I', array=arr3)
    col4 = fits.Column(name='Wavelengths', format='I', array=arr4)
    col5 = fits.Column(name='Zeros', format='I', array=arr5)
    col6 = fits.Column(name='Delta_RA', format='D', array=arr6)
    col7 = fits.Column(name='Delta_Dec', format='D', array=arr7)
    col8 = fits.Column(name='Fibre_OLD', format='I', array=arr8)
    
    # PyKoala Spaxels table
    py_koala_spaxels_table = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6,col7,col8])

    return py_koala_spaxels_table














