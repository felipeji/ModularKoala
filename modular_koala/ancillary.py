# =============================================================================
# Basics packages
# =============================================================================
import numpy as np
from scipy import interpolate

# =============================================================================
# Astropy and associated packages
# =============================================================================

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

