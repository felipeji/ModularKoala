# =============================================================================
# Basics packages
# =============================================================================
from scipy import interpolate
import numpy as np
import copy
import os


# =============================================================================
# Astropy and associated packages
# =============================================================================

# =============================================================================
# KOALA packages
# =============================================================================
# Modular
from modular_koala.ancillary import vprint

# =============================================================================








def apply_correction(rss,
                     airmass,
                     extinction_correction, 
                     extinction_file="",
                     verbose=True):
    """
    This corrects for extinction using a given extiction correction and an observatory file

    Parameters
    ----------
    extinction_correction: array
        array with the extiction correction derived using the airmass

    observatory_extinction_file : string, input data file (default = 'ssoextinct.dat')
        file containing the standard extinction curve at the observatory
    verbose : boolean (default = True)
        Print results
    """
    
    # Set verbose
    vprint.verbose = verbose
    
    
    vprint("  Intensities corrected for extinction stored in self.intensity_corrected")
    vprint("  Variance corrected for extinction stored in self.variance_corrected")

    rss.intensity_corrected *= extinction_correction[np.newaxis, :]
    rss.variance_corrected *= extinction_correction[np.newaxis, :]**2
    
    
    comment = ' '.join(["- Data corrected for extinction using file :",
                        extinction_file,
                        "Average airmass =",
                        str(airmass)])
    
    rss.log['extinction'] = comment 
    



    


def extinction(rss,
               airmass,
               extinction_file=None,
               verbose=True,
               # plot=False, Necesitamos redefinir los limite de wavelenght validos usando las mascaras
               ):
    """
    This task accounts and corrects for extinction due to gas and dusty
    between target and observer. It creates a extinction curve based off
    airmass input and observatory file data.

            Parameters
    ----------
    apply_extinction : boolean (default = True)
        Apply extinction curve to the data
    observatory_extinction_file : string, input data file (default = 'ssoextinct.dat')
        file containing the standard extinction curve at the observatory
    plot : boolean (default = True)
        Plot that generates extinction curve [extinction_curve_wavelengths,extinction_corrected_airmass]
    verbose : boolean (default = True)
        Print results
    """
    
    # Set verbose
    vprint.verbose = verbose

        
# =============================================================================
# Copy of the input RSS for containing the changes implemented by the task   
# =============================================================================    
    rss_out = copy.deepcopy(rss)
    
    
    vprint("\n> Computing extinction at given airmass...")
    vprint("  Airmass = ", np.round(airmass, 3))
    
    
    # Read data
    if extinction_file is None:
        
        dirname = os.path.dirname(__file__)
        extinction_file = os.path.join(dirname, 'input_data/observatory_extinction/ssoextinct.dat')

        

    data_observatory = np.loadtxt(extinction_file, unpack=True)
    extinction_curve_wavelenghts = data_observatory[0]
    extinction_curve = data_observatory[1]
    extinction_corrected_airmass = 10 ** (0.4 * airmass * extinction_curve)
    # Make fit
    tck = interpolate.splrep(extinction_curve_wavelenghts,
                             extinction_corrected_airmass, s=0)
    extinction_correction = interpolate.splev(rss_out.wavelength, tck, der=0)

    if verbose: print("  Observatory file with extinction curve :\n ", extinction_file)


    # if plot:
    #     cinco_por_ciento = 0.05 * (np.max(extinction_correction) - np.min(extinction_correction))
    #     plot_plot(extinction_curve_wavelenghts, extinction_corrected_airmass, xmin=np.min(rss_out.wavelength),
    #               xmax=np.max(rss_out.wavelength), ymin=np.min(extinction_correction) - cinco_por_ciento,
    #               ymax=np.max(extinction_correction) - cinco_por_ciento,
    #               vlines=[self.valid_wave_min, self.valid_wave_max],
    #               ptitle='Correction for extinction using airmass = ' + str(np.round(self.airmass, 3)),
    #               xlabel="Wavelength [$\mathrm{\AA}$]", ylabel="Flux correction", fig_size=fig_size,
    #               statistics=False)

# =============================================================================
# Apply_extinction
# =============================================================================
    apply_correction(rss = rss_out,
                     airmass = airmass,
                     extinction_correction = extinction_correction, 
                     extinction_file = extinction_file,
                     verbose = verbose
                     )
    


    return rss_out