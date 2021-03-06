# =============================================================================
# Basics packages
# =============================================================================
import numpy as np
import copy
from scipy.signal import medfilt

# =============================================================================
# Astropy and associated packages
# =============================================================================

# =============================================================================
# KOALA packages
# =============================================================================
# Modular
from modular_koala.ancillary import vprint

# Original
from koala.plot_plot import plot_plot

# =============================================================================




def relative_throughput(rss): 

    # TODO: New dummy function  
    #    * add histogram plot
    #    * Sigma clipping parameter
                            
# =============================================================================
# Copy input RSS for storage the changes implemented in the task   
# =============================================================================
    rss_out = copy.deepcopy(rss)
    
    mean_count = np.nanmean(rss.intensity_corrected, axis=1)
    perc50 = np.nanpercentile(mean_count, 50)
    # self.low_throughput = mean_count < perc5
    # self.high_throughput = mean_count > perc95
    
    rss_out.intensity_corrected = rss_out.intensity_corrected / perc50
    
    return rss_out



def get_from_sky_flat(rss_skyflat, 
                      plot=True,  
                      kernel_throughput=0,
                      verbose = False):
    """
    Get a 2D array with the throughput 2D using a COMBINED skyflat / domeflat.
    It is important that this is a COMBINED file, as we should not have any cosmics/nans left.
    A COMBINED flappy flat could be also used, but that is not as good as the dome / sky flats.

    Parameters
    ----------
    file_skyflat: string
        The fits file containing the skyflat/domeflat
    plot: boolean
        If True, it plots the results
    throughput_2D_file: string
        the name of the fits file to be created with the output throughput 2D
    no_nas: booleans
        If False, it indicates the mask will be built using the nan in the edges
    correct_ccd_deffects: boolean
        If True, it corrects for ccd defects when reading the skyflat fits file
    kernel_throughput: odd number
        If not 0, the 2D throughput will be smoothed with a this kernel
    """


    # Set verbose
    vprint.verbose = verbose
    
    
    n_spectra = rss_skyflat.intensity.shape[0]
    n_wave = rss_skyflat.intensity.shape[1]
    
    throughput_2D_ = np.zeros_like(rss_skyflat.intensity)
    
    
    
    
    vprint("\n> Getting the throughput per wavelength...")
    for i in range(n_wave):
        column = rss_skyflat.intensity_corrected[:, i]
        mcolumn = column / np.nanmedian(column)
        throughput_2D_[:, i] = mcolumn

    if kernel_throughput > 0:
        print("\n  - Applying smooth with kernel =", kernel_throughput)
        throughput_2D = np.zeros_like(throughput_2D_)
        for i in range(rss_skyflat.n_spectra):
            throughput_2D[i] = medfilt(throughput_2D_[i], kernel_throughput)
        #skyflat.RSS_image(image=throughput_2D, chigh=1.1, clow=0.9, cmap="binary_r")
        #skyflat.history.append('- Throughput 2D smoothed with kernel ' + np.str(kernel_throughput))
    else:
        throughput_2D = throughput_2D_


    if plot:
        x = np.arange(n_spectra)
        median_throughput = np.nanmedian(throughput_2D, axis=1)
        plot_plot(x, median_throughput, ymin=0.2, ymax=1.2, hlines=[1, 0.9, 1.1],
                  ptitle="Median value of the 2D throughput per fibre", xlabel="Fibre")
        #skyflat.RSS_image(image=throughput_2D, cmap="binary_r",title="\n ---- 2D throughput ----")

    print("\n> Throughput 2D obtained!")

    return throughput_2D





def apply_throughput(rss, 
                  throughput,
                  verbose=False,
                  plot=True):
    """
    Apply throughput_2D using the information of a variable or a fits file.
    """
    # Set print verbose
    vprint.verbose = verbose
# =============================================================================
# Copy input RSS for storage the changes implemented in the task   
# =============================================================================
    rss_out = copy.deepcopy(rss)
    try:
        vprint("\n> Applying 2D throughput correction")

    except:
        print("Warning: Incorrect throughput specified. Correction not applied")
    
    """
    if plot:
        print("\n> Plotting map BEFORE correcting throughput:")
        self.RSS_image()
    """
    

    rss_out.intensity_corrected = rss_out.intensity_corrected / throughput        
    rss_out.variance_corrected = rss_out.variance_corrected / throughput**2
    
    
    """
    if plot:
        print("  Plotting map AFTER correcting throughput:")
        self.RSS_image()
    """
    return rss_out
















