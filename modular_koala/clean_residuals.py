# =============================================================================
# Basics packages
# =============================================================================
import copy
import numpy as np

# =============================================================================
# Astropy and associated packages
# =============================================================================
from ccdproc import cosmicray_lacosmic

# =============================================================================
# KOALA packages
# =============================================================================
# Modular
from modular_koala.ancillary import vprint,interpolate_nan, nearest
from modular_koala.mask import mask_section
from modular_koala.rss import detect_edge



# =============================================================================




def fix_edges(rss,
              median_from=8800,
              median_to=6300,              
              verbose=False):
    """
    Fixing edges of a RSS file edges
    """
    # Set print verbose
    vprint.verbose = verbose
    
# =============================================================================
# Copy input RSS for storage the changes implemented in the task   
# =============================================================================
    rss_out = copy.deepcopy(rss)
    wavelength = rss.wavelength
    
    
    # Detect edges of the detector
    min_w, min_index,  max_w, max_index = detect_edge(rss)

    # median_to and median_from in index
    median_to_index = nearest(wavelength,median_to)
    median_from_index = nearest(wavelength,median_from)
    
# =============================================================================
# Blue edge
# =============================================================================
    lower = 0
    upper = min_index
    lower_mask = mask_section(rss, 
                          lower=lower,
                          upper=upper,
                          mask_value = 'nan',
                          verbose = verbose,
                          )
    
    rss_out.mask += (lower_mask * 2**1) # This mask has index 1 (2^1, i.e. maks value = 2)

    for i, f in enumerate(rss.intensity):
        median_value = np.nanmedian(f[min_index:median_to_index])
    
        rss_out.intensity_corrected[i][lower_mask[i]] = median_value
        
    comment = "\n - Blue edge has been checked up to " + str(min_w)+' Angstrom (pixel '+str(min_index)+')'
    rss_out.log['blue edge']['comment'] = comment
    vprint(comment)

        # if plot: plot_plot(w,[f,ff],vlines=[median_to,fix_to], xmax=median_to)
            
# =============================================================================
# Red edge
# =============================================================================
    lower = max_index
    upper = len(wavelength)
    upper_mask = mask_section(rss, 
                          lower=lower,
                          upper=upper,
                          mask_value = 'nan',
                          verbose = verbose,
                          )
    
    rss_out.mask += (upper_mask * 2**2) # This mask has index 2 (2^2, i.e. maks value = 4)

    for i, f in enumerate(rss.intensity):
        median_value = np.nanmedian(f[median_from_index:max_index])

        rss_out.intensity_corrected[i][upper_mask[i]] = median_value
        
        # if plot: plot_plot(w,[f,ff],vlines=[median_to,fix_to], xmax=median_to)

    comment = "\n - Red edge has been checked from " + str(max_w)+' Angstrom (pixel '+str(max_index)+')'
    rss_out.log['red edge']['comment'] = comment
    vprint(comment)
        

    return rss_out
    
"""    
    if plot:
        self.RSS_image(title=" - Before correcting edges") # PLOT BEFORE                                                                      
        self.RSS_image(title=" - After correcting edges") # PLOT AFTER


"""



def clean_nans(rss,
              verbose=False):
    # Set print verbose
    vprint.verbose = verbose
# =============================================================================
# Copy input RSS for storage the changes implemented in the task   
# =============================================================================
    rss_out = copy.deepcopy(rss)
    
    rss_out.mask = np.isnan(rss.intensity_corrected) * 8  # This mask has index 3 (2^3, i.e. maks value = 8)
    for index, fiber in enumerate(rss.intensity_corrected):
        rss_out.intensity_corrected[index] = interpolate_nan(fiber)

    return rss_out





def kill_cosmics(rss, 
                 verbose = False,
                 **kwargs
                 ):
    """
    """
    # Set print verbose
    vprint.verbose = verbose

# =============================================================================
# Copy input RSS for storage the changes implemented in the task   
# =============================================================================
    rss_out = copy.deepcopy(rss)
    
    
# =============================================================================
# Construct a sky model from the 20 faintest fibers for subtracting in the 
# cosmic ray detection proccess   
# =============================================================================
    # Find 20 brigthtest fibers for constructing the sky spectra template 
    
    masked_data = np.ma.masked_array(rss.intensity,rss.mask)
    
    integrated_fibre = np.nansum(masked_data,axis=1).data # appliying the mask

    # brightest_20_fibers = (-integrated_fibre).argsort()[:20]
    # brightest_line_pixel = np.nanargmax(median_spectrum) # Take 10 brigtest and then median ??
    
    faintest_20_fibers = integrated_fibre.argsort()[:20]
        
    sky_median_model = np.nanmedian(rss.intensity_corrected[faintest_20_fibers],0)
                                                         
    # 2D sky model by scaling the sky_median model to every fiber   
    sky_2D_model = np.zeros_like(rss.intensity_corrected)
    
    for index, fiber in enumerate(rss.intensity_corrected):
        proportion = sky_median_model / fiber
        scale_factor = np.nanmedian(proportion)
        sky_2D_model[index] = sky_median_model * scale_factor 
    
    
    corrected_data, cosmic_mask = cosmicray_lacosmic(rss.intensity_corrected,
                                       gain_apply=False,
                                       inbkg=sky_2D_model,
                                       **kwargs) 

    rss_out.intensity_corrected = corrected_data    
    rss_out.mask = cosmic_mask * 16  # This mask has index 4 (2^4, i.e. maks value = 16)
    
    return rss_out




def extreme_negatives(rss, 
                      fibre_list='all',
                      percentile_min=0.5,
                      plot=True,
                      verbose=True):
    """
    Remove pixels that have extreme negative values (that is below percentile_min) and replace for the median value

    Parameters
    ----------
    fibre_list : list of integers (default all)
        List of fibers to clean. The default is [], that means it will do everything.
    percentile_min : float, (default = 0.5)
        Minimum value accepted as good.
    plot : boolean (default = False)
        Display all the plots
    verbose: boolean (default = True)
        Print what is doing
    """
       
    # Set print verbose
    vprint.verbose = verbose
    

    # =============================================================================
    # Copy input RSS for storage the changes implemented in the task   
    # =============================================================================
    rss_out = copy.deepcopy(rss)

    if fibre_list == 'all':
        fibre_list = list(range(rss_out.intensity.shape[0]))
        vprint("\n> Correcting the extreme negatives in all fibres, making any pixel below")
    else:
        vprint("\n> Correcting the extreme negatives in given fibres, making any pixel below")


    minimo = np.nanpercentile(rss_out.intensity, percentile_min)

    vprint("  np.nanpercentile(intensity_corrected, ", percentile_min, ") = ", np.round(minimo, 2))
    vprint("  to have the median value of the fibre...")

    
     # Masking by fibre 
    for fibre in fibre_list:
        fibre_mask = rss.intensity[fibre]<minimo
        rss_out.intensity_corrected[fibre][fibre_mask] = np.nanmedian(rss_out.intensity[fibre])
        rss_out.mask[fibre][fibre_mask] += 2**5 # This mask has index 4 (2^5, i.e. maks value = 32)
        
        
    comment = "- Extreme negatives (values below percentile " + str(np.round(percentile_min, 3)) + " = " + str(np.round(minimo, 3)) + " ) cleaned"    
    rss_out.log['extreme negative']['comment'] = comment
                

    """
    if plot:
        correction_map = g / self.intensity_corrected  # / g
        self.RSS_image(image=correction_map, cmap="binary_r", title=" - Correction map")
    """

    return rss_out








