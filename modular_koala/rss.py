# =============================================================================
# Basics packages
# =============================================================================
import numpy as np
import copy
import os
import matplotlib.pyplot as plt

# =============================================================================
# Astropy and associated packages
# =============================================================================
from astropy.io import fits
from astropy.nddata import bitfield_to_boolean_mask
# =============================================================================
# KOALA packages
# =============================================================================
# Modular
from modular_koala.ancillary import vprint




# =============================================================================
# Anciliary Functions
# =============================================================================

def coord_range(rss_list):
    RA = [rss.RA_centre_deg + rss.offset_RA_arcsec / 3600. for rss in rss_list]
    RA_min = np.nanmin(RA)
    RA_max = np.nanmax(RA)
    DEC = [rss.DEC_centre_deg + rss.offset_DEC_arcsec / 3600. for rss in rss_list]
    DEC_min = np.nanmin(DEC)
    DEC_max = np.nanmax(DEC)
    return RA_min, RA_max, DEC_min, DEC_max


def detect_edge(rss):
    """
    This functions detects the edges of a RSS. Returs the wavelenght and index
    of the first and last valid data. 

    Parameters
    ----------
    rss : TYPE
        DESCRIPTION.

    Returns
    -------
    left : TYPE
        DESCRIPTION.
    rigth : TYPE
        DESCRIPTION.

    """
    collapsed = np.sum(rss.intensity,0)
    nans = np.isfinite(collapsed)
    wavelength = rss.wavelength
    
    
    min_w = wavelength[nans].min()
    min_index = wavelength.tolist().index(min_w)

    max_w = wavelength[nans].max()
    max_index = wavelength.tolist().index(max_w)
 
    
    return min_w, min_index,  max_w, max_index







blank_log = {'read':{'comment':None,'index':None},
         'mask from file':{'comment':None,'index':0},
         'blue edge':{'comment':None,'index':1},
         'red edge':{'comment':None,'index':2},
         'cosmic':{'comment':None,'index':3},
         'extreme negative':{'comment':None,'index':4},
         'wavelenght fix':{'comment':None,'index':None,'sol':[]},
         }


blank_header = fits.header.Header(cards=[], copy=False)



# Example how to add header value at the end
# blank_header.append(('DARKCORR', 'OMIT', 'Dark Image Subtraction'), end=True)



# =============================================================================
# RSS CLASS
# =============================================================================


class RSS(object):
    """
    Collection of row-stacked spectra (RSS).

    Attributes
    ----------
    wavelength: np.array(float)
        Wavelength, in Angstrom
    intensity: np.array(float)
        Intensity :math:`I_\lambda` per unit wavelength.
        Axis 0 corresponds to fiber ID
        Axis 1 Corresponds to spectral dimension
    corrected_intensity: np.array(float)
        Intensity with all the corresponding corrections applied.
    variance: np.array(float)
        Variance :math:`\sigma^2_\lambda` per unit wavelength
        (note the square in the definition of the variance).
    corrected_variance: np.array(float)
        Variance with all the corresponding corrections applied.
    """

    # -----------------------------------------------------------------------------

    def __init__(self, 
                 intensity,
                 wavelength,
                 variance,
                 mask,
                 intensity_corrected,
                 variance_corrected,
                 log,
                 header,
                 spaxels_table,
                 ):
        
        self.intensity = intensity        
        self.wavelength = wavelength
        self.variance = variance
        self.mask = mask
        self.log = log
        self.intensity_corrected = intensity_corrected
        self.variance_corrected = variance_corrected
        self.header = header
        self.spaxels_table = spaxels_table


    # =============================================================================
    # Mask layer
    # =============================================================================
    def mask_layer(self,index=-1):
        """
        Identify the layers in the binary mask layer.
        
        """
        mask = self.mask.astype(int)
        mask_value = 2**(index)
        ignore_flags = [1,2,4,8,16,32]
        try:
            ignore_flags.remove(mask_value)
        except:
            print('Warning: '+str(index)+' is not a valid index. All layers considered')
        return bitfield_to_boolean_mask(mask, ignore_flags=ignore_flags)


    # =============================================================================
    # Imshow data
    # =============================================================================
    def show(self,pmin=5,pmax=95,mask=False,**kwargs):
        

        
        if mask:
            data = np.ma.masked_array(self.intensity_corrected,self.mask)
        else:
            data = self.intensity_corrected
        
        vmin,vmax = np.nanpercentile(data,(pmin,pmax))

        
        plt.imshow(data,vmin=vmin,vmax=vmax,**kwargs)

    # =============================================================================
    # show/save formated log        
    # =============================================================================
    def formated_log(self,verbose=True, save=None):
        
        pretty_line = '-------------------------------------------------------------------------------'
        
        import textwrap
        
        if verbose:
            
            for procedure in self.log.keys():

                comment = self.log[procedure]['comment']
                index = self.log[procedure]['index']
                applied = isinstance(comment, str) 
                
                mask_index = {None:'',
                              0 :'Mask index: 0 (bit mask value 2e0)',
                              1 :'Mask index: 1 (bit mask value 2e1)',
                              2 :'Mask index: 2 (bit mask value 2e2)',
                              3 :'Mask index: 3 (bit mask value 2e3)',
                              4 :'Mask index: 4 (bit mask value 2e4)',
                              }
                    
                if applied:
                    print('\n'+pretty_line)
                    
                    print ("{:<49}{:<2}".format(procedure.capitalize(), mask_index[index]))
                    print(pretty_line)
                    
                    for i in textwrap.wrap(comment+'\n', 80):
                        print (i)
                    print('\n')
    
        # if save:
        #     file = open(save,"w")
            
        #     file.write(pretty_line+'\n')
        #     file.write('{:<30} {:<5}'.format('Procedure', 'Implemented')+'\n')
        #     file.write(pretty_line+'\n')
            
        #     for procedure in self.log.keys():
        #         comment = self.log[procedure]
        #         applied = isinstance(comment, str)
                
        #         file.write('{:<30} {:<5}'.format(procedure.capitalize(), str(applied))+'\n')
                
        #         if applied:
        #             file.write('\n' + comment+'\n')
        #         file.write(pretty_line+'\n')
        #     file.close()

# =============================================================================
# Save RSS in fits     
# =============================================================================

    def tofits(self, name, 
               layer='corrected',
               overwrite=False,
               ):
        data = {'corrected': self.intensity_corrected,
                  'mask':self.mask,
                  }

        primary_hdu = fits.PrimaryHDU(data = data[layer])
        primary_hdu.header = self.header
        primary_hdu.verify('fix') 
        
        if self.spaxels_table is not None:

            # TODO: Why add again this information again in the table header?
            self.spaxels_table.header['CENRA']  =  self.header['RACEN']  / ( 180/np.pi )   # Must be in radians 
            self.spaxels_table.header['CENDEC']  =  self.header['DECCEN'] / ( 180/np.pi )
            
            hdu_list = fits.HDUList([primary_hdu,self.spaxels_table]) 
            hdu_list.writeto(name, overwrite=True)     

        else:                
            # Write the fits
            primary_hdu.writeto(name=name,overwrite=overwrite)
        
        
        

# =============================================================================
# Reading rss from file 
# =============================================================================

def read_rss(file_path,
             wcs,
             intentity_axis=0,
             variance_axis=None, 
             bad_spaxels_list=[],
             instrument=None,
             verbose=False,
             log = blank_log,
             header = blank_header,
             spaxels_table = None,
             ):

    file_name = os.path.basename(file_path)
    
    vprint.verbose = verbose
    
    vprint("\n> Reading RSS file", file_name,"created with",instrument,"...")

   

   
    #  Open fits file
    rss_fits = fits.open(file_path)
    
    # Read intensity using rss_fits_file[0]
    all_intensities = rss_fits[intentity_axis].data
    intensity = np.delete(all_intensities,bad_spaxels_list,0)
    
    # Bad pixel verbose sumary 
    vprint("\n  Number of spectra in this RSS =", len(all_intensities), 
            ",  number of good spectra =",len(intensity), 
            " ,  number of bad spectra =", len(bad_spaxels_list))
    
    if bad_spaxels_list != []: vprint("  Bad fibres =", bad_spaxels_list)



    # Read errors if exist a dedicated axis
    if variance_axis is not None:
        all_variances = rss_fits[variance_axis].data
        variance = np.delete(all_variances,bad_spaxels_list,0)

    else:
        vprint("\n  WARNING! Variance extension not found in fits file!")
        variance = copy.deepcopy(intensity)
    
    # Close fits file
    rss_fits.close()

    # Create wavelength from wcs
    nrow,ncol = wcs.array_shape
    wavelength_index =  np.arange(ncol)
    wavelength = wcs.dropaxis(1).wcs_pix2world(wavelength_index, 0)[0]
    
    # log    
    comment = ' '.join(['- RSS readed from ',file_name])
    log['read']['comment'] = comment
   
    # First Header value added by the PyKoala routine
    
    header.append(('DARKCORR', 'OMIT', 'Dark Image Subtraction'), end=True)

   
    
    # Blank mask (all 0, i.e. makning nothing) of the same shape of the data
    mask = np.zeros_like(intensity)
     
    # Blank corrected intensity (i.e. a copy of the data)
    intensity_corrected = copy.deepcopy(intensity)
    
    # Blank corrected variance (i.e a copy of the variacnce)
    variance_corrected = copy.deepcopy(variance)
    
    
    
    
    return RSS(intensity=intensity,
               wavelength=wavelength,
               variance=variance,
               mask=mask,
               intensity_corrected=intensity_corrected,
               variance_corrected=variance_corrected,
               log=log,
               header = header, 
               spaxels_table = spaxels_table,
               )



