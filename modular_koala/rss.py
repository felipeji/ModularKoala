# =============================================================================
# Basics packages
# =============================================================================
import numpy as np
import copy
import os
import matplotlib.pyplot as plt
from matplotlib import colors

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
    Detect the edges of a RSS. Returns the minimum and maximum wavelength that 
    determine the maximum interval with valid (i.e. no masked) data in all the 
    spaxels.
    
    
    Parameters
    ----------
    rss : RSS object.

    Returns
    -------
    min_w : float
        The lowest value (in units of the RSS wavelength) with 
        valid data in all spaxels.
    min_index : int
        Index of min_w in the RSS wavelength variable.
    max_w : float
        The higher value (in units of the RSS wavelength) with 
        valid data in all spaxels.
    max_index : int
        Index of max_w in the RSS wavelength variable.

    """
    collapsed = np.sum(rss.intensity,0)
    nans = np.isfinite(collapsed)
    wavelength = rss.wavelength
    
    
    min_w = wavelength[nans].min()
    min_index = wavelength.tolist().index(min_w)

    max_w = wavelength[nans].max()
    max_index = wavelength.tolist().index(max_w)
 
    
    return min_w, min_index,  max_w, max_index





# Blank dictionary for the log 
blank_log = {'read':{'comment':None,'index':None},
         'mask from file':{'comment':None,'index':0},
         'blue edge':{'comment':None,'index':1},
         'red edge':{'comment':None,'index':2},
         'cosmic':{'comment':None,'index':3},
         'extreme negative':{'comment':None,'index':4},
         'wavelenght fix':{'comment':None,'index':None,'sol':[]},
         }

# Blank Astropy Header object for the RSS header
# Example how to add header value at the end
# blank_header.append(('DARKCORR', 'OMIT', 'Dark Image Subtraction'), end=True)
blank_header = fits.header.Header(cards=[], copy=False)






# =============================================================================
# RSS CLASS
# =============================================================================


class RSS(object):
    """
    Container class for row-stacked spectra (RSS).

    Attributes
    ----------
    intensity: numpy.ndarray(float)
        Intensity :math:`I_\lambda` per unit wavelength.
        Axis 0 corresponds to fiber ID
        Axis 1 Corresponds to spectral dimension
    wavelength: numpy.ndarray(float)
        Wavelength, in Angstrom
    variance: numpy.ndarray(float)
        Variance :math:`\sigma^2_\lambda` per unit wavelength
        (note the square in the definition of the variance).
    mask : numpy.ndarray(float)
        Bit mask that records the pixels with individual corrections performed 
        by the various processes:
            Mask value      Correction
            -----------     ----------------
            1               Readed from file
            2               Blue edge
            4               Red edge
            8               NaNs mask
            16              Cosmic rays
            32              Extreme negative
            
    intensity_corrected: numpy.ndarray(float)
        Intensity with all the corresponding corrections applied (see log).
    variance_corrected: numpy.ndarray(float)
        Variance with all the corresponding corrections applied (see log).
    log : dict
        Dictionary containing a log of the processes applied on the rss.
        
    header : astropy.io.fits.header.Header object 
        The header associated with data.
    spaxels_table : astropy.io.fits.hdu.table.BinTableHDU object
        Bin table containing spaxel metadata.
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
    def mask_layer(self,index=-1,verbose=False):
        """
        Filter the bit mask according the index value:
            
	     Index    Layer                Bit mask value    
         ------   ----------------     --------------              
         0        Readed from file     1                           
         1        Blue edge            2                                 
         2        Red edge             4                               
         3        NaNs mask            8                                  
         4        Cosmic rays          16                                 
         5        Extreme negative     32                               


        Parameters
        ----------
        index : int or list, optional
            Layer index. The default -1 means that all the existing layers are 
            returned.

        Returns
        -------
        numpy.ndarray(bool)
            Boolean mask with the layer (or layers) selected.

        """
        vprint.verbose = verbose
        mask = self.mask.astype(int)
        if type(index) is not list: index = [index]
        ignore_flags = [1,2,4,8,16,32]
        
        for i in index:
            mask_value = 2**(i)
    
            try:
                ignore_flags.remove(mask_value)
            except:
                print('Warning: '+str(index)+' is not a valid index.')
            
            vprint('Layers considered: ',ignore_flags)
            
        return bitfield_to_boolean_mask(mask, ignore_flags=ignore_flags)
    
        
    def show_mask(self,):
        cmap = colors.ListedColormap(['darkgreen','blue','red','darkviolet','black','orange'])
        plt.imshow(np.log2(self.mask),vmin=0,vmax=6,  cmap=cmap, interpolation='none')
        #cb = plt.colorbar(figure,location="bottom",aspect=20)
        
        #ticks_position = np.arange(6)+0.5
        #ticks_name = ['From file','Blue edge','Red edge','NaN\'s','Cosmics','Extreme \n negatives']
        #cb.set_ticks(ticks_position)
        #cb.set_ticklabels(ticks_name)



    # =============================================================================
    # Imshow data
    # =============================================================================
    def show(self,pmin=5,pmax=95,mask=False,**kwargs):
        """
        Simple "imshow" of the corrected data (i.e. self.intensity_corrected ).
        Accept all (matplotlib.pyplot) imshow parameters.
        
        In the future we will implement the PyKoala RSS display function here. 
        
        Parameters
        ----------
        pmin : float, optional
            Minimum percentile of the data range that the colormap will covers. 
            The default is 5.
        pmax : TYPE, optional
            Maximum percentile of the data range that the colormap will covers. 
            The default is 95.
        mask : bool, optional
            True show the image with correceted pixels masked. 
            The default is False.
        """

        
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
                              0 :'Mask index: 0 (bit mask value 2**0)',
                              1 :'Mask index: 1 (bit mask value 2**1)',
                              2 :'Mask index: 2 (bit mask value 2**2)',
                              3 :'Mask index: 3 (bit mask value 2**3)',
                              4 :'Mask index: 4 (bit mask value 2**4)',
                              5 :'Mask index: 5 (bit mask value 2**5)',
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



