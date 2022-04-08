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
from modular_koala.rss import detect_edge

# Original
from koala.plot_plot import plot_plot

# =============================================================================

        # 6. ---------------------------------------------------            
        # SKY SUBSTRACTION      sky_method                      (S)
        #          
        # Several options here: (1) "1D"      : Consider a single sky spectrum, scale it and substract it
        #                       (2) "2D"      : Consider a 2D sky. i.e., a sky image, scale it and substract it fibre by fibre
        #                       (3) "self"    : Obtain the sky spectrum using the n_sky lowest fibres in the RSS file (DEFAULT)
        #                       (4) "none"    : No sky substraction is performed (DEFAULT)
        #                       (5) "1Dfit"   : Using an external 1D sky spectrum, fits sky lines in both sky spectrum AND all the fibres 
        #                       (6) "selffit" : Using the n_sky lowest fibres, obtain an sky spectrum, then fits sky lines in both sky spectrum AND all the fibres.



    # def read_sky_spectrum(self, sky_spectrum_file, path="", verbose = True):
    #     """
    #     Reads a TXT file with a 1D spectrum 

    #     Parameters
    #     ----------
    #     sky_spectrum_file : string (default = None)
    #         Specify the name of sky spectrum file (including or not the path)
    #     path: string (default = "")
    #         path to the sky spectrum file
    #     verbose : Boolean (optional)
    #         Print what is doing. The default is True.

    #     Returns
    #     -------
    #     sky_spectrum : array
    #         1D sky spectrum
            
    #     It also adds the 1D sky spectrum to self.sky_spectrum
    #     """
        
    #     if path != "" : sky_spectrum_file=full_path(sky_spectrum_file,path)
        
    #     if verbose:
    #         print("\n> Reading file with a 1D sky spectrum :")
    #         print(" ", sky_spectrum_file)

    #     w_sky, sky_spectrum = read_table(sky_spectrum_file, ["f", "f"])
        
    #     self.sky_spectrum = sky_spectrum

    #     self.history.append('- 1D sky spectrum provided in file :')
    #     self.history.append('  ' + sky_spectrum_file)

    #     if np.nanmedian(self.wavelength - w_sky) != 0:
    #         if verbose or warnings: print("\n\n  WARNING !!!! The wavelengths provided on the sky file do not match the wavelengths on this RSS !!\n\n")
    #         self.history.append('  WARNING: The wavelengths provided on the sky file do not match the wavelengths on this RSS')
    #     return sky_spectrum








def substract_sky(rss,
                  sky_model,
                  plot=True,
                  verbose=True,
                  warnings=True,
                  
                  # correct_negative_sky=False,
                  # order_fit_negative_sky=3,
                  # kernel_negative_sky=51,
                  # exclude_wlm=[[0, 0]],
                  # individual_check=True,
                  # use_fit_for_negative_sky=False,
                  # low_fibres=10,
                  ):
    """
    Substracts the sky_model to all fibres in the rss

    Parameters
    ----------
    correct_negative_sky : boolean (default = True)
        If True, and if the integrated value of the median sky is negative, this is corrected
    plot : boolean (default = True)
        Plots results
    see task 'correcting_negative_sky()' for definition of the rest of the parameters
    """

    # Set print verbose
    vprint.verbose = verbose
    
    # Copy input RSS for storage the changes implemented in the task   
    rss_out = copy.deepcopy(rss)

    # Substract sky in all fibers
    rss_out.intensity_corrected -= sky_model[np.newaxis, :]
    

    
    # TODO: Warning!!! THIS IS CORRECT? SHOULDN WE NEED TO COMPUTE VAR_{SKY}?
    # self.variance_corrected -= self.sky_emission[np.newaxis, :]



    

    if len(self.sky_fibres) > 0: last_sky_fibre = self.sky_fibres[-1]
    median_sky_corrected = np.zeros(self.n_spectra)

    for i in range(self.n_spectra):
        median_sky_corrected[i] = np.nanmedian(self.intensity_corrected[i, self.valid_wave_min_index:self.valid_wave_max_index], axis=0)
    
    if len(self.sky_fibres) > 0: median_sky_per_fibre = np.nanmedian(median_sky_corrected[self.sky_fibres])

    if verbose:
        print("  Median flux all fibres          = ", np.round(np.nanmedian(median_sky_corrected), 3))
        if len(self.sky_fibres) > 0:
            print("  Median flux sky fibres          = ", np.round(median_sky_per_fibre, 3))
            print("  Median flux brightest sky fibre = ", np.round(median_sky_corrected[last_sky_fibre], 3))
            print("  Median flux faintest  sky fibre = ", np.round(median_sky_corrected[self.sky_fibres[0]], 3))



    # Plot median value of fibre vs. fibre
    # if plot:

    #     if len(self.sky_fibres) > 0:
    #         ymin = median_sky_corrected[self.sky_fibres[0]] - 1
    #         # ymax = np.nanpercentile(median_sky_corrected,90),
    #         hlines = [np.nanmedian(median_sky_corrected), median_sky_corrected[self.sky_fibres[0]],
    #                   median_sky_corrected[last_sky_fibre], median_sky_per_fibre]
    #         chlines = ["r", "k", "k", "g"]
    #         ptitle = "Median flux per fibre after sky substraction\n (red = median flux all fibres, green = median flux sky fibres, grey = median flux faintest/brightest sky fibre)"
    #     else:
    #         ymin = np.nanpercentile(median_sky_corrected, 1)
    #         # ymax=np.nanpercentile(self.sky_emission, 1)
    #         hlines = [np.nanmedian(median_sky_corrected), 0]
    #         chlines = ["r", "k"]
    #         ptitle = "Median flux per fibre after sky substraction (red = median flux all fibres)"

    #     plot_plot(list(range(self.n_spectra)), median_sky_corrected,
    #               ylabel="Median Flux [counts]", xlabel="Fibre",
    #               ymin=ymin, ymax=np.nanpercentile(median_sky_corrected, 90),
    #               hlines=hlines, chlines=chlines,
    #               ptitle=ptitle)

    if len(self.sky_fibres) > 0:
        if median_sky_corrected[self.sky_fibres[0]] < 0:
            if verbose or warnings: print(
                "  WARNING !  The integrated value of the sky fibre with the smallest value is negative!")

            # if correct_negative_sky:
            #     if verbose: print("  Fixing this, as 'correct_negative_sky' = True  ... ")
            #     self.correcting_negative_sky(plot=plot, low_fibres=low_fibres, exclude_wlm=exclude_wlm,
            #                                   kernel_negative_sky=kernel_negative_sky,
            #                                   use_fit_for_negative_sky=use_fit_for_negative_sky,
            #                                   order_fit_negative_sky=order_fit_negative_sky,
            #                                   individual_check=individual_check)

    vprint("  Intensities corrected for sky emission and stored in self.intensity_corrected !")
    #history.append("  Intensities corrected for the sky emission")

    
    return rss_out










    # # =============================================================================
    # def correcting_negative_sky(self, low_fibres=10, kernel_negative_sky=51, order_fit_negative_sky=3, edgelow=0,
    #                             clip_fit_negative_sky = 0.8,
    #                             edgehigh=0,  # step=11, weight_fit_median = 1, scale = 1.0,
    #                             use_fit_for_negative_sky=False, individual_check=True, force_sky_fibres_to_zero=True,
    #                             exclude_wlm=[[0, 0]],
    #                             show_fibres=[0, 450, 985], fig_size=12, plot=True, verbose=True):
    #     """
    #     Corrects negative sky with a median spectrum of the lowest intensity fibres

    #     Parameters
    #     ----------
    #     low_fibres : integer (default = 10)
    #         amount of fibres allocated to act as fibres with the lowest intensity
    #     kernel_negative_sky : odd integer (default = 51)
    #         kernel parameter for smooth median spectrum
    #     order_fit_negative_sky : integer (default = 3)
    #         order of polynomial used for smoothening and fitting the spectrum
    #     edgelow, edgehigh : integers (default = 0, 0)
    #         Minimum and maximum pixel number such that any pixel in between this range is only to be considered
    #     use_fit_for_negative_sky: boolean (default = False)
    #         Substract the order-order fit instead of the smoothed median spectrum
    #     individual_check: boolean (default = True)
    #         Check individual fibres and correct if integrated value is negative
    #     exclude_wlm : list
    #         exclusion command to prevent large absorption lines from being affected by negative sky correction : (lower wavelength, upper wavelength)
    #     show_fibres : list of integers (default = [0,450,985])
    #         List of fibres to show
    #     force_sky_fibres_to_zero : boolean (default = True)
    #         If True, fibres defined at the sky will get an integrated value = 0
    #     fig_size: float (default = 12)
    #         Size of the figure
    #     plot : boolean (default = False)
    #        Plot figure
    #     """

    #     # CHECK fit_smooth_spectrum and compare with medfilt
    #     w = self.wavelength
    #     # Set limits
    #     if edgelow == 0: edgelow = self.valid_wave_min_index
    #     if edgehigh == 0: edgehigh = np.int((self.n_wave - self.valid_wave_max_index) / 2)

    #     plot_this = False
    #     if len(show_fibres) > 0:
    #         show_fibres.append(self.integrated_fibre_sorted[-1])  # Adding the brightest fibre
    #         show_fibres.append(self.integrated_fibre_sorted[0])  # Adding the faintest fibre

    #     if individual_check:
    #         if verbose: print("\n> Individual correction of fibres with negative sky ... ")
    #         if force_sky_fibres_to_zero and verbose: print("  Also forcing integrated spectrum of sky_fibres = 0 ... ")
    #         corrected_not_sky_fibres = 0
    #         total_corrected = 0
    #         sky_fibres_to_zero = 0
    #         for fibre in range(self.n_spectra):
    #             corregir = False
    #             if fibre in show_fibres and plot:
    #                 print("\n - Checking fibre", fibre, "...")
    #                 plot_this = True
    #             else:
    #                 plot_this = False
    #             smooth, fit = fit_smooth_spectrum(w, self.intensity_corrected[fibre], 
    #                                               mask = [self.mask[0][fibre],self.mask[1][fibre]],
    #                                               edgelow=edgelow, edgehigh=edgehigh, #remove_nans=False,
    #                                               kernel_fit=kernel_negative_sky, 
    #                                               index_fit=order_fit_negative_sky, clip_fit = clip_fit_negative_sky,
    #                                               plot=plot_this, verbose=False, hlines=[0.], ptitle="",
    #                                               fcal=False)
    #             if np.nanpercentile(fit, 5) < 0:
    #                 if fibre not in self.sky_fibres: corrected_not_sky_fibres = corrected_not_sky_fibres + 1
    #                 corregir = True
    #             else:
    #                 if fibre in self.sky_fibres and force_sky_fibres_to_zero:
    #                     corregir == True
    #                     sky_fibres_to_zero = sky_fibres_to_zero + 1

    #             if corregir == True:
    #                 total_corrected = total_corrected + 1
    #                 if use_fit_for_negative_sky:
    #                     if fibre in show_fibres and verbose and plot: print(
    #                         "      Using fit to smooth spectrum for correcting the negative sky in fibre", fibre, " ...")
    #                     self.intensity_corrected[fibre] -= fit
    #                     # self.variance_corrected[fibre] -= fit
    #                 else:
    #                     if fibre in show_fibres and verbose and plot: print(
    #                         "      Using smooth spectrum for correcting the negative sky in fibre", fibre, " ...")
    #                     self.intensity_corrected[fibre] -= smooth
    #                     # self.variance_corrected[fibre] -= smooth
    #             else:
    #                 if fibre in show_fibres and verbose and plot: print("      Fibre", fibre,
    #                                                            "does not need to be corrected for negative sky ...")

    #         corrected_sky_fibres = total_corrected - corrected_not_sky_fibres
    #         if verbose:
    #             print("\n> Corrected {} fibres (not defined as sky) and {} out of {} sky fibres !".format(
    #                 corrected_not_sky_fibres, corrected_sky_fibres, len(self.sky_fibres)))
    #             if force_sky_fibres_to_zero:
    #                 print("  The integrated spectrum of", sky_fibres_to_zero, "sky fibres have been forced to 0.")
    #                 print("  The integrated spectrum of all sky_fibres have been set to 0.")
    #         self.history.append("- Individual correction of negative sky applied")
    #         self.history.append("  Corrected " + np.str(corrected_not_sky_fibres) + " not-sky fibres")
    #         if force_sky_fibres_to_zero:
    #             self.history.append("  All the " + np.str(len(self.sky_fibres)) + " sky fibres have been set to 0")
    #         else:
    #             self.history.append("  Corrected " + np.str(corrected_sky_fibres) + " out of " + np.str(
    #                 len(self.sky_fibres)) + " sky fibres")

    #     else:
    #         # Get integrated spectrum of n_low lowest fibres and use this for ALL FIBRES
    #         integrated_intensity_sorted = np.argsort(self.integrated_fibre)
    #         region = integrated_intensity_sorted[0:low_fibres]
    #         Ic = np.nanmedian(self.intensity_corrected[region], axis=0)

    #         if verbose:
    #             print("\n> Correcting negative sky using median spectrum combining the", low_fibres,
    #                   "fibres with the lowest integrated intensity")
    #             print("  which are :", region)
    #             print("  Obtaining smoothed spectrum using a {} kernel and fitting a {} order polynomium...".format(
    #                 kernel_negative_sky, order_fit_negative_sky))
    #         ptitle = self.object + " - " + str(low_fibres) + " fibres with lowest intensity - Fitting an order " + str(
    #             order_fit_negative_sky) + " polynomium to spectrum smoothed with a " + str(
    #             kernel_negative_sky) + " kernel window"
    #         smooth, fit = fit_smooth_spectrum(self.wavelength, Ic, kernel=kernel_negative_sky, edgelow=edgelow,
    #                                           edgehigh=edgehigh, verbose=False, #mask=self.mask[],
    #                                           order=order_fit_negative_sky, plot=plot, hlines=[0.], ptitle=ptitle,
    #                                           fcal=False)
    #         if use_fit_for_negative_sky:
    #             self.smooth_negative_sky = fit
    #             if verbose: print(
    #                 "  Sustracting fit to smoothed spectrum of {} low intensity fibres to all fibres ...".format(
    #                     low_fibres))
    #         else:
    #             self.smooth_negative_sky = smooth
    #             if verbose: print(
    #                 "  Sustracting smoothed spectrum of {} low intensity fibres to all fibres ...".format(low_fibres))

    #         for i in range(self.n_spectra):
    #             self.intensity_corrected[i, :] = self.intensity_corrected[i, :] - self.smooth_negative_sky
    #             # self.sky_emission = self.sky_emission - self.smooth_negative_sky

    #         # TODO: New implementation including variance
    #         # self.intensity_corrected -= self.smooth_negative_sky[np.newaxis, :]
    #         # self.variance_corrected -= self.smooth_negative_sky[np.newaxis, :]

    #         if verbose: print("  This smoothed spectrum is stored in self.smooth_negative_sky")
    #         self.history.append("- Correcting negative sky using smoothed spectrum of the")
    #         self.history.append("  " + np.str(low_fibres) + " fibres with the lowest integrated value")
    #     # -----------------------------------------------------------------------------









    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # def fit_and_substract_sky_spectrum(self, sky, w=1000, spectra=1000, rebin=False,
    #                                    # If rebin == True, it fits all wavelengths to be at the same wavelengths that SKY spectrum...
    #                                    brightest_line="Ha", brightest_line_wavelength=0,
    #                                    maxima_sigma=3.0, ymin=-50, ymax=600, wmin=0, wmax=0,
    #                                    auto_scale_sky=False,
    #                                    sky_lines_file="",
    #                                    warnings=False, verbose=False, plot=False, plot_step_fibres=True, step=100,
    #                                    fig_size=12, fibre=-1, max_flux_variation=15.,
    #                                    min_flux_ratio=-1, max_flux_ratio=-1):
    #     """
    #     Given a 1D sky spectrum, this task fits
    #     sky lines of each spectrum individually and substracts sky
    #     Needs the observed wavelength (brightest_line_wavelength) of the brightest emission line (brightest_line) .
    #     w is the wavelength
    #     spec the 2D spectra
    #     max_flux_variation = 15. is the % of the maximum value variation for the flux OBJ/SKY
    #                          A value of 15 will restrict fits to 0.85 < OBJ/SKY < 1.15
    #                          Similarly, we can use:  min_flux_ratio <  OBJ/SKY < max_flux_ratio
    #     Parameters
    #     ----------
    #     sky : list of floats
    #         Given sky spectrum
    #     w : list of floats  (default = 1000 )
    #         If given, sets this vector as wavelength
    #         It has to have the same len than sky.
    #     spectra : integer
    #         2D spectra
    #     rebin :  boolean
    #         wavelengths are fitted to the same wavelengths of the sky spectrum if True
    #     brightest_line : string (default "Ha")
    #         string name with the emission line that is expected to be the brightest in integrated spectrum
    #     brightest_line_wavelength : integer
    #         wavelength that corresponds to the brightest line
    #     maxima_sigma : float
    #         ## Sets maximum removal of FWHM for data
    #     ymin : integer
    #        sets the bottom edge of the bounding box (plotting)
    #     ymax : integer
    #        sets the top edge of the bounding box (plotting)
    #     wmin : integer
    #         sets the lower wavelength range (minimum), if 0, no range is set, specified wavelength (w) is investigated
    #     wmax : integer
    #         sets the upper wavelength range (maximum), if 0, no range is set, specified wavelength (w) is investigated
    #     auto_scale_sky : boolean
    #         scales sky spectrum for subtraction if True
    #     sky_lines_file :
    #         file containing list of sky lines to fit
    #     warnings : boolean
    #         disables warnings if set to False
    #     verbose : boolean (default = True)
    #         Print results
    #     plot : boolean (default = False)
    #         Plot results
    #     plot_step_fibres : boolean
    #        if True, plots ever odd fibre spectrum...
    #     step : integer (default = 50)
    #        step using for estimating the local medium value
    #     fig_size:
    #        Size of the figure (in x-axis), default: fig_size=10
    #     fibre: integer (default 0)
    #        If fibre is given, it identifies emission lines in the given fibre
    #     """
    #     if min_flux_ratio == -1: min_flux_ratio = 1. - max_flux_variation / 100.
    #     if max_flux_ratio == -1: max_flux_ratio = 1. + max_flux_variation / 100.

    #     template_path_prefix = './input_data/sky_lines/'
    #     # TODO: (FUTURE WORK) Include the possibility of using other files provided by the user
    #     if not sky_lines_file:
    #         sky_lines_file = template_path_prefix + "sky_lines_bright.dat"
    #         print(' > Using "sky_lines_bright.dat" as sky line template')
    #     if sky_lines_file == "ALL":
    #         sky_lines_file = template_path_prefix + "sky_lines.dat"
    #     if sky_lines_file == "BRIGHT":
    #         sky_lines_file = template_path_prefix + "sky_lines_bright.dat"
    #     if sky_lines_file == "IR":
    #         sky_lines_file = template_path_prefix + "sky_lines_IR.dat"
    #     if sky_lines_file in ["IRshort", "IRs", "IR_short"]:
    #         sky_lines_file = template_path_prefix + "sky_lines_IR_short.dat"

    #     self.history.append('  Skylines fitted following file:')
    #     self.history.append('  ' + sky_lines_file)

    #     print("\n> Fitting selected sky lines to both sky spectrum and object spectra ...\n")

    #     # TODO: It is easier to use dictionaries (also it is possible to
    #     # include a longer dictionary at constants module)
    #     common_bright_lines = {'Ha': 6562.82, 'O3b': 5006.84, 'Hb': 4861.33}
    #     if brightest_line in list(common_bright_lines.keys()):
    #         brightest_line_wavelength_rest = common_bright_lines[brightest_line]

    #     if brightest_line_wavelength != 0:
    #         print(
    #             "  - Using {} at rest wavelength {:6.2f} identified by the user at {:6.2f} to avoid fitting emission lines...".format(
    #                 brightest_line, brightest_line_wavelength_rest, brightest_line_wavelength))
    #     else:
    #         print(
    #             "  - No wavelength provided to 'brightest_line_wavelength', the object is NOT expected to have emission lines\n")

    #     redshift = brightest_line_wavelength / brightest_line_wavelength_rest - 1.

    #     if w == 1000:
    #         w = self.wavelength
    #     if spectra == 1000:
    #         spectra = copy.deepcopy(self.intensity_corrected)

    #     if wmin == 0: wmin = w[0]
    #     if wmax == 0: wmax = w[-1]

    #     print("  - Reading file with the list of sky lines to fit :")
    #     print("   ", sky_lines_file)

    #     # Read file with sky emission lines
    #     sl_center_, sl_name_, sl_fnl_, sl_lowlow_, sl_lowhigh_, sl_highlow_, sl_highhigh_, sl_lmin_, sl_lmax_ = read_table(
    #         sky_lines_file, ["f", "s", "f", "f", "f", "f", "f", "f", "f"])
    #     # number_sl = len(sl_center)

    #     # MOST IMPORTANT EMISSION LINES IN RED
    #     # 6300.30       [OI]  -0.263   30.0 15.0   20.0   40.0
    #     # 6312.10     [SIII]  -0.264   30.0 18.0    5.0   20.0
    #     # 6363.78       [OI]  -0.271   20.0  4.0    5.0   30.0
    #     # 6548.03      [NII]  -0.296   45.0 15.0   55.0   75.0
    #     # 6562.82         Ha  -0.298   50.0 25.0   35.0   60.0
    #     # 6583.41      [NII]  -0.300   62.0 42.0    7.0   35.0
    #     # 6678.15        HeI  -0.313   20.0  6.0    6.0   20.0
    #     # 6716.47      [SII]  -0.318   40.0 15.0   22.0   45.0
    #     # 6730.85      [SII]  -0.320   50.0 30.0    7.0   35.0
    #     # 7065.28        HeI  -0.364   30.0  7.0    7.0   30.0
    #     # 7135.78    [ArIII]  -0.374   25.0  6.0    6.0   25.0
    #     # 7318.39      [OII]  -0.398   30.0  6.0   20.0   45.0
    #     # 7329.66      [OII]  -0.400   40.0 16.0   10.0   35.0
    #     # 7751.10    [ArIII]  -0.455   30.0 15.0   15.0   30.0
    #     # 9068.90    [S-III]  -0.594   30.0 15.0   15.0   30.0

    #     el_list_no_z = [6300.3, 6312.10, 6363.78, 6548.03, 6562.82, 6583.41, 6678.15, 6716.47, 6730.85, 7065.28,
    #                     7135.78, 7318.39, 7329.66, 7751.1, 9068.9]
    #     el_list = (redshift + 1) * np.array(el_list_no_z)
    #     #  [OI]   [SIII]  [OI]   Ha+[NII]  HeI    [SII]     HeI   [ArIII]  [OII]  [ArIII]  [SIII]
    #     el_low_list_no_z = [6296.3, 6308.1, 6359.8, 6544.0, 6674.2, 6712.5, 7061.3, 7129., 7312., 7747.1, 9063.9]
    #     el_high_list_no_z = [6304.3, 6316.1, 6367.8, 6590.0, 6682.2, 6736.9, 7069.3, 7141., 7336., 7755.1, 9073.9]
    #     el_low_list = (redshift + 1) * np.array(el_low_list_no_z)
    #     el_high_list = (redshift + 1) * np.array(el_high_list_no_z)

    #     # Double Skylines
    #     dsky1_ = [6257.82, 6465.34, 6828.22, 6969.70, 7239.41, 7295.81, 7711.50, 7750.56, 7853.391, 7913.57, 7773.00,
    #               7870.05, 8280.94, 8344.613, 9152.2, 9092.7, 9216.5, 8827.112, 8761.2, 0]  # 8760.6, 0]#
    #     dsky2_ = [6265.50, 6470.91, 6832.70, 6978.45, 7244.43, 7303.92, 7715.50, 7759.89, 7860.662, 7921.02, 7780.43,
    #               7879.96, 8288.34, 8352.78, 9160.9, 9102.8, 9224.8, 8836.27, 8767.7, 0]  # 8767.2, 0] #

    #     # Be sure the lines we are using are in the requested wavelength range
    #     # print "  Checking the values of skylines in the file", sky_lines_file
    #     # for i in range(len(sl_center_)):
    #     #    print sl_center_[i],sl_fnl_[i],sl_lowlow_[i],sl_lowhigh_[i],sl_highlow_[i],sl_highhigh_[i],sl_lmin_[i],sl_lmax_[i]
    #     # print "  We only need skylines in the {} - {} range:".format(self.valid_wave_min, self.valid_wave_max)
    #     print("  - We only need sky lines in the {} - {} range ".format(np.round(self.wavelength[0], 2),
    #                                                                     np.round(self.wavelength[-1], 2)))

    #     # valid_skylines = np.where((sl_center_ < self.valid_wave_max) & (sl_center_ > self.valid_wave_min))

    #     valid_skylines = np.where((sl_center_ < self.wavelength[-1]) & (sl_center_ > self.wavelength[0]))

    #     sl_center = sl_center_[valid_skylines]
    #     sl_fnl = sl_fnl_[valid_skylines]
    #     sl_lowlow = sl_lowlow_[valid_skylines]
    #     sl_lowhigh = sl_lowhigh_[valid_skylines]
    #     sl_highlow = sl_highlow_[valid_skylines]
    #     sl_highhigh = sl_highhigh_[valid_skylines]
    #     sl_lmin = sl_lmin_[valid_skylines]
    #     sl_lmax = sl_lmax_[valid_skylines]
    #     number_sl = len(sl_center)

    #     dsky1 = []
    #     dsky2 = []
    #     for l in range(number_sl):
    #         if sl_center[l] in dsky1_:
    #             dsky1.append(dsky1_[dsky1_.index(sl_center[l])])
    #             dsky2.append(dsky2_[dsky1_.index(sl_center[l])])

    #     print("  - All sky lines: ", sl_center)
    #     print("  - Double sky lines: ", dsky1)
    #     print("  - Total number of skylines to fit =", len(sl_center))
    #     print("  - Valid values for OBJ / SKY Gauss ratio  = ( ", min_flux_ratio, ",", max_flux_ratio, ")")
    #     print("  - Maxima sigma to consider a valid fit  = ", maxima_sigma, " A\n")

    #     say_status = 0
    #     self.wavelength_offset_per_fibre = []
    #     self.sky_auto_scale = []
    #     f_new_ALL = []
    #     sky_sl_gaussian_fitted_ALL = []
    #     only_fibre = False
    #     if fibre != -1:
    #         f_i = fibre
    #         f_f = fibre + 1
    #         print("\n ----> Checking fibre ", fibre, " (only this fibre is corrected, use fibre = -1 for all)...")
    #         plot = True
    #         verbose = True
    #         warnings = True
    #         only_fibre = True
    #         say_status = fibre
    #     else:
    #         f_i = 0
    #         f_f = self.n_spectra

    #     # Check if skylines are located within the range of an emission line !
    #     skip_sl_fit = [False] * number_sl
    #     if verbose or fibre == -1: print("  - Checking skylines within emission line ranges...")
    #     for i in range(number_sl):
    #         for j in range(len(el_low_list)):
    #             if el_low_list[j] < sl_center[i] < el_high_list[j]:
    #                 skip_sl_fit[i] = True
    #                 if verbose or fibre == -1: print('  ------> SKY line', sl_center[i], 'in EMISSION LINE !  ',
    #                                                  el_low_list[j], sl_center[i], el_high_list[j])

    #                 # Gaussian fits to the sky spectrum
    #     sl_gaussian_flux = []
    #     sl_gaussian_sigma = []
    #     sl_gauss_center = []
    #     sky_sl_gaussian_fitted = copy.deepcopy(sky)
    #     if verbose or fibre == -1: print("  - Performing Gaussian fitting to sky lines in sky spectrum...")
    #     for i in range(number_sl):
    #         if sl_fnl[i] == 0:
    #             plot_fit = False
    #         else:
    #             plot_fit = True
    #         if sl_center[i] in dsky1:  # == dsky1[di] :
    #             if fibre == -1: print("    DOUBLE IN SKY: ", sl_center[i], dsky2[dsky1.index(sl_center[i])])
    #             warnings_ = False
    #             if sl_fnl[i] == 1:
    #                 warnings_ = True
    #                 if verbose: print("    Line ", sl_center[i], " blended with ", dsky2[dsky1.index(sl_center[i])])
    #             resultado = dfluxes(w, sky_sl_gaussian_fitted, sl_center[i], dsky2[dsky1.index(sl_center[i])],
    #                                 lowlow=sl_lowlow[i], lowhigh=sl_lowhigh[i], highlow=sl_highlow[i],
    #                                 highhigh=sl_highhigh[i], lmin=sl_lmin[i], lmax=sl_lmax[i], fmin=-20, fmax=0,
    #                                 broad1=2.1 * 2.355, broad2=2.1 * 2.355, plot=plot_fit, verbose=False,
    #                                 plot_sus=False,
    #                                 fcal=False, warnings=warnings_)  # Broad is FWHM for Gaussian sigm a= 1,

    #             sl_gaussian_flux.append(resultado[3])  # 15 is Gauss 1, 16 is Gauss 2, 3 is Total Gauss
    #             sl_gauss_center.append(resultado[1])
    #             sl_gaussian_sigma.append(resultado[5] / 2.355)
    #             # sl_gaussian_flux.append(resultado[16])
    #             # sl_gauss_center.append(resultado[12])
    #             # sl_gaussian_sigma.append(resultado[14]/2.355)
    #             # 12     13      14        15              16
    #             # fit[3], fit[4],fit[5], gaussian_flux_1, gaussian_flux_2 # KANAN

    #         else:
    #             resultado = fluxes(w, sky_sl_gaussian_fitted, sl_center[i], lowlow=sl_lowlow[i], lowhigh=sl_lowhigh[i],
    #                                highlow=sl_highlow[i], highhigh=sl_highhigh[i], lmin=sl_lmin[i], lmax=sl_lmax[i],
    #                                fmin=-20, fmax=0,
    #                                broad=2.1 * 2.355, plot=plot_fit, verbose=False, plot_sus=False, fcal=False,
    #                                warnings=warnings)  # Broad is FWHM for Gaussian sigm a= 1,

    #             sl_gaussian_flux.append(resultado[3])
    #             sl_gauss_center.append(resultado[1])
    #             sl_gaussian_sigma.append(resultado[5] / 2.355)

    #         if plot_fit:
    #             if verbose:  print("    Fitted wavelength for sky line ", sl_center[i], " : ", sl_gauss_center[-1],
    #                                "  sigma = ", sl_gaussian_sigma[-1])
    #             wmin = sl_lmin[i]
    #             wmax = sl_lmax[i]

    #         if skip_sl_fit[i] == False:
    #             sky_sl_gaussian_fitted = resultado[11]
    #         else:
    #             if verbose: print('  ------> SKY line', sl_center[i], 'in EMISSION LINE !')

    #     # Now Gaussian fits to fibres
    #     for fibre in range(f_i, f_f):  # (self.n_spectra):
    #         if fibre == say_status:
    #             if fibre == 0: print(" ")
    #             print("  - Checking fibre {:4} ...  ({:6.2f} % completed) ...".format(fibre,
    #                                                                                   fibre * 100. / self.n_spectra))
    #             say_status = say_status + step
    #             if plot_step_fibres: plot = True
    #         else:
    #             plot = False

    #         skip_el_fit = copy.deepcopy(skip_sl_fit)

    #         # Gaussian fit to object spectrum                       #BOBA
    #         object_sl_gaussian_flux = []
    #         object_sl_gaussian_sigma = []
    #         ratio_object_sky_sl_gaussian = []
    #         dif_center_obj_sky = []
    #         spec = spectra[fibre]
    #         object_sl_gaussian_fitted = copy.deepcopy(spec)
    #         object_sl_gaussian_center = []
    #         if verbose: print("\n  - Performing Gaussian fitting to sky lines in fibre", fibre, "of object data ...")

    #         for i in range(number_sl):
    #             if sl_fnl[i] == 0:
    #                 plot_fit = False
    #             else:
    #                 plot_fit = True
    #             if skip_el_fit[i]:
    #                 if verbose: print("    SKIPPING SKY LINE", sl_center[i],
    #                                   "as located within the range of an emission line!")
    #                 object_sl_gaussian_flux.append(float('nan'))  # The value of the SKY SPECTRUM
    #                 object_sl_gaussian_center.append(float('nan'))
    #                 object_sl_gaussian_sigma.append(float('nan'))
    #                 dif_center_obj_sky.append(float('nan'))
    #             else:
    #                 if sl_center[i] in dsky1:  # == dsky1[di] :
    #                     if fibre == -1: print("    DOUBLE IN SKY: ", sl_center[i], dsky2[dsky1.index(sl_center[i])])
    #                     warnings_ = False
    #                     if sl_fnl[i] == 1:
    #                         if fibre == -1:
    #                             warnings_ = True
    #                         if verbose: print("    Line ", sl_center[i], " blended with ",
    #                                           dsky2[dsky1.index(sl_center[i])])
    #                     resultado = dfluxes(w, object_sl_gaussian_fitted, sl_center[i],
    #                                         dsky2[dsky1.index(sl_center[i])], lowlow=sl_lowlow[i],
    #                                         lowhigh=sl_lowhigh[i], highlow=sl_highlow[i], highhigh=sl_highhigh[i],
    #                                         lmin=sl_lmin[i], lmax=sl_lmax[i], fmin=-20, fmax=0,
    #                                         broad1=sl_gaussian_sigma[i] * 2.355, broad2=sl_gaussian_sigma[i] * 2.355,
    #                                         plot=plot_fit, verbose=False, plot_sus=False, fcal=False,
    #                                         warnings=warnings_)
    #                     if verbose:
    #                         print(
    #                             "    line = {:.3f} : center = {:.3f}, gauss = {:.2f},  sigma = {:.2f}, flux = {:.2f}".format(
    #                                 sl_center[i], resultado[1], sl_gaussian_sigma[i], resultado[5] / 2.355,
    #                                 resultado[15]))
    #                         print(
    #                             "    line = {:.3f} : center = {:.3f}, gauss = {:.2f},  sigma = {:.2f}, flux = {:.2f}".format(
    #                                 dsky2[dsky1.index(sl_center[i])], resultado[12], sl_gaussian_sigma[i],
    #                                 resultado[14] / 2.355, resultado[16]))
    #                         print("    For skylines ", sl_center[i], "+", dsky2[dsky1.index(sl_center[i])],
    #                               " the total flux is ", np.round(sl_gaussian_flux[i], 3),
    #                               ",                     OBJ/SKY = ", np.round(resultado[3] / sl_gaussian_flux[i], 3))

    #                     if resultado[3] > 0 and resultado[5] / 2.355 < maxima_sigma and resultado[15] > 0 and resultado[
    #                         14] / 2.355 < maxima_sigma and resultado[3] / sl_gaussian_flux[i] > min_flux_ratio and \
    #                             resultado[3] / sl_gaussian_flux[
    #                         i] < max_flux_ratio * 1.25:  # and resultado[5] < maxima_sigma: # -100000.: #0:
    #                         object_sl_gaussian_fitted = resultado[11]

    #                         object_sl_gaussian_flux.append(resultado[15])
    #                         object_sl_gaussian_center.append(resultado[1])
    #                         object_sl_gaussian_sigma.append(resultado[5] / 2.355)
    #                         # object_sl_gaussian_flux.append(resultado[16])
    #                         # object_sl_gaussian_center.append(resultado[12])
    #                         # object_sl_gaussian_sigma.append(resultado[14]/2.355)

    #                         dif_center_obj_sky.append(object_sl_gaussian_center[i] - sl_gauss_center[i])
    #                     else:
    #                         if verbose: print("    Bad double fit for ", sl_center[i], "! trying single fit...")
    #                         average_wave = (sl_center[i] + dsky2[dsky1.index(sl_center[i])]) / 2
    #                         resultado = fluxes(w, object_sl_gaussian_fitted, average_wave, lowlow=sl_lowlow[i],
    #                                            lowhigh=sl_lowhigh[i], highlow=sl_highlow[i], highhigh=sl_highhigh[i],
    #                                            lmin=average_wave - 50, lmax=average_wave + 50, fmin=-20, fmax=0,
    #                                            broad=4.5, plot=plot_fit, verbose=False, plot_sus=False, fcal=False,
    #                                            warnings=warnings)  # Broad is FWHM for Gaussian sigma= 1,
    #                         if verbose: print(
    #                             "    line = {:.3f} : center = {:.3f}, gauss = {:.2f},  sigma = {:.2f}, flux = {:.2f},   OBJ/SKY = {:.3f}".format(
    #                                 sl_center[i], resultado[1], sl_gaussian_sigma[i], resultado[5] / 2.355,
    #                                 resultado[3], resultado[3] / sl_gaussian_flux[i]))
    #                         if resultado[3] > 0 and resultado[5] / 2.355 < maxima_sigma * 2. and resultado[3] / \
    #                                 sl_gaussian_flux[i] > min_flux_ratio and resultado[3] / sl_gaussian_flux[
    #                             i] < max_flux_ratio:  # and resultado[5] < maxima_sigma: # -100000.: #0:
    #                             object_sl_gaussian_flux.append(resultado[3])
    #                             object_sl_gaussian_fitted = resultado[11]
    #                             object_sl_gaussian_center.append(resultado[1])
    #                             object_sl_gaussian_sigma.append(resultado[5] / 2.355)
    #                             dif_center_obj_sky.append(object_sl_gaussian_center[i] - sl_gauss_center[i])
    #                         else:
    #                             if verbose: print("    -> Bad fit for ", sl_center[i], "! ignoring it...")
    #                             object_sl_gaussian_flux.append(float('nan'))
    #                             object_sl_gaussian_center.append(float('nan'))
    #                             object_sl_gaussian_sigma.append(float('nan'))
    #                             dif_center_obj_sky.append(float('nan'))
    #                             skip_el_fit[i] = True  # We don't substract this fit
    #                 else:
    #                     resultado = fluxes(w, object_sl_gaussian_fitted, sl_center[i], lowlow=sl_lowlow[i],
    #                                        lowhigh=sl_lowhigh[i], highlow=sl_highlow[i], highhigh=sl_highhigh[i],
    #                                        lmin=sl_lmin[i], lmax=sl_lmax[i], fmin=0, fmax=0,
    #                                        broad=sl_gaussian_sigma[i] * 2.355, plot=plot_fit, verbose=False,
    #                                        plot_sus=False, fcal=False,
    #                                        warnings=warnings)  # Broad is FWHM for Gaussian sigma= 1,
    #                     if verbose: print(
    #                         "    line = {:.3f} : center = {:.3f}, gauss = {:.2f},  sigma = {:.2f}, flux = {:.2f},   OBJ/SKY = {:.3f}".format(
    #                             sl_center[i], resultado[1], sl_gaussian_sigma[i], resultado[5] / 2.355, resultado[3],
    #                             resultado[3] / sl_gaussian_flux[i]))
    #                     if resultado[3] > 0 and resultado[5] / 2.355 < maxima_sigma and resultado[3] / sl_gaussian_flux[
    #                         i] > min_flux_ratio and resultado[3] / sl_gaussian_flux[
    #                         i] < max_flux_ratio:  # and resultado[5] < maxima_sigma: # -100000.: #0:
    #                         object_sl_gaussian_flux.append(resultado[3])
    #                         object_sl_gaussian_fitted = resultado[11]
    #                         object_sl_gaussian_center.append(resultado[1])
    #                         object_sl_gaussian_sigma.append(resultado[5] / 2.355)
    #                         dif_center_obj_sky.append(object_sl_gaussian_center[i] - sl_gauss_center[i])
    #                     else:
    #                         if verbose: print("    -> Bad fit for ", sl_center[i], "! ignoring it...")
    #                         object_sl_gaussian_flux.append(float('nan'))
    #                         object_sl_gaussian_center.append(float('nan'))
    #                         object_sl_gaussian_sigma.append(float('nan'))
    #                         dif_center_obj_sky.append(float('nan'))
    #                         skip_el_fit[i] = True  # We don't substract this fit

    #             try:
    #                 ratio_object_sky_sl_gaussian.append(object_sl_gaussian_flux[i] / sl_gaussian_flux[i])
    #             except Exception:
    #                 print("\n\n\n\n\n DIVISION FAILED in ", sl_center[i], "!!!!!   sl_gaussian_flux[i] = ",
    #                       sl_gaussian_flux[i], "\n\n\n\n")
    #                 ratio_object_sky_sl_gaussian.append(1.)

    #         # Scale sky lines that are located in emission lines or provided negative values in fit
    #         # reference_sl = 1 # Position in the file! Position 1 is sky line 6363.4
    #         # sl_ref_ratio = sl_gaussian_flux/sl_gaussian_flux[reference_sl]
    #         if verbose:
    #             print("  - Correcting skylines for which we couldn't get a Gaussian fit and are not in an emission line range...")
    #         for i in range(number_sl):
    #             if skip_el_fit[i] == True and skip_sl_fit[i] == False:  # Only those that are NOT in emission lines
    #                 # Use known center, sigma of the sky and peak
    #                 gauss_fix = sl_gaussian_sigma[i]
    #                 small_center_correction = 0.
    #                 # Check if center of previous sky line has a small difference in wavelength
    #                 small_center_correction = np.nanmedian(dif_center_obj_sky[0:i])
    #                 if verbose:
    #                     print("  - Small correction of center wavelength of sky line ", sl_center[i], "  :",
    #                           small_center_correction)

    #                 object_sl_gaussian_fitted = substract_given_gaussian(w, object_sl_gaussian_fitted,
    #                                                                      sl_center[i] + small_center_correction, peak=0,
    #                                                                      sigma=gauss_fix, flux=0, search_peak=True,
    #                                                                      lowlow=sl_lowlow[i], lowhigh=sl_lowhigh[i],
    #                                                                      highlow=sl_highlow[i], highhigh=sl_highhigh[i],
    #                                                                      lmin=sl_lmin[i], lmax=sl_lmax[i], plot=False,
    #                                                                      verbose=verbose)

    #                 # Substract second Gaussian if needed !!!!!
    #                 for di in range(len(dsky1)):
    #                     if sl_center[i] == dsky1[di]:
    #                         if verbose: print("    This was a double sky line, also substracting ",
    #                                           dsky2[dsky1.index(sl_center[i])], "  at ", np.round(
    #                                 np.array(dsky2[dsky1.index(sl_center[i])]) + small_center_correction, 2))
    #                         object_sl_gaussian_fitted = substract_given_gaussian(w, object_sl_gaussian_fitted, np.array(
    #                             dsky2[dsky1.index(sl_center[i])]) + small_center_correction, peak=0, sigma=gauss_fix,
    #                                                                              flux=0, search_peak=True,
    #                                                                              lowlow=sl_lowlow[i],
    #                                                                              lowhigh=sl_lowhigh[i],
    #                                                                              highlow=sl_highlow[i],
    #                                                                              highhigh=sl_highhigh[i],
    #                                                                              lmin=sl_lmin[i], lmax=sl_lmax[i],
    #                                                                              plot=False, verbose=verbose)
    #             else:
    #                 if skip_sl_fit[i] == True and skip_el_fit[i] == True:
    #                     if verbose: print("     - SKIPPING SKY LINE", sl_center[i],
    #                                       " as located within the range of an emission line!")

    #         offset = np.nanmedian(np.array(object_sl_gaussian_center) - np.array(sl_gauss_center))
    #         offset_std = np.nanstd(np.array(object_sl_gaussian_center) - np.array(sl_gauss_center))

    #         good_ratio_values = []
    #         for ratio in ratio_object_sky_sl_gaussian:
    #             if np.isnan(ratio) == False:
    #                 if ratio > min_flux_ratio and ratio < max_flux_ratio:
    #                     good_ratio_values.append(ratio)

    #         valid_median_flux = np.nanmedian(good_ratio_values)

    #         if verbose:
    #             print("  - Median center offset between OBJ and SKY :", np.round(offset, 3), " A ,    std = ",
    #                   np.round(offset_std, 3))
    #             print("    Median gauss for the OBJECT              :",
    #                   np.round(np.nanmedian(object_sl_gaussian_sigma), 3), " A ,    std = ",
    #                   np.round(np.nanstd(object_sl_gaussian_sigma), 3))
    #             print("    Median flux OBJECT / SKY                 :",
    #                   np.round(np.nanmedian(ratio_object_sky_sl_gaussian), 3), "   ,    std = ",
    #                   np.round(np.nanstd(ratio_object_sky_sl_gaussian), 3))
    #             print("    Median flux OBJECT / SKY VALID VALUES    :", np.round(valid_median_flux, 3),
    #                   "   ,    std = ", np.round(np.nanstd(good_ratio_values), 3))
    #             print("  - min and max flux OBJECT / SKY = ", np.round(np.nanmin(ratio_object_sky_sl_gaussian), 3), ",",
    #                   np.round(np.nanmax(ratio_object_sky_sl_gaussian), 3), "  -> That is a variation of ",
    #                   np.round(-100. * (np.nanmin(ratio_object_sky_sl_gaussian) - 1), 2), "% and ",
    #                   np.round(100. * (np.nanmax(ratio_object_sky_sl_gaussian) - 1), 2), "%")
    #             print("                                                        but only fits with < ",
    #                   max_flux_variation, "% have been considered")
    #         if plot == True and only_fibre == True:
    #             # for i in range(len(sl_gauss_center)):
    #             #    print i+1, sl_gauss_center[i],ratio_object_sky_sl_gaussian[i]
    #             plt.figure(figsize=(12, 5))
    #             plt.plot(sl_gauss_center, ratio_object_sky_sl_gaussian, "+", ms=12, mew=2)
    #             plt.axhline(y=np.nanmedian(ratio_object_sky_sl_gaussian), color="k", linestyle='--', alpha=0.3)
    #             plt.axhline(y=valid_median_flux, color="g", linestyle='-', alpha=0.3)
    #             plt.axhline(y=valid_median_flux + np.nanstd(good_ratio_values), color="c", linestyle=':', alpha=0.5)
    #             plt.axhline(y=valid_median_flux - np.nanstd(good_ratio_values), color="c", linestyle=':', alpha=0.5)
    #             plt.axhline(y=min_flux_ratio, color="r", linestyle='-', alpha=0.5)
    #             plt.axhline(y=max_flux_ratio, color="r", linestyle='-', alpha=0.5)
    #             # plt.ylim(0.7,1.3)
    #             ptitle = "Checking flux OBJECT / SKY for fitted skylines in fibre " + np.str(
    #                 fibre)  # +" with rms = "+np.str(rms[i])
    #             plt.title(ptitle)
    #             plt.xlabel("Wavelength [$\mathrm{\AA}$]")
    #             plt.ylabel("OBJECT / SKY ")
    #             # plt.legend(frameon=True, loc=2, ncol=6)
    #             plt.minorticks_on()
    #             plt.show()
    #             plt.close()

    #         self.wavelength_offset_per_fibre.append(offset)
    #         # self.sky_auto_scale.append(np.nanmedian(ratio_object_sky_sl_gaussian))
    #         self.sky_auto_scale.append(valid_median_flux)

    #         if auto_scale_sky:
    #             if verbose:  print("  - As requested, using this value to scale sky spectrum before substraction... ")
    #             auto_scale = np.nanmedian(ratio_object_sky_sl_gaussian)
    #         else:
    #             if verbose:  print(
    #                 "  - As requested, DO NOT using this value to scale sky spectrum before substraction... ")
    #             auto_scale = 1.0
    #         if rebin:
    #             if verbose:
    #                 print("\n> Rebinning the spectrum of fibre", fibre, "to match sky spectrum...")
    #             f = object_sl_gaussian_fitted
    #             f_new = rebin_spec_shift(w, f, offset)
    #         else:
    #             f_new = object_sl_gaussian_fitted

    #         # This must be corrected at then end to use the median auto_scale value
    #         # self.intensity_corrected[fibre] = f_new - auto_scale * sky_sl_gaussian_fitted
    #         f_new_ALL.append(f_new)
    #         sky_sl_gaussian_fitted_ALL.append(sky_sl_gaussian_fitted)

    #         if plot:
    #             plt.figure(figsize=(12, 5))
    #             plt.plot(w, spec, "purple", alpha=0.7, label="Obj")
    #             plt.plot(w, auto_scale * sky, "r", alpha=0.5, label="Scaled sky")
    #             plt.plot(w, auto_scale * sky_sl_gaussian_fitted, "lime", alpha=0.8, label="Scaled sky fit")
    #             plt.plot(w, object_sl_gaussian_fitted, "k", alpha=0.5, label="Obj - sky fit")
    #             plt.plot(w, spec - auto_scale * sky, "orange", alpha=0.4, label="Obj - scaled sky")
    #             plt.plot(w, object_sl_gaussian_fitted - sky_sl_gaussian_fitted, "b", alpha=0.9,
    #                      label="Obj - sky fit - scale * rest sky")

    #             plt.xlim(wmin, wmax)
    #             plt.ylim(ymin, ymax)
    #             ptitle = "Fibre " + np.str(fibre)  # +" with rms = "+np.str(rms[i])
    #             plt.title(ptitle)
    #             plt.xlabel("Wavelength [$\mathrm{\AA}$]")
    #             plt.ylabel("Flux [counts]")
    #             plt.legend(frameon=True, loc=2, ncol=6)
    #             plt.minorticks_on()
    #             for i in range(len(el_list)):
    #                 plt.axvline(x=el_list[i], color="k", linestyle='-', alpha=0.5)  # MARIO
    #             for i in range(number_sl):
    #                 if skip_sl_fit[i]:
    #                     alpha = 0.1
    #                 else:
    #                     alpha = 0.6
    #                 if sl_fnl[i] == 1:
    #                     plt.axvline(x=sl_center[i], color="brown", linestyle='-', alpha=alpha + 0.4)  # alpha=1)
    #                 else:
    #                     plt.axvline(x=sl_center[i], color="y", linestyle='--', alpha=alpha)
    #             for i in range(len(dsky2) - 1):
    #                 plt.axvline(x=dsky2[i], color="orange", linestyle='--', alpha=0.6)
    #             plt.show()
    #             plt.close()

    #         if only_fibre:
    #             ymax = np.nanpercentile(self.intensity_corrected[fibre], 99.5)
    #             ymin = np.nanpercentile(self.intensity_corrected[fibre], 0.1) - (
    #                         np.nanpercentile(self.intensity_corrected[fibre], 99.5) - np.nanpercentile(
    #                     self.intensity_corrected[fibre], 0.1)) / 15.
    #             self.intensity_corrected[fibre] = f_new - auto_scale * sky_sl_gaussian_fitted
    #             plot_plot(w, [self.intensity_corrected[fibre], self.intensity[fibre]], color=["b", "r"], ymin=ymin,
    #                       ymax=ymax,
    #                       ptitle="Comparison before (red) and after (blue) sky substraction using Gaussian fit to skylines")
    #             print("\n  Only fibre", fibre, " is corrected, use fibre = -1 for all...")

    #     if only_fibre == False:
    #         # To avoid bad auto scaling with bright fibres or weird fibres,
    #         # we fit a 2nd order polynomium to a filtered median value
    #         sas_m = medfilt(self.sky_auto_scale, 21)  ## Assuming odd_number = 21
    #         # fit=np.polyfit(range(self.n_spectra),sas_m,2)   # If everything is OK this should NOT be a fit, but a median
    #         fit = np.nanmedian(sas_m)
    #         # y=np.poly1d(fit)
    #         fity = [fit] * self.n_spectra
    #         # fity=y(range(self.n_spectra))
    #         if plot_step_fibres:
    #             # ptitle = "Fit to autoscale values:\n"+np.str(y)
    #             ptitle = "Checking autoscale values, median value = " + np.str(
    #                 np.round(fit, 2)) + " using median filter 21"
    #             ymin_ = np.nanmin(sas_m) - 0.1
    #             ymax_ = np.nanmax(sas_m) + 0.4
    #             plot_plot(list(range(self.n_spectra)), [sas_m, fity, self.sky_auto_scale, ],
    #                       color=["b", "g", "r"], alpha=[0.5, 0.5, 0.8], ptitle=ptitle, ymin=ymin_, ymax=ymax_,
    #                       xlabel="Fibre", ylabel="Flux ratio", label=["autoscale medfilt=21", "median", "autoscale"])
    #             # label=["autoscale med=21", "fit","autoscale"])

    #         self.sky_auto_scale_fit = fity
    #         if auto_scale_sky:
    #             # print "  Correcting all fluxes adding the autoscale value of the FIT above for each fibre..."
    #             print("  Correcting all fluxes adding the median autoscale value to each fibre (green line)...")
    #         else:
    #             # print "  Correcting all fluxes WITHOUT CONSIDERING the autoscale value of the FIT above for each fibre..."
    #             print("  Correcting all fluxes WITHOUT CONSIDERING the median autoscale value ...")

    #         for fibre in range(self.n_spectra):
    #             if auto_scale_sky:
    #                 self.intensity_corrected[fibre] = f_new_ALL[fibre] - self.sky_auto_scale_fit[fibre] * \
    #                                                   sky_sl_gaussian_fitted_ALL[fibre]
    #             else:
    #                 self.intensity_corrected[fibre] = f_new_ALL[fibre] - sky_sl_gaussian_fitted_ALL[fibre]
    #         print("\n  All fibres corrected for sky emission performing individual Gaussian fits to each fibre !")
    #         self.history.append(
    #             "  Intensities corrected for the sky emission performing individual Gaussian fits to each fibre")
    #         # self.variance_corrected += self.sky_variance # TODO: Check if telluric/ext corrections were applied before
    #         # DOES VARIANCE NEED SKY SUBSTRACTION?
    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
   
    
   
    
   
    
   
    # def find_sky_fibres(self, sky_wave_min=0, sky_wave_max=0, n_sky=200, plot=False, verbose=True, warnings=True):
    #     """
    #     Identify n_sky spaxels with the LOWEST INTEGRATED VALUES and store them in self.sky_fibres

    #     Parameters
    #     ----------
    #     sky_wave_min, sky_wave_max : float, float (default 0, 0)
    #         Consider the integrated flux in the range [sky_wave_min, sky_wave_max]
    #         If 0, they are set to self.valid_wave_min or self.valid_wave_max
    #     n_sky : integer (default = 200)
    #         number of spaxels used for identifying sky.
    #         200 is a good number for calibration stars
    #         for real objects, particularly extense objects, set n_sky = 30 - 50
    #     plot : boolean (default = False)
    #         plots a RSS map with sky positions
    #     """
    #     if sky_wave_min == 0: sky_wave_min = self.valid_wave_min
    #     if sky_wave_max == 0: sky_wave_max = self.valid_wave_max
    #     # Assuming cleaning of cosmics and CCD defects, we just use the spaxels with the LOWEST INTEGRATED VALUES
    #     self.compute_integrated_fibre(valid_wave_min=sky_wave_min, valid_wave_max=sky_wave_max, plot=False,
    #                                   verbose=verbose, warnings=warnings)
    #     sorted_by_flux = np.argsort(self.integrated_fibre)
    #     print("\n> Identifying sky spaxels using the lowest integrated values in the [", np.round(sky_wave_min, 2), ",",
    #           np.round(sky_wave_max, 2), "] range ...")
    #     print("  We use the lowest", n_sky, "fibres for getting sky. Their positions are:")
    #     # Compute sky spectrum and plot RSS map with sky positions if requested
    #     self.sky_fibres = sorted_by_flux[:n_sky]
    #     if plot: self.RSS_map(self.integrated_fibre, None, self.sky_fibres, title=" - Sky Spaxels")










    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # def find_sky_emission(self, intensidad=[0, 0],  n_sky=200,
    #                       sky_fibres=[], sky_wave_min=0, sky_wave_max=0,
    #                       win_sky=0, include_history=True,
    #                       log=True, gamma=0, plot=True):
    #     """
    #     Find the sky emission given fibre list or taking n_sky fibres with lowest integrated value.

    #     Parameters
    #     ----------
    #     intensidad :
    #         Matrix with intensities (self.intensity or self.intensity_corrected)
    #     plot : boolean (default = True)
    #         Plots results
    #     n_sky : integer (default = 200)
    #         number of lowest intensity fibres that will be used to create a SKY specrum
    #         200 is a good number for calibration stars
    #         for real objects, particularly extense objects, set n_sky = 30 - 50
    #     sky_fibres : list of floats (default = [1000])
    #        fibre or fibre range associated with sky
    #        If [1000], then the fibre list of sky will be computed automatically using n_sky
    #     sky_wave_min, sky_wave_max : float, float (default 0, 0)
    #         Only used when sky_fibres is [1000]
    #         Consider the integrated flux in the range [sky_wave_min, sky_wave_max]
    #         If 0, they are set to self.valid_wave_min or self.valid_wave_max
    #     log, gamma:
    #         Normalization scale, default is lineal scale.
    #         Lineal scale: norm=colors.Normalize().   log = False, gamma = 0
    #         Log scale:    norm=colors.LogNorm()      log = True, gamma = 0
    #         Power law:    norm=colors.PowerNorm(gamma=1./4.) when gamma != 0
    #         //
    #     win_sky : odd integer (default = 0)
    #         Width in fibres of a median filter applied to obtain sky spectrum
    #         If 0, it will not apply any median filter.
    #     include_history : boolean (default = True)
    #         If True, it includes RSS.history the basic information
    #     """
    #     if len(sky_fibres) == 0:
    #         if sky_wave_min == 0: sky_wave_min = self.valid_wave_min
    #         if sky_wave_max == 0: sky_wave_max = self.valid_wave_max
    #         self.find_sky_fibres(sky_wave_min=sky_wave_min, sky_wave_max=sky_wave_max, n_sky=n_sky)
    #     else:  # We provide a list with sky positions
    #         print("  We use the list provided to get the sky spectrum")
    #         print("  sky_fibres = ", sky_fibres)
    #         self.sky_fibres = np.array(sky_fibres)

    #     if plot: self.RSS_map(self.integrated_fibre, list_spectra=self.sky_fibres, log=log, gamma=gamma, title=" - Sky Spaxels")
    #     print("  List of fibres used for sky saved in self.sky_fibres")

    #     if include_history: self.history.append("- Obtaining the sky emission using " + np.str(n_sky) + " fibres")
    #     self.sky_emission = sky_spectrum_from_fibres(self, self.sky_fibres, win_sky=win_sky, plot=False,
    #                                                  include_history=include_history)

    #     if plot: plot_plot(self.wavelength, self.sky_emission, color="c",
    #                        ylabel="Relative flux [counts]", xlabel="Wavelength [$\mathrm{\AA}$]",
    #                        xmin=self.wavelength[0] - 10, xmax=self.wavelength[-1] + 10,
    #                        ymin=np.nanpercentile(self.sky_emission, 1), ymax=np.nanpercentile(self.sky_emission, 99),
    #                        vlines=[self.valid_wave_min, self.valid_wave_max],
    #                        ptitle="Combined sky spectrum using the requested fibres")
    #     print("  Sky spectrum obtained and stored in self.sky_emission !! ")

    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------





    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # def apply_self_sky(self, sky_fibres=[], sky_spectrum=[],  sky_spectrum_file="", path="",   
    #                    sky_wave_min=0, sky_wave_max=0, win_sky=0, scale_sky_1D=0,
    #                    brightest_line="Ha", brightest_line_wavelength=0, ranges_with_emission_lines=[0],
    #                    cut_red_end=0, low_fibres=10, use_fit_for_negative_sky=False, kernel_negative_sky=51,
    #                    order_fit_negative_sky=3, n_sky=50, correct_negative_sky=False,
    #                    individual_check=False, verbose=True, plot=True):
    #     """

    #     Apply sky correction using the specified number of lowest fibres in the RSS file to obtain the sky spectrum

    #     Parameters
    #     ----------
    #     sky_fibres : list of integers (default = none)
    #         Specify the fibres to use to obtain sky spectrum. Will automatically determine the best fibres if not specified
    #     sky_spectrum : list of floats (default = none)
    #         Specify the sky spectrum to be used for correction. If not specified, will derive it automatically
    #     sky_spectrum_file : string (default = None)
    #         Specify the name of sky spectrum file (including or not the path)
    #     path: string (default = "")
    #         path to the sky spectrum file    
    #     plot : boolean (default = True)
    #         Show the plots in the console
    #     sky_wave_min : float (default = 0)
    #         Specify the lower bound on wavelength range. If 0, it is set to self.valid_wave_min
    #     sky_wave_max : float (default = 0)
    #         Specify the upper bound on wavelength range. If 0, it is set to self.valid_wave_max
    #     win_sky : odd integer (default = 0)
    #         Width in fibres of a median filter applied to obtain sky spectrum, if 0, it will not apply any median filter
    #     scale_sky_1D : float (default = 0)
    #         Specify the scale between the sky emission and the object, if 0, will find it automatically
    #     brightest_line : string (default = "Ha")
    #         Specify the brightest emission line in the object spectrum, by default it is H-alpha
    #         Options: “O3”: [OIII] 5007, “O3b”: [OIII] 4959, “Ha”: H-alpha 6563, “Hb”: H-beta 4861.
    #     brightest_line_wavelength : float (default = 0)
    #         Wavelength of the brightest emission line, if 0, will take a stored value for emission line specified
    #     ranges_with_emission_lines = list of floats (default = [0])
    #         Specify ranges containing emission lines than needs to be corrected
    #     cut_red_end : float (default = 0)
    #         Apply mask to the red end of the spectrum. If 0, will proceed, if -1, will do nothing
    #     low_fibres : integer (default = 10)
    #         amount of fibres allocated to act as fibres with the lowest intensity
    #     use_fit_for_negative_sky: boolean (default = False)
    #         Substract the order-order fit instead of the smoothed median spectrum
    #     kernel_negative_sky : odd integer (default = 51)
    #         kernel parameter for smooth median spectrum
    #     order_fit_negative_sky : integer (default = 3)
    #         order of polynomial used for smoothening and fitting the spectrum
    #     verbose : boolean (default = True)
    #         Print detailed description of steps taken in console
    #     n_sky : integer (default = 50)
    #         Number of fibres to use for finding sky spectrum
    #     correct_negative_sky : boolean (default = True)
    #         If True, and if the integrated value of the median sky is negative, this is corrected
    #     individual_check: boolean (default = True)
    #         Check individual fibres and correct if integrated value is negative
    #     """

    #     self.history.append('- Sky sustraction using the self method')

    #     if len(sky_fibres) != 0:
    #         n_sky = len(sky_fibres)
    #         print("\n> 'sky_method = self', using list of", n_sky, "fibres to create a sky spectrum ...")
    #         self.history.append('  A list of ' + np.str(n_sky) + ' fibres was provided to create the sky spectrum')
    #         self.history.append(np.str(sky_fibres))
    #     else:
    #         print("\n> 'sky_method = self', hence using", n_sky, "lowest intensity fibres to create a sky spectrum ...")
    #         self.history.append(
    #             '  The ' + np.str(n_sky) + ' lowest intensity fibres were used to create the sky spectrum')
            
    #     if sky_spectrum_file != "":            
    #         sky_spectrum = self.read_sky_spectrum(sky_spectrum_file, path=path, verbose = verbose)

    #     if len(sky_spectrum) == 0:
    #         self.find_sky_emission(n_sky=n_sky, plot=plot, sky_fibres=sky_fibres,
    #                                sky_wave_min=sky_wave_min, sky_wave_max=sky_wave_max,
    #                                win_sky=win_sky, include_history=True)

    #     else:
    #         print("  Sky spectrum provided. Using this for replacing regions with bright emission lines...")

    #         self.find_sky_emission(n_sky=n_sky, plot=plot, sky_fibres=sky_fibres,
    #                                sky_wave_min=sky_wave_min, sky_wave_max=sky_wave_max,
    #                                win_sky=win_sky, include_history=False)

    #         sky_r_self = self.sky_emission

    #         self.sky_emission = replace_el_in_sky_spectrum(self, sky_r_self, sky_spectrum, scale_sky_1D=scale_sky_1D,
    #                                                        brightest_line=brightest_line,
    #                                                        brightest_line_wavelength=brightest_line_wavelength,
    #                                                        ranges_with_emission_lines=ranges_with_emission_lines,
    #                                                        cut_red_end=cut_red_end,
    #                                                        plot=plot)
    #         self.history.append('  Using sky spectrum provided for replacing regions with emission lines')

    #     self.substract_sky(plot=plot, low_fibres=low_fibres,
    #                        correct_negative_sky=correct_negative_sky, use_fit_for_negative_sky=use_fit_for_negative_sky,
    #                        kernel_negative_sky=kernel_negative_sky, order_fit_negative_sky=order_fit_negative_sky,
    #                        individual_check=individual_check)

    #     self.apply_mask(verbose=verbose, make_nans=True)
    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # def apply_1D_sky(self, sky_fibres=[], sky_spectrum=[], sky_wave_min=0, sky_wave_max=0,
    #                  win_sky=0, include_history=True,
    #                  scale_sky_1D=0, remove_5577=True, sky_spectrum_file="", path="",
    #                  plot=True, verbose=True, n_sky=50):
    #     """

    #     Apply sky correction using 1D spectrum provided

    #     Parameters
    #     ----------
    #     sky_fibres : list of integers (default = none)
    #         Specify the fibres to use to obtain sky spectrum. Will automatically determine the best fibres if not specified
    #     sky_spectrum : list of floats (default = none)
    #         Specify the sky spectrum to be used for correction. If not specified, will derive it automatically
    #     sky_wave_min : float (default = 0)
    #         Specify the lower bound on wavelength range. If 0, it is set to self.valid_wave_min
    #     sky_wave_max : float (default = 0)
    #         Specify the upper bound on wavelength range. If 0, it is set to self.valid_wave_max
    #     win_sky : odd integer (default = 0)
    #         Width in fibres of a median filter applied to obtain sky spectrum, if 0, it will not apply any median filter
    #     scale_sky_1D : float (default = 0)
    #         Specify the scale between the sky emission and the object, if 0, will find it automatically
    #     include_history : boolean (default = True)
    #         Include the task completion into the RSS object history
    #     remove_5577 : boolean (default = True)
    #         Remove the line 5577 from the data
    #     sky_spectrum_file : string (default = None)
    #         Specify the name of sky spectrum file (including or not the path)
    #     path: string (default = "")
    #         path to the sky spectrum file
    #     plot : boolean (default = True)
    #         Show the plots in the console
    #     verbose : boolean (default = True)
    #         Print detailed description of steps taken in console
    #     n_sky : integer (default = 50)
    #         Number of fibres to use for finding sky spectrum
    #     """

    #     self.history.append('- Sky sustraction using the 1D method')
        
    #     if sky_spectrum_file != "":            
    #         sky_spectrum = self.read_sky_spectrum(sky_spectrum_file, path=path, verbose = verbose)
  
    #     if verbose:
    #         print("\n> Sustracting the sky using the sky spectrum provided, checking the scale OBJ/SKY...")
    #     if scale_sky_1D == 0:
    #         if verbose:
    #             print("  No scale between 1D sky spectrum and object given, calculating...")

    #         # TODO !
    #         # Task "scale_sky_spectrum" uses sky lines, needs to be checked...
    #         # self.sky_emission,scale_sky_1D_auto=scale_sky_spectrum(self.wavelength, sky_spectrum, self.intensity_corrected,
    #         #                                     cut_sky=cut_sky, fmax=fmax, fmin=fmin, fibre_list=fibre_list)

    #         # Find self sky emission using only the lowest n_sky fibres (this should be small, 20-25)
    #         if n_sky == 50: n_sky = 20
    #         self.find_sky_emission(n_sky=n_sky, plot=plot, sky_fibres=sky_fibres,
    #                                sky_wave_min=sky_wave_min, sky_wave_max=sky_wave_max,
    #                                win_sky=win_sky, include_history=include_history)

    #         sky_r_self = self.sky_emission

    #         scale_sky_1D = auto_scale_two_spectra(self, sky_r_self, sky_spectrum, scale=[0.1, 1.01, 0.025],
    #                                               w_scale_min=self.valid_wave_min, w_scale_max=self.valid_wave_max,
    #                                               plot=plot, verbose=True)

    #     elif verbose:
    #         print("  As requested, we scale the given 1D sky spectrum by", scale_sky_1D)

    #     self.sky_emission = sky_spectrum * scale_sky_1D
    #     self.history.append('  1D sky spectrum scaled by =' + np.str(scale_sky_1D))

    #     if verbose: print("\n> Scaled sky spectrum stored in self.sky_emission, substracting to all fibres...")

    #     # For blue spectra, remove 5577 in the sky spectrum...
    #     if self.valid_wave_min < 5577 and remove_5577 == True:
    #         if verbose: print("  Removing sky line 5577.34 from the sky spectrum...")
    #         resultado = fluxes(self.wavelength, self.sky_emission, 5577.34, lowlow=30, lowhigh=10, highlow=10,
    #                            highhigh=30,
    #                            plot=False, verbose=False)  # fmin=-5.0E-17, fmax=2.0E-16,
    #         # resultado = [rms_cont, fit[0], fit_error[0], gaussian_flux, gaussian_flux_error, fwhm, fwhm_error, flux, flux_error, ew, ew_error, spectrum  ]
    #         self.sky_emission = resultado[11]
    #     else:
    #         if self.valid_wave_min < 5577 and verbose: print(
    #             "  Sky line 5577.34 is not removed from the sky spectrum...")

    #     # Remove 5577 in the object
    #     if self.valid_wave_min < 5577 and remove_5577 == True and scale_sky_1D == 0:  # and individual_sky_substraction == False:
    #         if verbose:
    #             print("  Removing sky line 5577.34 from the object...")
    #         self.history.append("  Sky line 5577.34 removed performing Gaussian fit")

    #         wlm = self.wavelength
    #         for i in range(self.n_spectra):
    #             s = self.intensity_corrected[i]
    #             # Removing Skyline 5577 using Gaussian fit if requested
    #             resultado = fluxes(wlm, s, 5577.34, lowlow=30, lowhigh=10, highlow=10, highhigh=30,
    #                                plot=False, verbose=False)  # fmin=-5.0E-17, fmax=2.0E-16,
    #             # resultado = [rms_cont, fit[0], fit_error[0], gaussian_flux, gaussian_flux_error, fwhm, fwhm_error, flux, flux_error, ew, ew_error, spectrum  ]
    #             self.intensity_corrected[i] = resultado[11]
    #     else:
    #         if self.valid_wave_min < 5577 and verbose:
    #             if scale_sky_1D == 0:
    #                 print("  Sky line 5577.34 is not removed from the object...")
    #             else:
    #                 print("  Sky line 5577.34 already removed in object during CCD cleaning...")

    #     self.substract_sky(plot=plot, verbose=verbose)

    #     if plot:
    #         text = "Sky spectrum (scaled using a factor " + np.str(scale_sky_1D) + " )"
    #         plot_plot(self.wavelength, self.sky_emission, hlines=[0], ptitle=text,
    #                   xmin=self.wavelength[0] - 10, xmax=self.wavelength[-1] + 10, color="c",
    #                   vlines=[self.valid_wave_min, self.valid_wave_max])
    #     if verbose:
    #         print("  Intensities corrected for sky emission and stored in self.intensity_corrected !")
    #     self.sky_emission = sky_spectrum  # Restore sky_emission to original sky_spectrum
    #     # self.apply_mask(verbose=verbose)














    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # def apply_1Dfit_sky(self, sky_spectrum=[], n_sky=50, sky_fibres=[], sky_spectrum_file="", path="",
    #                     sky_wave_min=0, sky_wave_max=0, win_sky=0, scale_sky_1D=0,
    #                     sky_lines_file="", brightest_line_wavelength=0,
    #                     brightest_line="Ha", maxima_sigma=3, auto_scale_sky=False,
    #                     plot=True, verbose=True, fig_size=12, fibre_p=-1, kernel_correct_ccd_defects=51):
    #     """
    #     Apply 1Dfit sky correction.

    #     Parameters
    #     ----------
    #     sky_spectrum : array or list of floats (default = none)
    #         Specify the sky spectrum to be used for correction. If not specified, will derive it automatically
    #     n_sky : integer (default = 50)
    #         Number of fibres to use for finding sky spectrum
    #     sky_fibres : list of integers (default = none)
    #         Specify the fibres to use to obtain sky spectrum. Will automatically determine the best fibres if not specified
    #     sky_spectrum_file : string (default = None)
    #         Specify the name of sky spectrum file (including or not the path)
    #     path: string (default = "")
    #         path to the sky spectrum file
    #     sky_wave_min : float (default = 0)
    #         Specify the lower bound on wavelength range. If 0, it is set to self.valid_wave_min
    #     sky_wave_max : float (default = 0)
    #         Specify the upper bound on wavelength range. If 0, it is set to self.valid_wave_max
    #     win_sky : odd integer (default = 0)
    #         Width in fibres of a median filter applied to obtain sky spectrum, if 0, it will not apply any median filter
    #     scale_sky_1D : float (default = 0)
    #         Specify the scale between the sky emission and the object, if 0, will find it automatically
    #     sky_lines_file : string (default = None)
    #         Specify the path and name of sky lines file
    #     brightest_line_wavelength : float (default = 0)
    #         Wavelength of the brightest emission line, if 0, will take a stored value for emission line specified
    #     brightest_line : string (default = "Ha")
    #         Specify the brightest emission line in the object spectrum, by default it is H-alpha
    #         Options: “O3”: [OIII] 5007, “O3b”: [OIII] 4959, “Ha”: H-alpha 6563, “Hb”: H-beta 4861.
    #     maxima_sigma : float (default = 3)
    #         Maximum allowed standard deviation for Gaussian fit
    #     auto_scale_sky : boolean (default = False)
    #         Scales sky spectrum for subtraction if True
    #     plot : boolean (default = True)
    #         Show the plots in the console
    #     verbose : boolean (default = True)
    #         Print detailed description of steps taken in console
    #     fig_size : integer (default = 12)
    #         Size of the image plotted
    #     fibre_p: integer (default = -1)
    #         if fibre_p=fibre only corrects that fibre and plots the corrections, if -1, applies correction to all fibres
    #     kernel_correct_ccd_defects : odd integer (default = 51)
    #         width used for the median filter
    #     """
    #     self.history.append('- Sky sustraction using the 1Dfit method')

    #     if sky_spectrum_file != "":            
    #         sky_spectrum = self.read_sky_spectrum(sky_spectrum_file, path=path, verbose = verbose)
            
    #     if verbose:
    #         print("\n> Fitting sky lines in both a provided sky spectrum AND all the fibres")
    #         print("  This process takes ~20 minutes for 385R if all skylines are considered!\n")
    #     if len(sky_spectrum) == 0:
    #         if verbose:
    #             print("  No sky spectrum provided, using", n_sky, "lowest intensity fibres to create a sky...")
    #         self.history.append('  ERROR! No sky spectrum provided, using self method with n_sky =' + np.str(n_sky))
    #         self.find_sky_emission(n_sky=n_sky, plot=plot, sky_fibres=sky_fibres,
    #                                sky_wave_min=sky_wave_min, sky_wave_max=sky_wave_max, win_sky=win_sky)
    #     else:
    #         if scale_sky_1D != 0:
    #             self.history.append('  1D sky spectrum scaled by =' + np.str(scale_sky_1D))
    #             if verbose:
    #                 print("  1D sky spectrum scaled by ", scale_sky_1D)
    #         else:
    #             if verbose:
    #                 print("  No scale between 1D sky spectrum and object given, calculating...")
    #             if n_sky == 50: n_sky = 20
    #             self.find_sky_emission(n_sky=n_sky, plot=plot, sky_fibres=sky_fibres,
    #                                    sky_wave_min=sky_wave_min, sky_wave_max=sky_wave_max,
    #                                    win_sky=win_sky, include_history=False)

    #             sky_r_self = self.sky_emission

    #             scale_sky_1D = auto_scale_two_spectra(self, sky_r_self, sky_spectrum, scale=[0.1, 1.01, 0.025],
    #                                                   w_scale_min=self.valid_wave_min, w_scale_max=self.valid_wave_max,
    #                                                   plot=plot, verbose=True)

    #         self.history.append('  1D sky spectrum scaled by =' + np.str(scale_sky_1D))

    #         self.sky_emission = np.array(sky_spectrum) * scale_sky_1D

    #     self.fit_and_substract_sky_spectrum(self.sky_emission, sky_lines_file=sky_lines_file,
    #                                         brightest_line_wavelength=brightest_line_wavelength,
    #                                         brightest_line=brightest_line,
    #                                         maxima_sigma=maxima_sigma, ymin=-50, ymax=600, wmin=0, wmax=0,
    #                                         auto_scale_sky=auto_scale_sky,
    #                                         warnings=False, verbose=False, plot=False, fig_size=fig_size, fibre=fibre_p)

    #     if fibre_p == -1:
    #         if verbose:
    #             print("\n> 1Dfit sky_method usually generates some nans, correcting ccd defects again...")
    #         self.correct_ccd_defects(kernel_correct_ccd_defects=kernel_correct_ccd_defects, verbose=verbose, plot=plot,
    #                                  only_nans=True)  # Not replacing values <0
  
    
  
    
  
    
  
    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # def apply_selffit_sky(self, sky_spectrum=[], n_sky=50,  sky_fibres=[] , sky_spectrum_file="", path ="",
    #                       sky_wave_min=0, sky_wave_max=0, win_sky=0, scale_sky_1D=0,
    #                       sky_lines_file="", brightest_line_wavelength=0,
    #                       ranges_with_emission_lines = [0],
    #                       cut_red_end = 0,
    #                       brightest_line="Ha", maxima_sigma=3, auto_scale_sky=False,
    #                       fibre_p=-1, kernel_correct_ccd_defects=51,
    #                       plot=True, verbose=True, fig_size=12):
    #     """
    #     Subtract sky using the selffit method.

    #     Parameters
    #     ----------
    #     sky_spectrum : TYPE, optional
    #         DESCRIPTION. The default is [].
    #     n_sky : TYPE, optional
    #         DESCRIPTION. The default is 50.
    #     sky_fibres : TYPE, optional
    #         DESCRIPTION. The default is [].
    #     sky_spectrum_file : TYPE, optional
    #         DESCRIPTION. The default is "".
    #     path : TYPE, optional
    #         DESCRIPTION. The default is "".
    #     sky_wave_min : TYPE, optional
    #         DESCRIPTION. The default is 0.
    #     sky_wave_max : TYPE, optional
    #         DESCRIPTION. The default is 0.
    #     win_sky : TYPE, optional
    #         DESCRIPTION. The default is 0.
    #     scale_sky_1D : TYPE, optional
    #         DESCRIPTION. The default is 0.
    #     sky_lines_file : TYPE, optional
    #         DESCRIPTION. The default is "".
    #     brightest_line_wavelength : TYPE, optional
    #         DESCRIPTION. The default is 0.
    #     ranges_with_emission_lines : TYPE, optional
    #         DESCRIPTION. The default is [0].
    #     cut_red_end : TYPE, optional
    #         DESCRIPTION. The default is 0.
    #     brightest_line : TYPE, optional
    #         DESCRIPTION. The default is "Ha".
    #     maxima_sigma : TYPE, optional
    #         DESCRIPTION. The default is 3.
    #     auto_scale_sky : TYPE, optional
    #         DESCRIPTION. The default is False.
    #     fibre_p : TYPE, optional
    #         DESCRIPTION. The default is -1.
    #     kernel_correct_ccd_defects : TYPE, optional
    #         DESCRIPTION. The default is 51.
    #     plot : TYPE, optional
    #         DESCRIPTION. The default is True.
    #     verbose : TYPE, optional
    #         DESCRIPTION. The default is True.
    #     fig_size : TYPE, optional
    #         DESCRIPTION. The default is 12.





    #     """

    #     self.history.append('- Sky sustraction using the selffit method')

    #     if verbose: print("\n> 'sky_method = selffit', hence using", n_sky,
    #                       "lowest intensity fibres to create a sky spectrum ...")

    #     self.find_sky_emission(n_sky=n_sky, plot=plot, sky_fibres=sky_fibres,
    #                            sky_wave_min=sky_wave_min, sky_wave_max=sky_wave_max,
    #                            win_sky=win_sky, include_history=True)
        
    #     if sky_spectrum_file != "":            
    #         sky_spectrum = self.read_sky_spectrum(sky_spectrum_file, path=path, verbose = verbose)
        
    #     if sky_spectrum[0] != -1 and np.nanmedian(sky_spectrum) != 0:
    #         if verbose: print(
    #             "\n> Additional sky spectrum provided. Using this for replacing regions with bright emission lines...")

    #         sky_r_self = self.sky_emission

    #         self.sky_emission = replace_el_in_sky_spectrum(self, sky_r_self, sky_spectrum,
    #                                                        scale_sky_1D=scale_sky_1D,
    #                                                        brightest_line=brightest_line,
    #                                                        brightest_line_wavelength=brightest_line_wavelength,
    #                                                        ranges_with_emission_lines=ranges_with_emission_lines,
    #                                                        cut_red_end=cut_red_end,
    #                                                        plot=plot)
    #         self.history.append('  Using sky spectrum provided for replacing regions with emission lines')

    #     self.fit_and_substract_sky_spectrum(self.sky_emission, sky_lines_file=sky_lines_file,
    #                                         brightest_line_wavelength=brightest_line_wavelength,
    #                                         brightest_line=brightest_line,
    #                                         maxima_sigma=maxima_sigma, ymin=-50, ymax=600, wmin=0, wmax=0,
    #                                         auto_scale_sky=auto_scale_sky,
    #                                         warnings=False, verbose=False, plot=False, fig_size=fig_size,
    #                                         fibre=fibre_p)

    #     if fibre_p == -1:
    #         if verbose: print("\n> 'selffit' sky_method usually generates some nans, correcting ccd defects again...")
    #         self.correct_ccd_defects(kernel_correct_ccd_defects=kernel_correct_ccd_defects, verbose=verbose,
    #                                  plot=plot, only_nans=True)  # not replacing values < 0
    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # def apply_2D_sky(self, sky_rss, scale_sky_rss=0,
    #                  plot=True, verbose=True, fig_size=12):
    #     """
    #     Task that uses a RSS file with a sky, scale and sustract it.
    #     #TODO: this method needs to be checked and use plot_plot in plots

    #     Parameters
    #     ----------
    #     sky_rss : OBJECT #TODO This needs to be also a file
    #         A RSS file with offset sky.
    #     scale_sky_rss : float, optional
    #         scale applied to the SKY RSS before sustracting
    #     plot : TYPE, optional
    #         DESCRIPTION. The default is True.
    #     verbose : TYPE, optional
    #         DESCRIPTION. The default is True.
    #     fig_size : TYPE, optional
    #         DESCRIPTION. The default is 12.
    #     """
    #     self.history.append('- Sky sustraction using the 2D method')

    #     if scale_sky_rss != 0:
    #         if verbose: print("\n> Using sky image provided to substract sky, considering a scale of",
    #                           scale_sky_rss, "...")
    #         self.sky_emission = scale_sky_rss * sky_rss.intensity_corrected
    #         self.intensity_corrected = self.intensity_corrected - self.sky_emission
    #     else:
    #         if verbose: print(
    #             "\n> Using sky image provided to substract sky, computing the scale using sky lines")
    #         # check scale fibre by fibre
    #         self.sky_emission = copy.deepcopy(sky_rss.intensity_corrected)
    #         scale_per_fibre = np.ones((self.n_spectra))
    #         scale_per_fibre_2 = np.ones((self.n_spectra))
    #         lowlow = 15
    #         lowhigh = 5
    #         highlow = 5
    #         highhigh = 15
    #         if self.grating == "580V":
    #             if verbose: print("  For 580V we use bright skyline at 5577 AA ...")
    #             sky_line = 5577
    #             sky_line_2 = 0
    #         if self.grating == "1000R":
    #             # print "  For 1000R we use skylines at 6300.5 and 6949.0 AA ..."   ### TWO LINES GIVE WORSE RESULTS THAN USING ONLY 1...
    #             if verbose: print("  For 1000R we use skyline at 6949.0 AA ...")
    #             sky_line = 6949.0  # 6300.5
    #             lowlow = 22  # for getting a good continuuem in 6949.0
    #             lowhigh = 12
    #             highlow = 36
    #             highhigh = 52
    #             sky_line_2 = 0  # 6949.0  #7276.5 fails
    #             lowlow_2 = 22  # for getting a good continuuem in 6949.0
    #             lowhigh_2 = 12
    #             highlow_2 = 36
    #             highhigh_2 = 52
    #         if sky_line_2 != 0 and verbose: print("  ... first checking", sky_line, "...")
    #         for fibre_sky in range(self.n_spectra):
    #             skyline_spec = fluxes(self.wavelength, self.intensity_corrected[fibre_sky], sky_line,
    #                                   plot=False, verbose=False, lowlow=lowlow, lowhigh=lowhigh,
    #                                   highlow=highlow, highhigh=highhigh)  # fmin=-5.0E-17, fmax=2.0E-16,
    #             # resultado = [rms_cont, fit[0], fit_error[0], gaussian_flux, gaussian_flux_error, fwhm, fwhm_error, flux, flux_error, ew, ew_error, spectrum  ]
    #             self.intensity_corrected[fibre_sky] = skyline_spec[11]

    #             skyline_sky = fluxes(self.wavelength, self.sky_emission[fibre_sky], sky_line, plot=False,
    #                                  verbose=False, lowlow=lowlow, lowhigh=lowhigh, highlow=highlow,
    #                                  highhigh=highhigh)  # fmin=-5.0E-17, fmax=2.0E-16,

    #             scale_per_fibre[fibre_sky] = skyline_spec[3] / skyline_sky[3]
    #             self.sky_emission[fibre_sky] = skyline_sky[11]

    #         if sky_line_2 != 0:
    #             if verbose: print("  ... now checking", sky_line_2, "...")
    #             for fibre_sky in range(self.n_spectra):
    #                 skyline_spec = fluxes(self.wavelength, self.intensity_corrected[fibre_sky], sky_line_2,
    #                                       plot=False, verbose=False, lowlow=lowlow_2, lowhigh=lowhigh_2,
    #                                       highlow=highlow_2,
    #                                       highhigh=highhigh_2)  # fmin=-5.0E-17, fmax=2.0E-16,
    #                 # resultado = [rms_cont, fit[0], fit_error[0], gaussian_flux, gaussian_flux_error, fwhm, fwhm_error, flux, flux_error, ew, ew_error, spectrum  ]
    #                 self.intensity_corrected[fibre_sky] = skyline_spec[11]

    #                 skyline_sky = fluxes(self.wavelength, self.sky_emission[fibre_sky], sky_line_2, plot=False,
    #                                      verbose=False, lowlow=lowlow_2, lowhigh=lowhigh_2, highlow=highlow_2,
    #                                      highhigh=highhigh_2)  # fmin=-5.0E-17, fmax=2.0E-16,

    #                 scale_per_fibre_2[fibre_sky] = skyline_spec[3] / skyline_sky[3]
    #                 self.sky_emission[fibre_sky] = skyline_sky[11]

    #                 # Median value of scale_per_fibre, and apply that value to all fibres
    #         if sky_line_2 == 0:
    #             scale_sky_rss = np.nanmedian(scale_per_fibre)
    #             self.sky_emission = self.sky_emission * scale_sky_rss
    #         else:
    #             scale_sky_rss = np.nanmedian((scale_per_fibre + scale_per_fibre_2) / 2)
    #             # Make linear fit
    #             scale_sky_rss_1 = np.nanmedian(scale_per_fibre)
    #             scale_sky_rss_2 = np.nanmedian(scale_per_fibre_2)
    #             if verbose:
    #                 print("  Median scale for line 1 :", scale_sky_rss_1, "range [", np.nanmin(scale_per_fibre),
    #                       ",", np.nanmax(scale_per_fibre), "]")
    #                 print("  Median scale for line 2 :", scale_sky_rss_2, "range [",
    #                       np.nanmin(scale_per_fibre_2), ",", np.nanmax(scale_per_fibre_2), "]")

    #             b = (scale_sky_rss_1 - scale_sky_rss_2) / (sky_line - sky_line_2)
    #             a = scale_sky_rss_1 - b * sky_line
    #             if verbose: print("  Appling linear fit with a =", a, "b =", b,
    #                               "to all fibres in sky image...")  # ,a+b*sky_line,a+b*sky_line_2

    #             for i in range(self.n_wave):
    #                 self.sky_emission[:, i] = self.sky_emission[:, i] * (a + b * self.wavelength[i])

    #         if plot:
    #             plt.figure(figsize=(fig_size, fig_size / 2.5))
    #             label1 = "$\lambda$" + np.str(sky_line)
    #             plt.plot(scale_per_fibre, alpha=0.5, label=label1)
    #             plt.minorticks_on()
    #             plt.ylim(np.nanmin(scale_per_fibre), np.nanmax(scale_per_fibre))
    #             plt.axhline(y=scale_sky_rss, color='k', linestyle='--')
    #             if sky_line_2 == 0:
    #                 text = "Scale OBJECT / SKY using sky line $\lambda$" + np.str(sky_line)
    #                 if verbose:
    #                     print("  Scale per fibre in the range [", np.nanmin(scale_per_fibre), ",",
    #                           np.nanmax(scale_per_fibre), "], median value is", scale_sky_rss)
    #                     print("  Using median value to scale sky emission provided...")
    #             if sky_line_2 != 0:
    #                 text = "Scale OBJECT / SKY using sky lines $\lambda$" + np.str(
    #                     sky_line) + " and $\lambda$" + np.str(sky_line_2)
    #                 label2 = "$\lambda$" + np.str(sky_line_2)
    #                 plt.plot(scale_per_fibre_2, alpha=0.5, label=label2)
    #                 plt.axhline(y=scale_sky_rss_1, color='k', linestyle=':')
    #                 plt.axhline(y=scale_sky_rss_2, color='k', linestyle=':')
    #                 plt.legend(frameon=False, loc=1, ncol=2)
    #             plt.title(text)
    #             plt.xlabel("Fibre")
    #             plt.show()
    #             plt.close()
    #         self.intensity_corrected = self.intensity_corrected - self.sky_emission
    #     self.apply_mask(verbose=verbose)
    #     self.history(" - 2D sky subtraction performed")
    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    # # -----------------------------------------------------------------------------
    
    
    
    
    
    
    
    # def is_sky(self, n_sky=50, win_sky=0, sky_fibres=[], sky_wave_min=0,
    #            sky_wave_max=0, plot=True, verbose=True):
    #     """
    #     If frame is an empty sky, apply median filter to it for further use in 2D sky emission fitting

    #     Parameters
    #     ----------
    #     n_sky : integer (default = 50)
    #         Number of fibres to use for finding sky spectrum
    #     sky_fibres : list of integers (default = none)
    #         Specify the fibres to use to obtain sky spectrum. Will automatically determine the best fibres if not specified
    #     sky_wave_min : float (default = 0)
    #         Specify the lower bound on wavelength range. If 0, it is set to self.valid_wave_min
    #     sky_wave_max : float (default = 0)
    #         Specify the upper bound on wavelength range. If 0, it is set to self.valid_wave_max
    #     win_sky : odd integer (default = 0)
    #         Width in fibres of a median filter applied to obtain sky spectrum, if 0, it will not apply any median filter
    #      plot : boolean (default = True)
    #         Show the plots in the console
    #     verbose : boolean (default = True)
    #         Print detailed description of steps taken in console
    #     """

    #     if verbose: print("\n> This RSS file is defined as SKY... identifying", n_sky,
    #                       " lowest fibres for getting 1D sky spectrum...")
    #     self.history.append('- This RSS file is defined as SKY:')
    #     self.find_sky_emission(n_sky=n_sky, plot=plot, sky_fibres=sky_fibres,
    #                            sky_wave_min=sky_wave_min, sky_wave_max=sky_wave_max, win_sky=0)
    #     # print "\n> This RSS file is defined as SKY... applying median filter with window",win_sky,"..."
    #     if win_sky == 0:  # Default when it is not a win_sky
    #         win_sky = 151
    #     print("\n  ... applying median filter with window", win_sky, "...\n")

    #     medfilt_sky = median_2D_filter(self.intensity_corrected, self.n_spectra, self.n_wave, win_sky=win_sky)
    #     self.intensity_corrected = copy.deepcopy(medfilt_sky)
    #     print("  Median filter applied, results stored in self.intensity_corrected !")
    #     self.history.append('  Median filter ' + np.str(win_sky) + ' applied to all fibres')

