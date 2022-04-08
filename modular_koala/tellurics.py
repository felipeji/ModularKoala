# =============================================================================
# Basics packages
# =============================================================================
import numpy as np
import copy
import matplotlib.pyplot as plt
from scipy.signal import medfilt

# =============================================================================
# Astropy and associated packages
# =============================================================================
from astropy.table import Table

# =============================================================================
# KOALA packages
# =============================================================================
# Modular
from modular_koala.ancillary import vprint
from modular_koala.rss import detect_edge

# Original
from koala.onedspec import smooth_spectrum
from koala.plot_plot import plot_plot
from koala.io import spectrum_to_text_file
# =============================================================================


# def obtain_telluric_correction(w, telluric_correction_list, plot=True, label_stars=[], scale=[]):
#     if len(scale) == 0:
#         for star in telluric_correction_list: scale.append(1.)

#     for i in range(len(telluric_correction_list)):
#         telluric_correction_list[i] = [1. if x * scale[i] < 1 else x * scale[i] for x in telluric_correction_list[i]]

#     telluric_correction = np.nanmedian(telluric_correction_list, axis=0)
#     if plot:
#         fig_size = 12
#         plt.figure(figsize=(fig_size, fig_size / 2.5))
#         plt.title("Telluric correction")
#         for i in range(len(telluric_correction_list)):
#             if len(label_stars) > 0:
#                 label = label_stars[i]
#             else:
#                 label = "star" + str(i + 1)
#             plt.plot(w, telluric_correction_list[i], alpha=0.3, label=label)
#         plt.plot(w, telluric_correction, alpha=0.5, color="k", label="Median")
#         plt.minorticks_on()
#         plt.legend(frameon=False, loc=2, ncol=1)
#         step_up = 1.15 * np.nanmax(telluric_correction)
#         plt.ylim(0.9, step_up)
#         plt.xlim(w[0] - 10, w[-1] + 10)
#         plt.show()
#         plt.close()

#     print("\n> Telluric correction = ", telluric_correction)
#     if np.nanmean(scale) != 1.: print("  Telluric correction scale provided : ", scale)
#     print("\n  Telluric correction obtained!")
#     return telluric_correction





class Tellurics():
    """

    """

    def __init__(self, rss,
                        n_fibres=10,
                        correct_from=6850.,
                        correct_to=10000.,
                        save_telluric_file="",
                        apply_tc=False,
                        step=10,
                        is_combined_cube=False,
                        weight_fit_median=0.5,
                        exclude_wlm=[[6450, 6700], [6850, 7050], [7130, 7380]],  # This is range for 1000R
                        wave_min=None,
                        wave_max=None,
                        plot=True,
                        fig_size=12,
                        verbose=False):
    
        
        # Set print verbose
        vprint.verbose = verbose
        vprint("\n> Obtaining telluric correction using spectrophotometric star...")
    
        # TODO: adapt the to analize RSS or a  new CUBE class containing the cubes
        if is_combined_cube:
            wlm = rss.combined_cube.wavelength
        else:
            wlm = rss.wavelength
        
        
        if wave_min is None: wave_min = wlm[0]
        if wave_max is None: wave_max = wlm[-1]
    
        if is_combined_cube:
            if rss.combined_cube.seeing == 0:
                rss.combined_cube.half_light_spectrum(5, plot=plot, min_wave=wave_min, max_wave=wave_max)
            estrella = rss.combined_cube.integrated_star_flux
        else:
            intensidad = rss.intensity_corrected
            integrated_fibre = np.nansum(rss.intensity_corrected,axis=1)
            # The n_fibres brightest fibre 
            brightest_fibres = (-integrated_fibre).argsort()[:n_fibres]
            # TODO: Wouldn't be better take median instead of the sum brightest_fibres ?
            estrella = np.nansum(intensidad[brightest_fibres], axis=0)
    
        smooth_med_star = smooth_spectrum(wlm, estrella, wave_min=wave_min, wave_max=wave_max, step=step,
                                          weight_fit_median=weight_fit_median,
                                          exclude_wlm=exclude_wlm, plot=plot, verbose=verbose)
    
        telluric_correction = np.ones(len(wlm))
    
        estrella_m = medfilt(estrella, 151)
        plot_plot(wlm, [estrella, smooth_med_star, estrella_m])
    
        # Avoid H-alpha absorption
        rango_ha = [0, 0]
        for rango in exclude_wlm:
            if rango[0] < 6563. and rango[1] > 6563.:  # H-alpha is here, skip
                print("  Skipping range with H-alpha...")
                rango_ha = rango
    
        # TODO: REMOVE LOOPS
        correct_from = 6000.
        for l in range(len(wlm)):
            if wlm[l] > correct_from and wlm[l] < correct_to:
    
                if wlm[l] > rango_ha[0] and wlm[l] < rango_ha[1]:
                    step = step + 0
                    # skipping Ha
                else:
                    telluric_correction[l] = smooth_med_star[l] / estrella[l]
    
        waves_for_tc_ = []
        for rango in exclude_wlm:
            if rango[0] < 6563. and rango[1] > 6563.:  # H-alpha is here, skip
                print("  Skipping range with H-alpha...")
            else:
                index_region = np.where((wlm >= rango[0]) & (wlm <= rango[1]))
                waves_for_tc_.append(index_region)
    
        waves_for_tc = []
        for rango in waves_for_tc_:
            waves_for_tc = np.concatenate((waves_for_tc, rango[0].tolist()), axis=None)
    
        # Now, change the value in telluric_correction
        for index in waves_for_tc:
            i = np.int(index)
            if smooth_med_star[i] / estrella[i] > 1.:
                telluric_correction[i] = smooth_med_star[i] / estrella[i]
    
        if plot:
            plt.figure(figsize=(fig_size, fig_size / 2.5))
            if is_combined_cube:
                print("  Telluric correction for this star (" + rss.combined_cube.object + ") :")
                plt.plot(wlm, estrella, color="b", alpha=0.3)
                plt.plot(wlm, estrella * telluric_correction, color="g", alpha=0.5)
                plt.ylim(np.nanmin(estrella), np.nanmax(estrella))
            else:
                print("  Example of telluric correction using fibres", brightest_fibres[0], " (blue) and ", brightest_fibres[1],
                      " (green):")
                plt.plot(wlm, intensidad[brightest_fibres[0]], color="b", alpha=0.3)
                plt.plot(wlm, intensidad[brightest_fibres[0]] * telluric_correction, color="g", alpha=0.5)
                plt.plot(wlm, intensidad[brightest_fibres[1]], color="b", alpha=0.3)
                plt.plot(wlm, intensidad[brightest_fibres[1]] * telluric_correction, color="g", alpha=0.5)
                plt.ylim(np.nanmin(intensidad[brightest_fibres[0]]),
                         np.nanmax(intensidad[brightest_fibres[0]]))  # CHECK THIS AUTOMATICALLY
            plt.axvline(x=wave_min, color='k', linestyle='--')
            plt.axvline(x=wave_max, color='k', linestyle='--')
            plt.xlim(wlm[0] - 10, wlm[-1] + 10)
            plt.xlabel("Wavelength [$\mathrm{\AA}$]")
            if exclude_wlm[0][0] != 0:
                for i in range(len(exclude_wlm)):
                    plt.axvspan(exclude_wlm[i][0], exclude_wlm[i][1], color='r', alpha=0.1)
            plt.minorticks_on()
            plt.show()
            plt.close()
        
        # Wave asignated self variables to the oabtained telluric correction
        self.wlm = wlm
        self.telluric_correction = telluric_correction



    def apply(self,rss,
              verbose=False,
              is_combined_cube=False):
        
        # Set print verbose
        vprint.verbose = verbose
        
        # Copy input RSS for storage the changes implemented in the task   
        rss_out = copy.deepcopy(rss)

        
        vprint("  Applying telluric correction to this star...")
        if is_combined_cube:
            rss.combined_cube.integrated_star_flux = rss.combined_cube.integrated_star_flux * self.telluric_correction
            for i in range(rss.combined_cube.n_rows):
                for j in range(rss.combined_cube.n_cols):
                    rss.combined_cube.data[:, i, j] = rss.combined_cube.data[:, i, j] * self.telluric_correction
        else:
            for i,flux in enumerate(rss.intensity_corrected):
                rss_out.intensity_corrected[i] = flux * self.telluric_correction

        """
        if is_combined_cube:
            rss.combined_cube.telluric_correction = telluric_correction
        else:
            self.telluric_correction = telluric_correction

            """
        return rss_out

        

    def save_to_txt(self,filename='telluric_correction.txt'):
        spectrum_to_text_file(self.wlm, self.telluric_correction, filename=filename)





def tellurics_from_file(rss, 
                        telluric_file,
                        plot=True,
                        fig_size=12,
                        verbose=True):
    """

    Apply telluric correction to the data from a 1D telluric correction model.

    Parameters
    ----------
    telluric_correction_file : string (default = none)
        Path to the file containing data necessary to apply telluric correction
    telluric_correction : list of floats (default = [0])
        Table data from the telluric correction file in format of Python list
    plot : boolean (default = True)
        Show the plots in the console
    fig_size : float (default = 12)
        Size of the plots
    verbose : boolean (default = True)
        Print detailed description of steps taken in console
    """
            
    # Set print verbose
    vprint.verbose = verbose
    
    # Copy input RSS for storage the changes implemented in the task   
    rss_out = copy.deepcopy(rss)

    
    vprint("\n> Reading file with the telluric correction: ")
    vprint(" ", telluric_file)
    telluric_correction  = Table.read(telluric_file,format='ascii')['col2']
 

    vprint("\n> Applying telluric correction...")
    rss_out.intensity_corrected *= telluric_correction[np.newaxis, :]
    rss_out.variance_corrected *= telluric_correction[np.newaxis, :]**2

    valid_wave_min, min_index, valid_wave_max, max_index = detect_edge(rss)

    if plot:
        plot_plot(rss.wavelength, telluric_correction, xmin=rss.wavelength[0] - 10,
                  xmax=rss.wavelength[-1] + 10, statistics=False,
                  ymin=0.9, ymax=2, ptitle="Telluric correction", xlabel="Wavelength [$\mathrm{\AA}$]",
                  vlines=[valid_wave_min, valid_wave_max])

        integrated_fibre = np.nansum(rss.intensity_corrected,axis=1)
        faintest = integrated_fibre.argsort()[0]    
        brightest = integrated_fibre.argsort()[1]

        
        
        print("  Example of telluric correction using faintest fibre", faintest, ":")
        ptitle = "Telluric correction in fibre " + str(faintest)
        plot_plot(rss.wavelength, [rss.intensity_corrected[faintest], rss_out.intensity_corrected[faintest]],
                  xmin=rss.wavelength[0] - 10, xmax=rss.wavelength[-1] + 10,
                  ymin=np.nanpercentile(rss.intensity_corrected[faintest], 1),
                  ymax=np.nanpercentile(rss.intensity_corrected[faintest], 99),
                  vlines=[valid_wave_min, valid_wave_max],
                  xlabel="Wavelength [$\mathrm{\AA}$]", ptitle=ptitle)
        
        print("  Example of telluric correction using brightest fibre", brightest, ":")
        ptitle = "Telluric correction in fibre " + str(brightest)
        plot_plot(rss.wavelength, [rss.intensity_corrected[brightest], rss_out.intensity_corrected[brightest]],
                  xmin=rss.wavelength[0] - 10, xmax=rss.wavelength[-1] + 10,
                  ymin=np.nanpercentile(rss.intensity_corrected[brightest], 1),
                  ymax=np.nanpercentile(rss.intensity_corrected[brightest], 99),
                  vlines=[valid_wave_min, valid_wave_max],
                  xlabel="Wavelength [$\mathrm{\AA}$]", ptitle=ptitle)

        # history.append("- Telluric correction applied reading from file:")
        # history.append("  " + telluric_correction_file)


    return rss_out





# def get_telluric_correction(rss,
#                             n_fibres=10,
#                             correct_from=6850.,
#                             correct_to=10000.,
#                             save_telluric_file="",
#                             apply_tc=False,
#                             step=10,
#                             is_combined_cube=False,
#                             weight_fit_median=0.5,
#                             exclude_wlm=[[6450, 6700], [6850, 7050], [7130, 7380]],  # This is range for 1000R
#                             wave_min=None,
#                             wave_max=None,
#                             plot=True,
#                             fig_size=12,
#                             verbose=False):
    
    
#     """
#     Get telluric correction using a spectrophotometric star

#     IMPORTANT: check tasks "telluric_correction_from_star" and "telluric_correction_using_bright_continuum_source"

#     Parameters
#     ----------
#     n_fibres: integer (default = 10)
#         number of fibers to add for obtaining spectrum
#     correct_from :  float (default = 6850)
#         wavelength from which telluric correction is applied
#     correct_to :  float (default = 10000)
#         last wavelength where telluric correction is applied
#     save_telluric_file : string (default = "")
#         If given, it saves the telluric correction in a file with that name
#     apply_tc : boolean (default = False)
#         apply telluric correction to data
#         Only do this when absolutely sure the correction is good!
#     step : integer (default = 10)
#        step using for estimating the local medium value
#     is_combined_cube : boolean (default = False)
#         Use True if the cube is a combined cube and needs to read from self.combined_cube
#     weight_fit_median : float between 0 and 1 (default = 0.5)
#         weight of the median value when calling task smooth_spectrum
#     exclude_wlm : list of [float, float] (default = [[6450,6700],[6850,7050], [7130,7380]] )
#         Wavelength ranges not considering for normalising stellar continuum
#         The default values are for 1000R grating
#     wave_min, wave_max : float (default = 0,0)
#         Wavelength range to consider, [wave_min, wave_max]
#         if 0, it uses  wave_min=wavelength[0] and wave_max=wavelength[-1]
#     plot : boolean (default = True)
#         Plot
#     fig_size: float (default = 12)
#        Size of the figure
#     verbose : boolean (default = True)
#         Print results

#     Example
#     ----------
#     telluric_correction_star1 = star1r.get_telluric_correction(n_fibres=15,
#                                 exclude_wlm= [ [6245,6390],[6450,6750],[6840,7000],
#                                 [7140,7400],[7550,7720],[8050,8450]])
#     """
    
#     # Set print verbose
#     vprint.verbose = verbose
    
# # =============================================================================
# # Copy input RSS for storage the changes implemented in the task   
# # =============================================================================
#     rss_out = copy.deepcopy(rss)
    
#     vprint("\n> Obtaining telluric correction using spectrophotometric star...")

#     # TODO: adapt the to analize RSS or a  new CUBE class containing the cubes
#     if is_combined_cube:
#         wlm = rss.combined_cube.wavelength
#     else:
#         wlm = rss.wavelength
    
    
#     if wave_min is None: wave_min = wlm[0]
#     if wave_max is None: wave_max = wlm[-1]

#     if is_combined_cube:
#         if rss.combined_cube.seeing == 0:
#             rss.combined_cube.half_light_spectrum(5, plot=plot, min_wave=wave_min, max_wave=wave_max)
#         estrella = rss.combined_cube.integrated_star_flux
#     else:
#         intensidad = rss.intensity_corrected
#         integrated_fibre = np.nansum(rss.intensity_corrected,axis=1)
#         # The n_fibres brightest fibre 
#         brightest_fibres = (-integrated_fibre).argsort()[:n_fibres]
#         # TODO: Wouldn't be better take median instead of the sum brightest_fibres ?
#         estrella = np.nansum(intensidad[brightest_fibres], axis=0)

#     smooth_med_star = smooth_spectrum(wlm, estrella, wave_min=wave_min, wave_max=wave_max, step=step,
#                                       weight_fit_median=weight_fit_median,
#                                       exclude_wlm=exclude_wlm, plot=plot, verbose=verbose)

#     telluric_correction = np.ones(len(wlm))

#     estrella_m = medfilt(estrella, 151)
#     plot_plot(wlm, [estrella, smooth_med_star, estrella_m])

#     # Avoid H-alpha absorption
#     rango_ha = [0, 0]
#     for rango in exclude_wlm:
#         if rango[0] < 6563. and rango[1] > 6563.:  # H-alpha is here, skip
#             print("  Skipping range with H-alpha...")
#             rango_ha = rango

#     # TODO: REMOVE LOOPS
#     correct_from = 6000.
#     for l in range(len(wlm)):
#         if wlm[l] > correct_from and wlm[l] < correct_to:

#             if wlm[l] > rango_ha[0] and wlm[l] < rango_ha[1]:
#                 step = step + 0
#                 # skipping Ha
#             else:
#                 telluric_correction[l] = smooth_med_star[l] / estrella[l]

#     waves_for_tc_ = []
#     for rango in exclude_wlm:
#         if rango[0] < 6563. and rango[1] > 6563.:  # H-alpha is here, skip
#             print("  Skipping range with H-alpha...")
#         else:
#             index_region = np.where((wlm >= rango[0]) & (wlm <= rango[1]))
#             waves_for_tc_.append(index_region)

#     waves_for_tc = []
#     for rango in waves_for_tc_:
#         waves_for_tc = np.concatenate((waves_for_tc, rango[0].tolist()), axis=None)

#     # Now, change the value in telluric_correction
#     for index in waves_for_tc:
#         i = np.int(index)
#         if smooth_med_star[i] / estrella[i] > 1.:
#             telluric_correction[i] = smooth_med_star[i] / estrella[i]

#     if plot:
#         plt.figure(figsize=(fig_size, fig_size / 2.5))
#         if is_combined_cube:
#             print("  Telluric correction for this star (" + rss.combined_cube.object + ") :")
#             plt.plot(wlm, estrella, color="b", alpha=0.3)
#             plt.plot(wlm, estrella * telluric_correction, color="g", alpha=0.5)
#             plt.ylim(np.nanmin(estrella), np.nanmax(estrella))
#         else:
#             print("  Example of telluric correction using fibres", brightest_fibres[0], " (blue) and ", brightest_fibres[1],
#                   " (green):")
#             plt.plot(wlm, intensidad[brightest_fibres[0]], color="b", alpha=0.3)
#             plt.plot(wlm, intensidad[brightest_fibres[0]] * telluric_correction, color="g", alpha=0.5)
#             plt.plot(wlm, intensidad[brightest_fibres[1]], color="b", alpha=0.3)
#             plt.plot(wlm, intensidad[brightest_fibres[1]] * telluric_correction, color="g", alpha=0.5)
#             plt.ylim(np.nanmin(intensidad[brightest_fibres[0]]),
#                      np.nanmax(intensidad[brightest_fibres[0]]))  # CHECK THIS AUTOMATICALLY
#         plt.axvline(x=wave_min, color='k', linestyle='--')
#         plt.axvline(x=wave_max, color='k', linestyle='--')
#         plt.xlim(wlm[0] - 10, wlm[-1] + 10)
#         plt.xlabel("Wavelength [$\mathrm{\AA}$]")
#         if exclude_wlm[0][0] != 0:
#             for i in range(len(exclude_wlm)):
#                 plt.axvspan(exclude_wlm[i][0], exclude_wlm[i][1], color='r', alpha=0.1)
#         plt.minorticks_on()
#         plt.show()
#         plt.close()

#     if apply_tc:  # Check this
#         print("  Applying telluric correction to this star...")
#         if is_combined_cube:
#             rss.combined_cube.integrated_star_flux = rss.combined_cube.integrated_star_flux * telluric_correction
#             for i in range(rss.combined_cube.n_rows):
#                 for j in range(rss.combined_cube.n_cols):
#                     rss.combined_cube.data[:, i, j] = rss.combined_cube.data[:, i, j] * telluric_correction
#         else:
#             for i,flux in enumerate(rss.intensity_corrected):
#                 rss_out.intensity_corrected[i] = flux * telluric_correction
#     else:
#         print("  As apply_tc = False , telluric correction is NOT applied...")

#     """
#     if is_combined_cube:
#         rss.combined_cube.telluric_correction = telluric_correction
#     else:
#         self.telluric_correction = telluric_correction

#     # save file if requested
#     if save_telluric_file != "":
#         spectrum_to_text_file(wlm, telluric_correction, filename=save_telluric_file)
#         """
#     return rss_out,telluric_correction

