from koala.constants import red_gratings
from koala.io import full_path, read_table, spectrum_to_text_file, list_fits_files_in_folder
from astropy.io import fits

#from koala.KOALA_RSS import KOALA_RSS
from koala.rss import RSS, get_throughput_2D
from koala.cube import plot_response,obtain_flux_calibration,obtain_telluric_correction
from koala.automatic_scripts.run_automatic_star import run_automatic_star
import numpy as np

import os




def is_calib_star(value):
    """
    This function checks whether the string <value> is associated with 
    any of the known calibration stars described in the calib_dict dictionary
    by searching if <value> exist in any of the lists of known names.
    
    
    If one of the  <value> is identified (exist in the name list), this function 
    returns the corresponding name of star  as interpreted by KOALA (key value
    of the dictyonaty). Otherwise returns None. 


    Parameters
    ----------
    value : str
        String with the name of the star.

    Returns
    -------
    found_star : str, False
        pretty name or False.

    """
    # =============================================================================
    # Dictionary describing the known calibration stars 
    # =============================================================================
    # key : preatty name (str) 
    # value: list of posible names ([str, str, str])
    
    calib_dict = { "Hilt600" : ["H600", "HILT600", "Hilt600", "Hiltner600", "HILTNER600"], 
                   "EG274" : ["EG274", "Eg274", "eg274", "eG274", "E274", "e274"],
                   "HD60753" : ["HD60753", "hd60753", "Hd60753", "HD60753FLUX"],
                   "HD49798" : ["HD49798", "hd49798", "Hd49798"],
                   "CD32d9927" : ["cd32d9927", "CD32d9927", "CD32D9927", "CD-32d9927", "cd-32d9927", "Cd-32d9927", "CD-32D9927", "cd-32D9927", "Cd-32D9927"  ],
                   "HR3454" : ["HR3454", "Hr3454", "hr3454"],
                   "HR718" : [ "HR718" ,"Hr718" , "hr718", "HR718FLUX","HR718auto" ,"Hr718auto" , "hr718auto", "HR718FLUXauto"], 
                 }
    
    found_star = False
    for star, names in calib_dict.items():  
            if value in names:
                found_star = star
    
    return found_star


def str_identify(string):
    """
    This funtions evaluates a string. Returns a number or a boolean if it 
    identifiyed, If not terurs the strigs

    Parameters
    ----------
    string : str
        Varaiable to be evaluated.

    Returns
    -------
        Returs the variable identifyied.

    """
    try:
        return eval(string)
    except:
        return string







def automatic_calibration_night(path,
                                CALIBRATION_NIGHT_FILE = None,
                                date = "",
                                grating = "",
                                pixel_size = 0,
                                kernel_size=0,
                                file_skyflat="",
                                throughput_2D_file="",
                                throughput_2D = 0,
                                skyflat=0,
                                do_skyflat = True, 
                                kernel_throughput =0,
                                correct_ccd_defects=True,
                                fix_wavelengths = False,
                                sol=[-1],
                                rss_star_file_for_sol = "",
                                plot=True,
                                CONFIG_FILE_path="",
                                CONFIG_FILE_list=[],
                                star_list=[],
                                abs_flux_scale=[],
                                flux_calibration_file="",
                                telluric_correction_file ="",
                                objects_auto=[],
                                auto = False,
                                rss_clean = False,
                                flux_calibration_name="flux_calibration_auto",
                                cal_from_calibrated_starcubes = False,
                                disable_stars=[],                      # stars in this list will not be used
                                skyflat_names = ["SKYFLAT", "skyflat", "Skyflat", "SkyFlat", "SKYFlat", "SkyFLAT"],
                                list_of_objects = [],
                                instrument = "KOALA"
                                ):
    """
    Use: 
        CALIBRATION_NIGHT_FILE = "./configuration_files/calibration_night.config"
        automatic_calibration_night(CALIBRATION_NIGHT_FILE)
    """
    
    # =============================================================================
    # Dictionary with the calibration parameters     
    # =============================================================================
    
    skyflat_variable = ""
    throughput_2D_variable = "" 

    params = {"path": path,  
            "file_skyflat": file_skyflat, 
            "rss_star_file_for_sol": rss_star_file_for_sol,
            "grating": grating,
            "pixel_size": pixel_size,
            "kernel_size": kernel_size,
            "date": date,
            "flux_calibration_file": flux_calibration_file,
            "telluric_correction_file": telluric_correction_file,
            "fix_wavelengths": fix_wavelengths,	 	
            "kernel_throughput": kernel_throughput,	 
           "throughput_2D_file": throughput_2D_file,	
            "throughput_2D_variable": throughput_2D_variable,	
            "CONFIG_FILE_path": CONFIG_FILE_path,	
            "CONFIG_FILE_list": CONFIG_FILE_list,	
            "plot": plot,	 		
            "cal_from_calibrated_starcubes": cal_from_calibrated_starcubes,
            "objects_auto": objects_auto,		
            "skyflat_variable": skyflat_variable,	
            "do_skyflat": do_skyflat,		
            "correct_ccd_defects": correct_ccd_defects, 	 
            "sol": sol, 			
            "abs_flux_scale": abs_flux_scale,
            "throughput_2D": throughput_2D,
            "skyflat": skyflat,
            }

    
    print(list_of_objects)
    
    
    
    
    
    # obj is a dict containing the resulting objects
    obj = {}

    w = []
    telluric_correction_list=[]  

    global flux_calibration_night
    global telluric_correction_night
    global throughput_2D_
    throughput_2D_ = [0]
         
    
    check_nothing_done = 0

    
    # =============================================================================
    # Do automatic calibration (auto = True) ?  
    # =============================================================================
    if auto:
        print("\n# =============================================================================")
        print("\n    COMPLETELY AUTOMATIC CALIBRATION OF THE NIGHT ")
        print("\n# =============================================================================")
        
        # Fix wavelengths in automatic calibration
        fix_wavelengths = True 
        
        # =============================================================================
        #  if no list_object is given the objects are automatically readed from <path>     
        # =============================================================================
        if list_of_objects == []:
            list_of_objects, list_of_files, list_of_exptimes, params["date"], params["grating"] = list_fits_files_in_folder(path, return_list=True)

        else:
            objects_auto = list_of_objects
            
            
            
    
        # =============================================================================


        # =============================================================================
        # Loop for identify objects types
        # =============================================================================
        
        list_of_files_of_stars = []
        
        for object_name, fits_files in zip(list_of_objects, list_of_files):
            
            # If any of the objects is a sky flatfield
            if object_name in skyflat_names:
                params["file_skyflat"] = fits_files[0]
                print ("  - SKYFLAT automatically identified")
            
            # If any of the objects is a know calibration star
            calib_name = is_calib_star(object_name)
            
            if calib_name and object_name not in disable_stars:
                print ("  - Calibration star "+ calib_name +" automatically identified")
                star_list.append(calib_name)                
                objects_auto.append(calib_name+"_"+grating+"_"+date)                                 
                params["rss_star_file_for_sol"] = fits_files[0]
                list_of_files_of_stars.append(fits_files)
                CONFIG_FILE_list.append("")            

        
            
        # =============================================================================

        # =============================================================================
        # throughput_2D_file ?
        # =============================================================================
        if throughput_2D_file != "":
            throughput_2D_file = full_path(throughput_2D_file,path) 
            do_skyflat = False
            print ("  - throughput_2D_file provided, no need of processing skyflat")          
            sol=[0,0,0]
            ftf = fits.open(throughput_2D_file)
            if ftf[0].data[0][0] == 1. :
                sol = [ftf[0].header["SOL0"],ftf[0].header["SOL1"],ftf[0].header["SOL2"]]

                print ("  - solution for fixing small wavelength shifts included in this file :\n    sol = ",sol)  

        
    else:
        list_of_files_of_stars=[[],[],[],[],[],[]]

    # =============================================================================

    
    # =============================================================================
    # is any CALIBRATION_NIGHT_FILE given ?        
    # =============================================================================
    if  CALIBRATION_NIGHT_FILE:
        # Generates a dictionary fqiet the reduccion parameters from 
        keys, values = read_table(CALIBRATION_NIGHT_FILE, ["s", "s"] )
        readed_params = dict(zip(keys,values))
        # Update params from CALIBRATION_NIGHT_FILE
        params.update(readed_params)
        
        print("\n> Reading configuration file ", CALIBRATION_NIGHT_FILE)
        print("  for performing the automatic calibration of the night...\n")
        
    else:
        print("\n# =============================================================================")
        print("\n> Using the values given in automatic_calibration_night()")
        print("  for performing the automatic calibration of the night...")
        print("\n# =============================================================================")
        # Takes parameter from user

        # Checkin  parameters needed are present        
        if not params['grating']:
            print ("  - No information about the grating provided")
        if not params['date']:
            print ("  - No information about the date provided")
        if pixel_size == 0:
            print ("  - No pixel size provided, considering pixel_size = 0.7")
            params["pixel_size"] = 0.7
        if kernel_size == 0:
            print ("  - No kernel size provided, considering kernel_size = 1.1")
            params["kernel_size"] = 1.1
        
        # label with the pixel and kernel size
        pixel_size_str = str(pixel_size).replace('.','p')    
        kernel_size_str = str(kernel_size).replace('.','k')
        pk = '_'+pixel_size_str+'_'+kernel_size_str
        
        #pk = "_"+str(int(pixel_size))+"p"+str(int((abs(pixel_size)-abs(int(pixel_size)))*10))+"_"+str(int(kernel_size))+"k"+str(int(abs(kernel_size*100))-int(kernel_size)*100)
        
        if sol[0] != 0 : fix_wavelengths = True
        if len(CONFIG_FILE_path) > 0:
            for i in range(len(CONFIG_FILE_list)):
                CONFIG_FILE_list[i] = full_path (CONFIG_FILE_list[i],CONFIG_FILE_path)                

 

    # Variable allocation from params dictionary
    date = params["date"]
    grating = params["grating"]
    pixel_size = float(params["pixel_size"])
    kernel_size = float(params['kernel_size'])
    path = params['path']

    params["throughput_2D_file"] = os.path.join(path,"throughput_2D_"+date+"_"+grating+".fits")
    params["flux_calibration_file"] = os.path.join(path,"flux_calibration_"+date+"_"+grating+pk+".dat") 

    if flux_calibration_name =="flux_calibration_auto" : flux_calibration_name = "flux_calibration_"+date+"_"+grating+pk 
    if grating in red_gratings:
        params["telluric_correction_file"] = os.path.join(path,"telluric_correction_"+date+"_"+grating+".dat" )
        telluric_correction_name = os.path.join("telluric_correction_"+date+"_"+grating)

    file_skyflat = full_path(params["file_skyflat"],path)
    rss_star_file_for_sol = full_path(params["rss_star_file_for_sol"],path)
    flux_calibration_file = full_path(params["flux_calibration_file"],path)
    telluric_correction_file = full_path(params["telluric_correction_file"],path)


    fix_wavelengths = str_identify(params[ "fix_wavelengths"])
    kernel_throughput = int(params["kernel_throughput"])  

    throughput_2D_file = full_path(params["throughput_2D_file"],path)         
    throughput_2D_variable = params["throughput_2D"]  

    CONFIG_FILE_path = params["CONFIG_FILE_path"]
    try:
        CONFIG_FILE_list.append(full_path(params["CONFIG_FILE"],CONFIG_FILE_path))
    except:
        pass
    plot = str_identify(params["plot"])
    cal_from_calibrated_starcubes = str_identify(params["cal_from_calibrated_starcubes"]) 
    try:
        objects_auto.append(params["object"])
    except:
        pass
    skyflat_variable = params["skyflat"] # Originally set as global variable
    do_skyflat = str_identify(params["do_skyflat"])
    correct_ccd_defects = str_identify(params["correct_ccd_defects"])    
        

    try:
        sol_ = params["sol"].strip('][').split(',')
        fix_wavelengths = True
        if float(sol_[0]) == -1:
                sol = [float(sol_[0])]
        else:
            if float(sol_[0]) != -0: sol = [float(sol_[0]),float(sol_[1]),float(sol_[2])]
    except:
        pass
        
    try:        
        abs_flux_scale_ = params["abs_flux_scale"].strip('][').split(',')
        for j in abs_flux_scale_:
            abs_flux_scale.append(float(j))
    except:
        pass


    # =============================================================================

    # =============================================================================
    #  flux_calibration_file and telluric_correction_file ?   
    # =============================================================================
    if flux_calibration_file == "":flux_calibration_file= os.path.join(path, "flux_calibration_file_"+date+"_"+grating+"_auto.dat")
    if telluric_correction_file == "":telluric_correction_file= os.path.join(path,"telluric_correction_file_"+date+"_"+grating+"_auto.dat")
    # =============================================================================

    # =============================================================================
    # abs_flux_scale ?     
    # =============================================================================
    if len(abs_flux_scale) == 0:
        for i in range(len(CONFIG_FILE_list)): abs_flux_scale.append(1.)
    # =============================================================================


    # =============================================================================
    # Print the summary of parameters readed
    # =============================================================================

    print("> Parameters for automatically processing the calibrations of the night:\n")
    print("  date                       = ",date)
    print("  grating                    = ",grating)
    print("  path                       = ",path)
    if cal_from_calibrated_starcubes == False:
        if do_skyflat:     
            print("  file_skyflat               = ",file_skyflat)
            if skyflat_variable != "" : print("  Python object with skyflat = ",skyflat_variable)
            print("  correct_ccd_defects        = ",correct_ccd_defects)            
            if fix_wavelengths:
                print("  fix_wavelengths            = ",fix_wavelengths)     
                if sol[0] != 0 and sol[0] != -1:
                    print("    sol                      = ",sol)
                else:
                    if rss_star_file_for_sol =="" :
                        print("    ---> However, no solution given! Setting fix_wavelength = False !")
                        fix_wavelengths = False
                    else:
                        print("    Star RSS file for getting small wavelength solution:\n   ",rss_star_file_for_sol)
        else:
            print("  throughput_2D_file         = ",throughput_2D_file)
            if throughput_2D_variable != "" : print("  throughput_2D variable     = ",throughput_2D_variable)
    
        print("  pixel_size                 = ", pixel_size )
        print("  kernel_size                = ", kernel_size )
    
        if len(CONFIG_FILE_list[0]) > 0:
    
            for config_file in range(len(CONFIG_FILE_list)):
                if config_file == 0 : 
                    if len(CONFIG_FILE_list) > 1:       
                        print("  CONFIG_FILE_LIST           =  [",CONFIG_FILE_list[config_file],",")
                    else:
                        print("  CONFIG_FILE_LIST           =  [",CONFIG_FILE_list[config_file],"]")
                else:
                    if config_file == len(CONFIG_FILE_list)-1:
                        print("                                 ",CONFIG_FILE_list[config_file]," ]")
                    else:        
                        print("                                 ",CONFIG_FILE_list[config_file],",")      
        else:
            print("  No configuration file list provided, using generic configuration file for all stars")

    else:
        print("\n> The calibration of the night will be obtained using these fully calibrated starcubes:\n")

    if len(objects_auto) != 0 :
        pprint = ""
        for i in range(len(objects_auto)):
            pprint=pprint+objects_auto[i]+ "  " 
        print("  Using stars in objects     = ",pprint)

    if len(abs_flux_scale) > 0 : print("  abs_flux_scale             = ",abs_flux_scale)
    print("  plot                       = ",plot)
    
    print("\n> Output files:\n")
    if do_skyflat:
        if throughput_2D_variable != "" : print("  throughput_2D variable     = ",throughput_2D_variable)   
        print("  throughput_2D_file         = ",throughput_2D_file)
    print("  flux_calibration_file      = ",flux_calibration_file)
    if grating in red_gratings:
        print("  telluric_correction_file   = ",telluric_correction_file)

    print("\n# =============================================================================")
    

    # =============================================================================
    # do_skyflat            
    # =============================================================================
    if do_skyflat:      
        if rss_star_file_for_sol != "" and sol[0] in [0,-1] :
            print("\n> Getting the small wavelength solution, sol, using star RSS file")
            print(" ",rss_star_file_for_sol,"...")                                  
            if grating in red_gratings :
                _rss_star_ = RSS(rss_star_file_for_sol, instrument=instrument)
                _rss_star_.process_rss(correct_ccd_defects = False, 
                                       fix_wavelengths=True, sol = sol,
                                       plot= plot, plot_final_rss = plot)
            if grating in ["580V"] :
                _rss_star_ = RSS(rss_star_file_for_sol, instrument=instrument)
                _rss_star_.process_rss(rss_star_file_for_sol, 
                                       correct_ccd_defects = True, remove_5577 = True,
                                       plot= plot, plot_final_rss = plot)               
            sol = _rss_star_.sol
            print("\n> Solution for the small wavelength variations:")
            print(" ",sol)
        
        throughput_2D_, skyflat_ =  get_throughput_2D(file_skyflat,
                                                      instrument =instrument,
                                                      plot = plot, 
                                                      plot_final_rss = plot,
                                                      also_return_skyflat = True,
                                                      correct_ccd_defects = correct_ccd_defects,
                                                      fix_wavelengths = fix_wavelengths, sol = sol,
                                                      throughput_2D_file = throughput_2D_file,
                                                      kernel_throughput = kernel_throughput)      
        
        if throughput_2D_variable != "":
            print("  Saving throughput 2D into Python variable:", throughput_2D_variable)
            throughput_2D_variable = throughput_2D_

        if skyflat_variable != "":
            print("  Saving skyflat into Python variable:", skyflat_variable)
            skyflat_variable = skyflat_

    else:
        if cal_from_calibrated_starcubes == False: print("\n> Skyflat will not be processed! Throughput 2D calibration already provided.\n")
        check_nothing_done = check_nothing_done + 1
    # =============================================================================

    good_CONFIG_FILE_list =[]
    good_star_names =[]
    stars=[]
    
    # =============================================================================
    #  cal_from_calibrated_starcubes ?
    # =============================================================================
    # This is for the case that we have individual star cubes ALREADY calibrated in flux
    if cal_from_calibrated_starcubes: 
        pprint = ""
        stars=[]
        good_star_names=[]
        
        for i in range(len(objects_auto)):
            pprint=pprint+objects_auto[i]+ "  "
            
            try: # This is for a combined cube
                stars.append(obj[objects_auto[i]].combined_cube)
                if grating in red_gratings:
                    telluric_correction_list.append(obj[objects_auto[i]].combined_cube.telluric_correction)
            except Exception: # This is when we read a cube from fits file
                stars.append(obj[objects_auto[i]]) 
                if grating in red_gratings:
                    telluric_correction_list.append(obj[objects_auto[i]].telluric_correction)  
            good_star_names.append(stars[i].object)
                
        print("\n> Fully calibrated star cubes provided :",pprint) 
        good_CONFIG_FILE_list = pprint
    
    
    else: 
        
        if len(CONFIG_FILE_list) > 0:
            for i in range(len(CONFIG_FILE_list)):
                run_star = True                
                if CONFIG_FILE_list[i] != "":    
                    print(objects_auto[i])
                    try:
                        config_property, config_value = read_table(CONFIG_FILE_list[i], ["s", "s"] )
                        if len(CONFIG_FILE_list) != len(objects_auto)  :               
                            for j in range (len(config_property)):
                                if config_property[j] == "obj_name" : running_star = config_value[j] 
                            if i < len(objects_auto) :
                                objects_auto[i] = running_star
                            else:    
                                objects_auto.append(running_star)                   
                        else:
                            running_star = objects_auto[i]
                    except Exception:
                        print("===================================================================================")
                        print("\n> ERROR! config file {} not found!".format(CONFIG_FILE_list[i]))
                        run_star = False
                else:
                    running_star = star_list[i]
                    

                if run_star:
                #try:
                    print("===================================================================================")   
                    
                    if CONFIG_FILE_list[i] == "":
                        print("\n> Running automatically calibration star",running_star, "using the standard configuration_file\n")
                    else:
                        print("\n> Running automatically calibration star",running_star, "in CONFIG_FILE:")
                        print(" ",CONFIG_FILE_list[i],"\n")
                    psol= sol
                    plist = list_of_files_of_stars[i]
                    
                    obj[objects_auto[i]] = run_automatic_star(CONFIG_FILE = CONFIG_FILE_list[i],
                                                              object_auto = objects_auto[i],
                                                              star=star_list[i],
                                                              sol = psol,
                                                              throughput_2D_file = throughput_2D_file,
                                                              rss_list = plist,
                                                              path_star = path, 
                                                              date = date, 
                                                              grating = grating,
                                                              pixel_size = pixel_size,
                                                              kernel_size =  kernel_size,
                                                              rss_clean = rss_clean,
                                                              plot = plot
                                                              )
                    
                    if CONFIG_FILE_list[i] == "":
                        print("\n> Running automatically calibration star using generic configuration file SUCCESSFUL !!\n")
                    else:
                        print("\n> Running automatically calibration star in CONFIG_FILE")
                        print("  ",CONFIG_FILE_list[i]," SUCCESSFUL !!\n")
                    
                    
                    good_CONFIG_FILE_list.append(CONFIG_FILE_list[i])
                    good_star_names.append(running_star)
    
    
                    try: # This is for a combined cube
                        stars.append(obj[objects_auto[i]].combined_cube)      
                        if grating in red_gratings:
                            telluric_correction_list.append(obj[objects_auto[i]].combined_cube.telluric_correction)
                    except Exception: # This is when we read a cube from fits file
                        #stars.append(objects_auto[i])
                        stars.append(obj[objects_auto[i]])     
                        if grating in red_gratings:
                            telluric_correction_list.append(obj[objects_auto[i]].telluric_correction)
                            
    # =============================================================================
            
  
    # =============================================================================
    #  CHECK AND GET THE FLUX CALIBRATION FOR THE NIGHT
    # =============================================================================

    if len(good_CONFIG_FILE_list) > 0:        
        
        print("===================================================================================")  
        if grating in red_gratings:
            print("\n> Obtaining the flux calibration and the telluric correction...")  
        else:
            print("\n> Obtaining the flux calibration...")  
        
        # Define in "stars" the cubes we are using, and plotting their responses to check  
        plot_response(stars, scale=abs_flux_scale, plot=plot)

        # We obtain the flux calibration applying:    
        flux_calibration_night = obtain_flux_calibration(stars)
        exec(flux_calibration_name + '= flux_calibration_night', globals())
        print("  Flux calibration saved in variable:", flux_calibration_name)
    
        # And we save this absolute flux calibration as a text file
        w= stars[0].wavelength
        spectrum_to_text_file(w, flux_calibration_night, filename=flux_calibration_file)

        # Similarly, provide a list with the telluric corrections and apply:            
        if grating in red_gratings:
            telluric_correction_night = obtain_telluric_correction(w,telluric_correction_list, label_stars=good_star_names, scale=abs_flux_scale, plot=plot)            
            telluric_correction_name = telluric_correction_night
            print("  Telluric calibration saved in variable:", telluric_correction_name)
    
 # Save this telluric correction to a file
            spectrum_to_text_file(w, telluric_correction_night, filename=telluric_correction_file)

    else:
        print("\n> No configuration files for stars available !")
        check_nothing_done = check_nothing_done + 1
    # =============================================================================

    # =============================================================================
    #  Print Summary
    # =============================================================================
    print("\n# =============================================================================")   
    if CALIBRATION_NIGHT_FILE :
        print("\n> SUMMARY of running configuration file", CALIBRATION_NIGHT_FILE,":\n") 
    else:
        print("\n> SUMMARY of running automatic_calibration_night() :\n") 


    if len(objects_auto) > 0 and cal_from_calibrated_starcubes == False:    
        pprint = ""
        for i in range(len(objects_auto)):
            pprint=pprint+objects_auto[i]+ "  " 
        print("  Created objects for calibration stars           :",pprint) 
    
        if len(CONFIG_FILE_list) > 0:    
            print("  Variable with the flux calibration              :",flux_calibration_name)
            if grating in red_gratings:
                print("  Variable with the telluric calibration          :",telluric_correction_name)
                print(" ")
        print("  throughput_2D_file        = ",throughput_2D_file)
        if throughput_2D_variable != "" : print("  throughput_2D variable    = ",throughput_2D_variable)
    
        if sol[0] != -1 and sol[0] != 0:
            print("  The throughput_2D information HAS BEEN corrected for small wavelength variations:")
            print("  sol                       =  ["+np.str(sol[0])+","+np.str(sol[1])+","+np.str(sol[2])+"]")
    
        if skyflat_variable != "" : print("  Python object created with skyflat = ",skyflat_variable)
        
        if len(CONFIG_FILE_list) > 0:
            print('  flux_calibration_file     = "'+flux_calibration_file+'"')
            if grating in red_gratings:
                print('  telluric_correction_file  = "'+telluric_correction_file+'"')
 
    if cal_from_calibrated_starcubes:
        print("  Variable with the flux calibration              :",flux_calibration_name)
        if grating in red_gratings:
                print("  Variable with the telluric calibration          :",telluric_correction_name)
                print(" ")       
        print('  flux_calibration_file     = "'+flux_calibration_file+'"')
        if grating in red_gratings:
            print('  telluric_correction_file  = "'+telluric_correction_file+'"')
        
 
    if check_nothing_done == 2:
        print("\n> NOTHING DONE!")
              
    print("\n# =============================================================================")
    # =============================================================================

    return obj





