#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a script to calculate and plot the bound fraction of a MAP
from tubulin co-sedimentation assay and extrapolate the Kd value.
Comparison between multiple conditions is enabled.

Author: Giulia Tonon
ZMBH Heidelberg, AG Pfeffer 
"""

import argparse
from operator import index
import sys
import os
import pathlib
import glob
import re
import inspect
from turtle import color
try:
    import matplotlib as mpl
    #mpl.use('pgf')   # pfg backend to export plots in a format that can be embedded in LaTeX
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from scipy.optimize import curve_fit

except:
	ImportError
	print("Please, make sure to have matplotlib, numpy, pandas and scipy installed. You might also need openpyxl")
	print("You might need to activate a conda environment (e.g. plotting)")

##########################################################################################

### DEFINE FUNCTIONS ###

def dir_path(string):
    """
    check that the passed folder does exist
    from https://stackoverflow.com/questions/38834378/path-to-a-directory-as-argparse-argument
    """
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def does_file_exist(f):
    """
    check that the input file exists
    This function is an extra-check but should not be necessary
    as argparse checks whether the input files are readable.
    """
    if not os.path.isfile(f):
        print(f + ' could not be found')
        print(f"Currently you are in' + {os.getcwd()}")
        sys.exit()
    else:
        return True


def is_xlsx(f):
    """
    Check that the file is a .xlsx
    """
    if not f.endswith('.xlsx'):
        print(f + ' has the wrong file extension.')
        sys.exit()
    else:
        return True


def select_out_folder(list_input_dir):
    """
    Function to select a folder where to save the output
    """
    cwd = os.getcwd()
    if len(list_input_dir) != 1:
        path_user = input(f"Please, enter the path where you want to save the output. Currently you are in {cwd}: ")
        if dir_path(str(path_user)):
            out_path = pathlib.Path(dir_path(str(path_user)))

    elif len(list_input_dir) == 1:
        out_path = pathlib.Path(list_input_dir[0])
    
    return out_path


def band_selection(band_list, reg):
    """
    Function to select which band to use for the plot
    """
    band_list = [reg.sub('',quantification_band) for quantification_band in band_list]  # keep only the band name for each quantification folder
    band_list = list(dict.fromkeys(band_list))   # remove duplicates
    if len(band_list) > 1:
        band = input("Please select which band you wish to use for your plot: ")   # ask the user to select a band to use for the plot for that construct
        while band not in band_list:
            band = input(f"Choose one of {','.join(band_list)}: ")
    else:
        band = band_list[0]
    return band


def get_filepaths_with_oslistdir(root_path: str, file_regex: str):
    """
    Function to use regex for matching a desidered folder
    edited from https://stackoverflow.com/questions/39293968/how-do-i-search-directories-and-find-files-that-match-regex
    """
    folder_paths = []
    pattern = re.compile(file_regex, re.IGNORECASE)
    list_dirs = os.listdir(root_path)
    for e_dir in list_dirs:
        if pattern.match(e_dir):
            folder_paths.append(e_dir)
    return folder_paths


def get_construct_def(stringPath):
    """
    Function to fetch construct specifications from folder names
    """
    e_path = pathlib.Path(stringPath)    # transform the string of the input sample folder into a path object
    os.chdir(e_path)
    construct_path = os.getcwd()
    parts_construct_path = os.path.normpath(construct_path).split(os.path.sep)
    construct_part1 = parts_construct_path[-2]
    construct_part2 = parts_construct_path[-1]
    print(f"Construct being analysed: {parts_construct_path[-2]} {parts_construct_path[-1]}")
    construct_def = (construct_part1, construct_part2)

    return construct_def


def get_specifications(pandas_xlsx, construct_defs):
    """
    Function to get construct name, colour for rendering and tubulin range 
    used in the experiment from excel file
    """
    construct_specifications = pandas_xlsx.get(construct_defs[0])
    construct_name = construct_specifications.loc['name',construct_defs[1]]
    tubulin_string = construct_specifications.loc['tubulin_conc', construct_defs[1]]
    tubulin_float_list = [float(tub) for tub in tubulin_string.split()]
    colour = construct_specifications.loc['colour',construct_defs[1]]

    return construct_name, colour, tubulin_float_list


def check_replica_quantification_folder(tubulin_vector):
    """
    Function to check how many replicas are present for the construct
    being analysed and how many folders with quantification files for each replica
    """
    tot_replicas = get_filepaths_with_oslistdir('.', r'.*[0-9]$')    #^(replica)|(#)')   # get a list of all the replica folders within the sample folder
    print(tot_replicas)
    
    if len(tot_replicas):

        empty_single_intensities_array = np.empty([len(tot_replicas),2*len(tubulin_vector)]) #initialize an empty vector to store the band intensity values for each replica
        quantification_bands = []   # initialize an empty list to store all possibilities of quantified bands performed for the construct being analysed
        quantification_directories = []   # intialize a list to store all directories containing the quantification files from Fiji for the construct under analysis

        for replica in tot_replicas:

            os.chdir(replica)
            quantification_directories.append(os.getcwd())
            quantification_folders = get_filepaths_with_oslistdir('.', r'^(quantification).*')
            #quantification_folders = glob.glob("quantification*")

            if len(quantification_folders):
                for quantification_folder in quantification_folders:
                    quantification_bands.append(quantification_folder)
                    
            else:
                print('It looks like you did not perform band quantification in Fiji for '+replica+'. Exiting now...')
                sys.exit()
                
            os.chdir('..')
            
        plain_quant_regex = re.compile(r'^(quantification)$', re.IGNORECASE)
        plain_quant_folders = [quant_band for quant_band in quantification_bands if plain_quant_regex.match(quant_band)]   # list of quantification folders that do not have an indicated band
        
        if len(plain_quant_folders):
            quant_bands_not_plain = [quant_band for quant_band in quantification_bands if quant_band not in plain_quant_folders]   # list of quantification folders that do have an indicated band
            if len(quant_bands_not_plain):
                print("It appears that you named the quantification folders in the wrong way")
                sys.exit()
            else:
                quantification_folders_use = quantification_bands

        else:
            pattern = re.compile(r'^(quantification)_', re.IGNORECASE)
            selected_band = band_selection(quantification_bands, pattern)   # select band to use for the plot
            quantification_folders_use = [quant_band for quant_band in quantification_bands if re.match(r'.*'+selected_band+r'$',quant_band)]   # retain quantification folders corresponding to teh selected band

        if len(quantification_folders_use) == len(quantification_directories):   # check that you are using only one quantiifcation folder for each replica
            quantification_directories = [os.path.join(quantification_directories[i], quantification_folders_use[i]) for i in range(len(quantification_directories))]   # add quantification folder to the single replica path
        else:
            print("It seems that the number of replica folders and of selected quantification folders do not match")
            sys.exit()
    
    else:
        print("It appears as if you don't have any replica for the selected construct. Make sure you select a valid construct. Exiting now...")
        sys.exit()
        
        
    return empty_single_intensities_array, quantification_directories


def values_wrangler(list_csv):
    """
    Function to get the band intensity values 
    from list of csv-files saved from Fiji,
    fill them with values corresponding to 0 
    according to the file name and sort them 
    according to the corresponding tubulin range
    """
    pos0_regex = re.compile(r'0pos(.*)_', re.IGNORECASE)
    low_regex = re.compile(r'.*_low', re.IGNORECASE)
    high_regex = re.compile(r'.*_high', re.IGNORECASE)
    
    if len(list_csv) == 2:
        for csv_file in list_csv:
        
            csvFile = pd.read_csv(csv_file, usecols = ["Area"])   # read csv-file
            values = csvFile.values[:,0]                          # numpy array with the band intensity values
            
            pos0_match = re.search(pos0_regex, csv_file)
            if pos0_match:
                pos0_match_string = str(pos0_match.group(1))
                if not pos0_match_string.startswith('-'):
                    pos0_list = pos0_match_string.split('-')        # get a list of the indexes of 0 values form the csv-file name
                    pos0_list = [int(pos0) for pos0 in pos0_list]   # convert the indexes into integers
                    for pos in pos0_list:
                        values = np.insert(values, pos-1, 0.0)      # add 0.0 where required in the values numpy.array obj
                
                # sort band intensity values according to the tubulin range (indicated in the file name)
                if re.search(low_regex, csv_file):
                    low_values = values
                elif re.search(high_regex, csv_file):
                    high_values = values
            
            else:
                raise ValueError("The file name is invalid.")
                sys.exit()
    
    else:
        print(f"{len(list_csv)} csv-files found. Make sure that you have one for low range of tubulin concentrations and one for high tubulin concentrations.")
        sys.exit()
    
    return low_values, high_values


def handle_csv(list_dirs, empty_np_array):
    """
    Function to read the cvs-files from the band quantification
    performed in Fiji for each construct and fill the array to 
    store the values for further manipulation
    """
    for e_dir in list_dirs:
    
        os.chdir(e_dir)
        print('Quantification folder being analysed: '+os.getcwd())
        csv_files = glob.glob("*.csv")   # get a list with all csv-files for the current folder

        if len(csv_files):
            list_values_array = values_wrangler(csv_files)   # get band intensity values sorted according to the tubulin range
            if len(list_values_array[0]) == len(list_values_array[1]) == empty_np_array.shape[1]/2:   # check that the length of the two arrays correspond to the tubulin range
                intensity_values = np.concatenate([list_values_array[0],list_values_array[1]])        # concatenate the values into a single band intensity array for each replica
                index = list_dirs.index(e_dir)
                empty_np_array[index] = intensity_values   # fill the values of the array storing the values for each replica

            else:
                print("The number of intensities does not match the tubulin range. Exiting now...")
                sys.exit()

        else:
            print("No quantification files detected. Exiting now...")
            sys.exit()

    return empty_np_array


def calculate_fraction_bound(replica_intensity_array):
    """
    Function to calculate the fraction of your MAP of interest
    bound to MTs based on the co-pelleted fraction for each replica
    of one construct
    
    """
    # Initialize empty numpy arrays to store the values for each replica
    replica_intensity_shape = np.shape(replica_intensity_array)
    intensity_pellet = np.empty([replica_intensity_shape[0], int(replica_intensity_shape[1]/2)])       # pellet intensity
    intensity_supernatant = np.empty([replica_intensity_shape[0], int(replica_intensity_shape[1]/2)])  # supernatant intensity
    intensity_tot = np.empty([replica_intensity_shape[0], int(replica_intensity_shape[1]/2)])          # total intensity (total MAP)
    frac_bound = np.empty([replica_intensity_shape[0], int(replica_intensity_shape[1]/2)])             # fraction of MAP bound to (=> co-pelleted with) MTs
    
    for i in range(0, replica_intensity_shape[0]):                                                   # loop through the arrays representing the intensity values for each replica
        intensity_pellet[i] = replica_intensity_array[i][1:len(replica_intensity_array[i]):2]        # select intensity values of pellet fractions (every 2nd element)
        intensity_supernatant[i] = replica_intensity_array[i][0:len(replica_intensity_array[i])-1:2] # select intensity values of supernatant (=> unbound) fractions
        if np.shape(intensity_pellet[i]) == np.shape(intensity_supernatant[i]):
            intensity_tot[i] = np.add(intensity_pellet[i], intensity_supernatant[i])                 # intensity of total protein
            frac_bound[i] = np.divide(intensity_pellet[i],intensity_tot[i])                          # fraction bound for that replica
        else:
            print("An error occurred. Exiting now...")
            sys.exit()

    return intensity_pellet, intensity_supernatant, intensity_tot, frac_bound


def calculate_mean_sttdev(np_array_2d):
    """
    Function to calculate the mean value of MAP bound to MTs
    across the replicas and the corresponding standard deviation
    """
    average_stddev = np.empty([2,np.shape(np_array_2d)[1]])   # initialize an empty np array to store the mean values (1st array) and the stddev values (2nd array)
    average_stddev[0] = np.mean(np_array_2d, axis = 0)        # calculate mean values for bound fraction across replicas
    average_stddev[1] = np.std(np_array_2d, axis = 0)         # calculate standard deviation across replicas for bound fraction
    return average_stddev


def quadr_eq(x,B_t,K_d):
    """
    Define the modified quadratic equation to fit
    the data and extrapolate the dissociation constant (K_d)
    Based on Ciferri et al (2008, Cell)
    """
    y = ((B_t + K_d + x) - np.sqrt((B_t+K_d+x)**2-(4*B_t*x)))/2
    return y


def fit_model(construct_list, x_data, y_data, eq):
    """
    Function to fit the model specified by eq to the experimental
    data and returns the optimal parameters with their standard
    errors, R-squared and the theoretical values given the values 
    of the independent variable
    """
    list_args_eq = inspect.signature(eq).parameters.values()   # get arguments of eq
    nr_params = len(list_args_eq) -1   # get how many parameters have to be optimized (accounts for just one independent variable)
    fit_results = {}                   # create an empty dictionary to store the results
    x_data_extremes = [func(l) for l in x_data for func in (np.min,np.max)]
    x_model = np.linspace(min(x_data_extremes),max(x_data_extremes),num = 1000)   # x-axis values for computing theoretical y-values
    
    # initialize empty arrays to store results
    opt_params = np.empty([len(construct_list),nr_params])   # for model parameters
    perr = np.empty([len(construct_list),nr_params])         # for standard errors of model parameters
    r_squared = np.empty(len(construct_list))                # R-squared
    y_model = np.empty([len(construct_list),len(x_model)])   # values of dependent variable according to the model (optimal y-values given our x-values)
    
    # fit model
    for construct in construct_list:
        index = construct_list.index(construct)
        y_data_i_array = np.asarray(y_data[index])   # convert "experimental" y_data into a format for the fit
        fit_i = curve_fit(eq,x_data[index],y_data_i_array, bounds = ([0,0],[1,np.inf]))  # fit model to the data, !!! bounds = ([0,0],[1,np.inf] )!!!
        opt_params[index] = fit_i[0]                 # optimal parameters for the construct being analysed
        
        # calculate the R-squared for the fit and the standard errors for the parameters
        # based on https://stackoverflow.com/questions/19189362/getting-the-r-squared-value-using-curve-fit
        residuals_i = y_data_i_array - eq(x_data[index],*fit_i[0])
        ss_res_i = np.sum(residuals_i**2)                                 # residuals sum of squares
        ss_tot_i = np.sum((y_data_i_array - np.mean(y_data_i_array))**2)  # total sum of squares
        r_squared[index] = 1 - (ss_res_i/ss_tot_i)                        # R-squared value
        perr[index] = np.sqrt(np.diag(fit_i[1]))                          # standard deviation errors on the parameters
        y_model[index] = eq(x_model, *fit_i[0])                           # theoretical y-values
    
    # Update and return results
    fit_results["Kd_[uM]"] = opt_params[:,1]
    fit_results["Kd_stddev"] = perr[:,1]
    fit_results["Bt (maximal MAP-tubulin complex)"] = opt_params[:,0]
    fit_results["Bt_stddev"] = perr[:,0]
    fit_results["r_squared"] = r_squared
    return x_model, y_model, fit_results


##########################################################################################

parser = argparse.ArgumentParser(description = 'This is a script to plot the bound fraction of a MAP \
                                assessed via tubulin co-sedimentation assay.\
                                Requirements: matplotlib, numpy, pandas and scipy. You might also need openpyxl')

parser.add_argument('--samplePath', nargs = '+', type = dir_path, required = True,
                    help = 'Folder with file from band intensity quantification carried out in Fiji')

parser.add_argument('--Kd', action = 'store_true',
                    help = 'Print the Kd values')

parser.add_argument('--Kd_round_digits', type = int, default = 1,
                    help = 'Number of digits to round the Kd')

parser.add_argument('--ext', type = str, choices = ['svg','png','pdf','tiff', 'jpg', 'jpeg'], default = 'png',
                    help = 'File extension for output plot')

parser.add_argument('--pathToExcel', type = argparse.FileType('r'), required = True, 
                    help = 'Path to Excel file containing one worksheet per construct,\
                        and each worksheet contains one column for each construct variant\
                            -> first row contains the string for the color, second row the string for the construct name for the plots,\
                            third row the tubulin concentrations used')

args = parser.parse_args()

##########################################################################################

start = os.getcwd()
out_folder = select_out_folder(args.samplePath)   # select the folder where to save the output

# colors and construct names from Excel input file
if does_file_exist(args.pathToExcel.name):
    if is_xlsx(args.pathToExcel.name):
        colors_constructs=pd.read_excel(open(args.pathToExcel.name, 'rb'), sheet_name=None,index_col=0, header=0)


########################################################################################################################

### Calculate data to plot ###

# Initialize empty arrays or lists to store relevant info for each construct
mean_across_replicas = []   # to store the mean value of fraction bound
stddev_across_replicas = []   # to store the standard deviation of fraction bound
construct_names = []    # to store the name of the constructs
color = []              # to store the colors
tubulin_conc = []       # list of lists to store tubulin (independent variable)


for construct in args.samplePath:

    index_construct = args.samplePath.index(construct)
    
    # Retrieve specifications for construct
    construct_parts = get_construct_def(construct)   # get the elements that define the construct under analysis from name
    construct_spec = get_specifications(colors_constructs, construct_parts)
    
    construct_names.append(construct_spec[0])   # construct variant name
    color.append(construct_spec[1])             # color for the plot

    tubulin_conc_i = construct_spec[2]
    tubulin_conc.append(tubulin_conc_i)   # list of lists
    
    quantification_folders = check_replica_quantification_folder(tubulin_conc_i)
    replicas_intensity_values = handle_csv(quantification_folders[1], quantification_folders[0])

    calculation_intensitites = calculate_fraction_bound(replicas_intensity_values)
    mean_stddev_across_replicas = calculate_mean_sttdev(calculation_intensitites[3])

    mean_across_replicas.append(mean_stddev_across_replicas[0])
    stddev_across_replicas.append(mean_stddev_across_replicas[1]) 
    
    os.chdir(start)


tubulin_conc_arrays = [np.asarray(tub_list) for tub_list in tubulin_conc]
os.chdir(out_folder)   # change folder inside directory in which to save the plot


### Plot the data ### 
# and save output in PGF format for importing in LaTeX (https://jwalton.info/Matplotlib-latex-PGF/) --> NOT ENABLED!

plt.rcParams['font.family'] = 'sans-serif'   # set the font type
fig, ax = plt.subplots()

# Set axis labels
ax.set_xlabel('Tubulin ['+r'$\mu$'+'M]', fontsize = '16')
ax.set_ylabel('Fraction bound', fontsize = '16')

# Adjust axis margins and limits
tubulin_up_lim = [max(tubulin) for tubulin in tubulin_conc]
ax.set_ylim(0,1.05)                       # set y axis range
ax.set_xlim(0,max(tubulin_up_lim)+0.5)   # set the x-axis range
plt.margins(0)     # resert the offset of the axis ot 0
ax.minorticks_on() # switch on minor ticks
ax.tick_params(axis='both', which='major', labelsize=12)   # set fontsize of major ticks for both axes

# Show only left and bottom frame
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(True)
ax.spines['bottom'].set_visible(True)

# Fit quadratic equation to data
fit = fit_model(args.samplePath, tubulin_conc_arrays, mean_across_replicas, quadr_eq)

df_fit = pd.DataFrame(fit[2])                   # create a dataframe to export the results of the fit
df_fit.insert(0,"Construts", construct_names)   # add construct identifiers

title_output = ('_').join(construct_names)   # base name for output files
spec_char = re.compile(r'(\W+)|((?<=\\)it)')   # regex to match all groups of non-alphanumeric characters and 'it' used to indicate italic
title_output = spec_char.sub('',title_output)   # cleaned base name for output files

df_fit.to_csv(title_output+'_fit.csv',index = False)   # export the results of the fit as csv-file


# define the size of the marker (both plot and legend), linewidth of the fit and opacity
marker_size = 15
line_width = 3
opacity = 1
marker_legend = 3

opacity_array = np.repeat(opacity, len(args.samplePath))   # initialize an array with opacity values for each construct
for construct in args.samplePath:
    index_construct = args.samplePath.index(construct)
    construct_id = construct_names[index_construct]
     
    
    Kd_construct = df_fit.loc[index_construct,'Kd_[uM]']
    r_squared = df_fit.loc[index_construct,'r_squared']
    if r_squared < 0.1:
        Kd_label_parts = ('K'+r'$_{d}$', '-')   # if R-squared < 0.1, do not indicate Kd in the legend
    else:
        Kd_label_parts = ('K'+r'$_{d}$', r'$\approx$', str(round(Kd_construct,args.Kd_round_digits)), r'$\mu$M')   
    
    # define the text of the legend
    if args.Kd:
        Kd_label = (' ').join(Kd_label_parts)
        label_text_parts = (construct_id, Kd_label)
        label_text = (', ').join(label_text_parts)
    else:
        label_text = construct_id
    
    plt.scatter(tubulin_conc_arrays[index_construct], mean_across_replicas[index_construct], marker = '.', s = marker_size, c = color[index_construct], 
                label = label_text, alpha=opacity_array[index_construct])
    plt.errorbar(tubulin_conc_arrays[index_construct],mean_across_replicas[index_construct], yerr = stddev_across_replicas[index_construct],
                c = color[index_construct], elinewidth = 0.5, capsize = 3, fmt = '.', alpha=opacity_array[index_construct])
    plt.plot(fit[0], fit[1][index_construct], c = color[index_construct], linewidth = line_width, alpha = opacity_array[index_construct])   # plot the fitted model
    
# Set legend
legend = ax.legend(loc = 'lower right', frameon = False, markerscale = marker_legend, fontsize = 'large')   # position the legend at the lower right and switch off the legend box

# Save figure
title_fig = title_output+'_pelleting_plot.'+args.ext
fig.savefig(title_fig, bbox_inches = 'tight', dpi = 700, transparent = True)


os.chdir(start)