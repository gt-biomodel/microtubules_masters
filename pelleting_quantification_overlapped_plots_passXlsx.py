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


def band_selection(band_list):
    """
    Function to select which band to use for the plot
    """
    if len(band_list) > 1:
        band_list = [quantification_band.lstrip("quantification_") for quantification_band in band_list]  # keep only the band name for each quantification folder
        band = input("Please select which band you wish to use for your plot: ")   # ask the user to select a band to use for the plot for that construct
        while band not in band_list:
            band = input(f"Choose one of {','.join(band_list)}: ")
        band = "quantification_"+band
    else:
        band = band_list[0]
    return band


def check_replica_quantification_folder(stringPath, tubulin_vector):
    """
    Function to check how many replicas are present for the construct
    being analysed and how many folders with quantification files for each replica
    """
    e_path = pathlib.Path(stringPath)    # transform the string of the input sample folder into a path object
    os.chdir(e_path)
    construct_path = os.getcwd()
    parts_construct_path = os.path.normpath(construct_path).split(os.path.sep)
    construct_part1 = parts_construct_path[-2]
    construct_part2 = parts_construct_path[-1]
    construct_def = (construct_part1, construct_part2)
    print(f"Construct being analysed: {parts_construct_path[-2]} {parts_construct_path[-1]}")
    tot_replicas = glob.glob("replica*") # get a list of all the replica folders within the sample folder
    
    if len(tot_replicas):

        empty_single_intensities_array = np.empty([len(tot_replicas),2*len(tubulin_vector)]) #initialize an empty vector to store the band intensity values for each replica
        quantification_bands = []   # initialize an empty list to store all possibilities of quantified bands performed for the construct being analysed
        quantification_directories = []   # intialize a list to store all directories containing the quantification files from Fiji for the construct under analysis

        for replica in tot_replicas:

            os.chdir(replica)
            quantification_directories.append(os.getcwd())
            quantification_folders = glob.glob("quantification_*")

            if len(quantification_folders):
                for quantification_folder in quantification_folders:
                    quantification_bands.append(quantification_folder)
                    
            else:
                print('It looks like you did not perform band quantification in Fiji for '+replica+'. Exiting now...')
                sys.exit()
                
            os.chdir('..')
            
        quantification_bands = list(dict.fromkeys(quantification_bands))   # remove the duplicates for band quantification choices
        selected_band = band_selection(quantification_bands)
        quantification_directories = [os.path.join(quantification_dir, selected_band) for quantification_dir in quantification_directories]   # add quantification folder to the single replica path
    
    else:
        print("It appears as if you don't have any replica for the selected construct. Make sure you select a valid construct. Exiting now...")
        sys.exit()
        
        
    return empty_single_intensities_array, quantification_directories, construct_def


def values_wrangler(list_csv):
    """
    Function to get the band intensity values 
    from list of csv-files saved from Fiji,
    fill them with values corresponding to 0 
    according to the file name and sort them 
    according to the corresponding tubulin range
    """
    pos0_regex = re.compile(r'0pos(.*)_')
    low_regex = re.compile(r'.*_low')
    high_regex = re.compile(r'.*_high')
    
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
    It accounts for x_data being a 1D array and being the same 
    for all constructs!!!
    """
    list_args_eq = inspect.signature(eq).parameters.values()   # get arguments of eq
    nr_params = len(list_args_eq) -1   # get how many parameters have to be optimized (accounts for just one independent variable)
    fit_results = {}                   # create an empty dictionary to store the results
    x_model = np.linspace(min(x_data),max(x_data),num = 1000)   # x-axis values for computing theoretical y-values
    
    # initialize empty arrays to store results
    opt_params = np.empty([len(construct_list),nr_params])   # for model parameters
    perr = np.empty([len(construct_list),nr_params])         # for standard errors of model parameters
    r_squared = np.empty(len(construct_list))                # R-squared
    y_model = np.empty([len(construct_list),len(x_model)])   # values of dependent variable according to the model (optimal y-values given our x-values)
    
    # fit model
    for construct in construct_list:
        index = construct_list.index(construct)
        y_data_i_array = np.asarray(y_data[index])   # convert "experimental" y_data into a format for the fit
        fit_i = curve_fit(eq,x_data,y_data_i_array)  # fit model to the data
        opt_params[index] = fit_i[0]                 # optimal parameters for the construct being analysed
        
        # calculate the R-squared for the fit and the standard errors for the parameters
        # based on https://stackoverflow.com/questions/19189362/getting-the-r-squared-value-using-curve-fit
        residuals_i = y_data_i_array - eq(x_data,*fit_i[0])
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
                                Requirements: matplotlib, numpy, pandas and scipy')

parser.add_argument('--samplePath', nargs = '+', type = dir_path, required = True,
                    help = 'Folder with file from band intensity quantification carried out in Fiji')

parser.add_argument('--Kd', action = 'store_true',
                    help = 'Print the Kd values')

parser.add_argument('--ext', type = str, choices = ['svg','png','pdf','tiff', 'jpg', 'jpeg'], default = 'png',
                    help = 'File extension for output plot')

parser.add_argument('--pathToExcel', type = argparse.FileType('r'), required = True, 
                    help = 'Path to Excel file containing one worksheet per construct,\
                        and each worksheet contains one column for each construct variant\
                            -> first row contains the string for the color, second row the string for the construct name for the plots')

args = parser.parse_args()

##########################################################################################

tubulin_conc = [0, 0.25, 0.75, 1.25, 1.75, 2.5, 5.0, 10.0]  # tubulin concentrations used in the experiment
tubulin_conc_array = np.asarray(tubulin_conc)               # "experimental" x values (independent variable) - required for fitting model


start = os.getcwd()
out_folder = select_out_folder(args.samplePath)   # select the folder where to save the output

# colors and construct names from Excel input file
if does_file_exist(args.pathToExcel.name):
    if is_xlsx(args.pathToExcel.name):
        colors_constructs=pd.read_excel(open(args.pathToExcel.name, 'rb'), sheet_name=None,index_col=0, header=0)


########################################################################################################################

### Calculate data to plot ###

# Initialize empty arrays or lists to store relevant info for each construct
mean_across_replicas = np.empty([len(args.samplePath),len(tubulin_conc)])   # to store the mean value of fraction bound
stddev_across_replicas = np.empty([len(args.samplePath),len(tubulin_conc)]) # to store the standard deviation of fraction bound
color = []              # to store the colors
construct_names = []    # to store the name of the constructs


for construct in args.samplePath:

    index_construct = args.samplePath.index(construct)
    
    quantification_folders = check_replica_quantification_folder(construct, tubulin_conc)
    construct_parts = quantification_folders[2]   # get the elements that define the construct under analysis
    replicas_intensity_values = handle_csv(quantification_folders[1], quantification_folders[0])

    calculation_intensitites = calculate_fraction_bound(replicas_intensity_values)
    mean_stddev_across_replicas = calculate_mean_sttdev(calculation_intensitites[3])

    mean_across_replicas[index_construct] = mean_stddev_across_replicas[0]
    stddev_across_replicas[index_construct] = mean_stddev_across_replicas[1]

    # Retrieve specifications for construct rendering
    spec = colors_constructs.get(construct_parts[0])   # get the data frame object corresponding to the construct of interest
    color.append(spec.loc['colour',construct_parts[1]])            # color
    construct_names.append(spec.loc['name', construct_parts[1]])   # construct variant name
    
    os.chdir(start)


os.chdir(out_folder)   # change folder inside directory in which to save the plot
print(color, construct_names)


### Plot the data ### 
# and save output in PGF format for importing in LaTeX (https://jwalton.info/Matplotlib-latex-PGF/) --> NOT ENABLED!

plt.rcParams['font.family'] = 'sans-serif'   # set the font type
fig, ax = plt.subplots()

# Set axis labels
ax.set_xlabel('Tubulin ['+r'$\mu$'+'M]', fontsize = '16')
ax.set_ylabel('Fraction bound', fontsize = '16')

# Adjust axis margins and limits
ax.set_ylim(0,1.05)                       # set y axis range
ax.set_xlim(0,max(tubulin_conc)+0.5)   # set the x-axis range
plt.margins(0)     # resert the offset of the axis ot 0
ax.minorticks_on() # switch on minor ticks
ax.tick_params(axis='both', which='major', labelsize=12)   # set fontsize of major ticks for both axes

# Show only left and bottom frame
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(True)
ax.spines['bottom'].set_visible(True)

# Fit quadratic equation to data
fit = fit_model(args.samplePath, tubulin_conc_array, mean_across_replicas, quadr_eq)

df_fit = pd.DataFrame(fit[2])                   # create a dataframe to export the results of the fit
df_fit.insert(0,"Construts", construct_names)   # add construct identifiers
df_fit.to_csv(('_').join(construct_names)+'_fit.csv',index = False)   # export the results of the fit as csv-file


# define the size of the marker (both plot and legend), linewidth of the fit and opacity
if len(args.samplePath) >1:
    marker_size = 15
    line_width = 3
    opacity = 1
    marker_legend = 3
else:
    marker_size = 15
    line_width = 3
    opacity = 1
    marker_legend = 3

opacity_array = np.repeat(opacity, len(args.samplePath))   # initialize an array woth opacity values for each construct -> trick to escape opacity for human wt
for construct in args.samplePath:
    index_construct = args.samplePath.index(construct)
    construct_id = construct_names[index_construct]
     
    
    Kd_construct = df_fit.loc[index_construct,'Kd_[uM]']
    r_squared = df_fit.loc[index_construct,'r_squared']
    if r_squared < 0.1:
        Kd_label_parts = ('K'+r'$_{d}$', '-')   # if R-squared < 0.1, do not indicate Kd in the legend
    else:
        Kd_label_parts = ('K'+r'$_{d}$', r'$\approx$', str(round(Kd_construct,1)), r'$\mu$M')   
    
    # define the text of the legend
    if args.Kd:
        Kd_label = (' ').join(Kd_label_parts)
        label_text_parts = (construct_id, Kd_label)
        label_text = (', ').join(label_text_parts)
    else:
        label_text = construct_id
    
    plt.scatter(tubulin_conc_array, mean_across_replicas[index_construct], marker = '.', s = marker_size, c = color[index_construct], 
                label = label_text, alpha=opacity_array[index_construct])
    plt.errorbar(tubulin_conc_array,mean_across_replicas[index_construct], yerr = stddev_across_replicas[index_construct],
                c = color[index_construct], elinewidth = 0.5, capsize = 3, fmt = '.', alpha=opacity_array[index_construct])
    plt.plot(fit[0], fit[1][index_construct], c = color[index_construct], linewidth = line_width, alpha = opacity_array[index_construct])   # plot the fitted model
    
# Set legend
legend = ax.legend(loc = 'lower right', frameon = False, markerscale = marker_legend, fontsize = 'large')   # position the legend at the lower right and switch off the legend box

# Save figure
title_fig = ('_').join(construct_names)+'_pelleting_plot.'+args.ext
fig.savefig(title_fig, bbox_inches = 'tight', dpi = 700, transparent = True)


os.chdir(start)