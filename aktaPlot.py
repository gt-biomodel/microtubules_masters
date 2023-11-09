#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a script to plot the results of AktaGo runs from Cytiva
controlled by Unicorn software v7.6
It can handle gel filtration and ion-exchange chromatography data.

USAGE EXAMPLES in Unix OS: 
    aktaPlot.py sup6increase*.csv sec -f svg -log_path "valid/path"    # output: plot of elution profile from SEC run saved as .svg-file, fractions highlighted, no conductivity
    aktaPlot.py sup6increase*.csv sec -f png -e -log_path "valid/path" # output: plot of elution profile from SEC run saved as .png-file, no fractions marked, no conductivity
    aktaPlot.py sup6increase*.csv sec -f svg -log_path "valid/path" -c # output: plot of elution profile from SEC run saved as .svg-file, fractions marked, conductivity plotted
    aktaPlot.py sup6increase*.csv sec -f pdf -e -log_path "valid/path" # output: plot of elution profile from SEC run saved as PDF, no fractions highlighted
USAGE EXAMPLE in Windows OS:
    python aktaPlot.py aktafile1.csv aktafile2.csv iex -f svg -log_path "valid/path"

INPUT: list of .csv files storing results of AktaGo runs
       Important: If you run the script on Windows Anaconda PowerShell, you need to specify the full filename without wildcards!

OUTPUT: file of specified format with plot of the elution profile (and concentration of high-salt buffer for IEX runs)

Programming language Python::3.9
Requirements: matplotlib, numpy, chardet

Author: Giulia Tonon
ZMBH Heidelberg, AG Pfeffer
"""

import argparse
import sys
import os
import pathlib
import datetime
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    from chardet import detect
except:
    ImportError
    print("Please make sure to have installed matplotlib, numpy, chardet")
    sys.exit()


parser = argparse.ArgumentParser(description = 'This is a script to plot the results from column chromatography runs performed with AktaGo controlled by Unicorn software.\
                                Python 3.9, requirements: matplotlib, numpy, chardet',
                                 epilog = 'Notes: By default, aktaPlot.py will plot data elution profile (UV absorbance) until the last fraction,\
                                 without conductivity and marking the elution fractions.')

parser.add_argument('list_csv_files', metavar='csv', nargs = '+', type = argparse.FileType('r'),
        help = 'list of csv files storing the results of the runs you wish to plot')

parser.add_argument('run_type', metavar='rt', type = str, choices = ['sec', 'iex'],
        help = 'specify the type of run: Either gel filtration (sec) or ion-exchange chromatography (iex)')

parser.add_argument('-f', '--format', type = str, choices = ['svg','png','pdf','tiff', 'jpg', 'jpeg'], required = True,
        help = 'specify the output file format')

parser.add_argument('-e', '--no_elutions', action = 'store_true',
        help = 'Disable plotting of fractions')

parser.add_argument('-c', '--conductivity', action = 'store_true',
        help = 'Enable plotting conductivity')

parser.add_argument('--xmin', type = float, default = None,
        help = 'Lower bound on the x-axis')

parser.add_argument('--xmax', type = float, default = None,
        help = 'Upper bound on the x-axis')

parser.add_argument('--dpi', nargs = '?', default = 300,
        help = 'Set DPI (dots per inch) of the output file for raster images.\
                Default is 300.')    # if no --dpi flag is passed, default value will be used (combination nargs='?' and default)

parser.add_argument('--color', type = str, choices = ['black','grey','blue1','blue2','green1', 'green2'],
        help = 'Choose color for elution profile. Default: black')

parser.add_argument('--no_legend', action = 'store_true',
        help = 'Disable legend')

parser.add_argument('--peak', nargs = 2, default = None,
        help = 'First and last elution fractions to color. Example: "C.10" "D.12"')

parser.add_argument('-lpath', '--log_path', type = pathlib.Path, required = True,
        help = 'Directory in which to write the output files and the log-file')

group_title = parser.add_mutually_exclusive_group()
group_title.add_argument('--no_title', action = 'store_true',
            help = 'Disable plot title')
group_title.add_argument('--title', type = str,
            help = 'Title of the chromatogram. Example: "Type of sample"')

args = parser.parse_args()

##############################################################################################

### DEFINE FUNCTIONS ###


def get_encoding_type(f):
    """
    get file encoding type
    taken from: https://stackoverflow.com/questions/191359/how-to-convert-a-file-to-utf-8-in-python
    """
    with open(f, 'rb') as e_file:
        rawdata = e_file.read()
    return detect(rawdata)['encoding']


def does_file_exist(f):
    """
    check that the input file exists
    This function is an extra-check but should not be necessary
    as argparse checks whether the input files are readable.
    """
    if not os.path.isfile(f):
        print(f + ' could not be found')
        sys.exit()


def is_csv(f):
    """
    Check that the file is a .csv
    """
    if not f.endswith('.csv'):
        print(f + ' has the wrong file extension.')
        sys.exit()
    else:
        return True


def is_valid_dir(path):
    """
    Check whether the output directory is valid
    Re-adapted from https://stackoverflow.com/questions/38834378/path-to-a-directory-as-argparse-argument
    """
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")
        sys.exit()


def prepare_data_file(f, path):
    """
    Define UTF8-encoded data file
    Return also the filename without extension
    """
    from_codec = get_encoding_type(f)    # get encoding type of the input file
    csv_file = pathlib.Path(f).stem      # get the filename without extension

    if from_codec != 'UTF-8' and from_codec != 'ascii':
        utf8_csv_file = os.path.join(path,csv_file+"_utf8.csv")  # utf-8 file counterpart

        # re-adapted from https://stackoverflow.com/questions/191359/how-to-convert-a-file-to-utf-8-in-python
        # add try: except block for reliability
        try: 
            with open(f, 'r', encoding=from_codec) as original_csv_file, open(utf8_csv_file, 'w', encoding='utf-8') as write_utf8:
                csv_content = original_csv_file.read() # for small files, for big use chunks
                write_utf8.write(csv_content)
                """
                other options:
                - one liner - NOT TESTED YET: python -c "from pathlib import Path; path = Path('yourfile.txt') ; path.write_text(path.read_text(encoding='utf16'), encoding='utf8')"
                - same but in a code snippet - NOT TESTED YET:
                ```
                path = pathlib.Path("yourfile.txt")
                path.write_text(path.read_text(encoding="utf16"), encoding="utf8")
                ```
                """
            return utf8_csv_file, csv_file
        except UnicodeDecodeError:
            print('Decode Error')
            sys.exit()
        except UnicodeEncodeError:
            print('Encode Error')
            sys.exit()

    else:
        return f, csv_file


def array_no_nan(array):
    """
    Get ndarray without missing values
    """
    final_array = np.delete(array, np.where(array == float(-999)))
    return final_array


def frac_ydelta(axes):
    """
    Define y-offset for fraction labels
    """
    y_lim = axes.get_ylim()
    bottom_ylim = y_lim[0]
    y_delta = bottom_ylim + 0.06*y_lim[1]
    return y_delta


def scaler(argmin,argmax,volume_abs_array):
    """
    Define the lower and upper boundaries of the x-axis
    """
    if argmin != None:
        x_min = argmin
    else:
        x_min = volume_abs_array[0]
    if args.xmax != None:
        x_max = argmax
    else:
        x_max = volume_abs_array[-1]
    return argmin, argmax


def frac_bound(peak, frac_ids, frac):
    """
    Define the boundaries of the peak area to color
    """
    if peak != None:
        index_min = frac_ids.index(peak[0])
        index_max = frac_ids.index(peak[1])
        if index_min < index_max:
            vol_min = frac[index_min]
            vol_max = frac[index_max] + frac[1]-frac[0]
            return vol_min, vol_max
        else:
            print("Please enter a valid range of fraction elutions")
            print("You should pass two fractions: The first and the last to color")
            print("########## Exiting now... ")
            sys.exit()
    else:
        return False

#############################################################################################

### CHECK STEPS ###

is_valid_dir(args.log_path)                                        # check whether the specified log-directory is valid
log_file = open(os.path.join(args.log_path,"aktaPlot.out"),'w')    # creates log-file

# Define header info
date = datetime.datetime.now()
python_version = sys.version
header = f'{date}\nThis script was tested in python::3.9\nCurrent python version: {python_version}'

log_file.write(header)    # write out the header on the log-file
print(header)             # print the header on the console

# Check the input files to be processed
for input_file in args.list_csv_files:
    if is_csv(input_file.name):
        does_file_exist(input_file.name)
print('\nAll input files have been checked...')
log_file.write('\n\nAll input files have been checked...')

# print out some information about user's request
dict_run_type = {"sec":"gel filtration", "iex":"ion-exchange chromatography"}
print("\nYou asked for plotting the results of ",dict_run_type[args.run_type], " Akta runs and save the results in ",args.format," format.\n")
log_file.write("\n\nYou asked for plotting the results of "+dict_run_type[args.run_type]+" Akta runs and save the results in "+args.format+" format.\n\n")


### PROCESS DATA AND PLOT ###

print("### -START- ###")
print("\nStart processing input files now...")
log_file.write("\n\n### -START- ###")
log_file.write("\n\nStart processing input files now...")

for input_file in args.list_csv_files:                                                                # loop through each file (one output file for each input file)

    prepare_data_file_out = prepare_data_file(input_file.name, args.log_path)
    csv_data_file = prepare_data_file_out[0]                                                          # ensure to have csv file in UTF-8 format ready to be processed

    csv_data = np.genfromtxt(csv_data_file, skip_header = 3, delimiter = '\t', filling_values = -999) # load data skipping header, converting to float type (default), "" is treated as -999
    
    ## Define data variables ##
    volume_uv_abs = array_no_nan(csv_data[:,0])
    mAU_profile = array_no_nan(csv_data[:,1])
    volume_cond = array_no_nan(csv_data[:,2])
    conductivity = array_no_nan(csv_data[:,3])
    volume_concB = array_no_nan(csv_data[:,4])
    concB = array_no_nan(csv_data[:,5])
    fractions_with_waste = array_no_nan(csv_data[:,10])
    fractions = fractions_with_waste[:len(fractions_with_waste)-2].tolist()                       # do not consider 'waste' fractions, convert ndarray into list
    fractions_ids_with_blank = np.genfromtxt(csv_data_file, dtype = str,
                                            usecols = (11), skip_header = 3, delimiter = '\t')    # get fraction IDs as string type
    fractions_ids = np.delete(fractions_ids_with_blank, 
                            np.where((fractions_ids_with_blank=="") | (fractions_ids_with_blank=='"Waste"'))).tolist() # discard the empty values and waste fractions, convert ndarray into list
    for fraction_id in fractions_ids:
        fraction_id_elements = fraction_id.split('.')
        fractions_ids[fractions_ids.index(fraction_id)] = '.'.join(fraction_id_elements[1:]).strip("\"")

    log_file.write("\n\nData for file "+input_file.name+" defined.")

    ## Plot ##

    fig, ax = plt.subplots()    # single plot (axes object)
    #font = {'size':14}
    #mpl.rc('font', **font)

    # define the elution profile colour
    palette = {'black': '#000000', 'grey': '#808080', 'blue1':'#00004d', 'blue2':'#1f78b4', 'green1':'#4d8964', 'green2':'#56b546ff'}
    if not args.color:
        colour = 'k'
    else:
        colour = palette[args.color]
    
    # define the x-axis boundaries
    scaler_out = scaler(args.xmin,args.xmax,volume_uv_abs)
    x_min = scaler_out[0]
    x_max = scaler_out[1]
    plt.margins(0)     # resert the offset of the axis to 0


    # Show only left and bottom frame
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)

    chromatogram = ax.plot(volume_uv_abs,mAU_profile, color = colour, label = "UV", linewidth = 2)    # plot UV absorbance (elution profile)

    if args.run_type == 'iex':                                                    # if ion-exchange
        concB_plot = ax.plot(volume_concB,concB, label = "ConcB %",               # add plot for concB on same axes as elution profile
                            color = '#00b300', lw = 1.2)

    # plot elution fractions
    if not args.no_elutions:
        frac_y_delta = frac_ydelta(ax)
        for fraction in fractions:
            ax.axvline(x=fraction, ymin=0.065, ymax=0.0, color='#cc0000', linewidth=0.85)
            ax.annotate(fractions_ids[fractions.index(fraction)], xy=(fraction + (fractions[1]-fractions[0]) * 0.55, frac_y_delta),
                         horizontalalignment='center', verticalalignment='bottom', size=2, rotation=90)

    # area of the peak to color
    fill = frac_bound(args.peak,fractions_ids,fractions)
    if fill:
        ax.fill_between(volume_uv_abs, 0, mAU_profile, color = colour, alpha = 0.2,
                        where=((volume_uv_abs>fill[0]) & (volume_uv_abs<fill[1])))

    ax.minorticks_on()                                       # initialize minorticks
    ax.yaxis.set_tick_params(which='minor', bottom=False)    # switch off minor ticks for y-axis
    ax.tick_params(axis='both', which='major', labelsize=12)

    ax.set_xlabel('Elution volume [ml]', fontsize=16)    # x-axis label
    ax.set_ylabel('Abs 280nm [mA.U]', fontsize=16)       # y-axis label
    
    ax.set_xlim(x_min,x_max)                # rescale x-axis
    ax.set_ylim(0,max(mAU_profile)+10)      # hard-coded limit y-axis for absorbance (not conc)

    # plot title
    if not args.no_title:
        if args.title:
            ax.set_title(args.title)
        else:
            ax.set_title("AKTA Chromatogram")

    # plot conductivity
    if args.conductivity:
        ax1 = ax.twinx()                                        # Instantiate a second axes that shares the same x-axis
        conduct = ax1.plot(volume_cond,conductivity, label = "Cond",
                            color = '#4d0000', alpha = 0.75, lw = 0.5)
        max_cond = np.max(conductivity)
        ax1.set_ylim(bottom = 0, top = max_cond + max_cond*0.5) # set lower and upper limit of y-axis for conduct
        ax1.set_ylabel('Conductivity [mS/cm]')
        ax1.yaxis.label.set_color('#4d0000')                    # set color of y-axis label
        
    if not args.no_legend:
        fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax.transAxes,
                    labelspacing = 0.3, fontsize = 12)             # set legend at upper right corner, include both ax and ax1 (if created)

    title_figure = prepare_data_file_out[1] + '_plot.' + args.format    # set figure title, no destination path!
    fig.savefig(os.path.join(args.log_path,title_figure), 
                bbox_inches = 'tight', dpi = int(args.dpi), transparent = True)             # save figure in the log-directory with the desired format and dpi
    plt.clf()                                                           # clear the current figure
    print("Plot for file "+input_file.name + " saved as "+title_figure)
    log_file.write("\nPlot for file "+input_file.name + " saved as "+title_figure)

print("\n### -END- ###")
print(f"\nPlots have been generated successfully. You may find them within {args.log_path} folder.")
log_file.write("\n\n### -END- ###")
log_file.write(f"\n\nPlots have been generated successfully. You may find them within {args.log_path} folder.")
log_file.close()