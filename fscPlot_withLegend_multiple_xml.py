#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a script to plot the FSC curve of cryoEM density maps from RELION post-processing jobs

INPUT: One or multiple postprocess_fsc.xlm files from RELION post-processing jobs. 
       You need to specify the output format and directory where to save your plot as well.
       You can pass a title and disable the gold-standard FSC
OUTPUT: Single plot of FSC curve for all files passed saved in the desired format.

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
    import pandas as pd
    from lxml import etree
except:
    ImportError
    print("Please make sure to have installed matplotlib, pandas and lxml modules.\nYou can do that via conda or pip")
    sys.exit()


parser = argparse.ArgumentParser(description = 'This is a script to plot the FSC curve of cryoEM density maps from RELION post-processing jobs.',
                                 epilog = 'Python 3.9, requirements: matplotlib, pandas.')

parser.add_argument('list_xml_files', metavar='xml', nargs = '+', type = argparse.FileType('r'),
        help = 'list of xml files storing the FSC values for the map(s) you wish to plot')

parser.add_argument('-f', '--format', type = str, choices = ['svg','png','pdf','tiff', 'jpg', 'jpeg'], required = True,
        help = 'specify the output file format')

parser.add_argument('--out_path', type = pathlib.Path, required = True,
        help = 'Directory in which to save the output plot')

parser.add_argument('--title', type = str,
            help = 'Title of the plot. Example: "Type of sample"')

parser.add_argument('--dpi', nargs = '?', default = 300,
        help = 'Set DPI (dots per inch) of the output file for raster images.\
                Default is 300.')    # if no --dpi flag is passed, default value will be used (combination nargs='?' and default)

parser.add_argument('-gs', '--no_goldSt', action = 'store_true',
        help = 'Disable highlighting FSC gold standard 0.143')

parser.add_argument('--legend_elements', nargs = '+',
                    help = 'List of entries for the legend (useful especially when plotting multiple curves')

args = parser.parse_args()

current = os.getcwd()

##############################################################################################

### DEFINE FUNCTIONS ###

def does_file_exist(f):
    """
    check that the input file exists
    This function is an extra-check but should not be necessary
    as argparse checks whether the input files are readable.
    """
    if not os.path.isfile(f):
        print(f + ' could not be found')
        print('Currently you are in' + current)
        sys.exit()


def is_xml(f):
    """
    Check that the file is a .xml
    """
    if not f.endswith('.xml'):
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

#########################################################################################################

is_valid_dir(args.out_path)

# Define header info
date = datetime.datetime.now()
python_version = sys.version
header = f'{date}\nThis script was tested in python::3.9\nCurrent python version: {python_version}'
print(header)


# Check the input files to be processed
for input_file in args.list_xml_files:
    if is_xml(input_file.name):
        does_file_exist(input_file.name)


### Prepare fig and axes objects ###

plt.style.use('seaborn-v0_8-colorblind')
fig, ax = plt.subplots()    # single plot (axes object)

ax.minorticks_on()          # set miinor ticks
ax.set_ylim(0,1.05)

if args.title:
    ax.set_title(args.title)

ax.set_xlabel('Resolution [1/'+r'$\AA$'+']', fontsize = 16)   # x-axis label
ax.set_ylabel('Fourier Shell Correlation', fontsize = 16)     # y-axis label
plt.margins(0)     # resert the offset of the axis to 0

ax.tick_params(axis='both', which='major', labelsize=12)   # set fontsize of major ticks for both axes



# Show only left and bottom frame
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(True)
ax.spines['bottom'].set_visible(True)


# Define the labels for the legend

labels = ['']*len(args.list_xml_files)   # initialize a list of empty labels

if args.legend_elements:
        if len(args.legend_elements) != len(args.list_xml_files):
            print('Please make sure to pass the correct number of labels for the legend')
            sys.exit()
        else:
            for label in args.legend_elements:
                labels[args.legend_elements.index(label)] = label


### Fetch data and plot them ###

for input_file in args.list_xml_files:
    if is_xml(input_file.name):
        does_file_exist(input_file.name)
        data = pd.read_xml(input_file)
        ax.plot(data['x'], data['y'], alpha = 0.75, linewidth = 3,label = labels[args.list_xml_files.index(input_file)])    #, color = '#00004d')

if not args.no_goldSt:
    ax.axhline(y = 0.143, color = 'k', linestyle = '--', lw = 0.5)
    #ax.annotate('FSC 0.143', (0.8,0.25), xycoords='axes fraction')

if args.legend_elements:
    fig.legend(loc = "upper right",frameon = False, bbox_to_anchor=(1,1), fontsize = 'large')

### Save figure ###

title_figure = 'FSC_plot.' + args.format    # set figure title, no destination path!
fig.savefig(os.path.join(args.out_path,title_figure),
                bbox_inches = 'tight', dpi = int(args.dpi), transparent = True)    # save figure in the log-directory with the desired format and dpi
plt.clf()                                                      # clear the current figure
print("FSC plot saved as "+title_figure)