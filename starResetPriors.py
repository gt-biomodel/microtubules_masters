#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a script to re-set psi, tilt or both angles to their priors.
Optionally, it can also print the max number of filaments per micrograph?. 
Mind that this is done by getting the highest filament ID number, so if the 
filaments underwent different processing steps and might have been fileterd, 
need still to check whether the IDs are redefined accordingly or not!!!
Tested for RELION 3.1
Requirements: numpy

Author: Giulia Tonon
ZMBH Heidelberg, AG Pfeffer
"""

import os
import sys
import argparse
import pathlib
try:
    import numpy as np
except:
    ImportError
    print("Please make sure to have installed numpy")
    sys.exit()


parser = argparse.ArgumentParser(description = 'This is a script to re-set psi, tilt or both angles of helical particles to their priors.\
                                                Optionally, you can also ask for printing the number of filaments in the console.',
                                epilog = 'Tested for RELION 3.1. Requirements: numpy')

parser.add_argument('input_starFile', metavar='PATHor', type = argparse.FileType('r'),
        help = 'Path to STAR-file of your helical particles.')

parser.add_argument('new_starFile', metavar = 'PATHnew', type = argparse.FileType('w'),
        help = 'Path to the new STAR-file.')

parser.add_argument('angle', nargs = '+', type = str,
        help = 'Angles that you wish to reset. You can pass psi and tilt')

parser.add_argument('-nf', '--nr_filaments', action = 'store_true',
        help = 'Ask for printing the max number of filaments per micrograph?')

args = parser.parse_args()


current = os.getcwd()        # current directory
index_angles_to_change = []  # list of indices columns of angles to reset
index_angles_for_change = [] # list of indices columns of corresponding priors

##############################################################################################

### DEFINE FUNCTIONS ###

def does_file_exist(f):
    """
    check that the input file exists
    This function is an extra-check but should not be necessary
    as argparse checks whether the input files are readable.
    """
    if not os.path.isfile(f):
        print('')
        print('-> I did not find ' + f + '!')
        print('-> Please make sure the file exists and has the right name.')
        print('-> It could also be that you run this script from the wrong directory.')
        print('-> Currently you are in' + current)
        sys.exit()


def is_star(f):
    """
    Check that the file is a .star
    """
    if not f.endswith('star'):
        print('--> Input is not a .star')
        print('     type help')
        print('')
        sys.exit()
    else:
        return True


def header_parser(input_file, output_file):
    """
    Reading input file
    Copy the header in the new STAR-file
    Count how many lines compose the header
    Get the list of data_particle IDs
    """
    open_input_file = open(input_file, 'r')
    open_output_file = open(output_file, 'w')
    counter_lines_header = 0
    linesstarted = False
    data_particles = False
    list_data_particles = []     # list of data_particles IDs

    for star_line in open_input_file.readlines():
    
        if data_particles == True and star_line[0] != "_" and star_line.rstrip().endswith("_") == False and star_line.rstrip() != "" and star_line != "\n":
            linesstarted = True
            break

        # fill list of data_particles IDs
        elif data_particles == True:
            if star_line[0] == "_":
                star_line_elements = star_line.split(' ')
                list_data_particles.append(star_line_elements[0].lstrip('_rln'))

        else:
            if star_line.rstrip() == 'data_particles':
                data_particles = True

        counter_lines_header = counter_lines_header+1
        open_output_file.write(star_line)

    open_output_file.close()
    return counter_lines_header, list_data_particles


##############################################################################################

# Check for valid input and write out the header
if is_star(args.input_starFile.name):
    does_file_exist(args.input_starFile.name)
    output_header_parser = header_parser(args.input_starFile.name, args.new_starFile.name)
    counter_lines_header = output_header_parser[0]
    list_data_particles = output_header_parser[1]

# Fill the lists of column indices to swap
for angle in args.angle:
    if angle == 'psi':
        index_angles_to_change.append(list_data_particles.index('AnglePsi'))
        index_angles_for_change.append(list_data_particles.index('AnglePsiPrior'))
    elif angle == 'tilt':
        index_angles_to_change.append(list_data_particles.index('AngleTilt'))
        index_angles_for_change.append(list_data_particles.index('AngleTiltPrior'))
    else:
        print("You passed an incorrect angle name.\nMind that you can only pass 'psi' or 'tilt' or both of them")
        sys.exit()

# Read block of particle data and reset angles
particle_data = np.loadtxt(args.input_starFile.name, skiprows=counter_lines_header, dtype = str)
if particle_data.shape[1] == len(list_data_particles):
    for col in index_angles_to_change:
        particle_data[:,col] = particle_data[:,index_angles_for_change[index_angles_to_change.index(col)]]

    # write out the new particle data
    new_starFile_open = open(args.new_starFile.name, 'a')
    for line in range(0,particle_data.shape[0]):
        new_line = '\t'.join(particle_data[line, :])
        new_starFile_open.write(new_line+'\n')
    new_starFile_open.close()

    # print the max number of filaments per micrograph?
    if args.nr_filaments:
        filament_IDs = particle_data[:, list_data_particles.index('HelicalTubeID')]
        nr_filaments = max([eval(i) for i in filament_IDs])
        print('The max number of filaments per micrograph in your dataset is '+ str(nr_filaments))
        print(filament_IDs)
else:
    raise Exception('The data of the particles do not match their identifier')
    sys.exit()

print('Script finished successfully.')