'''
Command line tool for setting up different versions of the
Hosono code. 
'''

import argparse
import os

# write help message description
desc = 'Command line tool for setting up different\n'+\
    'versions of Hosono code. This should be used \n'+\
    'before compilation in order to specify, for \n'+\
    'example, if the code should be compiled with \n'+\
    'MHD on or off.'

# create argument parser object
parser = argparse.ArgumentParser(prog='setup',
                                 description=desc,
                                 epilog='Written by Sean Riedel')

parser.add_argument('--mhd', action='store_true',
                    help='Setup with mhd enabled')
parser.add_argument('-b', '--build_dir', default='_build',
                    help='Directory where symlinks'
                    'will be generated and code will be compiled')

# parse the command line arguments to the script
args = parser.parse_args()

# set the build directory
if args.build_dir=='_build':
    if args.mhd:
        build_dir = '_build_mhd'
    else:
        build_dir = '_build_hydro'
else:
    build_dir = args.build_dir


# List of common source files (including Makefile)
commonF90list = 'main BCs constants eos getFluxes_mod grid '+\
    'ICs output readParamFile reconstruct simulation '+\
    'slopeLimiter update'
commonF90list = [f+'.F90' for f in commonF90list.split()]

# list of source files that are different depending
# on whether using hydro or MHD
diffntF90list = 'cfl riemannsolver convert eigen'
diffntF90list = [f+'.F90' for f in diffntF90list.split()]
diffntF90list.append('definition.h')
diffntF90list.append('Makefile') # include Makefile

# list of source files that are exclusive to MHD
MHDonlyF90list = 'divB'
MHDonlyF90list = [f+'.F90' for f in MHDonlyF90list.split()]

# get current working directory
cwd = os.getcwd()

# make build directory
os.mkdir(build_dir)

# change from hydro to MHD if requested
hydro_or_mhd = 'hydro'
if args.mhd:
    hydro_or_mhd = 'MHD'

# symlink the rest of the necessary files into the build directory
for f in commonF90list:
    os.symlink(cwd+'/source/'+f, cwd+"/"+build_dir+"/"+f)
for f in diffntF90list:
    os.symlink(cwd+'/source/'+hydro_or_mhd+"/"+f, cwd+"/"+build_dir+"/"+f)
for f in MHDonlyF90list:
    os.symlink(cwd+'/source/MHD/'+f, cwd+"/"+build_dir+"/"+f)