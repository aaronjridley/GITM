#!/usr/bin/env python

import sys
import datetime as dt
import numpy as np
from omniweb import *
from swmf_imf import *
import argparse
import matplotlib.pyplot as plt
import matplotlib.dates as dates

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(description = 'Plot AMIE files')
    parser.add_argument('start', metavar = 'start', nargs = 1, \
                        help = 'start date as YYYYMMDD')
    parser.add_argument('end', metavar = 'end', nargs = 1, \
                        help = 'end date as YYYYMMDD')
    parser.add_argument('-swmf', \
                        help='output swmf style file (imfYYMMDD.dat)', \
                        action="store_true")

    args = parser.parse_args()

    return args

def write_omni_file(lines, fileout):

    with open(fileout, "w") as file:
        file.writelines("%s" % l for l in lines)

    
#------------------------------------------------------------------------------
#SCRIPT USE
#example command line input: python omniweb_read.py 20110620 20110623 -all

args = parse_args()

#assuming first two args are the start/end dates, then info desired

start = args.start
if (not np.isscalar(start)):
    start = start[0]
end = args.end
if (not np.isscalar(end)):
    end = end[0]

results = download_omni_data(start, end, "-all")
omniDirty = parse_omni_data(results)
data = clean_omni(omniDirty)

if (args.swmf):
    fileout = data["times"][0].strftime('imf%Y%m%d.dat')
    message = "Data downloaded from OMNIWeb and processed by omniweb_read.py\n"
    write_swmf_imf_file(data, fileout, message)
else:
    fileout = data["times"][0].strftime('omni_%Y%m%d.txt')
    write_omni_file(results, fileout)
    
plot_imf(data)
