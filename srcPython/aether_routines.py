#!/usr/bin/env python

from datetime import datetime
from datetime import timedelta
import numpy as np
from netCDF4 import Dataset

#-----------------------------------------------------------------------------
# Get some ancillary information from the netcdf file.  This simply
# mirrors the gitm header info.
# Inputs:
#  - filelist: a list of netcdf files. Really, only the last file is used.
#              also, the number of files is recorded.
# Output:
#  - header: A dictionary containing information about the netcdf file, such
#            as nLons, nLons, nAlts, nVars, variable names, time(s)
#-----------------------------------------------------------------------------

def read_aether_header(filelist):

    header = {}
    header["nFiles"] = len(filelist)
    header["version"] = 0.1
    header["nLons"] = 0
    header["nLats"] = 0
    header["nAlts"] = 0
    header["nVars"] = 0
    header["vars"] = []
    header["time"] = []
    header["filename"] = []

    header["filename"].append(filelist[-1])
    
    ncfile = Dataset(filelist[-1],'r')
    header["nVars"] = 0
    for var in ncfile.variables.values():
        s = var.shape
        if (len(s) == 3):
            header["nLons"] = s[0]
            header["nLats"] = s[1]
            header["nAlts"] = s[2]
            header["vars"].append(var.name)
            header["nVars"] = header["nVars"] + 1

    base = datetime(1965,1,1,0,0,0)
    time = np.array(ncfile.variables['Time'])
    current = base + timedelta(seconds = time[0])
    header["time"].append(current)
            
    ncfile.close()

    return header
    
#-----------------------------------------------------------------------------
# Read in some variables from one file. Inputs:
#   - file: netcdf file to read
#   - vars: list of variable NAMES to read
# Output:
#   - data["time"]: datetime of the file
#   - data[NUMBER]: data that is read in.
#                   NUMBER does from 0 - number of vars read in (0-3 typical)
#-----------------------------------------------------------------------------

def read_aether_one_file(file, vars):

    data = {}
    
    ncfile = Dataset(file,'r')

    base = datetime(1965,1,1,0,0,0)
    time = np.array(ncfile.variables['Time'])
    current = base + timedelta(seconds = time[0])
    data["time"] = current

    iVar = 0
    for var in vars:
        print(iVar,var)
        data[iVar] = np.array(ncfile.variables[var])
        iVar = iVar + 1
    ncfile.close()

    return data
    


