#!/usr/bin/env python
#---------------------------------------------------------------------------
# read_files
#
# Author: Angeline G Burrell, UMichigan, Feb 2013
#
# Comments: Routines to support file processing in python
#
# Contains: loadASCII_data_header
#           loadASCII_data_hline
#           loadASCII_index_profile
#           loadnetCDF_data
#---------------------------------------------------------------------------

import string
import numpy as np
from os import path

module_name = "read_files"

#---------------------------------------------------------------------------
# loadASCII_data_header: A routine to open an ascii file with a single header
#                        line marked by a '#' containing the data column names.
#                        Additional header lines preceeding the data name line
#                        are allowed and imported to a list.
#
# Corrections: AGB 3/14/13 - added check for filesize, changed naming
#                            convention to PEP8 standards

def loadASCII_data_header(filename, miss=None, fill=np.nan, *args, **kwargs):
    '''
    Open an ascii data file and load it into a python numpy array.

    Input:
    filename = CINDI data file name
    miss     = string or list denoting missing value options (default=None)
    fill     = fill value (default = NaN)

    Output:
    header = a list containing the header strings without the '#'
    out    = a dict containing the data in np.arrays, the dict keys are
             specified by the header data line
    '''

    func_name = string.join([module_name, "loadASCII_data_header"], " ")

    #-----------------------------------------------------------------------
    # Test to ensure the file is small enough to read in.  Python can only
    # allocate 2GB of data.  If you load something larger, python will crash

    fsize  = path.getsize(filename)
    header = list()
    out    = dict()

    if(fsize > 2.0e9):
        print func_name, "WARNING: File size [", (fsize * 1e-9), "GB > 2 GB]"
        return(header, out)
    elif(fsize == 0):
        print func_name, "WARNING: empty file [", filename, "]"
        return(header, out)

    #----------------------------------------------
    # Open the datafile and read the header rows

    f = open(filename, "r")

    if not f:
        print func_name, "ERROR: unable to open input file [", filename, "]"
        return header, out

    line   = f.readline()
    check  = 0

    while line.find("#") >= 0:
        hline  = string.strip(line.replace("#", ""))
        line   = f.readline()
        check += 1

        if(len(hline) > 0):
            header.append(hline)

    if(check > 0):
        hline = hline.split()
    else:
        print func_name, "ERROR: no header in this file [", filename, "]"
        return(header, out)

    #-------------------------------------------
    # Open the datafile and read the data rows

    temp = np.genfromtxt(filename, comments="#", missing_values=miss,
                         filling_values=fill)

    if len(temp) > 0:
        #------------------------------------------
        # Create the output dictionary

        for num,name in enumerate(hline):
            out[name] = temp[:,num]

    del temp
    return(header, out)

# End loadASCII_data_header

#---------------------------------------------------------------------------
# loadASCII_data_hline: A routine to open an ascii file with a header
#                       of a specified length.  The last header line is assumed
#                       to contain the data column names.

def loadASCII_data_hline(filename, hlines, miss=None, fill=np.nan, *args,
                         **kwargs):
    '''
    Open an ascii data file and load it into a python numpy array.  File
    header may be any number of lines and does not have to be preceeded by a
    particular character.

    Input:
    filename = CINDI data file name
    hlines   = number of lines in header
    miss     = string or list denoting missing value options (default=None)
    fill     = fill value (default = NaN)

    Output:
    header = a list containing all specified header lines
    out = a dict containing the data in np.arrays, the dict keys are
          specified by the header data line
    '''

    func_name = string.join([module_name, "loadASCII_data_hline"])

    #-----------------------------------------------------------------------
    # Test to ensure the file is small enough to read in.  Python can only
    # allocate 2GB of data.  If you load something larger, python will crash

    fsize = path.getsize(filename)
    header = list()
    out = dict()

    if(fsize > 2.0e9):
        print func_name, "WARNING: File size [", (fsize * 1e-9), "GB > 2 GB]"
        return hearder, out
    elif(fsize == 0):
        print func_name, "WARNING: empty file [", filename, "]"
        return(header, out)

    #----------------------------------------------
    # Open the datafile and read the header rows

    f = open(filename, "r")

    if not f:
        print func_name, "ERROR: unable to open input file [", filename, "]"
        return out

    for h in range(hlines):
        header.append(f.readline())

    #-------------------------------------------
    # Open the datafile and read the data rows

    temp = np.genfromtxt(filename, skip_header=hlines, missing_values=miss,
                         filling_values=fill)

    if len(temp) > 0:
        #---------------------------------------------------------------------
        # Create the output dictionary, removing the point sign from any keys

        for num,name in enumerate(header[-1].split()):
            name = name.replace("#", "")

            if len(name) > 0:
                out[name] = temp[:,num]

    del temp
    return header, out

# End loadASCII_data_hline

def loadASCII_index_profile(filename, miss=None, fill=np.nan, *args, **kwargs):
    '''
    Open an ascii data file and load it into a python numpy array.  Assumes
    this file is seperated into index blocks, which should be maintained.

    Input:
    filename = CINDI data file name
    miss     = string or list denoting missing value options (default=None)
    fill     = fill value (default = NaN)

    Output:
    header  = a list containing the header strings without the '#'
    out     = a dict containing the data in np.arrays, the dict keys are
              specified by the header data line
    nblocks = number of indexed data blocks
    '''

    func_name = string.join([module_name, "loadASCII_data_header"])

    #-----------------------------------------------------------------------
    # Test to ensure the file is small enough to read in.  Python can only
    # allocate 2GB of data.  If you load something larger, python will crash

    fsize = path.getsize(filename)
    header = list()
    out = dict()
    nblocks = 0

    if(fsize > 2.0e9):
        print func_name, "WARNING: File size [", (fsize * 1e-9), "GB > 2 GB]"
        return(header, out, nblocks)
    elif(fsize == 0):
        print func_name, "WARNING: empty file [", filename, "]"
        return(header, out, nblocks)

    #----------------------------------------------
    # Open the datafile and read the header rows

    f = open(filename, "r")

    if not f:
        print func_name, "ERROR: unable to open input file [", filename, "]"
        return header, out, nblocks

    line = f.readline()
    check = 0

    while line.find("#") >= 0:
        hline = string.strip(line.replace("#", ""))
        line = f.readline()
        check += 1

        if(len(hline) > 0):
            header.append(hline)

    if(check > 0):
        hline = hline.split()
    else:
        print func_name, "ERROR: no header in this file [", filename, "]"
        return(header, out, nblocks)

    #-------------------------------------------------
    # Cycle through the data rows, identifying blocks

    while len(line) > 0:
        if (len(line) == 1 and line.find("\n") == 0):
            # A new block has been found.  Only incriment if data has been read
            if len(out) > 0:
                nblocks += 1
            print "TEST: at block", nblocks, len(line)

            # Cycle to new dataline
            while len(line) == 1 and line.find("\n") == 0:
                line = f.readline()

        # Load the dataline into the output structure
        dline = line.split()

        for num,name in enumerate(hline):
            if out.has_key(name):
                if len(out[name]) < nblocks:
                    out[name][nblocks].append(dline[num])
                else:
                    out[name].append([dline[num]])
            else:
                out[name] = [[dline[num]]]

    return(header, out, nblocks)

# End loadASCII_index_profile

#---------------------------------------------------------------------------
# load_multASCII_data: A routine to open multiple ascii files of the same
#                      format and load them into the same data array.

def load_multASCII_data(filelist, file_type, miss=None, fill=np.nan,
                        hlines=-1, *args, **kwargs):
    '''
    Load a list of ASCII files into a single data dictionary

    Input:
    filelist  = List of ASCII files to read
    file_type = types: index/header/hline
                index  -> files are seperated into index blocks
                header -> files have header lines signaled by '#'
                hline  -> files have same-length headers, specified in input
    miss      = string or list denoting missing value options (default=None)
    fill      = fill value (default = NaN)
    hlines    = number of header lines (only used if file_type is hline)
    '''

    import copy

    func_name = string.join([module_name, "load_multASCII_data"])

    out = dict()
    nfiles = 0

    #-----------------------------------
    # Cycle through the list of files

    for filename in filelist:
        if file_type.find("index") >= 0:
            (h, o, n) = loadASCII_index_profile(filename, miss, fill)
        elif file_type.find("header") >= 0:
            (h, o) = loadASCII_data_header(filename, miss, fill)
        elif file_type.find("hline") >= 0 and hlines >= 0:
            (h, o) = loadASCII_data_hline(filename, hlines, miss, fill)
        else:
            print func_name, "ERROR: uknown file type [", file_type, "]"
            return nfiles, out

        if o:
            nfiles += 1

            # Add data column for the filename
            fsplit = string.split(filename, "/")
            k = o.keys()
            o['filename'] = [fsplit[-1]] * len(o[k[0]])

            if out:
                out = combine_data_dictionaries([out, o])
            else:
                out = copy.deepcopy(o)

    return nfiles, out

# End load_multASCII_data

#---------------------------------------------------------------------------
# loadnetCDF_data: A routine to open a netCDF file and load the data into a
#                  dictonary.

def loadnetCDF_data(filename, mult=False, *args, **kwargs):
    '''
    Open one or many netCDF data files and load them into a python numpy array.

    Input:
    filename = netCDF filename or a string with wildcards that would include
               all desired netCDF files if used for an 'ls' command.  NetCDF4
               files are not supported. Ex: see__L3_*.ncdf
    mult = True for multiple files, False for single (default is False)

    Output:
    out   = a dict containing the data in np.arrays, the dict keys are
            specified by the header data line
    units = a dict containing the units for each data key
    desc  = a dict containing the descriptions for each data key
    '''

    import netCDF4 as cdf
    func_name = string.join([module_name, "loadnetCDF_data"], " ")

    # Open the file
    f = list()

    if mult:
        try:
            temp = cdf.MFDataset(filename)
            f.append(temp)
        except:
            # Probably has the wrong file format.  You'll need to cycle through
            # these individually now.
            import subprocess

            command = "ls %s" % (filename)
            pipe = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
            names = pipe.stdout.readlines()

            for name in names:
                print "Opening:", name[0:-1]
                f.append(cdf.Dataset(name[0:-1]))
    else:
        f.append(cdf.Dataset(filename))

    # Load the variables into a standard numpy array and save the description
    # as attributes
    out = dict()
    units = dict()
    desc = dict()
    unit_warning = False
    desc_warning = False

    for fdat in f:
        slist = []
        for fkeys in fdat.variables.keys():
            s = fdat.variables[fkeys].shape

            if(len(s) > 0):
                slist.append(s[0])
        dim1 = np.max(slist)

        for fkeys in fdat.variables.keys():
            # Ensure the data will be stored as an array
            if len(fdat.variables[fkeys].shape) < 1:
                farray = [fdat.variables[fkeys][0] for i in range(dim1)]
            else:
                farray = fdat.variables[fkeys][:]

            if out.has_key(fkeys):
                out[fkeys] = np.append(out[fkeys], farray)
            else:
                out[fkeys] = np.array(farray)

            try:
                if units.has_key(fkeys):
                    units[fkeys] = np.append(units[fkeys],
                                             fdat.variables[fkeys].Units)
                else:
                    units[fkeys] = fdat.variables[fkeys].Units
            except AttributeError:
                if not unit_warning:
                    print module_name, "ADVISEMENT: no unit attribute"
                    unit_warning = True
                
            try:
                if desc.has_key(fkeys):
                    desc[fkeys] = np.append(desc[fkeys],
                                            fdat.variables[fkeys].Description)
                else:
                    desc[fkeys] = fdat.variables[fkeys].Description
            except AttributeError:
                if not desc_warning:
                    print module_name, "ADVISEMENT: no description attribute"
                    desc_warning = True

        fdat.close()
    return(out, units, desc)

# END loadnetCDF_data

#---------------------------------------------------------------------------
# combine_data_dictionaries: A routine to combine a list of data dictionaries,
#                            extending each of the numpy arrays for overlapping
#                            keys.

def combine_data_dictionaries(data_dicts, *args, **kwargs):
    '''
    Ended arrays from matching keys in the list of data dictionaries

    Input:
    data_dicts = list of data dictionaries
    
    Output:
    combinded_dict = single dictionary containing np arrays for each key
                     found in every listed data dictionary
    '''

    rname = string.join(module_name, "combine_data_dictionaries")
    cdict = dict()
    start = True

    # Identify the common data keys
    for ddict in data_dicts:
        if start:
            ckeys = ddict.keys()
            start = False
        elif len(ckeys) == 0:
            print rname, "ERROR: no overlapping keys found"
            return(cdict)
        else:
            nkeys = list()

            for key in ckeys:
                if ddict.has_key(key):
                    nkeys.append(key)

            ckeys = nkeys
       
    if len(ckeys) == 0:
        print rname, "ERROR: no overlapping keys found"
        return(cdict)

    # Combine the different data arrays

    for key in ckeys:
        for ddict in data_dicts:
            if cdict.has_key(key):
                cdict[key] = np.append(cdict[key], ddict[key])
            else:
                cdict[key] = np.array(ddict[key])
                
    return(cdict)
