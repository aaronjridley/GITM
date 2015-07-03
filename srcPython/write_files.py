#!/usr/bin/env python
#---------------------------------------------------------------------------
# write_files
#
# Author: Angeline G Burrell, UMichigan, Feb 2014
#
# Comments: Routines to support file processing in python
#
# Contains: writeASCII_file
#---------------------------------------------------------------------------

import string
import numpy as np
from os import path

module_name = "write_files"

#---------------------------------------------------------------------------
# writeASCII_file: A routine to create an ascii file from a string or list
#                  of strings.  Will overwrite a file with the same name if
#                  such a file exists.

def writeASCII_file(filename, datalines, *args, **kwargs):
    '''
    A routine to create an ascii file from a string or list of strings.  Will
    overwrite a file with the same name if such a file exists

    Input:
    filename  = output file name
    datalines = string or list of strings to be written
    '''

    func_name = string.join([module_name, "writeASCII_file"], " ")

    #-----------------------------------------------------------------------
    # Open and test the file to ensure it can be written
    try:
        f = open(filename, 'w')

        # Print data after determinging the appropriate data type
        try:
            dlen = len(datalines)

            if type(datalines) is str:
                f.write(datalines)
            else:
                for line in datalines:
                    line = line+"\n"
                    f.write(line)
        except:
            f.write(datalines)

        f.close()
    except:
        print func_name, "ERROR: unable to open [", filename, "]"

# End writeASCII_file


#---------------------------------------------------------------------------
# writeASCII_data_w_sorttext: A routine to create an ascii file from a dict
#                             of dicts.  Will overwrite a file with the same
#                             name if such a file exists.

def writeASCII_data_w_sorttext(filename, data, dt_key=None, sort_key=None,
                               *args, **kwargs):
    '''
    A routine to create an ascii file from a dict of dicts (which contain 
    numpy arrays).  Will overwrite a file with the same name if such a file
    exists

    Input:
    filename = output file name
    data     = dict of dicts of numpy arrays
    dt_key   = Name of key containing datetime data (default=None)
    sort_key = Name of key to sort output lines by (default=None)
    '''

    func_name = string.join([module_name, "writeASCII_data_w_sorttext"], " ")

    #---------------------------------------------
    # Construct formatted strings for each line
    outlines = list()
    hkeys = data.keys()
    skeys = data[hkeys[0]].keys()

    if sort_key is not None:
        if data.has_key(sort_key):
            sort_dat = list()
        else:
            sort_key = None
            print func_name, "WARNING: sort key does not exist [", sort_key, "]"

    if dt_key is not None:
        if data.has_key(dt_key):
            hkeys.pop(hkeys.index(dt_key))
        else:
            dt_key = None

    # Format the header line
    hline = string.join([hkey for hkey in hkeys], " ")

    if dt_key is not None:
        hline = "#Date Time Key {:s}".format(hline)
    else:
        hline = "#Key {:s}".format(hline)

    # Format the data lines
    for skey in skeys:
        slen = len(data[hkeys[0]][skey])
        for i in range(slen):
            line = string.join(["{:}".format(data[hkey][skey][i])
                                for hkey in hkeys], " ")
            line = "{:s} {:s}".format(skey, line)

            if dt_key is not None:
                line = "{:} {:s}".format(data[dt_key][skey][i], line)

            outlines.append(line)

            if sort_key is not None:
                sort_dat.append(data[sort_key][skey][i])

    # Sort the output by time
    if sort_key is not None:
        from operator import itemgetter
        (sort_dat, outlines) = zip(*sorted(zip(sort_dat, outlines),
                                         key=itemgetter(0)))
        outlines = list(outlines)
        del sort_dat

    # Add the header line
    outlines.insert(0, hline)

    # Print the output file
    writeASCII_file(filename, outlines)
    return outlines

# End writeASCII_data_w_sorttext
