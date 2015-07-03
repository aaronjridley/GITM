#!/usr/bin/env python
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#---------------------------------------------------------------------------
# $Id: load_files.py,v 1.4 2013/12/04 17:02:33 agburr Exp $
#
# load_files
#
# Author: Angeline G Burrell, UMichigan, Sept 2013
#
# Comments: Routines to load data from different instruments
#
# Contains: loadCINDIorbit_ASCII
#           loadMadrigalVTEC_ASCII
#           loadGITMsat_ASCII: AGB added 2013
#           loadMadrigalVTEC_HDF5: AGB added 11/20/13
#---------------------------------------------------------------------------

import string
import read_files as rf
import datetime as dt
import numpy as np

module_name = "load_files"

#---------------------------------------------------------------------------
# loadCINDIorbit_ASCII: A routine to load one or many CINDI orbit-length data
#                       files (the text files that contain obtained from the
#                       website) and load them into a data dictionary.  If
#                       multiple files are loaded, the input is checked to
#                       ensure that multiple observations do not exist for a
#                       single universal time

def loadCINDIorbit_ASCII(input_file, *args, **kwargs):
    '''
    Open one or many CINDI orbit-length ASCII files and load the data into a
    python numpy array after ensuring that multiple observations do not exist
    at a single universal time.

    Input:
    input_file = CINDI data file name or a python list containing multiple
                 filenames

    Output:
    nfiles = number of files loaded successfully
    out    = a dict containing the data in np.arrays, the dict keys are
             specified by the header data line
    '''

    #-----------------------------------------------
    # Initialize the routine variables
    func_name = string.join([module_name, "loadCINDIorbit_ASCII"], " ")

    #-----------------------------------------------
    # Initialize the output
    nfiles = 0
    out = dict()

    #------------------------------------------------
    # Test to see if this is a single file or a list.

    if type(input_file) == list:
        nfiles, out = rf.load_multASCII_data(input_file, "hline", hlines=2)
    elif type(input_file) == str:
        h, out = rf.loadASCII_data_hline(input_file, 2)
        if out:
            nfiles = 1
    else:
        print func_name, "ERROR: unknown input file type", type(input_file)

    # If there is output data, construct a datetime data key and ensure
    # that there is no temporal overlap
    if nfiles > 0 and out:
        # When building the datetime value, test to ensure that the seconds
        # of day do not extend beyond 86400 and that the time is always
        # increasing.
        out['datetime'] = list()
        delindex = list()

        for i,t in enumerate(out['Time']):
            if t < 86400.0:
                this_date = dt.datetime.strptime("{:d} {:d} {:d} {:d} {:d}".format(int(out['Date'][i]/1000), int(out['Date'][i]-int(out['Date'][i]/1000)*1000), int(t / 3600), int((t - int(t / 3600)*3600) / 60), int(t - int((t - int(t / 3600)*3600) / 60) * 60 - int(t / 3600) * 3600)), "%Y %j %H %M %S")

                if i > 0:
                    time_diff = this_date - last_date
                    dsec = time_diff.total_seconds()
                else:
                    dsec = 1

                if dsec > 0:
                    out['datetime'].append(this_date)
                    last_date = this_date

                else:
                    out['datetime'].append(np.nan)
                    delindex.append(i)
            else:
                out['datetime'].append(np.nan)
                delindex.append(i)

        out['datetime'] = np.array(out['datetime'])

        # Remove any data lines where the time did not meet the requirements
        if len(delindex) > 0:
            print func_name, "ADVISEMENT: removing", len(delindex), "duplicate lines from structure"
            for k in out.keys():
                out[k] = np.delete(out[k], delindex)

        # If this is a good timepoint, test for fill values. These cannot be
        # removed with the file read options since a number is used instead of
        # a string
        if len(out['datetime']) > 0:
            largefill_keys = ["Vx", "Vy", "Vz", "Vzonal", "Vpara", "Vmerid",
                              "Ti", "Ni(cm^-3)"]
            smallfill_keys = ["FracO", "FrH", "FrHe", "FrNO"]

            for k in largefill_keys:
                out[k] = [np.nan if v == -99999.0 else v for v in out[k]]
            
            for k in smallfill_keys:
                out[k] = [np.nan if v == 9.99 else v for v in out[k]]

    return nfiles, out

# END loadCINDIorbit_ASCII

#---------------------------------------------------------------------------
# loadMadrigalVTEC_ASCII: A routine to load one or many VTEC data files
#                        (text files that contain obtained from MIT Haystack)
#                        and load them into a data dictionary.

def loadMadrigalVTEC_ASCII(input_file, *args, **kwargs):
    '''
    A routine to load one or many VTEC data files (text files that contain
    obtained from MIT Haystack) and load them into a data dictionary.
   
    Input:
    input_file = VTEC data file name or a python list containing multiple
                 filenames

    Output:
    nfiles = number of files loaded successfully
    out    = a dict containing the data in np.arrays, the dict keys are
             specified by the header data line
    '''

    #-----------------------------------------------
    # Initialize the routine variables
    func_name = string.join([module_name, "loadMadrigalVTEC_ASCII"], " ")
    miss_string = "missing"

    #-----------------------------------------------
    # Initialize the output
    nfiles = 0
    out = dict()

    #------------------------------------------------
    # Test to see if this is a single file or a list.

    if type(input_file) == list:
        nfiles, out = rf.load_multASCII_data(input_file, "hline", hlines=1,
                                             miss=miss_string)
    elif type(input_file) == str:
        h, out = rf.loadASCII_data_hline(input_file, 1, miss=miss_string)
        if out:
            nfiles = 1
    else:
        print func_name, "ERROR: unknown input file type", type(input_file)

    # If there is output data, construct a datetime data key and ensure
    # that there is no temporal overlap
    if nfiles > 0 and out:
        # When building the datetime value, test to ensure that the seconds
        # of day do not extend beyond 86400 and that the time is always
        # increasing.
        out['datetime'] = np.array([dt.datetime.strptime("{:.0f} {:.0f} {:.0f} {:.0f} {:.0f} {:.0f}".format(out['YEAR'][i], out['MONTH'][i], out['DAY'][i], h, out['MIN'][i], out['SEC'][i]), "%Y %m %d %H %M %S") for i,h in enumerate(out['HOUR'])])
        # Change the longitude range to go from 0-360 deg instead of +/- 180 deg
        out['GLON'] = np.array([l if l >= 0. else l+360. for l in out['GLON']])

    return nfiles, out

# END loadMadrigalVTEC_ASCII

#---------------------------------------------------------------------------
# loadGITMsat_ASCII: A routine to load one GITM satellite data input file
#                    into a data dictionary.  If data is available, instead of
#                    just position and time, it will also be included.

def loadGITMsat_ASCII(input_file, *args, **kwargs):
    '''
    Open one GITM satellite data input file into a data dictionary.  If data is
    available, instead of just position and time, it will also be included.

    Input:
    input_file = GITM satellite data file

    Output:
    out = a dict containing the data in np.arrays, the dict keys are specified
          by the GITM standards
    '''

    #-----------------------------------------------
    # Initialize the routine variables
    func_name = string.join([module_name, "loadGITMsat_ASCII"], " ")

    #-----------------------------------------------
    # Initialize the output
    out = dict()

    #-----------------------------------------------
    # Identify the number of header lines
    f = open(input_file, "r")
    line = f.readline()
    nhead = 1

    while line.find("#START") < 0:
        nhead += 1
        line = f.readline()

    #------------------------------------------------
    # Load text file
    temp = np.genfromtxt(input_file, skip_header=nhead)

    #------------------------------------------------
    # Load the data into the output dictionary
    if len(temp) > 0:
        # Assign the date, time, and position information
        out['year'] = temp[:,0]
        out['month'] = temp[:,1]
        out['day'] = temp[:,2]
        out['hour'] = temp[:,3]
        out['min'] = temp[:,4]
        out['sec'] = temp[:,5]
        out['ms'] = temp[:,6]
        out['lon'] = temp[:,7]
        out['lat'] = temp[:,8]
        out['alt'] = temp[:,9]

        # Assign any data
        len_temp = len(temp[0])
        lt = 0
        while len_temp - lt > 10:
            lt += 1
            out_key = "dat{:d}".format(lt)
            out[out_key] = temp[:,9+lt]
        
        # Construct a datetime array
        out['datetime'] = list()

        for t in temp[:]:
            this_date = dt.datetime.strptime("{:d} {:d} {:d} {:d} {:d} {:d} {:d}".format(int(t[0]), int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5]), int(t[6] * 1.0e3)), "%Y %m %d %H %M %S %f")
            out['datetime'].append(this_date)

        out['datetime'] = np.array(out['datetime'])

    return out

# END loadGITMsat_ASCII

#---------------------------------------------------------------------------
# loadMadrigalVTEC_HDF5: A routine to open a HDF5 file and load the data into a
#                        dictonary.

def loadMadrigalVTEC_HDF5(filename, *args, **kwargs):
    '''
    Open an Madrigal TEC HDF5 data files and loads the data into a dictionary
    of python numpy arrays.

    NOTE that HDF files can have datasets hidden under an arbitrary number
    of layers, so it's best to write a reader for each data format.

    Input:
    filename = Madrigal TEC HDF5 filename

    Output:
    out = a dict containing the data in np.arrays, the dict keys are
          specified by the header data line.  A datetime array is also added.
    '''
    import h5py
    func_name = string.join([module_name, "loadMadrigalVTEC_HDF5"], " ")

    # Open the file in read-only mode
    try:
        f = h5py.File(filename, 'r')
    except:
        print func_name, "ERROR: unable to open file [", filename, "]"
        return

    # Load the variables into a standard numpy array and save the description
    # as attributes
    out = dict()

    # Cycle through the file groups and locate the datasets
    try:
        dhandle = f['Data']
    except:
        print func_name, "ERROR: unable to locate Dataset key"
        f.close()
        return

    try:
        dhandle = dhandle['Table Layout']
    except:
        print func_name, "ERROR: unable to locate [Table Layout] key"
        f.close()
        return

    if isinstance(dhandle, h5py.Dataset):
        okeys = dhandle.dtype.fields

        # Load the data
        for ok in okeys.keys():
            out[ok.upper()] = dhandle[ok]

        # Construct the datetime array
        out['datetime'] = np.array([dt.datetime.strptime("{:d} {:d} {:d} {:d} {:d} {:d}".format(out['YEAR'][i], out['MONTH'][i], out['DAY'][i], h, out['MIN'][i], out['SEC'][i]), "%Y %m %d %H %M %S") for i,h in enumerate(out['HOUR'])])

        # Change the longitude range to go from 0-360 deg instead of +/- 180 deg
        out['GLON'] = np.array([l if l >= 0. else l+360. for l in out['GLON']])


    else:
        print func_name, "ERROR: dataset not where expected"

    # Close HDF5 file handle and return data
    f.close()
    return(out)

# END loadMadrigalVTEC_HDF5
