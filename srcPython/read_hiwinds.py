#!/usr/bin/env python

import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import argparse
import re

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args_hiwinds():

    parser = argparse.ArgumentParser(description =
                                     'Read and plot HI-WIND data file(s)')
    parser.add_argument('files',
                        nargs = '+',
                        help = 'List of files to read and plot')

    args = parser.parse_args()

    return args

def convert_time(yearday):
    year = int(yearday/1000.0)
    base = dt.datetime(year,1,1)
    doy = yearday - year * 1000.0
    time = base + dt.timedelta(days = doy-1)
    return time

def get_time_from_filename(filename):
    print(filename)
    m = re.match(r'.*_(\d*).asc', filename)
    if m:
        yearday = float(m.group(1))
        time = convert_time(yearday)
    else:
        time = dt.datetime(1970,1,1)
        print("I don't understand the filename, and can't get the date!!!")
    print("Setting start time to : ", time)
    return time

def read_hiwind_file(file):

    basetime = get_time_from_filename(file)
    
    fpin = open(file, 'r')

    header = fpin.readline()

    lats = []
    lons = []
    merid_time = []
    merid_wind = []
    merid_unc = []
    zonal_time = []
    zonal_wind = []
    zonal_unc = []

    for line in fpin:
        aline = line.split()
        
        if (len(aline) > 7):
        
            # check wind speeds:
            merid = float(aline[3])
            merid_err = float(aline[4])
            zonal = float(aline[6])
            zonal_err = float(aline[7])

            if ((np.abs(merid) < 500) and
                (np.abs(merid_err) < 500) and
                (np.abs(zonal) < 500) and
                (np.abs(zonal_err) < 500)):
        
                lats.append(float(aline[0]))
                lons.append(float(aline[1]))
                utm = basetime + dt.timedelta(seconds = float(aline[2])*3600.0)
        
                merid_time.append(utm)
                merid_wind.append(merid)
                merid_unc.append(merid_err)
                utz = basetime + dt.timedelta(seconds = float(aline[5])*3600.0)
                zonal_time.append(utz)
                zonal_wind.append(zonal)
                zonal_unc.append(zonal_err)

    data = {'lats' : np.array(lats),
            'lons' : np.array(lons),
            'merid_time' : np.array(merid_time),
            'merid_wind' : np.array(merid_wind),
            'merid_unc' : np.array(merid_unc),
            'zonal_time' : np.array(zonal_time),
            'zonal_wind' : np.array(zonal_wind),
            'zonal_unc' : np.array(zonal_unc)}

    return data

#------------------------------------------------------------------------------
# main code
#------------------------------------------------------------------------------

if __name__ == '__main__':  # main code block

    args = parse_args_hiwinds()

    filelist = args.files

    for file in filelist:
        print('Reading file: ',file)
        data = read_hiwind_file(file)

        mint = np.min([np.min(data['merid_time']), np.min(data['zonal_time'])])
        maxt = np.max([np.max(data['merid_time']), np.max(data['zonal_time'])])

        sTime = mint.strftime('%b %d, %Y %H:%M UT') + ' - ' + \
            maxt.strftime('%b %d, %Y %H:%M UT (Hours)')
        
        fig = plt.figure(figsize = (10,10))
        axm = fig.add_subplot(211)
        axm.plot(data['merid_time'], data['merid_wind'])
        axm.set_xlim(mint, maxt)
        axm.set_ylabel('Meridional Wind (m/s)')
        axm.axhline(y=0.0, color='r', linestyle=':')
        
        axz = fig.add_subplot(212)
        axz.plot(data['zonal_time'], data['zonal_wind'])
        axz.set_xlim(mint, maxt)
        axz.set_ylabel('Zonal Wind (m/s)')
        axz.set_xlabel(sTime)
        axz.axhline(y=0.0, color='r', linestyle=':')

        plotfile = data['merid_time'][0].strftime('hiwind%Y%m%d.png')
        print('writing : ',plotfile)    
        fig.savefig(plotfile)
        plt.close()

        
        
