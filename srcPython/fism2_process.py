#!/usr/bin/env python

import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import argparse
import os
import re

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args_fism():

    parser = argparse.ArgumentParser(description = 'Download and process FISM data')
    parser.add_argument('-start', \
                        help = 'start date as YYYYMMDD', default = '0')
    parser.add_argument('-end', \
                        help = 'end date as YYYYMMDD', default = '0')
    parser.add_argument('-euvfile', \
                        help='EUV file that provides wavelength bins',
                        default = 'euv_59.csv')
    parser.add_argument('-fismfile', \
                        help='FISM file to convert',
                        default = 'none')
    parser.add_argument('-gitm', \
                        help='Write out GITM-style file',
                        action='store_true')
    parser.add_argument('-flare', \
                        help='Download and process FISM flare data',
                        action='store_true')

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# Take string time YYYYMMDD.HHMM and convert to datetime
# ----------------------------------------------------------------------

def convert_ymdhm_to_dt(ymdhm):

    year = ymdhm[0:4]
    mo = ymdhm[4:6]
    da = ymdhm[6:8]

    if (len(ymdhm) >= 11):
        hr = ymdhm[9:11]
        if (len(ymdhm) >= 13):
            mi = ymdhm[11:13]
        else:
            mi = '00'
    else:
        hr = '00'
        mi = '00'
        
    time = dt.datetime(int(year), int(mo), int(da), int(hr), int(mi), 0)
    return time

# ----------------------------------------------------------------------
# Function to download FISM2 data.
# Need call to be in form:
# https://lasp.colorado.edu/lisird/latis/dap/fism_daily_hr.csv?&time>=2020-01-01T00:00:00.000Z&time<=2020-01-31T00:00:00.000Z
# ----------------------------------------------------------------------

def download_fism2(start, end, isFlare):

    sStart = start.isoformat()+'.000Z'
    sEnd = end.isoformat()+'.000Z'

    if (isFlare):
        site = 'https://lasp.colorado.edu/lisird/latis/dap/fism_flare_hr.csv'
    else:
        site = 'https://lasp.colorado.edu/lisird/latis/dap/fism_daily_hr.csv'

    url = site + '?&time>=' + sStart + '?&time<=' + sEnd
    
    ymdS = start.strftime('%Y%m%d')
    ymdE = end.strftime('%Y%m%d')
    filename = '.fism_raw_' + ymdS + '_to_' + ymdE + '.txt'

    command = 'curl "' + url + '" > ' + filename
    print("Running Command : ", command)
    os.system(command)
    
    return filename
    

#------------------------------------------------------------------------------
# Convert yearday (YYYYDDD) to date time
#   - center the day at 12 UT.
#------------------------------------------------------------------------------

def convert_time(yearday):
    year = int(yearday/1000.0)
    base = dt.datetime(year,1,1,12)
    doy = yearday - year * 1000.0
    time = base + dt.timedelta(days = doy-1)
    return time

#------------------------------------------------------------------------------
# Convert seconds since 1970 to date time
#   - center the day at 12 UT.
#------------------------------------------------------------------------------

def convert_time_seconds(seconds):
    base = dt.datetime(1970,1,1)
    time = base + dt.timedelta(seconds = seconds)
    return time

#------------------------------------------------------------------------------
# read euv file - this determines the wavelength bins
#------------------------------------------------------------------------------

def read_euv_csv_file(file):

    fpin = open(file, 'r')

    iFound = 0
    for line in fpin:
        aline = line.split(',')
        s = aline[-1].strip().split('.')[0]
        if (aline[0].strip() == "Short"):
            if (s.isnumeric()):
                short = np.asarray(aline[5:], dtype = float)
            else:
                short = np.asarray(aline[5:-1], dtype = float)
            iFound += 1
        if (aline[0].strip() == "Long"):
            if (s.isnumeric()):
                long = np.asarray(aline[5:], dtype = float)
            else:
                long = np.asarray(aline[5:-1], dtype = float)
            iFound += 1
        if (iFound == 2):
            break
    wavelengths = {'short': short,
                   'long': long}
    return wavelengths
    
#------------------------------------------------------------------------------
# read fism2 csv file
#------------------------------------------------------------------------------

def read_fism_csv_file(file):

    fpin = open(file, 'r')

    header = fpin.readline()
    vars = header.split(',')

    isSeconds = False

    m = re.match(r'.*seconds.*',vars[0])
    if m:
        isSeconds = True

    # Read in a 1d time array and 1d wavelength array
    # read in all of the irradiance and uncertainty data as
    # a large 1d array, then reshape it into a time vs wavelength 2d array
    # structure of the file is all wavelengths for one time, then
    # a new time and all wavelengths for that time.
    
    iRow = 0
    oldtime = 0.0
    nWaves = 0
    nTimes = 0
    times = []
    wavelengths = []
    irradiance = []
    uncertainty = []
    for line in fpin:
        aline = line.split(',')
        yearday = float(aline[0])
        irradiance.append(float(aline[2]))
        uncertainty.append(float(aline[3]))
        if (iRow == 0):
            oldtime = yearday

        if (yearday > oldtime):
            nWaves = 1
            if isSeconds:
                times.append(convert_time_seconds(oldtime))
            else:
                times.append(convert_time(oldtime))
            oldtime = yearday
            nTimes = nTimes + 1
        else:
            nWaves = nWaves + 1

        if (nTimes == 0):
            wavelengths.append(float(aline[1]))
            
        iRow = iRow+1

    if isSeconds:
        times.append(convert_time_seconds(oldtime))
    else:
        times.append(convert_time(oldtime))
    nTimes = nTimes + 1
    
    irradiance = np.array(irradiance).reshape((nTimes, nWaves))
    uncertainty = np.array(uncertainty).reshape((nTimes, nWaves))
    
    fism_data = {'vars': vars,
                 'time': times,
                 'nTimes': nTimes,
                 'wave': np.array(wavelengths)*10.0,  # convert to Angstroms
                 'irr': irradiance,
                 'unc': uncertainty}
            
    return fism_data

#------------------------------------------------------------------------------
# rebin FISM data into new wavelength bins
#  some bins are single wavelength, and some span lots of wavelenghts
#------------------------------------------------------------------------------

def rebin_fism(fism_waves, fism_vals, wavelengths):

    shorts = wavelengths['short']
    longs = wavelengths['long']
    nWaves = len(shorts)
    new_irr = np.zeros(nWaves)
    ave_wav = np.zeros(nWaves)

    # first go through all of the wavelengths that are singular
    for iWave, short in enumerate(shorts):
        long = longs[iWave]
        if (long == short):
            d = np.abs(fism_waves - short)
            i = np.argmin(d)
            new_irr[iWave] = fism_vals[i] * \
                    ((fism_waves[i+1] - fism_waves[i])/10.0)
            # zero out bin so we don't double count it.
            # fism_vals[i] = 0.0

    # then go through the ranges
    for iWave, short in enumerate(shorts):
        long = longs[iWave]
        ave_wav[iWave] = (short + long)/2.0
        if (long != short):
            d = np.abs(fism_waves - short)
            iStart = np.argmin(d)
            d = np.abs(fism_waves - long)
            iEnd = np.argmin(d)
            wave_int = 0.0
            for i in range(iStart+1, iEnd+1):
                new_irr[iWave] += fism_vals[i] * \
                    ((fism_waves[i+1] - fism_waves[i])/10.0)
                wave_int += (fism_waves[i+1] - fism_waves[i])/10.0
    return new_irr, ave_wav

#------------------------------------------------------------------------------
# main code:
#------------------------------------------------------------------------------

args = parse_args_fism()

if (args.fismfile == 'none'):

    if (args.start == '0'):
        print('Need to specify -start time! Use -h for help!')
        exit()
    
    start = convert_ymdhm_to_dt(args.start)

    if (args.flare):
        print('*** Downloading Flare data - Can only get 24 hours of data! ***')
        end = start + dt.timedelta(days = 1)
    else:
        if (args.end == '0'):
            print('Need to specify -end time! Use -h for help!')
            exit()
        end = convert_ymdhm_to_dt(args.end)

    fism_file = download_fism2(start, end, args.flare)

else:
    fism_file = args.fismfile

fism = read_fism_csv_file(fism_file)
wavelengths = read_euv_csv_file(args.euvfile)

nWaves = len(wavelengths['short'])

filetime = fism["time"][0].strftime('fism%Y%m')
filestart = filetime+'_nWaves_%03d' % nWaves

fig = plt.figure(figsize = (10,10))
ax = fig.add_subplot()

if (args.gitm):
    fileout = filestart + '_gitm.dat'
else:
    fileout = filestart + '.dat'
    
fp = open(fileout, 'wb')

if (args.gitm):
    fp.write('#START\n'.encode())

else:

    shortline = ' 0000 00 00 00 00 00 '
    for short in wavelengths['short']:
        shortline = shortline + "%8.1f" % short
    shortline = shortline + '\n'
    fp.write(shortline.encode())

    longline = ' 0000 00 00 00 00 00 '
    for long in wavelengths['long']:
        longline = longline + "%8.1f" % long
    longline = longline + '\n'
    fp.write(longline.encode())

for iTime, time in enumerate(fism['time']):
    new_irr, ave_wav = rebin_fism(fism['wave'], fism['irr'][iTime], wavelengths)

    if (args.gitm):
        ave_wav = np.flip(ave_wav)
        new_irr = np.flip(new_irr)
    ax.scatter(ave_wav, new_irr)

    sTime = time.strftime(' %Y %m %d %H %M %S')
    sData = ' '
    for irr in new_irr:
        sData = sData + "%15.8e" % irr
    line = sTime + sData + '\n'
    fp.write(line.encode())

fp.close()

ax.set_xlabel('Wavelength (A)')
ax.set_ylabel(fism['vars'][2])

plotfile = filestart + '.png'
print('writing : ',plotfile)    
fig.savefig(plotfile)
plt.close()
