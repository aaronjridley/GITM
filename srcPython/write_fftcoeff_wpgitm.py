#!/usr/bin/env python 

'''
PROGRAM:
    write_fftcoeff_wpgitm.py

PURPOSE:
    for a given time series of surface vertical displacement or vertical velocity,
    obtain FFT coefficients and write into a file for input to Wave Perturbation(WP)-GITM. 
    This is for specifying a multi-frequency wave perturbation to represent surface motions
    during tsunamis and earthquakes.    

USAGE: 
    write_fftcoeff_wpgitm.py -i <input_file> -b <time_begin> -e <time_end>
    format for time_begin and time_end: yyyymmddHHMMSS
    if time_begin and/or time_end not provided, perform FFT on the entire
    time series.
    example: 
    python write_fftcoeff_wpgitm.py -i test_pfrj.dat -b 20150916225505 -e 20150916225824
 
    Format of the time series file:
        # an arbitrary number of header lines
        year month day hour minute second (millisecond) data  


AUTHOR:
    Xing Meng (Xing.Meng@jpl.nasa.gov)

LICENSE:
    Copyright 2021 California Institute of Technology
    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0

'''

import sys, getopt
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt


def main(argv):

    # obtain input file name and begin/end time from command line
    input_file = ''
    time_begin = ''
    time_end = ''
    try:
        opts, args = getopt.getopt(argv,"hi:b:e:")
    except getopt.GetoptError:
        print 'write_fftcoeff_wpgitm.py -i <input_file> -b <time_begin> -e <time_end>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'write_fftcoeff_wpgitm.py -i <input_file> -b <time_begin> -e <time_end>'
            print 'format for time_begin and time_end: yyyymmddHHMMSS'         
            sys.exit()
        elif opt in ("-i"):
            input_file = arg
        elif opt in ("-b"):
            time_begin = datetime.strptime(arg,"%Y%m%d%H%M%S")
        elif opt in ("-e"):
            time_end = datetime.strptime(arg,"%Y%m%d%H%M%S")

    time = np.array([])
    signal = np.array([])
    f = open(input_file)
    for line in f:        
        if line[0] == '#':
            # skip header lines begin with '#'
            continue
        record = line.split()
        if len(record) == 8:
            # when the file contains milliseconds
            time_current = datetime.strptime(' '.join(record[0:7]),"%Y %m %d %H %M %S %f")
        elif len(record) == 7:
            time_current = datetime.strptime(' '.join(record[0:6]),"%Y %m %d %H %M %S")
        if time_begin == '' or time_current >= time_begin:
            if time_end != '' and time_current > time_end:
                break
            time = np.append(time, time_current)
            signal = np.append(signal, float(record[-1]))                        
    f.close()
    ndata = len(time)
    timestep = (time[1]-time[0]).seconds
    print 'peform FFT for time interval ', time[0], ' to ', time[-1]
    print 'signal temporal resolution = ', timestep, ' seconds' 
    print 'number of data points = ', ndata
    
    # peform FFT
    fftcoeff = np.fft.fft(signal)/ndata
    # natural frequency
    nu = np.fft.fftfreq(ndata,timestep)
    # write coefficients for the first half frequencies into file
    nfreq = ndata/2  
    output_file = 'FFTcoeff_'+input_file
    with open(output_file,'w') as coeff_file:
        coeff_file.write('FFT coefficients for '+input_file+'\n')
        coeff_file.write('Start time: ' + datetime.strftime(time[0],"%Y-%m-%d %H:%M:%S") \
                         +', End time: ' + datetime.strftime(time[-1],"%Y-%m-%d %H:%M:%S") \
                         + '\n') 
        coeff_file.write('Frequency(Hz)    Coeff_real   Coeff_imag\n')
        coeff_file.write('#START\n')
        for ifreq in range(0,nfreq):
            coeff_file.write('{:10.5f}{:18.8E}{:18.8E}\n'.\
                             format(nu[ifreq], np.real(fftcoeff[ifreq]),\
                                    np.imag(fftcoeff[ifreq])))

    # reconstruct the signal using sine and cosine waves to verify the signal
    # reconstruction in WP.
    time_recons = np.arange(ndata)*timestep
    signal_recons = np.real(fftcoeff[0])*np.cos(2*np.pi*0)
    for ifreq in range(1,nfreq):
        if nu[ifreq] <=  0.1:
            # WP only takes frequencies <= 0.1Hz
            signal_recons = signal_recons + \
                (np.real(fftcoeff[ifreq])*np.cos(-2*np.pi*nu[ifreq]*time_recons) \
                 + np.imag(fftcoeff[ifreq])*np.sin(-2*np.pi*nu[ifreq]*time_recons))*2

    # plot the original signal and reconstructed signal
    plt.figure()
    plt.plot(time, signal, 'k', label='original signal')
    plt.plot(time, signal_recons, 'r', label='reconstructed with frequencies <= 0.1Hz)')
    plt.legend()
    plt.xlabel('Time')
    plt.title('Time series reconstruction for '+input_file)
    plt.savefig('signal_reconstruct.pdf')
    plt.close()

if __name__ == "__main__":
   main(sys.argv[1:])

