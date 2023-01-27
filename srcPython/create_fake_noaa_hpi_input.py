#!/usr/bin/env python

# Author: Dogacan S. Ozturk

import os
import argparse
import datetime as dt
import numpy as np

try:
    from spacepy.pybats import kyoto
except ModuleNotFoundError:
    print('Spacepy not found. If using SME data, please pass it as a numpy array or a list to calculate_hp_from_ae.')

def calculate_hp_from_ae(ae):
    '''
    HP is in GW.
    AE is in nT.
    Formula taken from Wu et al, 2021. https://doi.org/10.1029/2020SW002629
    '''
    hp = 0.102*ae + 8.953
    return hp

def write_derived_hp(time_array, hp, output_filename = "empty"):

    if (output_filename == "empty"):
    
        savedir = './hp_from_ae'
        if not os.path.exists(savedir):
            os.mkdir(savedir)
        output_filename = savedir + \
            'power_from_ae_{0:%Y%m%d}'.format(time_array[0]) + \
            '_to_{0:%Y%m%d}.txt'.format(time_array[-1])

    output_file = open(output_filename, 'w')

    fmt_line = '{0:%Y-%m-%d} {0:%H:%M:%S} NOAA-17 (N)  6{1:7.2f}   0.75\n'
    ntimes = len(time_array)
    
    output_file.write(':Data_list: '+output_filename+'\n'.format(time_array[0]))
    output_file.write(':Created: {0:%a %b %d %H:%M:%S UTC %Y\n}'.format(dt.datetime.now()))
    output_file.write('# This file is created to replicate NOAA HPI files created before 2013.\n')
    output_file.write('# Please use with caution as Sat number, hemisphere, activity level, and normalizing factors are placeholders.\n')
    output_file.write('#\n')
    output_file.write('# Source: AE or SME index.\n')
    output_file.write('# Units: gigawatts\n\n')
    output_file.write('# Format:\n\n')
    output_file.write('# Each line is formatted as in this example:\n\n')
    output_file.write('# 2006-09-05 00:54:25 NOAA-16 (S)  7  29.67   0.82\n\n')
    output_file.write('# A19   Date and UT at the center of the polar pass as YYYY-MM-DD hh:mm:ss\n')
    output_file.write('# 1X    (Space)\n')
    output_file.write('# A7    NOAA POES Satellite number\n')
    output_file.write('# 1X    (Space)\n')
    output_file.write('# A3    (S) or (N) - hemisphere\n')
    output_file.write('# I3    Hemispheric Power Index (activity level)\n')
    output_file.write('# F7.2  Estimated Hemispheric Power in gigawatts\n')
    output_file.write('# F7.2  Normalizing factor\n\n')
    
    for i in range(ntimes):
        output_file.write(fmt_line.format(time_array[i], hp[i]))
        
    output_file.close()

#__________________________________________________________________________#
#    To use the code and generate a fake HPI file modify the dates below.  #
#__________________________________________________________________________#
if __name__ == '__main__':

    t_start = dt.datetime(2018,8,1)
    t_end = dt.datetime(2018,9,1)

    try:
        ae_data = kyoto.aefetch(t_start, t_end)
        ae = kyoto.KyotoAe(ae_data)
        hp = calculate_hp_from_ae(ae['ae'])
        write_derived_hp(ae['time'], hp)
        
    except AttributeError:
        print('Spacepy Kyoto library not found. For SME derived fake HPI, use standalone functions in this file.')



