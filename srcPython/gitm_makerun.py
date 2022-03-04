#!/usr/bin/env python

from gitm_routines import *
from omniweb import *
from supermag_api import *
from swmf_imf import *
import argparse
import re
from datetime import datetime
from datetime import timedelta
import os

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args(plotTypes):

    parser = argparse.ArgumentParser(description = 'Prepare GITM Run')

    parser.add_argument('-input', \
                        help = 'Use this file as the baseline UAM.in', \
                        default = 'UAM.in')
    
    parser.add_argument('-output', \
                        help = 'Output filename', \
                        default = 'none')
    
    parser.add_argument('start', metavar = 'start', \
                        help = 'start date as YYYYMMDD')
    parser.add_argument('end', metavar = 'end', \
                        help = 'end date as YYYYMMDD')

    parser.add_argument('-imf', \
                        help='set IMF file (find = download IMF from omniweb)')
    
    parser.add_argument('-fism', \
                        help='turn on FISM and point to FISM daily file', \
                        action="store_true")
    
    parser.add_argument('-dynamo', \
                        help='turn on DYNAMO', \
                        action="store_true")
    
    parser.add_argument('-fang', \
                        help='turn on Fang [2010] auroral energy deposition', \
                        action="store_true")

    parser.add_argument('-fta', \
                        help='turn on FTA auroral model', \
                        action="store_true")
    parser.add_argument('-newell', \
                        help='turn on Newell OVATION auroral model', \
                        action="store_true")

    parser.add_argument('-sme', \
                        help = 'SuperMAG AE/AU/AL file (find = cp file)')

    parser.add_argument('-hpi', \
                        help = 'NOAA Hemispheric Power file (find = cp file)')

    parser.add_argument('-datadir', \
                        help = 'Base directory to find data', \
                        default = '~/Data/')

    parser.add_argument('-conduction', nargs = 3, \
                        help = 'set thermal conduction')
    
    parser.add_argument('-nhe', \
                        help = 'set neutral heating efficiency')

    parser.add_argument('-phe', \
                        help = 'set photoelectron heating efficiency')

    # Grid stuff
    parser.add_argument('-nlats', \
                        help = 'set number of latitude blocks')
    parser.add_argument('-nlons', \
                        help = 'set number of longitude blocks')
    parser.add_argument('-latrange', nargs = 2, \
                        help = 'set [min, max] latitude')
    parser.add_argument('-lonrange', nargs = 2, \
                        help = 'set [min, max] longitude (0 0 for whole Earth)')
    
    # Plotting stuff
    for type in plotTypes:
        parser.add_argument('-' + type, \
                            help = 'Set time interval for ' +
                            type.upper() + ' plots (in sec.)')

    parser.add_argument('-3dmag', \
                        help = 'Output one 3DMAG file',
                        action = 'store_true')
    
    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
#
# ----------------------------------------------------------------------

def find_value_at_time(data, key, time):

    minTime = 1e32
    minDist = 1e32
    iPt = -1
    minValue = 1e32
    for i, t in enumerate(data['times']):
        d = np.abs((t-time).total_seconds())
        if (d < minDist):
            minDist = d
            minTime = t
            minValue = data[key][i]
            iPt = i

    return minValue, iPt

# ----------------------------------------------------------------------
#
# ----------------------------------------------------------------------

def find_ave_between_times(data, key, time, window):

    tStart = time + timedelta(seconds = -window/2)
    tEnd = time + timedelta(seconds = window/2)
    
    v, iStart = find_value_at_time(data, key, tStart)
    v, iEnd = find_value_at_time(data, key, tEnd)

    ave = np.mean(data[key][iStart:iEnd])

    return ave
        
#-----------------------------------------------------------------------------
# Reads a UAM.in file into a dictionary
#-----------------------------------------------------------------------------

def read_uam(file):

    fpin = open(file, 'r')
    lines = fpin.readlines()
    fpin.close()
    
    uam = {}

    IsInHash = 0
    iVal = 0
    for line in lines:
        
        line = line.strip('\n')
        smashed = line.replace(" ", "")
        if (len(smashed) == 0):
            IsInHash = 0
        else:
            if (IsInHash):
                values = line.split()
                quant = values[0]
                if (len(values) > 1):
                    desc = ""
                    for v in values[1:]:
                        desc = desc + " " + v
                else:
                    desc = " no descriptor"
                uam[hash][iVal] = [quant, desc]
                iVal = iVal + 1
            else:
                m = re.match(r'^(#.*)', line)
                if m:
                    hash = m.group(1)
                    if (hash.find('END') != 1):
                        IsInHash = 1
                        uam[hash] = {}
                        iVal = 0
        
    return uam

# ----------------------------------------------------------------------
#
# ----------------------------------------------------------------------

def write_uam(file, uam):

    fp = open(file, 'wb')
    
    fp.write("\n".encode())
    
    for key in uam.keys():
        fp.write(key.encode())
        fp.write("\n".encode())
        sub = uam[key]
        for k in sub.keys():
            spaces = '  '
            l = 15 - len(sub[k][0])
            for i in range(0,l):
                spaces = spaces + ' '
            s = sub[k][0]
            fp.write(s.encode())
            fp.write(spaces.encode())
            s = sub[k][1]
            fp.write(s.encode())
            fp.write("\n".encode())
        fp.write("\n".encode())
        
    fp.write("#END\n".encode())
    fp.write("\n".encode())
    fp.close()

    return

# ----------------------------------------------------------------------
#
# ----------------------------------------------------------------------

def download_sme_data(start, end):

    startStr = start.strftime('%Y-%m-%dT%M:%D')
    length = (end - start).total_seconds()
    print('--> Downloading AE data')
    userid = 'ridley'
    (status,idxdata) = SuperMAGGetIndices(userid, start, length, 'all')

    nPts = len(idxdata['tval'])
    print('  --> Found %d points' % nPts)

    if (nPts > 0):
        
        base = datetime(1970,1,1,0,0,0)
        current = base + timedelta(seconds = idxdata['tval'][0])

        ymd = current.strftime('%Y%m%d')
        fileout = 'ae_' + ymd + '.txt'
        print('  --> Writing file ' + fileout)
        
        l0 = 'Created by python code using SuperMAGGetIndices\n'
        l1 = '============================================================\n'
        l2 = '<year>  <month>  <day>  <hour>  <min>  <sec>  '
        l2 = l2 + '<SME (nT)>  <SML (nT)>  <SMU (nT)>\n'

        fp = open(fileout, 'wb')
        fp.write("\n".encode())
        fp.write(l0.encode())
        fp.write("\n".encode())
        fp.write(l1.encode())
        fp.write(l2.encode())

        for i, t in enumerate(idxdata['tval']):
            current = base + timedelta(seconds = t)
            ae = idxdata['SME'][i]
            al = idxdata['SML'][i]
            au = idxdata['SMU'][i]
            out = " %8.2f %8.2f %8.2f" % (ae, al, au)
            ymdhms = current.strftime('%Y  %m  %d  %H  %M  %S')
            line = ymdhms + out + "\n"
            fp.write(line.encode())

        fp.close()

    else:
        
        fileout = 'none.txt'
        print('  --> ERROR!! No AE data downloaded, no file written!!')

    return fileout

# ----------------------------------------------------------------------
#
# ----------------------------------------------------------------------

def read_f107(file):

    f107 = {
        'times': [],
        'f107': [] }
    
    fpin = open(file, 'r')
    lines = fpin.readlines()
    fpin.close()
    lines = lines[15:]

    for line in lines:
        tv = line.split(';')
        f107['f107'].append(float(tv[1]))
        dt = tv[0].split(' ')
        ymd = dt[0].split('-')
        hms = dt[1].split(':')
        y = int(int(ymd[0]))
        m = int(int(ymd[1]))
        d = int(int(ymd[2]))
        h = int(int(hms[0]))
        mi = int(int(hms[1]))
        s = int(int(hms[2]))
        f107['times'].append(datetime(y, m, d, h, mi, s))

    return f107
        
#------------------------------------------------------------------------------
# main code:
#------------------------------------------------------------------------------

plotTypes = ['3dall', '3danc', '3dthm', '3dneu', '3dion',
             '2danc', '2dmel', '2dgel']

args = parse_args(plotTypes)

start = args.start
end = args.end

file = args.input

if (os.path.exists(file)):
    print('--> Reading file : ', file)
    uam = read_uam(file)
else:
    print('Can not find input file ', file, '!')
    print('  Please run with -h if you need help!')
    exit()

# ----------------------------------------
# output types

nPlots = 0
savePlots = []
savePlotsDt = []
for type in plotTypes:
    if (getattr(args, type)):
        savePlots.append(type.upper())
        savePlotsDt.append('%d' % int(getattr(args, type)))
        nPlots = nPlots + 1
if (getattr(args, '3dmag')):
    savePlots.append('3DMAG')
    savePlotsDt.append('604800')
    nPlots = nPlots + 1
    
if (nPlots > 0):
    print('--> Setting save plots:')
    uam['#SAVEPLOTS'] = {
        0 : ['7200', ' dt for restarts'],
        1 : ['%d' % nPlots, ' number of plot types'] }
    print('    7200 -> restarts')
    for iPlot in range(nPlots):
        pt = [savePlots[iPlot], ' plot type']
        dt = [savePlotsDt[iPlot], ' dt for plot type ' + savePlots[iPlot]]
        uam['#SAVEPLOTS'][iPlot*2 + 2] = pt
        uam['#SAVEPLOTS'][iPlot*2 + 3] = dt
        print('    ' + savePlotsDt[iPlot] + ' -> ' + savePlots[iPlot])
        
# ----------------------------------------
# Figure out grid:

if '#GRID' in uam:
    print('--> Interpreting #GRID')
    nLons = int(uam['#GRID'][0][0])
    nLats = int(uam['#GRID'][1][0])
    minLat = float(uam['#GRID'][2][0])
    maxLat = float(uam['#GRID'][3][0])
    minLon = float(uam['#GRID'][4][0])
    maxLon = float(uam['#GRID'][5][0])

if (args.nlons):
    nLons = int(args.nlons)
    print(' --> Setting nLons: ', nLons)

if (args.nlats):
    nLats = int(args.nlats)
    print(' --> Setting nLats: ', nLats)

if (args.latrange):
    minLat = float(args.latrange[0])
    maxLat = float(args.latrange[1])
    print(' --> Setting latitude range : ', minLat, maxLat)
    
if (args.lonrange):
    minLon = float(args.lonrange[0])
    maxLon = float(args.lonrange[1])
    print(' --> Setting longitude range : ', minLon, maxLon)

uam['#GRID'] = {
    0: ['%d' % nLons, ' nBlocksLongitude'],
    1: ['%d' % nLats, ' nBlocksLatitude'],
    2: ['%f' % minLat, ' min latitude to model (deg)'],
    3: ['%f' % maxLat, ' max latitude to model (deg)'],
    4: ['%f' % minLon, ' min longitude to model (deg)'],
    5: ['%f' % maxLon, ' max longitude to model (0 = whole earth) (deg)'] }

# ----------------------------------------
# Grab IMF from OMNIWeb:

if (args.imf):
    print('--> Processing IMF file : ')

    if (args.imf.find('find') == 0):
        results = download_omni_data(start, end, "-all")
        omniDirty = parse_omni_data(results)
        data = clean_omni(omniDirty)

        imfFile = data["times"][0].strftime('imf%Y%m%d.dat')
        message = "Downloaded from OMNIWeb and processed by omniweb_read.py\n"
        write_swmf_imf_file(data, imfFile, message)
        print(' --> Downloaded IMF file : ', imfFile)
    else:
        imfFile = args.imf
        
    if (os.path.exists(imfFile)):
        print(' --> Found IMF file : ', imfFile)
    else:
        print(' --> ERROR!! Could not find IMF file : ', imfFile)
        print('     May not work right!  Check file!')
        
    uam['#MHD_INDICES'] = {
        0: [imfFile, ' imf input file'] }


# ----------------------------------------
# Set Times:

yStart = start[0:4]
m = start[4:6]
d = start[6:8]
start = datetime(int(yStart), int(m), int(d), 0, 0, 0)

uam['#TIMESTART'] = {
    0: [yStart, ' year'],
    1: [m, ' month'],
    2: [d, ' day'],
    3: ['00', ' hour'],
    4: ['00', ' minute'],
    5: ['00', ' second'] }

y = end[0:4]
m = end[4:6]
d = end[6:8]
end = datetime(int(y), int(m), int(d), 0, 0, 0)

uam['#TIMEEND'] = {
    0: [y, ' year'],
    1: [m, ' month'],
    2: [d, ' day'],
    3: ['00', ' hour'],
    4: ['00', ' minute'],
    5: ['00', ' second'] }

print('--> Setting Start Time : ', start)
print('--> Setting End Time : ', end)

# ----------------------------------------
# Deal with F107 stuff:

diff = (end - start).total_seconds()
mid = start + timedelta(seconds = diff/2)

f107file = 'UA/DataIn/f107_adjusted.txt'
f107data = read_f107(f107file)

f107, i = find_value_at_time(f107data, 'f107', mid)
f107a = find_ave_between_times(f107data, 'f107', mid, 81.0*86400.0)

f107s = ('%5.1f' % f107).strip()
f107as = ('%5.1f' % f107a).strip()
uam['#F107'] = {
    0: [f107s, ' data from '+f107file],
    1: [f107as, ' 81 day average of f107'] }

print('--> Setting F107, F107a : ' + f107s + ', ' + f107as)

# ----------------------------------------
# FISM

if (args.fism):
    fismFile = 'UA/DataIn/FISM/fismflux_daily_' + yStart + '.dat'
    if (os.path.exists(fismFile)):
        print('--> Setting fism file : ', fismFile)
        uam['#EUV_DATA'] = {
            0: ['T', ' Use FISM Solar Flux Data (daily)'],
            1: [fismFile, ' FISM Filename'] }
    else:
        print('--> ERROR : fism file does not exist : ', fismFile)
        print('    Setting #EUV_DATA to false!')
        uam['#EUV_DATA'] = {
            0: ['F', ' Use FISM Solar Flux Data (daily)'],
            1: [fismFile, ' FISM Filename'] }

# ----------------------------------------
# HPI file

if (args.hpi):
    # find hemispheric power file:
    if (args.hpi.find('find') == 0):
        hpiFile = args.datadir + '/Hpi/power_' + yStart + '.txt'
        if (os.path.exists(hpiFile)):
            command = 'cp -f ' + hpiFile + ' .'
            print('--> Hemispheric Power File Found')
            print(' --> running command : ' + command)
            os.system(command)
            hpiFile = 'power_' + yStart + '.txt'
    else:
        hpiFile = args.hpi
    if (os.path.exists(hpiFile)):
        print('--> Setting HPI file : ', hpiFile)
        uam['#NOAAHPI_INDICES'] = {
            0: [hpiFile, ' HPI Filename'] }
    else:
        print('--> ERROR : HPI file does not exist : ', hpiFile)
        print('  Taking no action! Check output file!!!')


# ----------------------------------------
# SME file

if (args.sme):
    # find hemispheric power file:
    if (args.sme.find('find') == 0):
        smeFile = download_sme_data(start, end)
        #smeFile = args.datadir + '/Ae/sme_' + yStart + '.txt'
        #if (os.path.exists(smeFile)):
        #    command = 'cp -f ' + smeFile + ' .'
        #    print('--> SuperMAG AE/AU/AL File Found')
        #    print(' --> running command : ' + command)
        #    os.system(command)
        #    smeFile = 'sme_' + yStart + '.txt'
    else:
        smeFile = args.sme
    if (os.path.exists(smeFile)):
        print('--> Setting SME file : ', smeFile)
        uam['#SME_INDICES'] = {
            0: [smeFile, ' SME Filename'],
            1: ['none', ' onset time delay file'] }
    else:
        print('--> ERROR : SME file does not exist : ', smeFile)
        print('  Taking no action! Check output file!!!')

        
# ----------------------------------------
# Thermal Conduction

if (args.conduction):
    print('--> Setting thermal conduction : ', args.conduction)
    uam['#THERMALCONDUCTION'][0][0] = args.conduction[0]
    uam['#THERMALCONDUCTION'][1][0] = args.conduction[1]
    uam['#THERMALCONDUCTION'][2][0] = args.conduction[2]
    
# ----------------------------------------
# Neutral Heating Efficiency

if (args.nhe):
    print('--> Setting neutral heating efficiency : ', args.nhe)
    uam['#NEUTRALHEATING'][0][0] = args.nhe
    
# ----------------------------------------
# Photoelectron Heating Efficiency

if (args.phe):
    print('--> Setting photoelectron heating efficiency : ', args.phe)
    uam['#PHOTOELECTRON'][0][0] = args.phe
    
# ----------------------------------------
# Dynamo

if (args.dynamo):
    # This is not a great way of doing this, but it is something
    if ((minLat < -80.0) and (maxLat > 80.0)):
        uam['#DYNAMO'] = {
            0 : ['T', ' UseDynamo'],
            1 : ['45.0', ' Dynamo magnetic latitude boundary'],
            2 : ['500', ' nIterations'], 
            3 : ['1.0', ' Max Residual'], 
            4 : ['F', ' Include Cowling Conductances'],
            5 : ['40.0', ' Longitude Averaging'] }
        print('--> Turning on Dynamo')
    else:
        print('--> ERROR : dynamo can not be turned on with partial Earth!')
        print('    Setting dynamo to false!')
        uam['#DYNAMO'] = {
            0 : ['F', ' UseDynamo'] }

# ----------------------------------------
# Fang

if (args.fang):
    # turn on fang auroral energy deposition
    uam['#FANGENERGY'] = {
        0 : ['T', ' Use Fang 2010 energy deposition'] }
    
    print('--> Turning on Fang auroral energy deposition')
    
# ----------------------------------------
# FTA model of the aurora

if (args.fta):
    # turn on fang auroral energy deposition
    uam['#FTAMODEL'] = {
        0 : ['T', ' Use FTA model of the aurora'] }
    print('--> Turning on FTA Model')

    # Turn off other aurora
    uam['#NEWELLAURORA'] = {
        0 : ['F', ' Use Newell Ovation model of the aurora'] }
    print(' --> Turning off Newell Auroral Model')

    if (not args.sme):
        print(' --> WARNING : sme not set, so may not work right! Careful!')

# ----------------------------------------
# Newell model of the aurora

if (args.newell):
    # turn on fang auroral energy deposition
    uam['#NEWELLAURORA'] = {
        0 : ['T', ' Use Newell Ovation model of the aurora'],
        1 : ['T', ' Use Diffuse Aurora'],
        2 : ['F', ' Use Monoenergetic Aurora'],
        3 : ['F', ' Use Wave driven Aurora'],
        4 : ['F', ' Remove spikes'],
        5 : ['F', ' Average patterns'] }   
    print('--> Turning on Newell Ovation Model')

    uam['#FTAMODEL'] = {
        0 : ['F', ' Use FTA model of the aurora'] }
    print(' --> Turning off FTA Model')

    if (not args.imf):
        print(' --> WARNING : imf not set, so may not work right! Careful!')

fileout = args.output
if (fileout.find('none') != 1):
    fileout = start.strftime('UAM_%Y%m%d.in')
    
print('--> Writing file : ', fileout)
write_uam(fileout, uam)

