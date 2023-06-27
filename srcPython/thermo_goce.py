#!/usr/bin/env python

from datetime import datetime
from datetime import timedelta
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from gitm_routines import *
from pylab import cm
import argparse

#-----------------------------------------------------------------------------
# get arguments from the user
#-----------------------------------------------------------------------------

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Compare GITM model results with GOCE or CHAMP')

    parser.add_argument('-winds',  \
                        action='store_true', default = False, \
                        help = 'plot winds instead of density')
    
    args = parser.parse_args()

    return args
    

#-----------------------------------------------------------------------------
# Plot a lat vs time contour
#-----------------------------------------------------------------------------

def plot_map(mapTime, mapLats, satMap, modelMap,
             mini, maxi, cmap,
             sTime,
             maxLat, StartTime,
             satellite, short_var,
             var, units, fac,
             nodeTitle, nodeValue):
        
    fig = plt.figure(figsize = (10,10))

    ax = fig.add_axes([0.075, 0.53, 1.0, 0.43])
    cax = ax.pcolor(mapTime, mapLats, satMap, \
                    vmin = mini, vmax = maxi, cmap = cmap)
    ax.set_xlabel(sTime)
    ax.set_ylabel('Latitude (deg)')
    ax.set_ylim([-maxLat, maxLat])
    ax.set_xlim([np.min(mapTime), np.max(mapTime)])
    title = satellite + " (top) vs GITM (bot) - "
    title = title + nodeTitle + ' node (LT  : %4.1f hours)' % nodeValue
    ax.set_title(title)
    cbar = fig.colorbar(cax, ax=ax, shrink = 0.75, pad=0.02)
    facString = '(%7.1e' % fac
    label = satellite + ' ' + var + facString + ')'
    cbar.set_label(label, rotation=90)

    ax2 = fig.add_axes([0.075, 0.05, 1.0, 0.43])
    cax2 = ax2.pcolor(mapTime, mapLats, modelMap, \
                      vmin = mini, vmax = maxi, cmap = cmap)
    ax2.set_xlabel(sTime)
    ax2.set_ylabel('Latitude (deg)')
    ax.set_ylim([-maxLat, maxLat])
    ax2.set_xlim([np.min(mapTime), np.max(mapTime)])
    cbar2 = fig.colorbar(cax2, ax=ax2, shrink = 0.75, pad=0.02)
    label = 'GITM' + ' ' + var + facString + ')'
    cbar2.set_label(label, rotation=90)

    outfile = StartTime.strftime(short_var + '_' + nodeTitle +'_%Y%m%d.png')
    print('--> Writing file : ' + outfile)
    plt.savefig(outfile)
    plt.close()
    
    return

#------------------------------------------------------------------------------
# Calculate the Wind direction given the flight path of the satellite
#------------------------------------------------------------------------------

def calc_wind_dir(lons, lats):

    dlats = lats[1:] - lats[:-1]
    dlats = np.concatenate((dlats, [dlats[-1]]))

    dlons = lons[1:] - lons[:-1]
    dlons = np.concatenate((dlons, [dlons[-1]]))

    # Longitude can go across the 0 - 360 or 360 - 0 boundary, so
    # we need to correct for this possibility:
    
    dlons[dlons > 180.0] = dlons[dlons > 180.0] - 360.0
    dlons[dlons < -180.0] = dlons[dlons < -180.0] + 360.0

    # Longitudes get closer together near the poles, so we need to
    # correct for that also:
    dlons = dlons * np.cos(lats * np.pi / 180.0)

    # Make a unit vector of the direction of travel:
    mag = np.sqrt(dlats**2 + dlons**2)
    unitn = dlats / mag 
    unite = dlons / mag

    # Rotate the unit vector, so that it points orthogonal to the
    # orbit plane.  This will be the actual wind vector direction:
    uniteRot = unitn
    unitnRot = -unite

    return uniteRot, unitnRot

#-----------------------------------------------------------------------------
# Read GOCE data
#-----------------------------------------------------------------------------

def read_goce(file):

    data = {}
    data["times"] = []
    data["alts"] = []
    data["lats"] = []
    data["lons"] = []
    data["lst"] = []
    data["density"] = []
    data["Ve"] = []
    data["Vn"] = []
    data["Vr"] = []
    data["densityError"] = []
    data["windError"] = []
    data["FlagOver"] = []
    data["FlagEclipse"] = []
    data["FlagAD"] = []
    data["FlagThuster"] = []

    f = open(file, 'r')

    for line in f:

        if (line.find('#') < 0):
            items = line.split()
            ymd = items[0].split('-')
            hms = items[1].split(':')
            s = float(hms[2])
            data["times"].append(datetime(int(ymd[0]),int(ymd[1]),int(ymd[2]),
                                         int(hms[0]),int(hms[1]),int(s)))
            data["alts"].append(float(items[3])/1000.0)
            data["lons"].append(float(items[4]))
            data["lats"].append(float(items[5]))
            data["lst"].append(float(items[6]))
            data["density"].append(float(items[8]))
            data["Ve"].append(float(items[9]))
            data["Vn"].append(float(items[10]))
            data["Vr"].append(float(items[11]))
            data["densityError"].append(float(items[12]))
            data["windError"].append(float(items[13]))
            data["FlagOver"].append(int(items[14]))
            data["FlagEclipse"].append(int(items[15]))
            data["FlagAD"].append(int(items[16]))
            data["FlagThuster"].append(int(items[17]))

    f.close()

    # Here we are calculating the direction of travel of the sat:

    uniteRot, unitnRot = calc_wind_dir(np.array(data["lons"]),
                                       np.array(data["lats"]))

    data["wind_e_dir"] = uniteRot
    data["wind_n_dir"] = unitnRot
    
    return data

#-----------------------------------------------------------------------------
# Read CHAMP data
#-----------------------------------------------------------------------------

def read_champ(file):

    data = {}
    data["times"] = []
    data["alts"] = []
    data["lats"] = []
    data["lons"] = []
    data["lst"] = []
    data["density"] = []
    data["Ue"] = []
    data["Un"] = []
    data["Ur"] = []
    data["densityError"] = []
    data["windError"] = []
    data["FlagOver"] = []
    data["FlagEclipse"] = []
    data["FlagAD"] = []
    data["FlagThuster"] = []

    f = open(file, 'r')

    for line in f:

        if (line.find('#') < 0):
            items = line.split()
            ymd = items[0].split('-')
            hms = items[1].split(':')
            s = float(hms[2])
            data["times"].append(datetime(int(ymd[0]),int(ymd[1]),int(ymd[2]),
                                         int(hms[0]),int(hms[1]),int(s)))
            data["alts"].append(float(items[3])/1000.0)
            data["lons"].append( (float(items[4])+360.0) % 360.0 )
            data["lats"].append(float(items[5]))
            data["lst"].append(float(items[6]))
            data["density"].append(float(items[8]))

    f.close()

    return data


#------------------------------------------------------------------------------
def read_champ_winds(file):

    data = {}
    data["times"] = []
    data["alts"] = []
    data["lats"] = []
    data["lons"] = []
    data["lst"] = []
    data["Ve"] = []
    data["Vn"] = []
    data["Vr"] = []
    data["density"] = []

    f = open(file, 'r')

    for line in f:

        if (line.find('#') < 0):
            items = line.split()
            ymd = items[0].split('-')
            hms = items[1].split(':')
            s = float(hms[2])
            data["times"].append(datetime(int(ymd[0]),int(ymd[1]),int(ymd[2]),
                                         int(hms[0]),int(hms[1]),int(s)))
            data["alts"].append(float(items[3])/1000.0)
            data["lons"].append((float(items[4]) + 360.0) % 360.0)
            data["lats"].append(float(items[5]))
            data["lst"].append(float(items[6]))
            data["Ve"].append(float(items[8]))
            data["Vn"].append(float(items[9]))

    f.close()

    # Here we are calculating the direction of travel of the sat:

    uniteRot, unitnRot = calc_wind_dir(np.array(data["lons"]),
                                       np.array(data["lats"]))

    data["wind_e_dir"] = uniteRot
    data["wind_n_dir"] = unitnRot

    return data


#-----------------------------------------------------------------------------
# finds the index of the time within a time array 
#-----------------------------------------------------------------------------

def find_index(time, t):

    if (t < time[0]):
        return 0
    if (t > time[-1]):
        return len(time)-1

    iLow = 0
    iHigh = len(time)
    iMid = int((iHigh + iLow)/2)
    while (iHigh - iLow > 1):
        if (time[iMid] == t):
            iHigh = iMid
            iLow = iMid
        else:
            if (t > time[iMid]):
                iLow = iMid
            else:
                iHigh = iMid
            iMid = int((iHigh + iLow)/2)
    return iMid
            
#-----------------------------------------------------------------------------
# smooths the data over a window of time
#-----------------------------------------------------------------------------

def smooth_data(time, data, window):

    smoothed = np.zeros(len(time))
    for i, t in enumerate(time):
        iMin = find_index(time, t-window/2)
        iMax = find_index(time, t+window/2)
        s = np.mean(data[iMin : iMax+1])
        smoothed[i] = s

    return smoothed
    

def make_map(ix, iy, values, lats, lsts, nX, nY):

    asc = 0.0
    ascN = 0
    des = 0.0
    desN = 0
    nPts = len(ix)
    dLat = lats[1:] - lats[:-1]

    mapA = np.zeros([nX, nY]) + np.nan
    mapD = np.zeros([nX, nY]) + np.nan
        
    for iPt in range(nPts-1):
        i = int(ix[iPt])
        j = int(iy[iPt])
        if (dLat[iPt] > 0):
            mapA[i][j] = values[iPt]
            if (np.abs(lats[iPt]) < 30.0):
                asc = asc + lsts[iPt]
                ascN = ascN + 1
        else:
            mapD[i][j] = values[iPt]
            if (np.abs(lats[iPt]) < 30.0):
                des = des + lsts[iPt]
                desN = desN + 1
        
    asc = asc / ascN
    des = des / desN

    return mapA, mapD, asc, des
    


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Main Code!
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

# define some constants
rtod =  180.0 / np.pi
dtor = 1.0 / rtod

# Get the input arguments
args = get_args()

headers = read_gitm_headers('3DALL')
nGitmFiles = len(headers["time"])

# get lons, lats, alts:
vars = [0,1,2]
data = read_gitm_one_file(headers["filename"][0], vars)
Alts = data[2][0][0]/1000.0;
Lons = data[0][:,0,0]*rtod;
Lats = data[1][0,:,0]*rtod;
[nLons, nLats, nAlts] = data[0].shape

dLon = Lons[1]-Lons[0]
dLat = Lats[1]-Lats[0]

goceDir_lap = '/Users/ridley/Data/Goce/timeseries_data/'
goceDir_zed = '/backup/Data/GOCE/'
champDir_mia = '/backup/Data/Champ/Originals/thermosphere.tudelft.nl/'

IsGoce = True
found = False
if (os.path.isdir(goceDir_lap)):
    goceDir = goceDir_lap
    found = True
if (os.path.isdir(goceDir_zed)):
    goceDir = goceDir_zed
    found = True
if (not found):
    print("Can't seem to find GOCE directory...")
    print("Checking for CHAMP...")
    if (os.path.isdir(champDir_mia)):
        goceDir = champDir_mia
        IsGoce = False
        found = True
    else:
        print("Nope. Shoot.")
        exit()

IsWind = args.winds

if (IsGoce):
    sat = "goce"
    goceFileFront = 'goce_denswind_ac082_v2_0_'
    goceYear = headers["time"][int(nGitmFiles/2)].strftime('%Y-%m.txt')
else:
    sat = "champ"
    if (IsWind):
        goceFileFront = 'CH_WND_ACC_'
    else:
        goceFileFront = 'CH_DNS_ACC_'
    goceYear = headers["time"][int(nGitmFiles/2)].strftime('%Y_%m_v01.txt')
    
goceFile = goceDir + goceFileFront + goceYear

if (not os.path.exists(goceFile)):
    print("Can not find GOCE file : ", goceFile)
    exit()

if (IsGoce):
    print("Reading GOCE file: ", goceFile)
    data = read_goce(goceFile)
else:
    if (IsWind):
        data = read_champ_winds(goceFile)
    else:
        data = read_champ(goceFile)
    
# extract winds -------->>>>
if (IsWind):
    Ve = np.array(data["Ve"])
    Vn = np.array(data["Vn"])
    eDir = data["wind_e_dir"]
    nDir = data["wind_n_dir"]
else:
    rho = np.array(data["density"])
    
goceAltsAll = np.array(data["alts"])
goceLonsAll = np.array(data["lons"])
goceLatsAll = np.array(data["lats"]) 
goceLstAll = np.array(data["lst"]) 
goceTimesAll = np.array(data["times"])

re = 6378.137
rp = 6356.75
diff = re - rp
rgitm = 6372.0
r = rp + diff * np.cos(goceLatsAll * dtor)
correction = r - rgitm
goceAltsPrime = goceAltsAll + correction

meanAlt = np.mean(goceAltsPrime)

goceRm = (r + goceAltsAll) * 1000.0
meanR = np.mean(goceRm)

G = 6.67259e-11
Me = 5.9722e24
mu = G * Me
v0 = np.sqrt(mu / meanR)
period = 2.0 * np.pi * meanR / v0

print("Period : ", period)

iBefore0 = -1
iAfter0 = -1
# here we just assume that we will start with the first set of files.
# this shouldn't matter...
iBefore = 0
iAfter = 1

if (IsWind):
    vars = [16, 17]
else:
    vars = [3]

nVars = len(vars)

AfterVals = np.zeros(nVars)
BeforeVals = np.zeros(nVars)

gitmRho = []
goceRho = []
goceTime = []
goceLats = []
goceLst = []

gitmVn = []
goceVn = []

gitmVe = []
goceVe = []

for i, time in enumerate(goceTimesAll):

    while (time > headers["time"][iAfter]):
        iAfter = iAfter+1
        if (iAfter >= nGitmFiles-1):
            break
        
    if (iAfter == nGitmFiles):
        break

    iBefore = iAfter-1

    if (iBefore != iBefore0):
        file = headers["filename"][iBefore]
        BeforeData = read_gitm_one_file(file, vars)
    if (iAfter != iAfter0):
        file = headers["filename"][iAfter]
        AfterData = read_gitm_one_file(file, vars)

    if (time >= headers["time"][iBefore]):
            
        dt = (headers["time"][iAfter] - \
              headers["time"][iBefore]).total_seconds()
        xt = (time - headers["time"][iBefore]).total_seconds() / dt

        lon = goceLonsAll[i]
        lat = goceLatsAll[i]
        alt = goceAltsPrime[i]
        if (IsWind):
            goceVe.append(Ve[i])
            goceVn.append(Vn[i])
        else:
            goceRho.append(rho[i])
        goceTime.append(time)
        goceLats.append(lat)
        goceLst.append(goceLstAll[i])
        
        xLon = (lon-Lons[0])/dLon
        iLon = int(xLon)
        xLon = xLon - iLon
        
        yLat = (lat-Lats[0])/dLat
        jLat = int(yLat)
        yLat = yLat - jLat

        kAlt = 0
        zAlt = 0.0
        if ((alt > Alts[0]) and (nAlts > 1)):
            if (alt > Alts[nAlts-1]):
                # above domain:
                kAlt = nAlts-2
                zAlt = 1.0
            else:
                while (Alts[kAlt] < alt):
                    kAlt = kAlt + 1
                kAlt = kAlt - 1
                zAlt = (alt - Alts[kAlt]) / (Alts[kAlt+1] - Alts[kAlt])
            kAltp1 = kAlt + 1
        else:
            kAltp1 = kAlt

        for i, v in enumerate(vars):
            BeforeVals[i] = \
                (1-xLon)*(1-yLat)*(1-zAlt)*BeforeData[v][iLon][jLat][kAlt]+\
                (  xLon)*(1-yLat)*(1-zAlt)*BeforeData[v][iLon+1][jLat][kAlt]+\
                (1-xLon)*(  yLat)*(1-zAlt)*BeforeData[v][iLon][jLat+1][kAlt]+\
                (  xLon)*(  yLat)*(1-zAlt)*BeforeData[v][iLon+1][jLat+1][kAlt]+\
                (1-xLon)*(1-yLat)*(  zAlt)*BeforeData[v][iLon][jLat][kAltp1]+\
                (  xLon)*(1-yLat)*(  zAlt)*BeforeData[v][iLon+1][jLat][kAltp1]+\
                (1-xLon)*(  yLat)*(  zAlt)*BeforeData[v][iLon][jLat+1][kAltp1]+\
                (  xLon)*(  yLat)*(  zAlt)*BeforeData[v][iLon+1][jLat+1][kAltp1]
            AfterVals[i] = \
                (1-xLon)*(1-yLat)*(1-zAlt)*AfterData[v][iLon][jLat][kAlt]+\
                (  xLon)*(1-yLat)*(1-zAlt)*AfterData[v][iLon+1][jLat][kAlt]+\
                (1-xLon)*(  yLat)*(1-zAlt)*AfterData[v][iLon][jLat+1][kAlt]+\
                (  xLon)*(  yLat)*(1-zAlt)*AfterData[v][iLon+1][jLat+1][kAlt]+\
                (1-xLon)*(1-yLat)*(  zAlt)*AfterData[v][iLon][jLat][kAltp1]+\
                (  xLon)*(1-yLat)*(  zAlt)*AfterData[v][iLon+1][jLat][kAltp1]+\
                (1-xLon)*(  yLat)*(  zAlt)*AfterData[v][iLon][jLat+1][kAltp1]+\
                (  xLon)*(  yLat)*(  zAlt)*AfterData[v][iLon+1][jLat+1][kAltp1]
            

        if (IsWind):
            i = 0
            ve = ((1-xt) * BeforeVals[i] + xt * AfterVals[i]) * eDir[i]
            i = 1
            vn = ((1-xt) * BeforeVals[i] + xt * AfterVals[i]) * nDir[i]
            gitmVe.append(ve)
            gitmVn.append(vn)
        else:
            i = 0
            gitmRho.append((1-xt) * BeforeVals[i] + xt * AfterVals[i])
        
        #date = time.strftime('%Y %m %d %H %M %S 00 ')
        #pos = '%7.2f %7.2f %8.2f' % (lon, lat, alt)
        #vals = ''
        #for i, v in enumerate(vars):
        #    v = (1-xt) * BeforeVals[i] + xt * AfterVals[i]
        #    vals = vals + '  %e' % v
        ##fpout.write(date+pos+vals+'\n')

    iBefore0 = iBefore
    iAfter0 = iAfter


dt = ((goceTime[-1] - goceTime[0]).total_seconds())/3600.0
if (dt < 2*86400.0):
    StartTime = datetime(goceTime[0].year, \
                         goceTime[0].month, \
                         goceTime[0].day, \
                         goceTime[0].hour)
else:
    StartTime = datetime(goceTime[0].year, \
                         goceTime[0].month, \
                         goceTime[0].day)
    
sTime = StartTime.strftime('%b %d, %Y %H:%M UT') + ' - ' + \
    goceTime[-1].strftime('%b %d, %Y %H:%M UT (Hours)')

tInS = []
for t in goceTime:
    tInS.append((t - StartTime).total_seconds())
tInS = np.array(tInS)

EndTimeS = tInS[-1]
if (IsGoce):
    dLat = 1.0
else:
    dLat = 2.5
nX = int(np.round(EndTimeS/period))+1
nY = int(180.0/dLat)

mapTime = []
mapTime1D = []
mapLats = np.zeros([nX+1, nY+1])

for i in range(nX+1):
    mapLats[i,:] = np.arange(nY+1) * dLat - 90.0 #+ dLat/2.0
    dt = i * period - period/2
    temp = [ dt ] * (nY+1)
    mapTime.append(temp)
    mapTime1D.append(StartTime + timedelta(seconds = dt))
    
nPts = len(goceTime)

goceLats = np.array(goceLats)
mapTime = np.array(mapTime)/3600.0

ix = np.round(tInS/period)
iy = np.round((goceLats + 90.0)/dLat)

dLat = goceLats[1:] - goceLats[0:-1]

if (IsWind):
    gitmRho = gitmVe
    goceRho = goceVe

all = [gitmRho, goceRho]
maxi = np.max(all)
mini = np.min(all)

if (mini >= 0):
    cmap = cm.plasma
    mini = 0.0
else:
    cmap = cm.bwr
    mini = -maxi

if (IsWind):
    fac = 1.0
else:
    fac = 10**(int(np.log10(maxi))-1)

gitmRho = np.array(gitmRho) / fac
goceRho = np.array(goceRho) / fac

mini = mini / fac
maxi = maxi / fac

norm = cm.colors.Normalize(vmax=mini, vmin=maxi)

if (IsGoce):
    maxLat = 82.0
else:
    maxLat = 90.0


if (IsWind):

    gitmVeMapa, gitmVeMapd, asc, dec = make_map(ix, iy, gitmVe, \
                                                  goceLats, goceLst, nX, nY)
    goceVeMapa, goceVeMapd, asc, dec = make_map(ix, iy, goceVe, \
                                                  goceLats, goceLst, nX, nY)

    gitmVnMapa, gitmVnMapd, asc, dec = make_map(ix, iy, gitmVn, \
                                                  goceLats, goceLst, nX, nY)
    goceVnMapa, goceVnMapd, asc, dec = make_map(ix, iy, goceVn, \
                                                  goceLats, goceLst, nX, nY)

    plot_map(mapTime, mapLats, goceVeMapa, gitmVeMapa, mini, maxi, cmap,
             sTime, maxLat, StartTime, sat, 've', 'Zonal Vel.', 'm/s', fac,
             'ascending', asc)
    plot_map(mapTime, mapLats, goceVeMapd, gitmVeMapd, mini, maxi, cmap,
             sTime, maxLat, StartTime, sat, 've', 'Zonal Vel.', 'm/s', fac,
             'descending', dec)
    plot_map(mapTime, mapLats, goceVnMapa, gitmVnMapa, mini, maxi, cmap,
             sTime, maxLat, StartTime, sat, 'vn', 'Merid. Vel.', 'm/s', fac,
             'ascending', asc)
    plot_map(mapTime, mapLats, goceVnMapd, gitmVnMapd, mini, maxi, cmap,
             sTime, maxLat, StartTime, sat, 'vn', 'Merid. Vel.', 'm/s', fac,
             'descending', dec)

else:    

    gitmRhoMapa, gitmRhoMapd, asc, dec = make_map(ix, iy, gitmRho, \
                                                  goceLats, goceLst, nX, nY)
    goceRhoMapa, goceRhoMapd, asc, dec = make_map(ix, iy, goceRho, \
                                                  goceLats, goceLst, nX, nY)

    plot_map(mapTime, mapLats, goceRhoMapa, gitmRhoMapa, mini, maxi, cmap,
             sTime, maxLat, StartTime, sat, 'rho', 'Rho', 'kg/m3', fac,
             'ascending', asc)
    plot_map(mapTime, mapLats, goceRhoMapd, gitmRhoMapd, mini, maxi, cmap,
             sTime, maxLat, StartTime, sat, 'rho', 'Rho', 'kg/m3', fac,
             'descending', dec)

    
    # ----------------------------------------------------------------------
    # Line plots

    goceRhoSmoothed = smooth_data(tInS, goceRho, period)
    gitmRhoSmoothed = smooth_data(tInS, gitmRho, period)

    package = {'times': goceTime,
               'alt': meanAlt,
               'GOCE_Rho': goceRho * fac,
               'GITM_Rho': gitmRho * fac,
               'GOCE_Rho_Smoothed': goceRhoSmoothed * fac,
               'GITM_Rho_Smoothed': gitmRhoSmoothed * fac}
    
    write_log(package, fileHeader = 'goce_rho',
              message = 'From thermo_goce.py')
    
    sTime = goceTime[0].strftime('%b %d, %Y %H:%M UT') + ' - ' + \
        goceTime[-1].strftime('%b %d, %Y %H:%M UT')

    fig = plt.figure(figsize = (10,10))
    ax = fig.add_axes([0.075, 0.53, 0.88, 0.43])
    ax.plot(goceTime, gitmRho, 'b', alpha = 0.1)
    ax.plot(goceTime, gitmRhoSmoothed, 'b', label = 'GITM')
    ax.plot(goceTime, goceRho, 'r', alpha = 0.1)
    ax.plot(goceTime, goceRhoSmoothed, 'r', label = 'GOCE')
    ax.set_ylabel('Rho (kg/m3)')
    ax.set_xlim(goceTime[0], goceTime[-1])
    ax.legend()

    ax2 = fig.add_axes([0.075, 0.05, 0.88, 0.43])
    rawDiff = (np.array(gitmRho) - np.array(goceRho)) / \
        np.array(goceRho) * 100.0
    smoothedDiff = (np.array(gitmRhoSmoothed) - np.array(goceRhoSmoothed)) / \
        np.array(goceRhoSmoothed) * 100.0
    ax2.plot(goceTime, rawDiff, 'k', alpha = 0.1)
    ax2.plot(goceTime, smoothedDiff, 'k')
    ax2.plot(goceTime, rawDiff*0.0, 'g:')
    ax2.set_ylabel('GITM - GOCE Rho Diff (%)')
    ax2.set_xlim(goceTime[0], goceTime[-1])
    ax2.set_xlabel(sTime)

    outfile = StartTime.strftime('goce_gitm_rho_line_%Y%m%d.png')
    print('--> Writing file : ' + outfile)
    plt.savefig(outfile)
    plt.close()

exit()


                                        
