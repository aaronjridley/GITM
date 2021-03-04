#!/usr/bin/env python

from glob import glob
from datetime import datetime
from datetime import timedelta
from struct import unpack
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from matplotlib.gridspec import GridSpec
from pylab import cm
from gitm_routines import *
from aether_routines import *
import re
import sys

rtod = 180.0/3.141592

#-----------------------------------------------------------------------------
# 
#-----------------------------------------------------------------------------

def get_args(argv):

    filelist = []
    IsLog = 0
    var = 15
    alt = 400.0
    lon = -100.0
    lat = -100.0
    tec = 0
    cut = 'alt'
    help = 0
    winds = 0
    diff = 0
    IsGitm = 1
    
    for arg in argv:

        IsFound = 0

        if (not IsFound):

            m = re.match(r'-var=(.*)',arg)
            if m:
                var = int(m.group(1))
                IsFound = 1

            m = re.match(r'-diff',arg)
            if m:
                diff = 1
                IsFound = 1

            m = re.match(r'-tec',arg)
            if m:
                var = 34
                tec = 1
                IsFound = 1

            m = re.match(r'-alt=(.*)',arg)
            if m:
                alt = int(m.group(1))
                IsFound = 1

            m = re.match(r'-lat=(.*)',arg)
            if m:
                lat = int(m.group(1))
                IsFound = 1

            m = re.match(r'-lon=(.*)',arg)
            if m:
                lon = int(m.group(1))
                IsFound = 1

            m = re.match(r'-cut=(.*)',arg)
            if m:
                cut = m.group(1)
                IsFound = 1

            m = re.match(r'-alog',arg)
            if m:
                IsLog = 1
                IsFound = 1

            m = re.match(r'-h',arg)
            if m:
                help = 1
                IsFound = 1

            m = re.match(r'-wind',arg)
            if m:
                winds = 1
                IsFound = 1

            if IsFound==0 and not(arg==argv[0]):
                filelist.append(arg)
                m = re.match(r'.bin',arg)
                if m:
                    IsGitm = 1
                else:
                    IsGitm = 0

                    
    args = {'filelist':filelist,
            'IsGitm':IsGitm,
            'var':var,
            'cut':cut,
            'diff':diff,
            'tec':tec,
            'help':help,
            'winds':winds,
            'alt':alt,
            'lat':lat,
            'lon':lon,
            'IsLog':IsLog}

    return args

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Main Code!
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

args = get_args(sys.argv)

IsGitm = args['IsGitm']

if (IsGitm):
    header = read_gitm_header(args["filelist"])
else:
    header = read_aether_header(args["filelist"])

if (args["var"] >= header["nVars"]):
    print('You asked for a variable that doesnt exist!!!')
    args["help"] = 1
    
if (args["help"]):

    print('Usage : ')
    print('gitm_plot_one_alt.py -var=N -tec -winds -cut=alt,lat,lon') 
    print('                     -alt=alt -lat=lat -lon=lon -alog ')
    print('                     -help [*.bin or a file]')
    print('   -help : print this message')
    print('   -var=number : number is variable to plot')
    print('   -cut=alt,lat,lon : which cut you would like')
    print('   -alt=altitude : can be either alt in km or grid number (closest)')
    print('   -lat=latitude : latitude in degrees (closest)')
    print('   -lon=longitude: longitude in degrees (closest)')
    print('   -alog : plot the log of the variable')
    print('   -winds: overplot winds')
    print('   At end, list the files you want to plot')
    print('   This code should work with GITM files (*.bin) and')
    print('   Aether netCDF files (*.nc)')
    
    iVar = 0
    for var in header["vars"]:
        print(iVar,var)
        iVar=iVar+1

    exit()

        
filelist = args["filelist"]
nFiles = len(filelist)

cut = args["cut"]

vars = [0,1,2]
vars.append(args["var"])

if (args["winds"]):
    if (cut=='alt'):
        iUx_ = 16
        iUy_ = 17
    if (cut=='lat'):
        iUx_ = 16
        iUy_ = 18
    if (cut=='lon'):
        iUx_ = 17
        iUy_ = 18
    vars.append(iUx_)
    vars.append(iUy_)
    AllWindsX = []
    AllWindsY = []

Var = header["vars"][args["var"]]

iVar_ = args["var"]

AllData2D = []
AllAlts = []
AllTimes = []

j = 0
for file in filelist:

    if (IsGitm):
        data = read_gitm_one_file(file, vars)
        iVar_ = args["var"]
    else:
        if (j == 0):
            VarList = []
            for v in vars:
                VarList.append(header["vars"][v])
        data = read_aether_one_file(file, VarList)
        iVar_ = 3
    if (j == 0):
        [nLons, nLats, nAlts] = data[0].shape
        Alts = data[2][0][0]/1000.0;
        Lons = data[0][:,0,0]*rtod;
        Lats = data[1][0,:,0]*rtod;
        if (cut == 'alt'):
            xPos = Lons
            yPos = Lats
            if (len(Alts) > 1):
                if (args["alt"] < 50):
                    iAlt = args["alt"]
                else:
                    if (args["alt"] > Alts[nAlts-3]):
                        iAlt = nAlts-3
                    else:
                        iAlt = 2
                        while (Alts[iAlt] < args["alt"]):
                            iAlt=iAlt+1
            else:
                iAlt = 0
            Alt = Alts[iAlt]
            
        if (cut == 'lat'):
            xPos = Lons
            yPos = Alts
            if (args["lat"] < Lats[1]):
                iLat = int(nLats/2)
            else:
                if (args["lat"] > Lats[nLats-2]):
                    iLat = int(nLats/2)
                else:
                    iLat = 2
                    while (Lats[iLat] < args["lat"]):
                        iLat=iLat+1
            Lat = Lats[iLat]
            
        if (cut == 'lon'):
            xPos = Lats
            yPos = Alts
            if (args["lon"] < Lons[1]):
                iLon = int(nLons/2)
            else:
                if (args["lon"] > Lons[nLons-2]):
                    iLon = int(nLons/2)
                else:
                    iLon = 2
                    while (Lons[iLon] < args["lon"]):
                        iLon=iLon+1
            Lon = Lons[iLon]
                        
    AllTimes.append(data["time"])
    
    if (args["tec"]):
        iAlt = 2
        tec = np.zeros((nLons, nLats))
        for Alt in Alts:
            if (iAlt > 0 and iAlt < nAlts-3):
                tec = tec + data[iVar_][:,:,iAlt] * (Alts[iAlt+1]-Alts[iAlt-1])/2 * 1000.0
            iAlt=iAlt+1
        AllData2D.append(tec/1e16)
    else:
        if (cut == 'alt'):
            AllData2D.append(data[iVar_][:,:,iAlt])
        if (cut == 'lat'):
            AllData2D.append(data[iVar_][:,iLat,:])
        if (cut == 'lon'):
            AllData2D.append(data[iVar_][iLon,:,:])
        if (args["winds"]):
            if (cut == 'alt'):
                AllWindsX.append(data[iUx_][:,:,iAlt])
                AllWindsY.append(data[iUy_][:,:,iAlt])
            if (cut == 'lat'):
                AllWindsX.append(data[iUx_][:,iLat,:])
                AllWindsY.append(data[iUy_][:,iLat,:])
            if (cut == 'lon'):
                AllWindsX.append(data[iUx_][iLon,:,:])
                AllWindsY.append(data[iUy_][iLon,:,:])
    j=j+1
    

AllData2D = np.array(AllData2D)
if (args['IsLog']):
    AllData2D = np.log10(AllData2D)
if (args["winds"]):
    AllWindsX = np.array(AllWindsX)
    AllWindsY = np.array(AllWindsY)

Negative = 0

maxi  = np.max(AllData2D)*1.01
mini  = np.min(AllData2D)*0.99

if (mini < 0):
    Negative = 1

if (Negative):
    maxi = np.max(abs(AllData2D))*1.05
    mini = -maxi

if (cut == 'alt'):
    maskNorth = ((yPos>45) & (yPos<90.0))
    maskSouth = ((yPos<-45) & (yPos>-90.0))
    DoPlotNorth = np.max(maskNorth)
    DoPlotSouth = np.max(maskSouth)
    if (DoPlotNorth):
        maxiN = np.max(abs(AllData2D[:,:,maskNorth]))*1.05
        if (Negative):
            miniN = -maxiN
        else:
            miniN = np.min(AllData2D[:,:,maskNorth])*0.95
    if (DoPlotSouth):
        maxiS = np.max(abs(AllData2D[:,:,maskSouth]))*1.05
        if (Negative):
            miniS = -maxiS
        else:
            miniS = np.min(AllData2D[:,:,maskSouth])*0.95
    
dr = (maxi-mini)/31
levels = np.arange(mini, maxi, dr)

i = 0

# Define plot range:
minX = (xPos[ 1] + xPos[ 2])/2
maxX = (xPos[-2] + xPos[-3])/2
minY = (yPos[ 1] + yPos[ 2])/2
maxY = (yPos[-2] + yPos[-3])/2

file = "var%2.2d_" % args["var"]
file = file+cut

for time in AllTimes:

    ut = time.hour + time.minute/60.0 + time.second/3600.0
    shift = ut * 15.0

    fig = plt.figure(constrained_layout=False,
                     tight_layout=True, figsize=(10, 8.5))

    gs1 = GridSpec(nrows=2, ncols=2, wspace=0.0, hspace=0)
    gs = GridSpec(nrows=2, ncols=2, wspace=0.0, left=0.0, right=0.9)

    norm = cm.colors.Normalize(vmax=mini, vmin=maxi)
    if (mini >= 0):
        cmap = cm.plasma
    else:
        cmap = cm.bwr

    d2d = np.transpose(AllData2D[i])
    if (args["winds"]):
        Ux2d = np.transpose(AllWindsX[i])
        Uy2d = np.transpose(AllWindsY[i])

    sTime = time.strftime('%y%m%d_%H%M%S')
    outfile = file+'_'+sTime+'.png'

    ax = fig.add_subplot(gs1[1, :2])

    cax = ax.pcolor(xPos, yPos, d2d, vmin=mini, vmax=maxi, cmap=cmap)

    if (args["winds"]):
        ax.quiver(xPos,yPos,Ux2d,Uy2d)
    ax.set_ylim([minY,maxY])
    ax.set_xlim([minX,maxX])

    if (cut == 'alt'):
        ax.set_ylabel('Latitude (deg)')
        ax.set_xlabel('Longitude (deg)')
        title = time.strftime('%b %d, %Y %H:%M:%S')+'; Alt : '+"%.2f" % Alt + ' km'
        ax.set_aspect(1.0)

    if (cut == 'lat'):
        ax.set_xlabel('Longitude (deg)')
        ax.set_ylabel('Altitude (km)')
        title = time.strftime('%b %d, %Y %H:%M:%S')+'; Lat : '+"%.2f" % Lat + ' km'

    if (cut == 'lon'):
        ax.set_xlabel('Latitude (deg)')
        ax.set_ylabel('Altitude (km)')
        title = time.strftime('%b %d, %Y %H:%M:%S')+'; Lon : '+"%.2f" % Lon + ' km'

    ax.set_title(title)
    cbar = fig.colorbar(cax, ax=ax, shrink = 0.75, pad=0.02)
    cbar.set_label(Var,rotation=90)

    if (cut == 'alt'):
        
        if (DoPlotNorth):
            # Top Left Graph Northern Hemisphere
            ax2 = fig.add_subplot(gs[0, 0],projection='polar')
            r, theta = np.meshgrid(90.0-yPos[maskNorth], \
                                   (xPos+shift-90.0)*3.14159/180.0)
            cax2 = ax2.pcolor(theta, r, AllData2D[i][:,maskNorth], \
                              vmin=miniN, vmax=maxiN, cmap=cmap)
            xlabels = ['', '12', '18', '00']
            ylabels = ['80', '70', '60', '50']
            ax2.set_xticklabels(xlabels)
            ax2.set_yticklabels(ylabels)
            cbar2 = fig.colorbar(cax2, ax=ax2, shrink = 0.5, pad=0.01)
            ax2.grid(linestyle=':', color='black')
            pi = 3.14159
            ax2.set_xticks(np.arange(0,2*pi,pi/2))
            ax2.set_yticks(np.arange(10,50,10))
            
        if (DoPlotSouth):
            # Top Right Graph Southern Hemisphere
            r, theta = np.meshgrid(90.0+yPos[maskSouth], \
                                   (xPos+shift-90.0)*3.14159/180.0)
            ax3 = fig.add_subplot(gs[0, 1],projection='polar')
            cax3 = ax3.pcolor(theta, r, AllData2D[i][:,maskSouth], \
                              vmin=miniS, vmax=maxiS, cmap=cmap)
            xlabels = ['', '12', '18', '00']
            ylabels = ['80', '70', '60', '50']
            ax3.set_xticklabels(xlabels)
            ax3.set_yticklabels(ylabels)
            cbar3 = fig.colorbar(cax3, ax=ax3, shrink = 0.5, pad=0.01)
            ax3.grid(linestyle=':', color='black')
            pi = 3.14159
            ax3.set_xticks(np.arange(0,2*pi,pi/2))
            ax3.set_yticks(np.arange(10,50,10))

    print("Writing file : "+outfile)
    fig.savefig(outfile)
    plt.close()

    i=i+1


