#!/usr/bin/env python

from glob import glob
from datetime import datetime
from datetime import timedelta
from struct import unpack
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from pylab import cm
from gitm_routines import *
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

    for arg in argv:

        IsFound = 0

        if (not IsFound):

            m = re.match(r'-var=(.*)',arg)
            if m:
                var = int(m.group(1))
                IsFound = 1

            m = re.match(r'-alt=(.*)',arg)
            if m:
                alt = int(m.group(1))
                IsFound = 1

            m = re.match(r'-alog',arg)
            if m:
                IsLog = 1
                IsFound = 1

            if IsFound==0 and not(arg==argv[0]):
                filelist.append(arg)

    args = {'filelist':filelist,
            'var':var,
            'alt':alt,
            'IsLog':IsLog}

    return args

#-----------------------------------------------------------------------------
# 
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Main Code!
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

header = read_gitm_header()

args = get_args(sys.argv)
filelist = args["filelist"]

nFiles = len(filelist)

vars = [0,1,2]
vars.append(args["var"])

Var = header["vars"][args["var"]].replace(" ","")

AllData = []
AllData2D = []
AllAlts = []
AllTimes = []

j = 0
for file in filelist:

    data = read_gitm_one_file(file, vars)
    if (j == 0):
        [nLons, nLats, nAlts] = data[0].shape
        Alts = data[2][0][0]/1000.0;
        Lons = data[0][:,0,0]*rtod;
        Lats = data[1][0,:,0]*rtod;
        if (args["alt"] < 50):
            iAlt = args["alt"]
        else:
            if (args["alt"] > Alts[nAlts-1]):
                iAlt = nAlts-3
            else:
                iAlt = 2
                while (Alts[iAlt] < args["alt"]):
                    iAlt=iAlt+1
        Alt = Alts[iAlt]
    AllTimes.append(data["time"])
    AllData.append(data[args["var"]][nLons/2][nLats/2])
    AllData2D.append(data[args["var"]][:,:,iAlt])
    j=j+1
    

AllData2D = np.array(AllData2D)

maxi  = np.max(AllData2D)*1.05
mini  = np.min(AllData2D)*0.95
dr = (maxi-mini)/31
levels = np.arange(mini, maxi, dr)

i = 0
for time in AllTimes:

    fig = plt.figure()
    ax = fig.add_subplot(111)

    norm = cm.colors.Normalize(vmax=np.max(AllData2D), vmin=np.min(AllData2D))
    cmap = cm.rainbow

    d2d = AllData2D[i]

    outfile = time.strftime('test_%y%m%d_%H%M%S.png')

    cax = ax.pcolor(Lons, Lats, np.transpose(d2d), vmin=mini, vmax=maxi)
    ax.set_ylim([-90,90])
    ax.set_xlim([0,360])
    title = time.strftime('%b %d, %Y %H:%M:%S')+'; Alt : '+"%.2f" % Alt + ' km'
    ax.set_title(title)
    ax.set_xlabel('Latitude (deg)')
    ax.set_ylabel('Longitude (deg)')

    cbar = fig.colorbar(cax)
    cbar.set_label(Var,rotation=90)
    fig.savefig(outfile)

    i=i+1


