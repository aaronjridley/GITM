#!/usr/bin/env python

from datetime import datetime
from datetime import timedelta
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from gitm_routines import *

def thermo_get_points(lats_in, lons_in, alts_in, times_in, vars_in):

    rtod =  180.0 / np.pi
    dtor = 1.0 / rtod

    print("Reading gitm headers")
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
        
    print("Looping through times...")

    times_out = []
    values_out = []

    iBefore0 = -1
    iAfter0 = -1
    # here we just assume that we will start with the first set of files.
    # this shouldn't matter...
    iBefore = 0
    iAfter = 1
    nVars = len(vars_in)
    AfterVals = np.zeros(nVars)
    BeforeVals = np.zeros(nVars)

    for i, time in enumerate(times_in):

        while (time > headers["time"][iAfter]):
            iAfter = iAfter+1
            if (iAfter >= nGitmFiles-1):
                break
        
        if (iAfter == nGitmFiles):
            break

        iBefore = iAfter-1

        if (iBefore != iBefore0):
            file = headers["filename"][iBefore]
            BeforeData = read_gitm_one_file(file, vars_in)
        if (iAfter != iAfter0):
            file = headers["filename"][iAfter]
            AfterData = read_gitm_one_file(file, vars_in)

        if (time >= headers["time"][iBefore]):

            times_out.append(time)
            
            dt = (headers["time"][iAfter] - \
                  headers["time"][iBefore]).total_seconds()
            xt = (time - headers["time"][iBefore]).total_seconds() / dt

            lon = lons_in[i]
            lat = lats_in[i]
            alt = alts_in[i]
        
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

            for i, v in enumerate(vars_in):
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
            
            vals = []
            for i, v in enumerate(vars_in):
                vals.append( (1-xt) * BeforeVals[i] + xt * AfterVals[i] )

            values_out.append(vals)

        iBefore0 = iBefore
        iAfter0 = iAfter

    return np.array(values_out), np.array(times_out)
