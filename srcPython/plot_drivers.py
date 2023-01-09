#!/usr/bin/env python

import re
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates

from swmf_imf import *
from read_ae import *

#------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------

def plot_day_boundaries(ax, times):

    year = times[0].year
    month = times[0].month
    day = times[0].day
    nDays = int((times[-1] - times[0]).total_seconds()/86400.0) + 1

    for iDay in range(nDays):
        d = dt.datetime(year, month, day) + dt.timedelta(days=iDay)
        ax.axvline(d, color = 'k', ls = '--')

    ax.set_xlim(times[0], times[-1])
    ax.xaxis.set_major_formatter(dates.DateFormatter('%m-%d'))
        
    return

imfFile = 'imf.dat'
aeFile = 'ae_20180624.txt'

imf = read_swmf_file(imfFile)
ae = read_ae(aeFile)

fig = plt.figure(figsize = (7,10))
plt.rcParams.update({'font.size': 14})

zeros = np.array(imf["bz"]) * 0.0

nY = 4
ax = []

xStart = 0.15
xSize = 0.81
ySize = 0.19
yGap = 0.04
yEnd = 0.95

iY = 1
ax.append(fig.add_subplot(nY*100 + 10 + iY))
yStart = yEnd - ySize * iY - yGap * (iY-1)
ax[-1].set_position([xStart, yStart, xSize, ySize])

ax[-1].plot(imf["times"], imf["bz"])
ax[-1].plot(imf["times"], zeros, 'k:')
ax[-1].set_ylabel('(a) IMF Bz (nT)')
ax[-1].set_title('IMF, Solar Wind, and AE Drivers for GITM')
plot_day_boundaries(ax[-1], imf["times"])

iY = 2
ax.append(fig.add_subplot(nY*100 + 10 + iY))
yStart = yEnd - ySize * iY - yGap * (iY-1)
ax[-1].set_position([xStart, yStart, xSize, ySize])
ax[-1].plot(imf["times"], imf["by"])
ax[-1].plot(imf["times"], zeros, 'k:')
ax[-1].set_ylabel('(b) IMF By (nT)')
plot_day_boundaries(ax[-1], imf["times"])

iY = 3
ax.append(fig.add_subplot(nY*100 + 10 + iY))
yStart = yEnd - ySize * iY - yGap * (iY-1)
ax[-1].set_position([xStart, yStart, xSize, ySize])
ax[-1].plot(imf["times"], imf["vx"])
ax[-1].set_ylabel('(c) SW Vx (km/s)')
plot_day_boundaries(ax[-1], imf["times"])

iY = 4
ax.append(fig.add_subplot(nY*100 + 10 + iY))
yStart = yEnd - ySize * iY - yGap * (iY-1)
ax[-1].set_position([xStart, yStart, xSize, ySize])
ax[-1].plot(ae["time"], ae["ae"])
ax[-1].set_ylabel('(d) AE (nT)')
plot_day_boundaries(ax[-1], imf["times"])

mint = np.min(imf["times"])
maxt = np.max(imf["times"])

sTime = mint.strftime('%b %d, %Y %H:%M UT') + ' - ' + \
    maxt.strftime('%b %d, %Y %H:%M UT')
ax[-1].set_xlabel(sTime)

plotfile = imf["times"][0].strftime('drivers%Y%m%d.png')
fig.savefig(plotfile)
plt.close()

