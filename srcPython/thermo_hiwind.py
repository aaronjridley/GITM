#!/usr/bin/env python

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

from thermo_get_points import *
from read_hiwinds import *

def parse_args_hiwinds():

    parser = argparse.ArgumentParser(description =
                                     'Compare GITM to HI-WIND files')
    parser.add_argument('files',
                        nargs = '+',
                        help = 'List of files to read and plot')

    args = parser.parse_args()

    return args


def local_time(times, lons):

    uts = [t.hour + t.minute/60.0 + t.second/3600.0 for t in times]
    localtime = (lons/15.0 + uts) % 24.0
    return localtime

#------------------------------------------------------------------------------
# main code
#------------------------------------------------------------------------------

args = parse_args_hiwinds()
filelist = args.files

file = filelist[0]

print('Reading file: ',file)
data = read_hiwind_file(file)

lats = data['lats']
lons = data['lons']
alts = lats * 0.0 + 250.0

mtime = data['merid_time']

localtimes = local_time(mtime, lons)
d = np.abs(localtimes - 6)
tDawn = mtime[np.argmin(d)]

d = np.abs(localtimes - 12)
tNoon = mtime[np.argmin(d)]

d = np.abs(localtimes - 18)
tDusk = mtime[np.argmin(d)]

vars = [17] # n/s wind in GITM
merid_gitm, mtime_gitm = thermo_get_points(lats, lons, alts, mtime, vars)
merid_gitm.flatten()
mtime_gitm.flatten()

ztime = data['zonal_time']
vars = [16] # e/w wind in GITM
zonal_gitm, ztime_gitm = thermo_get_points(lats, lons, alts, mtime, vars)
zonal_gitm.flatten()
ztime_gitm.flatten()

mint = np.min([np.min(mtime), np.min(ztime)])
maxt = np.max([np.max(mtime), np.max(ztime)])

sTime = mint.strftime('%b %d, %Y %H:%M UT') + ' - ' + \
    maxt.strftime('%b %d, %Y %H:%M UT (Hours)')
sDate = mtime[int(len(mtime)/2)].strftime('%b %d, %Y')
        
fig = plt.figure(figsize = (10,10))

axm = fig.add_subplot(211)
axm.plot(mtime, data['merid_wind'], 'r', label = 'HiWind')
axm.plot(mtime_gitm, merid_gitm, 'b', label = 'GITM')
axm.set_xlim(mint, maxt)
axm.set_ylabel('Meridional Wind (m/s)')
axm.axhline(y=0.0, color='r', linestyle=':')
axm.axvline(x = tDawn, color = 'g', linestyle = 'dashed')
axm.axvline(x = tNoon, color = 'k', linestyle = 'solid')
axm.text(tNoon, np.min(merid_gitm), ' Noon')
axm.axvline(x = tDusk, color = 'g', linestyle = 'dashed')
axm.legend()
axm.set_title('GITM Comparison to HIWIND on ' + sDate)

axz = fig.add_subplot(212)
axz.plot(ztime, data['zonal_wind'], 'r', label = 'HiWind')
axz.plot(ztime_gitm, zonal_gitm, 'b', label = 'GITM')
axz.set_xlim(mint, maxt)
axz.set_ylabel('Zonal Wind (m/s)')
axz.set_xlabel(sTime)
axz.axhline(y=0.0, color='r', linestyle=':')
axz.axvline(x = tDawn, color = 'g', linestyle = 'dashed')
axz.axvline(x = tNoon, color = 'k', linestyle = 'solid')
axz.text(tNoon, np.min(zonal_gitm), ' Noon')
axz.axvline(x = tDusk, color = 'g', linestyle = 'dashed')
axz.legend()

plotfile = mtime[int(len(mtime)/2)].strftime('gitm_hiwind_%Y%m%d.png')
print('writing : ',plotfile)    
fig.savefig(plotfile)
plt.close()

