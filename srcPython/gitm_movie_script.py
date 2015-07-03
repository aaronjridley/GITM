#!/usr/bin/env python
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#----------------------------------------------------------------------------
# $Id: gitm_movie_script.py,v 1.2 2013/10/12 04:01:00 kopmanis Exp $
#
# GITM_Movie_Script
#
# Author: Angeline G. Burrell, UMich, July 2013
#
# Comments: A script to make a movie of GITM 2D or 3D parameters
#-----------------------------------------------------------------------------

import matplotlib as mpl
mpl.use('Agg') # This needs to be fixed
import string
import subprocess # To issue commands to the OS
import spacepy
import gitm
import gitm_3D_global_plots as gplots
import gitm_time as gt
import gitm_plot_rout as gpr
import numpy as np

#------------------------------------
# Initialize the data lists

linc = 6
zflag = True
zkeys = list()
zalts = list()
zmaxs = list()
zmins = list()
zunits = list()
zindexes = list()
earths = list()
plots = list()
nlimits = list()
slimits = list()

#--------------------------------------
# Get the desired input values

gitmlist = raw_input('Ordered list of GITM binary files: ')

while zflag:
    zplot=raw_input('GITM plot type (rectangular, polar, nspolar, snapshot): ')
    plots.append(zplot)
    zkey = raw_input('GITM key to plot on z axis (eg Temperature): ')
    zkeys.append(zkey)
    zalt = float(raw_input('Altitude to plot z value at (eg 250): '))
    zalts.append(zalt)
    zunit = raw_input('Units of altitude (km or m): ')
    zunits.append(zunit)
    earth = bool(raw_input('Use map of Earth? (empty for False): '))
    earths.append(earth)

    if zplot == "nspolar":
        nlat = float(raw_input('Polar latitude limit (degrees): '))
        nlimits.append(nlat)
        slat = float(raw_input('Equatorial latitude limit (degrees): '))
        slimits.append(slat)
    elif zplot == "snapshot":
        nlat = float(raw_input('Polar latitude limit (degrees): '))
        nlimits.append(nlat)
        slimits.append(0.0)
    else:
        nlat = float(raw_input('Northern latitude limit (degrees): '))
        nlimits.append(nlat)
        slat = float(raw_input('Southern latitude limit (degrees): '))
        slimits.append(slat)

    zflag = bool(raw_input('Load another z axis key? (empty for False): '))

#------------------------------------------
# Load the GITM binaries

gitmbins = gt.load_multiple_gitm_bin(gitmlist)

#---------------------------------------------
# Determine the limits of each z parameter.  To do this we need to open
# one of the gitm files

for z,zkey in enumerate(zkeys):
    # Find the index for the specified altitude.  This will be the same for
    # all files, longitudes, and latititudes.

    aindex = gpr.find_alt_index(gitmbins[0], 1, 1, zalts[z], zunits[z])
    zindexes.append(aindex)

    # Find the indexes for the specified latitude limits.  The latitude will
    # be the same for all files and longitudes.  The northern limit will be
    # the same for all plot types.

    (nlon,nlat) = gpr.find_lon_lat_index(gitmbins[0], -1.0 , nlimits[z],
                                         "degrees")

    if zplot == "nspolar" or zplot == "snapshot":
        slimit = -1.0 * nlimits[z]
    else:
        slimit = slimits[z]

    (slon,slat) = gpr.find_lon_lat_index(gitmbins[0], 361.0, slimit, "degrees")
    
    # Test to ensure the validity of the indexes

    if nlat == slat:
        slat += 1
    elif nlat < slat:
        temp = nlat
        nlat = slat
        slat = temp

    # Find the z parameter limits within the desired lat/lon/alt range

    (zmin, zmax) = gpr.find_data_limits_irange(gitmbins, zkey, slat, nlat, nlon,
                                               slon, aindex, aindex+1)

    zmins.append(zmin)
    zmaxs.append(zmax)

#------------------------------------------
# For each key, create the desired images

plotfiles = list()
nplots = len(gitmbins)
nmag = np.floor(np.log10(nplots)) + 1
i = 1

while gitmbins:
    # Get the daet from the base filename
    sfile = gitmbins[0].attrs['file'].split("/")
    sbase = sfile[-1].split(".")
    sdate = sbase[0].split("_")
    fdate = "%s/%s/%s-%s:%s:%s" % (sdate[1][1:3], sdate[1][3:5], sdate[1][5:8],
                                   sdate[2][0:2], sdate[2][2:4], sdate[2][4:7])

    # Construct the plot number

    sbase = str(int(10**nmag + i))

    # Cycle through the z keys

    for z,zkey in enumerate(zkeys):
        title = "%s %s" % (gitmbins[0][zkey].attrs['name'], fdate)
        fkey = string.join(zkey.split(" "), "")

        if plots[z] == "polar":
            # Polar plot
            if i == 1:
                plotfiles.append("Pole_%s_%d%s" % (fkey, zalts[z], zunits[z]))

            figname = "%s_%s.png" % (plotfiles[z], sbase[1:len(sbase)])
            gplots.plot_single_3D_image("polar", zkey, gitmbins[0], title,
                                        figname, False, zindexes[z], nlimits[z],
                                        slimits[z], linc, earths[z], 90,
                                        zmaxs[z], zmins[z])
        elif plots[z] == "rectangular":
            # Global plot
            if i == 1:
                plotfiles.append("Globe_%s_%d%s" % (fkey, zalts[z], zunits[z]))
            figname = "%s_%s.png" % (plotfiles[z], sbase[1:len(sbase)])
            gplots.plot_single_3D_image("rectangular", zkey, gitmbins[0], title,
                                        figname, False, zindexes[z], nlimits[z],
                                        slimits[z], linc, earths[z], 90,
                                        zmaxs[z], zmins[z])
        elif plots[z] == "nspolar":
            # Global plot with polar circles
            if i == 1:
                plotfiles.append("NSG_%s_%d%s" % (fkey, zalts[z], zunits[z]))
            figname = "%s_%s.png" % (plotfiles[z], sbase[1:len(sbase)])
            gplots.plot_single_nsglobal_3D_image(zkey, gitmbins[0], title,
                                                 figname, False, zindexes[z],
                                                 nlimits[z], slimits[z], 3, 90,
                                                 earths[z], zmaxs[z], zmins[z])
        else:
            # Global snapshot
            if i == 1:
                plotfiles.append("Gsnap_%s_%d%s" % (fkey, zalts[z], zunits[z]))
            figname = "%s_%s.png" % (plotfiles[z], sbase[1:len(sbase)])
            gplots.plot_global_3D_snapshot(zkey, gitmbins[0], title, figname,
                                           False, zindexes[z], 90, earths[z],
                                           zmaxs[z], zmins[z])
    del gitmbins[0]
    i += 1

#---------------------------------------------------
# Create movies from each of the series of figures

for plotbase in plotfiles:
    # Construct the ffmpeg command
    video = ('ffmpeg', '-f', 'image2', '-i', "%s_%%%dd.png" % (plotbase, nmag),
             "%s.mpg" % (plotbase))

    # Execute the movie command using ffmpeg
    try:
        subprocess.call(video)

        print "Successfully created video: %s.png" % (plotbase) 
    except:
        com = string.join(video, " ")
        print "Unable to execute video command, is ffmpeg installed? ", com

# End
