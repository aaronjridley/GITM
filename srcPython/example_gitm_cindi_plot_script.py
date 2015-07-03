#!/usr/bin/env python
#-----------------------------------------------------------------------------
# example_gitm_plot_script
#
# Author: Angeline G. Burrell, UMichigan, May 2013
#
# Comments: Script routine that creates a figure comparing CINDI and GITM data
# from the command line
#----------------------------------------------------------------------------

import string
from os import path
import spacepy
import gitm
import gitm_time
import gitm_plot_rout as gpr
import gitm_loc_rout as glr
import gitm_comparison_plots as gcp
import load_files as lf

#------------------------------------------------------------------------------
# Set defaults for the input parameters
errkey = False

#------------------------------------------------------------------------------
# Get the desired input values FIX
filename = raw_input('Output image file (eg Test.png): ')
gitmlist = raw_input('List of input GITM binaries: ')
gitmtype = raw_input('GITM binary type (3D, 2D, 1D): ')
gitmmag = raw_input('Name of GITM 3DMAG or 3DION file: ')
cindilist = raw_input('List of input CINDI text files: ')
match = raw_input('Satellite-Model match method (nearest, average, median): ')
ykey = raw_input('GITM key of the dependent variable: ')
title = raw_input('Plot title: ')

max_tdelt = 300.0
max_ldelt = 100000.0

if match.find("nearest") >= 0:
    max_tdelt = float(raw_input('Maximum allowable seperation (in seconds) between nearest neighbor points: '))
    max_ldelt = float(raw_input('Maximum allowable seperation (in meters) between nearest neighbor points: '))
elif len(match) > 0:
    max_alt = float(raw_input('Width of alt boxcar for ave/median (in km): '))
    max_lat = float(raw_input('Width of lat boxcar for ave/median (in deg): '))
    max_lon = float(raw_input('Width of lon boxcar for ave/median (in deg): '))
    max_ldelt = [max_lon, max_lat, max_alt]
    max_tdelt = float(raw_input('Width of UT boxcar for ave/median (in sec): '))

#------------------------------------------------------------------------------
# Test the desired input values
if gitmtype.find("3D")<0 and gitmtype.find("2D")<0 and gitmtype.find("1D")<0:
    print "ERROR: unknown GITM binary type [", gitmtype, "]"
    errkey = True

if(match.find("nearest")<0 and match.find("average")<0
   and match.find("median")<0):
    print "ERROR: unknown match method [", match, "]"
    errkey = True

if(filename.find(".png") < 1):
    print "WARNING: unknown image name [%s], no file will be output" % filename
    filename = None

lsize = path.getsize(gitmlist)
if lsize > 2.0e9:
    print "ERROR: GITM file list is too long"
    errkey = True
elif lsize == 0:
    print "ERROR: GITM file list is empty"
    errkey = True

lsize = path.getsize(cindilist)
if lsize > 2.0e9:
    print "ERROR: CINDI file list is too long"
    errkey = True
elif lsize == 0:
    print "ERROR: CINDI file list is empty"
    errkey = True

cindi_ykey = gpr.match_cindi_key(ykey)

if not cindi_ykey:
    errkey = True

#---------------------------------
# Load the CINDI data
if not errkey:
    f = open(cindilist, "r")
    cindi_list = f.readlines()
    for i,u in enumerate(cindi_list):
        cindi_list[i] = u[0:len(u)-1]

    nfiles, cdata = lf.loadCINDIorbit_ASCII(cindi_list)

    if nfiles < 1:
        print "ERROR: unable to load CINDI data from files in [", cindilist, "]"
        errkey = True

#---------------------------------
# Match the GITM and CINDI data
if not errkey:
    f = open(gitmlist, "r")
    gitm_list = f.readlines()
    for i,u in enumerate(gitm_list):
        gitm_list[i] = u[0:len(u)-1]

    magfile = None
    if len(gitmmag) > 0:
        magfile = gitmmag

    ckeys = gpr.match_cindi_key("ALL", "GITM")
    ctrack = glr.gitm_inst_loc(cdata['datetime'], cdata['GLAT'], cdata['GLONG'],
                               cdata['Alt.'], ckeys, "satellite", gitm_list,
                               gitmtype, magfile)

    if ctrack:
        # For certian types of data, the CINDI and GITM results require some
        # processing.  Do that here to the cindi data.
        if cindi_ykey.find("Fr") == 0:
            cdata[cindi_ykey] *= cdata['Ni(cm^-3)']*1.0e6 # cm^-3 -> m^-3

        if match.find("nearest") >= 0:
            obs_key = ctrack.appendobs([cdata[cindi_ykey]], [cindi_ykey], match,
                                       cdata['datetime'], cdata['GLONG'],
                                       cdata['GLAT'], cdata['Alt.'], max_tdelt,
                                       max_ldelt)
        else:
            obs_key = ctrack.appendobs([cdata[cindi_ykey]], [cindi_ykey], match,
                                       cdata['datetime'], cdata['GLONG'],
                                       cdata['GLAT'], cdata['Alt.'],
                                       max_locdelt=max_ldelt,
                                       boxcar_sec=max_tdelt)

        if not ctrack.has_key(obs_key[0]):
            errkey = True
    else:
        errkey = True

#---------------------------------
# Plot the data
if not errkey:
    tkey = None
    bkey = None
    if match.find("average") >= 0:
        tkey = "{:s}_sig".format(obs_key[0])
        bkey = "{:s}_sig".format(obs_key[0])
    elif match.find("median") >= 0:
        tkey = "{:s}_terr".format(obs_key[0])
        bkey = "{:s}_berr".format(obs_key[0])

    # Set UT limits
    tmin = ctrack['time'][0]
    tmax = ctrack['time'][-1]
    # Set percentage limits
    pmin = -100.0
    pmax = 100.0

    gcp.plot_sat_gitm_comp(cdata['datetime'], cdata[cindi_ykey], ctrack, ykey,
                           obs_key[0], tkey, bkey, title, tmin=tmin, tmax=tmax,
                           pmin=pmin, pmax=pmax, figname=filename, draw=False)
                           

# End
