#!/usr/bin/env python
#-----------------------------------------------------------------------------
# $Id: gitm_3D_global_plots.py,v 1.16 2014/02/17 17:30:33 agburr Exp $
# gitm_3D_global_plots
#
# Author: Angeline G. Burrell, UMichigan, Jan 2013
#
# Comments: Routine to make 3D color plots of thermospheric and ionospheric
#           quantities that are either at a single altitude or are constant
#           with altitude as functions of geographic latitude and longitude.
#           Map contours are also availabel if you have the mpl_toolkit
#           baseline installed (available on fink, macports, or online at
#           http://matplotlib.org/basemap).
#
# AGB: 10/9/13: Redid plots to use subroutines that allow any class of data
#              (not just GitmBin) be used.  Added option for scatter plots in
#               addition to contour plots.  Added options for solar terminator
#               and geomagnetic equator.  Added option to define the high-lat
#               limit in the snapshot plots.
#
# Includes: gitm_single_3D_image          - plots a single rectangular or polar
#                                           alt slice as a contour or scatter
#           gitm_single_nsglobal_3D_image - plots northern and southern polar
#                                           altitude slice as a contour or
#                                           scatter
#           gitm_global_3D_snapshot       - plots northern and southern polar
#                                           and a midlatitude rectangular
#                                           altitude slice as a contour or
#                                           scatter
#           gitm_mult_3D_slices           - plot multiple polar or rectangular
#                                           altitude slices  as a contour or
#                                           scatter
#----------------------------------------------------------------------------

'''
Plot GITM data in lat/lon contours or color-coded scatter plots to facilitate
model evaluation
'''

# Import modules
import math
import numpy as np
from spacepy.pybats import gitm
import gitm # Temporary until the new GITM is incorporated into spacepy
import matplotlib.pyplot as plt
import gitm_plot_rout as gpr
import plot_3D_global as p3g

def gitm_single_3D_image(plot_type, zkey, gdata, title=None, figname=None,
                         draw=True, aindex=-1, nlat=90, slat=-90, linc=6,
                         earth=False, tlon=90, zmax=None, zmin=None,
                         zcolor=None, data_type="contour", faspect=True,
                         meq=False, terminator=False, m=None, *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range.
    Input: plot_type  = key to determine plot type (rectangular, polar)
           zkey       = key for z variable (ie 'Vertical TEC')
           gData      = gitm bin structure
           title      = plot title
           figname    = file name to save figure as (default is none)
           draw       = draw to screen? (default is True)
           aindex     = altitude index (default -1 if it is a 2D parameter)
           nlat       = northern latitude limit (degrees North, default 90)
           slat       = southern latitude limit (degrees North, defalut -90)
           linc       = number of latitude tick incriments (default 6)
           earth      = include continent outlines for Earth (default False)
           tlon       = longitude at the top of a polar dial (degrees east,
                        default 90)
           zmax       = Maximum z range (default None)
           zmin       = Minimum z range (default None)
           zcolor     = Color map for the z variable.  If none, will be chosen
                        based on the z range (default=None)
           data_type  = scatter or contour (default=scatter)
           faspect    = Fix the aspect of Earth if using outlines (default=True)
           meq        = Include the geomagnetic equator?  (default=False)
           terminator = Include the solar terminator? (default=False)
           m          = Handle for earth map (default=None)

    Output: f = figure handle
            m = map handle (or None)
    '''
    # Set the altitude and latitude limits.  For polar plots, latitudes
    # Above and below 90 degrees will cause the routine to fail
    ialt = 0
    if aindex > 0:
        ialt = aindex

    if plot_type.find("polar") >= 0:
        (i, imin) = gpr.find_lon_lat_index(gdata, 0.0, -90.0, "degrees")
        (i, imax) = gpr.find_lon_lat_index(gdata, 0.0, 90.0, "degrees")
        imax += 1
    else:
        imin = 0
        imax = gdata.attrs['nLat']

    # If including the solar terminator, extract the UT date
    tdt = None
    if terminator:
        tdt = gdata['time']

    # Format title
    spec_title = "{:} UT".format(gdata['time'])
    if aindex >= 0:
        spec_title = "{:s} slice at {:.2f} km".format(spec_title, 1.0e-3 *
                                                      gdata['Altitude'][0,0,
                                                                        ialt])

    if title:
        title = "{:s}\n{:s}".format(spec_title, title)
    else:
        title = spec_title

    # Output the figure
    fm = p3g.plot_single_3D_image(plot_type,
                                  np.array(gdata['dLat'][:,imin:imax,ialt]),
                                  np.array(gdata['dLon'][:,imin:imax,ialt]),
                                  np.array(gdata[zkey][:,imin:imax,ialt]),
                                  gdata[zkey].attrs['name'],
                                  gdata[zkey].attrs['scale'],
                                  gdata[zkey].attrs['units'], zmax=zmax,
                                  zmin=zmin, zcolor=zcolor, title=title,
                                  figname=figname, draw=draw, nlat=nlat,
                                  slat=slat, linc=linc, tlon=tlon,
                                  data_type=data_type, meq=meq, earth=earth,
                                  m=m, faspect=faspect, term_datetime=tdt)
    return fm
# End gitm_single_3D_image

def gitm_single_nsglobal_3D_image(zkey, gdata, title=None, figname=None,
                                  draw=True, aindex=-1, plat=90, elat=0, linc=3,
                                  tlon=90, earth=False, zmax=None, zmin=None,
                                  zcolor=None, data_type="contour",
                                  terminator=False, mn=None, ms=None,
                                  *args, **kwargs):
    '''
    Creates a figure with two polar map projections for the northern and 
    southern ends of a specified latitude range.
    Input: zkey       = key for z variable (ie 'Vertical TEC')
           gdata      = gitm bin structure
           title      = plot title
           figname    = file name to save figure as (default is none)
           draw       = draw to screen? (default is True)
           aindex     = altitude index (default -1 if it is a 2D parameter)
           plat       = polar latitude limit (degrees North, default +/-90)
           elat       = equatorial latitude limit (degrees North, defalut 0)
           linc       = number of latitude tick incriments (default 6)
           tlon       = longitude to place on the polar dial top (degrees east,
                        default 90)
           earth      = include Earth continent outlines (default False)
           zmax       = maximum z range (default None)
           zmin       = mininimum z range (default None)
           zcolor     = Color scale for plotting the z data.  If not specified,
                        this will be determined by the z range (default=None)
           data_type  = Type of plot to make scatter/contour (default=scatter)
           terminator = Include the solar terminator by shading the night
                        time regions? (default=False)
           mn         = Northern latitude map handle (default=None)
           ms         = Southern latitude map handle (default=None)

    Output: f  = figure handle
            mn = Northern latitude map handle
            ms = Southern latitude map handle
    '''
    # Set the altitude and latitude limits.  For polar plots, latitudes
    # Above and below 90 degrees will cause the routine to fail
    ialt = 0
    if aindex > 0:
        ialt = aindex

    (i, imin) = gpr.find_lon_lat_index(gdata, 0.0, -90.0, "degrees")
    (i, imax) = gpr.find_lon_lat_index(gdata, 0.0, 90.0, "degrees")
    imax += 1

    # If including the solar terminator, extract the UT date
    tdt = None
    if terminator:
        tdt = gdata['time']

    # Format title
    spec_title = "{:} UT".format(gdata['time'])
    if aindex >= 0:
        spec_title = "{:s} slice at {:.2f} km".format(spec_title, 1.0e-3 *
                                                      gdata['Altitude'][0,0,
                                                                        ialt])
    if title:
        title = "{:s}\n{:s}".format(spec_title, title)
    else:
        title = spec_title

    # Output figure
    fm = p3g.plot_single_nsglobal_3D_image(np.array(gdata['dLat'][:,imin:imax,ialt]), np.array(gdata['dLon'][:,imin:imax,ialt]), np.array(gdata[zkey][:,imin:imax,ialt]),
                                           gdata[zkey].attrs['name'],
                                           gdata[zkey].attrs['scale'],
                                           gdata[zkey].attrs['units'],
                                           zmax=zmax, zmin=zmin, zcolor=zcolor,
                                           title=title, figname=figname,
                                           draw=draw, plat=plat, elat=elat,
                                           linc=linc, tlon=tlon, earth=earth,
                                           mn=mn, ms=ms, data_type=data_type,
                                           term_datetime=tdt)

    return fm
# End gitm_single_nsglobal_3D_image

def gitm_global_3D_snapshot(zkey, gdata, title=None, figname=None, draw=True,
                            aindex=-1, tlon=90, polar_blat=45, rect_blat=45,
                            earth=False, zmax=None, zmin=None,
                            zcolor=None, meq=False, data_type="contour",
                            terminator=False, ml=None, mn=None, ms=None,
                            *args, **kwargs):
    '''
    Creates a map projection plot for the entire globe, seperating the polar
    and central latitude regions.
    Input: zkey       = key for z variable (ie 'Vertical TEC')
           gData      = gitm bin structure
           title      = plot title
           figname    = file name to save figure as (default is none)
           draw       = output a screen image? (default is True)
           aindex     = altitude index (default -1 if it is a 2D parameter)
           tlon       = longitude at the top of the polar dial (degrees East,
                        default 90)
           polar_blat = co-latitude of the lower boundary of the polar dials
                        (default 45)
           rect_blat  = Upper bounding co-latitude of the rectangular map
                        (default 45)
           earth      = include Earth continent outlines (default False)
           zmax       = maximum z limit (default None)
           zmin       = minimum z limit (default None)
           zcolor     = Color scale for z variable.  If not specified, will be
                        determined by the z range (default=None)
           meq        = Add a line for the geomagnetic equator? (default=False)
           data_type  = Type of plot to make scatter/contour (default=contour)
           terminator = Include the solar terminator by shading the night
                        time regions? (default=False)
           ml         = Low latitude map handle (default=None)
           mn         = Northern latitude map handle (default=None)
           ms         = Southern latitude map handle (default=None)

    Output: f  = Figure handle
            ml = Low latitude map handle
            mn = Northern latitude map handle
            ms = Southern latitude map handle
    '''
    # Set the altitude and latitude limits.  For polar plots, latitudes
    # Above and below 90 degrees will cause the routine to fail
    ialt = 0
    if aindex > 0:
        ialt = aindex

    (i, imin) = gpr.find_lon_lat_index(gdata, 0.0, -90.0, "degrees")
    (i, imax) = gpr.find_lon_lat_index(gdata, 0.0, 90.0, "degrees")
    imax += 1

    # If including the solar terminator, extract the UT date
    tdt = None
    if terminator:
        tdt = gdata['time']

    # Format title
    spec_title = "{:} UT".format(gdata['time'])
    if aindex >= 0:
        spec_title = "{:s} slice at {:.2f} km".format(spec_title, 1.0e-3 *
                                                      gdata['Altitude'][0,0,
                                                                        ialt])

    if title:
        title = "{:s}\n{:s}".format(spec_title, title)
    else:
        title = spec_title

    # Output figure
    fm = p3g.plot_global_3D_snapshot(np.array(gdata['dLat'][:,imin:imax,ialt]),
                                     np.array(gdata['dLon'][:,imin:imax,ialt]),
                                     np.array(gdata[zkey][:,imin:imax,ialt]),
                                     gdata[zkey].attrs['name'],
                                     gdata[zkey].attrs['scale'],
                                     gdata[zkey].attrs['units'], zmax=zmax,
                                     zmin=zmin, zcolor=zcolor, title=title,
                                     figname=figname, draw=draw, tlon=tlon,
                                     polar_blat=polar_blat, rect_blat=rect_blat,
                                     meq=meq, earth=earth, ml=ml, mn=mn, ms=ms,
                                     data_type=data_type, term_datetime=tdt)
    return fm
# End gitm_global_3D_snapshot

def gitm_mult_3D_slices(plot_type, zkey, gdata, aindex, title=None,
                        figname=None, draw=True, nlat=90, slat=-90, linc=6,
                        earth=False, tlon=90, zmax=None, zmin=None, 
                        zcolor=None, data_type="contour", meq=False,
                        faspect=True, terminator=False, *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range.
    Input: plot_type  = key to determine plot type (rectangular, polar)
           zkey       = key for z variable (ie 'Vertical TEC')
           gdata      = gitm bin structure
           aindex     = list of altitude indices
           title      = plot title
           figname    = file name to save figure as (default is none)
           draw       = draw to screen? (default is True)
           nlat       = northern latitude limit (degrees North, default 90)
           slat       = southern latitude limit (degrees North, defalut -90)
           linc       = number of latitude tick incriments (default 6)
           earth      = include Earth continent outlines (default False)
           tlon       = longitude on the top, for polar plots (degrees East,
                        default 90)
           zmax       = maximum z range (default None)
           zmin       = minimum z range (default None)
           zcolor     = Color spectrum for z data.  If not specified, will be
                        determined by the z range. (default=None)
           data_type  = Contour or scatter plot? (default=scatter)
           meq        = Add a line for the geomagnetic equator? (default=False)
           earth      = include Earth continent outlines (default False)
           faspect    = Fix aspect ratio if using continents (default=True)
           terminator = Include the solar terminator by shading the night
                        time regions? (default=False)
    '''
    # Set the latitude limits.  For polar plots, latitudes above and below
    # 90 degrees will cause the routine to fail
    if plot_type.find("polar") >= 0:
        (i, imin) = gpr.find_lon_lat_index(gdata, 0.0, -90.0, "degrees")
        (i, imax) = gpr.find_lon_lat_index(gdata, 0.0, 90.0, "degrees")
        imax += 1
    else:
        imin = 0
        imax = gdata.attrs['nLat']

    # If including the solar terminator, extract the UT date
    tdt = None
    if terminator:
        tdt = gdata['time']

    # Format title
    spec_title = "{:} UT slice".format(gdata['time'])

    if title:
        title = "{:s}\n{:s}".format(spec_title, title)
    else:
        title = spec_title

    # Output the figure
    f = p3g.plot_mult_3D_slices(plot_type, 2, aindex,
                                np.array(gdata['dLat'][:,imin:imax,:]),
                                np.array(gdata['dLon'][:,imin:imax,:]),
                                np.array(gdata[zkey][:,imin:imax,:]),
                                gdata[zkey].attrs['name'],
                                gdata[zkey].attrs['scale'],
                                gdata[zkey].attrs['units'], zmax=zmax,
                                zmin=zmin, zcolor=zcolor, title=title,
                                figname=figname, draw=draw, nlat=nlat,
                                slat=slat, linc=linc, tlon=tlon,
                                data_type=data_type, meq=meq, earth=earth,
                                faspect=faspect, term_datetime=tdt)
    return f
# End gitm_mult_3D_slices
