#!/usr/bin/env python
#-----------------------------------------------------------------------------
# $Id: plot_3D_global.py,v 1.14 2014/03/27 13:56:40 agburr Exp $
# plot_3D_global
#
# Author: Angeline G. Burrell, UMichigan, Oct 2013
#
# Comments: Routines to plot lat/lon data as scatter or contours.  Earth
#           continental outlines, solar terminator shading, and the geomagnetic
#           equator may also be added.  The mpl_toolkit baseline must also
#           be installed (available on fink, macports, or online at
#           http://matplotlib.org/basemap).
#
# Includes: plot_single_3D_image          - plots a single rectangular or polar
#                                           alt slice as a contour or scatter
#           plot_single_nsglobal_3D_image - plots northern and southern polar
#                                           altitude slice as a contour or
#                                           scatter
#           plot_global_3D_snapshot       - plots northern and southern polar
#                                           and a midlatitude rectangular
#                                           altitude slice as a contour or
#                                           scatter
#           plot_mult_3D_slices           - plot multiple polar or rectangular
#                                           altitude slices  as a contour or
#                                           scatter
#           -----------------------------------------------------------------
#           plot_snapshot_subfigure       - plots the northern and southern
#                                           polar caps along with a midlat
#                                           rectuangular section
#           plot_nsglobal_subfigure       - plots the northern and southern
#                                           polar caps
#           -----------------------------------------------------------------
#           plot_rectangular_3D_global    - plot a rectangular geographic
#                                           contour with or without a map as a
#                                           contour or scatter
#           plot_polar_3D_global          - plot a polar geographic contour
#                                           with or without a map as a contour
#                                           or scatter
#----------------------------------------------------------------------------

'''
Plot data in lat/lon contours or color-coded scatter plots to facilitate
GITM/data comparisons with networks of ground-based observations
'''

# Import modules
from mpl_toolkits.basemap import Basemap
import string 
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FormatStrFormatter
import datetime as dt
import gitm_plot_rout as gpr

def plot_single_3D_image(plot_type, lat_data, lon_data, z_data, zname, zscale,
                         zunits, zmax=None, zmin=None, zcolor=None, title=None,
                         figname=None, draw=True, nlat=90, slat=-90, linc=6,
                         tlon=90, data_type="scatter", bcolor="#D7DBE0",
                         meq=False, earth=False, m=None, faspect=True,
                         term_datetime=False, extra_line=False,
                         *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range.
    Input: plot_type = key to determine plot type (rectangular, polar)
           lat_data  = numpy array containing the latitudes (in degrees North)
                       This can be 1D or 2D for scatter plots and must be 2D
                       for contours
           lon_data  = numpy array contaiing the longitudes (in degrees East)
                       This can be 1D or 2D for scatter plots and must be 2D
                       for contours
           z_data    = numpy array containing the z variable data
                       This should be 1D for scatter plots and 2D for contours
           zname     = Descriptive name of z variable data
           zscale    = Scaling type (linear/exponential) for z variable data
           zunits    = Descriptive units for z variable data
           zmin      = minimum value for z variable (default=None)
           zmax      = maximum value for z variable (default=None)
           zcolor    = Color map for the z variable.  If not specified, will be
                        determined by the z range. (default=None)
           title     = plot title
           figname   = file name to save figure as (default is none)
           draw      = draw to screen? (default is True)
           nlat      = northern latitude limit (degrees North, default 90)
           slat      = southern latitude limit (degrees North, defalut -90)
           linc      = number of latitude tick incriments (default 6)
           tlon      = longitude at the top of a polar dial (degrees east,
                       default 90)
           data_type     = scatter or contour (default=scatter)
           bcolor        = Background color, or None (default="#D7DBE0")
           meq           = Add the geomagnetic equator? (default=False)
           earth         = include continent outlines for Earth (default False)
           m             = Earth map handle (default=None)
           faspect       = Fix the aspect of Earth (default=True)
           term_datetime = Include the solar terminator by shading the night
                           time regions?  If so, include a datetime object
                           with the UT for this map. (default=False)
           extra_line    = Plot a specified line in (default=False) eg:
                           [x np.array, y np.array, style string (eg 'k-')]
                           where x is in degrees longitude and y is in 
                           degrees latitude
                           
    Output: f = figure handle
            m = map handle
    '''

    # Initialize the z variable limits
    if(zmin is None):
        zmin = np.nanmin(z_data)
    if(zmax is None):
        zmax = np.nanmax(z_data)

    zran = round((zmax-zmin)/6.0)
    if(zran != 0.0 and zscale.find("exp") < 0):
        zmin = math.floor(float("{:.14f}".format(zmin / zran))) * zran
        zmax = math.ceil(float("{:.14f}".format(zmax / zran))) * zran

    if zmax == zmin:
        # Exit gracefully
        print "plot_single_3D_image ERROR: z variable at a constant value", zmax
        return

    # Initialize the color scheme, if desired
    if zcolor is None:
        if zmin < 0.0 and zmax > 0.0:
            zcolor = "seismic_r"
            # Change zrange to be centered, making white = 0.0
            if abs(zmin) > zmax:
                zmax = abs(zmin)
            else:
                zmin = -zmax
        else:
            zcolor = "Spectral_r"

    # Initialize the new figure
    gf = True
    pf = False

    if(string.lower(plot_type)=="polar" and not earth):
        pf = True

    f = plt.figure()
    ax = f.add_subplot(111, polar=pf)

    if(string.lower(plot_type)=="rectangular"):
        (con, m) = plot_rectangular_3D_global(ax, lat_data, lon_data, z_data,
                                              zname, zscale, zunits, zmin, zmax,
                                              zcolor, nlat=nlat, slat=slat,
                                              linc=linc, title=title,
                                              bcolor=bcolor, meq=meq,
                                              earth=earth, m=m, faspect=faspect,
                                              data_type=data_type,
                                              term_datetime=term_datetime)
    elif(string.lower(plot_type)=="polar"):
        (con, m) = plot_polar_3D_global(ax, 1, lat_data, lon_data, z_data,
                                        zname, zscale, zunits, zmin, zmax,
                                        zcolor, center_lat=nlat, edge_lat=slat,
                                        linc=linc, top_lon=tlon, title=title,
                                        bcolor=bcolor, earth=earth, m=m,
                                        data_type=data_type, meq=meq,
                                        term_datetime=term_datetime,
                                        extra_line=extra_line)
    else:
        print "ERROR: unknown input type [", plot_type, "]\n"
        return

    if draw:
        # Draw to screen.
        if plt.isinteractive():
            plt.draw() #In interactive mode, you just "draw".
        else:
            # W/o interactive mode, "show" stops the user from typing more 
            # at the terminal until plots are drawn.
            plt.show()

    # Save output file
    if figname is not None:
        plt.savefig(figname)

    return f, ax, m

def plot_single_nsglobal_3D_image(lat_data, lon_data, z_data, zname, zscale,
                                  zunits, zmax=None, zmin=None,
                                  zcolor=None, title=None, figname=None,
                                  draw=True, plat=90, elat=0, linc=3,
                                  tlon=90, bcolor="#D7DBE0", earth=False,
                                  mn=None, ms=None, data_type="scatter",
                                  term_datetime=False, extra_line_n=False,
                                  extra_line_s=False, *args, **kwargs):
    '''
    Creates a figure with two polar map projections for the northern and 
    southern ends of a specified latitude range.
    Input: lat_data      = numpy array containing the latitudes (in degrees
                           North) This can be 1D or 2D for scatter plots and
                           must be 2D for contours
           lon_data      = numpy array contaiing the longitudes (in degrees
                           East). This can be 1D or 2D for scatter plots and
                           must be 2D for contours
           z_data        = numpy array containing the z variable data.  This
                           should be 1D for scatter plots and 2D for contours
           zname         = Descriptive name of z variable data
           zscale        = Scaling type (linear/exponential) for z variable data
           zunits        = Descriptive units for z variable data
           zmax          = maximum z range (default None)
           zmin          = mininimum z range (default None)
           zcolor        = Color scale for plotting the z data.  If not
                           specified, will be determined by the z range.
                           (default=None)
           title         = plot title
           figname       = file name to save figure as (default is none)
           draw          = draw to screen? (default is True)
           plat          = polar latitude limit (degrees North, default +/-90)
           elat          = equatorial latitude limit (degrees North, defalut 0)
           linc          = number of latitude tick incriments (default 6)
           tlon          = longitude to place on the polar dial top (degrees
                           east, default 90)
           bcolor        = Background color, or None (default="#D7DBE0")
           earth         = include Earth continent outlines (default False)
           mn            = Northern latitude map handle (default=None)
           ms            = Southern latitude map handle (default=None_
           data_type     = Make scatter/contour plot (default=scatter)
           term_datetime = Include the solar terminator by shading the night
                           time regions?  If so, include a datetime object
                           with the UT for this map. (defualt=False)
           extra_line_n  = Plot a specified line in north (default=False) eg:
                           [x np.array, y np.array, style string (eg 'k-')]
                           where x is in degrees longitude and y is in 
                           degrees latitude
           extra_line_s  = Plot a specified line in south (default=False) eg:
                           [x np.array, y np.array, style string (eg 'k-')]
                           where x is in degrees longitude and y is in 
                           degrees latitude

    Output: f         = figure handle
            handle[1] = Northern latitude map handle
            handle[3] = Southern latitude map handle
    '''

    # Initialize the new figure
    f  = plt.figure()

    if title:
        f.suptitle(title, size="medium")

    # Add the northern and southern polar dials
    handle = plot_nsglobal_subfigure(f, 1, 0, lat_data, lon_data, z_data, zname,
                                     zscale, zunits, zmax, zmin, zcolor,
                                     title=True, cb=True, plat=plat, elat=elat,
                                     linc=linc, tlon=tlon, bcolor=bcolor,
                                     earth=earth, mn=mn, ms=ms,
                                     data_type=data_type,
                                     term_datetime=term_datetime,
                                     extra_line_n=extra_line_n,
                                     extra_line_s=extra_line_s)

    if handle:
        # Output the plot
        if draw:
            # Draw to screen.
            if plt.isinteractive():
                plt.draw() #In interactive mode, you just "draw".
            else:
                # W/o interactive mode, "show" stops the user from typing more 
                # at the terminal until plots are drawn.
                plt.show()

        # Save output file
        if figname is not None:
            plt.savefig(figname)

        return(f, handle[1], handle[3])
    else:
        return
# End plot_single_nsglobal_3D_image

def plot_global_3D_snapshot(lat_data, lon_data, z_data, zname, zscale, zunits,
                            zmax=None, zmin=None, zcolor=None, title=None,
                            figname=None, draw=True, tlon=90, polar_blat=45,
                            rect_blat=45, bcolor="#D7DBE0", meq=False,
                            earth=False, ml=None, mn=None, ms=None,
                            data_type="scatter", term_datetime=False,
                            *args, **kwargs):
    '''
    Creates a map projection plot for the entire globe, seperating the polar
    and central latitude regions.
    Input: lat_data      = numpy array containing the latitudes (in degrees
                           North). This can be 1D or 2D for scatter plots and 
                           must be 2D for contours.  This should be 1D for
                           scatter plots and 2D for contours
           lon_data      = numpy array contaiing the longitudes (in degrees
                           East) This can be 1D or 2D for scatter plots and
                           must be 2D for contours
           z_data        = numpy array containing the z variable data.  This
                           should be 1D for scatter plots and 2D for contours
           zname         = Descriptive name of z variable data
           zscale        = Scaling type (linear/exponential) for z variable
           zunits        = Descriptive units for z variable data
           zmax          = maximum z limit (default None)
           zmin          = minimum z limit (default None)
           zcolor        = Color scale for z variable.  If not specified, will
                           be determined by the z range. (default=None)
           title         = plot title
           figname       = file name to save figure as (default is none)
           draw          = output a screen image? (default is True)
           nlat          = northern latitude limit (degrees North, default 90)
           slat          = southern latitude limit (degrees North, defalut 90)
           tlon          = longitude at the top of the polar dial (degrees East,
                           default 90)
           polar_blat    = co-latitude of the lower boundary of the polar dials
                           (default 45)
           rect_blat     = Upper bounding co-latitude of the rectangular map
                           (default 45)
           bcolor        = Background color, or None (default="#D7DBE0")
           meq           = Add the geomagnetic equator? (default=False)
           earth         = include Earth continent outlines (default False)
           ml            = Low latitude map handle (default=None)
           mn            = Norhtern latitude map handle (default=None)
           ms            = Southern latitude map handle (default=None)
           data_type     = Make a scatter/contour plot (default=scatter)
           term_datetime = Include the solar terminator by shading the night
                           time regions?  If so, include a datetime object
                           with the UT for this map.  (default=False)
    Output: f         = figure handle
            handle[1] = Lower latitude map handle
            handle[3] = Northern latitude map handle
            handle[5] = Southern latitude map handle
    '''

    # Initialize the new figure, starting with the mid- and low-latitudes
    f = plt.figure()

    # Call the snapshot subroutine for a subplot number of one
    handles = plot_snapshot_subfigure(f, 1, 0, lat_data, lon_data, z_data,
                                      zname, zscale, zunits, zmax=zmax,
                                      zmin=zmin, zcolor=zcolor, cb=True,
                                      tlon=tlon, polar_blat=polar_blat,
                                      rect_blat=rect_blat, bcolor=bcolor,
                                      meq=meq, earth=earth, ml=ml, mn=mn, ms=ms,
                                      data_type=data_type,
                                      term_datetime=term_datetime)

    if handles:
        if title:
            f.suptitle(title, size="medium")

        if draw:
            # Draw to screen.
            if plt.isinteractive():
                plt.draw() #In interactive mode, you just "draw".
            else:
                # W/o interactive mode, "show" stops the user from typing more 
                # at the terminal until plots are drawn.
                plt.show()

        # Save output file
        if figname is not None:
            plt.savefig(figname)

        return(f, handles[1], handles[3], handles[5])
    else:
        return

# End snapshot

def plot_mult_3D_slices(plot_type, isub, subindex, lat_data, lon_data, z_data,
                        zname, zscale, zunits, zmax=None, zmin=None,
                        zcolor=None, title=None, figname=None, draw=True,
                        nlat=90, slat=-90, linc=6, tlon=90, data_type="scatter",
                        bcolor="#D7DBE0", meq=False, earth=False, faspect=True,
                        term_datetime=False, extra_line=False, *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range.
    Input: plot_type = key to determine plot type (rectangular, polar)
           isub      = index to interate over for subplots (in the data arrays)
           sub_index = list containing the indexes to include in subplots
           lat_data  = numpy array containing the latitudes (in degrees North)
                       for each of the desired subplots. This can be 2D or 3D
                       for scatter plots and must be 3D for contours.  isub
                       indicates which dimension to hold constant when making
                       each subplot.
           lon_data  = numpy array contaiing the longitudes (in degrees East)
                       for each of the desired subplots. This can be 2D or 3D
                       for scatter plots and must be 3D for contours.  isub
                       indicates which dimension to hold constant when making
                       each subplot.
           z_data    = numpy array containing the z variable data for each of
                       the desired subplots. This can be 2D or 3D for scatter
                       plots and must be 3D for contours.  isub indicates which
                       dimension to hold constant when making each subplot.
           zname     = Descriptive name of z variable data
           zscale    = Scaling type (linear/exponential) for z variable data
           zunits    = Descriptive units for z variable data
           zmax      = maximum z range (default None)
           zmin      = minimum z range (default None)
           zcolor    = Color spectrum for z data.  If not specified, will be
                       determined by the z range. (default=None)
           title     = plot title
           figname   = file name to save figure as (default is none)
           draw      = draw to screen? (default is True)
           nlat      = northern latitude limit (degrees North, default 90)
           slat      = southern latitude limit (degrees North, defalut -90)
           linc      = number of latitude tick incriments (default 6)
           tlon          = longitude on the top, for polar plots (degrees East,
                           default 90)
           data_type     = Contour or scatter plot? (default=scatter)
           bcolor        = Background color, or None (default="#D7DBE0")
           meq           = Add the geomagnetic equator? (default=False)
           earth         = include Earth continent outlines (default False)
           faspect       = Fix aspect ratio if using continents (default=True)
           term_datetime = Include the solar terminator by shading the night
                           time regions?  If so, include a datetime object
                           with the UT for this map.  (default=False)
           extra_line    = Plot a specified line in (default=False) eg:
                           [x np.array, y np.array, style string (eg 'k-')]
                           where x is in degrees longitude and y is in 
                           degrees latitude
    '''

    m = None

    # Initialize the z variable limits
    if(zmin is None):
        zmin = np.nanmin(z_data)
    if(zmax is None):
        zmax = np.nanmax(z_data)

    zran = round((zmax-zmin)/6.0)
    if(zran != 0.0 and zscale.find("exp") < 0):
        zmin = math.floor(float("{:.14f}".format(zmin / zran))) * zran
        zmax = math.ceil(float("{:.14f}".format(zmax / zran))) * zran

    if zmax == zmin:
        # Exit gracefully
        print "plot_mult_3D_slices ERROR: z variable at a constant value", zmax
        return

    # Initialize the color scheme, if desired
    if zcolor is None:
        if zmin < 0.0 and zmax > 0.0:
            zcolor = "seismic_r"
            # Change zrange to be centered, making white = 0.0
            if abs(zmin) > zmax:
                zmax = abs(zmin)
            else:
                zmin = -zmax
        else:
            zcolor = "Spectral_r"

    # Initialize the new figure
    pf = False
    sn = len(subindex)

    if(sn < 1):
        print "plot_mult_3D_slices ERROR: no altitude slices specified"
        return

    if(string.lower(plot_type)=="polar" and not earth):
        pf = True

    f  = plt.figure()
    ax = list()
    tl = " "

    if title:
        f.suptitle(title, size="medium")

    # Adjust the figure height to accomadate the number of subplots

    if(sn > 2):
        fheight = f.get_figheight()
        f.set_figheight(fheight * 0.5 * sn)

    if(sn > 1 and string.lower(plot_type) == "polar"):
        fwidth = f.get_figwidth()
        f.set_figwidth(fwidth / 1.5)

    for snum in reversed(range(0, sn)):
        cl = False
        xl = False
        yl = False
        fnum = (sn * 100) + 11 + snum
        si = subindex[snum]
        ax = f.add_subplot(fnum, polar=pf)

        if(sn == snum + 1):
            xl = True

        if(math.floor(sn * 0.5) == snum):
            yl = True

        if(snum == sn - 1):
            cl = True

        if isub == 0:
            latdat = np.array(lat_data[si])
            londat = np.array(lon_data[si])
            zdat = np.array(z_data[si])
        elif isub == 1:
            latdat = np.array(lat_data[:,si])
            londat = np.array(lon_data[:,si])
            zdat = np.array(z_data[:,si])
        elif isub == 2:
            latdat = np.array(lat_data[:,:,si])
            londat = np.array(lon_data[:,:,si])
            zdat = np.array(z_data[:,:,si])
        else:
            print "plot_mult_3D_slices ERROR: subplot index out of range"
            return

        if(string.lower(plot_type)=="rectangular"):
            (con, m) = plot_rectangular_3D_global(ax, latdat, londat, zdat,
                                                  zname, zscale, zunits, zmin,
                                                  zmax, zcolor=zcolor,
                                                  nlat=nlat, slat=slat,
                                                  linc=linc, cb=False, cloc="t",
                                                  title=tl, tloc="r", xl=xl,
                                                  yl=yl, bcolor=bcolor, meq=meq,
                                                  earth=earth, m=m,
                                                  data_type=data_type,
                                                  faspect=faspect,
                                                  term_datetime=term_datetime)
        elif(string.lower(plot_type)=="polar"):
            # Because the polar plots take up half the width as the rectangular
            # plots, decrease the level of latitude incriments by half
            # and adjust the y label output pad
            (con, m) = plot_polar_3D_global(ax, 2, latdat, londat, zdat, zname,
                                            zscale, zunits, zmin, zmax, zcolor,
                                            center_lat=nlat, edge_lat=slat,
                                            linc=(linc/2), top_lon=tlon,
                                            cb=False, cloc="t", title=tl,
                                            tloc="l", tl=xl, rl=yl,
                                            bcolor=bcolor, earth=earth, m=m,
                                            data_type=data_type,
                                            term_datetime=term_datetime,
                                            extra_line=extra_line)
        else:
            print "ERROR: unknown input type [", plot_type, "]\n"
            return

        if data_type.find("scatter") >= 0:
            cax = con.axes
        else:
            cax = con.ax


        cpr = list(cax.get_position().bounds)
        if(cl == False):
            cpr[1] = cpr[1] + ydiff
            cax.set_position(cpr)
        else:
            # Add and adjust colorbar
            cbar = gpr.add_colorbar(con, zmin, zmax, 6, "horizontal",
                                    zscale, zname, zunits)

            bp = list(cbar.ax.get_position().bounds)
            cp = list(cax.get_position().bounds)

            cp[1] = cp[1] - bp[3]
            ydiff = cp[1] - cpr[1]
            cp[3] = cpr[3]
            bp[1] = 4.0 * bp[3]
            bp[3] = 0.5 * bp[3]

            cbar.ax.set_position(bp)
            cax.set_position(cp)

    if draw:
        # Draw to screen.
        if plt.isinteractive():
            plt.draw() #In interactive mode, you just "draw".
        else:
            # W/o interactive mode, "show" stops the user from typing more 
            # at the terminal until plots are drawn.
            plt.show()

    # Save output file
    if figname is not None:
        plt.savefig(figname)

    return f
# End plot_mult_3D_slices

def plot_nsglobal_subfigure(f, nsub, isub, lat_data, lon_data, z_data, zname,
                            zscale, zunits, zmax=None, zmin=None,
                            zcolor=None, title=True, cb=True, plat=90, elat=0,
                            linc=3, tlon=90, rl=True, tl=True, bcolor="#D7DBE0",
                            earth=False, mn=None, ms=None, data_type="scatter",
                            term_datetime=False, extra_line_n=False,
                            extra_line_s=False, *args, **kwargs):
    '''
    Creates a figure with two polar map projections for the northern and 
    southern ends of a specified latitude range.
    Input: f             = figure handle
           nsub          = number of subplots to include
           isub          = number of this subplot (zero offset)
           lat_data      = numpy array containing the latitudes (in degrees
                           North) This can be 1D or 2D for scatter plots and
                           must be 2D for contours
           lon_data      = numpy array contaiing the longitudes (degrees East)
                           This can be 1D or 2D for scatter plots and must be 2D
                           for contours
           z_data        = numpy array containing the z variable data.  This
                           should be 1D for scatter plots and 2D for contours
           zname         = Descriptive name of z variable data
           zscale        = Scaling type (linear/exponential) for z variable data
           zunits        = Descriptive units for z variable data
           zmax          = maximum z range (default None)
           zmin          = mininimum z range (default None)
           zcolor        = Color scale for plotting the z data.  If not
                           specified, will be determined by the z range.
                           (default=None)
           title         = Include titles indicating which polar dial is North
                           and South? (default=True)
           cb            = Include colorbar? (default=True)
           plat          = polar latitude limit (degrees North, default +/-90)
           elat          = equatorial latitude limit (degrees North, defalut 0)
           linc          = number of latitude tick incriments (default 6)
           tlon          = longitude to place on the polar dial top (degrees
                           east, default 90)
           bcolor        = Background color, or None (default="#D7DBE0")
           earth         = include Earth continent outlines (default False)
           mn            = Northern latitude map handle (default=None)
           ms            = Southern latitude map handle (default=None)
           data_type     = Plot type: scatter/contour (default=scatter)
           term_datetime = Include the solar terminator by shading the night
                           time regions?  If so, include a datetime object
                           with the UT for this map.  (default=False)
           extra_line_n  = Plot a specified line in north (default=False) eg:
                           [x np.array, y np.array, style string (eg 'k-')]
                           where x is in degrees longitude and y is in 
                           degrees latitude
           extra_line_s  = Plot a specified line in south (default=False) eg:
                           [x np.array, y np.array, style string (eg 'k-')]
                           where x is in degrees longitude and y is in 
                           degrees latitude
    '''
    # Initialize the z variable limits
    if(zmin is None):
        zmin = np.nanmin(z_data)
    if(zmax is None):
        zmax = np.nanmax(z_data)

    zran = round((zmax-zmin)/6.0)
    if(zran != 0.0 and zscale.find("exp") < 0):
        zmin = math.floor(float("{:.14f}".format(zmin / zran))) * zran
        zmax = math.ceil(float("{:.14f}".format(zmax / zran))) * zran

    if zmax == zmin:
        # Exit gracefully
        print "plot_nsglobal_subfigure ERROR: z variable at a const value", zmax
        return

    # Initialize the color scheme, if desired
    if zcolor is None:
        if zmin < 0.0 and zmax > 0.0:
            zcolor = "seismic_r"
            # Change zrange to be centered, making white = 0.0
            if abs(zmin) > zmax:
                zmax = abs(zmin)
            else:
                zmin = -zmax
        else:
            zcolor = "Spectral_r"

    # Initialize the flags for the new subfigure
    pf = True
    nc = False
    sc = True

    if earth:
        pf = False
        nc = True
        sc = False

    if not rl:
        nc = False
        sc = False

    rsub = 2 * isub + 1

    ntitle = False
    stitle = False

    if title:
        ntitle = "North"
        stitle = "South"

    # Northern Plot
    axn = f.add_subplot(nsub,2,rsub, polar=pf)
    (con, mn) = plot_polar_3D_global(axn, 2, lat_data, lon_data, z_data, zname,
                                     zscale, zunits, zmin, zmax, zcolor,
                                     center_lat=plat, edge_lat=elat, linc=linc,
                                     top_lon=tlon, cb=cb, cloc="l",
                                     title=ntitle, rl=nc, tl=tl, bcolor=bcolor,
                                     earth=earth, m=mn, data_type=data_type,
                                     term_datetime=term_datetime,
                                     extra_line=extra_line_n)
    psn = list(axn.get_position().bounds)

    # Southern Plot
    rsub += 1
    plat = -plat
    if(elat != 0.0):
        elat = -elat

    axs = f.add_subplot(nsub,2,rsub, polar=pf)
    (con, ms) = plot_polar_3D_global(axs, 2, lat_data, lon_data, z_data, zname,
                                     zscale, zunits, zmin, zmax, zcolor,
                                     center_lat=plat, edge_lat=elat, linc=linc,
                                     top_lon=tlon, cb=False, title=stitle,
                                     rl=sc, tl=tl, bcolor=bcolor, earth=earth,
                                     m=ms, data_type=data_type,
                                     term_datetime=term_datetime,
                                     extra_line=extra_line_s)
    pss = list(axs.get_position().bounds)

    # Adjust plot sizes
    pwidth = 0.5 * (psn[2] + pss[2])
    if earth:
        psn[0] = psn[0] + (pwidth - psn[2]) * 0.75 - 0.05
    else:
        psn[0] = psn[0] - (pwidth - psn[2]) * 0.25
    psn[2] = pwidth
    pss[0] = psn[0] + psn[2] + pwidth * 0.2
    pss[2] = pwidth

    if earth:
        pss[0] = pss[0] + 0.05

    axs.set_position(pss)
    axn.set_position(psn)    

    return(axn, mn, axs, ms)
# End plot_nsglobal_subfigure

def plot_snapshot_subfigure(f, nsub, isub, lat_data, lon_data, z_data, zname,
                            zscale, zunits, zmax=None, zmin=None, zcolor=None,
                            cb=True, tlon=90, polar_blat=45, rect_blat=45,
                            title=True, xl=True, yl=True, xt=True, yt=True, 
                            bcolor="#D7DBE0", meq=False, earth=False, ml=None,
                            mn=None, ms=None, data_type="scatter",
                            term_datetime=False, *args, **kwargs):
    '''
    Creates a map projection plot for the entire globe, seperating the polar
    and central latitude regions.
    Input: f             = figure handle
           nsub          = number of subplots to include
           isub          = number of this subplot (zero offset)
           lat_data      = numpy array containing the latitudes (in degrees
                           North) This can be 1D or 2D for scatter plots and
                           must be 2D for contours
           lon_data      = numpy array contaiing the longitudes (in degrees
                           East) This can be 1D or 2D for scatter plots and
                           must be 2D for contours
           z_data        = numpy array containing the z variable data.  This
                           should be 1D for scatter plots and 2D for contours
           zname         = Descriptive name of z variable data
           zscale        = Scaling type (linear/exponential) for z variable data
           zunits        = Descriptive units for z variable data
           zmax          = maximum z limit (default None)
           zmin          = minimum z limit (default None)
           zcolor        = Color scale for z data.  If not specified, will be
                           determined by the z range. (default=None)
           cb            = Include colorbar? (default=True)
           tlon          = longitude at the top of the polar dial (degrees East,
                           default 90)
           polar_blat    = co-latitude of the lower boundary of the polar dials
                           (default 45)
           rect_blat     = Upper bounding co-latitude of the rectangular map
                           (default 45)
           title         = Include North/South polar titles? (default=True)
           xl            = Include x (Longitude) label? (default=True)
           yl            = Include y (Latitude) label? (default=True)
           xt            = Include x (Longitude) ticks? (default=True)
           yt            = Include y (Latitude) ticks? (default=True)
           bcolor        = Background color, or None (default="#D7DBE0")
           meq           = Add the geomagnetic equator? (default=False)
           earth         = include Earth continent outlines (default False)
           ml            = Map handle for low latitudes (default=None)
           mn            = Map handle for northern latiutdes (default=None)
           ms            = Map handle for southern latitudes (default=None)
           term_datetime = Include the solar terminator by shading the night
                           time regions?  If so, include a datetime object
                           with the UT for this map.  (default=False)

           Output: axl = Plot handle for low latitudes
                   ml  = Map handle for low latitudes
                   axn = Plot handle for northern latitudes
                   mn  = Map handle for northern latitudes
                   axs = Plot handle for southern latitudes
                   ms  = Map handle for southern latitudes
           
    '''

    # Initialize the z variable limits
    if(zmin is None):
        zmin = np.nanmin(z_data)
    if(zmax is None):
        zmax = np.nanmax(z_data)

    zran = round((zmax-zmin)/6.0)
    if(zran != 0.0 and zscale.find("exp") < 0):
        zmin = math.floor(float("{:.14f}".format(zmin / zran))) * zran
        zmax = math.ceil(float("{:.14f}".format(zmax / zran))) * zran

    if zmax == zmin:
        # Exit gracefully
        print "plot_snapshot_subfigure ERROR: z variable at a const value", zmax
        return

    # Initialize the color scheme, if desired
    if zcolor is None:
        if zmin < 0.0 and zmax > 0.0:
            zcolor = "seismic_r"
            # Change zrange to be centered, making white = 0.0
            if abs(zmin) > zmax:
                zmax = abs(zmin)
            else:
                zmin = -zmax
        else:
            zcolor = "Spectral_r"

    # Initialize the subplot numbers. Snapshot uses two subplot spaces per plot
    lsub = 2 * nsub
    rsub = 2 * (isub + 1)

    # Add the mid- and low-latitude rectuangular plot
    axl = f.add_subplot(lsub,1,rsub)
    (con, ml) = plot_rectangular_3D_global(axl, lat_data, lon_data, z_data,
                                           zname, zscale, zunits, zmin, zmax,
                                           zcolor, nlat=rect_blat,
                                           slat=-1.0*rect_blat, cb=False, xl=xl,
                                           xt=xt, yl=yl, yt=yt, bcolor=bcolor,
                                           meq=meq, earth=earth, m=ml,
                                           data_type=data_type,
                                           term_datetime=term_datetime)
    psl = list(axl.get_position().bounds)

    # Add the North pole
    pf = True
    if earth:
        pf = False

    rsub = 4 * isub + 1 # polar plots exist in 2 nsub x 2 space

    axn = f.add_subplot(lsub,2,rsub, polar=pf)
    if title:
        t="North"
    else:
        t=None

    if nsub > 1:
        xl = False
        yl = False

    (con, mn) = plot_polar_3D_global(axn, 2, lat_data, lon_data, z_data, zname,
                                     zscale, zunits, zmin, zmax, zcolor,
                                     center_lat=90, edge_lat=polar_blat, linc=3,
                                     top_lon=tlon, cb=False, title=t,
                                     tl=False, rl=False, bcolor=bcolor,
                                     earth=earth, m=mn, data_type=data_type,
                                     term_datetime=term_datetime)
    psn = list(axn.get_position().bounds)

    # Add the South pole
    rsub += 1
    axs = f.add_subplot(lsub,2,rsub, polar=pf)
    if title:
        t="South"
    (con, ms) = plot_polar_3D_global(axs, 2, lat_data, lon_data, z_data, zname,
                                     zscale, zunits, zmin, zmax, zcolor,
                                     center_lat=-90, edge_lat=-1.0*polar_blat,
                                     linc=3, top_lon=tlon, cb=False, title=t,
                                     tl=False, rl=False, bcolor=bcolor,
                                     earth=earth, m=ms, data_type=data_type,
                                     term_datetime=term_datetime)
    pss = list(axs.get_position().bounds)

    # Add a colorbar for the entire plot if desired
    if cb:
        cb = gpr.add_colorbar(con, zmin, zmax, 6, 'vertical', zscale, zname,
                              zunits)
        bp = list(cb.ax.get_position().bounds)

        # Ridley Adjustment of plots and colorbar
        bp[0] = .82
        bp[1] = psl[1] - 0.02
        bp[3] = pss[3] + .85
        bp[2] = bp[2] / 1.8 

        if earth is True:
            pss[0] = pss[0] - .08
        else:
            pss[0] = pss[0] - .12
        pss[1] = pss[1] - .02

        psn[0] = psn[0] - .05
        psn[1] = psn[1] - .02

        #AGB Adjustment of Ridley alignment
        if earth:
            pss[0] = pss[0] - .05
            if rect_blat > 70.0: psl[0] = psl[0] - .09
        else:
            psl[0] = psl[0] - .02

        cb.ax.set_position(bp)
    else:
    # AGB plot alignments without a colorbar
        pss[0] = pss[0] - 0.08
        psn[0] = psn[0] - 0.08
        psl[0] = psl[0] - 0.02

    # Adjust the plot alignments
    axs.set_position(pss)
    axl.set_position(psl)
    axn.set_position(psn)

    return(axl,ml,axn,mn,axs,ms)
# End snapshot

def plot_rectangular_3D_global(ax, lat_data, lon_data, z_data, zname, zscale,
                               zunits, zmin, zmax, zcolor=None, zinc=6,
                               nlat=90, slat=-90, linc=6, cb=True, cloc="r",
                               title=None, tloc="t", xl=True, xt=True, yl=True,
                               yt=True, bcolor="#D7DBE0", meq=False,
                               earth=False, m=None, data_type="scatter",
                               faspect=True, term_datetime=False,
                               *args, **kwargs):
    '''
    Creates a rectangular map projection plot for a specified latitude range.
    Input: ax        = axis handle
           lat_data  = numpy array containing the latitudes (in degrees North)
                       This can be 1D or 2D for scatter plots and must be 2D
                       for contours
           lon_data  = numpy array contaiing the longitudes (in degrees East)
                       This can be 1D or 2D for scatter plots and must be 2D
                       for contours
           z_data    = numpy array containing the z variable data
                       This should be 1D for scatter plots and 2D for contours
           zname     = Descriptive name of z variable data
           zscale    = Scaling type (linear/exponential) for z variable data
           zunits    = Descriptive units for z variable data
           zmin      = minimum value for z variable
           zmax      = maximum value for z variable
           zcolor    = Color map for the z variable.  If not specified, will be
                       determined by the z range. (default=None) 
           zinc      = number of tick incriments for z variable (default 6)
           nlat      = northern latitude limit (degrees North, default 90)
           slat      = southern latitude limit (degrees North, defalut 90)
           linc      = number of tick incriments in latitude (default 6)
           cb        = Add a colorbar (default is True)
           cloc      = Specify the colorbar location (t=top, r=right, l=left,
                       b=bottom, default is right)
           title     = plot title (default is None)
           tloc      = Specify the title location (t=top, r=right, l=left,
                       b=bottom, default is top)
           xl        = Include x (longitude) label (default is True)
           xt        = Include x tics (default=True)
           yl        = Include y (latitude) label (default is True)
           yt        = Include y tics (default=True)
           bcolor    = Background color, or None (default="#D7DBE0")
           meq       = Add a line for the geomagnetic equator? (default=False)
           earth     = Include Earth continent outlines (default is False)
           m         = Earth map handle (default=None)
           data_type = scatter or contour (default=scatter)
           faspect   = Fix the aspect of Earth if using outlines (default=True)
           term_datetime = Include the solar terminator by shading the night
                           time regions?  If so, include a datetime object
                           with the UT for this map.  (default=False)

    Output: con = Contour handle
            m   = Earth map handle
    '''

    # Exit gracefully, if necessary
    rout_name = "plot_rectangular_3D_global"

    if(nlat == slat):
        print rout_name, "ERROR: no latitude range"
        return

    if zmax == zmin:
        print rout_name, "ERROR: z variable at a constant value", zmax
        return

    # Initialize the color scheme, if desired
    if zcolor is None:
        if zmin < 0.0 and zmax > 0.0:
            zcolor = "seismic_r"
            # Change zrange to be centered, making white = 0.0
            if abs(zmin) > zmax:
                zmax = abs(zmin)
            else:
                zmin = -zmax
        else:
            zcolor = "Spectral_r"

    # Set latitude range values
    yrange = nlat - slat
    ywidth = yrange / linc

    # If desired, map the Earth using the Equidistant Cylindrical Projection
    if earth:
        if not m:
            m = Basemap(projection='cyl',llcrnrlat=slat,urcrnrlat=nlat,
                        llcrnrlon=0,urcrnrlon=360,resolution='c',
                        fix_aspect=faspect)
        m.drawcoastlines(linewidth=0.5)
        if bcolor is not None:
            m.drawmapboundary(fill_color=bcolor)
        m.drawparallels(np.arange(slat, nlat+1, ywidth))
        m.drawmeridians(np.arange(0,361,60))

        if type(term_datetime) is dt.datetime:
            m.nightshade(term_datetime, alpha=.2)
    elif bcolor is not None:
        # Change the background color
        ax.patch.set_facecolor(bcolor)

    # Determine whether the z scale is linear or exponential
    if zscale.find("exp") >= 0:
        v = np.logspace(math.log10(zmin), math.log10(zmax), zinc*10,
                        endpoint=True)
        norm = LogNorm(vmin=zmin, vmax=zmax)
    else:
        v = np.linspace(zmin, zmax, zinc*10, endpoint=True)
        norm = None

    # Plot the data
    if data_type.find("con") >= 0:
        con = ax.contourf(lon_data, lat_data, z_data, v, cmap=get_cmap(zcolor),
                          norm=norm, vmin=zmin, vmax=zmax)
        cax = con.ax
    else:
        con = ax.scatter(lon_data, lat_data, c=z_data, cmap=get_cmap(zcolor),
                         norm=norm, vmin=zmin, vmax=zmax, edgecolors="none",
                         s=15)
        cax = con.axes

    if not earth and type(term_datetime) is dt.datetime:
        # Add solar terminator, if desired
        try:
            gpr.add_solar_terminator(term_datetime, ax=ax)
            gpr.add_subsolar_point(term_datetime, ax=ax, style="*")
        except:
            print "Unable to add solar terminator without PySolar"


    # Configure axis
    ytics  = MultipleLocator(ywidth)
    ax.yaxis.set_major_locator(ytics)
    if not yt:
        ax.yaxis.set_major_formatter(FormatStrFormatter(""))
    else:
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f"))

    if yl:
        ax.set_ylabel('Latitude ($degrees$)')
    plt.ylim(slat, nlat)

    xtics = MultipleLocator(60)
    ax.xaxis.set_major_locator(xtics)
    if not xt:
        ax.xaxis.set_major_formatter(FormatStrFormatter(""))
    if xl:
        ax.set_xlabel('Longitude ($degrees$)')
    plt.xlim(0., 360.)

    # Add a line for the geomagnetic equator, if desired
    if meq:
        gpr.add_geomagnetic_equator(ax)

    # Set the title
    if title:
        rot  = 'horizontal'
        yloc = 1.05
        xloc = .5

        if(tloc == "l" or tloc == "r"):
            xloc = -.1
            yloc = .5
            rot  = 'vertical'

            if(tloc == "r"):
                xloc = 1.1

        if(tloc == "b"):
            yloc = -.1
            
        if title:
            ax.set_title(title, size='medium', rotation=rot, y=yloc, x=xloc)

    # Add a colorbar
    if cb:
        cp = list(cax.get_position().bounds)
        orient = 'vertical'

        if(cloc == 't' or cloc == 'b'):
            orient = 'horizontal'

        cbar = gpr.add_colorbar(con, zmin, zmax, zinc, orient, zscale, zname,
                                zunits)
        bp = list(cbar.ax.get_position().bounds)

        if(cloc == 't'):
            cp[1] = bp[1]
            bp[1] = cp[1] + cp[3] + 0.07
        elif(cloc == 'l'):
            bp[0] = 0.125
            cp[0] = bp[0] + 0.1 + bp[2]
        elif(cloc == 'b'):
            bp[1] = cp[1] - 0.01
            cp[1] = cp[1] + 0.075
        else:
            if earth and faspect:
                cp[0] = bp[0] - cp[2] + 0.005
                bp[0] = bp[0] + 0.04
            else:
                cp[2] = bp[0] - 0.01 - cp[0]

        cax.set_position(cp)
        cbar.ax.set_position(bp)

    return(con, m)

#End plot_rectangular_3D_global

def plot_polar_3D_global(ax, nsub, lat_data, lon_data, z_data, zname, zscale,
                         zunits, zmin, zmax, zcolor=None, zinc=6,
                         center_lat=90, edge_lat=0, linc=6, top_lon=90, cb=True,
                         cloc="r", title = None, tloc="t", tl = True, rl = True,
                         bcolor="#D7DBE0", earth=False, m=None,
                         data_type="scatter", meq=False, term_datetime=False,
                         extra_line=False, *args, **kwargs):
    '''
    Creates a single polar projection, with the latitude center and range
    determined by the input.
    Input: ax         = axis handle
           nsub       = number of vertically stacked subplots
           lat_data   = numpy array containing the latitudes (in degrees North)
                        This can be 1D or 2D for scatter plots and must be 2D
                        for contours.  Latitudes > 90 or < -90 will cause the
                        contour option to crash.
           lon_data   = numpy array contaiing the longitudes (in degrees East)
                        This can be 1D or 2D for scatter plots and must be 2D
                        for contours
           z_data     = numpy array containing the z variable data
                        This should be 1D for scatter plots and 2D for contours
           zname      = Descriptive name of z variable data
           zscale     = Scaling type (linear/exponential) for z variable data
           zunits     = Descriptive units for z variable data
           zmin       = minimum value for z variable
           zmax       = maximum value for z variable
           zcolor     = Color map for the z variable.  If not specified, will be
                        determined by the z range. (default=None)
           zinc       = number of tick incriments for z variable (default 6)
           center_lat = upper (center) latitude limit (degrees N, default 90)
           edge_lat   = lower (edge) latitude limit (degrees N, default 0)
           linc       = number of tick incriments in latitude (default 6)
           top_lon    = longitude to place on top (degrees E, default 90)
           cb         = Add a colorbar (default is True)
           cloc       = Colorbar location (t=top, r=right, l=left, b=bottom, 
                        default is right)
           title      = plot title (default is none)
           tloc       = title location (t=top, r=right, l=left, b=bottom,
                        default is top)
           tl         = include theta (longitude) label (default is True)
           rl         = include r (latitude) label (default is True)
           bcolor     = Background color, or None (default="#D7DBE0")
           earth      = include Earth continent outlines (default is False)
           m          = Earth map handle.  (default=None)
           data_type  = scatter or contour (default=scatter)
           meq           = Include the geomagnetic equator?
           term_datetime = Include the solar terminator by shading the night
                           time regions?  If so, include a datetime object
                           with the UT for this map.  (default=False)
           extra_line    = Plot a specified line (default=False) eg:
                           [x np.array, y np.array, style string (eg 'k-')]
                           where x is in degrees longitude and y is in 
                           degrees latitude

    Output: con = Contour plot handle
            m   = Earth map handle (or None if earth=False and not specified)
    
    '''

    rout_name = "plot_polar_3D_global"

    if zmax == zmin:
        print rout_name, "ERROR: z variable at a constant value", zmax
        return

    # Initialize the color scheme, if desired
    if zcolor is None:
        if zmin < 0.0 and zmax > 0.0:
            zcolor = "seismic_r"
            # Change zrange to be centered, making white = 0.0
            if abs(zmin) > zmax:
                zmax = abs(zmin)
            else:
                zmin = -zmax
        else:
            zcolor = "Spectral_r"

    # Assign the Longitude, Latitude, and Z data structures.  For the polar
    # plots we do not want to include co-latitudes above 90 degrees
    csign = np.sign(center_lat)
    if csign == 0:
        csign = 1.0

    # Set range values
    lon = top_lon - 67.5
    rrange = center_lat - edge_lat
    rwidth = (rrange / linc)

    # Plot the polar contours
    if earth:
        blon = top_lon - 180
        r = lat_data
        theta = lon_data
        # If desired, map the Earth using a Polar Azimuthal Equidistant
        # Projection

        if center_lat < 0.0:
            if not m:
                m = Basemap(projection='spaeqd', lon_0=blon, lat_0=center_lat,
                            boundinglat=edge_lat, round=True, resolution='c')
            lats = list(np.arange(edge_lat, center_lat-1, rwidth))
            lon = lon + 135
        else:
            if not m:
                m = Basemap(projection='npaeqd', lon_0=blon, lat_0=center_lat,
                            boundinglat=edge_lat, round=True, resolution='c')
            lats = list(np.arange(edge_lat, center_lat+1, rwidth))
            lats.reverse()

        m.drawcoastlines(linewidth=0.5)
        if bcolor is not None:
            m.drawmapboundary(fill_color=bcolor)
        llab = m.drawparallels(lats)
        m.drawmeridians(np.arange(0,360,45), labels=[1,1,0,0], labelstyle="+/-")

        if type(term_datetime) is dt.datetime:
            try:
                m.nightshade(term_datetime, alpha=.2)
            except AttributeError:
                print rout_name, "AttributeError in Terminator shading"

        # Determine whether the z scale is linear or exponential
        if zscale.find("exp") >= 0:
            v = np.logspace(math.log10(zmin), math.log10(zmax), zinc*10,
                            endpoint=True)
            norm = LogNorm(vmin=zmin, vmax=zmax)
        else:
            v = np.linspace(zmin, zmax, zinc*10, endpoint=True)
            norm = None

        # Plot the data on top of the map
        if data_type.find("con") >= 0:
            con = m.contourf(theta, r, z_data, v, cmap=get_cmap(zcolor),
                             norm=norm, vmin=zmin, vmax=zmax, latlon=True)
            cax = con.ax
        else:
            con = m.scatter(theta, r, c=z_data, cmap=get_cmap(zcolor),
                            norm=norm, vmin=zmin, vmax=zmax, edgecolors="none",
                            s=15, latlon=True)
            cax = con.axes
        lpad = 20
        x, y = m(list(lon for i in lats), lats)

        for i, label in enumerate(lats):
            ax.text(x[i], y[i], "%.0f$^\circ$" % (label))
 
        # Add magnetic equator, if desired
        if meq is True:
            meq_lon, meq_lat = gpr.add_geomagnetic_equator()
            m.plot(meq_lon, meq_lat, "-k", latlon=True) 

        # Add a specified line, if desired
        if type(extra_line) is list:
            if len(extra_line) >= 2:
                stylestring="k:"
                if len(extra_line) > 2:
                    stylestring = extra_line[2]

                m.plot(extra_line[0], extra_line[1], stylestring, latlon=True)
            else:
                print "Unable to plot extra line, missing data"


    else:
        # Change the background color
        if bcolor is not None:
            ax.patch.set_facecolor(bcolor)

        # Set the contour
        rwidth *= csign
        toff = (450.0 - top_lon) * np.pi / 180.0
        r = csign * lat_data
        theta = lon_data * np.pi / 180.0

        # Determine whether the z scale is linear or exponential
        if zscale.find("exp") >= 0:
            v = np.logspace(math.log10(zmin), math.log10(zmax), zinc*10,
                            endpoint=True)
            norm = LogNorm(vmin=zmin, vmax=zmax)
        else:
            v = np.linspace(zmin, zmax, zinc*10, endpoint=True)
            norm = None

        # plot the data as a scatter or contour plot
        if data_type.find("con") >= 0:
            con = ax.contourf(theta,
                              np.array(gpr.center_polar_cap(abs(center_lat),
                                                            abs(edge_lat),r)),
                              z_data, v, cmap=get_cmap(zcolor), norm=norm,
                              vmin=zmin, vmax=zmax)
            cax = con.ax
        else:
            con = ax.scatter(theta,
                             np.array(gpr.center_polar_cap(abs(center_lat),
                                                           abs(edge_lat),r)),
                             c=z_data, cmap=get_cmap(zcolor), norm=norm,
                             vmin=zmin, vmax=zmax, edgecolors="none", s=15)
            cax = con.axes

        # get minimum value of r for tics
        rmin = ax.get_rmin()

        ax.set_theta_offset(toff)

        if(edge_lat < center_lat):
            ledge = edge_lat
        else:
            ledge = center_lat

        rtics = [rmin + rwidth * (x) for x in range(linc+1)]
        rlabels = map(str, (map(int, [ledge + rwidth * (x)
                                      for x in range(linc+1)])))

        if(center_lat > edge_lat):
            rlabels.reverse()

        if(min(rtics) == 0.0 or max(rtics) == 0.0):
            # set_rgrids cannot handle zero.  Set minimum value to a fraction
            # with the same sign as the maximum value
            i = rtics.index(0.0)
            if(min(rtics) == 0.0):
                j = rtics.index(max(rtics))
            else:
                j = rtics.index(min(rtics))
            rtics[i] = .1 * (rtics[j] / abs(rtics[j]))

        if(min(rtics) > 0.0):
            ax.set_rgrids(rtics, labels=rlabels, angle=lon)
        else:
            rtics = range(-edge_lat, -center_lat, -rwidth)
            ax.set_rgrids(rtics, labels=rlabels, angle=lon)
            
        # Add magnetic equator, if desired
        if meq is True:
            meq_lon, meq_lat = gpr.add_geomagnetic_equator()
            meq_theta = np.array([l * np.pi / 180.0 for l in meq_lon])
            meq_r = np.array([l * csign for l in meq_lat])
            meq_r = np.array(gpr.center_polar_cap(abs(center_lat),
                                                  abs(edge_lat), meq_r))
            ax.plot(meq_theta, meq_r, "-k") 

        # Add solar terminator, if desired
        if type(term_datetime) is dt.datetime:
            try:
                # Add the solar terminator
                (tlon, tlat) = gpr.add_solar_terminator(term_datetime)
                tlon = np.radians(tlon)
                tlat *= csign
                tlat = np.array(gpr.center_polar_cap(abs(center_lat),
                                                     abs(edge_lat), tlat))
                ax.plot(tlon, tlat, "-k")

                # Add the subsolar point
                (sslon, sslat) = gpr.add_subsolar_point(term_datetime)
                sslon *= (np.pi / 180.0)
                sslat = gpr.center_polar_cap(abs(center_lat), abs(edge_lat),
                                             csign * sslat)
                ax.plot([sslon], [sslat], marker="*", ms=5,
                        markerfacecolor='none', markeredgecolor='k')
            except:
                print "Unable to add solar terminator without PySolar"

        if type(extra_line) is list:
            if len(extra_line) >= 2:
                x = np.radians(extra_line[0])
                y = np.array(gpr.center_polar_cap(abs(center_lat),
                                                  abs(edge_lat),
                                                  (csign*extra_line[1])))
                stylestring="k:"
                if len(extra_line) > 2:
                    stylestring = extra_line[2]

                ax.plot(x, y, stylestring)
            else:
                print "Unable to plot extra line, missing data"

        ax.set_rmax(max(rtics))
        ax.set_rmin(min(rtics))
        ax.set_rticks(rtics)

    # Configure axis.
    if tl:
        ax.set_xlabel('Longitude')
        if nsub == 1:
            ax.xaxis.set_label_coords(0.5, 0.0)
        else:
            ax.xaxis.set_label_coords(0.5, -0.1)

    if rl:
        label_string = ""
        if not earth:
            label_string = " ($degrees$)"

        ax.set_ylabel('Latitude{:s}'.format(label_string))
        if earth:
            if nsub == 1:
                ax.yaxis.set_label_coords(1.125, 0.5)
            else:
                ax.yaxis.set_label_coords(1.25, 0.5)
        else:
            ax.yaxis.set_label_coords(1.2, 0.5)

    # Set the title
    if title:
        rot  = 'horizontal'
        yloc = 1.07
        xloc = 0.5

        if tloc == "b":
            yloc = -.1
        elif tloc != "t":
            rot  = 'vertical'
            yloc = 0.5
            xloc = -.2

            if earth:
                xloc = -.3

            if tloc == "r":
                xloc = 1.1
                
                if earth:
                    xloc = 1.2

        # Set the plot title, if desired
        if title:
            ax.set_title(title, size='medium', y=yloc, x=xloc, rotation=rot)

    # Add a colorbar
    if cb:
        cp = list(cax.get_position().bounds)
        orient = 'vertical'

        if(cloc == 't' or cloc == 'b'):
            orient = 'horizontal'

        cbar = gpr.add_colorbar(con, zmin, zmax, zinc, orient,
                                zscale, zname, zunits)
        bp = list(cbar.ax.get_position().bounds)
        cbp = list(cax.get_position().bounds)
        
        if(cloc == 't'):
            cp[1] = bp[1]
            cp[3] = cbp[3]
            bp[1] = cp[1] + cp[3] + 3.0*bp[3]
        elif(cloc == 'l'):
            bp[0] = 0.005
            cp[0] = bp[0] + 0.125 + bp[2]
        elif(cloc == 'r'):
            cp[0] = cp[0] - 0.2
            bp[0] = bp[0] - 0.025

        if(cloc != 'b'):
            cbar.ax.set_position(bp)
            cax.set_position(cp)

    return (con, m)
#End
