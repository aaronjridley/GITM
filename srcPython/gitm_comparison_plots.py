#!/usr/bin/env python
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#-----------------------------------------------------------------------------
# gitm_comparison_plots
#
# Author: Angeline G. Burrell, UMichigan, Oct 2013
#
# Comments: Routine to make GITM/data comparision plots for data from a network
#           of ground-based instruments.  The plots use lat/lon contours or
#           scatter plots as appropriate with data in the third dimension.
#           Map contours are also availabel if you have the mpl_toolkit
#           baseline installed (available on fink, macports, or online at
#           http://matplotlib.org/basemap).
#
# Includes: extract_gitm_time_arrays - Extract non-filler data from GitmTime
#                                      class and construct numpy arrays
#           plot_net_gitm_comp       - plots ground data, gitm data, and the
#                                      difference between the two
#           plot_sat_gitm_comp       - Plots satellite data along with GITM
#                                      data along the satellite track.  The
#                                      difference and percent difference between
#                                      the observations and models is also shown
#----------------------------------------------------------------------------

'''
Process and plot data to facilitate GITM/data comparisons
'''

# Import modules
import math
import numpy as np
from spacepy import pybats
from spacepy.pybats import gitm
import gitm # Temporary until the new GITM is incorporated into spacepy
import gitm_time
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter, FuncFormatter
import matplotlib.pyplot as plt
import datetime as dt
import gitm_plot_rout as gpr
import plot_3D_global as p3g

def extract_data_matched_arrays(dkey, data_dict, bad_value=np.nan):
    '''
    Extract points from matched data arrays that do not have the specified
    bad values.

    Input: dkey      = key to check for bad values
           data_dict = dictionary of lists or numpy arrays with the same
                       number of data values
           bad_value = bad value (default=[np.nan])

    Output: good_dict = dictionary of lists with only good, matched values
    '''
    # Extract the data, removing any instances where the specified data key has
    # the bad values

    data_key = data_dict.keys()
    if np.isnan(bad_value):
        good_index = [i for i,l in enumerate(data_dict[dkey])
                      if not np.isnan(l)]
    elif type(bad_value) is bool:
        good_index = [i for i,l in enumerate(data_dict[dkey])
                      if not bad_value]
    else:
        good_index = [i for i,l in enumerate(data_dict[dkey])
                      if l != bad_value]

    # If not all keys have the same amount of data, remove trailing elements
    data_len = [len(data_dict[k]) for k in data_key]
    min_len = np.min(data_len)

    # Cycle through the rest of the data keys to see if they hold bad values
    if len(data_key) > 1:
        bad_index = list()
        dkey = data_key[1:len(data_key)]

        for i,igood in enumerate(good_index):
	    if igood >= min_len:
	        bad_index.append(i)
	    else:
                if np.isnan(bad_value):
                    badlist = [k for k in dkey if np.isnan(data_dict[k][igood])]
                elif type(bad_value) is bool:
                    badlist = [k for k in dkey if data_dict[k][igood]
                               is bad_value]
                else:
                    badlist = [k for k in dkey if data_dict[k][igood]
                               == bad_value]

                if len(badlist) > 0:
                    bad_index.append(i)
                
        if len(bad_index) > 0:
            good_index = list(np.delete(np.array(good_index), bad_index))

    # Save the good data
    good_data = dict()
    if len(good_index) > 0:
        for k in data_key:
            good_data[k] = data_dict[k][good_index]

    return(good_data)

def extract_gitm_time_arrays(gtime, data_key, ut_index=-1,
                             lon_index=-1, lat_index=-1, alt_index=-1,
                             *args, **kwargs):
    '''
    Extract 1D, 2D, or 3D numpy arrays from a GitmTime structure for specified
    UT, lon, lat, and/or altitudes.

    Input: gtime     = GitmTime structure
           data_key  = List of key to extract
           ut_index  = UT Index to hold constant (or -1, default)
           lon_index = Longitude index to hold constant (or -1, default)
           lat_index = Latitude index to hold constant (or -1, default)
           alt_index = Altitude index to hold constant (or -1, default)

    Output: data_list = list of numpy arrays
    '''
    rout_name = "extract_gitm_time_arrays"
    data_list = list()

    # Set up the ut, lon, lat, and alt index limits
    dim = list(gtime['Longitude'].shape)

    utmin = 0
    utmax = dim[0]
    if ut_index >= 0:
        utmin = ut_index
        utmax = ut_index + 1
        dim[0] = 1

    lonmin = 0
    lonmax = dim[1]
    if lon_index >= 0:
        lonmin = lon_index
        lonmax = lon_index + 1
        dim[1] = 1

    latmin = 0
    latmax = dim[2]
    if lat_index >= 0:
        latmin = lat_index
        latmax = lat_index + 1
        dim[2] = 1

    altmin = 0
    altmax = dim[3]
    if alt_index >= 0:
        altmin = alt_index
        altmax = alt_index + 1
        dim[3] = 1

    # Extract the data, removing any instances where any data key is np.nan
    good_index = [i for i,l in enumerate(gtime[data_key[0]][utmin:utmax,lonmin:lonmax,latmin:latmax,altmin:altmax].flatten()) if not np.isnan(l)]

    def get_unflattened_index(i):
        # Convert from a flattened index to an unflattened index
        if ut_index >= 0:
            iut = ut_index
        else:
            iut = int(i / (dim[1] * dim[2] * dim[3]))

        if lon_index >= 0:
            ilon = lon_index
        else:
            ilon = int((i - iut*dim[1]*dim[2]*dim[3]) / (dim[2]*dim[3]))

        if lat_index >= 0:
            ilat = lat_index
        else:
            ilat = int((i - (ilon + iut*dim[1])*dim[2]*dim[3]) / dim[3])

        if alt_index >= 0:
            ialt = alt_index
        else:
            ialt = int(i - ((iut*dim[1] + ilon)*dim[2] + ilat)*dim[3])

        return(iut, ilon, ilat, ialt)
    # End get_unflattened_index

    if len(data_key) > 1:
        bad_index = list()
        dkey = data_key[1:len(data_key)]

        for i,igood in enumerate(good_index):
            iut, ilon, ilat, ialt = get_unflattened_index(igood)
            nanlist = [gtime[k][iut,ilon,ilat,ialt] for k in dkey
                       if np.isnan(gtime[k][iut,ilon,ilat,ialt])]
            if len(nanlist) > 0:
                bad_index.append(i)

        if len(bad_index) > 0:
            good_index = list(np.delete(np.array(good_index), bad_index))

    # Flatten the data arrays
    data_list = [list() for k in data_key]
    for igood in good_index:
        iut, ilon, ilat, ialt = get_unflattened_index(igood)
        for i,k in enumerate(data_key):
            data_list[i].append(gtime[k][iut,ilon,ilat,ialt])

    # Recast the data as numpy arrays
    last_len = -1
    for i,dlist in enumerate(data_list):
        if last_len > -1 and last_len != len(dlist):
            print rout_name, "ERROR: unequal lengths of vetted arrays"
            return
        else:
            last_len = len(dlist)

        data_list[i] = np.array(dlist)

    return(data_list)
# END


def plot_net_gitm_comp(plot_type, lon_data, lat_data, obs_data, obs_name,
                       obs_scale, obs_units, diff_data, diff_name, diff_scale,
                       diff_units, gitm_key, gitm_alt, gdata, gitm_name,
                       diff_max=None, zmax=None, zmin=None, title=None,
                       color=True, bcolor='#747679', data_coff=False,
                       diff_coff=True, figname=None, draw=True, latlim1=90,
                       latlim2=-90, linc=6, tlon=90, meq=False, earth=False,
                       map_list=[], faspect=True, term_datetime=None,
                       extra_lines=False, *args, **kwargs):
    '''
    Creates three plots of a specified type, one showing the observations, one
    showing the GITM data, and one showing the difference between the two.

    Input: plot_type  = key to determine plot type (rectangular, polar,
                        nsglobal, or snapshot)
           lon_data   = Numpy array with longitude data for matching model-obs
                        points
           lat_data   = Numpy array with latitude data for matching model-obs
                        points
           obs_data   = Numpy array with observational data for matching
                        model-obs points
           obs_name   = Name portion of the observational data label
           obs_scale  = Scale (linear/exponential) for plotting obs. data
           obs_units  = Unit portion of the observational data label
           diff_data  = Numpy array with differences for matching model-obs
                        points
           gitm_key   = Key for the GITM data
           gitm_alt   = Altitude in km to plot the GITM data at.  For a 2D
                        variable like hmF2 or TEC, use 0.0 km.
           gdata      = GitmBin structure with model observations.
           gitm_name  = Name portion of the GITM data label
           diff_max   = Maximum value for the difference (absolute value),
                        if None, will be determined in script (default=None)
           zmin       = minimum value for z variable (default=None)
           zmax       = maximum value for z variable (default=None)
           title      = Plot title (default=None)
           color      = Color (True, default) or black and white (False)?
           bcolor     = Background color (default=)
           data_coff  = Center the data color scale about zero (False, default)?
           diff_coff  = Center the diff color scale about zero (True, default)?
           figname    = Output figure name with a .png suffix (default=None)
           draw       = Draw to screen? (default=True)
           latlim1    = First latitude limit (degrees North, default=90).
                        Purpose varies depending on plot type.  For rectangular,
                        this is the northern latitude limit.  For polar, this
                        is the latitude at the center of the dial.  For
                        snapshot, this is the lower boundary of polar dials.
                        It is not used for nsglobal.
           latlim2    = Second latitude limit (degrees North, default=-90).
                        Purpose varies depending on plot type.  For rectangular,
                        this is the southern latitude limit.  For polar, this
                        is the latitude at the edge of the dial.  This option is
                        not used with the snapshot or nsglobal option.
           linc       = Number of latitude tick incriments (default=6)
           tlon       = Longitude on top of the polar dial (degrees East,
                        default=90)
           meq            = Add a line for the geomagnetic equator?
                            (default=False)
           earth         = Include continent outlines for Earth (default=False)
           map_list      = List of map handles for the specified plot_type
                           (default=empty list)
           faspect       = Keep a true aspect ratio for maps? (default=True)
           term_datetime = Include the solar terminator by shading the night
                           time regions?  If so, include a datetime object
                           with the UT for this map.  Only used if earth=True.
           extra_lines   = Plot a specified lines (good for showing regional
                           boundaries) (default=False).  Provide a list of lists
                           which have the format shown:
                           [x np.array, y np.array, style string (eg 'k-')]
                           where x is in degrees longitude and y is in 
                           degrees latitude

    Output: f = handle to figure
    '''
    rout_name = "plot_net_gitm_comp"

    # Get the desired color bars
    data_color = gpr.choose_contour_map(color, data_coff)
    diff_color = gpr.choose_contour_map(color, diff_coff)

    # Get the altitude index
    ialt = 0
    if gitm_alt > 0.0:
        ialt = gpr.find_alt_index(gdata, 0, 0, alt, units="km")

    # Initialize the z variables, if desired.  GITM and Observational data
    # should share the same scale.
    if(zmin is None):
        obsmin = np.nanmin(obs_data)
        gitmin = np.nanmin(gdata[gitm_key][:,:,ialt])
        zmin = min(obsmin,gitmin)
    if(zmax is None):
        obsmax = np.nanmax(obs_data)
        gitmax = np.nanmax(gdata[gitm_key][:,:,ialt])
        zmax = max(obsmax, gitmax)

    zran = round((zmax-zmin)/6.0)
    if(zran != 0.0):
        zmin = math.floor(float("{:.14f}".format(zmin / zran))) * zran
        zmax = math.ceil(float("{:.14f}".format(zmax / zran))) * zran

    # Set the difference max/min limits, if desired
    if diff_max is None:
        diff_max = max(np.nanmax(diff_data), abs(np.nanmin(diff_data)))

    diff_min = -1.0 * diff_max

    # Initialize the figure, setting the height for a 3 subfigure stack
    fwidth = 6
    fheight = 12
    if(plot_type.find("global") > 0):
        fwidth *= 1.5
    if(plot_type.find("shot") > 0):
        fwidth *= 1.5
        fheight *= 1.5

    f = plt.figure(figsize=(fwidth,fheight))

    # Plot the three datasets using the desired format
    if plot_type.find("shot") > 0:
        if len(map_list) == 3:
            ml = map_list[0]
            mn = map_list[1]
            ms = map_list[2]
        else:
            ml = None
            mn = None
            ms = None
        # Output the observations as a scatter plot
        axl,ml,axn,mn,axs,ms = p3g.plot_snapshot_subfigure(f, 3, 0, lat_data,
                                                           lon_data, obs_data,
                                                           obs_name, obs_scale,
                                                           obs_units, zmax,
                                                           zmin, data_color,
                                                           tlon=tlon,
                                                           blat=latlim1,
                                                           xl=False, yl=False,
                                                           xt=False,
                                                           bcolor=bcolor,
                                                           meq=meq, earth=earth,
                                                           ml=ml, mn=mn, ms=ms,
                                                           faspect=faspect,
                                                           term_datetime=term_datetime)
        # Output the gitm data as a contour after ensuring that the GITM array
        # isn't padded to include unrealistic latitudes
        (i, imin) = gpr.find_lon_lat_index(gdata, 0.0, -90.0, "degrees")
        (i, imax) = gpr.find_lon_lat_index(gdata, 0.0, 90.0, "degrees")
        imax += 1 
        p3g.plot_snapshot_subfigure(f, 3, 1,
                                    np.array(gdata['dLat'][:,imin:imax,ialt]),
                                    np.array(gdata['dLon'][:,imin:imax,ialt]),
                                    np.array(gdata[gitm_key][:,imin:imax,ialt]),
                                    gitm_name, gdata[gitm_key].attrs["scale"],
                                    gdata[gitm_key].attrs["units"], zmax, zmin,
                                    data_color, cb=True, cloc="r", tlon=tlon,
                                    blat=latlim1, title=False, xl=False,
                                    xt=False, bcolor=bcolor, meq=meq,
                                    earth=earth, ml=ml, mn=mn, ms=ms,
                                    faspect=faspect, data_type="contour",
                                    term_datetime=term_datetime)
        # Output the differences as a scatter plot
        p3g.plot_snapshot_subfigure(f, 3, 2, lat_data, lon_data, diff_data,
                                    diff_name, diff_scale, diff_units, diff_max,
                                    diff_min, diff_color, tlon=tlon,
                                    blat=latlim1, title=False, yl=False,
                                    bcolor=bcolor, meq=meq, earth=earth, ml=ml,
                                    mn=mn, ms=ms, faspect=faspect,
                                    term_datetime=term_datetime)
        map_list = list([ml, mn, ms])
    elif plot_type.find("nsglobal") >= 0:
        if len(map_list) == 2:
            mn = map_list[0]
            ms = map_list[1]
        else:
            mn = None
            ms = None

        # Check for boundary lines to plot
        eline_north = False
        eline_south = False
        if type(extra_lines) is list:
            if len(extra_lines) >= 1:
                eline_north = extra_lines[0]

                if len(extra_lines) >= 2:
                    eline_south = extra_lines[1]
                else:
                    print "Only one boundary provided, plotting in north"
            else:
                print "No boundaries provided, better to declare as False"

        # Output the observations as a scatter plot
        axn1,mn,axs1,ms = p3g.plot_nsglobal_subfigure(f, 3, 0, lat_data,
                                                      lon_data, obs_data,
                                                      obs_name, obs_scale,
                                                      obs_units, zmax, zmin,
                                                      data_color, title=True, cb=True,
                                                      elat=latlim1, tlon=tlon, rl=False,
                                                      tl=False, bcolor=bcolor,
                                                      earth=earth, mn=mn, ms=ms,
                                                      faspect=faspect, term_datetime=term_datetime,
                                                      extra_line_n=eline_north,
                                                      extra_line_s=eline_south)
        # Output the gitm data as a contour after ensuring that the GITM array
        # isn't padded to include unrealistic latitudes 
        (i, imin) = gpr.find_lon_lat_index(gdata, 0.0, -90.0, "degrees")
        (i, imax) = gpr.find_lon_lat_index(gdata, 0.0, 90.0, "degrees")
        imax += 1 

        axn2,mn,axs2,ms = p3g.plot_nsglobal_subfigure(f, 3, 1, np.array(gdata['dLat'][:,imin:imax,ialt]), np.array(gdata['dLon'][:,imin:imax,ialt]), np.array(gdata[gitm_key][:,imin:imax,ialt]), gitm_name, gdata[gitm_key].attrs["scale"], gdata[gitm_key].attrs["units"], zmax, zmin, data_color, title=False, cb=True, elat=latlim1, tlon=tlon, tl=False, bcolor=bcolor, earth=earth, mn=mn, ms=ms, data_type="contour", term_datetime=term_datetime, extra_line_n=eline_north, extra_line_s=eline_south)

        # Output the differences as a scatter plot
        p3g.plot_nsglobal_subfigure(f, 3, 2, lat_data, lon_data, diff_data,
                                    diff_name, diff_scale, diff_units, diff_max,
                                    diff_min, diff_color, title=False, cb=True,
                                    elat=latlim1, tlon=tlon, rl=False, bcolor=bcolor,
                                    earth=earth, mn=mn, ms=ms, faspect=faspect,
                                    term_datetime=term_datetime,
                                    extra_line_n=eline_north,
                                    extra_line_s=eline_south)
        map_list = list([mn, ms])
    elif plot_type.find("rect") >= 0:
        if len(map_list) == 1:
            m = map_list[0]
        else:
            m = None
        # Output the observations as a scatter plot
        ax = f.add_subplot(3,1,1)
        con1, m = p3g.plot_rectangular_3D_global(ax, lat_data, lon_data,
                                                 obs_data, obs_name, obs_scale,
                                                 obs_units, zmin, zmax,
                                                 data_color, nlat=latlim1,
                                                 slat=latlim2, linc=linc,
                                                 cloc="r", xl=False, xt=False,
                                                 yl=False, meq=meq,
                                                 bcolor=bcolor, earth=earth,
                                                 m=m, faspect=faspect,
                                                 term_datetime=term_datetime)
        # Output the gitm data as a contour
        ax = f.add_subplot(3,1,2)
        con2, m = p3g.plot_rectangular_3D_global(ax, np.array(gdata['dLat'][:,:,ialt]), np.array(gdata['dLon'][:,:,ialt]), np.array(gdata[gitm_key][:,:,ialt]),
                                                 gitm_name,
                                                 gdata[gitm_key].attrs["scale"],
                                                 gdata[gitm_key].attrs["units"],
                                                 zmin, zmax, data_color,
                                                 nlat=latlim1, slat=latlim2,
                                                 linc=linc, cb=True, cloc="r",
                                                 xl=False, xt=False,
                                                 bcolor=bcolor, meq=meq,
                                                 earth=earth, m=m,
                                                 faspect=faspect,
                                                 data_type="contour",
                                                 term_datetime=term_datetime)
        # Adjust plot dimensions if necessary
        if not earth:
            con1_dim = list(con1.axes.get_position().bounds)
            con2_dim = list(con2.ax.get_position().bounds)
            con2_dim[2] = con1_dim[2]
            con2.ax.set_position(con2_dim)

        # Output the differences as a scatter plot
        ax = f.add_subplot(3,1,3)
        p3g.plot_rectangular_3D_global(ax, lat_data, lon_data, diff_data,
                                       diff_name, diff_scale, diff_units,
                                       diff_min, diff_max, diff_color,
                                       nlat=latlim1, slat=latlim2, linc=linc,
                                       cloc="r", yl=False, bcolor=bcolor,
                                       meq=meq, earth=earth, m=m,
                                       faspect=faspect,
                                       term_datetime=term_datetime)
        map_list = list([m])
    elif plot_type.find("polar") >= 0:
        if len(map_list) == 1:
            m = map_list[0]
        else:
            m = None

        pf = True
        if earth:
            pf = False
        # Output the observations as a scatter plot
        ax = f.add_subplot(3,1,1, polar=pf)
        con1,m = p3g.plot_polar_3D_global(ax, 3, lat_data, lon_data, obs_data,
                                          obs_name, obs_scale, obs_units, zmin,
                                          zmax, data_color, center_lat=latlim1,
                                          edge_lat=latlim2, linc=linc,
                                          top_lon=tlon, cloc="r", tl=False,
                                          rl=False, bcolor=bcolor, earth=earth,
                                          m=m, faspect=faspect,
                                          term_datetime=term_datetime)
        # Output the gitm data as a contour after ensuring that the GITM
        # array isn't padded to include unrealistic latitudes
        ax = f.add_subplot(3,1,2, polar=pf)
        (i, imin) = gpr.find_lon_lat_index(gdata, 0.0, -90.0, "degrees")
        (i, imax) = gpr.find_lon_lat_index(gdata, 0.0, 90.0, "degrees")
        imax += 1 
        con2,m = p3g.plot_polar_3D_global(ax, 3, np.array(gdata['dLat'][:,imin:imax,ialt]), np.array(gdata['dLon'][:,imin:imax,ialt]), np.array(gdata[gitm_key][:,imin:imax,ialt]),
                                          gitm_name,
                                          gdata[gitm_key].attrs["scale"],
                                          gdata[gitm_key].attrs["units"], zmin,
                                          zmax, data_color, center_lat=latlim1,
                                          edge_lat=latlim2, linc=linc,
                                          top_lon=tlon, cb=True, cloc="r",
                                          tl=False, bcolor=bcolor, earth=earth,
                                          m=m, faspect=faspect,
                                          data_type="contour",
                                          term_datetime=term_datetime)

        con1_dim = list(con1.axes.get_position().bounds)
        con2_dim = list(con2.ax.get_position().bounds)
        con2_dim[0] = con2_dim[0] - 0.05
        con2_dim[2] = con1_dim[2]
        con2.ax.set_position(con2_dim)

        # Output the differences as a scatter plot
        ax = f.add_subplot(3,1,3, polar=pf)
        p3g.plot_polar_3D_global(ax, 3, lat_data, lon_data, diff_data,
                                 diff_name, diff_scale, diff_units, diff_min,
                                 diff_max, diff_color, center_lat=latlim1,
                                 edge_lat=latlim2, linc=linc, top_lon=tlon,
                                 cloc="r", rl=False, bcolor=bcolor, earth=earth,
                                 m=m, faspect=faspect,
                                 term_datetime=term_datetime)
        map_list = list([m])
    else:
        print rout_name, "ERROR: uknown plot type [", plot_type, "]"
        return

    if title:
        f.suptitle(title, size="medium")

    # Adjust subplot locations
    if plot_type.find("rect") >= 0 or plot_type.find("polar") >= 0:
        plt.subplots_adjust(left=.15)

    # Draw to screen if desired
    if draw:
        if plt.isinteractive():
            plt.draw() #In interactive mode, you just "draw".
        else:
            # W/o interactive mode, "show" stops the user from typing more 
            # at the terminal until plots are drawn.
            plt.show()

    # Save output file
    if figname is not None:
        plt.savefig(figname)

    return(f, map_list)

# END plot_net_gitm_comp

def plot_sat_gitm_comp(sat_datetime, sat_data, gtrack, gkey, skey, tkey, bkey,
                       title, tmin=None, tmax=None, ymin=None, ymax=None,
                       dmin=None, dmax=None, pmin=None, pmax=None, figname=None,
                       draw=True, fsize=14, sc='#0039A6', gc='#ffcc00',
                       *args, **kwargs):
    '''
    A routine to plot satellite and GITM data to show how a single physical
    quantity varies over the orbit.  Four panels are included; the top panel
    shows the raw satellite data and the GITM data along the track.  The second
    panel shows the matched GITM/satellite data.  The third panel shows the
    difference between the satellite and GITM data.  The fourth panel shows the
    percent difference 100*(sat-GITM)/sat.

    Input:
    sat_datetime = Satellite datetime numpy array: dim(nsat,)
    sat_data     = Satelite data numpy array: dim(nsat,)
    gtrack       = GitmTime structure with GITM data along satellite
                   track and the matching satellite measurements
    gkey         = GITM data key
    skey         = Matched satellite data key
    tkey         = Matched satellite top errorbar key (or None)
    bkey         = Matched satellite bottom errorbar key (or None)
    title        = Plot title

    tmin = UT Minimum (default None)
    tmax = UT Maximum (default None)
    ymin = Dependent variable minimum (default None)
    ymax = Dependent variable maximum (default None)
    dmin = Difference minimum (default None)
    dmax = Difference maximum (default None)
    pmin = Percent difference minimum (default None)
    pmax = Percent difference maximum (default None)

    figname = Output file name (must be .png, default None)
    draw    = Draw to screen? (default is True)

    fsize = Font size (defualt=14)
    sc    = satellite color (default Michigan Blue)
    gc    = GITM color (default Michigan Maize)
    '''
    # Change fontsize, marker size, and line width
    mpl.rc('xtick', labelsize=fsize)
    mpl.rc('ytick', labelsize=fsize)
    mpl.rc('font', size=fsize)
    ms = 8
    lw = 3
    
    # Initialize the figure
    f = plt.figure(figsize=(12,12))
    if not tmin:
        tmin = gtrack['time'][0]
    if not tmax:
        tmax = gtrack['time'][-1]

    # Create the first subplot, showing the raw satellite data and model data
    ax1 = f.add_subplot(411)
    con1 = ax1.plot_date(sat_datetime, sat_data, color=sc,fmt="-")
    con1 = ax1.plot_date(gtrack['time'],gtrack[gkey][:,0,0,0], fmt="o",
                         color=gc, markersize=ms)
    plt.xlim(tmin, tmax)
    if type(ymin) is float and type(ymax) is float:
        plt.ylim(ymin, ymax)

    # Create the second subplot, showing the matched satellite and model data
    ax2 = f.add_subplot(412)
    con2 = ax2.plot_date(gtrack['time'], gtrack[gkey][:,0,0,0], color=gc,
                         fmt="o", markersize=ms)
    if gtrack.has_key(tkey) and gtrack.has_key(bkey):
        con2 = ax2.errorbar(gtrack['time'], gtrack[skey][:,0,0,0],
                            yerr=[gtrack[bkey][:,0,0,0], gtrack[tkey][:,0,0,0]],
                            color=sc, fmt="+", markersize=ms, linewidth=lw)
    else:
        con2 = ax2.plot_date(gtrack['time'], gtrack[skey][:,0,0,0], color=sc,
                             marker="+", markersize=ms, linewidth=lw)
    plt.xlim(tmin,tmax)
    if type(ymin) is float and type(ymax) is float:
        plt.ylim(ymin,ymax)

    # Create the third subplot, showing the difference between the satellite
    # and model data
    ax3 = f.add_subplot(413)
    sdiff = gtrack[skey][:,0,0,0] - gtrack[gkey][:,0,0,0]
    szero = np.zeros(shape=gtrack['time'].shape)

    if gtrack.has_key(tkey) and gtrack.has_key(bkey):
        con3 = ax3.errorbar(gtrack['time'], sdiff,
                            yerr=[gtrack[bkey][:,0,0,0],gtrack[tkey][:,0,0,0]],
                            color=sc, fmt="o", markersize=ms, linewidth=lw)
    else:
        con3 = ax3.plot_date(gtrack['time'], sdiff, color=sc, fmt="o",
                             markersize=ms, linewidth=lw)
    con3 = ax3.plot_date(gtrack['time'], szero, "k-")
    plt.xlim(tmin,tmax)
    if type(dmin) is float and type(dmax) is float:
        plt.ylim(dmin,dmax)

    # Create the fourth subplot, showing the difference between the satellite
    # and model data
    ax4 = f.add_subplot(414)
    sperc = 100.0 * sdiff / gtrack[skey][:,0,0,0]

    if gtrack.has_key(tkey) and gtrack.has_key(bkey):
        bperc = 100.0 * gtrack[bkey][:,0,0,0] / gtrack[skey][:,0,0,0]
        tperc = 100.0 * gtrack[tkey][:,0,0,0] / gtrack[skey][:,0,0,0]
        con4 = ax4.errorbar(gtrack['time'], sperc, yerr=[bperc,tperc],
                            color=sc, fmt="o", markersize=ms, linewidth=lw)
    else:
        con4 = ax4.plot_date(gtrack['time'], sperc, color=sc, fmt="o",
                             markersize=ms, linewidth=lw)
    con4 = ax4.plot_date(gtrack['time'], szero, "k-")
    plt.xlim(tmin,tmax)
    if type(pmin) is float and type(pmax) is float:
        plt.ylim(pmin,pmax)

    # Set the labels, title, and tick intervals
    f.suptitle(title)
    ylabel1 = "{:s} (${:s}$)".format(gtrack[gkey].attrs['name'],
                                     gtrack[gkey].attrs['units'])
    ax1.set_ylabel(ylabel1)
    ax2.set_ylabel(ylabel1)
    ax3.set_ylabel("Obs$-$GITM {:s}".format(ylabel1))
    ax4.set_ylabel(r"% Diff {:s}".format(ylabel1))

    dtics, dmtics, dfmt = pybats.smart_timeticks([gtrack['time'][0],
                                                  gtrack['time'][-1]])
    ax1.xaxis.set_major_locator(dtics)
    ax2.xaxis.set_major_locator(dtics)
    ax3.xaxis.set_major_locator(dtics)
    ax4.xaxis.set_major_locator(dtics)
    ax1.xaxis.set_minor_locator(dmtics)
    ax2.xaxis.set_minor_locator(dmtics)
    ax3.xaxis.set_minor_locator(dmtics)
    ax4.xaxis.set_minor_locator(dmtics)

    xfmt1 = FormatStrFormatter("")
    xfmt2 = FuncFormatter(gtrack.sat_dateloc_ticks)
    ax1.xaxis.set_major_formatter(xfmt1)
    ax2.xaxis.set_major_formatter(xfmt1)
    ax3.xaxis.set_major_formatter(xfmt1)
    ax4.xaxis.set_major_formatter(xfmt2)

    gitm_time.set_sat_dateloc_label(ax4, gtrack.has_key("Magnetic Latitude"))
    plt.subplots_adjust(hspace=0.1, top=.94, bottom=.2, right=.9)

    # Output the plots to the screen and/or file as desired
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
# END plot_sat_gitm_comp
