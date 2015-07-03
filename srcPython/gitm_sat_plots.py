#!/usr/bin/env python
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

'''
Open a GITM list of satellite files and a satellite path file, print an output
 file with this data, and create plots of the atmospheric variables as a
 function of local time.
'''

# Global Imports
import string
import numpy as np
import datetime as dt
from datetime import timedelta
from spacepy import pybats
import gitm
import gitm_time
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, FuncFormatter

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

# Interesting plot, needs to be cleaned and inputs fixed
#def linear_eventplot(self, ykey, fstart, fend, altmin, altmax, color):
#    '''
#    Make a linear plot where the x axis it time in datetime format
#    '''
#    from spacepy import pybats
#    import matplotlib
#    from matplotlib import pyplot
#    from matplotlib import dates
#    import math
    # Open the figure
#    matplotlib.pyplot.figure(num=1, figsize=(12,18))
#    matplotlib.pyplot.subplots_adjust(hspace=.075)
#    matplotlib.pyplot.suptitle('Climatology and Event GITM Runs along Satellite Track')
    # Open the first subfigure
#    matplotlib.pyplot.subplot(212)
    # Prepare the date format using pybats smart date system
#    sdates = matplotlib.dates.date2num(self['time'][0][0])
#    ax     = matplotlib.pyplot.gca()
#    xfmt   = matplotlib.dates.DateFormatter("%H:%M\n%j/%Y")
#    xtics  = matplotlib.dates.MinuteLocator(interval=30)
#    xmtics = matplotlib.dates.MinuteLocator(interval=5)
#    xmin   = fstart + timedelta(hours=-1)
#    xmax   = fend + timedelta(hours=3)    
    # Find the index range that match the desired time range
#    imin = self['time'][0][0].searchsorted(xmin)
#    imax = self['time'][0][0].searchsorted(xmax)
    # Set the date axis
#    ax.xaxis.set_major_locator(xtics)
#    ax.xaxis.set_minor_locator(xmtics)
#    ax.xaxis.set_major_formatter(xfmt)
#    matplotlib.pyplot.xlim(xmin, xmax)
#    matplotlib.pyplot.xlabel("Universal Time")
    # Set the y axis
#    ylabel = self[ykey].attrs['label'] if self[ykey].attrs.has_key('label') else ykey
#    if self[ykey].attrs.has_key('units'):
#        ylabel = string.join([ylabel, " (", self[ykey].attrs['units'], ")"], "")
#    matplotlib.pyplot.ylabel(ylabel)
#    ecolor = color
#    color  = string.join([color, "o:"], "")
    # Use a nonlinear scale, if indicated
#    if self[ykey].attrs.has_key('scale'):
#        ax.set_yscale(self[ykey].attrs['scale'])
    # Plot the data in the first subplot
#    matplotlib.pyplot.plot_date(sdates[imin:imax], self[ykey][0][0][imin:imax], fmt=color, xdate=True, ydate=False, linewidth=.5, markeredgecolor=ecolor,markerfacecolor='none')
    # Add vertical lines flanking the event
#    ymin, ymax = matplotlib.pyplot.ylim()
#    matplotlib.pyplot.vlines(fstart,ymin,ymax,color='k',linestyles='solid')
#    matplotlib.pyplot.vlines(fend,ymin,ymax,color='k',linestyles='solid')
    # Open the second subfigure
#    matplotlib.pyplot.subplot(412)
    # Format the x axis to match the previous subfigure
#    ax   = matplotlib.pyplot.gca()
#    xfmt = matplotlib.dates.DateFormatter("")
#    ax.xaxis.set_major_locator(xtics)
#    ax.xaxis.set_minor_locator(xmtics)
#    ax.xaxis.set_major_formatter(xfmt)
#    matplotlib.pyplot.xlim(xmin, xmax)
    # Format the y axis
#    ytics = matplotlib.ticker.MultipleLocator(6)
#    yfmt  = matplotlib.ticker.FormatStrFormatter('%d')
#    ax.yaxis.set_major_locator(ytics)
#    ax.yaxis.set_major_formatter(yfmt)
#    matplotlib.pyplot.ylabel("Local Time (h)")
#    matplotlib.pyplot.ylim(0, 24)
    # Plot the local time data for the second subplot
#    matplotlib.pyplot.plot_date(sdates[imin:imax], self["LT"][0][0][imin:imax], fmt='ko', xdate=True, ydate=False, linewidth=.5, ms=5)
    # Add vertical lines flanking the event
#    matplotlib.pyplot.vlines(fstart, 0, ymax, color='k',linestyles='solid')
#    matplotlib.pyplot.vlines(fend, 0, ymax, color='k',linestyles='solid')
    # Find the index range that match the event time range
#    fmin = self['time'][0][0].searchsorted(fstart)
#    fmax = self['time'][0][0].searchsorted(fend)
    # Add hash lines to event points
#    matplotlib.pyplot.plot_date(sdates[fmin:fmax], self["LT"][0][0][fmin:fmax], fmt='k|', xdate=True, ydate=False, ms=15)
    # Mark the starting and ending points
#    matplotlib.pyplot.plot_date(sdates[imin], self['LT'][0][0][imin], fmt='kx', ms=15)
#    matplotlib.pyplot.plot_date(sdates[imax-1],self['LT'][0][0][imax-1], fmt='ko', ms=15, mfc='none', mec='k')
    # Mark every 15 points
#    matplotlib.pyplot.plot_date(sdates[imin+15:imax:15], self["LT"][0][0][imin+15:imax:15], fmt='k|', xdate=True, ydate=False, ms=15, lw=3)
    # Open the third subfigure
#    matplotlib.pyplot.subplot(411)
    # Format the x axis
#    ax    = matplotlib.pyplot.gca()
#    xtics = matplotlib.ticker.MultipleLocator(60)
#    xfmt  = matplotlib.ticker.FormatStrFormatter('%d')
#    ax.xaxis.set_ticks_position("top")
#    ax.xaxis.set_label_position("top")
#    ax.xaxis.set_major_locator(xtics)
#    ax.xaxis.set_major_formatter(xfmt)
#    matplotlib.pyplot.xlim(0, 360)
#    matplotlib.pyplot.xlabel("Longitude (degrees)")
    # Format the y axis
#    ytics = matplotlib.ticker.MultipleLocator(5)
#    yfmt  = matplotlib.ticker.FormatStrFormatter('%d')
#    ymin  = math.floor(dmarray.min(self['dLat'][0][0])/5)*5
#    ymax  = math.ceil(dmarray.max(self['dLat'][0][0])/5)*5
#    ax.yaxis.set_major_locator(ytics)
#    ax.yaxis.set_major_formatter(yfmt)
#    matplotlib.pyplot.ylabel("Latitude (degrees)")
#    matplotlib.pyplot.ylim(ymin, ymax)
    # Plot the local time data for the second subplot
#    cm = matplotlib.pyplot.get_cmap('YlGnBu_r')
#    sc = matplotlib.pyplot.scatter(self['dLon'][0][0][imin:imax], self['dLat'][0][0][imin:imax], c=self['Altitude'][0][0][imin:imax]/1000, marker='o', edgecolors='none', vmin=altmin, vmax=altmax, cmap=cm)
#    pos   = ax.get_position().bounds
#    bar   = matplotlib.pyplot.colorbar(sc, pad=0.02)
#    bp    = list(bar.ax.get_position().bounds)
#    bp[0] = .91
#    bar.set_label('Altitude (km)')
#    ax.set_position(pos)
#    bar.ax.set_position(bp)
    # Add hash lines to event points
#    matplotlib.pyplot.scatter(self['dLon'][0][0][fmin:fmax], self['dLat'][0][0][fmin:fmax], marker='|', edgecolors='k', s=150)
    # Mark the starting and ending points
#    matplotlib.pyplot.scatter(self['dLon'][0][0][imin], self['dLat'][0][0][imin], marker='x', edgecolors='k', s=100)
#    matplotlib.pyplot.scatter(self['dLon'][0][0][imax-1], self['dLat'][0][0][imax-1], marker='o', edgecolors='k', facecolors='none', s=100)
    # Mark every 15 points
#    matplotlib.pyplot.scatter(self['dLon'][0][0][imin+15:imax:15], self['dLat'][0][0][imin+15:imax:15], marker='|', edgecolor='k', s=150, linewidth=3)
# END linear_eventplot
