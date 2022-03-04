#!/usr/bin/env python
""" Standard model visualization routines
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import re
from glob import glob
import argparse
import sys

from aetherpy.io import read_routines
from aetherpy.utils import inputs, time_conversion
from aetherpy.plot import data_prep, movie_routines

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Plot Aether / GITM model results')
    
    parser.add_argument('-list',  \
                        action='store_true', default = False, \
                        help = 'list variables in file')

    parser.add_argument('-timeplot',  \
                        action='store_true', default = False, \
                        help = 'Plot integrated (or mean) value vs. time')

    parser.add_argument('-mean',  \
                        action='store_true', default = False, \
                        help = 'Plot mean value instead of integrated value')
    
    parser.add_argument('-var',  \
                        default = 3, type = int, \
                        help = 'variable to plot (number)')
    parser.add_argument('-cut', metavar = 'cut',  default ='alt', \
                        choices = ['alt', 'lat', 'lon'], 
                        help = 'alt,lat,lon : which cut you would like')
    parser.add_argument('-ext',  default ='png', \
                        choices = ['png', 'jpg', 'pdf'], 
                        help = 'plot type file extention')
    parser.add_argument('-winds', default = False,\
                        help='overplot winds', \
                        action="store_true")
    parser.add_argument('-alt', metavar = 'alt', default =400.0, type = int, \
                        help = 'altitude :  alt in km (closest)')
    parser.add_argument('-lat', metavar = 'lat',  default =-100.0, \
                        help = 'latitude : latitude in degrees (closest)')
    parser.add_argument('-lon', metavar = 'lon',  default =-100.0,\
                        help = 'longitude in degrees (closest)')
    parser.add_argument('-alog',  default = False,
                        action="store_true",
                        help = 'plot the log of the variable')
    parser.add_argument('-IsLog', default =False,
                        help='plot the log of the variable', 
                        action="store_true")    
    parser.add_argument('-diff', default = False, 
                        action = 'store_true',
                        help = 'plot difference of files (2 files needed)')
    parser.add_argument('-mkv',
                        action = 'store_true',
                        default = False,
                        help = 'movie format = mkv')
    parser.add_argument('-mp4',
                        action = 'store_true',
                        default = True,
                        help = 'movie format = mp4')
    parser.add_argument('-gif',
                        action='store_true',
                        default = False, 
                        help = 'movie format = gif')
    parser.add_argument('-movie',  default = False,\
                        action='store_true',
                        help = 'Make a movie out of results')
    parser.add_argument('-tec',  default = False, \
                        action='store_true',
                        help = 'plot total electron content (TEC)')
    parser.add_argument('-rate',  default =30,\
                        help = 'framerate for movie')
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')
    
    args = parser.parse_args()

    return args


## ----------------------------------------------------------------------------
## Define the support routines
#
#def get_help(file_vars=None):
#    """ Provide string explaining how to run the command line interface
#
#    Parameters
#    ----------
#    file_vars : list or NoneType
#        List of file variables or None to exclude this output (default=None)
#
#    Returns
#    -------
#    help_str : str
#        String with formatted help statement
#
#    """
#
#    mname = os.path.join(
#        os.path.commonpath([inputs.__file__, data_prep.__file__]),
#        'run_plot_model_results.py') if __name__ == '__main__' else __name__
#
#    help_str = 'Usage:\n{:s} -[flags] [filenames]\n'.format(mname)
#    help_str += 'Flags:\n'
#    help_str += '       -help : print this message, include filename for '
#    help_str += 'variable names and indices\n'
#    help_str += '       -var=number : index of variable to plot\n'
#    help_str += '       -cut=alt, lat, or lon : which cut you would like\n'
#    help_str += '       -alt=number : alt in km or grid number (closest)\n'
#    help_str += '       -lat=number : latitude in degrees (closest)\n'
#    help_str += '       -lon=number: longitude in degrees (closest)\n'
#    help_str += '       -log : plot the log of the variable\n'
#    help_str += '       -winds : overplot winds\n'
#    help_str += '       -tec : plot the TEC variable\n'
#    help_str += '       -movie=number : provide a positive frame rate to '
#    help_str += 'create a movie\n'
#    help_str += '       -ext=str : figure or movie extension\n'
#    help_str += 'At end, list the files you want to plot. This code should '
#    help_str += 'work with either GITM files (*.bin) or Aether netCDF files '
#    help_str += '(*.nc)'
#
#    if file_vars is not None:
#        help_str += "File Variables (index, name):\n"
#        for ivar, var in enumerate(file_vars):
#            help_str += "               ({:d}, {:s})\n".format(ivar, var)
#
#    return help_str
#
#
#def get_command_line_args(argv):
#    """ Parse the arguements and set to a dictionary
#
#    Parameters
#    ----------
#    argv : list
#        List of arguments fed on the command line
#
#    Returns
#    -------
#    args : dict
#        A dictionary containing information about arguements, including:
#        filelist (list of filenames), gitm (flag that is true for GITM input,
#        determined by examining filelist naming convention),
#        var (variable index to plot), cut (coordinate to hold constant),
#        diff (difference with other plots),
#        movie (framerate for movie, which is > 0 if a movie is desired),
#        ext (output extension), winds (flag to plot with winds),
#        alt (to plot), lat (to plot), lon (to plot),
#        log (flag to use log scale), and help (flag to display help)
#
#    """
#    # Initialize the arguments to their default values
#    args = {'filelist': [], 'log': False, 'var': 15, 'alt': 400, 'tec': False,
#            'lon': np.nan, 'lat': np.nan, 'cut': 'alt', 'winds': False,
#            'diff': False, 'IsGitm': False, 'HasHeader': False, 'movie': 0,
#            'ext': 'png'}
#
#    arg_type = {'filelist': list, 'log': bool, 'var': int, 'alt': int,
#                'tec': bool,
#                'lon': float, 'lat': float, 'cut': str, 'help': bool,
#                'winds': bool, 'diff': bool, 'IsGitm': bool, 'HasHeader': bool,
#                'tec': bool,
#                'movie': int, 'ext': str}
#
#    # If there is input, set default help to False
#    args['help'] = False if len(argv) > 0 else True
#
#    # Cycle through all arguments except the first, saving input
#    for arg in argv:
#        # Treat the file list and formatting seperately
#        if arg.find('-') == 0:
#            # This is not a filename, remove the dash to get the key
#            split_arg = arg.split('=')
#            akey = split_arg[0][1:]
#
#            # Get the argument value as the desired type
#            if akey not in arg_type.keys():
#                raise ValueError(''.join(['unknown command line input, ',
#                                          arg, ', try -help for details']))
#
#            if len(split_arg) == 1:
#                if arg_type[akey] == bool:
#                    arg_val = True
#                else:
#                    raise ValueError('expected equality after flag {:}'.format(
#                        akey))
#            else:
#                if arg_type[akey] == int:
#                    arg_val = int(split_arg[1])
#                elif arg_type[akey] == float:
#                    arg_val = float(split_arg[1])
#                elif arg_type[akey] == str:
#                    arg_val = split_arg[1]
#                else:
#                    # This is boolean input
#                    arg_val = inputs.bool_string(split_arg[1])
#
#            # Assign the output
#            if akey.find('tec') == 0:
#                args['var'] = 34
#            else:
#                args[akey] = arg_val
#        else:
#            # Save the filenames
#            args['filelist'].append(arg)
#
#            m = re.match(r'(.*)bin',arg)
#            if m:
#                args['IsGitm'] = 1
#                args['HasHeader'] = 0
#                # check for a header file:
#                checkFile = glob(m.group(1)+"header")
#                if (len(checkFile) > 0):
#                    if (len(checkFile[0]) > 1):
#                        args['HasHeader'] = 1
#            else:
#                args['IsGitm'] = 0
#
#    # Update default movie extention for POSIX systems
#    if args['movie'] > 0 and args['ext'] == 'png':
#        if (os.name == "posix"):
#            args['ext'] = "mkv"
#        else:
#            args['ext'] = "mp4"
#
#    return args


def determine_file_type(file):

    IsGitm = False
    HasHeader = False
    m = re.match(r'(.*)bin', file)
    if m:
        IsGitm = True
        # check for a header file:
        checkFile = glob(m.group(1)+"header")
        if (len(checkFile) > 0):
            if (len(checkFile[0]) > 1):
                HasHeader = True

    return IsGitm, HasHeader

def fix_vars(vars):
    newvars = []
    for v in vars:
        nv = re.sub('!U', '', v)
        nv = re.sub('!N', '', nv)
        nv = re.sub('!D', '', nv)
        newvars.append(nv)

    return newvars

# ----------------------------------------------------------------------------
# Define the main plotting routine

def plot_model_results():

    # Get the input arguments
    args = get_args()

    # determine what kind of files we are dealing with
    IsGitm, HasHeader = determine_file_type(args.filelist[0])
    
    if ((IsGitm) and (not HasHeader)):
        header = read_routines.read_gitm_headers(args.filelist, finds = 0)
    else:
        if (HasHeader):
            header = read_routines.read_aether_ascii_header(args.filelist)
            IsGitm = 0
        else:
            header = read_routines.read_aether_header(args.filelist)

    header['vars'] = fix_vars(header['vars'])
    
    if (args.list):
        for k, v in header.items():
            if (k != 'vars'):
                print(k, '-> ', v)
            else:
                print('vars : ')
                for i, var in enumerate(v):
                    print(i, var)
        exit()
        
    if (args.var >= len(header["vars"])):
        raise ValueError("requested variable doesn't exist: {:d}>{:d}".format(
            args.var, len(header["vars"])))

    # Define the plotting inputs
    plot_vars = [0, 1, 2, args.var]

    # Update plotting variables to include the wind, if desired
    if args.winds:
        plot_vars.append(16 if args.cut in ['alt', 'lat'] else 17)
        plot_vars.append(18 if args.cut in ['lat', 'lon'] else 17)
        all_winds_x = []
        all_winds_y = []

    # Prepare to load the desired file data
    all_2dim_data = []
    all_times = []
    all_int_data = []

    for j, filename in enumerate(args.filelist):
        # Read in the data file
        if IsGitm:
            data = read_routines.read_gitm_file(filename, plot_vars)
            ivar = args.var
        else:
            if j == 0:
                var_list = []
                for pvar in plot_vars:
                    var_list.append(header["vars"][pvar])
            if (HasHeader):
                data = read_routines.read_aether_one_binary_file(header, j, plot_vars)
                ivar = args.var
            else:
                data = read_routines.read_aether_file(filename, var_list)
                ivar = 3

        # For the first file, initialize the necessary plotting data
        if j == 0:
            # Get 1D arrays for the coordinates
            alts = data[2][0][0] / 1000.0  # Convert from m to km
            lons = np.degrees(data[0][:, 0, 0])  # Convert from rad to deg
            lats = np.degrees(data[1][0, :, 0])  # Convert from rad to deg
            # Find the desired index to cut along to get a 2D slice
            isgrid = False
            if (args.cut == 'alt'):
                pos = args.alt
                if (len(alts) == 1):
                    print("Only one alt found, setting alt pos = 0");
                    pos = 0
                    isgrid = True
                lat2d = data[1][:, :, 0]  # Convert from rad to deg
                dlon = data[0][1, 0, 0] - data[0][0, 0, 0]
                dlat = data[1][0, 1, 0] - data[1][0, 0, 0]
                area = np.cos(lat2d) * dlon * dlat *((6372.0 + 100.0)*1000.0)**2
                int_area = np.sum(area)
            if (args.cut == 'lon'):
                pos = args.lon
            if (args.cut == 'lat'):
                pos = args.lat
                
            icut, cut_data, x_pos, y_pos, z_val = data_prep.get_cut_index(
                lons, lats, alts, pos, isgrid, args.cut)

        if (args.cut == 'alt'):
            int_data = data[ivar][cut_data] * area
            if (args.mean):
                int_data = int_data / int_area
            all_int_data.append(np.sum(int_data))
                
        # Save the time data
        all_times.append(data["time"])

        # Save the z-axis data
        if args.tec:
            all_2dim_data.append(data_prep.calc_tec(alts, data[ivar], 2, -4))
        else:
            all_2dim_data.append(data[ivar][cut_data])

            if (args.winds):
                all_winds_x.append(data[plot_vars[-1]][cut_data])
                all_winds_y.append(data[plot_vars[-1]][cut_data])

    # Convert data list to a numpy array
    all_2dim_data = np.array(all_2dim_data)
    
    if args.winds:
        all_winds_x = np.array(all_winds_x)
        all_winds_y = np.array(all_winds_y)

    # If desired, take the log of the data
    if args.alog:
        all_2dim_data = np.log10(all_2dim_data)

    # Define plotting limits
    symmetric = False
    cmap = mpl.cm.plasma
    
    maxi = all_2dim_data.max() * 1.01
    mini = all_2dim_data.min() * 0.99

    factorString = ''
    if ((mini < 0.0) and (not args.alog)):
        symmetric = True
        cmap = mpl.cm.bwr
        maxi = abs(all_2dim_data).max() * 1.05
        mini = -maxi
    else:
        if (not args.alog):
            if ((np.log10(maxi) > 5.0) or (np.log10(maxi) < -5.0)):
                factor = 10**float(int(np.log10(maxi)))
                all_2dim_data = all_2dim_data / factor
                maxi = maxi / factor
                mini = mini / factor
                factorString = '(x%7.1e)' % factor

    if args.cut == 'alt':

        mask_north = ((y_pos >= 40) & (y_pos <= 90.0))
        mask_south = ((y_pos <= -40) & (y_pos >= -90.0))
        plot_north = mask_north.max()
        plot_south = mask_south.max()

        if plot_north:
            if symmetric:
                maxi_north = abs(all_2dim_data[:, :, mask_north]).max() * 1.05
                mini_north = -maxi_north
            else:
                maxi_north = all_2dim_data[:, :, mask_north].max() * 1.05
                mini_north = all_2dim_data[:, :, mask_north].min() * 0.95

        if plot_south:
            if symmetric:
                maxi_south = abs(all_2dim_data[:, :, mask_south]).max() * 1.05
                mini_south = -maxi_south
            else:
                maxi_south = all_2dim_data[:, :, mask_south].max() * 1.05
                mini_south = all_2dim_data[:, :, mask_south].min() * 0.95

    # Define plot range
    minx = (x_pos[1] + x_pos[2]) / 2.0
    maxx = (x_pos[-2] + x_pos[-3]) / 2.0
    miny = (y_pos[1] + y_pos[2]) / 2.0
    maxy = (y_pos[-2] + y_pos[-3]) / 2.0

    # Prepare the output filename
    filename = "var{:02d}_{:s}{:03d}".format(args.var, args.cut, icut)

    if args.movie > 0:
        img_file_fmt = movie_routines.setup_movie_dir(filename)
    else:
        img_file_fmt = filename+'_{:}.'+args.ext

    if (args.timeplot):
        fig = plt.figure(figsize=(10, 8.5))
        ax = fig.add_subplot(111)
        ax.plot(all_times, all_int_data)

        start = all_times[0].strftime('%b %d, %Y %H:%M')
        end = all_times[-1].strftime('%b %d, %Y %H:%M')
        ax.set_xlabel(start + ' to ' + end)
        if (args.mean):
            type = 'mean'
        else:
            type = 'integral'
        ax.set_ylabel('Global '+type+' (' + header["vars"][args.var] + ')')

        title = 'Global '+type+' of ' + header["vars"][args.var]
        title = title + ' at {:.2f} km'.format(z_val)
        ax.set_title(title)
        
        stime = all_times[0].strftime('%y%m%d')
        fig.savefig(filename+'_'+stime+'.'+args.ext)
        exit()

    # Create a plot for each time
    for itime, utime in enumerate(all_times):
        # Initialize the figure
        fig = plt.figure(constrained_layout=False, figsize=(10, 8.5))

        gs1 = mpl.gridspec.GridSpec(nrows=2, ncols=2, wspace=0.0, hspace=0)
        gs = mpl.gridspec.GridSpec(nrows=2, ncols=2, wspace=-0.05,
                                   left=0.02, right=0.95,
                                   top = 0.99, bottom = 0.05)
        #ax = fig.add_subplot(gs1[1, 0:2])
        ax = fig.add_axes([0.07, 0.06, 0.98, 0.48])

        # Plot the global data set (square plot at bottom if three plots):

        dx = (x_pos[1] - x_pos[0])/2.0
        xp = np.append(x_pos - dx, x_pos[-1:]+dx)
        dy = (y_pos[1] - y_pos[0])/2.0
        yp = np.append(y_pos - dy, y_pos[-1]+dy)
        con = ax.pcolormesh(xp, yp, all_2dim_data[itime].transpose(),
                        vmin=mini, vmax=maxi, cmap=cmap, shading='auto')

        # Add the winds, if desired
        if args.winds:
            ax.quiver(x_pos, y_pos, all_winds_x[itime].transpose(),
                      all_winds_y[itime].transpose())
        ax.set_ylim([miny, maxy])
        ax.set_xlim([minx, maxx])

        # Set the labels and aspect ratio
        ax.set_title("{:s}; {:s}: {:.2f} {:s}".format(
            utime.strftime("%d %b %Y %H:%M:%S UT"), args.cut, z_val,
            'km' if args.cut == 'alt' else r'$^\circ$'))
        ax.set_xlabel(r'Latitude ($^\circ$)' if args.cut == 'lon'
                      else r'Longitude ($^\circ$)')
        ax.set_ylabel(r'Latitude ($^\circ$)' if args.cut == 'alt'
                      else r'Altitude (km)')
        if args.cut == 'alt':
            ax.set_aspect(1.0)

        # Set the colorbar
        cbar = fig.colorbar(con, ax=ax, shrink=0.75, pad=0.02)
        cbar.set_label(header["vars"][args.var] + factorString, rotation=90)

        # If this is an altitude slice, add polar dials
        if args.cut == 'alt' and (plot_north or plot_south):
            # Set the common inputs
            shift = time_conversion.calc_time_shift(utime)

            #xlabels = ['12', '18', '00']
            #xlabelpos = [np.pi/2, np.pi, 3*np.pi/2]
            xlabels = []
            xlabelpos = []
            ylabels = [r'80$^\circ$', r'70$^\circ$', r'60$^\circ$',
                       r'50$^\circ$']

            ylabelpos = [10.0, 20.0, 30.0, 40.0]
            xticks = np.arange(0, 2 * np.pi, np.pi / 2.0)
            yticks = np.arange(10, 50, 10)

            if plot_north:
                # Top Left Graph Northern Hemisphere
                #ax2 = fig.add_subplot(gs[0, 0], projection='polar')
                ax2 = fig.add_axes([0.06, 0.55, 0.43, 0.43], projection='polar')
                yp = 90.0 - y_pos[mask_north]
                dy = (int(100.0*(yp[1]-yp[0]))/100.0)/2.0
                yp = np.append(yp - dy, yp[-1] + dy)
                xp = np.radians(x_pos + shift - 90.0)
                dx = (xp[1] - xp[0])/2
                xp = np.append(xp - dx, xp[-1] + dx)
                z = all_2dim_data[itime][:, mask_north].transpose()
                ax2.grid(False)
                conn = ax2.pcolormesh(xp, yp,
                                      z,
                                      shading = 'auto',
                                      vmin=mini_north, vmax=maxi_north,
                                      cmap=cmap)
                ax2.set_xticks(xlabelpos)
                ax2.set_xticklabels(xlabels)
                ax2.text(-np.pi/2, 45.0, '00 LT',
                         verticalalignment='top',
                         horizontalalignment='center')
                ax2.text(np.pi/2, 45.0, '12 LT',
                         verticalalignment='bottom',
                         horizontalalignment='center')
                ax2.text(-np.pi, 47.0, '18 LT',
                         verticalalignment='center',
                         horizontalalignment='center',
                         rotation = 90)
                ax2.text(3*np.pi/4, 45.0, 'North',
                         verticalalignment='bottom',
                         horizontalalignment='center',
                         rotation = 45)
                ax2.set_yticks(ylabelpos)
                ax2.set_yticklabels(ylabels)
                ax2.grid(linestyle=':', color='black')
                ax2.set_xticks(xticks)
                ax2.set_yticks(yticks)
                ax2.set_ylim([0, 45])
                cbar2 = fig.colorbar(conn, ax=ax2, shrink=0.5, pad=0.01)
                cbar2.set_label(header["vars"][args.var] + factorString, rotation=90)                

            if plot_south:
                # Top Right Graph Southern Hemisphere
                rad, theta = np.meshgrid(90.0 + y_pos[mask_south],
                                         np.radians(x_pos + shift - 90.0))
                #ax3 = fig.add_subplot(gs[0, 1], projection='polar')
                ax3 = fig.add_axes([0.54, 0.55, 0.43, 0.43], projection='polar')

                yp = 90.0 + y_pos[mask_south]
                dy = (int(100.0*(yp[1]-yp[0]))/100.0)/2.0
                yp = np.append(yp - dy, yp[-1] + dy)
                xp = np.radians(x_pos + shift - 90.0)
                dx = (xp[1]-xp[0])/2.0
                xp = np.append(xp - dx, xp[-1] + dx)
                z = all_2dim_data[itime][:, mask_south].transpose()
                ax3.grid(False)
                cons = ax3.pcolormesh(xp, yp, z,
                                  shading = 'auto',
                                  vmin=mini_south, vmax=maxi_south, cmap=cmap)
                ax3.set_xticks(xlabelpos)
                ax3.set_xticklabels(xlabels)
                ax3.text(-np.pi/2, 45.0, '00 LT',
                         verticalalignment='top',
                         horizontalalignment='center')
                ax3.text(np.pi/2, 45.0, '12 LT',
                         verticalalignment='bottom',
                         horizontalalignment='center')
                ax3.text(-np.pi, 47.0, '18 LT',
                         verticalalignment='center',
                         horizontalalignment='center',
                         rotation = 90)
                ax3.text(3*np.pi/4, 45.0, 'South',
                         verticalalignment='bottom',
                         horizontalalignment='center',
                         rotation = 45)
                ax3.set_yticks(ylabelpos)
                ax3.set_yticklabels(ylabels)
                ax3.grid(linestyle=':', color='black')
                ax3.set_xticks(xticks)
                ax3.set_yticks(yticks)
                ax3.set_ylim([0, 45])
                cbar3 = fig.colorbar(cons, ax=ax3, shrink=0.5, pad=0.01)
                cbar3.set_label(header["vars"][args.var] + factorString, rotation=90)                

        # Format the output filename
        if args.movie > 0:
            fmt_input = itime
        else:
            fmt_input = utime.strftime('%y%m%d_%H%M%S')
        outfile = img_file_fmt.format(fmt_input)

        # Save the output file
        print("Writing file : ", outfile)
        fig.savefig(outfile)
        plt.close(fig)

    # Create a movie, if desired
    if args.movie > 0:
        movie_routines.save_movie(filename, ext=args.ext,
                                  rate=args.rate)
        
    return


# Needed to run main script as the default executable from the command line
if __name__ == '__main__':
    plot_model_results()
