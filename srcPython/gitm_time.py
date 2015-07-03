#!/usr/bin/env python
#-----------------------------------------------------------------------------
# $Id: gitm_time.py,v 1.6 2013/12/23 15:46:33 agburr Exp $
#
# GITM_Time.py, Angeline Burrell (AGB), UMich, June 2013
#
# Comments: Defines a class to hold data from multiple GITM binary output files,
#           allowing UT dependence to be explored
#
# AGB 10/17/13: Added routines to assist in plotting satellite data and matching
#               observational in GITM data
#
# Contains: class GitmTime - The class for the GITM binary, which will contain
#                            data from multiple GITM output binaries
#           def load_multiple_gitm_bin - A routine to load multiple GITM binary
#                                        files, returning a list of the
#                                        GitmBin data structures
#           def set_sat_dateloc_label - A routine to add a label to a plot
#                                       using GitmTime satellite ticks
#------------------------------------------------------------------------------

'''
PyBats submodel for handling input/output for the Global Ionosphere-Thermosphere
Model (GITM), providing data structures with universal time dependence.
'''

# Global Imports
import numpy as np
from spacepy.pybats import PbData
from spacepy.datamodel import dmarray
from spacepy.pybats import gitm
import string
from copy import deepcopy as dc

# Temporary
import gitm

class GitmTime(PbData):
    '''
    Object containing GITM data from multiple GTIM binaries, providing a 
    data structure with UT time dependence.  Requires a list of GitmBin
    data structures as input
    '''

    def __init__(self, gitmlist, *args, **kwargs):

        super(GitmTime, self).__init__(*args, **kwargs) # Init as PbData
        self.attrs['nFiles'] = 0

        self._appendgitm(gitmlist)
        
    def __repr__(self):
        return 'File with list of GITM binaries: %s' % (self.attrs['file'])

    def _appendgitm(self, gitmlist):
        '''
        Append GITM binary file into a single data structure.  Requires a
        list of GitmBin data structures as input.  If the GitmTime structure
        already exists, the new data will be appended
        '''

        # Import local packages
        from datetime import datetime

        # Update the number of files in the attributes
        self.attrs['nFiles'] += len(gitmlist)

        #Initialize and fill the new data structure
        for i,gData in enumerate(gitmlist):

            # Determine which keys are new and which already exist
            newkeys = gData.keys()
            oldkeys = list()

            if i > 0:
                oldkeys = self.keys()

                for old in oldkeys:
                    try:
                        n = newkeys.index(old)
                        newkeys.pop(n)
                    except ValueError:
                        if old != 'file' and old != 'magfile':
                            print "ADVISEMENT: file [",i+1,"] is missing [",old,"]"
            else:
                # Add the filename as a key
                if self.has_key('file'):
                    self['file'] = gitm.dmarray(np.append(self['file'],
                                                          [gData.attrs['file']],
                                                          0),
                                                attrs=self['file'].attrs)
                else:
                    self['file'] = gitm.dmarray(np.array([gData.attrs['file']]),
                                                attrs={"name":
                                                           "GITM Bin Filename"})

                # Add the magnetic field filename as a key
                magfile = None
                if gData.attrs.has_key('magfile'):
                    magfile = gData.attrs['magfile']

                if self.has_key('magfile'):
                    self['magfile'] = gitm.dmarray(np.append(self['magfile'],
                                                             [magfile], 0),
                                                   attrs=self['magfile'].attrs)
                else:
                    self['magfile'] = gitm.dmarray(np.array([magfile]), attrs={"name":"GITM 3DMAG Filename"})

            # Initialize the first instance of each key and append the new data
            # if it is not
            for k in newkeys:
                if type(gData[k]) is dmarray:
                    scale = 0.0
                    for j in reversed(range(i+1)):
                        if j == 0:
                            scale = 1.0
                        if j == i:
                            data = np.array([gData[k]*scale])
                        else:
                            data = np.append(data, [gData[k]*scale], 0)
                elif i != 0:
                    print "WARNING: key [",k,"] not temporally aligned"
                else:
                    data = np.array([gData[k]])

                try:
                    self[k] = gitm.dmarray(data, attrs=gData[k].attrs)
                except AttributeError:
                    if k.find('time') >= 0:
                        self[k] = gitm.dmarray(data,
                                               attrs={"name":"Universal Time",
                                                      "units":"date",
                                                      "scale":"date"})

            tempk = gData.keys()
            for k in oldkeys:
                if gData.has_key(k):
                    self[k] = gitm.dmarray(np.append(self[k], [gData[k]], 0),
                                           attrs=self[k].attrs)
                elif gData.attrs.has_key(k):
                    self[k] = gitm.dmarray(np.append(self[k], [gData.attrs[k]],
                                                     0), attrs=self[k].attrs)
                elif k.find("magfile") >= 0:
                    self[k] = gitm.dmarray(np.append(self[k], [None], 0),
                                           attrs=self[k].attrs)
                else:
                    self[k] = gitm.dmarray(np.append(self[k],
                                                     [gData[tempk[0]] * np.nan],
                                                     0), attrs=self[k].attrs)

    def appendobs(self, nplist, keylist, match_method="none", dtdata=None,
                  londata=None, latdata=None, altdata=None, max_timedelt=300.0,
                  max_locdelt=100000.0, boxcar_sec=300.0, lat_unit="degrees",
                  lon_unit="degrees", alt_unit="km", *args, **kwargs):
        '''
        A routine to add observational data to a GitmTime data structure.
        There are several available methods to align the observational and model
        data, including none (assumes data is already aligned and structured),
        nearest neighbor, running average, and running median.  If any of the
        matching methods are used, the latitude, longitude, altitude, and
        datetime information must be supplied.

        Input:
        nplist       = list of numpy arrays containing the observational data
        keylist      = list of key names corresponding to the numpy arrays
        match_method = none/nearest/average/median (default = none)  No matching
                       appends the observational data without aligning the
                       times.  Nearest matches the observational data with the
                       GitmLoc data using the datetime array to find the closest
                       points.  A maximum time difference can be specified (in
                       seconds).  Average and Median computes a running average
                       or median using a specified temperal width for a boxcar
                       size and outputs points at every GitmLoc point.  If no
                       data is available, np.nan is inserted.

        Inputs if match_method is not equal to none:
        dtdata       = numpy array for the obs. datetime data (default is None)
        latdata      = numpy array for the obs. latitude (default is None)
        londata      = numpy array for the obs. longitude (default is None)
        altdata      = numpy array for the obs. altitude (default is None)
        max_timedelt = Maximum allowable seperation (in seconds) between nearest
                       neighbor points. (default = 300.0 s)
        max_locdelt  = Maximum allowable location seperation.  For nearest
                       neighbor method (default), this value is spedified in
                       units of altitude (m). For the average or median methods,
                       this should be a 3D list containing the longitude
                       seperation, latitude seperation, and altitude seperation
                       in units of the observable data
        boxcar_sec   = For the average or median match methods, a boxcar window
                       (in seconds) must be specified. (default = 300.0 s)
        lat_unit     = latitude units (default degrees)
        lon_unit     = longitude units (default degrees)
        alt_unit     = altitude units (default km)
        '''
        import gitm_loc_rout as glr
        import datetime as dt

        # Initialize the unit scaling and data dimension
        rlon = 1.0
        rlat = 1.0
        malt = 1.0
        dims = self['dLon'].shape

        if lon_unit.find("deg") >= 0:
            rlon = np.pi / 180.0
        if lat_unit.find("deg") >= 0:
            rlat = np.pi / 180.0
        if alt_unit.find("km") >= 0:
            malt = 1000.0

        # Initialize the output keys
        obs_keylist = list()
        for i,k in enumerate(keylist):
            obs_keylist.append("obs_{:s}".format(k))

        # Ensure the match method is in lowercase
        match_method = match_method.lower()

        if match_method.find("none") >= 0:
            # Append the desired data to the GitmLoc structure
            for i,k in enumerate(obs_keylist):
                self[k] = dc(nplist[i])        

        elif match_method.find("nearest") >= 0:
            # Match the nearest observational and GitmLoc data.  Begin by
            # initializing the desired arrays
            selfdelt = list([[max_timedelt + 1.0,
                              max_locdelt + 1.0, -1, -1]]) * len(self['time'])
            for k in obs_keylist:
                self[k] = np.ndarray(shape=dims) * np.nan

            # Cycle through the observations, aligning times and then
            # locations
            for i,t in enumerate(dtdata):
                timedelt, itime = glr.find_nearest_datetime(self['time'], t)

                if abs(timedelt) <= max_timedelt:
                    # If this point is acceptable, find the nearest location
                    # if location comparison indexes have been provided
                    if(latdata is not None or londata is not None or
                       altdata is not None):
                        x_list = list()
                        y_list = list()
                        z_list = list()

                        for j in range(dims[1] * dims[2] * dims[3]):
                            ilon = int(j / (dims[2] * dims[3]))
                            ilat = int((j - ilon * dims[2] * dims[3]) / dims[3])
                            ialt = int(j - (ilon * dims[2] + ilat) * dims[3])

                            if(not np.isnan(self['Longitude'][itime,ilon,ilat,
                                                              ialt])
                               and not np.isnan(self['Latitude'][itime, ilon,
                                                                 ilat, ialt])
                               and not np.isnan(self['Altitude'][itime,ilon,
                                                                 ilat, ialt])):
                                if altdata:
                                    x_list.append(self['Altitude'][itime, ilon,
                                                                   ilat, ialt])
                                else:
                                    x_list.append(1.0)

                                y_list.append(self['Longitude'][itime,ilon,ilat,
                                                                ialt])
                                z_list.append(self['Latitude'][itime,ilon,ilat,
                                                               ialt])
                        if altdata:
                            obsloc = [altdata[i]*malt, londata[i]*rlon,
                                      latdata[i]*rlat]
                        else:
                            obsloc = [1.0, londata[i]*rlon, latdata[i]*rlat]

                        (locdelt, j) = glr.find_nearest_location(x_list, y_list,
                                                                 z_list, obsloc,
                                                                 "sph")
                        if abs(locdelt) <= max_locdelt:
                            # If this pairing is within the desired range,
                            # test to see if it is closer than another pairing
                            # for the same self data point.  Prioritize time
                            # over location

                            if(selfdelt[itime][0] > abs(timedelt) or
                               (selfdelt[itime][0] == abs(timedelt) and
                                selfdelt[itime][1] > abs(locdelt))):
                                selfdelt[itime]=[abs(timedelt),abs(locdelt),i,j]
                    else:
                        # If this pairing is within the desired range,
                        # test to see if it is closer than another pairing
                        # for the same self data point.

                        if(selfdelt[itime][0] > abs(timedelt)):
                            selfdelt[itime]=[abs(timedelt),0.0,i,0]
                                
            for itime,delt in enumerate(selfdelt):
                if delt[2] >= 0:
                    # Then the time and location are appropriate,
                    # save the observational data at the appropriate index
                    ilon = int(delt[3] / (dims[2] * dims[3]))
                    ilat = int((delt[3] - ilon * dims[2] * dims[3]) / dims[3])
                    ialt = int(delt[3] - (ilon * dims[2] + ilat) * dims[3])
                    
                    for ik,k in enumerate(obs_keylist):
                        self[k][itime,ilon,ilat,ialt] = nplist[ik][delt[2]]

        elif match_method.find("average")<0 and match_method.find("method")<0:
            # This is an unknown matching method
            print "append_obs ERROR: unknown matching method", match_method

        else:
            # Compute a running median or average at the GitmLoc times.  Begin
            # by constructing the location lists
            nploc = np.array([np.array([londata[i]*rlon,latdata[i]*rlat,a*malt])
                              for i,a in enumerate(altdata)])

            # Convert max_locdelt from deg/km to rad/m
            max_locdelt[0] *= rlon
            max_locdelt[1] *= rlat
            max_locdelt[2] *= malt

            def get_itime(j):
                return int(j / (dims[1] * dims[2] * dims[3]))
            def get_ilon(j):
                d = dims[1] * dims[2] * dims[3]
                return int((j - int(j/d) * d) / (dims[2] * dims[3]))
            def get_ilat(j):
                d = dims[1] * dims[2] * dims[3]
                k = int((j - int(j/d)*d)/ (dims[2]*dims[3])) * dims[2] * dims[3]
                return int((j- int(j/d) * d - k)/dims[3])
            def get_ialt(j):
                d = dims[1] * dims[2] * dims[3]
                k = int((j - int(j/d)*d)/ (dims[2]*dims[3])) * dims[2] * dims[3]
                l = int((j - k - int(j/d) * d) / dims[3]) * dims[3]
                return int(j - int(j/d) * d - k - l)

            d = dims[0] * dims[1] * dims[2] * dims[3]
            selfloc = np.array([[self['Longitude'][get_itime(j),get_ilon(j),
                                                   get_ilat(j),get_ialt(j)],
                                 self['Latitude'][get_itime(j),get_ilon(j),
                                                  get_ilat(j),get_ialt(j)],
                                 self['Altitude'][get_itime(j),get_ilon(j),
                                                  get_ilat(j),get_ialt(j)]]
                                for j in range(d)])
            selftime = np.array([self['time'][get_itime(j)] for j in range(d)])
            locwrap = [2.0*np.pi, 0.0, 0.0]

            if match_method.find("average") >= 0:
                yout,ystd,nout = glr.match_running_average(nplist, dtdata,
                                                           nploc, selftime,
                                                           selfloc, boxcar_sec,
                                                           max_locdelt, 0,
                                                           locwrap)
            else:
                yout,yup,ydown,nout = glr.match_running_average(nplist, dtdata,
                                                                nploc, selftime,
                                                                selfloc,
                                                                boxcar_sec,
                                                                max_locdelt, 0,
                                                                locwrap)
            
            # Assign output
            for i,k in enumerate(obs_keylist):
                self[k] = np.ndarray(shape=dims, buffer=yout[i])
                okey = "{:s}_nout".format(k)
                self[okey] = np.ndarray(shape=dims, dtype=int, buffer=nout[i])

                if match_method.find("average") >= 0:
                    okey = "{:s}_sig".format(k)
                    self[okey] = np.ndarray(shape=dims, buffer=ystd[i])
                else:
                    okey = "{:s}_terr".format(k)
                    self[okey] = np.ndarray(shape=dims, buffer=yup[i])
                    okey = "{:s}_berr".format(k)
                    self[okey] = np.ndarray(shape=dims, buffer=ydown[i])

        return obs_keylist

    def sat_dateloc_ticks(self, x, pos):
        '''
        Define ticks to include all the information necessary to know where
        measurements lie in spacetime.  This is most useful for data output
        along a satellite track.
        '''
        from matplotlib.dates import num2date
        import datetime as dt
        
        # Format the UT date
        nt = num2date(x, tz=None)
        nowtime = dt.datetime(nt.year, nt.month, nt.day, nt.hour, nt.minute, 
                              nt.second, tzinfo=None)
        nowstring = str(nowtime)
        fmtstring = "{:s}".format(nowstring.replace(" ", "\n "))

        # Get the time index 
        deltime = (self['time']-nowtime)
        mintime = min(abs(deltime))
        try:
            # Get the index corresponding to this time
            index = (deltime.tolist()).index(mintime)
        except ValueError:
            # Only the UT is available at this time.
            return(fmtstring)

        # Get the desired coordinates
        lth = int(self['LT'][index,0,0,0])
        ltm = int((self['LT'][index,0,0,0] - lth) * 60.0)
        lat = self['dLat'][index,0,0,0]
        lon = self['dLon'][index,0,0,0]
        alt = self['Altitude'][index,0,0,0] / 1000.0

        # Build the format string
        fmtstring = "{:s}\n {:02d}:{:02d}\n{: 3.1f}\n{:4.1f}\n{:4.0f}".format(fmtstring, lth, ltm, lat, lon, alt)

        # Add magnetic latitude if possible
        if self.has_key('Magnetic Latitude'):
            fmtstring = "{:s}\n{:.1f}".format(fmtstring, self['Magnetic Latitude'][index,0,0,0])

        return(fmtstring)
#End Class

def load_multiple_gitm_bin(filelist, magfile=None, *args, **kwargs):
    '''
    Loads a list of GITM binary files into their own GitmBin data structures.
    The list may be an ascii file containing a list of files or a list object.
    A list of the data structures is returned.  A 3DION or 3DMAG file may
    be specified for the entire list using the keyword arguement "magfile".

    Input:
    filelist = python list of file names or an ASCII file containing a list
               of filenames
    magfile = 3DMAG or 3DION file (default=None)
    '''
    # import local packages
    from os import path
    outlist = list()

    if type(filelist) == str:
        func_name = "load_multiple_gitm_bin"
        fsize = path.getsize(filelist)

        if(fsize > 2.0e9):
            print func_name,"ERROR: File list size [",(fsize*1e-9),"GB > 2 GB]"
            return outlist
        elif(fsize == 0.0):
            print func_name, "ERROR: empty file list [", filelist, "]"
            return outlist

        # Read in the list of files
        f = open(filelist, "r")
        namelist = f.readlines()
    else:
        namelist = filelist

    for name in namelist:
        name = string.strip(name)
        outlist.append(gitm.GitmBin(name, magfile=magfile))

        if outlist[-1].attrs['nAlt'] > 1:
            outlist[-1].calc_2dion()

    return(outlist)
# End load_multiple_gitm_bin

def set_sat_dateloc_label(ax, mlat=True):
    '''
    Create a label for the GitmTime.sat_dateloc_ticks formatted Ticks

    Input:
    ax   = subplot axis handle
    mlat = Is magnetic latitude included? default is True
    '''
    slabel = "Date\nTime\nSLT\nGLat (deg)\nGlon (deg)\nAlt (km)"
    yoff = -0.68

    if mlat:
        slabel = "{:s}\nMLat (deg)".format(slabel)
        yoff -= 0.1

    ax.text(1.01, yoff, slabel, transform=ax.transAxes)
# END set_sat_dateloc_label
