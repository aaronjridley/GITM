#!/usr/bin/env python

from glob import glob
from datetime import datetime
from datetime import timedelta
from struct import unpack
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from pylab import cm
import os


def write_line(fp, sLine):
    sLineN = sLine + '\n'
    fp.write(sLineN.encode())
    return

#-----------------------------------------------------------------------------
# Here is a useful routine for writing out data to a log file
#-----------------------------------------------------------------------------

def write_log(data, fileHeader = 'log', message = ''):

    sTime = data['times'][0].strftime('%Y%m%d')
    fileout = fileHeader + '_' + sTime + '.txt'
    print('--> Writing file : ' + fileout)
    fp = open(fileout, 'wb')

    write_line(fp, '')
    write_line(fp, message)
    
    pwd = os.getcwd()
    write_line(fp, '')
    write_line(fp, '#DIRECTORY')
    write_line(fp, pwd)

    if ('alt' in data):
        if (np.isscalar(data['alt'])):
            alt = data['alt']
        else:
            alt = mean(np.array(data['alt']))
        sLine = '%f' % alt
        write_line(fp, '')
        write_line(fp, '#ALTITUDE')
        write_line(fp, sLine)
    
    write_line(fp, '')
    write_line(fp, '#VARIABLES')
    write_line(fp, 'Year')
    write_line(fp, 'Month')
    write_line(fp, 'Day')
    write_line(fp, 'Hour')
    write_line(fp, 'Minute')
    write_line(fp, 'Second')
    for key in data.keys():
        if ((key != 'times') and (key != 'alt')):
            write_line(fp, key)

    write_line(fp, '')
    write_line(fp, '#START')

    for i, t in enumerate(data['times']):

        sLine = t.strftime(' %Y %m %d %H %M %S')

        for key in data.keys():
            if ((key != 'times') and (key != 'alt')):
                sLine = sLine + ' %e' % data[key][i]
        write_line(fp, sLine)
        
    fp.close()
    
    return



#-----------------------------------------------------------------------------

def read_gitm_header(file):
    r""" Grab ancillary information from the GITM file

    Parameters
    ----------
    file - name of the file to read and get header from

    Returns
    -------
    header: A dictionary containing information about the netcdf file, such
            as nLons, nLons, nAlts, nVars, variable names, time(s)

    """
    
    if (len(file) == 0):

        filelist = sorted(glob('./3DALL*.bin'))

        if (len(filelist) == 0):
            print("No 3DALL files found. Checking for 1DALL.")
            filelist = glob('./1DALL*.bin')
            if (len(filelist) == 0):
                print("No 1DALL files found. Stopping.")
                exit()
            file = filelist[0]

    else:
        filelist = glob(file[0])
        file = filelist[0]
            
    header = {"nFiles": len(filelist), \
              "version": 0.1, \
              "nLons": 0, \
              "nLats": 0, \
              "nAlts": 0, \
              "nVars": 0, \
              "vars": [], \
              "time": [], \
              "filename": [file] }

    f=open(file, 'rb')

    # This is all reading header stuff:

    endChar='>'
    rawRecLen=f.read(4)
    recLen=(unpack(endChar+'l',rawRecLen))[0]
    if (recLen>10000)or(recLen<0):
        # Ridiculous record length implies wrong endian.
        endChar='<'
        recLen=(unpack(endChar+'l',rawRecLen))[0]

    # Read version; read fortran footer+header.
    header["version"] = unpack(endChar+'d',f.read(recLen))[0]

    (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

    # Read grid size information.
    (header["nLons"],header["nLats"],header["nAlts"]) = \
        unpack(endChar+'lll',f.read(recLen))
    (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

    # Read number of variables.
    header["nVars"]=unpack(endChar+'l',f.read(recLen))[0]
    (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

    # Collect variable names.
    for i in range(header["nVars"]):
        v = unpack(endChar+'%is'%(recLen),f.read(recLen))[0]
        header["vars"].append(v.decode('utf-8').replace(" ",""))
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

    print(recLen)
    
        
    # Extract time. 
    (yy,mm,dd,hh,mn,ss,ms)=unpack(endChar+'lllllll',f.read(recLen))
    header["time"].append(datetime(yy,mm,dd,hh,mn,ss,ms*1000))
    # print(header["time"][-1])

    f.close()

    return header

#-----------------------------------------------------------------------------

def read_gitm_headers(pre='./3DALL'):
    r""" Grab ancillary information from the GITM files

    Parameters
    ----------
    None - does a glob to find files.

    Returns
    -------
    header: A dictionary containing information about the netcdf file, such
            as nLons, nLons, nAlts, nVars, variable names, time(s)

    """

    filelist = sorted(glob(pre+'*.bin'))
    print("Found ", len(filelist), "files")
    
    header = {"nFiles": len(filelist), \
              "version": 0.1, \
              "nLons": 0, \
              "nLats": 0, \
              "nAlts": 0, \
              "nVars": 0, \
              "vars": [], \
              "time": [], \
              "filename": [] }

    for file in filelist:

        header["filename"].append(file)

        f=open(file, 'rb')

        # This is all reading header stuff:

        endChar='>'
        rawRecLen=f.read(4)
        recLen=(unpack(endChar+'l',rawRecLen))[0]
        if (recLen>10000)or(recLen<0):
            # Ridiculous record length implies wrong endian.
            endChar='<'
            recLen=(unpack(endChar+'l',rawRecLen))[0]

        # Read version; read fortran footer+header.
        header["version"] = unpack(endChar+'d',f.read(recLen))[0]

        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Read grid size information.
        (header["nLons"],header["nLats"],header["nAlts"]) = \
            unpack(endChar+'lll',f.read(recLen))
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Read number of variables.
        header["nVars"]=unpack(endChar+'l',f.read(recLen))[0]
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Collect variable names.
        for i in range(header["nVars"]):
            v = unpack(endChar+'%is'%(recLen),f.read(recLen))[0]
            if (file == filelist[0]):
                header["vars"].append(v.decode('utf-8').replace(" ",""))
            (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Extract time. 
        (yy,mm,dd,hh,mn,ss,ms)=unpack(endChar+'lllllll',f.read(recLen))
        header["time"].append(datetime(yy,mm,dd,hh,mn,ss,ms*1000))

        f.close()

    return header

#-----------------------------------------------------------------------------

def read_gitm_one_file(file_to_read, vars_to_read=-1):
    r""" Read list of variables from one GITM file

    Parameters
    ----------
    file_to_read: GITM file to read
    vars_to_read: list of variable NUMBERS to read

    Returns
    -------
    data["time"]: datetime of the file
    data[NUMBER]: data that is read in.
                  NUMBER goes from 0 - number of vars read in (0-3 typical)
    (Also include header information, as described above)
    """

    print("Reading file : "+file_to_read)

    data = {"version": 0, \
            "nLons": 0, \
            "nLats": 0, \
            "nAlts": 0, \
            "nVars": 0, \
            "time": 0, \
            "vars": []}

    f=open(file_to_read, 'rb')

    # This is all reading header stuff:

    endChar='>'
    rawRecLen=f.read(4)
    recLen=(unpack(endChar+'l',rawRecLen))[0]
    if (recLen>10000)or(recLen<0):
        # Ridiculous record length implies wrong endian.
        endChar='<'
        recLen=(unpack(endChar+'l',rawRecLen))[0]

    # Read version; read fortran footer+data.
    data["version"] = unpack(endChar+'d',f.read(recLen))[0]

    (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

    # Read grid size information.
    (data["nLons"],data["nLats"],data["nAlts"]) = \
        unpack(endChar+'lll',f.read(recLen))
    (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

    # Read number of variables.
    data["nVars"]=unpack(endChar+'l',f.read(recLen))[0]
    (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

    if (vars_to_read[0] == -1):
        vars_to_read = np.arange[nVars]

    # Collect variable names.
    for i in range(data["nVars"]):
        data["vars"].append(unpack(endChar+'%is'%(recLen),f.read(recLen))[0])
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

    # Extract time. 
    (yy,mm,dd,hh,mn,ss,ms)=unpack(endChar+'lllllll',f.read(recLen))
    data["time"] = datetime(yy,mm,dd,hh,mn,ss,ms*1000)
    #print(data["time"])

    # Header is this length:
    # Version + start/stop byte
    # nLons, nLats, nAlts + start/stop byte
    # nVars + start/stop byte
    # Variable Names + start/stop byte
    # time + start/stop byte

    iHeaderLength = 8 + 4+4 + 3*4 + 4+4 + 4 + 4+4 + \
        data["nVars"]*40 + data["nVars"]*(4+4) + 7*4 + 4+4

    nTotal = data["nLons"]*data["nLats"]*data["nAlts"]
    iDataLength = nTotal*8 + 4+4

    for iVar in vars_to_read:
        f.seek(iHeaderLength+iVar*iDataLength)
        s=unpack(endChar+'l',f.read(4))[0]
        data[iVar] = np.array(unpack(endChar+'%id'%(nTotal),f.read(s)))
        data[iVar] = data[iVar].reshape( 
            (data["nLons"],data["nLats"],data["nAlts"]),order="F")

    f.close()

    return data


#-----------------------------------------------------------------------------
# locations should contain :
#   "lats" -> latitudes in degrees
#   "lons" -> longitudes in degrees
#   "alts" -> altitudes in km

def gitm_extract_data(locations, times, variables, dir = './'):

    headers = read_gitm_headers(dir + '3DALL')

    file = headers["filename"][0]

    outdata = {}
    for v in variables:
        outdata[v] = []
    
    # get lon/lat/alt:
    vars = [0,1,2]
    print("Reading first GITM file to get location info")
    data = read_gitm_one_file(file, vars)

    Alts = data[2][0][0]/1000.0
    Lons = np.rad2deg(data[0][:,0,0])
    Lats = np.rad2deg(data[1][0,:,0])

    minLat = np.min(Lats)
    minLon = np.min(Lons)

    dLat = Lats[1]-Lats[0]
    dLon = Lons[1]-Lons[0]

    nFiles = len(headers["time"])
    nPts = len(locations["lats"])

    SatLats = locations["lats"]
    SatLons = locations["lons"]
    SatAlts = locations["alts"]

    nAlts = len(Alts)

    iSatLats = (SatLats-minLat)/dLat
    rSatLats = np.array(iSatLats - np.floor(iSatLats))
    iSatLats = (np.array((np.floor(iSatLats)))).astype(int)

    iSatLons = (SatLons-minLon)/dLon
    rSatLons = iSatLons - np.floor(iSatLons)
    iSatLons = (np.array((np.floor(iSatLons)))).astype(int)

    # We need to figure out Alts, which are not uniform:
    
    SatMinAlt = np.min(SatAlts)
    iSatMinAlt = 0
    while (SatMinAlt > Alts[iSatMinAlt]):
        iSatMinAlt = iSatMinAlt + 1
    if (iSatMinAlt > 0):
        iSatMinAlt = iSatMinAlt-1

    iSatAlts = []
    rSatAlts = []
    for iPt in range(nPts):

        if (SatAlts[iPt] < Alts[0]):
            iSatAlts.append(0)
            rSatAlts.append(0.0)
        else:
            if (SatAlts[iPt] > Alts[nAlts-1]):
                iSatAlts.append(nAlts-2)
                rSatAlts.append(1.0)
            else:
                i = iSatMinAlt
                while (SatAlts[iPt] > Alts[i]):
                    i = i + 1
                i = i - 1
                r = (SatAlts[iPt] - Alts[i])/(Alts[i+1]-Alts[i])
                iSatAlts.append(i)
                rSatAlts.append(r)


    print("Reading files to extract data...")
    for iPt in range(nPts):

        iFileSave = -1
        for iFile in range(nFiles):
            if (times[iPt] >= headers["time"][iFile]):
                iFileSave = iFile
        if (iFileSave == nFiles-1 or iFileSave == -1):
            i = nFiles-2
            r = 1.0
        else:
            i = iFileSave
            r = (times[iPt] - headers["time"][i]).total_seconds() / (headers["time"][i+1] - headers["time"][i]).total_seconds()

        if (iPt == 0):
            iFileLeft = i
            iFileRight = i+1
            fileLeft = headers["filename"][iFileLeft]
            dataLeft = read_gitm_one_file(fileLeft, variables)
            fileRight = headers["filename"][iFileRight]
            dataRight = read_gitm_one_file(fileRight, variables)

        if (iFileLeft != i):
            if (iFileRight == i):
                # Move right to left and grab new data:
                fileLeft = fileRight
                dataLeft = dataRight
                iFileLeft = iFileRight
                iFileRight = i+1
                fileRight = headers["filename"][iFileRight]
                dataRight = read_gitm_one_file(fileRight, variables)
            else:
                # may have skipped a few files or something, so reset:
                iFileLeft = i
                iFileRight = i+1
                fileRight = headers["filename"][iFileRight]
                dataRight = read_gitm_one_file(fileRight, variables)
                fileLeft = headers["filename"][iFileLeft]
                dataLeft = read_gitm_one_file(fileLeft, variables)
                print("Skipping to a new set of files : ", fileLeft, fileRight)

        # At this point, the data should be good, so we do interpolation:
        
        # First get lat, lon, alt for the two data sets:
        
        iLon = iSatLons[iPt]
        rLon = rSatLons[iPt]
        iLat = iSatLats[iPt]
        rLat = rSatLats[iPt]

        iAlt = iSatAlts[iPt]
        rAlt = rSatAlts[iPt]

        data = []

        for iVar in variables:
            dataLeft0 = ( (1.0-rLon) * (1.0-rLat) * (1.0-rAlt) * dataLeft[iVar][iLon  ][iLat  ][iAlt] +
                          (1.0-rLon) * (1.0-rLat) * (    rAlt) * dataLeft[iVar][iLon  ][iLat  ][iAlt+1] +
                          (1.0-rLon) * (    rLat) * (1.0-rAlt) * dataLeft[iVar][iLon  ][iLat+1][iAlt] +
                          (1.0-rLon) * (    rLat) * (    rAlt) * dataLeft[iVar][iLon  ][iLat+1][iAlt+1] +
                          (    rLon) * (1.0-rLat) * (1.0-rAlt) * dataLeft[iVar][iLon+1][iLat  ][iAlt] +
                          (    rLon) * (1.0-rLat) * (    rAlt) * dataLeft[iVar][iLon+1][iLat  ][iAlt+1] +
                          (    rLon) * (    rLat) * (1.0-rAlt) * dataLeft[iVar][iLon+1][iLat+1][iAlt] +
                          (    rLon) * (    rLat) * (    rAlt) * dataLeft[iVar][iLon+1][iLat+1][iAlt+1] )

            dataRight0 = ( (1.0-rLon) * (1.0-rLat) * (1.0-rAlt) * dataRight[iVar][iLon  ][iLat  ][iAlt] +
                           (1.0-rLon) * (1.0-rLat) * (    rAlt) * dataRight[iVar][iLon  ][iLat  ][iAlt+1] +
                           (1.0-rLon) * (    rLat) * (1.0-rAlt) * dataRight[iVar][iLon  ][iLat+1][iAlt] +
                           (1.0-rLon) * (    rLat) * (    rAlt) * dataRight[iVar][iLon  ][iLat+1][iAlt+1] +
                           (    rLon) * (1.0-rLat) * (1.0-rAlt) * dataRight[iVar][iLon+1][iLat  ][iAlt] +
                           (    rLon) * (1.0-rLat) * (    rAlt) * dataRight[iVar][iLon+1][iLat  ][iAlt+1] +
                           (    rLon) * (    rLat) * (1.0-rAlt) * dataRight[iVar][iLon+1][iLat+1][iAlt] +
                           (    rLon) * (    rLat) * (    rAlt) * dataRight[iVar][iLon+1][iLat+1][iAlt+1] )

            outdata[iVar].append((1.0-r) * dataLeft0 + r * dataRight0)
            if (iVar < 2):
                outdata[iVar][-1] = outdata[iVar][-1] * rtod
            if (iVar == 2):
                outdata[iVar][-1] = outdata[iVar][-1] / 1000.0


    return outdata

