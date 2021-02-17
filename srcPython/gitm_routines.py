#!/usr/bin/env python

from glob import glob
from datetime import datetime
from datetime import timedelta
from struct import unpack
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from pylab import cm


#-----------------------------------------------------------------------------
# 
#-----------------------------------------------------------------------------

def read_gitm_header(file):

    if (len(file) == 0):

        filelist = glob('./3DALL*.bin')

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
            
    header = {}
    header["nFiles"] = len(filelist)
    header["version"] = 0
    header["nLons"] = 0
    header["nLats"] = 0
    header["nAlts"] = 0
    header["nVars"] = 0
    header["vars"] = []
    header["time"] = []
    header["filename"] = []

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
    (header["nLons"],header["nLats"],header["nAlts"]) = unpack(endChar+'lll',f.read(recLen))
    (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

    # Read number of variables.
    header["nVars"]=unpack(endChar+'l',f.read(recLen))[0]
    (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

    # Collect variable names.
    for i in range(header["nVars"]):
        v = unpack(endChar+'%is'%(recLen),f.read(recLen))[0]
        header["vars"].append(v.decode('utf-8').replace(" ",""))
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

    # Extract time. 
    (yy,mm,dd,hh,mn,ss,ms)=unpack(endChar+'lllllll',f.read(recLen))
    header["time"].append(datetime(yy,mm,dd,hh,mn,ss,ms*1000))
    # print(header["time"][-1])

    f.close()

    return header

#-----------------------------------------------------------------------------
# 
#-----------------------------------------------------------------------------

def read_gitm_headers():

    filelist = sorted(glob('./3DALL*.bin'))

    header = {}
    header["nFiles"] = len(filelist)
    header["version"] = 0
    header["nLons"] = 0
    header["nLats"] = 0
    header["nAlts"] = 0
    header["nVars"] = 0
    header["vars"] = []
    header["time"] = []
    header["filename"] = []

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
        (header["nLons"],header["nLats"],header["nAlts"]) = unpack(endChar+'lll',f.read(recLen))
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
        # print(header["time"][-1])

        f.close()

    return header

#-----------------------------------------------------------------------------
# 
#-----------------------------------------------------------------------------

def read_gitm_one_file(file_to_read, vars_to_read=-1):

    print("Reading file : "+file_to_read)

    data = {}
    data["version"] = 0
    data["nLons"] = 0
    data["nLats"] = 0
    data["nAlts"] = 0
    data["nVars"] = 0
    data["time"] = 0
    data["vars"] = []

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
    (data["nLons"],data["nLats"],data["nAlts"]) = unpack(endChar+'lll',f.read(recLen))
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

    iHeaderLength = 8 + 4+4 + 3*4 + 4+4 + 4 + 4+4 + data["nVars"]*40 + data["nVars"]*(4+4) + 7*4 + 4+4

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
# 
#-----------------------------------------------------------------------------

