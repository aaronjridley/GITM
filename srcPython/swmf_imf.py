#!/usr/bin/env python

import re
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates

#------------------------------------------------------------------------------
# read swmf file
#------------------------------------------------------------------------------

def read_swmf_file(file):

    fpin = open(file, 'r')

    Vars = ["bx", "by", "bz", "vx", "vy", "vz", "n", "t"]
    
    data = {"Vars" : Vars,
            "nVars" : len(Vars),
            "times" : []}

    for v in Vars:
        data[v] = []
    
    for line in fpin:
        m = re.match(r'#START',line)
        if m:
            break

    for line in fpin:
        aline = line.split()
        year = int(aline[0])
        month = int(aline[1])
        day = int(aline[2])
        hour = int(aline[3])
        minute = int(aline[4])
        second = int(aline[5])
        t = dt.datetime(year, month, day, hour, minute, second)
        
        data["times"].append(t)
        for i,v in enumerate(data["Vars"]):
            data[v].append(float(aline[i + 7]))
    
    return data

#------------------------------------------------------------------------------
# write swmf file
#------------------------------------------------------------------------------

def write_swmf_imf_file(data, fileout,
                        message = 'output from write_swmf_imf_file\n'):

    fp = open(fileout, 'wb')

    fp.write("\n".encode())
    fp.write(message.encode())
    fp.write("\n".encode())
    fp.write("No AlphaRatio was considered\n".encode())
    fp.write("\n".encode())
    fp.write("#TIMEDELAY\n".encode())
    fp.write("0.0\n".encode())
    fp.write("\n".encode())
    fp.write("#START\n".encode())
    
    i = 0
    while (np.abs(data["vx"][i]) > 10000.0):
        i += 1
    
    vx0 = data["vx"][i]
    vy0 = data["vy"][i]
    vz0 = data["vz"][i]
    n0 = data["n"][i]
    temp0 = data["t"][i]
    
    for i, t in enumerate(data["times"]):
        bx = data["bx"][i]
        by = data["by"][i]
        bz = data["bz"][i]
        vx = data["vx"][i]
        vy = data["vy"][i]
        vz = data["vz"][i]
        n = data["n"][i]
        temp = data["t"][i]

        if (np.abs(vx) > 10000):
            vx = vx0
            vy = vy0
            vz = vz0
            n = n0
            temp = temp0
        else:
            vx0 = vx
            vy0 = vy
            vz0 = vz
            n0 = n
            temp0 = temp
            
        if (np.abs(bx) < 1000):
            sTime = t.strftime(' %Y %m %d %H %M %S 000')
            sImf = "%8.2f %8.2f %8.2f" % (bx, by, bz)
            sSWV = "%9.2f %9.2f %9.2f" % (vx, vy, vz)
            sSWnt = "%8.2f %11.1f" % (n, temp)

        line = sTime + sImf + sSWV + sSWnt + "\n"

        fp.write(line.encode())

    fp.close()

#------------------------------------------------------------------------------
# plot swmf file
#------------------------------------------------------------------------------

def plot_imf(data):

    fig = plt.figure(figsize = (10,10))
    zeros = np.array(data["bz"]) * 0.0

    ax = fig.add_subplot(511)
    ax.plot(data["times"], data["bz"])
    ax.plot(data["times"], zeros, 'k:')
    ax.set_ylabel('IMF Bz (nT)')
    ax.set_xlim(data["times"][0],data["times"][-1])

    ax = fig.add_subplot(512)
    ax.plot(data["times"], data["by"])
    ax.plot(data["times"], zeros, 'k:')
    ax.set_ylabel('IMF By (nT)')
    ax.set_xlim(data["times"][0],data["times"][-1])

    ax = fig.add_subplot(513)
    ax.plot(data["times"], data["vx"])
    ax.set_ylabel('SW Vx (km/s)')
    ax.set_xlim(data["times"][0], data["times"][-1])

    ax = fig.add_subplot(514)
    ax.plot(data["times"], data["n"])
    ax.set_ylabel('SW N (/cm3)')
    ax.set_xlim(data["times"][0], data["times"][-1])

    vx = np.array(data["vx"])
    vy = np.array(data["vy"])
    vz = np.array(data["vz"])

    v = np.sqrt(vx*vx + vy*vy + vz*vz) * 1000.0

    bx = np.array(data["bx"])
    by = np.array(data["by"])
    bz = np.array(data["bz"])

    b = np.sqrt(bx*bx + by*by + bz*bz) / 1.0e9

    n = np.array(data["n"]) * 1e6

    mp = 1.67e-27

    mu0 = 1.256e-6

    ca = b / np.sqrt(mu0 * n * mp)
    ma = v / ca

    ax = fig.add_subplot(515)
    ax.plot(data["times"], ma)
    ax.set_ylabel('SW Ma')
    ax.set_xlim(data["times"][0], data["times"][-1])
    ax.set_ylim(0.0, 20.0)

    plotfile = data["times"][0].strftime('imf%Y%m%d.png')
    fig.savefig(plotfile)
    plt.close()
    
