#!/usr/bin/env python

import requests
import re
import datetime as dt
import numpy as np

# This is derived from https://omniweb.gsfc.nasa.gov/html/omni_min_data.html
    
#inputs: dates and info desired
#date1 and date2 format example: "20110620" means June 20, 2011
def download_omni_data(date1, date2, info):

    url_i = "https://omniweb.gsfc.nasa.gov/cgi/nx1.cgi"
    params_i = "activity=retrieve&res=min&spacecraft=omni_min&start_date=" + date1 + \
             "&end_date=" + date2
    if info == "-imf":
        params_i += "&vars=14&vars=17&vars=18"
    elif info == "-ae":
        params_i += "&vars=37&vars=38&vars=39"
    elif info == "-sw":
        params_i += "&vars=22&vars=23&vars=24&vars=25&vars=26"
    elif info == "-all":
        params_i += "&vars=14&vars=17&vars=18&vars=22&vars=23&vars=24&vars=25&vars=26&vars=37&vars=38&vars=39"
    else:   #if nothing is inputted for desired info
        params_i += "&vars=14&vars=17&vars=18&vars=22&vars=23&vars=24&vars=25&vars=26&vars=37&vars=38&vars=39"
    params_i += "&back="

    res = requests.get(url = url_i, params = params_i)
    return res.text     #this returns the info file from the get request

def parse_omni_data(results):

    names = {'BX, nT (GSE, GSM)' : 'bx', 
             'BY, nT (GSM)' : 'by',
             'BZ, nT (GSM)' : 'bz',
             'Vx Velocity,km/s' : 'vx',
             'Vy Velocity, km/s' : 'vy',
             'Vz Velocity, km/s' : 'vz',
             'Proton Density, n/cc' : 'n',
             'Temperature, K' : 't',
             'Proton Temperature, K' : 't',
             'AE-index, nT' : 'ae',
             'AL-index, nT' : 'al',
             'AU-index, nT' : 'au'}
    
    lines = results.splitlines()

    IsFound = 0
    iLine = 0
    while (not IsFound):
        line = lines[iLine]
        m = re.match(r'<B>.*',line)
        if m:
            IsFound = 1
        iLine += 1
    
    data = {}
    data["Vars"] = []
    data["nVars"] = 0
    data["times"] = []
    
    IsFound = 0
    while (not IsFound):
        line = lines[iLine]
        if (len(line) < 2):
            IsFound = 1
        else:
            v = names[line[3:]]
            data["Vars"].append(v)
            data[v] = []
            data["nVars"] += 1
        iLine += 1

    # Skip over YYYY line:
    iLine += 1

    nV = data["nVars"] + 4
    IsFound = 0
    while (not IsFound):
        aline = lines[iLine].split()
        if (len(aline) < nV):
            IsFound = 1
        else:
            year = int(aline[0])
            doy = int(aline[1])-1
            hour = int(aline[2])
            minute = int(aline[3])
            base = dt.datetime(year,1,1,hour,minute,0)
            actual = base + dt.timedelta(days=doy)
            data["times"].append(actual)
            for i,v in enumerate(data["Vars"]):
                data[v].append(float(aline[i + 4]))
        iLine += 1

    return data

def clean_omni(data):

    newdata = {"Vars" : data["Vars"],
               "nVars" : data["nVars"],
               "times" : [], 
               "bx" : [], 
               "by" : [], 
               "bz" : [], 
               "vx" : [], 
               "vy" : [], 
               "vz" : [], 
               "n" : [], 
               "t" : [], 
               "ae" : [], 
               "au" : [], 
               "al" : []}
    
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
            newdata["times"].append(t)
            newdata["bx"].append(bx)
            newdata["by"].append(by)
            newdata["bz"].append(bz)
            newdata["vx"].append(vx)
            newdata["vy"].append(vy)
            newdata["vz"].append(vz)
            newdata["n"].append(n)
            newdata["t"].append(temp)
            newdata["ae"].append(data["ae"][i])
            newdata["au"].append(data["au"][i])
            newdata["al"].append(data["al"][i])

    return newdata
