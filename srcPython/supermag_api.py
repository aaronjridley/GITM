# This is not written by UM - this is from APL / SuperMag project!

import urllib.request

# the 'certifi' library is required at APL and other sites that
# require SSL certs for web fetches.  If you need this, install certifi
# (pip install certifi)
import importlib
certspec = importlib.util.find_spec("certifi")
found = certspec is not None
if found: import certifi

import pandas as pd  # dataframes and also to_datetime
import json
import re
import datetime

"""
;supermag-api.py
; ================
; Author S. Antunes, based on supermag-api.pro by R.J.Barnes


; (c) 2021  The Johns Hopkins University Applied Physics Laboratory
;LLC.  All Rights Reserved. 

;This material may be only be used, modified, or reproduced by or for
;the U.S. Government pursuant to the license rights granted under the 
;clauses at DFARS 252.227-7013/7014 or FAR 52.227-14. For any other
;permission, 
;please contact the Office of Technology Transfer at JHU/APL.

; NO WARRANTY, NO LIABILITY. THIS MATERIAL IS PROVIDED "AS IS."
; JHU/APL MAKES NO REPRESENTATION OR WARRANTY WITH RESPECT TO THE
; PERFORMANCE OF THE MATERIALS, INCLUDING THEIR SAFETY, EFFECTIVENESS, 
; OR COMMERCIAL VIABILITY, AND DISCLAIMS ALL WARRANTIES IN THE
; MATERIAL, WHETHER EXPRESS OR IMPLIED, INCLUDING (BUT NOT LIMITED TO)
; ANY AND ALL IMPLIED WARRANTIES OF PERFORMANCE, MERCHANTABILITY,
; FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT OF
; INTELLECTUAL PROPERTY OR OTHER THIRD PARTY RIGHTS. ANY USER OF THE
; MATERIAL ASSUMES THE ENTIRERISK AND LIABILITY FOR USING THE
; MATERIAL. IN NO EVENT SHALL JHU/APL BE LIABLE TO ANY USER OF THE
; MATERIAL FOR ANY ACTUAL, INDIRECT, CONSEQUENTIAL, SPECIAL OR OTHER
; DAMAGES ARISING FROM THE USE OF, OR INABILITY TO USE, THE MATERIAL,
; INCLUDING, BUT NOT LIMITED TO, ANY DAMAGES FOR LOST PROFITS. 
"""
# Sample URLs, type into browser if you want to compare the data vs python
#https://supermag.jhuapl.edu/services/data-api.php?fmt=json&logon=YOURNAME&start=2019-10-15T10:40&extent=3600&all&station=NCK
#https://supermag.jhuapl.edu/services/indices.php?fmt=json&logon=YOURNAME&start=2019-10-15T10:40&extent=3600&all
#https://supermag.jhuapl.edu/services/inventory.php?fmt=json&logon=YOURNAME&start=2019-10-15T10:40&extent=3600


def sm_coreurl(page,logon,start,extent):
  # internal helper function
  baseurl = "https://supermag.jhuapl.edu/"

  mytime = sm_parsestart(start)
  urlstr = baseurl + 'services/'+page+'?python&nohead'
  urlstr+='&start='+mytime
  urlstr += '&logon='+logon
  urlstr+='&extent='+ ("%12.12d" % extent)

  #print("debug:",urlstr)

  return(urlstr)

# handy helper function when using complicated CVS encoding of lists
def sm_csvitem_to_list(myarr):
  # converts entity of form ['HOP', 'NVS', 'IRT'] to an actual list of HOP, NVS, IRT
  mylist=[]
  for myline in myarr:
    myline=re.sub("'","",myline[1:-1])
    mylist.append(myline.split(", "))
  return(mylist)


# handy helper function when using complicated CVS encoding of dicts
def sm_csvitem_to_dict(myarr,**kwargs):
  # converts entity of form {'X': -12.213, 'Y': -5.5, 'Z': 1.2} to an actual dict of var.X, var.Y, etc
  # items are presumed strings by default
  # optional argument "convert=num" will convert them to doubles
  mylist=[]
  for myline in myarr:
    myline2=re.sub(" ","",myline[1:-1]) # scrub out extra spaces
    myline2=re.sub("'","",myline2)
    elements = dict(x.split(":") for x in myline2.split(","))
    # little sanity check to make sure float subitems remain floats
    try: elements = {item: float(value) for (item, value) in elements.items()}
    except: pass
    
    #type(elements)
    mylist.append(elements)
  return(mylist)

def sm_parsestart(start):
  # internal helper function
  # takes either list of [yyyy, mo, dd, hh, mm, opt_ss]
  # or string of a normal datetime 'YYYY-MM-DD hh-mm' (optional ss)
  # or the SuperMAG-ready 'YYYY-MM-DDThh-mm-ss'

  if isinstance(start,list):
    timestring = "%4.4d-%2.2d-%2.2dT%2.2d:%2.2d" % tuple(start[0:5])
  elif isinstance(start,datetime.date):
    # good to go, TBD
    timestring=start.strftime("%Y-%m-%dT%H:%M")
  else:
    # is a string, reparse, TBD
    timestring=start
  
  return(timestring)
  
  
def sm_DateToYMDHMS(tval,yr,mo,dy,hr,mt,sc):
  # not used but here as an example of date conversion
  julday=(tval/86400.0)+2440587.5
  datestr=pd.to_datetime(julday,unit='D',origin='julian')  # format YYYY-MM-DD HH:MM:SS.ssssss
  return(datestr)

def sm_keycheck_data(flagstring):
  # internal helper function
  toggles=['mlt','mag','geo','decl','sza','delta=start','baseline=yearly','baseline=none']

  myflags=''
  flags=[x.strip() for x in flagstring.split(',')]

  for i in range(0,len(flags)):
    chk=flags[i]
    chk=chk.lower()
    # check for the '*all', also individual keys, and assemble url flags
    if chk == 'all': myflags += '&mlt&mag&geo&decl&sza'
    for ikey in range(0,len(toggles)):
      if chk == toggles[ikey]: myflags += '&'+toggles[ikey]

  return(myflags)


def sm_keycheck_indices(flagstring):
  # internal helper function
  # For category='indices', always returns:
  #        tval
  # additional flags to return data include:
  #        indicesall (or its alias: all)
  #  (or any of)
  #        baseall, sunall, darkall, regionalall, plusall
  #  (or specify individual items to include, from the sets below)
  #        
  basekeys=["sme","sml","smu","mlat","mlt","glat","glon","stid","num"]
  # sunkeys: alias allowed of SUN___ -> ___s
  sunkeys=["smes","smls","smus","mlats","mlts","glats","glons","stids","nums"]
  # darkkeys: alias allowed of DARK___ -> ___d
  darkkeys=["smed","smld","smud","mlatd","mltd","glatd","glond","stidd","num"]
  # regkeys: alias allowed of REGIONAL___ -> ___r
  regkeys=["smer","smlr","smur","mlatr","mltr","glatr","glonr","stidr","numr"]
  pluskeys=["smr","ltsmr","ltnum","nsmr"]
  indiceskeys = basekeys + sunkeys + darkkeys + regkeys + pluskeys
  # 'all' means all the above                                                 

  imfkeys=["bgse","bgsm","vgse","vgsm"] # or imfall for all these            
  swikeys=["pdyn","epsilon","newell","clockgse","clockgsm","density"] # % or swiall for all these                                                             
  myflags=''
  indices='&indices='
  swi='&swi='
  imf='&imf='

  flags=[x.strip() for x in flagstring.split(',')]

  for i in range(0,len(flags)):
    chk=flags[i]
    chk=chk.lower()
    
    # check for the '*all', also individual keys, and assemble url flags
    if chk == 'all': indices += 'all,'
    if chk == 'indicesall': indices += 'all,'
    if chk == 'imfall': imf += 'all,'
    if chk == 'swiall': swi += 'all,'
    # available keywords, we allow both the url version and the
    # aliases of "SUN___ -> ___s", "DARK___ -> ___d", "REGIONAL___ -> ___r"

    for ikey in range(0,len(indiceskeys)):
      mykey=indiceskeys[ikey]
      sunkey="sun"+mykey # allow alias
      darkkey="dark"+mykey # allow alias
      regkey1="regional"+mykey # allow alias
      regkey2="reg"+mykey # allow alias
      if chk == mykey:
        indices += mykey+','  # base key is correct
      elif sunkey == mykey:
        indices += mykey+'s,'  # alias, so base key + 's'
      elif darkkey == mykey:
        indices += mykey+'d,'  # alias, so base key + 'd'
      elif regkey1 == mykey or regkey2 == mykey:
        indices += mykey+'r,'  # alias, so base key + 'r'

    for ikey in range(0,len(swikeys)):
      if chk == swikeys[ikey]: swi += swikeys[ikey] + ','
     
    for ikey in range(0,len(imfkeys)):
      if chk == imfkeys[ikey]: imf += imfkeys[ikey] + ','
  
    # more aliases to the user
    if chk == 'baseall': indices += ','.join(basekeys)
    if chk == 'sunall': indices += ','.join(sunkeys)
    if chk == 'darkall': indices += ','.join(darkkeys)
    if chk == 'regionalall' or chk == 'regall': indices += ','.join(regkeys)
    if chk == 'plusall': indices += ','.join(pluskeys)

  # clean it up a bit by removing extraneous tags/characters
  if indices == "&indices=": indices=""
  if swi == "&swi=": swi=""
  if imf == "&imf=": imf=""
  # add them together
  myflags = indices + swi + imf
  # a little more cleaning for tidiness, removes extraneous commas
  myflags = re.sub(',&','&',myflags)
  myflags = re.sub(',$','',myflags)

  return(myflags)

def sm_GetUrl(fetchurl,fetch='raw'):
  # internal helper function
  # returned data choices are 'raw' or 'json', default is 'raw'
  # converts an http bytestream into a python list (raw) or list of dict (json)
  # 'stations' should be 'raw' and returns a list
  # 'data' should be 'json', returns a list (which converts to a dataframe)
  # 'indices' should be 'json', returns a list (which converts to a dataframe)
  
  success = 0 # gets changed to 1 good data is fetched
  longstring=b"ERROR: Unknown error" # prepare for the worst
  
  #print("debug: url trying ",fetch,"is",fetchurl)
  # If the url object throws an error it will be caught here
  try:
    #with urllib.request.urlopen(fetchurl,cafile=certifi.where()) as response:
    # note that 'cafile' is deprecated past python 3.5 but we keep it here
    # to have stronger backward compatability with earlier versions
    try:
      cafile=certifi.where()
    except:
      cafile=''
    with urllib.request.urlopen(fetchurl,cafile=cafile) as response:
      longstring = response.read()

      if fetch == 'json':
        if len(longstring) > 3:
          #mydata = json.loads(longstring[3:]) # skipping initial OK
          mydata = json.loads(longstring)
        else:
          mydata=[] # just the word 'OK', no data, so return no data
        success=1
      else:
        # default is raw byte strings, which we split into a list
        mydata = (longstring.decode('UTF-8')).split('\n')
        success = 1 # it worked
        if re.search(r'ERROR',mydata[0]): success=0 # legit return, but of an err

  except urllib.error.URLError as e:
    #print, !ERROR_STATE.msg
    mydata=['ERROR:HTTP error',e.reason]
  except:
    longstring = longstring.decode('UTF-8')
    mydata=[longstring] # catch-all if nothing below works
    
  #print("debug: function return type is:",type(mydata),".")
  # returns a list for 'raw' or a list of dictionaries for 'json'
  return(success,mydata)


# Gets a list of stations.
# Return value is either '1' plus list of stations, or
#                        '0' plus a string with the error message
# Sample usage:
# 'extent' is how long a window, in seconds, to grab. (86400 sec = 1 day)
# (status, stations) =supermaggetinventory('myname',2019,11,2,20,24,00,86400)
# In this case, 'status'=1 and 'stations' is a list of 184 stations

def SuperMAGGetInventory(logon,start,extent):
  # One of the core 3 functions
  
  iarr=""
  errstr=""

  # construct URL       
  urlstr = sm_coreurl('inventory.php',logon,start,extent)

  # get the string array of stations
  (success,stations)=sm_GetUrl(urlstr,'raw')
  
  # if the inventory is valid extract the stations to an array
  # if an error occurs set the the ERROR keyword to be the error string

  if success:
    # first data item is how many stations were found
    numstations = int(stations[0])
    #print("debug: found",numstations,"stations")
    if numstations > 0: stations=stations[1:-1] # remove OK,#stations,''
    else: stations=[]  # empty list because no stations found
    #print('debug, got back:',stations)
  # success: return true (1) if the call was successful otherwise false (0) 
  # stations: return list with true/1 plus data, or false/0 plus errstr
  return(success,stations)


def SuperMAGGetIndices(logon,start,extent,flagstring='',**kwargs):
  # One of the core 3 functions

  urlstr = sm_coreurl('indices.php',logon,start,extent)
  indices = sm_keycheck_indices(flagstring)
  urlstr += indices
  
  # get the string array of JSON data     
  (status,data_list)=sm_GetUrl(urlstr,'json')

  # default is to return a dataframe, but can also return an array
  if (kwargs.get('FORMAT','none')).lower() == 'list':
    return(status,data_list)
  else:
    # default, converts the json 'list of dictionaries' into a dataframe
    data_df = pd.DataFrame(data_list)
    return(status,data_df)


def SuperMAGGetData(logon,start,extent,flagstring,station,**kwargs):
  # One of the core 3 functions

  # optional options for 'data':
  # ALL=&mlt&mag&geo&decl&sza
  # MLT=&mlt,MAG=&mag,GEO=&geo,DECL=&decl,SZA=&sza,
  # DELTA='start',BASELINE='none/yearly'
  # e.g. can pass  MLT=1,MAG=1  and they will be evaluated.  Full set checked: ALL, MLT, MAG, GEO, DECL, SZA, also values for DELTA, BASELINE
  # also arg FORMAT='list', otherwise defaults to FORMAT='dataframe'  NOT YET DONE!!!!

  # default FORMAT='dataframe', alt is FORMAT='list'
  
  urlstr = sm_coreurl('data-api.php',logon,start,extent)
  indices = sm_keycheck_data(flagstring)
  urlstr += indices
  urlstr += '&station='+station.upper()
  
  (status,data_list)=sm_GetUrl(urlstr,'json')

  # default is to return a dataframe, but can also return an array
  if (kwargs.get('FORMAT','none')).lower() == 'list':
    return(status,data_list)
  else:
    # default, converts the json 'list of dictionaries' into a dataframe
    data_df = pd.DataFrame(data_list)
    return(status,data_df)

def sm_grabme(dataf,key,subkey):
  # syntactical sugar to grab nested subitems from a dataframe
  data = dataf[key]
  subdata = [temp[subkey] for temp in data]
  return(subdata)
  
# Unlike IDL, which returns as Array or Struct,
# we return as List (of dictionaries) or DataFrame

def sm_microtest(choice,userid):
  # 3 simple unit tests to verify the core fetches work
  import matplotlib.pyplot as plt

  start=[2019,11,15,10,40,00] # alt: start='2019-11-15T10:40'

  if choice == 1 or choice == 4:
    (status,stations) = SuperMAGGetInventory(userid,start,3600)
    print(status)
    print(stations)

  if choice == 2 or choice == 4:
    (status,data) = SuperMAGGetData(userid,start,3600,'all,delta=start,baseline=yearly','HBK')
    print(status)
    print(data)
    print(data.keys())
    
    tval=data.tval
    mlt=data.mlt
    ### Python way
    N_nez = [temp['nez'] for temp in data.N]
    N_geo = [temp['geo'] for temp in data.N]
    ### or, supermag helper shorthand way
    N_nez = sm_grabme(data,'N','nez')
    N_geo = sm_grabme(data,'N','geo')
    #
    plt.plot(tval,N_nez)
    plt.plot(tval,N_geo)
    plt.ylabel('N_geo vs N_nez')
    plt.xlabel('date')
    plt.show()

  if choice == 3 or choice == 4:
    (status,idxdata) = SuperMAGGetIndices(userid,start,3600,'swiall,density,darkall,regall,smes')
    #print(status)
    #print(idxdata)
    idxdata.keys()
    tval=idxdata.tval
    hours=list(range(24))
    y=idxdata.SMLr
    for i in range(len(tval)-1):
      plt.plot( hours, y[i] )
      plt.ylabel('SMLr')
      plt.xlabel('hour')
      plt.title('SMLr variation by hour, for successive days')
    plt.show()

def supermag_testing(userid):

  start=[2019,11,15,10,40,00] # alt: start='2019-11-15T10:40'

  (status,stations) = SuperMAGGetInventory(userid,start,3600)


  # DATA fetches
  # BARE CALL, dataframe returned
  (status,mydata1a) = SuperMAGGetData(userid,start,3600,'','HBK')
  mydata1a        # is 1440 rows x 6 columns dataframe
  mydata1a.keys() # Index(['tval', 'ext', 'iaga', 'N', 'E', 'Z'], dtype='object')

  # CALL with ALLINDICES, dataframe returned
  (status,mydata1a) = SuperMAGGetData(userid,start,3600,'all','HBK')
  mydata1a        # is 1440 rows x 12 columns dataframe
  mydata1a.keys() # Index(['tval', 'ext', 'iaga', 'glon', 'glat', 'mlt', 'mcolat', 'decl', 'sza', 'N', 'E', 'Z'], dtype='object')

  # BARE CALL, list returned
  (status,mydata1b) = SuperMAGGetData(userid,start,3600,'','HBK',FORMAT='list')
  len(mydata1b)  # is 1440 rows of dicts (key-value pairs)
  mydata1b[0:1]  # {'tval': 1572726240.0, 'ext': 60.0, 'iaga': 'DOB', 'N': {'nez': -3.942651, 'geo': -5.964826}, 'E': {'nez': 4.492887, 'geo': 0.389075}, 'Z': {'nez': 7.608168, 'geo': 7.608168}}

  # CALL with ALLINDICES, list returned
  (status,mydata1b) = SuperMAGGetData(userid,start,3600,'all','HBK',FORMAT='list')
  mydata1b        # is 1440 rows of dicts (key-value pairs)
  mydata1b[0:1]  # {'tval': 1572726240.0, 'ext': 60.0, 'iaga': 'DOB', 'glon': 9.11, 'glat': 62.07, 'mlt': 21.694675, 'mcolat': 30.361519, 'decl': 3.067929, 'sza': 124.698227, 'N': {'nez': -3.942651, 'geo': -5.964826}, 'E': {'nez': 4.492887, 'geo': 0.389075}, 'Z': {'nez': 7.608168, 'geo': 7.608168}}
  
  ####################
  # INDICES fetches
  (status,idxdata) = SuperMAGGetIndices(userid,start,3600)
  idxdata  # empty!

  (status,idxdata) = SuperMAGGetIndices(userid,start,3600,'all,swiall,imfall')
  idxdata  # 1440 rows x 77 columns dataframe
  idxdata.keys() # Index(['tval', 'SME', 'SML', 'SMLmlat', 'SMLmlt', 'SMLglat', 'SMLglon', 'SMLstid', 'SMU', 'SMUmlat', 'SMUmlt', 'SMUglat', 'SMUglon', 'SMUstid', 'SMEnum', 'SMEs', 'SMLs', 'SMLsmlat', 'SMLsmlt', 'SMLsglat', 'SMLsglon', 'SMLsstid', 'SMUs', 'SMUsmlat', 'SMUsmlt', 'SMUsglat', 'SMUsglon', 'SMUsstid', 'SMEsnum', 'SMEd', 'SMLd', 'SMLdmlat', 'SMLdmlt', 'SMLdglat', 'SMLdglon', 'SMLdstid', 'SMUd', 'SMUdmlat', 'SMUdmlt', 'SMUdglat', 'SMUdglon', 'SMUdstid', 'SMEdnum', 'SMEr', 'SMLr', 'SMLrmlat', 'SMLrmlt', 'SMLrglat', 'SMLrglon', 'SMLrstid', 'SMUr', 'SMUrmlat', 'SMUrmlt', 'SMUrglat', 'SMUrglon', 'SMUrstid', 'SMErnum', 'smr', 'smr00', 'smr06', 'smr12', 'smr18', 'smrnum', 'smrnum00', 'smrnum06', 'smrnum12', 'smrnum18', 'bgse', 'bgsm', 'vgse', 'vgsm', 'clockgse', 'clockgsm', 'density', 'dynpres', 'epsilon', 'newell'], dtype='object')
  #
  # just INDICESALL = 67 columns, above 'tval' through 'smrnum18'
  # just IMFALL = 5 columns, Index(['tval', 'bgse', 'bgsm', 'vgse', 'vgsm'], dtype='object')
  # just SWIALL = 7 columns, Index(['tval', 'clockgse', 'clockgsm', 'density', 'dynpres', 'epsilon', 'newell'], dtype='object')
  #
  # Dataframes are awesome!  To manipulate, just pull out what you need
  import pandas as pd  # call once at the top of your code if you are using dataframes
  tval = idxdata.tval
  density = idxdata.density
  vgse = idxdata.vgse
  # or all as 1 line of code
  tval, density, vgse = idxdata.tval, idxdata.density, idxdata.vgse
  # note that vgse is itself a dictionary of values for X/Y/Z, so you can get subitems from it like this
  vgse_x = [d.get('X') for d in idxdata.vgse]

  # to save the data, there are many formats.  Here is how to save as csv
  idxdata.to_csv('mydata.csv')

  # to read it back in later
  import pandas as pd
  import re
  mydata2b=pd.read_csv('mydata.csv',index_col=0) # you can read it into any variable name, we just used 'mydata2b' as an example
  # now you can do all the above items again, with one exception: each line of the CVS file got split into a dict (key-value pairs) but items like 'vsge' are part of the pandas structure
  # the 'd.get()' approach will _not_ work once read from csv
  stationlist = mydata2b.SMLrstid # item is a pandas series (not python list)
  print(stationlist[0]) # prints a list of stations as a string, but cannot easily access a single item because it is a pandas series
  # so you can convert a pandas series to a list
  stationlist2=sm_csvitem_to_list(mydata2b.SMLrstid) # goal is a list of stations
  slist = stationlist2[0] # grabs a list of stations for row 0
  s1 = stationlist2[0][0] # grabs the first station for row 0

  vgse=sm_csvitem_to_dict(mydata2b.vgse) # goal is a dict of coords or other values
  x = vgse[0]['X'] # grab just the 'X' value for the 1st row of data
  vgse_x = [mydat['X'] for mydat in vgse] # grab all the 'X' values as a new list
  vgse_xyz = [(mydat['X'],mydat['Y'],mydat['Z']) for mydat in vgse] # grab all 3

  # We also offer a list format, for users who prefer to work in python lists
  (status,mydata2c) = SuperMAGGetIndices(userid,start,3600,'all,swiall,imfall',FORMAT='list')
  len(mydata2c)  # is 1440 rows of dicts (key-value pairs)
  mydata2c[0:1] # {'tval': 1572726240.0, 'SME': 58.887299, 'SML': -27.709055, 'SMLmlat': 73.529922, 'SMLmlt': 23.321493, 'SMLglat': 76.510002, 'SMLglon': 25.01, 'SMLstid': 'HOP', 'SMU': 31.178246, 'SMUmlat': 74.702339, 'SMUmlt': 2.090216, 'SMUglat': 79.480003, 'SMUglon': 76.980003, 'SMUstid': 'VIZ', 'SMEnum': 118, 'SMEs': 34.451469, 'SMLs': -16.599854, 'SMLsmlat': 62.368008, 'SMLsmlt': 9.399416, 'SMLsglat': 62.299999, 'SMLsglon': 209.800003, 'SMLsstid': 'T39', 'SMUs': 17.851616, 'SMUsmlat': 73.989975, 'SMUsmlt': 18.228394, 'SMUsglat': 67.93, 'SMUsglon': 306.429993, 'SMUsstid': 'ATU', 'SMEsnum': 54, 'SMEd': 58.887299, 'SMLd': -27.709055, 'SMLdmlat': 73.529922, 'SMLdmlt': 23.321493, 'SMLdglat': 76.510002, 'SMLdglon': 25.01, 'SMLdstid': 'HOP', 'SMUd': 31.178246, 'SMUdmlat': 74.702339, 'SMUdmlt': 2.090216, 'SMUdglat': 79.480003, 'SMUdglon': 76.980003, 'SMUdstid': 'VIZ', 'SMEdnum': 64, 'SMEr': [29.685059, 29.857538, 31.387127, 41.707573, 10.320444, 10.885443, 9.604616, 13.479583, 15.471248, 15.471248, 15.714731, 5.434914, 12.13654, 11.156847, 9.62884, 14.752592, 14.752592, 24.204388, 21.41181, 21.41181, 27.121195, 46.345322, 51.403328, 51.403328], 'SMLr': [-27.709055, 1.320708, -0.208882, -10.529325, -10.529325, -10.529325, -9.248499, -13.123466, -16.599854, -16.599854, -16.599854, -5.449972, -5.449972, -4.470279, -2.942272, -6.352773, -6.352773, -6.352773, -3.560194, -3.560194, -7.514064, -22.651047, -27.709055, -27.709055], 'SMLrmlat': [73.529922, 51.264774, 47.791527, 66.696564, 66.696564, 66.696564, 41.771515, 70.602707, 62.368008, 62.368008, 62.368008, 67.471809, 67.471809, 60.639145, 68.500282, 72.20977, 72.20977, 72.20977, 75.762718, 75.762718, 77.33667, 71.889503, 73.529922, 73.529922], 'SMLrmlt': [23.321493, 2.119074, 3.578985, 4.929673, 4.929673, 4.929673, 5.414416, 8.57761, 9.399416, 9.399416, 9.399416, 11.35623, 11.35623, 12.266475, 13.977451, 16.720993, 16.720993, 16.720993, 19.65963, 19.65963, 21.307804, 22.863134, 23.321493, 23.321493], 'SMLrglat': [76.510002, 55.029999, 52.169998, 71.580002, 71.580002, 71.580002, 47.799999, 71.300003, 62.299999, 62.299999, 62.299999, 61.756001, 61.756001, 53.351002, 58.763, 63.75, 63.75, 63.75, 72.300003, 72.300003, 76.769997, 74.5, 76.510002, 76.510002], 'SMLrglon': [25.01, 82.900002, 104.449997, 129.0, 129.0, 129.0, 132.414001, 203.25, 209.800003, 209.800003, 209.800003, 238.770004, 238.770004, 247.026001, 265.920013, 291.480011, 291.480011, 291.480011, 321.700012, 321.700012, 341.369995, 19.200001, 25.01, 25.01], 'SMLrstid': ['HOP', 'NVS', 'IRT', 'TIK', 'TIK', 'TIK', 'BRN', 'BRW', 'T39', 'T39', 'T39', 'FSP', 'FSP', 'C06', 'FCC', 'IQA', 'IQA', 'IQA', 'SUM', 'SUM', 'DMH', 'BJN', 'HOP', 'HOP'], 'SMUr': [1.976003, 31.178246, 31.178246, 31.178246, -0.208882, 0.356117, 0.356117, 0.356117, -1.128606, -1.128606, -0.885122, -0.015059, 6.686568, 6.686568, 6.686568, 8.399819, 8.399819, 17.851616, 17.851616, 17.851616, 19.60713, 23.694275, 23.694275, 23.694275], 'SMUrmlat': [52.904049, 74.702339, 74.702339, 74.702339, 47.791527, 54.29908, 54.29908, 54.29908, 66.244217, 66.244217, 57.76614, 54.597057, 55.715378, 55.715378, 55.715378, 57.829525, 57.829525, 73.989975, 73.989975, 73.989975, 70.473801, 68.194489, 68.194489, 68.194489], 'SMUrmlt': [0.510692, 2.090216, 2.090216, 2.090216, 3.578985, 6.394085, 6.394085, 6.394085, 9.99274, 9.99274, 11.729218, 12.269058, 13.969843, 13.969843, 13.969843, 16.160952, 16.160952, 18.228394, 18.228394, 18.228394, 21.200783, 22.967857, 22.967857, 22.967857], 'SMUrglat': [56.432999, 79.480003, 79.480003, 79.480003, 52.169998, 59.970001, 59.970001, 59.970001, 64.047997, 64.047997, 51.882999, 47.664001, 45.870998, 45.870998, 45.870998, 48.650002, 48.650002, 67.93, 67.93, 67.93, 70.900002, 71.089996, 71.089996, 71.089996], 'SMUrglon': [58.567001, 76.980003, 76.980003, 76.980003, 104.449997, 150.860001, 150.860001, 150.860001, 220.889999, 220.889999, 239.973999, 245.791, 264.916992, 264.916992, 264.916992, 287.549988, 287.549988, 306.429993, 306.429993, 306.429993, 351.299988, 25.790001, 25.790001, 25.790001], 'SMUrstid': ['ARS', 'VIZ', 'VIZ', 'VIZ', 'IRT', 'MGD', 'MGD', 'MGD', 'DAW', 'DAW', 'C13', 'C10', 'C08', 'C08', 'C08', 'T50', 'T50', 'ATU', 'ATU', 'ATU', 'JAN', 'NOR', 'NOR', 'NOR'], 'SMErnum': [5, 3, 3, 4, 5, 6, 6, 4, 8, 9, 12, 13, 20, 17, 17, 11, 12, 14, 12, 14, 22, 51, 51, 35], 'smr': 0.252399, 'smr00': -0.531382, 'smr06': 0.885406, 'smr12': 1.051192, 'smr18': -0.395618, 'smrnum': 72, 'smrnum00': 26, 'smrnum06': 23, 'smrnum12': 6, 'smrnum18': 17, 'bgse': {'X': 1.07, 'Y': -3.75, 'Z': -0.74}, 'bgsm': {'X': 1.07, 'Y': -3.82, 'Z': -0.06}, 'vgse': {'X': -351.100006, 'Y': -5.5, 'Z': -4.0}, 'vgsm': {'X': 351.100006, 'Y': 6.128625, 'Z': -2.947879}, 'clockgse': 258.340698, 'clockgsm': 268.664337, 'density': 5.03, 'dynpres': 1.25, 'epsilon': 29.468521, 'newell': 2504.155029}
  # sample accessing
  print(mydata2c[0]['tval'],mydata2c[0]['density'])  # single element
  result=[ (myeach['tval'],myeach['density']) for myeach in mydata2c] # pull out pairs e.g. 'tval, density')
  # two-line method for extracting any variable set from this
  pairsets= [ (myeach['tval'],myeach['density'],myeach['vgse']) for myeach in mydata2c] # same, pull out pairs, only assign e.g. x=tval, y=density
  tval, density, vgse = [ [z[i] for z in pairsets] for i in (0,1,2)]
  # since 'vgse' is itself an dict of 3 values X/Y/Z, you can pull out nested items like this
  pairsets= [ (myeach['tval'],myeach['density'],myeach['vgse']['X']) for myeach in mydata2c] # same, pull out pairs, only assign e.g. x=tval, y=density
  tval, density, vgse_x = [ [z[i] for z in pairsets] for i in (0,1,2)]
  # the above methods are extensible to any number of variables, just update the (0,1,2) to reflect now many you have
  

#  Uncomment to run quick sample tests
# userid=YOUR_SUPERMAG_USER_ID
#sm_microtest(1,userid)   # sample stations fetch
#sm_microtest(2,userid)   # sample data fetch, with plotting
#sm_microtest(3,userid)   # sample indices fetch, with plotting
