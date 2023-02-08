#!/usr/bin/env python

from supermag_api import *
import argparse
from datetime import datetime
from datetime import timedelta

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args_sm():

    parser = argparse.ArgumentParser(description = 'Download SuperMAG AE file')

    parser.add_argument('start', metavar = 'start', \
                        help = 'start date as YYYYMMDD[.HHMM]')
    parser.add_argument('end', metavar = 'end', \
                        help = 'end date as YYYYMMDD[.HHMM]')

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
#
# ----------------------------------------------------------------------

def write_sme_file(data, message = "none"):

    ymd = data['times'][0].strftime('%Y%m%d')
    fileout = 'ae_' + ymd + '.txt'
    print('  --> Writing file ' + fileout)
        
    l0 = 'File created by python code using write_sme_file\n'
    l1 = '============================================================\n'
    l2 = '<year>  <month>  <day>  <hour>  <min>  <sec>  '
    l2 = l2 + '<SME (nT)>  <SML (nT)>  <SMU (nT)>\n'

    fp = open(fileout, 'wb')
    fp.write(l0.encode())
    fp.write("\n".encode())
    if (message != "none"):
        m = message + "\n"
        fp.write(m.encode())
    fp.write(l1.encode())
    fp.write(l2.encode())

    for i, t in enumerate(data['times']):
        ae = data['ae'][i]
        al = data['al'][i]
        au = data['au'][i]
        out = " %8.2f %8.2f %8.2f" % (ae, al, au)
        ymdhms = t.strftime('%Y  %m  %d  %H  %M  %S')
        line = ymdhms + out + "\n"
        fp.write(line.encode())

    fp.close()

    return fileout

# ----------------------------------------------------------------------
#
# ----------------------------------------------------------------------

def download_sme_data(start, end):

    startStr = start.strftime('%Y-%m-%dT%M:%D')
    length = (end - start).total_seconds()

    if (length < 60):
        print('Length requested: ', length)
        print('This is too short. Stopping')
        exit()
    
    if (length > 31.0*86400.0):
        print('Length requested: ', length)
        print('This is too long. Stopping')
        exit()
    
        
    print('--> Downloading AE data')
    userid = 'ridley'
    (status,idxdata) = SuperMAGGetIndices(userid, start, length, 'all')

    nPts = len(idxdata['tval'])
    print('  --> Found %d points' % nPts)

    if (nPts > 0):
        
        base = datetime(1970,1,1,0,0,0)
        current = base + timedelta(seconds = idxdata['tval'][0])

        ymd = current.strftime('%Y%m%d')
        fileout = 'ae_' + ymd + '.txt'
        print('  --> Writing file ' + fileout)
        
        l0 = 'File created by python code using SuperMAGGetIndices\n'
        l1 = '============================================================\n'
        l2 = '<year>  <month>  <day>  <hour>  <min>  <sec>  '
        l2 = l2 + '<SME (nT)>  <SML (nT)>  <SMU (nT)>\n'

        fp = open(fileout, 'wb')
        fp.write(l0.encode())
        fp.write("\n".encode())
        fp.write(l1.encode())
        fp.write(l2.encode())

        for i, t in enumerate(idxdata['tval']):
            current = base + timedelta(seconds = t)
            ae = idxdata['SME'][i]
            al = idxdata['SML'][i]
            au = idxdata['SMU'][i]
            out = " %8.2f %8.2f %8.2f" % (ae, al, au)
            ymdhms = current.strftime('%Y  %m  %d  %H  %M  %S')
            line = ymdhms + out + "\n"
            fp.write(line.encode())

        fp.close()

    else:
        
        fileout = 'none.txt'
        print('  --> ERROR!! No AE data downloaded, no file written!!')

    return fileout

#------------------------------------------------------------------------------
# main code:
#------------------------------------------------------------------------------

if __name__ == '__main__':
    
    args = parse_args_sm()

    # ----------------------------------------
    # Set Times:

    start = args.start
    end = args.end

    yStart = start[0:4]
    mo = start[4:6]
    da = start[6:8]

    if (len(start) >= 11):
        hr = start[9:11]
        if (len(start) >= 13):
            mi = start[11:13]
        else:
            mi = '00'
    else:
        hr = '00'
        mi = '00'

    start = datetime(int(yStart), int(mo), int(da), int(hr), int(mi), 0)

    yr = end[0:4]
    mo = end[4:6]
    da = end[6:8]

    if (len(end) >= 11):
        hr = end[9:11]
        if (len(end) >= 13):
            mi = end[11:13]
        else:
            mi = '00'
    else:
        hr = '00'
        mi = '00'

    end = datetime(int(yr), int(mo), int(da), int(hr), int(mi), 0)

    smeFile = download_sme_data(start, end)
