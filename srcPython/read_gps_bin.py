#!/usr/bin/env python
'''
Read TEC binary file that is converted from the Madrigal TEC files
'''

import datetime as dt
from struct import unpack
import numpy as np

nBytesLong = 4
nBytesFloat = 4

class GpsFile:

    def __init__(self, file):
        self.FileName = file
        self._read_header()

    def _read_header(self):

        f = open(self.FileName,'rb')
        self.nTimes = unpack('<l',f.read(nBytesLong))[0]
        #print(self.nTimes)

        self.time=[]
        self.offsets=[]
        self.nPts=[]
        for i in range(self.nTimes):
            self.offsets.append(unpack('<l',f.read(nBytesLong))[0])
            (yy,mm,dd,hh,mn,ss)=unpack('<llllll',f.read(nBytesLong*6))
            self.time.append(dt.datetime(yy,mm,dd,hh,mn,ss))
            if (i>1):
                self.nPts.append(self.offsets[i]-self.offsets[i-1])
            #print(i)
        LastPoint = unpack('<l',f.read(nBytesLong))[0]
        self.nPts.append(LastPoint-self.offsets[self.nTimes-1])

        nBytesHeader = \
            1*nBytesLong + \
            (1+6)*nBytesLong*self.nTimes + \
            1*nBytesLong
        nBytesPerLine = 4*nBytesFloat

        for i in range(self.nTimes):
            self.offsets[i] = self.offsets[i]*nBytesPerLine + \
                nBytesHeader + i*(nBytesLong+6*nBytesLong)

        f.close()

    def read_time(self, InTime):

        f = open(self.FileName,'rb')

        iPt = 0
        if (InTime > self.time[self.nTimes-1]):
            iPt = self.nTimes-1
        else:
            while (self.time[iPt] < InTime):
                iPt=iPt+1

        print("iPt to seek : ",iPt)

        f.seek(self.offsets[iPt])
        iPtRead = unpack('<l',f.read(nBytesLong))[0]
        (yy,mm,dd,hh,mn,ss)=unpack('<llllll',f.read(nBytesLong*6))
        print("iPt read  : ",iPtRead)
        print("Time read : ",yy,mm,dd,hh,mn,ss)

        self.lon  = [0]*self.nPts[iPt]
        self.lat  = [0]*self.nPts[iPt]
        self.Tec  = [0]*self.nPts[iPt]
        self.dTec = [0]*self.nPts[iPt]

        iCount=0;
        for i in range(self.nPts[iPt]):
            self.lon[i]  = unpack('<f',f.read(nBytesFloat))
            self.lat[i]  = unpack('<f',f.read(nBytesFloat))
            self.Tec[i]  = unpack('<f',f.read(nBytesFloat))
            self.dTec[i] = unpack('<f',f.read(nBytesFloat))
            iCount=iCount+1

        print(iCount," times read. Expected ",self.nPts[iPt],".")

        f.close()

        # Recast output as a 1D numpy array
        self.lon = np.array(self.lon).flatten()
        self.lat = np.array(self.lat).flatten()
        self.Tec = np.array(self.Tec).flatten()
        self.dTec = np.array(self.dTec).flatten()





