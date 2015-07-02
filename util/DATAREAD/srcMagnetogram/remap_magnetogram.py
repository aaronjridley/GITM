#!/usr/bin/env python

#this magnetogram remapping can either be run as a script from the unix command line
#  or imported into Python
#accompanying files (must be in same directory):
#    remap_magnetogram.py.README.txt
#         by  Richard A. Frazin July 2014 - February 2015

from astropy.io import fits
from scipy import interpolate
from scipy import integrate
import numpy as np
import sys
import time
import argparse
import pdb


def remap(inputfile, outputfile, nlat = -1, nlong = -1, out_grid = 'unspecified'):
    """
    Flux-conserving magnetogram remapping tool.
    inputfile - FITS file containing original magnetogram (include path)
    outputfile - contains result in customized output format (include path)
    nlat  (opitonal) = desired number of latitude  points in output,
         if not specified, will be set to the same as the input file
    nlong (optional) = desired number of longitude points in output
        if not specified, will be set to the same as the input file
    out_grid (optional), choices are 'sin(lat)' or 'uniform'
        if not specified, the output grid will be the same type
    If nlat, nlong and out_grid are ALL left out, no remapping is done,
        and the code simply reformats.
    Note that ADAPT files may have multiple maps.  Only the 1st is utilized.
    My MATLAB code remap_mag.m uses the same algorithm but runs much faster.
       by Richard Frazin, July 2014 - Feb 2015
    """
    pi = 3.141592653589793

    if ( (out_grid != 'sin(lat)') and (out_grid != 'uniform') and (out_grid != 'unspecified') ):
        print "Unknown output grid type.  Choices are blank, 'unspecified', 'uniform' and 'sin(lat)' "
        return(-1)
    
    cc =  FITS_RECOGNIZE(inputfile)
    if cc == -1:
        print "Input file not recognized."
        return(-1)
    else:
        magtype = cc[0]
        grid_type = cc[1]

    if out_grid == 'unspecified':
        out_grid = grid_type

    #what kind of transformation are we doing?
    if grid_type == out_grid:
        transformation = 'rebin' #no change in grid type, so just rebin
    elif ( (grid_type == 'uniform') and (out_grid == 'sin(lat)') ):
        transformation = 'reg2sin'
    elif ( (grid_type == 'sin(lat)') and (out_grid == 'uniform') ):
        transformation = 'sin2reg'
    else:
        print "Unknown transformation type."
        return(-1)

    g = fits.open(inputfile)
    d = g[0].data
    if magtype == 'ADAPT Synchronic':
        nim = g[0].header['NAXIS3'] # number of images
        imdex = 0  #which of the 12 maps do you want?
        print 'This file contains ', str(nim), ' images.  Using only number ', str(imdex)
        if nim > 1:  #just keep one of them for now
            d = d[imdex,:,:]

    nlo = g[0].header['NAXIS1'] # number of longitude points
    nla = g[0].header['NAXIS2'] #           latitude
    
    if nlat == -1:
        nlat = nla
    elif nlat < 1:
        print "nlat has to be -1 or a positive integer."
        return(-1)

    if nlong == -1:
        nlong = nlo
    elif nlong < 1:
        print "nlong has to be -1 or a positive integer."
        return(-1)
    
    
    g[0].header['NAXIS1'] = nlong # new number of longitude points
    g[0].header['NAXIS2'] = nlat  #               latitude
    
    try:
        g[0].header['GRID'] = out_grid #change the existing header value
    except KeyError, er:
        g[0].header.set('GRID',out_grid) #create FITS header keyword

    try:
        g[0].header['CTYPE2'] = out_grid #change the existing header value
    except KeyError, er:
        g[0].header.set('CTYPE2',out_grid) #create FITS header keyword

    if out_grid == 'sin(lat)':
        newlat = (180/pi)*np.arcsin(np.linspace(-1. + 1./2/nlat,1. - 1./2/nlat,nlat))
    elif out_grid == 'uniform':
        newlat = (180/pi)*np.linspace(-pi/2 + pi/2/nlat,pi/2 - pi/2/nlat,nlat)
    else:
        print "out_grid incorrectly set."
        return(-1)


    if ( (nlo == nlong) and (nla == nlat) and (grid_type == out_grid) ):
        newmap = d  #no remapping
    else:
        #first make a hybrid map that is (nla X nlong) by using the rebin
        #    and add alg. in the longitude direction.  If nlo = nlong, hybrid --> d
        hybrid = np.zeros([nla,nlong]) 
        crap = np.arange(nlo+1)
        for pf in crap[nlo+1:0:-1]: #pf will be the greatest common factor of nlong and nlo
            if ( (np.mod(nlo,pf) == 0) and (np.mod(nlong,pf) == 0)): #common factor test
                nlo_fac   = nlo/pf
                nlong_fac = nlong/pf
                break
        for k in np.arange(nla):
            w = np.kron(d[k,:],np.ones(nlong_fac)) #this array has length pf*nlo_fac*nlong_fac
            for l in np.arange(nlong):  #take the average over nlo_fac bins of w
                hybrid[k,l] = np.sum(w[l*nlo_fac:(l+1)*nlo_fac])/nlo_fac

        newmap = np.zeros([nlat,nlong]) #output map                           
        if transformation == 'rebin':  #do rebin and add in the latitude direction, if nlo = nlat, newmap --> d
            crap = np.arange(nla+1)
            for pf in crap[nla+1:0:-1]: #pf will be the greatest common factor of nla and nlat
                if ( (np.mod(nla,pf) == 0) and (np.mod(nlat,pf) == 0) ): # common factor test
                    nla_fac  = nla/pf
                    nlat_fac = nlat/pf
                    break
                
            for k in np.arange(nlong):
                w = np.kron(hybrid[:,k].T,np.ones(nlat_fac))#length is pf*nla_fac*nlat_fac
                for l in np.arange(nlat):
                    newmap[l,k] = np.sum(w[l*nla_fac:(l+1)*nla_fac])/nla_fac

        elif transformation == 'reg2sin':
            oldlat =  np.linspace(-pi/2 + pi/2/nla,pi/2 - pi/2/nla,nla) #old latitude grid
            oldlat = np.hstack((-pi/2-1.e-9,oldlat,pi/2+1.e-9)) #for the interpolator
            bin_boundary = np.arcsin(np.linspace(-1.,1.,nlat+1))   #boundaries of new sin(latitude) grid
            for k in np.arange(nlong):   #the magnetic field value assigned is the flux divided by the area.  
                u = np.hstack((hybrid[0,k],hybrid[:,k],hybrid[nla-1,k]))
                crap = interpolate.interp1d(oldlat,u,kind='linear') #magnetic field interpolator
                fcn = lambda x : crap(x)*np.cos(x)  #this is B(theta)*cos(theta)
                for l in np.arange(nlat):
                    result = integrate.quad(fcn,bin_boundary[l],bin_boundary[l+1],epsabs=1.e-3,epsrel=1.e-3)/(np.sin(bin_boundary[l+1]) - np.sin(bin_boundary[l]))
                    newmap[l,k] = result[0]

        elif transformation == 'sin2reg':
            oldlat = np.arcsin(np.linspace(-1. + 1./2/nla,1. - 1./2/nla,nla)) #arcsin(old sin(latitude) grid)
            oldlat = np.hstack((-pi/2-1.e-9,oldlat,pi/2+1.e-9)) #for the interpolator
            bin_boundary = np.linspace(-pi/2,pi/2,nlat+1) #boundaries of new latitude grid
            #pdb.set_trace()
            for k in np.arange(nlong):   #the magnetic field value assigned is the flux divided by the area.  
                u = np.hstack((hybrid[0,k],hybrid[:,k],hybrid[nla-1,k]))
                crap = interpolate.interp1d(oldlat,u,kind='linear') #magnetic field interpolator
                fcn = lambda x : crap(x)*np.cos(x)  #this is B(theta)*cos(theta)
                for l in np.arange(nlat):
                    result = integrate.quad(fcn,bin_boundary[l],bin_boundary[l+1],epsabs=1.e-3,epsrel=1.e-3)/(np.sin(bin_boundary[l+1]) - np.sin(bin_boundary[l]))
                    newmap[l,k] = result[0]

        else:
            print "Unknown transformation type."
            return(-1)


    #test for flux conservation in the transformation        
    test_flux = False 
    if test_flux:
        if grid_type == 'uniform':
            latt =  np.cos(np.linspace(-pi/2 + pi/2/nla,pi/2 - pi/2/nla,nla))
            cosgrid = np.kron(latt,np.ones((nlo,1))).T
            oldflux = np.sum(np.multiply(cosgrid,d))*2.*pi*pi/nlo/nla
        elif grid_type == 'sin(lat)':
            oldflux = np.sum(d)*4.*pi/nlo/nla
        else:
            print "Bad grid_type."
            return(-1)
        if out_grid == 'uniform':
            latt =  np.cos(np.linspace(-pi/2 + pi/2/nlat,pi/2 - pi/2/nlat,nlat))
            cosgrid = np.kron(latt,np.ones((nlong,1))).T
            newflux = np.sum(np.multiply(cosgrid,newmap))*2.*pi*pi/nlong/nlat
        elif out_grid == 'sin(lat)':
            newflux = np.sum(newmap)*4.*pi/nlong/nlat
        else:
            print "Bad out_grid."
            return(-1)
        print "original flux =",str(oldflux),", new flux =",str(newflux)
    

    #try to get some context information from the FITS file to include in output
    try:
        CRnumber = str(g[0].header['CRROTEDG'])  #works for ADAPT
    except KeyError, er:
        CRnumber = '0'

    if CRnumber == '0':    
        try :
            CRnumber = str(g[0].header['CAR_ROT']) #works on GONG and MDI
        except KeyError,er:
            CRnumber = '0'

    if magtype.find('GONG') > -1:
        mapdate = g[0].header['DATE'] #works for GONG
    else:                          
        try:
            mapdate = g[0].header['MAPTIME']  #works for ADAPT
        except KeyError, er:
            mapdate = '0000-00-00T00:00:00'

        if mapdate == '0000-00-00T00:00:00':    
            try :
                mapdate = g[0].header['MAP_DATE']  #works for Hathaway
            except KeyError, er:
                mapdate = '0000-00-00T00:00:00'

        if mapdate == '0000-00-00T00:00:00':    
            try :
                mapdate = g[0].header['T_OBS']  #works for MDI
            except KeyError, er:
                mapdate = '0000-00-00T00:00:00'
    try:
        bunit = g[0].header['BUNIT']  #works on GONG, MDI
    except KeyError, er:   #Hathaway and ADAPT don't list units
        bunit = 'Gauss'  #assume it's Gauss if you don't know
                
    try :
        long0 = g[0].header['LONG0'] #works on GONG 
    except KeyError, er:
        long0 = 0

    #ascii output file, Gabor format, the first line is arbitary
    fid = open(outputfile,'w')
    
    line0 = 'magnetogram type = '+magtype+', grid_type = sin(lat), CR'+CRnumber+ ', MapDate = '+mapdate+', units: ['+bunit+'], created at: '+time.ctime()+'\n' 
    fid.write(line0)
    line0 = '       0      0.00000       2       1       1 \n'
    fid.write(line0)
    fid.write('      '+str(nlong)+'     '+str(nlat)+'\n')
    fid.write(str(long0 + 0.5*360./nlong) + ' \n') #longitude shift (important for GONG Hourly)
    fid.write('Longitude Latitude Br LongitudeShift \n')
    
    for k in np.arange(nlat):
         for l in np.arange(nlong):
             line0 = str(l*360./nlong) + ' ' + str(newlat[k]) + ' ' + str(newmap[k,l]) + ' \n'
             fid.write(line0)

    #old output format
    ## fid.write('#CR\n')
    ## try:
    ##     fid.write(str(g[0].header['CRROTEDG']) + '\n')
    ## except KeyError, er:
    ##     fid.write('-1')
    ## fid.write('#nMax\n')
    ## fid.write(str(-1) + '\n')
    ## fid.write('#ARRAYSIZE\n')
    ## fid.write(str(nlong) + '\n')
    ## fid.write(str(nlat) + '\n')
    ## fid.write('#START\n')
    ## for k in np.arange(nlat):
    ##     for l in np.arange(nlong):
    ##         fid.write(str(newmap[k,l]) + '\n')

    g.close()
    fid.close()
    return(newmap,d)
    

def FITS_RECOGNIZE(inputfile):
    """
    This function opens inputfile and tries to determine what type of magnetogram
    it is as well as the type of grid on which the datatype is represented.  The
    magnetogram types and grid types are: 
      Hathaway Synchronic, regular
      ADAPT Synchronic, regular 
      GONG Synoptic, sin(lat)
      GONG Hourly updated, sin(lat)
      MDI Synoptic, sin(lat)
    This function returns a tuple with this information.  The output tuple is:
      (magnetogram_type,grid_type)
    """

    magnetogram_type = 'unknown'
    grid_type = 'unknown'
    g = fits.open(inputfile)

    try:
        telescope = g[0].header['TELESCOP'] #works for MDI, GONG
    except KeyError, er:
        telescope = 'unknown'
        
    try:
        inst = g[0].header['INSTRUME'] #works for MDI
    except KeyError,er:
        inst = 'unknown'

    try:
        ctyp = g[0].header['CTYPE2'] #works for MDI, GONG
    except KeyError, er:
        ctyp = 'unknown'

    try:
        model = g[0].header['MODEL'] #works for ADAPT
    except KeyError, err:
        model = 'unknown'

    try:
        sft = g[0].header['SFT_TYP'] #works for Hathaway
    except KeyError,er:
        sft = 'unknown'

    nlo = g[0].header['NAXIS1'] # number of longitude points
    nla = g[0].header['NAXIS2'] #           latitude
        
    g.close()

    if telescope.find('NSO-GONG') > -1 :
        magnetogram_type = 'NSO-GONG Synoptic'
        try:
            long0 = g[0].header['LONG0']
            if float(long0) > 0.:
                magnetogram_type = 'NSO-GONG Hourly'
        except KeyError, er:
            long0 = - 1
        if ctyp.find('CRLT-CEA') > -1:
            grid_type = 'sin(lat)'
        else:
            print "unknown NSO-GONG magnetogram type"
            return(-1)

    if telescope.find('SOHO') > -1:
        if ( (inst.find('MDI') > -1) & (ctyp.find('Sine Latitude') > -1) ):
            magnetogram_type = 'MDI Synoptic'
            grid_type = 'sin(lat)'
        else :
            print "unknown SOHO magnetogram type"
            return(-1)

    if model.find('ADAPT') > -1:
        magnetogram_type = 'ADAPT Synchronic'
        try:
            adapt_grid = g[0].header['GRID']
        except KeyError, er:
            adapt_grid = -1.
        if adapt_grid == 1.:
            grid_type = 'uniform'
        else:
            print "unknown ADAPT magnetogram type"
            return(-1)

    if sft.find('Baseline / Assimilation') > -1:
        magnetogram_type = 'Hathaway Synchronic'
        grid_type = 'uniform'

    if  ( (magnetogram_type == 'unknown') or (grid_type == 'unknown') ):
        print "I don't recognize the type of this magnetogram."
        return(-1)
                
    print "I think this is a",magnetogram_type,"magnetogram on a",str(nla),"X",str(nlo),grid_type,"grid."
    return( (magnetogram_type, grid_type) )

    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="""
    remap_magnetogram.py pre-processes the FITS format magnetograms 
    into ASCII files that can by read by FDIPS.exe, BATSRUS.exe and SWMF.exe
    and IDL macros. The script can read the following types of magnetograms:

       Hathaway Synchronic
       ADAPT Synchronic
       GONG Synoptic
       GONG Hourly updated
       MDI Synoptic

    The code opens the .fits file and automatically recognizes the type of
    map it is, which determines whether it is on a sin(latitude) or regular
    spherical grid.  The output can be any desired resolution, on either a
    sin(latitude) or regular spherical grid.  If the output grid type is not
    specified, it will be the same as the original .fits file.  If the
    resolution is not specified, it will be the same as the original .fits
    file. The calling syntax from the command line is shown above. Some examples:

    ./remap_magnetogram.py test.fits test.out
    ./remap_magnetogram.py test.fits test.out 180 360
    ./remap_magnetogram.py test.fits test.out -grid=uniform
 
    Within Python, the remapping is done with the remap function contained
    in this file.
    
    The script uses the scipy and astropy packages that can be installed, 
    for example, with MacPorts.
    """)
    parser.add_argument('inputfile', help='Input FITS file name including path')
    parser.add_argument('outputfile', help='Output magnetogram file name including path')
    parser.add_argument('nlat', nargs='?', type=int, default=-1, help='Number of latitude points in output. Default is same as input.')
    parser.add_argument('nlon', nargs='?', type=int, default=-1, help='Number of longitude points in output. Default is same as input.')
    parser.add_argument('-grid',choices=['uniform','sinlat'],help="type of latitude grid in the output. Default is same as input.")

    args = parser.parse_args()

    if args.nlat < -1:
        print "nlat must be -1 or a postive integer.  No output."
        quit()
    if args.nlon < -1:
        print "nlon must be -1 or a postive integer.  No output."
        quit()

    grid_type = 'unspecified'
    if args.grid == 'sinlat':
        grid_type = 'sin(lat)'
    elif args.grid == 'uniform':
        grid_type = 'uniform'

    remap(args.inputfile, args.outputfile, args.nlat, args.nlon, grid_type )

    

        

    
