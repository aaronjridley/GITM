!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine init_iri

  use EUA_ModIri90, ONLY: iri90

  use ModGITM
  use ModInputs
  use ModConstants
  use ModPlanet
  use ModTime

  implicit none

  ! iri variables

  integer :: jmag,nzkm
  real, dimension(1:11, 1:nAlts) :: outf
  real, dimension(1:30) :: oarr
  logical, dimension(1:12) :: jf
  !  character, dimension(1:80) :: ccirnm, ursinm

  integer :: iBlock, iAlt, iLat, iLon, iIon
  real :: geo_lat, geo_lon, geo_alt(1), geo_lst

  call report("init_iri",0)

  !--------------------------------------------------------------------------
  !
  !  From the iri90 library:
  !
  !          MODEL INPUTS (see also comments in iri90.f)
  !            JF    = Option Flags .TRUE. (or .FALSE.)
  !                    JF(1) = electron density is (not) calculated
  !                    JF(2) = temperatures are (not) calculated
  !                    JF(3) = ion composition is (not) calculated
  !                    JF(4) = B0 from table (from Gulyeava 1987)
  !                    JF(5) = F2 peak from CCIR (from URSI)
  !                    JF(6) = ion composition standard(Danilov-Yaichnikov-1985
  !                    JF(7) = standard IRI topside (IRI-79)
  !                    JF(8) = NmF2 peak model (input value).
  !                          = .F. => set OARR(1) = Nmax(m-3) or FoF2 (MHz).
  !                    JF(9) = HmF2 peak model (input value).
  !                          = .F. => set OARR(2) = Hmax (km)
  !                    JF(10)= Te model (Te-Ne model with Ne input).
  !                          = .F. => set OARR(3) = Ne (m-3) at 300 km
  !                                and/or OARR(4) = Ne (m-3) at 400 km
  !                                and/or OARR(5) = Ne (m-3) at 600 km
  !                          (If only specifying Ne 1 ht, set the others to 0)
  !                    JF(11)= Ne standard (LAY-functions version)
  !                    JF(12)= message are written to unit=12 (=6)
  !
  !            JMAG  = Lat/lon are geoditic (0) or geomagnetic (1)
  !            ALATI = Lat (degress, positive North)
  !            ALONG = Lon (degrees, positive East)
  !            RZ12  = Sunspot number (12-month running mean); or equivalent
  !                    F10.7 Solar Flux entered as a negative.
  !            MMDD  = Month and day; or day number of year entered as negative
  !            DHOUR = Hour of day, LT (or UT + 25) in decimal hours.
  !            ZKM   = heights
  !            NZKM  = no. heights
  !            CCIRNM= name of ccir coefficient file(ccir11.asc ... ccir22.asc)
  !            URSINM= name of ursi coefficient file(ursi11.asc ... ursi22.asc)
  !          MODEL OUTPUTS
  !            OUTF = Profiles
  !                   OUTF(1,*) = electron density (m-3)
  !                   OUTF(2,*) = neutral temperature (K)
  !                   OUTF(3,*) = ion temperature (K)
  !                   OUTF(4,*) = electron temperature (K)
  !                   OUTF(5,*) = percent O+ ions
  !                   OUTF(6,*) = percent H+ ions
  !                   OUTF(7,*) = percent He+ ions
  !                   OUTF(8,*) = percent O2+ ions
  !                   OUTF(9,*) = percent No+ ions
  !                   OUTF(10,*)= percent cluster ions when JF(6) = .F.
  !                   OUTF(11,*)= percent N+ ions      when JF(6) = .F.
  !            OARR = Additional output parameters
  !                   OARR(1)  = NmF2 (m-3)        OARR(2)  = HmF2 (km)
  !                   OARR(3)  = NmF1 (m-3)        OARR(4)  = HmF1 (km)
  !                   OARR(5)  = NmE  (m-3)        OARR(6)  = HmE  (km)
  !                   OARR(7)  = NmD  (m-3)        OARR(8)  = HmD  (km)
  !                   OARR(9)  = Hhalf (km)        OARR(10) = B0   (km)
  !                   OARR(11) = Valley-base (m-3) OARR(12) = Valley-top (km)
  !                   OARR(13) = Te-peak (K)     OARR(14) = Te-peak height (km)
  !                   OARR(15) = Te-mod(300km)   OARR(16) = Te-mod(400km) (K)
  !                   OARR(17) = Te-mod(600km)   OARR(18) = Te-mod(1400km) (K)
  !                   OARR(19) = Te-mod(3000km)  OARR(20) = Te(120km)=Tn=Ti (K)
  !                   OARR(21) = Ti-mod(430km)   OARR(22) = x (KM), where Te=Ti
  !                   OARR(23) = Solar zenith angle (degrees)
  !                   OARR(24) = Solar declination (degrees)
  !                   OARR(25) = Dip
  !                   OARR(26) = Dip latitude
  !                   OARR(27) = modified dip latitude
  !
  !*** THE ALTITUDE LIMITS ARE:  LOWER (DAY/NIGHT)  UPPER        ***
  !***     ELECTRON DENSITY         60/80 KM       1000 KM       ***
  !***     TEMPERATURES              120 KM        3000 KM       ***
  !***     ION DENSITIES             100 KM        1000 KM       ***
  !--------------------------------------------------------------------------
  !  Set recommended defaults for logicals (except set JF(6)=.F. to get N+)
  jf(:) = .true.
  jf(4) = .false.
  jf(5) = .false.
  jf(6) = .false.
  ! jf(12) = .false. means all print-out goes to 6, the console, 
  !          (instead of 12)
  !  jf(12) = .false.
  !  Assume lat/lon are geodetic (jmag=0)
  jmag = 0
  !  Get only 1 height at a time, so nzkm = 1
  nzkm = 1
  !  Set names of CCIR and URSI coefficient files
  !  ccirnm = 'ccir.cofcnts'
  !  ursinm = 'ursi.cofcnts'

  do iBlock = 1, nBlocks
     iTemperature(:,:,:,iBlock) = 200.0
     eTemperature(:,:,:,iBlock) = 200.0
     do iAlt = -1, nAlts+2
        do iLon=-1,nLons+2
           do iLat=-1,nLats+2

              geo_lat = Latitude(iLat,iBlock)*180.0/pi
              geo_lon = Longitude(iLon,iBlock)*180.0/pi

              geo_alt(1) = Altitude_GB(iLon, iLat, iAlt, iBlock)/1000.0
              geo_lst    = mod(utime/3600.0+geo_lon/15.0,24.0)

              if (geo_lat < -90.0) then
                 geo_lat = -180.0 - geo_lat
                 geo_lon = mod(geo_lon + 180.0,360.0)
              elseif (geo_lat > 90.0) then
                 geo_lat = 180.0 - geo_lat
                 geo_lon = mod(geo_lon + 180.0,360.0)
              endif

              !  write(*,*) "iri : ", geo_lat, geo_lon, -f107, -iJulianDay, &
              !             utime/3600.+25.,geo_alt,nzkm

              ! zero out densities that are not set below
              IRIDensity(iLon,iLat,iAlt,:,iBlock) = 1.0

              call iri90 (jf,jmag,geo_lat,geo_lon,-f107,-iJulianDay, &
                   utime/3600.+25.,geo_alt,nzkm,'UA/DataIn/ccir.cofcnts', &
                   'UA/DataIn/ursi.cofcnts',outf,oarr,iProc)

              IRIDensity(iLon,iLat,iAlt,ie_,iBlock)    = abs(outf(1,1))
              IRIDensity(iLon,iLat,iAlt,iO_4SP_,iBlock) = &
                   abs(outf(5,1)*outf(1,1))/100.0+1
              IRIDensity(iLon,iLat,iAlt,iO2P_,iBlock)  = &
                   abs(outf(8,1)*outf(1,1))/100.0+1
              IRIDensity(iLon,iLat,iAlt,iNOP_,iBlock)  = &
                   abs(outf(9,1)*outf(1,1))/100.0+1
              IRIDensity(iLon,iLat,iAlt,iHP_,iBlock)   = &
                   1.0 ! abs(outf(6,1)*outf(1,1))/100.0+1
              IRIDensity(iLon,iLat,iAlt,iHeP_,iBlock)  = &
                   1.0 ! abs(outf(7,1)*outf(1,1))/100.0+1

              ! IRI Ne from 60(day)/80(night) - 1000 km (-1 missing), 
              !  ions 100-1000 km
              ! IRI temperatures from 120-3000 km (-1 is missing)

              if (outf(4,1) < 0.) then
                 eTemperature(iLon,iLat,iAlt,iBlock) =  200.0
                 ITemperature(iLon,iLat,iAlt,iBlock) =  200.0
              else
                 eTemperature(iLon,iLat,iAlt,iBlock) = outf(4,1)
                 ITemperature(iLon,iLat,iAlt,iBlock) = outf(3,1)
              endif

              IDensityS(iLon,iLat,iAlt,:,iBlock) = &
                   IRIDensity(iLon,iLat,iAlt,:,iBlock)

              IDensityS(iLon,iLat,iAlt,nIons,iBlock) = 0.0
              do iIon = 1, nIons-1
                 IDensityS(iLon,iLat,iAlt,nIons,iBlock) = &
                      IDensityS(iLon,iLat,iAlt,nIons,iBlock) + &
                      IDensityS(iLon,iLat,iAlt,iIon,iBlock)
              enddo

           enddo
        enddo
     enddo

  enddo

end subroutine init_iri
