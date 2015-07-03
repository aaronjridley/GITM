C Subroutine SOLZEN
C
C This software is part of the GLOW model.  Use is governed by the Open Source
C Academic Research License Agreement contained in the file glowlicense.txt.
C For more information see the file glow.txt.
C
C Stan Solomon, 1988
C Temporary Y2K fix-up, SCS, 2005.
C
C Returns Solar Zenith Angle SZA in degrees for specified date in form
C yyddd, universal time in seconds, geographic latitude and longitude
C in degrees.
C
C
      SUBROUTINE SOLZEN (IDATE, UT, GLAT, GLONG, SZA)
C
      DATA PI/3.1415926536/
C
      RLAT = GLAT * PI/180.
      RLONG = GLONG * PI/180.
      CALL SUNCOR (IDATE, UT, SDEC, SRASN, GST)
      RH = SRASN - (GST+RLONG)
      COSSZA = SIN(SDEC)*SIN(RLAT) + COS(SDEC)*COS(RLAT)*COS(RH)
      SZA = ACOS(COSSZA) * 180./PI
      RETURN
      END
C
C
C
C
C Subroutine SUNCOR returns the declination SDEC and right ascension
C SRASN of the sun in GEI coordinates, radians, for a given date IDATE
C in yyddd format and universal time UT in seconds.  Greenwich Sidereal
C Time GST in radians is also returned.  Reference:  C.T. Russell,
C Geophysical Coordinate Transforms.
C
      SUBROUTINE SUNCOR (IDATE, UT, SDEC, SRASN, GST)
      DATA PI/3.1415926536/
C
      FDAY=UT/86400.
      IYR=IDATE/1000
      IDAY=IDATE-IYR*1000
C
C Temporary Y2K fix-up:
C Should work with either yyddd or yyyyddd format from 1950 to 2050.
C Note deteriorating accuracy after ~2050 anyway.
C Won't work after 2100 due to lack of a leap year.
      IF (IYR .GE. 1900) IYR=IYR-1900
      IF (IYR .LT. 50) IYR=IYR+100
C
      DJ=365*IYR+(IYR-1)/4+IDAY+FDAY-0.5
      T=DJ/36525.
      VL=DMOD(279.696678+.9856473354*DJ,360.)
      GST=DMOD(279.696678+.9856473354*DJ+360.*FDAY+180.,360.) * PI/180.
      G=DMOD(358.475845+.985600267*DJ,360.) * PI/180.
      SLONG=VL+(1.91946-.004789*T)*SIN(G)+.020094*SIN(2.*G)
      OBLIQ=(23.45229-0.0130125*T) *PI/180.
      SLP=(SLONG-.005686) * PI/180.
      SIND=SIN(OBLIQ)*SIN(SLP)
      COSD=SQRT(1.-SIND**2)
      SDEC=ATAN(SIND/COSD)
      SRASN=3.14159-ATAN2(1./TAN(OBLIQ)*SIND/COSD,-COS(SLP)/COSD)
      RETURN
      END
