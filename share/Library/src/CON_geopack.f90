!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module CON_geopack
  use ModNumConst
  implicit none

  !Contains some subroutine of the geopack code (by N.V.Tsyganenko), 
  !rewritten as the .f90 module procedures. 
  !Added procedures: JulianDay(A.Ridley)and a  computation for 
  !the coordinate transformation like HGI=>other systems

  real,dimension(3,3)::GeiGse_DD,HgiGse_DD,GeiGsm_DD,GsmGse_DD
  real,dimension(3)::AxisMagGeo_D

  ! Offset longitude angle for HGI in degrees and in radians
  ! Note: dLongitudeHgiDeg is the input value, which can be negative
  !       with a special meaning. For calculations use dLongitudeHgi only!!
  real :: dLongitudeHgiDeg = 0.0, dLongitudeHgi = 0.0

  !  SunE(arth)M(oon)B(arycenter) - The distance from the Sun to
  !                                 the Earth-and-Moon barycentre
  real::SunEMBDistance
contains
  !----------------------------------------------------------------------
  integer function JulianDay(iYear,iMonth,iDay)
    !Coded by A.Ridley
    integer,intent(in)::iYear,iMonth,iDay
    integer, dimension(1:12),parameter :: nDayInMonth_I = (/ &
         31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    JulianDay=iDay
    if(iMonth>1)JulianDay=JulianDay+sum(nDayInMonth_I(1:iMonth-1))
    if(iMonth>2.and.mod(iYear,4)==0)JulianDay=JulianDay+1
  end function JulianDay
  !----------------------------------------------------------------------
  subroutine CON_sun(iYear,jDay,iHour,iMin,iSec,&
       GSTime,SunLongitude,Obliq)
    !
    !  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE 
    !  TRANSFORMATIONS
    !  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME 
    !  AND SEASON)
    !  From geopack.f by N.V.Tsyganenko
    use ModCoordTransform
    integer,intent(in)::iYear,jDay,iHour,iMin,iSec
    !-------  INPUT PARAMETERS:
    !  iYear,jDay,iHour,iMin,iSec -  YEAR, DAY, AND UNIVERSAL TIME 
    !  IN HOURS, MINUTES,
    !    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
    !
    real,intent(out):: GSTime,SunLongitude,Obliq
    !-------  OUTPUT PARAMETERS:
    !  GSTime - GREENWICH MEAN SIDEREAL TIME
    !  SunLongitude - The Sun Longitude
    !  Obliq        - The inclination of the Earth equatorial plane
    !  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
    !  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
    !     LAST MODIFICATION:  JAN 5, 2001 (NO ESSENTIAL CHANGES, BUT
    !     SOME REDUNDANT STATEMENTS TAKEN OUT FROM THE PREVIOUS VERSION)
    !
    !     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead

    real*8:: DJ,FDAY
    real::Century,VL,G !Miscellaneous
    real,parameter::cDegToRadHere=1.0/57.295779513
    !----------------------------------------------------------------------
    IF(iYear.LT.1901.OR.iYear.GT.2099)then
       write(*,*)'CON_geopack ERROR: No ephemers data for the year of ',iYear
       call CON_stop('CON_geopack ERROR')
    end IF
    FDAY=dble(IHOUR*3600+iMIN*60+ISEC)/86400.E0
    DJ=365*(IYear-1900)+(IYear-1901)/4+jDAY-0.5E0+FDAY
    Century=DJ/36525.0
    VL=MOD(279.696678+0.9856473354*DJ,360.D0)
    GSTime=&
         MOD(279.690983+0.9856473354*DJ+360.0*FDAY+180.,360.D0)&
         *cDegToRadHere
    G=MOD(358.475845+0.985600267*DJ,360.D0)*cDegToRadHere
    SunLongitude=(VL+(1.91946-0.004789*Century)*SIN(G)+&
         0.020094*SIN(2.*G))*cDegToRadHere

    !The following formula is for the distance from the Sun to the
    !Earth+Moon barycentre. See Eq.(36) from
    !Heliospheric Coordinate Systems by M.Franz and D.Harper
    !Planetary and Space Scieance, V.50, 217ff(2002),
    !Also see the corrected version in 
    !http://www.space-plasma.qmul.ac.uk/heliocoords/
    !Added by I.Sokolov&I.Roussev, 08.17.03
    SunEMBDistance=1.000140-0.016710*cos(G)-0.000140*cos(G+G)

    IF(SunLongitude.GT.6.2831853)SunLongitude=SunLongitude-6.2831853
    IF(SunLongitude.LT.0.)SunLongitude=SunLongitude+6.2831853
    Obliq=(23.45229-0.0130125*Century)*cDegToRadHere
    GeiGse_DD=&
         matmul(rot_matrix_x(Obliq),rot_matrix_z(SunLongitude-9.924E-5))
  end subroutine CON_sun
  !---------------------------------------------------------------------------
  subroutine CON_mag_axis(iYear,jDay)
    !This is a part of the RECALC subroutine form geopack.f by Tsyganenko
    ! NOTE: Code modified 1/26/2010 by DDZ to extend to 2015
    !       with updated IGRF-11 coefficients
    integer,intent(in)::iYear,jDay
    !-----INPUT PARAMETERS:
    !
    !     IYear   -  YEAR NUMBER (FOUR DIGITS)
    !     jDAY  -  DAY OF YEAR (DAY 1 = JAN 1)
    integer::IYE=0,IDE=0,IPR=0
    integer::iY
    real::F1,F2,H11,G10,G11,SQR,DT
    IF (IYear.EQ.IYE.AND.jDAY.EQ.IDE) return
    !
    !   IYE AND IDE ARE THE CURRENT VALUES OF YEAR AND DAY NUMBER
    !
    IY=IYear
    IDE=jDAY
    IF(IY.LT.1965) THEN
       write(*,*)'No IGRF coefficients, year set to 1965 from ',IY
       IY=1965
    ENDIF
    IF(IY.GT.2015) THEN
       write(*,*)'No IGRF coefficients, year set to 2015 from ',IY
       write(*,*)' Is it time to upgrade IGRF coefficients in CON_mag_axis?'
       IY=2015
    ENDIF
    !
    !  WE ARE RESTRICTED BY THE INTERVAL 1965-2005,
    !  FOR WHICH THE IGRF COEFFICIENTS ARE KNOWN; IF IYR IS OUTSIDE THIS INTERVAL
    !  THE SUBROUTINE PRINTS A WARNING (BUT DOES NOT REPEAT IT AT NEXT INVOCATIONS)
    !
    !!! IF(IY.NE.IYear.AND.IPR.EQ.0)&
    !!!     WRITE (*,*) 'No Igrf Coefficients are availble for the year ',&
    !!!     IYear,'We use date for year ',IY
    IF(IY.NE.IYear) IPR=1
    IYE=IY   
    !
    !  LINEAR INTERPOLATION OF THE GEODIPOLE MOMENT COMPONENTS BETWEEN THE
    !  VALUES FOR THE NEAREST EPOCHS:
    !
    IF (IY.LT.1970) THEN                            !1965-1970
       F2=(FLOAT(IY)+FLOAT(jDAY)/365.-1965.)/5.
       F1=1.E0-F2
       G10=30334.*F1+30220.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-2119.*F1-2068.*F2
       H11=5776.*F1+5737.*F2
    ELSEIF (IY.LT.1975) THEN                        !1970-1975
       F2=(FLOAT(IY)+FLOAT(jDAY)/365.-1970.)/5.
       F1=1.E0-F2
       G10=30220.*F1+30100.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-2068.*F1-2013.*F2
       H11=5737.*F1+5675.*F2
    ELSEIF (IY.LT.1980) THEN                        !1975-1980
       F2=(dble(IY)+dble(jDAY)/365.-1975.)/5.
       F1=1.E0-F2
       G10=30100.*F1+29992.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-2013.*F1-1956.*F2
       H11=5675.*F1+5604.*F2
    ELSEIF (IY.LT.1985) THEN                        !1980-1985
       F2=(FLOAT(IY)+FLOAT(jDAY)/365.-1980.)/5.
       F1=1.E0-F2
       G10=29992.*F1+29873.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-1956.*F1-1905.*F2
       H11=5604.*F1+5500.*F2
    ELSEIF (IY.LT.1990) THEN                        !1985-1990
       F2=(FLOAT(IY)+FLOAT(jDAY)/365.-1985.)/5.
       F1=1.E0-F2
       G10=29873.*F1+29775.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-1905.*F1-1848.*F2
       H11=5500.*F1+5406.*F2
    ELSEIF (IY.LT.1995) THEN                        !1990-1995
       F2=(FLOAT(IY)+FLOAT(jDAY)/365.-1990.)/5.
       F1=1.E0-F2
       G10=29775.*F1+29692.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-1848.*F1-1784.*F2
       H11=5406.*F1+5306.*F2
    ELSEIF (IY.LT.2000) THEN                        !1995-2000
       F2=(FLOAT(IY)+FLOAT(jDAY)/365.-1995.)/5.
       F1=1.E0-F2
       G10=29692.*F1+29619.4*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-1784.*F1-1728.2*F2
       H11=5306.*F1+5186.1*F2
    ELSEIF (IY.LT.2005) THEN                        !2000-2005
       F2=(FLOAT(IY)+FLOAT(jDAY)/365.-2000.)/5.
       F1=1.E0-F2
       G10=29619.4*F1+29554.63*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-1728.2*F1-1669.05*F2
       H11=5186.1*F1+5077.99*F2
    ELSEIF (IY.LT.2010) THEN                        !2005-2010
       F2=(FLOAT(IY)+FLOAT(jDAY)/365.-2005.)/5.
       F1=1.E0-F2
       G10=29554.63*F1+29496.5*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-1669.05*F1-1585.9*F2
       H11=5077.99*F1+4945.1*F2
    ELSE                                            !2010-2015
       !   LINEAR EXTRAPOLATION BEYOND 2010:
       !     Use coefficient of F2 above and DT*(Coeff_F2-Coeff_F1)/5
       !
       DT  = IY + jDAY/365.0 - 2010
       G10 = 29496.5 - 11.6*DT
       G11 = -1585.9 + 16.6*DT
       H11 =  4945.1 - 26.6*DT
    ENDIF
    !
    !  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EzMAG IN GEO COORD.SYSTEM:
    !   SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
    !         ST0 * CL0                ST0 * SL0                CT0
    !
    SQR=G11**2+H11**2+G10**2
    AxisMagGeo_D(3)=G10/SQR
    AxisMagGeo_D(1)=-G11/SQR
    AxisMagGeo_D(2)=-H11/SQR
  end subroutine CON_mag_axis
  !------------------------------------------------------------------------
  subroutine CON_recalc(iYear,iMonth,iDay,iHour,iMin,iSec)
    use ModCoordTransform,ONLY:rot_matrix_z,rot_matrix_x

    !Updates matrices for the coordinate transformations
    !Computations for GeiGse_DD and GeiGsm_DD are from the subroutine
    !RECALC of geopack.f by N.V.Tsyganenko
    !Computations for GeiHgi_DD - i.Roussev and I.Sokolov,
    !igorsok@umich.edu, phone (734)647-4705

    ! 3/9/2005: G.Toth - corrected HgiGse_DD calculation,
    !                    which was 180 degrees off. 
    !                    NOTE: the GeiHgi_DD is only defined in the test.

    integer,intent(in)::iYear,iMonth,iDay,iHour,iMin,iSec
    integer::jDay
    real::AxisMagGei_D(3),GSTime,SunLongitude,Obliq
    real,parameter :: cLongAscNodeSolEquator = 75.77*cDegToRad
    ! Inclination of the solar equator on the ecliptic of date
    real,parameter :: cInclinationSolEquator = 7.25*cDegToRad
    integer,parameter::x_=1,y_=2,z_=3
    !-------------------------------------------------------------------
    jDay=JulianDay(iYear,iMonth,iDay)
    call CON_mag_axis(iYear,jDay)
    call CON_sun(iYear,jDay,iHour,iMin,iSec,GSTime,SunLongitude,Obliq)

    GeiGse_DD=&
         matmul(rot_matrix_x(Obliq),rot_matrix_z(SunLongitude-9.924E-5))
    !
    !   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION  
    !   DUE TO THE ORBITAL MOTION OF THE EARTH   

    HgiGse_DD = matmul( &
         rot_matrix_x(-cInclinationSolEquator),&
         rot_matrix_z( SunLongitude - cLongAscNodeSolEquator ))

    ! Offset the HGI coordinate system if required
    if(dLongitudeHgi > 0.0) &
         HgiGse_DD = matmul( rot_matrix_z(-dLongitudeHgi), HgiGse_DD)

    !   THE COMPONENTS OF THE UNIT VECTOR EXGSM=EXGSE IN THE
    !   SYSTEM GEI POINTING FROM THE EARTH'S CENTER TO THE SUN:
    GeiGsm_DD(:,x_)=GeiGse_DD(:,x_)


    !   THE COMPONENTS OF THE UNIT VECTOR EZSM=EZMAG
    !   IN THE SYSTEM GEI:
    AxisMagGei_D=matmul(rot_matrix_z(GSTime),AxisMagGeo_D)

    !
    !  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EYGSM 
    !  IN THE SYSTEM GEI BY TAKING THE VECTOR PRODUCT
    !   D x S AND NORMALIZING IT TO UNIT LENGTH:

    GeiGsm_DD(1,y_)=AxisMagGei_D(2)*GeiGsm_DD(3,x_)-&
         AxisMagGei_D(3)*GeiGsm_DD(2,x_)
    GeiGsm_DD(2,y_)=AxisMagGei_D(3)*GeiGsm_DD(1,x_)-&
         AxisMagGei_D(1)*GeiGsm_DD(3,x_)
    GeiGsm_DD(3,y_)=AxisMagGei_D(1)*GeiGsm_DD(2,x_)-&
         AxisMagGei_D(2)*GeiGsm_DD(1,x_)
    GeiGsm_DD(:,y_)=GeiGsm_DD(:,y_)/&
         sqrt(dot_product(GeiGsm_DD(:,y_),GeiGsm_DD(:,y_)))

    !
    !   THEN IN THE GEI SYSTEM THE UNIT VECTOR 
    !   Z = EZGSM = EXGSM x EYGSM = S x Y
    !   HAS THE COMPONENTS:
    GeiGsm_DD(1,z_)=GeiGsm_DD(2,x_)*GeiGsm_DD(3,y_)-&
         GeiGsm_DD(3,x_)*GeiGsm_DD(2,y_)
    GeiGsm_DD(2,z_)=GeiGsm_DD(3,x_)*GeiGsm_DD(1,y_)-&
         GeiGsm_DD(1,x_)*GeiGsm_DD(3,y_)
    GeiGsm_DD(3,z_)=GeiGsm_DD(1,x_)*GeiGsm_DD(2,y_)-&
         GeiGsm_DD(2,x_)*GeiGsm_DD(1,y_)

    GsmGse_DD=matmul(transpose(GeiGsm_DD),GeiGse_DD)

  end subroutine CON_recalc
  !----------------------------------------------------------------
  subroutine CON_test_geopack
    ! Coded by I Sokolov, and I Roussev, 08.16.2003
    ! The test compares the mean position of the pole of the solar
    ! eqautor (see SUN,2001, p.C3) with the unity vector n_z of
    ! the Heliographic inertial coordinates with respect to GEI  
    ! coordinate system
    real::GeiHgi_DD(3,3)
    real,parameter::RightAssention=286.13*cDegToRad,&
         Declination=63.87*cDegToRad
    integer::iYear=2000,iMonth,iDay,iHour,iMin=0,iSec=0
    !For perihelion
    iMonth=1;iDay=3;iHour=5
    call CON_recalc(iYear,iMonth,iDay,iHour,iMin,iSec)
    write(*,'(a,f14.10,a)')'SunEMBDistance=',SunEMBDistance,&
         ', should be 0.98329'
    GeiHgi_DD=matmul(GeiGse_DD,transpose(HgiGse_DD))
    write(*,'(a,3es16.8)')&
         'Solar rotation axis vector calculated as GeiHgi_DD(:,3)',&
         GeiHgi_DD(:,3)
    write(*,'(a,3es16.8)')&
         'The vector calculated in terms of RightAss=286.13,Declin=63.87',&
         cos(RightAssention)*cos(Declination),&
         sin(RightAssention)*cos(Declination),&
         sin(Declination)
    !For aphelion
    iMonth=7;iDay=4;iHour=0
    call CON_recalc(iYear,iMonth,iDay,iHour,iMin,iSec)
    write(*,'(a,f14.10,a)')'SunEMBDistance=',SunEMBDistance,&
         ', should be 1.01671'
    GeiHgi_DD=matmul(GeiGse_DD,transpose(HgiGse_DD))
    write(*,'(a,3es16.8)')&
         'Solar rotation axis vector calculated as GeiHgi_DD(:,3)',&
         GeiHgi_DD(:,3)
    write(*,'(a,3es16.8)')&
         'The vector calculated in terms of RightAss=286.13,Declin=63.87',&
         cos(RightAssention)*cos(Declination),&
         sin(RightAssention)*cos(Declination),&
         sin(Declination)

  end subroutine CON_test_geopack

end Module CON_geopack
