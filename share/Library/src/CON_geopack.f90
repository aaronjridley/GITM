!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_geopack
  use ModNumConst, ONLY: cDegToRad, cHalfPi, cTwoPi
  implicit none

  ! Contains some subroutine of the geopack code (by N.V.Tsyganenko), 
  ! rewritten as the .f90 module procedures. 
  ! Added procedures: JulianDay(A.Ridley) and a computation for 
  ! the coordinate transformation like HGI=>other systems

  real:: GeiGse_DD(3,3), HgiGse_DD(3,3), GeiGsm_DD(3,3), GsmGse_DD(3,3)
  real:: DipoleStrengthGeopack, AxisMagGeo_D(3)

  ! For Earth we can use Geopack to get rotation axis too
  real:: RotAxisPhiGeopack   = -1.0
  real:: RotAxisThetaGeopack = -1.0

  ! Offset longitude angle for HGI in degrees and in radians
  ! Note: dLongitudeHgiDeg is the input value, which can be negative
  !       with a special meaning. For calculations use dLongitudeHgi only!!
  real:: dLongitudeHgiDeg = 0.0, dLongitudeHgi = 0.0

  !  SunE(arth)M(oon)B(arycenter) - The distance from the Sun to
  !                                 the Earth-and-Moon barycentre
  real:: SunEMBDistance

  public:: geopack_recalc ! recalculate all quantities as needed
  public:: test_geopack   ! unit test
  
  interface geopack_recalc
     module procedure geopack_recalc_array, geopack_recalc
  end interface

contains
  !=======================================================================
  integer function JulianDay(iYear,iMonth,iDay)
    ! Coded by A.Ridley
    integer,intent(in)::iYear,iMonth,iDay
    integer, dimension(1:12),parameter :: nDayInMonth_I = (/ &
         31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    !-------------------------------------------------------------------
    JulianDay=iDay
    if(iMonth>1)JulianDay=JulianDay+sum(nDayInMonth_I(1:iMonth-1))
    if(iMonth>2.and.mod(iYear,4)==0)JulianDay=JulianDay+1

  end function JulianDay
  !=======================================================================
  subroutine geopack_sun(iYear, jDay, iHour, iMin, iSec,&
       GSTime, SunLongitude, Obliq)
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

    double precision:: DJ,FDAY
    real::Century,VL,G !Miscellaneous
    real, parameter:: cDegToRadHere=1.0/57.295779513

    character(len=*), parameter:: NameSub = 'geopack_sun'
    !----------------------------------------------------------------------
    if(iYear < 1901 .or. iYear > 2099)then
       write(*,*) NameSub,' ERROR: No ephemeris data for the year of ',iYear
       call CON_stop('CON_geopack ERROR')
    end if
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

    ! The following formula is for the distance from the Sun to the
    ! Earth+Moon barycentre. See Eq.(36) from
    ! Heliospheric Coordinate Systems by M.Franz and D.Harper
    ! Planetary and Space Scieance, V.50, 217ff(2002),
    ! Also see the corrected version in 
    ! http://www.space-plasma.qmul.ac.uk/heliocoords/
    ! Added by I.Sokolov&I.Roussev, 08.17.03
    SunEMBDistance=1.000140-0.016710*cos(G)-0.000140*cos(2*G)

    IF(SunLongitude > 6.2831853) SunLongitude=SunLongitude - 6.2831853
    IF(SunLongitude < 0.) SunLongitude = SunLongitude + 6.2831853
    Obliq = (23.45229-0.0130125*Century)*cDegToRadHere
    GeiGse_DD=&
         matmul(rot_matrix_x(Obliq),rot_matrix_z(SunLongitude-9.924E-5))

  end subroutine geopack_sun
  !===========================================================================
  subroutine geopack_mag_axis(iYearIn, iDayIn)

    ! This was part of the RECALC subroutine from geopack.f by Tsyganenko
    !
    ! 1/26/2010: Darren De Zeeuw extend to 2015 with updated 
    !            IGRF-11 coefficients.
    !
    ! 01/22/2016: G.Toth updated to IGRF-12 coefficients to 2020.
    !             Fixed JDAY/365 to (JDAY-1)/365.25 (as in GEOPACK-2008).
    !             Rewrote the whole thing.

    integer, intent(in):: iYearIn  ! Year number (four digits) 
    integer, intent(in):: iDayIn   ! Day of year (day 1 = JAN 1)

    integer:: iYearLast=0, iDayLast=0  ! Store year and day of last call
    logical:: DoWarnSmallYear = .true. ! Only warn once iYearIn < MinYear
    logical:: DoWarnLargeYear = .true. ! Only warn once iYearIn > MaxYear

    integer:: iYear ! limited year value
    integer:: iDay  ! limited day value

    ! Year range
    integer, parameter:: MinYear = 1965
    integer, parameter:: MaxYear = 2015

    ! IGRF coefficients are given every 5 year
    integer, parameter:: DnYear = 5
    integer, parameter:: nEpoch  = (MaxYear - MinYear)/DnYear + 1

    ! Array of IGRF coefficients from MinYear to MaxYear
    ! The 1st column is -1* 3rd element of the IGRF G coefficients (G11)
    ! The 2nd column is -1* 2nd element of the IGRF H coefficients (H11)
    ! The 3rd column is +1* 2nd element of the IGRF G coefficients (G10)
    ! These provide the X, Y and Z components of the magnetic dipole in nT.

    real, parameter:: Dipole_DI(3,nEpoch) = reshape( (/ &
         -2119.,   +5776.,   -30334.,   & ! 1965
         -2068.,   +5737.,   -30220.,   & ! 1970
         -2013.,   +5675.,   -30100.,   & ! 1975
         -1956.,   +5604.,   -29992.,   & ! 1980
         -1905.,   +5500.,   -29873.,   & ! 1985
         -1848.,   +5406.,   -29775.,   & ! 1990
         -1784.,   +5306.,   -29692.,   & ! 1995
         -1728.2,  +5186.1,  -29619.4,  & ! 2000
         -1669.05, +5077.99, -29554.63, & ! 2005
         -1586.42, +4944.26, -29496.57, & ! 2010
         -1501.0,  +4797.1,  -29442.0   & ! 2015
         /), (/3, nEpoch/) )

    integer:: iEpoch        ! Index of the 5 year "epoch"
    real:: Weight1, Weight2 ! interpolation weights
    real:: Dipole_D(3)      ! interpolated dipole strength

    character(len=*), parameter:: NameSub = 'geopack_mag_axis'
    !-------------------------------------------------------------------------
    iYear = iYearIn
    iDay  = iDayIn

    if(iYearIn < MinYear) then
       if(DoWarnSmallYear)then
          write(*,*) NameSub, ' WARNING: no IGRF coefficients before year ', &
               MinYear
          write(*,*)NameSub,': setting iYear=', MinYear,' iDay=1'
          DoWarnSmallYear = .false.
       end if
       iYear = MinYear
       iDay  = 1
    endif
    if(iYearIn > MaxYear+DnYear) then
       if(DoWarnLargeYear)then
          write(*,*) 'WARNING!!! Update IGRF coefficients in ',NameSub, &
               ' in share/Library/src/CON_geopack !!!'
          write(*,*) NameSub, ': no IGRF coefficients beyond year ', MaxYear
          write(*,*)NameSub,': setting iYear=', MaxYear+DnYear,' iDay=365'
          DoWarnLargeYear = .false.
       end if
       iYear = MaxYear+DnYear
       iDay  = 365
    endif

    ! No need to recalculate if the call uses the same values as last time
    if (iYear == iYearLast .and. iDay == iDayLast) RETURN
    iYearLast = iYear
    iDayLast  = iDay

    ! Find the epoch index for iYear
    iEpoch  = min(nEpoch-1, (iYear - MinYear)/DnYear + 1)

    ! Calculate time relative to epoch start normalized to epoch length
    Weight2 = (iYear + (iDay-1)/365.25 - MinYear - (iEpoch-1)*DnYear)/DnYear
    Weight1 = 1.0 - Weight2

    Dipole_D = Weight1*Dipole_DI(:,iEpoch) + Weight2*Dipole_DI(:,iEpoch+1)

    ! Take the dipole strength with negative magnitude (conventional)
    DipoleStrengthGeopack = -sqrt(sum(Dipole_D**2))

    ! The negative sign is because the dipole strength is negative
    AxisMagGeo_D = Dipole_D/DipoleStrengthGeopack

    ! Convert from nT to Tesla (SI units)
    DipoleStrengthGeopack = DipoleStrengthGeopack*1e-9

  end subroutine geopack_mag_axis
  !===========================================================================
  subroutine geopack_recalc_array(iTimeIn_I)

    ! Allow calling geopack_recalc with an array argument.
    ! Elements of iTimeIn_I are year, month, day, hour, min, sec, ...
    ! Elements beyond second are ignored. Missing elements are set to 0.

    integer, intent(in):: iTimeIn_I(:)
    integer:: iTime_I(6), n
    !------------------------------------------------------------------------
    iTime_I=0
    n = min(6,size(iTimeIn_I))
    iTime_I(1:n) = iTimeIn_I(1:6)
    call geopack_recalc(iTime_I(1), iTime_I(2), iTime_I(3), &
         iTime_I(4), iTime_I(5), iTime_I(6))

  end subroutine geopack_recalc_array
  !===========================================================================
  subroutine geopack_recalc(iYear,iMonth,iDay,iHour,iMin,iSec)

    use ModCoordTransform, ONLY: rot_matrix_z, rot_matrix_x

    ! Updates matrices for the coordinate transformations
    ! Computations for GeiGse_DD and GeiGsm_DD are from the subroutine
    ! RECALC of geopack.f by N.V.Tsyganenko
    ! Computations for GeiHgi_DD - i.Roussev and I.Sokolov,
    ! igorsok@umich.edu, phone (734)647-4705

    ! 3/9/2005: G.Toth - corrected HgiGse_DD calculation,
    !                    which was 180 degrees off. 
    !                    NOTE: the GeiHgi_DD is only defined in the test.

    integer,intent(in):: iYear, iMonth, iDay, iHour, iMin, iSec

    integer::jDay
    real::AxisMagGei_D(3),GSTime,SunLongitude,Obliq
    real,parameter :: cLongAscNodeSolEquator = 75.77*cDegToRad
    ! Inclination of the solar equator on the ecliptic of date
    real,parameter :: cInclinationSolEquator = 7.25*cDegToRad
    integer,parameter::x_=1,y_=2,z_=3
    !-------------------------------------------------------------------
    jDay=JulianDay(iYear,iMonth,iDay)
    call geopack_mag_axis(iYear,jDay)
    call geopack_sun(iYear,jDay,iHour,iMin,iSec,GSTime,SunLongitude,Obliq)

    RotAxisPhiGeoPack   = modulo(cHalfPi - (SunLongitude-9.924E-5), cTwoPi)
    RotAxisThetaGeopack = Obliq

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
         sqrt(sum(GeiGsm_DD(:,y_)**2))

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

  end subroutine geopack_recalc
  !===========================================================================
  subroutine test_geopack

    ! Coded by I Sokolov, and I Roussev, 08.16.2003
    ! The test compares the mean position of the pole of the solar
    ! eqautor (see SUN,2001, p.C3) with the unity vector n_z of
    ! the Heliographic inertial coordinates with respect to GEI  
    ! coordinate system

    real:: GeiHgi_DD(3,3)
    real, parameter::RightAscension=286.13*cDegToRad,&
         Declination=63.87*cDegToRad
    integer:: iYear=2000, iMonth, iDay, iHour, iMin=0, iSec=0
    !----------------------------------------------------------------------
    ! For perihelion
    iMonth=1;iDay=3;iHour=5
    call geopack_recalc(iYear,iMonth,iDay,iHour,iMin,iSec)
    write(*,'(a,f14.4)')  'DipoleStrength=', DipoleStrengthGeopack*1e9
    write(*,'(a,3f14.10)')'AxisMagGeo_D=', AxisMagGeo_D
    write(*,'(a,/,3f14.10,/,3f14.10,/,3f14.10)')'GsmGse_DD=', GsmGse_DD
    
    write(*,'(a,f14.10,a)')'SunEMBDistance=',SunEMBDistance,&
         ', should be 0.98329'
    GeiHgi_DD=matmul(GeiGse_DD,transpose(HgiGse_DD))
    write(*,'(a,3es16.8)')&
         'Solar rotation axis vector calculated as GeiHgi_DD(:,3)',&
         GeiHgi_DD(:,3)
    write(*,'(a,3es16.8)')&
         'The vector calculated in terms of RightAsc=286.13,Declin=63.87',&
         cos(RightAscension)*cos(Declination),&
         sin(RightAscension)*cos(Declination),&
         sin(Declination)

    ! For aphelion
    iMonth=7;iDay=4;iHour=0
    call geopack_recalc(iYear,iMonth,iDay,iHour,iMin,iSec)
    write(*,'(a,f14.10,a)')'SunEMBDistance=',SunEMBDistance,&
         ', should be 1.01671'
    GeiHgi_DD=matmul(GeiGse_DD,transpose(HgiGse_DD))
    write(*,'(a,3es16.8)')&
         'Solar rotation axis vector calculated as GeiHgi_DD(:,3)',&
         GeiHgi_DD(:,3)
    write(*,'(a,3es16.8)')&
         'The vector calculated in terms of RightAsc=286.13,Declin=63.87',&
         cos(RightAscension)*cos(Declination),&
         sin(RightAscension)*cos(Declination),&
         sin(Declination)

  end subroutine test_geopack

end module CON_geopack
