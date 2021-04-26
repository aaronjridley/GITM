!BP (8/24/2020): Compute the time scales for eddy diffusion, horizontal winds 
!and vertical motion
subroutine timescales

  use ModConstants, only:  Boltzmanns_Constant
  use ModPlanet
  use ModSizeGITM, only: nAlts  
  use ModGITM
  use ModSources, only: KappaEddyDiffusion
  use ModTime 

  integer :: iLon, iLat, iAlt
  integer ::  iBlock = -1
  real, dimension(nAlts) :: tauEddy
  real, dimension(nAlts) :: tauHorizontal, tauVertical
  real :: surfaceGravity = 8.87 !m/s^2
  real :: scaleHeight, SSLon, SSLat
  real :: RP_hours = RP_Venus/3600.0
  logical :: exist   

  call get_subsolar(CurrentTime, VernalTime, SSLon, SSLat) 
  if (SSLon == 0) then
    SSLon = 1e-2
  endif
  
  if (SSLat == 0) then
    SSLat = 1e-2
  endif
  
  call LocationIndex(SSLon, SSLat, iBlock, iLon, iLat, rLon, rLat)
!  if (iBlock .eq. 1) then
!    write(*,*) "Trying to find:", SSLon*180.0/pi, SSLat*180.0/pi
!    write(*,*) "Found:", Longitude(iLon,iBlock)*180.0/pi, &
!                         Latitude(iLat, iBlock)*180.0/pi
!  endif

  if (iBlock .eq. 1) then
    !Check for 12 LST
    geo_lon = SSLon*180.0/pi
    geo_lst = mod(UTime/3600.0 + SSLon*180.0/(360*pi/RP_hours), &
                         RP_hours)

    !Convert to a 0 - 24 scale instead of the 0 - VenusHoursPerDay                     
    geo_lst = 24*geo_lst/RP_hours  

    do iAlt = 1, nAlts
      gravity = surfaceGravity * (R_Venus /(R_Venus + Altitude_GB(iLon,iLat,iAlt,iBlock)))**2
 
      SSTemperature = &
        rLon*rLat*Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt) + &
        (1-rLon)*rLat*Temperature(iLon+1,iLat,iAlt,iBlock)*TempUnit(iLon+1,iLat,iAlt) + &
        rLon*(1-rLat)*Temperature(iLon,iLat+1,iAlt,iBlock)*TempUnit(iLon,iLat+1,iAlt) + &
        (1-rLon)*(1-rLat)*Temperature(iLon+1,iLat+1,iAlt,iBlock)*TempUnit(iLon+1,iLat+1,iAlt)

      scaleHeight = Boltzmanns_Constant*SSTemperature &
                    /(Mass(iCO2_)*gravity)

      tauEddy(iAlt) = scaleHeight**2/KappaEddyDiffusion(iLon,iLat,iAlt,iBlock)       

      !Horizontal time calculation
      
      tauHorizontal(iAlt) = pi*(R_Venus + Altitude_GB(iLon,iLat,iAlt,iBlock))/&
                            Velocity(iLon,iLat,iAlt,iEast_,iBlock)
      tauVertical(iAlt) = scaleHeight/Velocity(iLon,iLat,iAlt,iNorth_,iBlock)

    enddo

    inquire(file = "eddy.txt", exist = exist)
    if (exist .eqv. .False.) then
      open(32, file = "eddy.txt", status = "new", action = "write")
      write(32, *) tauEddy
      close(32)
    end if

    inquire(file = "horizontal.txt", exist = exist)
    if (exist .eqv. .False.) then
      open(33, file = "horizontal.txt", status = "new", action = "write")
      write(33, *) tauHorizontal
      close(33)
    end if
    
    inquire(file = "vertical.txt", exist = exist)
    if (exist .eqv. .False.) then
      open(33, file = "vertical.txt", status = "new", action = "write")
      write(33, *) tauVertical
      close(33)
    end if
  endif
end subroutine timescales

!BP (5/29/2020) Read in the IR table needed in calc_sources.f90
!to include IR heating from Crips & Roland Venus IR heating rates

subroutine readVenusIRTable

  use ModIoUnit, only : UnitTmp_
  use ModSources, only: qIR_NLTE_table, diurnalHeating, semiDiurnalHeating
  real, dimension(40) :: qAlts
  logical :: exist 
  integer :: iiAlt, iSza

  inquire(file = 'DataIn/VenusIRHeatingTable.txt', exist = exist)
  if (.not. exist) then
    call stop_gitm("Can't read VenusIRHeatingTable.txt")
  endif
  open(UNIT = UnitTmp_, FILE = 'DataIn/VenusIRHeatingTable.txt', &
       STATUS='OLD', ACTION = 'READ')

  !Skip first line of table headings
  read(UnitTmp_, *)

  do iiAlt = 1,40
     read(UnitTmp_,*) qAlts(iiAlt), qIR_NLTE_table(iiAlt, :)
  enddo
  close(Unit = UnitTmp_)

  !Convert from K/day to K/sec
  qIR_NLTE_table = qIR_NLTE_table/86400

   !Reading the Crisp 1986 data                                                                    
    !In a latitude vs. altitude format with the heating rate K/day as the contour                   
    !  There are two files here containing Figure 16 and Figure 17 in the paper                     
    !  which are the diurnal and semidiurnal components respectively                                
    open(UNIT = UnitTmp_, FILE = 'DataIn/crispDiurnalData_2km5deg.txt', &
         STATUS='OLD', ACTION = 'READ')

    

    !Skip two lines of header stuff                                                                 
    read(UnitTmp_, *)
    read(UnitTmp_, *)

    !0 - 90 deg, in 5 deg increments                                                                
    do iiLat = 1,19
       read(UnitTmp_,*) diurnalHeating(iiLat, :)
    enddo
    close(Unit = UnitTmp_)

    !Read the semidiurnal file                                                                      
    open(UNIT = 97, FILE = 'DataIn/crispSemidiurnalData_2km5deg.txt', &
         STATUS='OLD', ACTION = 'READ')

    !Skip two lines of header stuff                                           
    read(97, *)
    read(97, *)

    do iiLat = 1,19
       read(97,*) semiDiurnalHeating(iiLat, :)
    enddo
    close(Unit = 97)
  
    diurnalHeating = diurnalHeating/86400.0
    semiDiurnalHeating = semiDiurnalHeating/86400.0
    
end subroutine readVenusIRTable

subroutine fill_photo

! CVS_new_code:  Dec. 19, 2011 (DP additions)
! -- Timing Fix, photoabs(58,iCO2_) = 0.0
! -- E(EUV) = 0.21 (off) or 0.18 (on)
! -- TOTAL(1) = 0.0; TOTAL(2) = 0.0
! -- populate 1-D and 3-D fields for diagnostics from RT code
!    (RadCoolingRate,LowAtmosRadRate subroutines)
! -- ALS = constant for Amanda Brech = 0.2 (for Comparison Studies Studes Only!)
! CVS_new_code:  Dec. 20, 2011 (DP additions)
! -- ALS = array from 2-D table (standard input for all cases)
! Adjustment:  May 2012
! -- KMAX = 1000., KMIN = 500.
! Adjustments:  November 2012  (standard)
! -- ALS = array from 2-D table (standard input for all cases)
! Adjustment:  June 2017
! -- KMAX = 1500., KMIN = 500.
! -- KMAX = 2000., KMIN = 500.
! Adjustment:  June-July 2017
! O-CO2 cooling coefficient enhancement
! -- k20xc = 4.e-12 * rfvto3p
! -- k20xc = 3.e-12 * rfvto3p
! New nlte_tcool code from FGG and MLV: November 2017
! Modified by S. W. Bougher : November 2017
! -- nlte_setup routine (1)
! -- nlte_tcool.F (major code)
! -- supporting subroutines (5)
! -- supporting parameters and inputs (1)
! -- supporting array declarations (real, integer) (1)

  use ModPlanet
  use ModEUV

  use ModGITM
  
  implicit none

  integer :: i, iSpecies, iIon, NWH, iWave

  ! --------------------------------------------------------------------
  ! Secondary PE enhancement factors (simplified; check against Fox & Sung)

  real,parameter ::  sfn2p = 1.40, sfop = 1.60, sfco2p = 1.50
  ! --------------------------------------------------------------------
  photoabs           = 0.0
  photoion            = 0.0
  photodis            = 0.0

  NWH = Num_WaveLengths_High

  photoabs               = 0.0
  photoabs(:,iCO2_)  = PhotoAbs_CO2
  photoabs(:,iCO_)    = PhotoAbs_CO
  photoabs(:,iO_)      = PhotoAbs_O
  photoabs(:,iN2_)    = PhotoAbs_N2
  photoabs(:,iO2_)    = PhotoAbs_O2
  photoabs(58,iCO2_) = 0.0! timing fix (DP: Nov. 2011)

 ! ---------------------------------------------------------------------
  !  Specific Photoionization Cross Sections (nIons-1)
  !  Total Ionization Cross Sections * Specific Branching Ratios
  !  Need:  O+, O2+, CO2+, N2+ Productions
  ! ---------------------------------------------------------------------

  photoion(1:NWH,iOP_)        = PhotoIon_OPlus4S(1:NWH)*sfop
  photoion(1:NWH,iO2P_)      = PhotoIon_O2(1:NWH)
  photoion(1:NWH,iCO2P_)    = PhotoIon_CO2(1:NWH)* &
       BranchingRatio_CO2_to_CO2Plus(1:NWH)*sfco2p
  photoion(1:NWH,iN2P_)      = PhotoIon_N2(1:NWH)*  &
       BranchingRatio_N2_to_N2Plus(1:NWH)*sfn2p
  photoion(1:NWH,iNOP_)     = PhotoIon_CO2(1:NWH)* &
       BranchingRatio_CO2_to_OPlus(1:NWH)*sfco2p

  PhotoIonFrom(iOP_) = iO_
  PhotoIonFrom(iO2P_) = iO2_
  PhotoIonFrom(iCO2P_) = iCO2_
  PhotoIonFrom(iN2P_) = iN2_
  PhotoIonFrom(iNOP_)   = iNO_

  ! ---------------------------------------------------------------------
  !  Specific PhotoiDissociations (Absorption-Total_Ionization) (nSpecies)
  !  Need:  CO2, O2, and N2 dissociations (only)
  ! ---------------------------------------------------------------------

  photodis(1:NWH,iCO2_)  = PhotoAbs_CO2(1:NWH)-PhotoIon_CO2(1:NWH)
  photodis(1:NWH,iN2_)   = PhotoAbs_N2(1:NWH)-PhotoIon_N2(1:NWH)
  photodis(1:NWH,iO2_)   = PhotoAbs_O2(1:NWH)-PhotoIon_O2(1:NWH)

  do iWave = 1, nWavelengths
    if (photodis(iWave,iCO2_) < 0.0) then
       photodis(iWave,iCO2_) = 0.0
    else if (photodis(iWave,iN2_) < 0.0) then
       photodis(iWave, iN2_) = 0.0
    else if (photodis(iWave,iO2_) < 0.0) then
       photodis(iWave, iO2_) = 0.0
    endif
  enddo

  PhotoElecIon = 0.0
  PhotoElecDiss = 0.0

  !Brandon Ponder (3/8/2021)
  !Write out the cross-sections
  open(unit=61, file='gitm_cross_sections.txt', action='write', status="unknown")

  do i = 1,nWavelengths
     write(61,*) (shortWavelengths(i) + longWavelengths(i))/2, &
                 photoabsorptioncrosssection(i, iCO2_), &
                 photodissociationcrosssection(i,iCO2_), &
                 photoionizationcrosssection(i,iCO2P_) !use ion index
  end do
  
  close(61)
end subroutine fill_photo

!---------------------------------------------------------------------
! Initialize Heating Rates
!---------------------------------------------------------------------

subroutine init_heating_efficiency

  use ModEUV, only: HeatingEfficiency_CB, eHeatingEfficiency_CB

  implicit none

  ! HeatingEfficiency_CB  = 0.18
  ! HeatingEfficiency_CB  = 0.19
  HeatingEfficiency_CB  = 0.20
  eHeatingEfficiency_CB = 0.0

  !  call init_radcool
  call init_nlte_setup

end subroutine init_heating_efficiency

!---------------------------------------------------------+
!
!---------------------------------------------------------+
subroutine calc_planet_sources(iBlock)

  !  All new source code in here

  use ModInputs
  use ModSources
  use ModGITM
  use ModPlanet
  use ModEUV, only:SunPlanetDistance,AveCosSza
  use ModTime
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock
  integer :: nw, iLat, iLon, iAlt
  integer :: L_LAYERS, L_LEVELS, L_NLAYRAD,L_NLEVRAD
  integer, dimension(4) :: locationOfMin
  real :: tmp2(nLons, nLats, nAlts)
  real :: tmp3(nLons, nLats, nAlts)

  ! Ls Variables
  real :: deltat
  real :: em
  real :: alpha_pms
  real :: pbs
  real :: deltat_h
  integer :: i

  call start_timing("calc_planet_sources")
  ! New sources specificly for Mars include:
  ! (1) calc_radcooling(iBlock):  added 1/31/07 (BOUGHER)
  ! (2) calc_radcode(iBlock)   :  to be added later (NELLI)
  ! (3) OCooling code : added 4/14/12 (BOUGHER)

  if (useGravityWave)  call calc_GW(iBlock)

  !\ -------------------------------------------------------------------
  ! CO2 NLTE Cooling Formulation from Miguel Lopez-Valverde (2017)
  !/
  !/ Cooling ON
  if (useRadCooling) call nlte_tcool(iBlock)

  !\
  ! RadCoolingRate is in K/s (from nlte_tcool.F routine)
  ! --------------------------------------------------
  ! RadCoolingRate = K/s
  ! RadCoolingRate/TempUnit = <T>/s
  ! TempUnit = MeanMajorMass/kb = K*kg/J = K/(m/s)^2
  ! --------------------------------------------------
  ! cp is in J/kg/K
  ! rho is in kg/m^3
  ! cp*rho is in J/K/m^3
  ! 1.0/(cp*rho) is in (K*m^3)/J
  ! --------------------------------------------------
  !/ Cooling ON
  if (useRadCooling) then
    RadCooling(1:nLons,1:nLats,1:nAlts,iBlock) = &
         RadCoolingRate(1:nLons,1:nLats,1:nAlts,iBlock)/&
         (TempUnit(1:nLons,1:nLats,1:nAlts))

    UserData2d(1:nLons,1:nLats,1,1,iBlock) = &
         RadCoolingRate(1:nLons,1:nlats,1,iBlock)
  endif

  !/ Cooling OFF (zeroed out)
  !RadCooling = 0.0
  !RadCoolingRate = 0.0

  !\ -------------------------------------------------------------------
  ! O(63 micron) Cooling Formulation from Kockarts (1970)
  !/
  !write(*,*) "UseOCooling", UseOCooling
  if (UseOCooling) then

     ! [O] cooling
     ! Initial Reference: Kockarts, G., P. Sp.Sci., Vol. 18, pp. 271-285, 1970
     ! Exact LTE (opt. thin) equations found in Banks and Kockarts [1973] pg. 22.
     ! a. Where does tau = 1.0, above which optically thin approx is valid?
     ! b. We reduce the LTE 63-um cooling rate by a factor of 2 for
     !    the non-LTE effects.[Roble et.al,1987, JGR, 92, 8745]

     tmp2 = exp(-228./(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
          TempUnit(1:nLons,1:nLats,1:nAlts)))
     tmp3 = exp(-326./(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
          TempUnit(1:nLons,1:nLats,1:nAlts)))

     ! In erg/cm3/s
!  Factor of two applied:
! ---------------------------------------------------------------------
     OCooling = 0.5* (1.69e-18*tmp2 + 4.59e-20*tmp3) * &
          (NDensityS(1:nLons,1:nLats,1:nAlts,iO_,iBlock)/1.0e6) / &
          (1.0 + 0.6*tmp2 + 0.2*tmp3)
!  Utilized without factor of two in Earth GITM:
! ---------------------------------------------------------------------
!    OCooling = (1.69e-18*tmp2 + 4.59e-20*tmp3) * &
!         (NDensityS(1:nLons,1:nLats,1:nAlts,iO_,iBlock)/1.0e6) / &
!         (1.0 + 0.6*tmp2 + 0.2*tmp3)
     ! In w/m3/3
     OCooling = OCooling/10.0
     ! In our special units:
     OCooling = OCooling/ TempUnit(1:nLons,1:nLats,1:nAlts) / &
          (Rho(1:nLons,1:nLats,1:nAlts,iBlock)*cp(:,:,1:nAlts,iBlock))

  else

     OCooling = 0.0

  endif

  RadCooling(1:nLons,1:nLats,1:nAlts,iBlock) = &
         RadCooling(1:nLons,1:nLats,1:nAlts,iBlock) + OCooling
  !\
  ! ---------------------------------------------------------------
  ! This calls the lower atmosphere radiation code
  ! Steven Nelli May 7, 2007
  ! ---------------------------------------------------------------
  !/
  ! --------------------------------------------------
  ! xLowAtmosRadRate = K/s
  ! xLowAtmosRadRate/TempUnit = <T>/s
  ! TempUnit = MeanMajorMass/kb = K*kg/J = K/(m/s)^2
  ! --------------------------------------------------
  !
  !     Mars GITM ground temperature (on its grid) -------------------

  !      grtemp(1:nLons,1:nLats,iBlock)=280.0 !240.0

  ! Calculating the solar flux at Mars

  ! SunPlanetDistance = sqrt(2.64236)
  do NW=1,L_NSPECTV
     SOL(nw) = SOLARF(NW)/(SunPlanetDistance**2)
     !SOL(nw) = SOLARF(NW)/(2.64236)
  end do

  !##############################################################
  ! MAIN LOOP, for MarsGITM horizontal grid (Latitude,Longitude)

  if( floor((tSimulation - dT)/DtLTERadiation) == &
       floor(tSimulation/DtLTERadiation)) then
     ! Do nothing
  else

     ! Calculate Ls
     deltat_h = iTimeArray(4)+iTimeArray(5)/60.0+iTimeArray(6)/3600.0

     deltat = 367*iTimeArray(1)-7*(iTimeArray(1)+(iTimeArray(2)+9)/12)/4 &
          +275*iTimeArray(2)/9+iTimeArray(3)-730531.5+deltat_h/24.0

     em = (19.3870+0.52402075*deltat)*pi/180.0

     alpha_pms = 270.3863+0.52403840*deltat

     ell_s = alpha_pms+(10.691+3.d-7*deltat)*dsin(em)+ &
          0.623*dsin(2.0*em)+0.050*dsin(3.0*em)+ &
          0.005*dsin(4.0*em)+0.0005*dsin(5.0*em)

     pbs = 0.0
     do i=1,7
        pbs = pbs+Ls_a(i)*dcos(pi/180.0*(0.985626*deltat/ &
             Ls_tau(i)+Ls_phi(i)))
     enddo
     ell_s = ell_s+pbs

     if (ell_s > 0.0) then
        ell_s = 360.0*(ell_s/360.0-floor(ell_s/360.0))
     else
        ell_s = 360.0+360.0*(ell_s/360.0-ceiling(ell_s/360.0))
     endif

!     	 print*, ell_s

     LowAtmosRadRate(1:nLons,1:nLats,1:nAlts,iBlock)=0.0

     SurfaceTemp(1:nLons,1:nLats,iBlock) = SurfaceTemp(1:nLons,1:nLats,iBlock)+&
          dSurfaceTemp(1:nLons,1:nLats,iBlock)*DtLTERadiation

     SubsurfaceTemp(1:nLons,1:nLats,iBlock) = SubsurfaceTemp(1:nLons,1:nLats,iBlock)+&
          dSubsurfaceTemp(1:nLons,1:nLats,iBlock)*DtLTERadiation

     do iLat = 1, nLats
        do iLon = 1, nLons

           !     Mars GITM ground temperature (on its grid) -------------------
           if (minval(Altitude_GB(:,:,0,iBlock)) .lt. 0 .or. UseTopography) then

              Temperature(iLon,iLat,0,iBlock)=SurfaceTemp(iLon,iLat,iBlock)/&
                   TempUnit(iLon,iLat,0)
              Temperature(iLon,iLat,-1,iBlock)=SurfaceTemp(iLon,iLat,iBlock)/&
                   TempUnit(iLon,iLat,-1)

           endif

           ! Determining the top of the lower atmosphere radiative zone
           do iAlt=1,nAlts
              if(pressure(iLon,iLat,iAlt,iBlock)*0.01.gt.prad) L_LAYERS = iAlt
           end do

           !        L_LAYERS = 20

           ! setting up the vertical fields
           L_LEVELS  = 2*L_LAYERS+3
           L_NLAYRAD  = L_LAYERS+1
           L_NLEVRAD  = L_LAYERS+2

           !BP: turn off until figure out what this is
           if (.false.) then
             call calc_lowatmosrad(iblock,iLat,iLon,L_LAYERS,L_LEVELS,&
                  L_NLAYRAD,L_NLEVRAD)
           endif
        end do !Latitude Loop
     end do    !Longitude Loop


!!$     xLowAtmosRadRate(1:nLons,1:nLats,1:nAlts,iBlock) = &
!!$          xLowAtmosRadRate(1:nLons,1:nLats,1:nAlts,iBlock)/&
!!$          (TempUnit(1:nLons,1:nLats,1:nAlts))

  endif

  !################# END MAIN COMPUTATIONAL LOOP#################

  call end_timing("calc_planet_sources")
  
  if (iTimeArray(3) .eq. 21 .and. iTimeArray(4) &
      .eq. 12 .and. iTimeArray(5) .eq. 1) then 
    call timescales
  endif
end subroutine calc_planet_sources

!---------------------------------------------------------+
!
!---------------------------------------------------------+

subroutine init_topography

  use ModGITM
  use ModInputs

  implicit None

!  real, intent(out) :: altzero2(nLons,nLats,nBlocks)

  integer, parameter :: nMOLALons = 1440 , nMOLALats = 720 !1/4 degree resolution

  real, dimension(nMolaLons,nMOLALats, 3) :: SurfaceAltitude
  integer :: ilon, ilat, iilon, iilat, jlon, jlat, iBlock
  real :: rlat, rlon,latfind,lonfind


  open(unit=iInputUnit_, file='DataIn/Mars_MOLA_topo.dat', action='read', status="old")
  if (iDebugLevel > 4) write(*,*) "=====> Reading Topography"

  do iLat = 1, nMOLALats
     do iLon = 1, nMOLALons

        read(iInputUnit_,*) SurfaceAltitude(iLon,iLat,iNorth_), &
             SurfaceAltitude(iLon,iLat,iEast_), &
             SurfaceAltitude(iLon,iLat,iUp_)

     enddo
  enddo
  close(iInputUnit_)

  do iBlock = 1, nBlocks
     do iLon = 1, nLons
        do iLat = 1, nLats

           LonFind = Longitude(iLon,iBlock)*180/pi
           LatFind = latitude(iLat,iBlock)*180/pi

           do jLon = 1, nMOLALons-1

              if (SurfaceAltitude(jLon,1,iEast_) <= LonFind .and. &
                   SurfaceAltitude(jLon+1,1,iEast_) >= LonFind) then
                 iiLon = jLon
                rLon = 1.0 - (LonFind -  SurfaceAltitude(jLon,1,iEast_))/ &
                      (SurfaceAltitude(jLon+1,1,iEast_)-SurfaceAltitude(jLon,1,iEast_))

              endif
           enddo

           do jLat = 1, nMOLALats-1

              if (SurfaceAltitude(1,jLat,iNorth_) <= LatFind .and. &
                   SurfaceAltitude(1,jLat+1,iNorth_) >= LatFind) then
                 iiLat = jLat
                 rLat = 1.0 - (LatFind -  SurfaceAltitude(1,jLat,iNorth_))/ &
                      (SurfaceAltitude(1,jLat+1,iNorth_)-SurfaceAltitude(1,jLat,iNorth_))
              endif
           enddo

           altzero(iLon,iLat,iBlock) =  (rLon)*(rLat)*SurfaceAltitude(iiLon,iiLat,iUp_) + &
                (1-rLon)*(  rLat)*SurfaceAltitude(iiLon +1,iiLat,iUp_) + &
                (rLon)*(1-rLat)*SurfaceAltitude(iiLon,iiLat+1,iUp_) + &
                (1-rLon)*(1-rLat)*SurfaceAltitude(iiLon+1,iiLat+1,iUp_)

        enddo
     enddo

  enddo

end subroutine init_topography
!---------------------------------------------------------+
!
!---------------------------------------------------------+
!------------------------------------------------------------------
subroutine nlte_tcool(iBlock)

  !  ---------------------------------------------------------------------
  !  Purpose: to calculate NLTE CO2 15-micron cooling in Mars upper atmosphere
  !           (above 80 km, where NLTE physics becomes important)
  !  Source: nltecool.F from mtgcm16
  !  Coders: M. Lopez-Valverde (2012)
  !          S. W. Bougher (2017-MGITM)
  !  Inputs: from MarsGITM  (Pressure, Temperature, VMRCO2, VMRO, VMRN2, VMRCO)
  !          all on MarsGITM grid (nLons,nLats,nAlts)
  !  Output: COOLTOT (K/sec) and RadCoolingRate(1:nLons,1:nLats,1:nAlts,iBlock)
  !  =======================================================================
  !  From Miguel Lopez-Valverde: Main code for modification and calling
  !
  !    subroutine nlte_tcool(ngridgcm,n_gcm,
  !  $     p_gcm, t_gcm, z_gcm,
  !  $     co2vmr_gcm, n2vmr_gcm, covmr_gcm, o3pvmr_gcm,
  !  $     q15umco2_gcm , ierr, varerr)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Fast scheme for NLTE cooling rates at 15um by CO2 in a Martian GCM !
  !                 Version dlvr11_03. 2012.                           !
  ! Software written and provided by IAA/CSIC, Granada, Spain,         !
  ! under ESA contract "Mars Climate Database and Physical Models"     !
  ! Person of contact: Miguel Angel Lopez Valverde  valverde@iaa.es    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! November 2017   Adapt to the MGITM model (3D version); S. Bougher
  !                                                        & F. G.Galindo
  !  -----------------------------------------------------------------

  use ModInputs
  use ModSources, only: RadCoolingRate
  use ModPlanet
  use ModGITM
  use ModUserGITM
  use ModConstants, only:  Boltzmanns_Constant, Speed_Light, &
       Planck_Constant, Avogadros_Number
  use ModIndicesInterfaces

  implicit none

  ! Common Blocks (from LMD-MGCM) !   Added to ModPlanet.f90
  ! include 'nlte_paramdef.h'
  ! include 'nlte_commons.h'

  ! INPUTS:----------------------------------------------------------

  integer,intent(in):: iBlock

  !  -----------------------------------------------------------------
  ! Input variables

  integer, parameter :: n_gcm = nAlts ! no. of atmospheric layers
  integer, parameter :: nlon = nLons ! no. of long gridpoints
  integer, parameter :: nlat = nLats ! no. of latitude gridpoints
  integer, parameter :: ngridgcm = nlat*nlon ! no. of horizontal gridpoints

  !logical :: found

  !  -----------------------------------------------------------------

  ! Local Variables and Constants

  integer :: iLat, iLon, iAlt
  integer :: i,ig, l, indice, nl_cts_real, nzy_cts_real
  real*8 :: q15umco2_nltot(nltot), zld(nltot)
  real*8 :: hr110CTS(nl_cts)
  real :: xx,factor
  real :: auxgcm(n_gcm)
  real*8 :: auxgcmd(n_gcm), aux2gcmd(n_gcm)
  real :: zmin_gcm
  integer :: ierr
  real*8 :: varerr

  real :: refCooling = 0.0
  integer :: jAlt
  real :: alpha
  real :: pMerge = 0.5 !Pa


  !  -----------------------------------------------------------------
  ! Repository for MGITM fields into 2-D arrays (not needed)
  ! real,dimension(ngridgcm,n_gcm) :: p_gcm, t_gcm, z_gcm, co2vmr_gcm, &
  !    n2vmr_gcm, covmr_gcm, o3pvmr_gcm, q15umco2_gcm
  !
  !  -----------------------------------------------------------------
  ! Repository for MGITM fields into 1-D arrays

  real,dimension(n_gcm) :: p_ig,z_ig,t_ig, &
      co2_ig,n2_ig,co_ig,o3p_ig,mmean_ig, cpnew_ig

  !  -----------------------------------------------------------------
  ! Internal fields recast from MarsGITM

  real,dimension(1:nLons,1:nLats,1:nAlts) ::    &
       TN2,P,vmrco2,vmro,vmrn2,vmrco,q15umco2_gcm,mnd,mmean,cpm,zht

  !  Mars GITM real temperature (on its grid)  ----------------------

  TN2(1:nLons,1:nLats,1:nAlts) = &
       Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*  &
       TempUnit(1:nLons,1:nLats,1:nAlts)

  !     Mars GITM pressure in pascals  (on its grid) -------------------
  !     (convert to ATM units for nlte_tcool subroutine)
  !     (1.0 mb = 100.0 Pa;  1013.25 mb = 1.01325e+05 Pa = 1.0 ATM)
  !     (1.0 Pa = 9.869e-6 ATM)

  P(1:nLons,1:nLats,1:nAlts) = &
       Pressure(1:nLons,1:nLats,1:nAlts,iBlock)*9.869e-06

  !     Mars GITM heights in meters  (on its grid) -------------------

  zht(1:nLons,1:nLats,1:nAlts) = &
        Altitude_GB(1:nLons,1:nLats,1:nAlts,iBlock)

  !     Mars GITM cp in J/kg/K (on its grid) -------------------
  !     (same units as LMD-MGCM cooling code, and will vary with altitude)
  !     (cpco2 = 8.4e+06 erg/gm/K = 0.84 J/gm/K = 8.4e-04 J/kg/K)

  cpm(1:nLons,1:nLats,1:nAlts) = &
       cp(1:nLons,1:nLats,1:nAlts,iBlock)

  !    Volume Mixing Ratio amd MWT Calculations

  mnd(1:nLons,1:nLats,1:nAlts) = &
       NDensity(1:nLons,1:nLats,1:nAlts,iBlock)+1.0
  vmro(1:nLons,1:nLats,1:nAlts)  = &
       NdensityS(1:nLons,1:nLats,1:nAlts,iO_,iBlock)&
       /mnd(1:nLons,1:nLats,1:nAlts)
  vmrco(1:nLons,1:nLats,1:nAlts) = &
       NdensityS(1:nLons,1:nLats,1:nAlts,iCO_,iBlock)&
       /mnd(1:nLons,1:nLats,1:nAlts)
  vmrn2(1:nLons,1:nLats,1:nAlts) = &
       NdensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock)&
       /mnd(1:nLons,1:nLats,1:nAlts)
  vmrco2(1:nLons,1:nLats,1:nAlts) = &
       NdensityS(1:nLons,1:nLats,1:nAlts,iCO2_,iBlock)&
       /mnd(1:nLons,1:nLats,1:nAlts)

!  mean molecular mass (gm/molecule) = mean molecular weight(gm/mole)/
!                                      Avogardos number(molecules/mole)

!  (gm/molecule)
! mmean(1:nLons,1:nLats,1:nAlts) = &
!      (vmro(1:nLons,1:nLats,1:nAlts)*16.+ &
!      (vmrco(1:nLons,1:nLats,1:nAlts)+vmrn2(1:nLons,1:nLats,1:nAlts))&
!      *28. + vmrco2(1:nLons,1:nLats,1:nAlts)*44.)/  &
!      Avogadros_Number
!  (gm/mole) = mmwt
  mmean(1:nLons,1:nLats,1:nAlts) = &
       (vmro(1:nLons,1:nLats,1:nAlts)*16.+ &
       (vmrco(1:nLons,1:nLats,1:nAlts)+vmrn2(1:nLons,1:nLats,1:nAlts))&
       *28. + vmrco2(1:nLons,1:nLats,1:nAlts)*44.)

  !  -----------------------------------------------------------------
  !
  ! ****  Initialize to zero
  q15umco2_gcm(1:nLons,1:nLats,1:nAlts) = 0.0

  !-------------------------------------------------------------

  ! MAIN LOOP, for MarsGITM horizontal grid (Latitude,Longitude)
  ! (Corresponds to the ig = 1,ngridgcm loop from MLV)

  do iLat = 1, nlat
     do iLon = 1, nlon
         ierr = 0
         nl_cts_real = 0
         nzy_cts_real = 0
  ! *****  populate 1-D arrays(n_gcm) for calc_radcooling for MarsGITM
        do  iAlt=1,n_gcm
            p_ig(iALT)= P(iLon,iLat,iAlt)
            t_ig(iALT)= TN2(iLon,iLat,iAlt)
            co2_ig(iALT)=vmrco2(iLon,iLat,iAlt)
            n2_ig(iALT)=vmrn2(iLon,iLat,iAlt)
            o3p_ig(iALT)=vmro(iLon,iLat,iAlt)
            co_ig(iALT)=vmrco(iLon,iLat,iAlt)
            z_ig(iALT)=zht(iLon,iLat,iAlt)/1000.
            mmean_ig(iALT)=mmean(iLon,iLat,iAlt)
            cpnew_ig(iALT)=cpm(iLon,iLat,iAlt)
        enddo


  ! ****** From GCM's grid to NLTE's grid

         call NLTEdlvr11_ZGRID (n_gcm,  &
              p_ig, t_ig, z_ig,         &
              co2_ig, n2_ig, co_ig, o3p_ig, &
              mmean_ig,cpnew_ig, &
              nl_cts_real, nzy_cts_real )

  ! *****  Isotopic Tstar & VC at the NLTE grid
         call interdp_ESCTVCISO

  ! *****   Tstar para NLTE-CTS
  !      call MZESC110 ( ig,nl_cts_real, nzy_cts_real,ierr,varerr )
         call MZESC110 (nl_cts_real, nzy_cts_real,ierr,varerr )
         if (ierr .gt. 0) call ERRORS (ierr,varerr)

         ! 626FB C.M.
         call leetvt
         c110(1:nl,1:nl)=0.d0
!         call zerom (c110, nl)
         call zero2v (vc110,taustar11, nl)
         call MZTUD110 ( ierr, varerr )
         if (ierr .gt. 0) call ERRORS (ierr,varerr)

         input_cza = 0
         call NLTEdlvr11_CZALU(ierr,varerr)
         if (ierr .gt. 0) call ERRORS (ierr,varerr)

         input_cza = 1
         call NLTEdlvr11_CZALU(ierr,varerr)
         if (ierr .gt. 0) call ERRORS (ierr,varerr)

                                !  call NLTEdlvr11_FB626CTS
                                ! Falta un merging del hr110CTS con el HR110


                                ! NLTE-CTS
         call NLTEdlvr11_FB626CTS ( hr110CTS , nl_cts_real )



                                ! total TCR
         do i = 1, nl
            q15umco2_nltot(i) =hr110(i) + hr210(i) + hr310(i) + hr410(i) &
                 + hr121(i)
         enddo


                                ! Merging con / actualizacion del HR_total
                                !   Eliminamos el ultimo pto de hrTotal, y en el penultimo
                                !   (que coincide con i=1 en el grid nl_cts)
                                !   hacemos la media entre hrTotal y hr110CTS :
         i=nl-1
         q15umco2_nltot(i) = 0.5d0*( q15umco2_nltot(i) + hr110CTS(1) )
         do i=2,nl_cts_real
            indice = (nl-2) + i
            q15umco2_nltot(indice) = hr110CTS(i)
         enddo
         do i=nl_cts_real+1,nl_cts
            indice = (nl-2) + i
	    q15umco2_nltot(indice) = 0.0d0
         enddo

                                ! Interpol to original Pgrid
                                !
                                ! Primero, la parte conocida ([1,nl_cts_real])
         do i=1,nl
            zld(i) = - dble ( alog(pl(i)) )
                                !write (*,*) i, zld(i), q15umco2_nltot(i)
         enddo
         do i=3,nl_cts_real
            indice = (nl-2) + i
            zld(indice) = - dble ( alog(pl_cts(i)) )
                                !write (*,*) indice, zld(indice), q15umco2_nltot(indice)
         enddo
                                ! En caso que nl_cts_real < nl_cts , extrapolo el grid alegremente
         factor = pl_cts(nl_cts_real)/pl_cts(nl_cts_real-1)
         xx = pl_cts(nl_cts_real)
         do i=nl_cts_real+1,nl_cts
            indice = (nl-2) + i
            xx = xx * factor
            zld(indice) = - dble ( alog(xx) )
         enddo

         do i=1,n_gcm
            auxgcmd(i) = - dble( alog(p_ig(i)))
         enddo
!         call zerov( aux2gcmd, n_gcm )
         aux2gcmd(1:n_gcm)=0.d0
         call interdp_limits (aux2gcmd, auxgcmd, n_gcm,  &
               jlowerboundary,jtopCTS,  &
              q15umco2_nltot, zld, nltot, 1,  nltot, 1)

                                ! Smoothing
         call suaviza ( aux2gcmd, n_gcm, 1, auxgcmd )

         do i=1,n_gcm
            q15umco2_gcm(iLon,iLat,i) = sngl( aux2gcmd(i) )
         enddo

     enddo  !-------- END OF MAIN LONGITUDE LOOP
  enddo  !-------- END OF MAIN LATITUDE LOOP

  !-------------------------------------------------------------

! CO2 Cooling Rate  (K/Earth_day):  Use Positive sign from Loop above
!    Retain dynamical and radiative heat balance terms in K/Earth_day
!    Apply negative sign to q15umco2_gcm to make work inside MGITM
!    Convert to K/sec for Internal Usage inside M-GITM code
  RadCoolingRate(1:nLons,1:nLats,1:nAlts,iBlock) = &
       -q15umco2_gcm(1:nLons,1:nLats,1:nAlts)/86400.

!  Diagnostics
     !-----------------------------------------------------------------
     ! S. W. BOUGHER defined:  11-11-01 UserData1D
     ! S. W. BOUGHER defined:  11-11-01 UserData3D
     ! Pressure (ubar units for MTGCM compatibility)
     !-----------------------------------------------------------------
!    UserData3D(:,:,:,8,iBlock) = 0.0
!    UserData3D(1:nLons, 1:nLats, 1:nAlts, 8, iBlock) =  &
!              Pressure(1:nLons, 1:nLats, 1:nAlts, iBlock)*10.0
!    UserData3D(:,:,:,9,iBlock) = 0.0
!    UserData3D(1:nLons, 1:nLats, 1:nAlts, 9, iBlock) =  &
!              Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*  &
!              TempUnit(1:nLons,1:nLats,1:nAlts)
!    UserData3D(:,:,:,10,iBlock) = 0.0
!    UserData3D(1:nLons, 1:nLats, 1:nAlts, 10, iBlock) =  &
!              vmro(1:nLons,1:nLats,1:nAlts)
!    UserData3D(:,:,:,11,iBlock) = 0.0
!    UserData3D(1:nLons, 1:nLats, 1:nAlts, 11, iBlock) =  &
!              vmrco2(1:nLons,1:nLats,1:nAlts)
!
     !-----------------------------------------------------------------
     UserData1D(1,1,:,4) = 0.0
     UserData1D(1, 1, 1:nAlts, 4) =    &
               RadCoolingRate(1,1,1:nAlts,iBlock)
     UserData1D(1,1,:,8) = 0.0
     UserData1D(1, 1, 1:nAlts, 8) =  Pressure(1, 1, 1:nAlts, iBlock)*10.0
     UserData1D(1,1,:,9) = 0.0
     UserData1D(1, 1, 1:nAlts, 9) =   &
               Temperature(1,1,1:nAlts,iBlock)*  &
               TempUnit(1,1,1:nAlts)
     UserData1D(1,1,:,10) = 0.0
     UserData1D(1, 1, 1:nAlts, 10) = vmro(1,1,1:nAlts)
     UserData1D(1,1,:,11) = 0.0
     UserData1D(1, 1, 1:nAlts, 11) =  vmrco2(1,1,1:nAlts)

  !-------------------------------------------------------------

  !-------------------------------------------------------------
  
  !BP (7/7/2020)
  !Hotfix to Venus RadCooling that produces positive "cooling rates (~250,000 K/day)
  !below 85 km. Making it zero until we implement a better cooling scheme in the 
  !LTE region like Newtonian cooling or something else
  do iAlt = 1, nAlts
    if (Altitude_GB(1,1,iAlt,iBlock)/1000.0 < 85.0) then
      RadCoolingRate(1:nLons,1:nLats,iAlt,iBlock) = 0.0
    endif
  enddo

  !Newtonian cooling
  do iAlt = 1, nAlts
    do iLon = 1, nLons
      do iLat = 1, nLats
        if (RadCoolingRate(iLon,iLat,iAlt,iBlock) .le. 0.0) then
          RadCoolingRate(iLon,iLat,iAlt,iBlock) = 0.0
        endif
      !Altitudes above 0.01 Pa, NLTE effects dominate                                   
      !Blend below these altitudes according to the Gilli 2016 paper 
        if (Pressure(iLon,iLat,iAlt,iBlock) > 0.01) then
          !Find the cooling amount
          do jAlt = 1,31
            if (Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0 .le. referenceAltitude2(jAlt+1) .and. &
                Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0 .ge. referenceAltitude2(jAlt)) then
              refCooling = -newtonianCoolingRate(jAlt)
              exit
            endif
          enddo
            
          !Scale the cooling appropriately
          do jAlt = 1, 56
            if (Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0 .le. referenceAltitude1(jAlt+1) .and. &
                Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0 .ge. referenceAltitude1(jAlt)) then
              !Altitudes above 0.01 Pa, NLTE effects dominate
              !Blend below these altitudes according to the Gilli 2016 paper

              alpha = 1/(1 + Pressure(iLon,iLat,iAlt,iBlock)/pMerge)**4
              RadCoolingRate(iLon,iLat,iAlt,iBlock) = &
                  alpha*RadCoolingRate(iLon,iLat,iAlt,iBlock) + &
                  (1 - alpha)*refCooling*Temperature(iLon,iLat,iAlt,iBlock)*&
                  TempUnit(iLon,iLat,iAlt)/ReferenceTemperature(jAlt)
              exit
            endif
          enddo
        endif
      enddo
    enddo
  enddo

end subroutine nlte_tcool
!---------------------------------------------------------

subroutine init_isochem
  use ModPlanet, only: ialtminiono
  use ModInputs, only: altminiono
  use ModGITM, only: Altitude_GB
  use ModSizeGITM


  implicit none

  integer :: iBlock, iLon, iLat, iAlt

  iAltMinIono = 1

  do iBlock = 1, nBlocks
     do iLon = 1, nLons
        do iLat = 1, nLats
           do ialt = 1, nAlts
              if (Altitude_GB(iLon,iLat,ialt,iBlock)/1000.0 .le. AltMinIono) &
                   iAltMinIono(iLon,iLat,iBlock) = iAlt
              if (ialtminiono(ilon,ilat,iblock) .lt. 1)  ialtminiono(ilon,ilat,iblock) = 1
           enddo
        enddo
     enddo
  enddo

end subroutine init_isochem

   !---------------------------------------------------------------------
   ! Calculate Eddy Diffusion Coefficient
   !---------------------------------------------------------------------

   subroutine calc_eddy_diffusion_coefficient(iBlock)

     use ModSizeGITM
     use ModGITM, only: pressure, NDensity
     use ModInputs, only: EddyDiffusionPressure0,EddyDiffusionPressure1, &
          EddyDiffusionCoef
     use ModSources, only: KappaEddyDiffusion

     implicit none

     integer, intent(in) :: iBlock
     integer :: iAlt, iLat, iLon
     integer :: First
     real :: PEddyMax
     real :: NEddyMax(nLons,nLats)
     real :: EddyProfile(nLons,nLats,-1:nAlts+2)

     real :: KMax
     real :: KMin

     KappaEddyDiffusion(:,:,:,iBlock) = 0.0
! KStandard
!    KMax = 1000.0
!    KMin = 500.0

! KHigh
!     KMax = 2000.0
!     KMin = 500.0

! KModerate
!     KMax = 1500.0
      KMax = 2000.0
      KMin = 500.0

! KLow
!    KMax = 1200.0
!    KMin = 500.0

     ! \
     ! First, find the altitude level corresponding to the asymptotic
     ! upper bound for Eddy Diffusion.
     ! Call this upper limit, NEddyMax

     ! This upper limit is set to 1.26e-09 bars
     ! conversion to pascals -> 1 bar = 1e+05 pascals
     ! Thus, PEddyMax -> 1.26e-04 pascals


     PEddyMax = 1.26e-04  ! Pascals (SI Units)

     do iLat = 1, nLats
        do iLon = 1, nLons

           First = 0

           do iAlt = 1, nAlts

              if (Pressure(iLon,iLat,iAlt,iBlock) > PEddyMax) then
                 cycle
              else
                 if (First == 0) then
                    NEddyMax(iLon,iLat) = NDensity(iLon,iLat,iAlt,iBlock)
                    First = 1
                 endif
              endif

           enddo

        enddo
     enddo

     ! Now we have all the trigger densities as a function of Longitude and Latitude

     ! Next, extend the profile downward as 1/sqrt(N) and put a lower bound of 100 m^2/s

     do iAlt = -1, nAlts+2
        do iLat = 1, nLats
           do iLon = 1, nLons

              ! Krasnopolsky et al. [2005] Helium modeling
!              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) =  &
!                   (1.0e-04)*(1.8e+13)/sqrt( (1.0e-06)*NDensity(iLon,iLat,iAlt,iBlock))

              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) =  &
                   KMax* sqrt( NEddyMax(iLon,iLat) / NDensity(iLon,iLat,iAlt,iBlock))
              !
              !! \
              !! This gives an upper bound of Kmax
              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = &
                   min(KMax, KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) )

             KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = &
                  max(KMin, KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) )
              !
              !! This gives an lower bound of Kmin
!             KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = &
!                  max(100.0e+02, KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) )
             !

           enddo
        enddo

     enddo


   end subroutine calc_eddy_diffusion_coefficient


   ! ----------------------------------------------------------------------

   subroutine calc_gw(iBlock)

     use ModInputs
     use ModSources
     use ModPlanet
     use ModGITM
     use ModConstants
     use ModTime

     integer, parameter :: NLAY=nAlts
!    real, parameter :: KAP = 0.000025
     real, parameter :: KAP = 0.0000025
     !  declare  GWAccel  for compilation testing purposes only; this will
     !  be commented out when the  Mod..  files are available

!    real :: u(nlay),v(nlay),t(nlay),rhogw(nlay),rholev(nlay), &
!         geot(nlay),pz(nlay),du(nlay),uav(nlay),plev(nlay), &
!         tavg(nlay),theta(nlay),thtavg(nlay), &
!         dtheta(nlay),stress(0:nlay),dstress(nlay), &
!         dudt(nlay),cpco2
     real :: u(nlay),v(nlay),tgw(nlay),rhogw(nlay),rholev(nlay), &
          geot(nlay),pz(nlay),du(nlay),uav(nlay),plev(nlay), &
          tavg(nlay),theta(nlay),thtavg(nlay), &
          dtheta(nlay),stress(0:nlay),dstress(nlay), &
          dudt(nlay),cpco2

     real :: N(nlay),dh(nlay),dhsat(nlay),H(nlay), &
          dz(nlay),Ri(nlay),Rimin(nlay),eps(nlay)
     real :: Nsfc,k,C, VAR,TS,xRi,rhosfc,psfc,RCO2,gravv,gravsfc

     integer :: ilon,ilat,ialt,nlgw

     PSFC = 400.0
     RCO2 = 189.0
     VAR  = 400.0**2
!    MKS Units here
!    cpco2=3.5*RGAS*AMU/Mass(iCO2_)
     cpco2=3.5*RGAS*1.0E-07*AMU/Mass(iCO2_)*1000.

     !   MAIN DO LOOPS... OVER LONGITUDE AND LATITUDE
     GWAccel(:,:,:,:) = 0.0

     if (CurrentTime-StartTime .lt. 3600) then
        return
     endif
     do ilon=1,nLons
        do ilat=1,nLats

           gravsfc = -Gravity_GB(ilon,ilat,1,iBlock)


!   FILL LOCAL VERTICAL VECTORS FROM GITM ARRAY VALUES

           U = Velocity(iLon,iLat,1:nAlts,iEast_,iBlock)+1.0E-4
           rhogw = Rho(ilon,ilat,1:nAlts,iBlock) ! Is this correct ??

           geot = Altitude_GB(ilon,ilat,1:nAlts,iBlock) ! Is this correct ?
           TGW = Temperature(ilon,ilat,1:nAlts,iBlock)*TempUnit(ilon,ilat,1:nAlts)
           pz = rhogw * RCO2 * tgw
           THETA = tgw * (psfc/pz)**(RCO2/CpCO2)

           eps = 0.0
           dhsat = 0.0
           dh = 0.0
           Ri=0.0
           Rimin = 0.0
           stress = 0.0
           dstress = 0.0
           dudt = 0.0

           do nlgw=1,nlay-1
              !    Calculate the geometric distance between adjacent layer midpoints;
              !   units are  meters

              dz(nlgw) = geot(nlgw+1) - geot(nlgw)

              !    Calculate the average wind speed at the top of each model layer
              !   (but not the top layer), using a simple averaging of the zonal wind
              !   speeds at the two neighboring layer midpoints

              uav(nlgw) = 0.5 * (u(nlgw)+u(nlgw+1))

              !    Calculate the potential temperature at layer boundaries, via a simple
              !   averaging of the potential temperatures at the neighboring layer midpoints
              !    This can be imrpoved by using the layer boundary pressure (known
              !   via  psfc  and  sigma) to directly calculate the pressure term
              !   in the potential temperature calculation.. need a temperature at
              !   the layer boundary

              thtavg(nlgw) = 0.5 * (theta(nlgw)+theta(nlgw+1))

              !    Calculate the average temperature at layer boundaries, via a simple
              !   averaging of the temperatures at the neighboring layer midpoints

              tavg(nlgw) = 0.5 * (tgw(nlgw)+tgw(nlgw+1))

              !   Calculate the difference in zonal wind speed across a layer boundary

              du(nlgw) =   u(nlgw+1) - u(nlgw)
              !
              !   Calculate the density values at layer boundaries using the average of
              !  the pressures at the two neighboring layer midpoints and the average
              !  temperature at that same layer boundary;  this is a poor way of calculating
              !  this average density and will be improved upon using independently
              !  determined values of pressure and temperature at layer boundaries

              plev(nlgw) = 0.5*(pz(nlgw)+pz(nlgw+1))/2.0
              rholev(nlgw) = 0.5*(pz(nlgw)+pz(nlgw+1))/(RCO2*Tavg(nlgw))

              !  Difference in potential temperature values across layer boundaries,
              ! calcultes by calculating the difference between the two neighboring
              ! layer midpoint potential temperatures;  this should be  (nlgw-1)-(nlgw)

              dtheta(nlgw) = theta(nlgw+1) - theta(nlgw)

           enddo



           !      Calculate the surface drag magnitude, using lowest layer midpoint
           !     density and wind speed values, and N (Brunt-Vaisala frequency) value
           !     at the interface between the bottom two atmosphere layers

           rhosfc = psfc / (RCO2 * TGW(1))

          Nsfc = (gravsfc/theta(1))*(dtheta(1)/dz(1))**0.5

           TS = KAP * rhosfc * Nsfc * (U(1)) * VAR
           !write(88,*) ts,kap,rhosfc,nsfc,theta(1), &
            !    dtheta(1),dz(1),u(1)
           !print *, TS, KAP,Nsfc,U(1),VAR

           !     Set STRESS(NLAY+1) value equal to the surface drag value
           stress(0) = TS

          ! write(14,99) 0, stress(0), psfc, 1.0*msol, 1.0*mbin

           !  THIS IS THE PRIMARY  LOOP  OVER VERTICAL LAYERS, WITHIN
           !  WHICH WAVE VERTICAL PROPAGATION AND BREAKING ARE DIAGNOSED


           do  22 nlgw=1,nlay-1

              gravv = -Gravity_GB(ilon,ilat,nlgw,iBlock)

              !    DO NOT IMPLEMENT WAVE PROPAGATION (more properly, 'BREAKING')
              !    WITHIN THE BOTTOM ONE MODEL LAYER(S)

              IF(NLGW.LE.1) THEN
                 STRESS(NLGW) = STRESS(NLGW-1)
                 DSTRESS(NLGW) = 0.0

                 !   The  ELSE below accounts for model layers above the bottom three
                 !  layers... in these upper layers wave breaking is permitted
              ELSE

                 !  Determine if a CRITICAL LEVEL situation is occurring, and if it is
                 ! deposit all of the stress within the offending layer (the layer within
                 ! which the sign of the zonal wind speed is opposite the sign of the
                 ! stress magnitude at the bottom of that layer)

                 !  DIAGNOSE CRITICAL LEVEL OCCURRENCE BASED UPON EITHER LAYER MIDPOINT
                 !  ZONAL WIND (u(nlgw)) OR AVERAGE ZONAL WIND AT THE TOP OF THE LAYER
                 !  (uav(nlgw)) BEING OPPOSITE IN SIGN TO THE STRESS VALUE AT THE BOTTOM
                 !  OF THE LAYER

                 if(u(nlgw)*stress(nlgw-1) .lt. 0.0 .or. &
                      uav(nlgw)*stress(nlgw-1) .lt. 0.0) then
                    !     stress at the layer top is set equal to zero
                    stress(nlgw) = 0.0
                    !     dstress at the layer top is set equal to the stress value
                    !   at the bottom of the layer

                    dstress(nlgw) = stress(nlgw-1)

                    !   for this layer within which a CRITICAL LEVEL has arisen and
                    !   wave momentum has been deposited, calculate the wind acceleration
                    !   (units of m/s/s) by dividing the  dstress value by the mass (kg)
                    !   in the layer; a NEGATIVE sign indicates a westward acceleration,
                    !   and a positive value indicates an eastward accaleration

                    dudt(nlgw)=-dstress(nlgw)/((plev(nlgw)-plev(nlgw+1))/gravv)
                    GWAccel(ilon,ilat,nlgw,iEast_) = dudt(nlgw) * Dt

                    !       set dh(nlgw) and dhsat(nlgw) values to -9.99 in this CRITICAL LAYER
                    !     occurrence to easily help identify such layers
                    dh(nlgw) = -9.99
                    dhsat(nlgw) = -9.99
                    !     the  goto  statement below sends the code to the next layer,
                    !   but in actuality there is no further need to continue upwards
                    !   in this column since there is no wave energy to consider...
                    !   THIS CAN BE IMPROVED AND MADE MORE EFFICIENT

                    !   WRITE to unit 17 to permit statistical look at wave breaking level

!                    write(17,*) msol,mbin,nlgw,dudt(nlgw),1.25+(nlgw-1)*2.5
                    !    below it would be more appropriate to  goto  statement  6  since
                    !    at this point this column is completed
                    goto 22
                 endif
                    !    the above  endif  ends the treatment of CRITICAL LEVELS

                    !  Calculate the Brunt-Vaisala frequency at the top boundary of
                    ! layer  nlgw

                    N(nlgw) = (gravv/thtavg(nlgw)) * &
                         (theta(nlgw+1)-theta(nlgw))/ &
                         dz(nlgw)

                    !  THERE IS CURRENTLY NO 'CHECK' FOR NEGATIVE N(nlgw) VALUES, WHICH
                    !  would arise in the presence of a superadiabatic lapse rate.. which
                    !  will not occur in teh Ames model output but could, I believe, arise
                    !  the no-hydrostatic GITM..

                    N(nlgw) = sqrt(N(nlgw))
!                    write(88,*) nlgw,N(nlgw),theta(nlgw+1),theta(nlgw),dz(nlgw)

                    !  Estimate the isentropic vertical displacement at the top of layer  nlgw
                    ! using the calculated stress value at the value at the lnect lower layer
                    ! boundary and the density, Brunt-Vaisala frequency and average zonal wind
                    ! speed at the layer boundary of interest

                    !    use ABS(stress(nlgw+1)) and ABS(uav(nlgw))

                    dh(nlgw) = sqrt(abs(stress(nlgw-1)/(rhogw(nlgw)*kap*N(nlgw)*uav(nlgw))))

!                    print * , '**',stress(nlgw-1),rhogw(nlgw),kap,n(nlgw),uav(nlgw),nlgw

                    !  Calculate the Richardson Number at the top boundary of Layer  nlgw

                    !   MAKE Ri  and Rimin  vectors with NLAY elements

                    Ri(nlgw)= (gravv/thtavg(nlgw))* (dtheta(nlgw)/dz(nlgw))/ &
                         ((du(nlgw)/dz(nlgw))**2)

!                    write(33,*) nlgw, Ri(nlgw),gravv,thtavg(nlgw),dtheta(nlgw), &
!                         dz(nlgw),du(nlgw),dz(nlgw)

                    !  Calculate the 'minimum' Richardson Number at the top of Layer  nlgw
                    ! using the Brunt-Vaisala frequency, estimated isentropic vertical displacement,
                    ! average zonal wind speed, and Richardson Number already calculated at that
                    ! layer boundary of interest

                    !    implement  ABS(uav(nlgw))
                    Rimin(nlgw) = Ri(nlgw) * (1.0-(N(nlgw)*dh(nlgw)/ABS(uav(nlgw)))) / &
                         (1.0 + ((Ri(nlgw)**0.5)*N(nlgw)*dh(nlgw)/abs(uav(nlgw))))**2.0

                    !  If Rimin is greater than 0.25, than there is no wave breaking within
                    ! layer  nlgw  and thus the stress value at the top of Layer  nlgw  is the same
                    ! as the stress value at the bottom of Layre  nlgw  and the dstress value in
                    ! Layer  nlgw  is zero, and there is no wind acceleration within that layer

                    if(Rimin(nlgw).gt.0.25) then
                       stress(nlgw) = stress(nlgw-1)
                       dstress(nlgw) = 0.0
                    else

                       !   if Rimin is less than 0.25, than wave breaking is occurring and some stress
                       !  (acceleration.. really deceleration for westerly flow) is being applied to
                       !  the zonal wind in Layer  nlgw;

                       !  with the below inclusion of xRi I am NOT permitting
                       !  the value of  Ri  to be less than 0.25 when I calculate
                       !  the value of  epsILON... I am not sure if this is the
                       !  proper way to be proceeding here

                       xRi = max (Ri(nlgw),0.25000)

                       eps(nlgw)=(xRi**(-0.5))*(1.0+2.0*xRi**0.5)* &
                            ((2*(xRi**0.250)*(1.0+2.0*(xRi**0.50))**(-0.50))-1.0)
!                       write(15,99) nlgw,xRi,Ri(nlgw),Rimin(nlgw),eps(nlgw),2*(xRi**0.250), &
!                            (1.0+2.0*(xRi**0.50))**(-0.50),  &
!                            ((2*(xRi**0.250)*(1.0+2.0*(xRi**0.50))**(-0.50))-1.0)

                       !  eps  will not be less than zero, since if it did have a
                       !  negative value it would result in a negative  dh  value

                       eps(nlgw) = max(eps(nlgw),0.0)

                       !  Calculate a more representative value of the isentropic displacement
                       ! at the top of Layer  nlgw, and then use that value to calculate a new value
                       ! of the stress at the top of Layer  nlgw, and then a value (dstress) that
                       ! is the change in the stress across the layer.  it is the magnitude of this
                       ! dstress value that determines the magnitude of the wind acceleration that
                       ! results

                       !  implement  ABS(uav(nlgw))

                       dhsat(nlgw)= eps(nlgw) / (N(nlgw)/ABS(uav(nlgw)))

                       !   Determine the stress value at a layer's top boundary if wave breaking
                       !  is diagnosed within that layer

                       stress(nlgw) =eps(nlgw)*eps(nlgw)*kap*rhogw(nlgw)*(uav(nlgw)**3)/ &
                            N(nlgw)

                       !  If the stress at the top of layer  nlgw  is greater than the stress value
                       ! at the bottom of Layer  nlgw, then the stress value at the top of Layer nlgw
                       ! is set equal to zero and the  dstress value for layer  nlgw  is set equal
                       ! to the stress value at the bottom of Layer  nlgw

                       !  IS THIS CORRECT??
                       if(stress(nlgw).gt.stress(nlgw-1)) then
                          dstress(nlgw) = stress(nlgw-1)
                          stress(nlgw) = 0.0
                       else

                          !  If the stress value at the top of Layer  nlgw  is greater than
                          ! zero but less than the stress value at the bottom of layer  nlgw,
                          ! the stress value at the top of Layer  nlgw  is as calculated before
                          ! entering this  if  statement and the dstress value is set equal to
                          ! stress(nlgw) - stress(nlgw+1)

                          dstress(nlgw) = stress(nlgw-1) - stress(nlgw)
                       endif
                    endif
                    !     Endif NLGW.LE.NLAY-3

              endif
              !    Calculate the layer midoint wind accaleration value based upon
              !   the layer's value of  dstress  and the mass within the layer

              dudt(nlgw) = -dstress(nlgw) / ((pz(nlgw)-pz(nlgw-1))/gravv)

              !    Now, end the primary DO 22 loop

              !      write out to look at distribution of breaking vs layer
              if(dudt(nlgw).ne.0.0) then
                 !     write(17,*) msol,mbin,nlgw,dudt(nlgw),1.25+(NLAY-nlgw)*2.5
!                 write(17,*) msol,mbin,nlgw,dudt(nlgw),1.25+(nlgw-1)*2.5
              endif

              !    Populate the GWAccel array
              !   Be sure that the GWAccel array is properly zeroed upon each entry in to
              !  subroutine GW... the  "goto 22" statement would bypass the GWAccel
              !  assignment below, as would a   "goto 6"  statement

              GWAccel(ilon,ilat,nlgw,iEast_) = dudt(nlgw) * Dt

22            continue

           !  The above   22  continue statement ends the primary loop over
           ! vertical layers (1 to NLAY-1)


           !  Now, deal with the top model layer; if the stress value at the
           ! bottom of this top layer is zero, the layer's  dstress  value is
           ! set equal to the stress value at the bottom of that layer, and an
           ! acceleration is then calculated for the top layer;  no wave energy
           ! or momentum flux leaks out of the top of the model
           !
           stress(nlay) = 0.0
           dstress(nlay) = 0.0
           dudt(nlay) = 0.0

           if(abs(stress(nlay-1)).gt.1.0E-13) then
              dstress(nlay) = stress(nlay-1)
              dudt(nlay) = -dstress(nlay) / ((plev(nlay)-plev(nlay-1))/gravv)
!              write(17,*) msol,mbin,nlay,dudt(nlay),1.25+nlay*2.5
           endif

           !  Populate the top active layer of the GWAccel  array
           GWAccel(ilon,ilat,nlay,iEast_) = dudt(nlay) * Dt
           !  Now, write out results...
!           do nlgw=1,nlay
!              print *, u(nlgw),pz(nlgw),tgw(nlgw),theta(nlgw),dz(nlgw), &
!                   uav(nlgw),thtavg(nlgw),N(nlgw),dh(nlgw),Ri(nlgw),Rimin(nlgw),stress(nlgw), &
!                   dstress(nlgw),ts,dhsat(nlgw)
!              write(13,*) nlgw,stress(nlgw),dstress(nlgw),dudt,dz(nlgw),pz(nlgw),rhogw(nlgw)
!              write(14,99) nlgw,stress(nlgw),dstress(nlgw),dudt(nlgw),dz(nlgw),pz(nlgw), &
!                   Ri(nlgw),Rimin(nlgw),u(nlgw),t(nlgw),dh(nlgw),dhsat(nlgw),uav(nlgw), &
!                   eps(nlgw),1.250+(nlgw-1)*2.5
!              write(24,98) nlgw,stress(nlgw),dstress(nlgw),dudt(nlgw),dz(nlgw),pz(nlgw), &
!                   Ri(nlgw),Rimin(nlgw),u(nlgw),t(nlgw),dh(nlgw),dhsat(nlgw),uav(nlgw), &
!                   eps(nlgw),1.250+(nlgw-1)*2.5
!           enddo

99         format(i3,14(2x,1pe10.3))
98         format(i3,14(3x,1pe10.3))


        enddo
     enddo

     !  Set all GWAccel values to zero.. so there will be
     ! no effect upon the zonal wind speed, but all components
     ! of this subroutine will have been exercised
      GWAccel = 0.0
!     write(*,*) GWAccel(1,1,:,1)
   END subroutine CALC_GW

  !---------------------------------------------------------+
  subroutine  calc_lowatmosrad(iblock,iLat,iLon,L_LAYERS,L_LEVELS,&
       L_NLAYRAD,L_NLEVRAD)

    !  ---------------------------------------------------------------------
    !  Purpose: to calculate LTE CO2 15-micron heating/cooling in Mars
    !  lower atmosphere(below 80 km, where LTE physics dominates)
    !  Source: correlated k radiation code in Ames GCM v2.0
    !  Coders: Bob Haberle and Jim Schaeffer (2001)
    !          Converted to Mars GITM by S. M. Nelli (2007)
    !  Inputs: from MarsGITM  (Pressure, Temperature)
    !          all on MarsGITM grid (nLons,nLats,nAlts)
    !  Output: LowAtmosRadRate(1:nLons,1:nLats,1:nAlts,iBlock)
    !        **lowatmosradrate(iLon,iLat,L_NLAYRAD-L1,iBlock)=TOTAL(L1)
    !        : qnirtot(1:nLons,1:nLats,1:nAlts,iBlock)
    !        **qnirtot(iLon,iLat,L_NLAYRAD-L1,iBlock)= HEATING(L1)/XLTECORRECTION(L1)
    !        : qnirlte(1:nLons,1:nLats,1:nAlts,iBlock)
    !        **qnirlte(iLon,iLat,L_NLAYRAD-L1,iBlock)= HEATING(L1)
    !        : cirlte(1:nLons,1:nLats,1:nAlts,iBlock)
    !        **cirlte(iLon,iLat,L_NLAYRAD-L1,iBlock)= GREENHOUSE(L1)*COOLCORRECTION(L1)
    !        : all K/sec
    !        : Last Revised :  11-10-28 S. W. Bougher
    !  =======================================================================
    !
    use ModInputs
    use ModSources, only: LowAtmosRadRate,QnirLTE,QnirTOT,CirLTE
    use ModPlanet
    use ModGITM
    use ModEUV, only: AveCosSza
    use ModTime

    implicit none
    ! ----------------------------------------------------------------------

    integer,intent(in):: iBlock

    !C     Number of lower atmosphere layers
    integer,intent(in) :: L_LAYERS

    !C     Number of lower atmosphere levels:   2 * L_LAYERS + 3
    integer,intent(in) :: L_LEVELS

    !C     L_NLEVRAD corresponds to the surface - i.e., the GCM Level that
    !C     is at the surface.  PLEV(L_NLEVRAD) = P(J,I)+PTROP,
    !C     PLEV(2) = PTROP, PLEV(1) = ptop

    !C     L_NLAYRAD is the number of radiation code layers
    !C     L_NLEVRAD is the number of radiation code levels.  Level N is the
    !C               top of layer N.

    integer, intent(in) :: L_NLAYRAD, L_NLEVRAD, ilat,ilon

    LOGICAL :: isFirstDust(nBlocksMax) = .true.
    real(Real8_) :: DustTime(nDustTimes),ConrathTime(nConrathTimes),invDDiff,ctDiff(nDustTimes)
    REAL :: PLEV(LL_LEVELS+1), TLEV(LL_LEVELS)
    REAL :: PMID(LL_LEVELS), TMID(LL_LEVELS)
    REAL :: TAUREF(LL_LEVELS+1),TAUCUM(LL_LEVELS)
    REAL :: PTROP,ALS,ALBI
    REAL :: QH2O(LL_LEVELS)

    integer :: nw, n, nx, k, L,L1,MALT

    !C  VISUAL

    real :: DTAUV(LL_NLAYRAD,L_NSPECTV,L_NGAUSS)
    real :: TAUV(LL_NLEVRAD,L_NSPECTV,L_NGAUSS)
    real :: TAUCUMV(LL_LEVELS,L_NSPECTV,L_NGAUSS)
    real :: COSBV(LL_NLAYRAD,L_NSPECTV,L_NGAUSS)
    real :: WBARV(LL_NLAYRAD,L_NSPECTV,L_NGAUSS)
    real :: FMNETV(LL_NLAYRAD), diffvt, tDiff(nDustTimes)
    real :: conrnu(nDustAlts),tautot(nDustAlts)
    real :: fluxupv(LL_NLAYRAD), fluxdnv(LL_NLAYRAD), NFLUXTOPV

    real :: fluxvd(LL_LAYERS),HEATING(LL_LAYERS),TOTAL(LL_LAYERS),rtime,conrathrtime
    real :: fluxid(LL_LAYERS),GREENHOUSE(LL_LAYERS),XLTECORRECTION(LL_LAYERS)
    real :: timefactor(nDustAlts)

    integer :: ngwv(L_NSPECTV),ialt,imin(1),cmin(1)

    !C  IR

    real :: DTAUI(LL_NLAYRAD,L_NSPECTI,L_NGAUSS)
    real :: TAUCUMI(LL_LEVELS,L_NSPECTI,L_NGAUSS)
    real :: COSBI(LL_NLAYRAD,L_NSPECTI,L_NGAUSS)
    real :: WBARI(LL_NLAYRAD,L_NSPECTI,L_NGAUSS)
    real :: FMNETI(LL_NLAYRAD)
    real :: fluxupi(LL_NLAYRAD), fluxdni(LL_NLAYRAD), NFLUXTOPI
    real :: COOLCORRECTION(LL_LAYERS) !zero's out longwave cooling where
    !radcool is working
    integer ngwi(L_NSPECTI)

    ! Internal fields recast from MarsGITM

    real,dimension(nAlts) ::    &
         T, P

    !     Mars GITM real temperature (on its grid)  ----------------------

    T(1:nAlts) = &
         Temperature(iLon,iLat,1:nAlts,iBlock)*  &
         TempUnit(iLon,iLat,1:nAlts)

    !     Mars GITM pressure in mbars  (on its grid) -------------------

    P(1:nAlts) = &
         Pressure(iLon,iLat,1:nAlts,iBlock)*0.01

    !C              RADIATIVE CALCULATIONS.

    !C  Fill the new radiation code variables.
    !C  PLEV and TLEV are the pressure and temperatures on a vertical grid
    !C  that the new radiation code uses.

    CALL FILLPT(P,T,L_LEVELS,L_LAYERS,&
         SurfaceTemp(iLon,iLat,iBlock),altmin,&
         altitude_GB(iLon,iLat,1:nAlts,iBlock),&
         PLEV,TLEV,PMID,TMID,pressure(iLon,iLat,0,iBlock),&
         altitude_GB(iLon,iLat,0,iBlock))
    !

    ! Determining the top of the lower atmosphere radiative zone
    ptrop=plev(3)


    !C     Fill cumulative dust optical depth arrays (cum. dust optical
    !C     depth from the top of the lower atmosphere to the bottom of level K).

!write(*,*) plev(l_levels),l_levels
!stop

      !C  Fill QPI with water information

    !C  QH2O is the mixing ratio.  The GCM computes it as mass mixing ratio.
    !C  The radiation code wants number mixing ratio.  The TOTAL mass in each
    !C  layer is just the mass of CO2, not CO2+H2O.  If we change the total
    !C  mass, then the expresion for MWRATIO would change to read
    !C  MWRATIO = (MWCO2+MWH2O)/MWH2O.

    !   THIS HOOK IS FOR FUTURE POSSIBILITIES CONCERNING WATER VAPOR

    QH2O   = 1.0D-7

!!$          DO  L = 1, L_LAYERS
!!$            K = 2*L+2
!!$           QH2O(K)   = MWRATIO*QTRACE(JCMN,ICMN,L,M)
!!$           QH2O(K+1) = QH2O(K)
!!$          END DO

    !C  Set up, and solve for, the solar (visual) fluxes, if the sun
    !C  is up

    if(AveCosSza(iLon,iLat,iBlock).ge.1.0e-5) then

       !C  Check for ground ice.  Change albedo if there is any ice.

!      Scalar albedo for 1-D Tests with Brecht
!      ALS   = 0.24
!      2-D Array albedo for Production Simulations
       ALS   = SurfaceAlbedo(iLon,iLat,iBlock)

       ! FOR FUTURE CONSIDERATION OF CO2 GROUND ICE

!!$             IF(SurfaceTemp(iLon,iLat,iBlock).LE.150.0) THEN
!!$                IF(latitude(iLat,iBlock)*180.0/pi.LT.0.0) THEN
!!$                   ALS = ALICES
!!$                ELSE
!!$                   ALS = ALICEN
!!$                END IF
!!$             ENDIF

       !C  Get the optical depth (due to all sources) in the optical.
!write(*,*) ilon,ilat
       call optcv(DTAUV,TAUV,TAUCUMV,TLEV,PLEV,L_LAYERS,&
            L_LEVELS,L_NLAYRAD,L_NLEVRAD,WBARV,COSBV,&
            TAUREF,TMID,PMID,NGWV,QH2O)

       !C  Calculate the fluxes in the visual

       call sfluxv(DTAUV,TAUV,TAUCUMV,ALS,WBARV,COSBV,&
            AVECOSSZA(iLon,iLat,iBlock),NFLUXTOPV,FMNETV,&
            fluxupv,fluxdnv,diffvt,ngwv,&
            L_LAYERS,L_LEVELS,L_NLAYRAD,L_NLEVRAD,ilon,ilat)

       !C  If the sun is down, no solar flux, nor downward flux. . .

    else
       NFLUXTOPV = 0.0
       FLUXUPV = 0.0
       FLUXDNV = 0.0
       FMNETV = 0.0
    endif

    !C  Check for ground ice.  Change the IR albedo if there is any ice.

    ALBI   = 0.04

    ! FOR FUTURE CONSIDERATION OF CO2 GROUND ICE

!!$               IF(GT(JCMN,ICMN).LE.0.0) THEN
!!$                 IF(JCMN.LT.JEQUATOR) THEN
!!$                   ALBI = 1.0D0-EGOCO2S
!!$                 ELSE
!!$                   ALBI = 1.0D0-EGOCO2N
!!$                 END IF
!!$               ENDIF

    !C  Get the optical depth (due to all sources) in the infrared.

    call optci(DTAUI,TAUCUMI,TLEV,PLEV,L_LAYERS,&
         L_LEVELS,L_NLAYRAD,L_NLEVRAD,QextREF,COSBI,WBARI,&
         TAUREF,TMID,PMID,NGWI,QH2O)
!if (ilon .eq. 9 .and. ilat .eq. 9) stop
    !C  Calculate the fluxes in the IR.

    call sfluxi(PLEV,TLEV,DTAUI,TAUCUMI,ALBI,L_LAYERS,L_LEVELS,&
         L_NLAYRAD,L_NLEVRAD,COSBI,WBARI,NFLUXTOPI,FMNETI,&
         FLUXUPI,FLUXDNI,NGWI)
!write(*,*) fmneti
!stop

    do L=2,l_nlayrad
       L1=L-1
       IF(P(L_NLAYRAD-L1).GT.PLONG*0.01) THEN
          COOLCORRECTION(L1) = 1.0
       ELSE
          COOLCORRECTION(L1) = 0.0
       ENDIF
       IF(P(L_NLAYRAD-L1).GT.XLTEPRESSURE(1)) THEN
          XLTECORRECTION(L1) = 1.0
       ELSE
          DO N = 1,L_NLTE-1
             IF (P(L_NLAYRAD-L1).GT.XLTEPRESSURE(N)) THEN
                GOTO 199
             ELSE
                MALT = N
             ENDIF
          ENDDO

199       CONTINUE
          XLTECORRECTION(L1) = XLTEFACTOR(MALT) + &
               (XLTEFACTOR(MALT+1)-XLTEFACTOR(MALT))*&
               DLOG(P(L_NLAYRAD-L1)/XLTEPRESSURE(MALT))/&
               DLOG(XLTEPRESSURE(MALT+1)/XLTEPRESSURE(MALT))
       ENDIF

       fluxid(L1)  = FMNETI(L)-FMNETI(L-1)
       fluxvd(L1)  = FMNETV(L)-FMNETV(L-1)
       HEATING(L1) = fluxvd(L-1)*(-gravity_GB(iLon,iLat,L_NLAYRAD-L1,iBlock))/&
            (HeatCapacityCO2*100.0*(PLEV(2*L+1)-PLEV(2*L-1)))
!if (iproc .eq. 21 .and. ilon .eq. 1 .and. ilat .eq. 4) then
!   write(*,*)"LATM: ", fluxvd(l1),fmnetv(l)
!endif

       GREENHOUSE(L1) = fluxid(L-1)*(-gravity_GB(iLon,iLat,L_NLAYRAD-L1,iBlock))/&
            (HeatCapacityCO2*100.0*(PLEV(2*L+1)-PLEV(2*L-1)))

       TOTAL(L1)=(HEATING(L1)/XLTECORRECTION(L1))+&
            GREENHOUSE(L1)*COOLCORRECTION(L1)



       !Buffer region, Lopez-Valverde found a buffer needed for appropriate cooling rates
       !near the top of the atmosphere
       TOTAL(1) = 0.0
       TOTAL(2) = 0.0
       GREENHOUSE(1) = 0.0
       GREENHOUSE(2) = 0.0
       !Buffer region, Bougher found a buffer needed for appropriate near-IR rates
       !near the top of the atmosphere, since NLTE correction breaks down.
       HEATING(1) = 0.0
       HEATING(2) = 0.0

       lowatmosradrate(iLon,iLat,L_NLAYRAD-L1,iBlock)=TOTAL(L1)
!      qnirlte(iLon,iLat,L_NLAYRAD-L1,iBlock)= XLTECORRECTION(L1)
       qnirlte(iLon,iLat,L_NLAYRAD-L1,iBlock)= HEATING(L1)
       qnirtot(iLon,iLat,L_NLAYRAD-L1,iBlock)= HEATING(L1)/XLTECORRECTION(L1)
       cirlte(iLon,iLat,L_NLAYRAD-L1,iBlock)= GREENHOUSE(L1)*COOLCORRECTION(L1)


       dSubsurfaceTemp(iLon,iLat,iBlock) = -2.0*PI/sqrt(Pa*Pd)*&
            (SubsurfaceTemp(iLon,iLat,iBlock)-&
            SurfaceTemp(iLon,iLat,iBlock))+2.0*PI/Pa*&
            (CoreTemp-SubsurfaceTemp(iLon,iLat,iBlock))

       dSurfaceTemp(iLon,iLat,iBlock) = 1.0/(0.5*tinertia(iLon,iLat,iBlock)*&
            sqrt(Pd/PI))*(fluxdnv(L_NLAYRAD)*(1.0-SurfaceAlbedo(iLon,iLat,iBlock))+&
            fluxdni(L_NLAYRAD)-SBconstant*SurfaceTemp(iLon,iLat,iBlock)**4.0)+&
            2.0*PI/Pd*(SubsurfaceTemp(iLon,iLat,iBlock)-SurfaceTemp(iLon,iLat,iBlock))

    end do
!stop
    fir(iLon,iLat,iBlock) = fluxdni(L_NLAYRAD)
    fvis(iLon,iLat,iBlock) = fluxdnv(L_NLAYRAD)
    Tbot(iLon,iLat,iBlock) = T(1)
    TopL(iLon,iLat,iBlock) = L_NLAYRAD
    Psurf(iLon,iLat,iBlock)= PLEV(L_LEVELS)
    P125(iLon,iLat,iBlock) = P(50)

  end subroutine calc_lowatmosrad
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------

  subroutine fillpt(p,t,L_LEVELS,L_LAYERS,tg,altmin,altmid,plev,tlev,&
       pmid,tmid,pbot,altbot)

    !C  Put the T & P GCM arrays onto the NRC grid:  PLEV, PMID, TLEV, TMID
    !C  Sept 2002
    !C
    !C  PMID and TMID are the pressure and temperature at the GCM layer
    !C  mid-points.  PLEV and TLEV are the pressures and temperatures at
    !C  the GCM layer boundaries, i.e. at GCM levels.

!!!!!!Why is altmin and altbot different? !!!!!!!!!!!!!!!!!!!1

    use ModPlanet
    use ModGITM, only: iproc
    implicit none

    integer :: K, L, NK,L_LEVELS,L_LAYERS
    real :: PLEV(LL_LEVELS+1), PMID(LL_LEVELS)
    real :: TLEV(LL_LEVELS), TMID(LL_LEVELS)
    real :: TG, altmin,pbot,altbot
    real,dimension(1:nAlts) :: T, P,altmid
    real,dimension(0:ll_layers+1) :: altbound
    real,dimension(0:ll_levels) :: altlevels

    !C======================================================================C

    PBOT = PBOT*0.01 !convert bottom pressure to mbars
!altmin = 2429.59375
    ! Calculate boundary altitudes

!if (altmid(1) .gt. 0) then
!     ALTBOUND(0)=(ALTMIN+altmid(1))/2.
!  else
!     ALTBound(0) = Altmid(1) - (altmid(2) - altmid(1))/2.
!  endif
ALTBOUND(0) = altbot

!write(*,*) altbound(0)
!   write(*,*) altmid(1)
    do L = 1, L_Layers+1
       altbound(L) = (altmid(l) + altmid(l+1)) /2.
!       write(*,*) "alts: ",altbound(l), altmid(l),altmid(l+1)
    enddo
!    ALTBOUND(1)=ALTMID(1) + (ALTMID(1)-ALTMIN)
 !    DO L=2,L_LAYERS+1
!       ALTBOUND(L)=ALTMID(L) + (ALTMID(L)-ALTBOUND(L-1))
!       write(*,*) altbound(l),altmid(l), altbound(l-1)
!    END DO

    !C  Fill the new radiation code variables.
    !C  PLEV and TLEV are the pressure and tempertures on a vertical grid
    !C  that the new radiation code uses.

    ALTLEVELS(0) = ALTMID(L_LAYERS+2)
    ALTLEVELS(1) = ALTBOUND(L_LAYERS+1)
    ALTLEVELS(2) = ALTMID(L_LAYERS+1)
    ALTLEVELS(3) = ALTBOUND(L_LAYERS)
    ALTLEVELS(L_LEVELS+1) = (ALTBOUND(0)+Altmin)/2.

    PLEV(1) = EXP(DLOG(P(L_LAYERS+1)) + (DLOG(P(L_LAYERS+2))-&
         DLOG(P(L_LAYERS+1)))*&
         (ALTLEVELS(1)-ALTLEVELS(2))/&
         (ALTLEVELS(0)-ALTLEVELS(2)))
    PLEV(2) = P(L_LAYERS+1)
    PLEV(L_LEVELS) = pbot
    PLEV(L_LEVELS+1) = pbot
   DO K=4,L_LEVELS-1,2
       NK=L_LAYERS -(K/2 - 2)
       PLEV(K) = P(NK)
       ALTLEVELS(K) = ALTMID(NK)
       ALTLEVELS(K+1) = ALTBOUND(NK-1)

!       write(*,*) k,nk,altbound(nk-1)
!write(*,*) plev(k),k
!write(*,*) altmid(nk),altbound(nk),altlevels(k)-altlevels(k+1)
!if (iproc .eq. 21) write(*,*) k,nk,plev(k),P(nk),pbot
!write(*,*)
    END DO

!write(*,*) l_levels,l_layers
!if (iproc .eq. 21) stop
!!!!!!!!!!!!!!!!!! This is the issue (maybe) !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!stop
    DO K=3,L_LEVELS-2,2
       PLEV(K) = EXP(DLOG(PLEV(K+1)) + (DLOG(PLEV(K-1))-DLOG(PLEV(K+1)))*&
            (ALTLEVELS(K)-ALTLEVELS(K+1))/&
            (ALTLEVELS(K-1)-ALTLEVELS(K+1)))
!write(*,*) altlevels(k),k
!write(*,*) plev(k),k
    END DO

!    if (altbot .gt. altmin) then
       PLEV(L_LEVELS) = PBOT!EXP(DLOG(PBOT) - (DLOG(PLEV(L_LEVELS-1))-DLOG(PBOT))*&
!            (ALTBOT-ALTMIN)/(altbot-ALTLEVELS(L_LEVELS-1)))
!    else
!       PLEV(L_LEVELS) = EXP(DLOG(PBOT) + (DLOG(PLEV(L_LEVELS-1))-DLOG(PBOT))*&
!            (ALTBOT-ALTMIN)/(ALTLEVELS(L_LEVELS-1)-ALTBOT))
!    endif
    PLEV(L_LEVELS+1) = PLEV(L_LEVELS)


!    if (plev(l_levels) .lt. plev(l_levels -1)) then
!       write(*,*) "p: ",altbot,altmin,pbot,altbound(0),altmid(0),altmid(1)
!       do l = 1,l_levels+1
          !if (plev(l) .lt. plev(l-1)) then
!          write(*,*) l,plev(l),pbot,altlevels(l),altbot
!
!       enddo
!          stop
!    endif


!stop
    TLEV(1) = T(L_LAYERS+1) + (T(L_LAYERS+2)-T(L_LAYERS+1))*&
         DLOG(PLEV(1)/PLEV(2))/&
         DLOG(P(L_LAYERS+2)/PLEV(2))
    TLEV(2) = T(L_LAYERS+1)
    DO K=4,L_LEVELS-1,2
       NK=L_LAYERS -(K/2 - 2)
       TLEV(K) = T(NK)
    END DO

    DO K=3,L_LEVELS-2,2
       TLEV(K) = TLEV(K+1) + (TLEV(K-1)-TLEV(K+1))*&
            DLOG(PLEV(K)/PLEV(K+1))/&
            DLOG(PLEV(K-1)/PLEV(K+1))
    END DO

    !C  Temperature of the bottom level is the ground temperature.

    TLEV(L_LEVELS) = TG

    !C  Fill the PMID & TMID arrays used by OPTCI and OPTCV subroutines.
    !C  TMID and PMID used to get the index for CO2 k-coefficient interpolation.

    TMID(1) = TLEV(2)
    TMID(2) = TLEV(2)
    PMID(1) = PLEV(1)
    PMID(2) = PLEV(2)

    DO L=1,L_LAYERS
       TMID(2*L+1) = TLEV(2*L+1)
       TMID(2*L+2) = TLEV(2*L+1)
       PMID(2*L+1) = PLEV(2*L+1)
       PMID(2*L+2) = PLEV(2*L+1)
    END DO

    TMID(L_LEVELS) = TLEV(L_LEVELS)
    PMID(L_LEVELS) = PLEV(L_LEVELS)


  END SUBROUTINE FILLPT

  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  function jsrchgt(N,SX,INC,TARGET)

    !C  March 2002   c-grid
    !C  GCM2.0  Sept 2002
    !C
    !C      Find the first array element that is greater than the
    !C      TARGET value.  If N < 1, then 0 is returned.  If no
    !C      value is found, N+1 is returned.  Replaces Cray version
    !C      on the workstation.
    !C      Bisection search implemented 09/29/93 to improve speed.
    !C
    !C      Started: 08/23/93
    !C
    !C      Input:
    !C        N      - Number of elements in array SX.
    !C        SX     - Array of numbers to be searched.  Assumed to be an
    !C                 ordered array.
    !C        INC    - Increment between elements of searched array.  Kept
    !C                 for compatibility with Cray call.
    !C        TARGET - Value searched for in array SX.
    !C
    !C      Output:
    !C        JSRCHGT - location in array SX where value of SX is first
    !C                   greater than the TARGET value.
    !C----------------------------------------------------------------------C

    REAL :: SX(N), TARGET

    !C======================================================================C

    if(N.lt.1) then
       ians = 0
    elseif(TARGET.gt.SX(N)) then
       ians = N+1
    elseif(TARGET.lt.SX(1)) then
       ians = 1
    else

       JL = 1
       JH = N

10     CONTINUE
       if(JH-JL.gt.1) then
          JM = (JL+JH)/2

          if(TARGET.GT.SX(JM)) then
             JL = JM
             JM = (JL+JH)/2
          else
             JH = JM
             JM = (JL+JH)/2
          end if

          GOTO 10
       end if

       if(TARGET.EQ.SX(JH)) JH = JH+1

       ians = JH
    end if

    JSRCHGT = ians

  END FUNCTION


  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------

  SUBROUTINE OPTCV(DTAUV,TAUV,TAUCUMV,TLEV,PLEV,L_LAYERS,&
       L_LEVELS,L_NLAYRAD,L_NLEVRAD,WBARV,COSBV,&
       TAUREF,TMID,PMID,NGWV,QH2O)

!!$C  GCM2.0  Feb 2003
!!$C
!!$C THIS SUBROUTINE SETS THE OPTICAL CONSTANTS IN THE VISUAL
!!$C IT CALCUALTES FOR EACH LAYER, FOR EACH SPECRAL INTERVAL IN THE VISUAL
!!$C LAYER: WBAR, DTAU, COSBAR
!!$C LEVEL: TAU
!!$C
!!$C TAUV(L,NW,NG) is the cumulative optical depth at the top of radiation code
!!$C layer L. NW is spectral wavelength interval, ng the Gauss point index.
!!$C
!!$C     TLEV(L) - Temperature at the layer boundary
!!$C     PLEV(L) - Pressure at the layer boundary (i.e. level)
!!$C     CO2V(NT,NPS,NW,NG) - Visual CO2 k-coefficients
!!$C
!!$C----------------------------------------------------------------------C
    use ModPlanet
    use ModGITM, only: iproc
    implicit none

    integer :: L_LAYERS,L_LEVELS,L_NLAYRAD,L_NLEVRAD

    real :: DTAUV(LL_NLAYRAD,L_NSPECTV,L_NGAUSS)
    real :: DTAUKV(LL_LEVELS+1,L_NSPECTV,L_NGAUSS)
    real :: TAUV(LL_NLEVRAD,L_NSPECTV,L_NGAUSS)
    real :: TAUCUMV(LL_LEVELS,L_NSPECTV,L_NGAUSS)
    real :: TLEV(LL_LEVELS), PLEV(LL_LEVELS+1)
    real :: TMID(LL_LEVELS), PMID(LL_LEVELS)
    real :: COSBV(LL_NLAYRAD,L_NSPECTV,L_NGAUSS)
    real :: WBARV(LL_NLAYRAD,L_NSPECTV,L_NGAUSS)
    real :: TAUREF(LL_LEVELS+1)

    integer :: L, NW, NG, K, NG1(L_NSPECTV), LK
    integer :: MT(LL_LEVELS), MP(LL_LEVELS)
    real ::  ANS, TAUREFL, TAUGAS
    real ::  TRAY(LL_LEVELS,L_NSPECTV)
    real ::  TAEROS(LL_LEVELS,L_NSPECTV)
    real ::  DPR(LL_LEVELS), U(LL_LEVELS), TAUCLD
    real ::  LCOEF(4), LKCOEF(LL_LEVELS,4)

    real :: taugsurf(L_NSPECTV,L_NGAUSS-1), TRAYAER

    !C  Tlimit:  If the CO2 optical depth (top to the surface) is less than
    !C  this value, we place that Gauss-point into the "zeros" channel.
    !C  Set in driver, passed via common


    integer :: ngwv(L_NSPECTV)

    !C  Reference wavelength is (now) bin #6 - put into qextref

    !C  Water mixing ratio stuff

    real :: QH2O(LL_LEVELS), WRATIO(LL_LEVELS)
    real :: KCOEF(4)
    integer :: nh2o(LL_LEVELS)

    !C======================================================================C

    QextREF = Qextv(L_NREFV)

    !C  Determine the total gas opacity throughout the column, for each
    !C  spectral interval, NW, and each Gauss point, NG.
    !C  Calculate the continuum opacities, i.e., those that do not depend on
    !C  NG, the Gauss index.

    DO NG=1,L_NGAUSS-1
       do NW=1,L_NSPECTV
          TAUGSURF(NW,NG) = 0.0D0
       end do
    end do

    do K=2,L_LEVELS
       DPR(k) = PLEV(K)-PLEV(K-1)
       U(k)   = Cmk*DPR(k)
!       if (u(k) .lt. 0) then
!          write(*,*) k,plev(k),plev(k-1)
!          stop
!       endif
       call tpindex(PMID(K),TMID(K),QH2O(K),pfgasref,tgasref,&
            LCOEF,MT(K),MP(K),NH2O(K),WRATIO(K))

       do LK=1,4
          LKCOEF(K,LK) = LCOEF(LK)
       end do

       DO NW=1,L_NSPECTV
          TRAY(K,NW)   = TAURAY(NW)*DPR(K)
          TAEROS(K,NW) = TAUREF(K)*Qextv(NW)/QextREF
       END DO
    end do

    !C  TAUCLD = is cloud opacity, zero until further notice
    !C  TRAYAER is Tau RAYleigh scattering, plus AERosol opacity

    TAUCLD = 0.0D0
    do NW=1,L_NSPECTV
       ngwv(nw) = L_NGAUSS

       !C  Now fill in the "clear" part of the spectrum (NG = L_NGAUSS)
       !C  Which holds continuum opacity only

       do K=2,L_LEVELS
          DTAUKV(K,nw,L_NGAUSS) = TAEROS(K,NW)+TRAY(K,NW) + TAUCLD
       end do

       do ng=L_NGAUSS-1,1,-1
          do K=2,L_LEVELS

             !C           NOW COMPUTE TAUGAS

             !C  Interpolate between water mixing ratios
             !C  WRATIO = 0.0 if the requested water amount is equal to, or outside the
             !C  the range of water amount data.

             KCOEF(1) = CO2V(MT(K),MP(K),NH2O(K),NW,NG) + WRATIO(K)*&
                  (CO2V(MT(K),MP(K),NH2O(K)+1,NW,NG) - &
                  CO2V(MT(K),MP(K),NH2O(K),NW,NG))

             KCOEF(2) = CO2V(MT(K),MP(K)+1,NH2O(K),NW,NG) + WRATIO(K)*&
                  (CO2V(MT(K),MP(K)+1,NH2O(K)+1,NW,NG) - &
                  CO2V(MT(K),MP(K)+1,NH2O(K),NW,NG))

             KCOEF(3) = CO2V(MT(K)+1,MP(K)+1,NH2O(K),NW,NG) + WRATIO(K)*&
                  (CO2V(MT(K)+1,MP(K)+1,NH2O(K)+1,NW,NG) - &
                  CO2V(MT(K)+1,MP(K)+1,NH2O(K),NW,NG))

             KCOEF(4) = CO2V(MT(K)+1,MP(K),NH2O(K),NW,NG) + WRATIO(K)*&
                  (CO2V(MT(K)+1,MP(K),NH2O(K)+1,NW,NG) - &
                  CO2V(MT(K)+1,MP(K),NH2O(K),NW,NG))

             !C  Interpolate the CO2 k-coefficients to the requested T,P


             ANS = LKCOEF(K,1)*KCOEF(1) + LKCOEF(K,2)*KCOEF(2) +&
                  LKCOEF(K,3)*KCOEF(3) + LKCOEF(K,4)*KCOEF(4)

             TAUGAS          = U(k)*ANS
             TAUGSURF(NW,NG) = TAUGSURF(NW,NG) + TAUGAS
             DTAUKV(K,nw,ng) = TAUGAS + TAUCLD + TRAY(K,NW) +&
                  TAEROS(K,NW)
!             if (u(k) .lt. 0) write(*,*) u(k),k,nw,ng
!             if (dtaukv(k,nw,ng) .lt. -700000) then
!                write(*,*) dtaukv(k,nw,ng),taugas,u(k),ans
!                stop
!endif

          end do
          if(TAUGSURF(NW,NG) .LT. TLIMIT) THEN
             goto 10
          else
             NGWV(NW) = NG
          end if

       end do
10     continue
    end do

    !C  Now the full treatment for the layers, where besides the opacity
    !C  we need to calculate the scattering albedo and asymmetry factors
    !C  for each layer

    DO NW=1,L_NSPECTV

       !C  First, the special "clear" channel

       NG = L_NGAUSS
       DO L=1,L_LAYERS
          K              = 2*L+1
          TAUREFL        = (TAUREF(K)+TAUREF(K+1))/QextREF
          DTAUV(L,nw,ng) = DTAUKV(K,NW,NG)+DTAUKV(K+1,NW,NG)
          COSBV(L,NW,NG) = (GV(NW)*Qscatv(NW)*TAUREFL)/&
               (TRAY(K,NW)+TRAY(K+1,NW) + &
               QSCATV(NW)*TAUREFL)
          WBARV(L,nw,ng) = (QSCATV(NW)*TAUREFL + &
               (TRAY(K,NW)+TRAY(K+1,NW))*0.9999)/&
               DTAUV(L,nw,ng)
       END DO

       !C  Special bottom layer

       L              = L_NLAYRAD
       K              = 2*L+1
       TAUREFL        = TAUREF(K)/QextREF
       DTAUV(L,nw,ng) = DTAUKV(K,NW,NG)
       COSBV(L,NW,NG) = (GV(NW)*Qscatv(NW)*TAUREFL)/&
            (TRAY(K,NW) + QSCATV(NW)*TAUREFL)
       WBARV(L,nw,ng) = (QSCATV(NW)*TAUREFL + TRAY(K,NW)*0.9999)/&
            DTAUV(L,nw,ng)

    END DO

    !C  . . .Now the other Gauss points, if needed.

    do NW=1,L_NSPECTV
       DO NG=L_NGAUSS-1,NGWV(NW),-1
          DO L=1,L_LAYERS
             K              = 2*L+1
             TAUREFL        = (TAUREF(K)+TAUREF(K+1))/QextREF
             DTAUV(L,nw,ng) = DTAUKV(K,NW,NG)+DTAUKV(K+1,NW,NG)
             COSBV(L,NW,NG) = COSBV(L,NW,L_NGAUSS)
             WBARV(L,nw,ng) = (QSCATV(NW)*TAUREFL + &
                  (TRAY(K,NW)+TRAY(K+1,NW))*0.9999)/&
                  DTAUV(L,nw,ng)
          END DO

          !C  Special bottom layer

          L              = L_NLAYRAD
          K              = 2*L+1
          TAUREFL        = TAUREF(K)/QextREF
          DTAUV(L,nw,ng) = DTAUKV(K,NW,NG)
          COSBV(L,NW,NG) = (GV(NW)*Qscatv(NW)*TAUREFL)/&
               (TRAY(K,NW) + QSCATV(NW)*TAUREFL)
          WBARV(L,nw,ng) = (QSCATV(NW)*TAUREFL + TRAY(K,NW)*0.9999)/&
               DTAUV(L,nw,ng)
!          if (dtauv(l,nw,ng) .lt. -700000) then
!             write(*,*) dtauv(l,nw,ng),dtaukv(k,nw,ng),l,nw,ng
!             stop
!             endif
       END DO

    END DO     ! NW spectral loop

    !C     TOTAL EXTINCTION OPTICAL DEPTHS

    DO NW=1,L_NSPECTV
       NG = L_NGAUSS
       TAUV(1,NW,NG) = 0.0D0
       DO L=1,L_NLAYRAD
          TAUV(L+1,NW,NG) = TAUV(L,NW,NG)+DTAUV(L,NW,NG)
       END DO

       TAUCUMV(1,NW,NG)=0.0D0
       DO K=2,L_LEVELS
          TAUCUMV(K,NW,NG)=TAUCUMV(K-1,NW,NG)+DTAUKV(K,NW,NG)
       END DO

       DO NG=L_NGAUSS-1,NGWV(NW),-1
          TAUV(1,NW,NG)=0.0D0
          DO L=1,L_NLAYRAD
             TAUV(L+1,NW,NG)=TAUV(L,NW,NG)+DTAUV(L,NW,NG)
          END DO

          TAUCUMV(1,NW,NG)=0.0D0
          DO K=2,L_LEVELS
             TAUCUMV(K,NW,NG)=TAUCUMV(K-1,NW,NG)+DTAUKV(K,NW,NG)
          END DO
       END DO
    END DO


  END SUBROUTINE OPTCV
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------

  subroutine tpindex(pw,tw,qh2o,pref,tref,LCOEF,MT,MP,&
       NH2O,wratio)

!!$C  GCM2.0  Feb 2003
!!$C
!!$C
!!$C  PURPOSE
!!$C    Get the TI, UI values for a 2-dimensional interpolation
!!$C    based on the following (The interpolation is done in interpco2):
!!$C    Interpolate the CO2 K-coefficients to the current P,T values.
!!$C    The CO2 coefficients are given on a P,T grid:
!!$C    P = {1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1, 1E+1, 1E+2, 1E+3, 1E+4},
!!$C    T = {50, 100, 150, 200, 250, 300, 350}.
!!$C
!!$C    The interpolation is the usual interpolation in 2-dimensions given
!!$C    in "Numerical Recipes", where the "X" are P, the "Y" are
!!$C    T, and the F(X,Y) are the CO2 K-coefficients.
!!$C
!!$C     The interpolating box is designated as follows:
!!$C
!!$C           (PL,TU)                        (PRR,TU)
!!$C
!!$C                          (TW,PW)
!!$C
!!$C
!!$C           (PL,TL)                        (PRR,TL)
!!$C
!!$C     PL  - Pressure left
!!$C     PRR - Pressure right
!!$C     TL  - Temperature lower
!!$C     TU  - Temperature upper
!!$C     PW  - Pressure wanted
!!$C     TW  - Temperature wanted
!!$C
!!$C
!!$C  AUTHOR
!!$C      JIM SCHAEFFER     STERLING SOFTWARE  TASK 404  MAR 1998
!!$C
!!$C  FOR
!!$C      BOB HABERLE       PART OF GCM II/RADTRAN UPDATES
!!$C
!!$C  ENVIRONMENT
!!$C      SUN ULTRA-2       SOLARIS 2.5.1   FORTRAN 77
!!$C
!!$C  REVISION HISTORY
!!$C     3-98  JRS  ORIGINAL
!!$C
!!$C  INPUT PARAMETERS
!!$C    PW                 - The pressure to interpolate to
!!$C    TW                 - The temperature to interpolate to
!!$C    Pref(NP)           - The pressure grid array.
!!$C    Tref(NT)           - The temperature grid array.
!!$C
!!$C  OUTPUT PARAMETERS
!!$C    TI                 - Interpolation term (pressure)
!!$C    UI                 - Interpolation term (temperature)
!!$C    MT                 - Temperature index (bottom left Temperature)
!!$C                         of bounding box
!!$C    MP                 - Pressure index (bottom left pressure)
!!$C                         of bounding box
!!$C
!!$C  CALLED BY
!!$C    SETRAD
!!$C
!!$C  SUBROUTINES CALLED
!!$C      NONE
!!$C
!!$C----------------------------------------------------------------------C

    use ModPlanet

    implicit none


    real :: Tref(L_NTREF)
    real :: pref(L_PINT)

    integer :: MT, MP, N, M, NH2O
    real ::  PW, TW, Qh2o, wratio
    real ::  PWL, LCOEF(4), T, U

    !C======================================================================C

    !C     Get the upper and lower Temperature-grid indicies that bound the
    !C     requested temperature.  If the requested temperature is outside
    !C     the T-grid, set up to extrapolate from the appropriate end.

    IF(TW.LE.TREF(1)) THEN
       MT = 1
       U  = 0.0D0
    ELSE
       do n=1,L_NTREF-1
          if(tw.gt.Tref(n) .and. TW.LE.TREF(N+1)) then
             MT = n
             U = (TW-TREF(MT))/(TREF(MT+1)-TREF(MT))
             goto 10
          end if
       end do

       MT = L_NTREF-1
       U  = 1.0D0

10     continue
    END IF


    !C     Get the upper and lower Pressure-grid indicies that bound the
    !C     requested pressure.  If the requested pressure is outside
    !C     the P-grid, set up to extrapolate from the appropiate end.

    pwl = log10(pw)

    if(pwl.le.Pref(1)) then
       MP = 1
       T  = 0.0D0
    else
       do n=2,L_PINT-1
          if(pwl.le.Pref(n)) then
             MP = n-1
             T = (PWL-PREF(MP))/(PREF(MP+1)-PREF(MP))
             goto 20
          end if
       end do

       MP = L_PINT-1
       T  = 1.0D0

20     continue
    end if

    !C  Fill the interpolation coeficients:

    LCOEF(1) = (1.0-T)*(1.0-U)
    LCOEF(2) = T*(1.0-U)
    LCOEF(3) = T*U
    LCOEF(4) = (1.0-T)*U

    !C  Get the indicies for water abundance.  There are 10 sets of
    !C  k-coefficients with differing amounts of water vs. CO2.

    IF(Qh2o.le.WREFH2O(1)) then
       NH2O   = 1
       WRATIO = 0.0D0
    ELSEIF(Qh2o.ge.WREFH2O(L_REFH2O)) then
       NH2O   = L_REFH2O - 1
       WRATIO = 1.0D0
    ELSE
       DO N=2,L_REFH2O
          IF(QH2O.GE.WREFH2O(N-1) .and. QH2O.lt.WREFH2O(N)) then
             NH2O   = N-1
             WRATIO = (QH2O - WREFH2O(N-1))/(WREFH2O(N) - WREFH2O(N-1))
             GOTO 30
          END IF
       END DO
    END IF

30  CONTINUE


  end subroutine tpindex

  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------

  SUBROUTINE SFLUXV(DTAUV,TAUV,TAUCUMV,RSFV,WBARV,COSBV,&
       UBAR0,NFLUXTOPV,FMNETV,&
       FLUXUPV,FLUXDNV,DIFFVT,NGWV,&
       L_LAYERS,L_LEVELS,L_NLAYRAD,L_NLEVRAD,ilon,ilat)

    !C  GCM2.0  Feb 2003

    use ModPlanet
    use ModGITM, only: iproc
    implicit none

    integer :: L_LAYERS,L_LEVELS,L_NLAYRAD,L_NLEVRAD,ilon,ilat

    real :: FMNETV(LL_NLAYRAD)
    real :: TAUCUMV(LL_LEVELS,L_NSPECTV,L_NGAUSS)
    real :: TAUV(LL_NLEVRAD,L_NSPECTV,L_NGAUSS)
    real :: DTAUV(LL_NLAYRAD,L_NSPECTV,L_NGAUSS)
    real :: FMUPV(LL_NLAYRAD), FMDV(LL_NLAYRAD)
    real :: COSBV(LL_NLAYRAD,L_NSPECTV,L_NGAUSS)
    real :: WBARV(LL_NLAYRAD,L_NSPECTV,L_NGAUSS)
    real :: FLUXUPV(LL_NLAYRAD), FLUXDNV(LL_NLAYRAD)
    real :: NFLUXTOPV, FLUXUP, FLUXDN

    integer :: L, NG, NW, NG1
    real ::  rsfv, ubar0, f0pi, btop, bsurf, taumax, eterm

    real :: DIFFV, DIFFVT

    integer :: ngwv(L_NSPECTV)

    real :: fzero

    !C======================================================================C

    TAUMAX = L_TAUMAX

    !C     ZERO THE NET FLUXES

    NFLUXTOPV = 0.0

    DO L=1,L_NLAYRAD
       FMNETV(L)  = 0.0
       FLUXUPV(L) = 0.0
       FLUXDNV(L) = 0.0
    END DO

    DIFFVT = 0.0

    !C     WE NOW ENTER A MAJOR LOOP OVER SPECTRAL INTERVALS IN THE VISIBLE
    !C     TO CALCULATE THE NET FLUX IN EACH SPECTRAL INTERVAL

    DO 500 NW=1,L_NSPECTV

       F0PI = SOL(NW)

       FZERO = FZEROV(NW)
       IF(FZERO.ge.0.99) goto 40

       !C  Skip over the gauss points with minimal (opacity < TLIMIT) total
       !C  gas opacity

       DO NG=1,NGWV(NW)-1
          fzero = fzero + (1.0-FZEROV(NW))*GWEIGHT(NG)
       END DO

       !C  This loop includes only those Gauss points with sufficient gas opacity

       DO NG=NGWV(NW),L_NGAUSS-1

          !C         SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE VISIBLE

          BTOP = 0.0

          !C         LOOP OVER THE NTERMS BEGINNING HERE

          ETERM = MIN(TAUV(L_NLEVRAD,NW,NG),TAUMAX)
          BSURF = RSFV*UBAR0*SOL(NW)*EXP(-ETERM/UBAR0)

          !C         WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
          !C         CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
          !C         WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER
          !C
          !C         FUW AND FDW ARE WORKING FLUX ARRAYS THAT WILL BE USED TO
          !C         RETURN FLUXES FOR A GIVEN NT

          CALL GFLUXV(DTAUV(1,NW,NG),TAUV(1,NW,NG),TAUCUMV(1,NW,NG),&
               WBARV(1,NW,NG),COSBV(1,NW,NG),UBAR0,F0PI,RSFV,&
               BTOP,BSURF,FMUPV,FMDV,DIFFV,FLUXUP,FLUXDN,&
               L_LAYERS,L_LEVELS,L_NLAYRAD,L_NLEVRAD)

          !C         NOW CALCULATE THE CUMULATIVE VISIBLE NET FLUX

          NFLUXTOPV = NFLUXTOPV+(FLUXUP-FLUXDN)*GWEIGHT(NG)*&
               (1.0-FZEROV(NW))
          DO L=1,L_NLAYRAD
             FMNETV(L)=FMNETV(L)+( FMUPV(L)-FMDV(L) )*&
                  GWEIGHT(NG)*(1.0-FZEROV(NW))
             FLUXUPV(L) = FLUXUPV(L) + FMUPV(L)*GWEIGHT(NG)*&
                  (1.0-FZEROV(NW))
             FLUXDNV(L) = FLUXDNV(L) + FMDV(L)*GWEIGHT(NG)*&
                  (1.0-FZEROV(NW))

          END DO

          !C         THE DIFFUSE COMPONENT OF THE DOWNWARD SOLAR FLUX

          DIFFVT = DIFFVT + DIFFV*GWEIGHT(NG)*(1.0-FZEROV(NW))

       END DO   ! the Gauss loop

40     continue

       !C       Special 17th Gauss point

       NG = L_NGAUSS

       !C       SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE VISIBLE

       BTOP = 0.0

       !C       LOOP OVER THE NTERMS BEGINNING HERE

       ETERM = MIN(TAUV(L_NLEVRAD,NW,NG),TAUMAX)
       BSURF = RSFV*UBAR0*SOL(NW)*EXP(-ETERM/UBAR0)

       !C       WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
       !C       CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
       !C       WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER
       !C
       !C       FUW AND FDW ARE WORKING FLUX ARRAYS THAT WILL BE USED TO
       !C       RETURN FLUXES FOR A GIVEN NT

       CALL GFLUXV(DTAUV(1,NW,NG),TAUV(1,NW,NG),TAUCUMV(1,NW,NG),&
            WBARV(1,NW,NG),COSBV(1,NW,NG),UBAR0,F0PI,RSFV,&
            BTOP,BSURF,FMUPV,FMDV,DIFFV,FLUXUP,FLUXDN,&
            L_LAYERS,L_LEVELS,L_NLAYRAD,L_NLEVRAD)

       !C       NOW CALCULATE THE CUMULATIVE VISIBLE NET FLUX

       NFLUXTOPV = NFLUXTOPV+(FLUXUP-FLUXDN)*FZERO
       DO L=1,L_NLAYRAD
          FMNETV(L)=FMNETV(L)+( FMUPV(L)-FMDV(L) )*FZERO
          FLUXUPV(L) = FLUXUPV(L) + FMUPV(L)*FZERO
          FLUXDNV(L) = FLUXDNV(L) + FMDV(L)*FZERO


       END DO
       !C       THE DIFFUSE COMPONENT OF THE DOWNWARD SOLAR FLUX

       DIFFVT = DIFFVT + DIFFV*FZERO

500    CONTINUE

       !C     *** END OF MAJOR SPECTRAL INTERVAL LOOP IN THE VISIBLE*****


     END SUBROUTINE SFLUXV

     ! ----------------------------------------------------------------------

     ! ----------------------------------------------------------------------

     SUBROUTINE GFLUXV(DTDEL,TDEL,TAUCUM,WDEL,CDEL,UBAR0,F0PI,RSF,&
          BTOP,BSURF,FMIDP,FMIDM,DIFFV,FLUXUP,FLUXDN,&
          L_LAYERS,L_LEVELS,L_NLAYRAD,L_NLEVRAD)

!!$C  GCM2.0  Feb 2003
!!$C
!!$C ?? COMBINE INTO ONE GLFLUX ROUTINE, NOT VIS AND IR
!!$C
!!$C  THIS SUBROUTINE TAKES THE OPTICAL CONSTANTS AND BOUNDARY CONDITIONS
!!$C  FOR THE VISIBLE  FLUX AT ONE WAVELENGTH AND SOLVES FOR THE FLUXES AT
!!$C  THE LEVELS. THIS VERSION IS SET UP TO WORK WITH LAYER OPTICAL DEPTHS
!!$C  MEASURED FROM THE TOP OF EACH LAYER.  (DTAU) TOP OF EACH LAYER HAS
!!$C  OPTICAL DEPTH TAU(N).IN THIS SUB LEVEL N IS ABOVE LAYER N. THAT IS LAYER N
!!$C  HAS LEVEL N ON TOP AND LEVEL N+1 ON BOTTOM. OPTICAL DEPTH INCREASES
!!$C  FROM TOP TO BOTTOM. SEE C.P. MCKAY, TGM NOTES.
!!$C THIS SUBROUTINE DIFFERS FROM ITS IR COUNTERPART IN THAT HERE WE SOLVE FOR
!!$C THE FLUXES DIRECTLY USING THE GENERALIZED NOTATION OF MEADOR AND WEAVOR
!!$C J.A.S., 37, 630-642, 1980.
!!$C THE TRI-DIAGONAL MATRIX SOLVER IS DSOLVER AND IS DOUBLE PRECISION SO MANY
!!$C VARIABLES ARE PASSED AS SINGLE THEN BECOME DOUBLE IN DSOLVER
!!$C
!!$C NLL           = NUMBER OF LEVELS (NAYER + 1) THAT WILL BE SOLVED
!!$C NAYER         = NUMBER OF LAYERS (NOTE DIFFERENT SPELLING HERE)
!!$C WAVEN         = WAVELENGTH FOR THE COMPUTATION
!!$C DTDEL(NLAYER) = ARRAY OPTICAL DEPTH OF THE LAYERS
!!$C TDEL(NLL)     = ARRAY COLUMN OPTICAL DEPTH AT THE LEVELS
!!$C WDEL(NLEVEL)  = SINGLE SCATTERING ALBEDO
!!$C CDEL(NLL)     = ASYMMETRY FACTORS, 0=ISOTROPIC
!!$C UBARV         = AVERAGE ANGLE,
!!$C UBAR0         = SOLAR ZENITH ANGLE
!!$C F0PI          = INCIDENT SOLAR DIRECT BEAM FLUX
!!$C RSF           = SURFACE REFLECTANCE
!!$C BTOP          = UPPER BOUNDARY CONDITION ON DIFFUSE FLUX
!!$C BSURF         = REFLECTED DIRECT BEAM = (1-RSFI)*F0PI*EDP-TAU/U
!!$C FP(NLEVEL)    = UPWARD FLUX AT LEVELS
!!$C FM(NLEVEL)    = DOWNWARD FLUX AT LEVELS
!!$C FMIDP(NLAYER) = UPWARD FLUX AT LAYER MIDPOINTS
!!$C FMIDM(NLAYER) = DOWNWARD FLUX AT LAYER MIDPOINTS
!!$C added Dec 2002
!!$C DIFFV         = downward diffuse solar flux at the surface
!!$C

       use ModPlanet
       use ModGITM, only: iProc

       implicit none
       integer :: L_LAYERS,L_LEVELS,L_NLAYRAD,L_NLEVRAD
!      integer,PARAMETER :: NL=101  ! MUST BE LARGER THAN NLL
       integer,PARAMETER :: NL3=101  ! MUST BE LARGER THAN NLL

       REAL :: EM, EP
       REAL :: W0(LL_NLAYRAD), COSBAR(LL_NLAYRAD), DTAU(LL_NLAYRAD)
       REAL :: TAU(LL_NLEVRAD), WDEL(LL_NLAYRAD), CDEL(LL_NLAYRAD)
       REAL :: DTDEL(LL_NLAYRAD), TDEL(LL_NLEVRAD)
       REAL :: FMIDP(LL_NLAYRAD), FMIDM(LL_NLAYRAD)
       REAL :: LAMDA(NL3), ALPHA(NL3),XK1(NL3),XK2(NL3)
       REAL :: G1(NL3), G2(NL3), G3(NL3), GAMA(NL3), CP(NL3), CM(NL3), CPM1(NL3)
       REAL :: CMM1(NL3), E1(NL3), E2(NL3), E3(NL3), E4(NL3), EXPTRM(NL3)
       REAL :: FLUXUP, FLUXDN
       REAL :: TAUCUM(LL_LEVELS), taucump(LL_LEVELS+2)

       integer :: NAYER, L, K
       real ::  ubar0, f0pi, rsf, btop, bsurf, g4, denom, am, ap
       real ::  taumax, taumid, cpmid, cmmid
       real ::  diffv

       !C======================================================================C

       NAYER  = L_NLAYRAD
       TAUMAX = L_TAUMAX    !Default is 35.0

!!$C  Delta function updates.  Use scaled values, calculated below.
!!$c     DO L=1,L_NLAYRAD
!!$c       W0(L)     = WDEL(L)
!!$c       COSBAR(L) = CDEL(L)
!!$c       DTAU(L)   = DTDEL(L)
!!$c       TAU(L)    = TDEL(L)
!!$c     END DO
!!$c     TAU(L_NLEVRAD) = TDEL(L_NLEVRAD)
!!$
!!$C  Scaled (Delta function) values

       TAU(1)=TDEL(1)*(1.-WDEL(1)*CDEL(1)**2)
       TAUCUMP(1) = 0.0
       TAUCUMP(2) = TAUCUM(2)*(1.-WDEL(1)*CDEL(1)**2)
       TAUCUMP(3) = TAUCUM(2)+ (TAUCUM(3)-TAUCUM(2))*&
            (1.-WDEL(1)*CDEL(1)**2)
       DO L=1,L_NLAYRAD-1
          W0(L)        = WDEL(L)*(1.-CDEL(L)**2)/&
               (1.-WDEL(L)*CDEL(L)**2)
          COSBAR(L)    = CDEL(L)/(1.+CDEL(L))
          DTAU(L)      = DTDEL(L)*(1.-WDEL(L)*CDEL(L)**2)
          TAU(L+1)     = TAU(L)+DTAU(L)
          K            = 2*(L+1)
          TAUCUMP(K)   = TAU(L+1)
          TAUCUMP(K+1) = TAUCUMP(k)+(TAUCUM(K+1) - TAUCUM(k))*&
               (1.-WDEL(L)*CDEL(L)**2)


       END DO

       !C  Bottom layer

       L=L_NLAYRAD
       W0(L)      = WDEL(L)*(1.-CDEL(L)**2)/(1.-WDEL(L)*CDEL(L)**2)
       COSBAR(L)  = CDEL(L)/(1.+CDEL(L))
       DTAU(L)    = DTDEL(L)*(1.-WDEL(L)*CDEL(L)**2)
       TAU(L+1)   = TAU(L)+DTAU(L)
       K          = 2*(L+1)
       TAUCUMP(K) = TAU(L+1)
!       if  (iproc .eq. 21 .and. dtau(47) .lt. -700000) then
!          write(*,*) "dtau: ",DTAU(47),DTDEL(47),WDEL(47),CDEL(47),l_nlayrad
!       endif
       !C     WE GO WITH THE QUADRATURE APPROACH HERE.  THE "SQRT(3)" factors
       !C     ARE THE UBARV TERM.

       DO L=1,L_NLAYRAD

          ALPHA(L)=SQRT( (1.0-W0(L))/(1.0-W0(L)*COSBAR(L)) )

          !C       SET OF CONSTANTS DETERMINED BY DOM

          G1(L)    = (SQRT(3.0)*0.5)*(2.0- W0(L)*(1.0+COSBAR(L)))
          G2(L)    = (SQRT(3.0)*W0(L)*0.5)*(1.0-COSBAR(L))
          G3(L)    = 0.5*(1.0-SQRT(3.0)*COSBAR(L)*UBAR0)
          LAMDA(L) = SQRT(G1(L)**2 - G2(L)**2)
          GAMA(L)  = (G1(L)-LAMDA(L))/G2(L)
       END DO

       DO L=1,L_NLAYRAD
          G4    = 1.0-G3(L)
          DENOM = LAMDA(L)**2 - 1./UBAR0**2

          !C       THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
          !C       THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN
          !C       THE SCATTERING GOES TO ZERO
          !C       PREVENT THIS WITH AN IF STATEMENT

          IF ( DENOM .EQ. 0.) THEN
             DENOM=1.E-10
          END IF

          AM = F0PI*W0(L)*(G4   *(G1(L)+1./UBAR0) +G2(L)*G3(L) )/DENOM
          AP = F0PI*W0(L)*(G3(L)*(G1(L)-1./UBAR0) +G2(L)*G4    )/DENOM

          !C       CPM1 AND CMM1 ARE THE CPLUS AND CMINUS TERMS EVALUATED
          !C       AT THE TOP OF THE LAYER, THAT IS LOWER   OPTICAL DEPTH TAU(L)

          CPM1(L) = AP*EXP(-TAU(L)/UBAR0)
          CMM1(L) = AM*EXP(-TAU(L)/UBAR0)

          !C       CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
          !C       BOTTOM OF THE LAYER.  THAT IS AT HIGHER OPTICAL DEPTH TAU(L+1)

          CP(L) = AP*EXP(-TAU(L+1)/UBAR0)
          CM(L) = AM*EXP(-TAU(L+1)/UBAR0)

       END DO

       !C     NOW CALCULATE THE EXPONENTIAL TERMS NEEDED
       !C     FOR THE TRIDIAGONAL ROTATED LAYERED METHOD

       DO L=1,L_NLAYRAD
          EXPTRM(L) = MIN(TAUMAX,LAMDA(L)*DTAU(L))  ! CLIPPED EXPONENTIAL
          EP = EXP(EXPTRM(L))

          EM        = 1.0/EP
          E1(L)     = EP+GAMA(L)*EM
          E2(L)     = EP-GAMA(L)*EM
          E3(L)     = GAMA(L)*EP+EM
          E4(L)     = GAMA(L)*EP-EM
!          if (E2(L) .ne. e2(L)) then
!             write(*,*) "e2: ",exptrm(L),taumax,lamda(l),dtau(l),iproc,l
!          endif
       END DO

       CALL DSOLVER(NAYER,GAMA,CP,CM,CPM1,CMM1,E1,E2,E3,E4,BTOP,&
            BSURF,RSF,XK1,XK2)

       !C     NOW WE CALCULATE THE FLUXES AT THE MIDPOINTS OF THE LAYERS.

       !C  For original, unscaled, version, use TAUCUM, not TAUCUMP
       DO L=1,L_NLAYRAD-1
          !c     EXPTRM(L) = MIN(TAUMAX,LAMDA(L)*(TAUCUM(2*L+1)-TAUCUM(2*L)))
          EXPTRM(L) = MIN(TAUMAX,LAMDA(L)*(TAUCUMP(2*L+1)-TAUCUMP(2*L)))

          EP = EXP(EXPTRM(L))

          EM    = 1.0/EP
          G4    = 1.0-G3(L)
          DENOM = LAMDA(L)**2 - 1./UBAR0**2

          !C       THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
          !C       THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN
          !C       THE SCATTERING GOES TO ZERO
          !C       PREVENT THIS WITH A IF STATEMENT

          IF ( DENOM .EQ. 0.) THEN
             DENOM=1.E-10
          END IF

          AM = F0PI*W0(L)*(G4   *(G1(L)+1./UBAR0) +G2(L)*G3(L) )/DENOM
          AP = F0PI*W0(L)*(G3(L)*(G1(L)-1./UBAR0) +G2(L)*G4    )/DENOM

          !C       CPMID AND CMMID  ARE THE CPLUS AND CMINUS TERMS EVALUATED
          !C       AT THE MIDDLE OF THE LAYER.

          TAUMID   = TAUCUMP(2*L+1)
          !c       TAUMID   = TAUCUM(2*L+1)

          CPMID = AP*EXP(-TAUMID/UBAR0)
          CMMID = AM*EXP(-TAUMID/UBAR0)

          FMIDP(L) = XK1(L)*EP + GAMA(L)*XK2(L)*EM + CPMID
          FMIDM(L) = XK1(L)*EP*GAMA(L) + XK2(L)*EM + CMMID

          !C       ADD THE DIRECT FLUX TO THE DOWNWELLING TERM

          FMIDM(L)= FMIDM(L)+UBAR0*F0PI*EXP(-MIN(TAUMID,TAUMAX)/UBAR0)

       END DO

       !C     FLUX AT THE Ptop layer

       EP    = 1.0
       EM    = 1.0
       G4    = 1.0-G3(1)
       DENOM = LAMDA(1)**2 - 1./UBAR0**2

       !C     THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
       !C     THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN
       !C     THE SCATTERING GOES TO ZERO
       !C     PREVENT THIS WITH A IF STATEMENT

       IF ( DENOM .EQ. 0.) THEN
          DENOM=1.E-10
       END IF

       AM = F0PI*W0(1)*(G4   *(G1(1)+1./UBAR0) +G2(1)*G3(1) )/DENOM
       AP = F0PI*W0(1)*(G3(1)*(G1(1)-1./UBAR0) +G2(1)*G4    )/DENOM

       !C     CPMID AND CMMID  ARE THE CPLUS AND CMINUS TERMS EVALUATED
       !C     AT THE MIDDLE OF THE LAYER.

       CPMID  = AP
       CMMID  = AM

       FLUXUP = XK1(1)*EP + GAMA(1)*XK2(1)*EM + CPMID
       FLUXDN = XK1(1)*EP*GAMA(1) + XK2(1)*EM + CMMID

       !C     ADD THE DIRECT FLUX TO THE DOWNWELLING TERM

       fluxdn = fluxdn+UBAR0*F0PI*EXP(-MIN(TAUCUMP(2),TAUMAX)/UBAR0)
       !c     fluxdn = fluxdn+UBAR0*F0PI*EXP(-MIN(TAUCUM(1),TAUMAX)/UBAR0)

       !C     This is for the "special" bottom layer, where we take
       !C     DTAU instead of DTAU/2.

       L     = L_NLAYRAD
       EXPTRM(L) = MIN(TAUMAX,LAMDA(L)*(TAUCUMP(L_LEVELS)-&
            TAUCUMP(L_LEVELS-1)))
       !c     EXPTRM(L) = MIN(TAUMAX,LAMDA(L)*(TAUCUM(L_LEVELS)-&
       !c                                      TAUCUM(L_LEVELS-1)))

       EP    = EXP(EXPTRM(L))
       EM    = 1.0/EP
       G4    = 1.0-G3(L)
       DENOM = LAMDA(L)**2 - 1./UBAR0**2

       !C     THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
       !C     THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN
       !C     THE SCATTERING GOES TO ZERO
       !C     PREVENT THIS WITH A IF STATEMENT

       IF ( DENOM .EQ. 0.) THEN
          DENOM=1.E-10
       END IF

       AM = F0PI*W0(L)*(G4   *(G1(L)+1./UBAR0) +G2(L)*G3(L) )/DENOM
       AP = F0PI*W0(L)*(G3(L)*(G1(L)-1./UBAR0) +G2(L)*G4    )/DENOM

       !C     CPMID AND CMMID  ARE THE CPLUS AND CMINUS TERMS EVALUATED
       !C     AT THE MIDDLE OF THE LAYER.

       TAUMID   = MIN(TAUCUMP(L_LEVELS),TAUMAX)
       !c     TAUMID   = MIN(TAUCUM(L_LEVELS),TAUMAX)
       CPMID    = AP*EXP(-MIN(TAUMID,TAUMAX)/UBAR0)
       CMMID    = AM*EXP(-MIN(TAUMID,TAUMAX)/UBAR0)

       FMIDP(L) = XK1(L)*EP + GAMA(L)*XK2(L)*EM + CPMID
       FMIDM(L) = XK1(L)*EP*GAMA(L) + XK2(L)*EM + CMMID
!if (iproc .eq. 21 .and. ilon .eq. 1 .and. ilat .eq. 4)
!   if (fmidp(l) .ne. fmidp(l)) then
!      write(*,*) fmidp(l),fmidm(l),ep, gama(l),xk2(l),em,cpmid
!endif
       !C  Save the diffuse downward flux for TEMPGR calculations

       DIFFV = FMIDM(L)

       !C     ADD THE DIRECT FLUX TO THE DOWNWELLING TERM

       FMIDM(L)= FMIDM(L)+UBAR0*F0PI*EXP(-MIN(TAUMID,TAUMAX)/UBAR0)


     END SUBROUTINE GFLUXV
     ! ----------------------------------------------------------------------

     ! ----------------------------------------------------------------------

     SUBROUTINE DSOLVER(NLX,GAMA,CP,CM,CPM1,CMM1,E1,E2,E3,E4,BTOP,&
          BSURF,RSF,XK1,XK2)

       !C  GCM2.0  Feb 2003
       !C
       !C DOUBLE PRECISION VERSION OF SOLVER

       implicit none

       integer,PARAMETER :: NMAX=201
!      integer :: NL, L, N, LM1, LM2, I
       integer :: NLX, L, N, LM1, LM2, I
!      real :: GAMA(NL),CP(NL),CM(NL),CPM1(NL),CMM1(NL),XK1(NL)
!      real :: XK2(NL),E1(NL),E2(NL),E3(NL),E4(NL),RSF,BSURF,BTOP
       real :: GAMA(NLX),CP(NLX),CM(NLX),CPM1(NLX),CMM1(NLX),XK1(NLX)
       real :: XK2(NLX),E1(NLX),E2(NLX),E3(NLX),E4(NLX),RSF,BSURF,BTOP
       real :: AF(NMAX),BF(NMAX),CF(NMAX),DF(NMAX),XK(NMAX)
!!$C*********************************************************
!!$C* THIS SUBROUTINE SOLVES FOR THE COEFFICIENTS OF THE    *
!!$C* TWO STREAM SOLUTION FOR GENERAL BOUNDARY CONDITIONS   *
!!$C* NO ASSUMPTION OF THE DEPENDENCE ON OPTICAL DEPTH OF   *
!!$C* C-PLUS OR C-MINUS HAS BEEN MADE.                      *
!!$C* NLX     = NUMBER OF LAYERS IN THE MODEL                *
!!$C* CP     = C-PLUS EVALUATED AT TAO=0 (TOP)              *
!!$C* CM     = C-MINUS EVALUATED AT TAO=0 (TOP)             *
!!$C* CPM1   = C-PLUS  EVALUATED AT TAOSTAR (BOTTOM)        *
!!$C* CMM1   = C-MINUS EVALUATED AT TAOSTAR (BOTTOM)        *
!!$C* EP     = EXP(LAMDA*DTAU)                              *
!!$C* EM     = 1/EP                                         *
!!$C* E1     = EP + GAMA *EM                                *
!!$C* E2     = EP - GAMA *EM                                *
!!$C* E3     = GAMA*EP + EM                                 *
!!$C* E4     = GAMA*EP - EM                                 *
!!$C* BTOP   = THE DIFFUSE RADIATION INTO THE MODEL AT TOP  *
!!$C* BSURF  = THE DIFFUSE RADIATION INTO THE MODEL AT      *
!!$C*          THE BOTTOM: INCLUDES EMMISION AND REFLECTION *
!!$C*          OF THE UNATTENUATED PORTION OF THE DIRECT    *
!!$C*          BEAM. BSTAR+RSF*FO*EXP(-TAOSTAR/U0)          *
!!$C* RSF    = REFLECTIVITY OF THE SURFACE                  *
!!$C* XK1    = COEFFICIENT OF THE POSITIVE EXP TERM         *
!!$C* XK2    = COEFFICIENT OF THE NEGATIVE EXP TERM         *
!!$C*********************************************************

       !C======================================================================C

!      L=2*NL
       L=2*NLX

       !C     ************MIXED COEFFICENTS**********
       !C     THIS VERSION AVOIDS SINGULARITIES ASSOC.
       !C     WITH W0=0 BY SOLVING FOR XK1+XK2, AND XK1-XK2.

       AF(1) = 0.0
       BF(1) = GAMA(1)+1.
       CF(1) = GAMA(1)-1.
       DF(1) = BTOP-CMM1(1)
       N     = 0
       LM2   = L-2

       !C     EVEN TERMS

       DO I=2,LM2,2
          N     = N+1
          AF(I) = (E1(N)+E3(N))*(GAMA(N+1)-1.)
          BF(I) = (E2(N)+E4(N))*(GAMA(N+1)-1.)
          CF(I) = 2.0*(1.-GAMA(N+1)**2)
          DF(I) = (GAMA(N+1)-1.) * (CPM1(N+1) - CP(N)) +&
               (1.-GAMA(N+1))* (CM(N)-CMM1(N+1))
       END DO

       N   = 0
       LM1 = L-1
       DO I=3,LM1,2
          N     = N+1
          AF(I) = 2.0*(1.-GAMA(N)**2)
          BF(I) = (E1(N)-E3(N))*(1.+GAMA(N+1))
          CF(I) = (E1(N)+E3(N))*(GAMA(N+1)-1.)
          DF(I) = E3(N)*(CPM1(N+1) - CP(N)) + E1(N)*(CM(N) - CMM1(N+1))
       END DO

!      AF(L) = E1(NL)-RSF*E3(NL)
!      BF(L) = E2(NL)-RSF*E4(NL)
       AF(L) = E1(NLX)-RSF*E3(NLX)
       BF(L) = E2(NLX)-RSF*E4(NLX)
       CF(L) = 0.0
!      DF(L) = BSURF-CP(NL)+RSF*CM(NL)
       DF(L) = BSURF-CP(NLX)+RSF*CM(NLX)
!if (BF(L) .ne. bf(L)) then
!   write(*,*) "BF: ",bf(L),e2(NlX),rsf,e4(NLX)
!endif
       CALL DTRIDGL(L,AF,BF,CF,DF,XK)

       !C     ***UNMIX THE COEFFICIENTS****

!      DO 28 N=1,NL
       DO 28 N=1,NLX
          XK1(N) = XK(2*N-1)+XK(2*N)
          XK2(N) = XK(2*N-1)-XK(2*N)

          !C       NOW TEST TO SEE IF XK2 IS REALLY ZERO TO THE LIMIT OF THE
          !C       MACHINE ACCURACY  = 1 .E -30
          !C       XK2 IS THE COEFFICEINT OF THE GROWING EXPONENTIAL AND MUST
          !C       BE TREATED CAREFULLY

          IF(XK2(N) .EQ. 0.0) GO TO 28

!!!! next line added by ridley and pawlowski
          XK(2*N-1) = max(abs(XK(2*n-1)),1.0e-30)*sign(1.0,XK(2*n-1))


          IF (ABS (XK2(N)/XK(2*N-1)) .LT. 1.E-30) XK2(N)=0.0

28        CONTINUE

        END SUBROUTINE DSOLVER
        ! ----------------------------------------------------------------------

        ! ----------------------------------------------------------------------


        SUBROUTINE DTRIDGL(L,AF,BF,CF,DF,XK)

          !C  GCM2.0  Feb 2003

          !C     DOUBLE PRESCISION VERSION OF TRIDGL

          INTEGER,PARAMETER :: NMAX=201
          INTEGER :: L,I
          REAL :: AF(L),BF(L),CF(L),DF(L),XK(L)
          REAL :: AS(NMAX),DS(NMAX),X,XKB

!!$C*    THIS SUBROUTINE SOLVES A SYSTEM OF TRIDIAGIONAL MATRIX
!!$C*    EQUATIONS. THE FORM OF THE EQUATIONS ARE:
!!$C*    A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = D(I)
!!$C*    WHERE I=1,L  LESS THAN 103.
!!$C* ..............REVIEWED -CP........

          !C======================================================================C

          AS(L) = AF(L)/BF(L)
          DS(L) = DF(L)/BF(L)
!          if (AS(L) .ne. AS(L) .or. ds(L) .ne. ds(l)) then
!             write(*,*) "ne: ", as(l),ds(l),BF(L)
!endif
          DO I=2,L
             X         = 1./(BF(L+1-I) - CF(L+1-I)*AS(L+2-I))
             AS(L+1-I) = AF(L+1-I)*X
             DS(L+1-I) = (DF(L+1-I)-CF(L+1-I)*DS(L+2-I))*X
!             if (ds(L+1-I) .ne. ds(L+1-I))then
!                write(*,*)"dtr: ",I, df(L+1-I),DF(L+1-I),DS(L+2-I),X,BF(L+1-I),cf(L+1-I),AS(L+2-I),af(L+1-I)
!             endif
          END DO

          XK(1)=DS(1)
          DO I=2,L
             XKB   = XK(I-1)
             XK(I) = DS(I)-AS(I)*XKB

END DO


        END SUBROUTINE DTRIDGL

        ! ----------------------------------------------------------------------

        ! ----------------------------------------------------------------------

        SUBROUTINE OPTCI(DTAUI,TAUCUMI,TLEV,PLEV,L_LAYERS,&
             L_LEVELS,L_NLAYRAD,L_NLEVRAD,QREFV,COSBI,WBARI,TAUREF,&
             TMID,PMID,NGWI,QH2O)

!!$C  GCM2.0  Feb 2003
!!$C
!!$C THIS SUBROUTINE SETS THE OPTICAL CONSTANTS IN THE INFRARED
!!$C IT CALCUALTES FOR EACH LAYER, FOR EACH SPECRAL INTERVAL IN THE IR
!!$C LAYER: WBAR, DTAU, COSBAR
!!$C LEVEL: TAU
!!$C
!!$C Qrefv is the extinction coefficient at the reference (visible)
!!$C wavelength - 0.67 microns.
!!$C
!!$C TAUI(L,LW) is the cumulative optical depth at level L (or alternatively
!!$C at the *bottom* of layer L), LW is the spectral wavelength interval.
!!$C
!!$C     TLEV(L) - Temperature at the layer boundary (i.e. level)
!!$C     PLEV(L) - Pressure at the layer boundary (i.e. level)
!!$C     CO2_KI(NT,NP,NW,NG) - IR CO2 k-coefficients
!!$C                           CO2_K(temp,Pres,Waveln,gauss)
!!$C                           currently: CO2_K(7,11,5,17)
!!$C
!!$C----------------------------------------------------------------------C

          use ModPlanet

          implicit none

          integer :: L_LAYERS,L_LEVELS,L_NLAYRAD,L_NLEVRAD

          real :: DTAUI(LL_NLAYRAD,L_NSPECTI,L_NGAUSS)
          real :: DTAUKI(LL_LEVELS+1,L_NSPECTI,L_NGAUSS)
          real :: TAUI(LL_NLEVRAD,L_NSPECTI,L_NGAUSS)
          real :: TAUCUMI(LL_LEVELS,L_NSPECTI,L_NGAUSS)
          real :: TAUGAS, Qrefv
          real :: TLEV(LL_LEVELS), PLEV(LL_LEVELS+1)
          real :: TMID(LL_LEVELS), PMID(LL_LEVELS)
          real :: COSBI(LL_NLAYRAD,L_NSPECTI,L_NGAUSS)
          real :: WBARI(LL_NLAYRAD,L_NSPECTI,L_NGAUSS)
          real :: TAUREF(LL_LEVELS+1)

          integer :: L, NW, NG, K, LK
          integer :: MT(LL_LEVELS), MP(LL_LEVELS)
          real ::  ANS, TAUREFL
          real ::  TAEROS(LL_LEVELS,L_NSPECTI)
          real ::  DPR(LL_LEVELS), U(LL_LEVELS), TAUCLD, TAUAC
          real ::  LCOEF(4), LKCOEF(LL_LEVELS,4)

          ! fraction of zeros in each spectral interval, as a function of T, P

          real :: dt, tt
          real :: taugsurf(L_NSPECTI,L_NGAUSS-1)

          !C  Tlimit:  If the CO2 optical depth (top to the surface) is less than
          !C  this value, we place that Gauss-point into the "zeros" channel.
          !C  Set in driver, passed via common


          integer :: NGWI(L_NSPECTI)

          !C  Water mixing ratio variables

          real ::  QH2O(LL_LEVELS),  WRATIO(LL_LEVELS)
          real ::  KCOEF(4)
          integer :: NH2O(LL_LEVELS)

          !C======================================================================C

          !C  Determine the total gas opacity throughout the column, for each
          !C  spectral interval, NW, and each Gauss point, NG.

          DTAUKI = 0.0

          DO NG=1,L_NGAUSS-1
             do NW=1,L_NSPECTI
                TAUGSURF(NW,NG) = 0.0D0
             end do
          end do

          do K=2,L_LEVELS
             DPR(k) = PLEV(K)-PLEV(K-1)
             U(k)   = Cmk*DPR(k)
!             if (U(k) .lt. 0) then
!                write(*,*) u(k),dpr(k),plev(k), plev(k-1)

!endif
             call tpindex(PMID(K),TMID(K),QH2O(K),pfgasref,tgasref,&
                  LCOEF,MT(K),MP(K),NH2O(K),WRATIO(K))

             do LK=1,4
                LKCOEF(K,LK) = LCOEF(LK)
             end do

             DO NW=1,L_NSPECTI
                TAEROS(K,NW) = TAUREF(K)*Qexti(NW)/Qrefv
             END DO
          end do
!stop
          !C  TAUCLD = is cloud opacity, zero until further notice

          TAUCLD = 0.0

          do NW=1,L_NSPECTI
             ngwi(NW) = L_NGAUSS

             !C  Now fill in the "clear" part of the spectrum (NG = L_NGAUSS)
             !C  Which holds continuum opacity only

             do K=2,L_LEVELS
                DTAUKI(K,NW,L_NGAUSS) = TAEROS(K,NW) + TAUCLD
             end do

             do ng=L_NGAUSS-1,1,-1
                do K=2,L_LEVELS

                   !C           NOW COMPUTE TAUGAS

                   !C  Interpolate between water mixing ratios
                   !C  WRATIO = 0.0 if the requested water amount is equal to, or outside the
                   !C  the range of water amount data.

                   KCOEF(1) = CO2I(MT(K),MP(K),NH2O(K),NW,NG) + WRATIO(K)*&
                        (CO2I(MT(K),MP(K),NH2O(K)+1,NW,NG) -&
                        CO2I(MT(K),MP(K),NH2O(K),NW,NG))

                   KCOEF(2) = CO2I(MT(K),MP(K)+1,NH2O(K),NW,NG) + WRATIO(K)*&
                        (CO2I(MT(K),MP(K)+1,NH2O(K)+1,NW,NG) -&
                        CO2I(MT(K),MP(K)+1,NH2O(K),NW,NG))

                   KCOEF(3) = CO2I(MT(K)+1,MP(K)+1,NH2O(K),NW,NG) + WRATIO(K)*&
                        (CO2I(MT(K)+1,MP(K)+1,NH2O(K)+1,NW,NG) -&
                        CO2I(MT(K)+1,MP(K)+1,NH2O(K),NW,NG))

                   KCOEF(4) = CO2I(MT(K)+1,MP(K),NH2O(K),NW,NG) + WRATIO(K)*&
                        (CO2I(MT(K)+1,MP(K),NH2O(K)+1,NW,NG) -&
                        CO2I(MT(K)+1,MP(K),NH2O(K),NW,NG))

                   !C  Interpolate the CO2 k-coefficients to the requested T,P


                   ANS = LKCOEF(K,1)*KCOEF(1) + LKCOEF(K,2)*KCOEF(2) +&
                        LKCOEF(K,3)*KCOEF(3) + LKCOEF(K,4)*KCOEF(4)


                   TAUGAS          = U(k)*ANS
                   TAUGSURF(NW,NG) = TAUGSURF(NW,NG) + TAUGAS
                   DTAUKI(K,nw,ng) = TAUGAS+TAEROS(K,NW)+TAUCLD
!                   if (nw .eq. 3 .and. ng .eq. 15) then1
!                      write(*,*) dtauki(k,nw,ng),taugas,u(k),ans
!                   endif
    end do

                if(TAUGSURF(NW,NG) .LT. TLIMIT) THEN
                   goto 10
                else
                   NGWI(NW) = NG
                end if

             end do
10           continue

          end do
!stop
          !C  Now the full treatment for the layers, where besides the opacity
          !C  we need to calculate the scattering albedo and asymmetry factors
          !C  for each layer

          DO NW=1,L_NSPECTI

             !C  First, the special "clear" channel

             NG = L_NGAUSS

             DO L=1,L_NLAYRAD
                K              = 2*L+1
                TAUREFL        = (TAUREF(K)+TAUREF(K+1))/Qrefv
                DTAUI(L,nw,ng) = DTAUKI(K,NW,NG)+DTAUKI(K+1,NW,NG)
                if(DTAUI(L,NW,NG) .GT. 1.0E-9) then
                   WBARI(L,nw,ng) = (QSCATI(NW)*TAUREFL)/DTAUI(L,NW,NG)
                else
                   WBARI(L,nw,ng) = 0.0D0
                   DTAUI(L,NW,NG) = 1.0E-9
                endif

                TAUAC = TAEROS(K,NW) + TAUCLD
                if(TAUAC .GT. 0.0) then
                   cosbi(L,NW,NG) = GI(NW)           !change formula to add clouds
                else
                   cosbi(L,NW,NG) = 0.0D0
                end if

             END DO

          END DO

          !C  . . .Now the other Gauss points, if needed.

          DO NW=1,L_NSPECTI
             DO NG=L_NGAUSS-1,NGWI(NW),-1

                DO L=1,L_NLAYRAD
                   K              = 2*L+1
                   TAUREFL        = (TAUREF(K)+TAUREF(K+1))/Qrefv
                   DTAUI(L,nw,ng) = DTAUKI(K,NW,NG)+DTAUKI(K+1,NW,NG)
                   if(DTAUI(L,NW,NG) .GT. 1.0E-9) then
                      WBARI(L,nw,ng) = (QSCATI(NW)*TAUREFL)/DTAUI(L,NW,NG)
                   else
                      WBARI(L,nw,ng) = 0.0D0
                      DTAUI(L,NW,NG) = 1.0E-9
                   endif

                   cosbi(L,NW,NG) = cosbi(L,NW,L_NGAUSS)
                END DO

             END DO

          END DO     ! NW spectral loop

          !C     TOTAL EXTINCTION OPTICAL DEPTHS

          DO NW=1,L_NSPECTI
             NG = L_NGAUSS
             TAUI(1,NW,NG) = 0.0D0
             DO L=1,L_NLAYRAD
                TAUI(L+1,NW,NG) = TAUI(L,NW,NG)+DTAUI(L,NW,NG)
             END DO

             TAUCUMI(1,NW,NG)=0.0D0
             DO K=2,L_LEVELS
                TAUCUMI(K,NW,NG)=TAUCUMI(K-1,NW,NG)+DTAUKI(K,NW,NG)
             END DO

             DO NG=L_NGAUSS-1,NGWI(NW),-1

                TAUI(1,NW,NG)=0.0D0
                DO L=1,L_NLAYRAD
                   TAUI(L+1,NW,NG)=TAUI(L,NW,NG)+DTAUI(L,NW,NG)
                END DO

                TAUCUMI(1,NW,NG)=0.0D0
                DO K=2,L_LEVELS
                   TAUCUMI(K,NW,NG)=TAUCUMI(K-1,NW,NG)+DTAUKI(K,NW,NG)
                END DO

             END DO
          END DO

        END SUBROUTINE OPTCI


        ! ----------------------------------------------------------------------

        ! ----------------------------------------------------------------------

        SUBROUTINE SFLUXI(PLEV,TLEV,DTAUI,TAUCUMI,RSFI,L_LAYERS,L_LEVELS,&
             L_NLAYRAD,L_NLEVRAD,COSBI,WBARI,NFLUXTOPI,FMNETI,&
             fluxupi,fluxdni,NGWI)

          !C  GCM2.0  Feb 2003

          use ModPlanet

          implicit none

          integer :: L_LAYERS,L_LEVELS,L_NLAYRAD,L_NLEVRAD
          integer :: NLEVRAD, L, NW, NG, NTS, NTT

          real :: TLEV(LL_LEVELS), PLEV(LL_LEVELS+1)
          real :: TAUCUMI(LL_LEVELS,L_NSPECTI,L_NGAUSS)
          real :: FMNETI(LL_NLAYRAD)

          real :: DTAUI(LL_NLAYRAD,L_NSPECTI,L_NGAUSS)
          real :: FMUPI(LL_NLAYRAD), FMDI(LL_NLAYRAD)
          real :: COSBI(LL_NLAYRAD,L_NSPECTI,L_NGAUSS)
          real :: WBARI(LL_NLAYRAD,L_NSPECTI,L_NGAUSS)
          real :: NFLUXTOPI
          real :: FTOPUP

          real :: RSFI, TSURF, BSURF, TTOP, BTOP, TAUTOP
          real :: PLANCK, PLTOP
          real :: fluxupi(LL_NLAYRAD), fluxdni(LL_NLAYRAD)

          real :: fzero
          integer :: NGWI(L_NSPECTI)


          !C======================================================================C

          NLEVRAD = L_NLEVRAD

          !C     ZERO THE NET FLUXES

          NFLUXTOPI = 0.0

          DO L=1,L_NLAYRAD
             FMNETI(L)  = 0.0
             FLUXUPI(L) = 0.0
             FLUXDNI(L) = 0.0
          END DO

          !C     WE NOW ENTER A MAJOR LOOP OVER SPECTRAL INTERVALS IN THE INFRARED
          !C     TO CALCULATE THE NET FLUX IN EACH SPECTRAL INTERVAL

          TTOP  = TLEV(2)
          TSURF = TLEV(L_LEVELS)

          NTS   = TSURF*10.0D0-499
          NTT   = TTOP *10.0D0-499

          DO 501 NW=1,L_NSPECTI

             !C       SURFACE EMISSIONS - INDEPENDENT OF GAUSS POINTS

             BSURF = (1.-RSFI)*PLANCKIR(NW,NTS)
             PLTOP = PLANCKIR(NW,NTT)

             !C  If FZEROI(NW) = 1, then the k-coefficients are zero - skip to the
             !C  special Gauss point at the end.

             FZERO = FZEROI(NW)
             IF(FZERO.ge.0.99) goto 40

             DO NG=1,NGWI(NW)-1
                fzero = fzero + (1.0-FZEROI(NW))*GWEIGHT(NG)
             END DO

             DO NG=NGWI(NW),L_NGAUSS-1

                !C         SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE IR
                !C         CALCULATE THE DOWNWELLING RADIATION AT THE TOP OF THE MODEL
                !C         OR THE TOP LAYER WILL COOL TO SPACE UNPHYSICALLY

                TAUTOP = DTAUI(1,NW,NG)*PLEV(2)/(PLEV(4)-PLEV(2))
                BTOP   = (1.0-EXP(-TAUTOP/UBARI))*PLTOP

                !C         WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
                !C         CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
                !C         WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER

                CALL GFLUXI(NLEVRAD,TLEV,NW,DWNI(NW),DTAUI(1,NW,NG),&
                     TAUCUMI(1,NW,NG),&
                     WBARI(1,NW,NG),COSBI(1,NW,NG),RSFI,BTOP,&
                     BSURF,FTOPUP,FMUPI,FMDI,&
                     L_LAYERS,L_LEVELS,L_NLAYRAD,L_NLEVRAD)
!if (ng .eq. 15 .and. nw .eq. 3) stop
                !C         NOW CALCULATE THE CUMULATIVE IR NET FLUX

                NFLUXTOPI = NFLUXTOPI+FTOPUP*DWNI(NW)*GWEIGHT(NG)*&
                     (1.0-FZEROI(NW))

                DO L=1,L_NLEVRAD-1

                   !C           CORRECT FOR THE WAVENUMBER INTERVALS
                   FMNETI(L)  = FMNETI(L)+(FMUPI(L)-FMDI(L))*DWNI(NW)*&
                        GWEIGHT(NG)*(1.0-FZEROI(NW))
                   FLUXUPI(L) = FLUXUPI(L) + FMUPI(L)*DWNI(NW)*GWEIGHT(NG)*&
                        (1.0-FZEROI(NW))
                   FLUXDNI(L) = FLUXDNI(L) + FMDI(L)*DWNI(NW)*GWEIGHT(NG)*&
                        (1.0-FZEROI(NW))

                END DO

             END DO       !End NGAUSS LOOP

40           CONTINUE

             !C      SPECIAL 17th Gauss point

             NG     = L_NGAUSS

             TAUTOP = DTAUI(1,NW,NG)*PLEV(2)/(PLEV(4)-PLEV(2))
             BTOP   = (1.0-EXP(-TAUTOP/UBARI))*PLTOP

             !C      WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
             !C      CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
             !C      WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER

             CALL GFLUXI(NLEVRAD,TLEV,NW,DWNI(NW),DTAUI(1,NW,NG),&
                  TAUCUMI(1,NW,NG),&
                  WBARI(1,NW,NG),COSBI(1,NW,NG),RSFI,BTOP,&
                  BSURF,FTOPUP,FMUPI,FMDI,&
                  L_LAYERS,L_LEVELS,L_NLAYRAD,L_NLEVRAD)

             !C      NOW CALCULATE THE CUMULATIVE IR NET FLUX

             NFLUXTOPI = NFLUXTOPI+FTOPUP*DWNI(NW)*FZERO

             DO L=1,L_NLEVRAD-1

                !C        CORRECT FOR THE WAVENUMBER INTERVALS

                FMNETI(L)  = FMNETI(L)+(FMUPI(L)-FMDI(L))*DWNI(NW)*FZERO
                FLUXUPI(L) = FLUXUPI(L) + FMUPI(L)*DWNI(NW)*FZERO
                FLUXDNI(L) = FLUXDNI(L) + FMDI(L)*DWNI(NW)*FZERO
             END DO


501          CONTINUE      !End Spectral Interval LOOP

             !C *** END OF MAJOR SPECTRAL INTERVAL LOOP IN THE INFRARED****


           END SUBROUTINE SFLUXI


           ! ----------------------------------------------------------------------

           ! ----------------------------------------------------------------------

           SUBROUTINE GFLUXI(NLL,TLEV,NW,DW,DTAU,TAUCUM,W0,COSBAR,&
                RSF,BTOP,BSURF,FTOPUP,FMIDP,FMIDM,&
                L_LAYERS,L_LEVELS,L_NLAYRAD,L_NLEVRAD)

!!$C  GCM2.0  Feb 2003
!!$C
!!$C
!!$C ?? COMBINE ONTO ONE GLFLUX ROUTINE, NOT VIS AND IR
!!$C  THIS SUBROUTINE TAKES THE OPTICAL CONSTANTS AND BOUNDARY CONDITIONS
!!$C  FOR THE INFRARED FLUX AT ONE WAVELENGTH AND SOLVES FOR THE FLUXES AT
!!$C  THE LEVELS. THIS VERSION IS SET UP TO WORK WITH LAYER OPTICAL DEPTHS
!!$C  MEASURED FROM THE TOP OF EACH LAYER.  THE TOP OF EACH LAYER HAS
!!$C  OPTICAL DEPTH ZERO.  IN THIS SUB LEVEL N IS ABOVE LAYER N. THAT IS LAYER N
!!$C  HAS LEVEL N ON TOP AND LEVEL N+1 ON BOTTOM. OPTICAL DEPTH INCREASES
!!$C  FROM TOP TO BOTTOM. SEE C.P. MCKAY, TGM NOTES.
!!$C THE TRI-DIAGONAL MATRIX SOLVER IS DSOLVER AND IS DOUBLE PRECISION SO MANY
!!$C VARIABLES ARE PASSED AS SINGLE THEN BECOME DOUBLE IN DSOLVER
!!$C
!!$C NLL            = NUMBER OF LEVELS (NLAYERS + 1) MUST BE LESS THAT NL (101)
!!$C TLEV(L_LEVELS) = ARRAY OF TEMPERATURES AT GCM LEVELS
!!$C WAVEN          = WAVELENGTH FOR THE COMPUTATION
!!$C DW             = WAVENUMBER INTERVAL
!!$C DTAU(NLAYER)   = ARRAY OPTICAL DEPTH OF THE LAYERS
!!$C W0(NLEVEL)     = SINGLE SCATTERING ALBEDO
!!$C COSBAR(NLEVEL) = ASYMMETRY FACTORS, 0=ISOTROPIC
!!$C UBARI          = AVERAGE ANGLE, MUST BE EQUAL TO 0.5 IN IR
!!$C RSF            = SURFACE REFLECTANCE
!!$C BTOP           = UPPER BOUNDARY CONDITION ON IR INTENSITY (NOT FLUX)
!!$C BSURF          = SURFACE EMISSION = (1-RSFI)*PLANCK, INTENSITY (NOT FLUX)
!!$C FP(NLEVEL)     = UPWARD FLUX AT LEVELS
!!$C FM(NLEVEL)     = DOWNWARD FLUX AT LEVELS
!!$C FMIDP(NLAYER)  = UPWARD FLUX AT LAYER MIDPOINTS
!!$C FMIDM(NLAYER)  = DOWNWARD FLUX AT LAYER MIDPOINTS
!!$C
!!$C----------------------------------------------------------------------C

             use ModPlanet

             implicit none

             INTEGER :: L_LAYERS,L_LEVELS,L_NLAYRAD,L_NLEVRAD
!            integer,PARAMETER :: NL=101 ! MUST BE LARGER THAN NLEVEL
             integer,PARAMETER :: NL3=101 ! MUST BE LARGER THAN NLEVEL


!            INTEGER :: NLL, NLAYER, L, NW, NT, NT2
             INTEGER :: NLL, NLAYER, L, NW, NT3, NT2
             REAL ::  TERM, CPMID, CMMID
             REAL ::  EM,EP
             REAL ::  COSBAR(LL_NLAYRAD), W0(LL_NLAYRAD), DTAU(LL_NLAYRAD)
             REAL ::  TAUCUM(LL_LEVELS), DTAUK
             REAL ::  TLEV(LL_LEVELS)
             REAL ::  WAVEN, DW, RSF
             REAL ::  BTOP, BSURF, FMIDP(LL_NLAYRAD), FMIDM(LL_NLAYRAD)
             REAL ::  B0(NL3),B1(NL3),ALPHA(NL3),LAMDA(NL3),XK1(NL3),XK2(NL3)
             REAL ::  GAMA(NL3),CP(NL3),CM(NL3),CPM1(NL3),CMM1(NL3),E1(NL3),E2(NL3)
             REAL ::  E3(NL3),E4(NL3)
             REAL ::  TAUMAX

             REAL ::  FTOPUP, FLUXUP, FLUXDN


             DATA TAUMAX / L_TAUMAX /

             !C======================================================================C

             !C     WE GO WITH THE HEMISPHERIC CONSTANT APPROACH IN THE INFRARED
    
             IF (NLL .GT. NL3) STOP 'PARAMETER NL TOO SMALL IN GLUFV'

             NLAYER = L_NLAYRAD

             DO L=1,L_NLAYRAD-1
                ALPHA(L) = SQRT( (1.0-W0(L))/(1.0-W0(L)*COSBAR(L)) )
                LAMDA(L) = ALPHA(L)*(1.0-W0(L)*COSBAR(L))/UBARI

                NT2   = TLEV(2*L+2)*10.0D0-499
!               NT    = TLEV(2*L)*10.0D0-499
                NT3   = TLEV(2*L)*10.0D0-499

!               B1(L) = (PLANCKIR(NW,NT2)-PLANCKIR(NW,NT))/DTAU(L)
                B1(L) = (PLANCKIR(NW,NT2)-PLANCKIR(NW,NT3))/DTAU(L)
!               B0(L) = PLANCKIR(NW,NT)
                B0(L) = PLANCKIR(NW,NT3)
             END DO

             !C     Take care of special lower layer

             L        = L_NLAYRAD
             ALPHA(L) = SQRT( (1.0-W0(L))/(1.0-W0(L)*COSBAR(L)) )
             LAMDA(L) = ALPHA(L)*(1.0-W0(L)*COSBAR(L))/UBARI

!            NT    = TLEV(2*L+1)*10.0D0-499
             NT3   = TLEV(2*L+1)*10.0D0-499
             NT2   = TLEV(2*L)*10.0D0-499
!            B1(L) = (PLANCKIR(NW,NT)-PLANCKIR(NW,NT2))/DTAU(L)
             B1(L) = (PLANCKIR(NW,NT3)-PLANCKIR(NW,NT2))/DTAU(L)
             B0(L) = PLANCKIR(NW,NT2)

             DO L=1,L_NLAYRAD
                GAMA(L) = (1.0-ALPHA(L))/(1.0+ALPHA(L))
                TERM    = UBARI/(1.0-W0(L)*COSBAR(L))

                !C       CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
                !C       BOTTOM OF THE LAYER.  THAT IS AT DTAU OPTICAL DEPTH

                CP(L) = B0(L)+B1(L)*DTAU(L) +B1(L)*TERM
                CM(L) = B0(L)+B1(L)*DTAU(L) -B1(L)*TERM

                !C       CPM1 AND CMM1 ARE THE CPLUS AND CMINUS TERMS EVALUATED
                !C       AT THE TOP OF THE LAYER, THAT IS ZERO OPTICAL DEPTH

                CPM1(L) = B0(L)+B1(L)*TERM
                CMM1(L) = B0(L)-B1(L)*TERM
             END DO

             !C     NOW CALCULATE THE EXPONENTIAL TERMS NEEDED
             !C     FOR THE TRIDIAGONAL ROTATED LAYERED METHOD
             !C     WARNING IF DTAU(J) IS GREATER THAN ABOUT 35 (VAX)
             !C     WE CLIP IT TO AVOID OVERFLOW.

             DO L=1,L_NLAYRAD

                !C       CLIP THE EXPONENTIAL HERE.

                EP    = EXP( MIN((LAMDA(L)*DTAU(L)),TAUMAX))
                EM    = 1.0/EP
                E1(L) = EP+GAMA(L)*EM
                E2(L) = EP-GAMA(L)*EM
                E3(L) = GAMA(L)*EP+EM
                E4(L) = GAMA(L)*EP-EM
             END DO

             !c     B81=BTOP  ! RENAME BEFORE CALLING DSOLVER - used to be to set
             !c     B82=BSURF ! them to real*8 - but now everything is real*8
             !c     R81=RSF   ! so this may not be necessary

             !C     Double precision tridiagonal solver

             CALL DSOLVER(NLAYER,GAMA,CP,CM,CPM1,CMM1,E1,E2,E3,E4,BTOP,&
                  BSURF,RSF,XK1,XK2)

             !C     NOW WE CALCULATE THE FLUXES AT THE MIDPOINTS OF THE LAYERS.

             DO L=1,L_NLAYRAD-1
                DTAUK = TAUCUM(2*L+1)-TAUCUM(2*L)
                EP    = EXP(MIN(LAMDA(L)*DTAUK,TAUMAX)) ! CLIPPED EXPONENTIAL
                EM    = 1.0/EP
                TERM  = UBARI/(1.-W0(L)*COSBAR(L))

                !C       CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
                !C       BOTTOM OF THE LAYER.  THAT IS AT DTAU  OPTICAL DEPTH

                CPMID    = B0(L)+B1(L)*DTAUK +B1(L)*TERM
                CMMID    = B0(L)+B1(L)*DTAUK -B1(L)*TERM
                FMIDP(L) = XK1(L)*EP + GAMA(L)*XK2(L)*EM + CPMID
                FMIDM(L) = XK1(L)*EP*GAMA(L) + XK2(L)*EM + CMMID

                !C       FOR FLUX WE INTEGRATE OVER THE HEMISPHERE TREATING INTENSITY CONSTANT

                FMIDP(L) = FMIDP(L)*PI
                FMIDM(L) = FMIDM(L)*PI

             END DO

             !C     And now, for the special bottom layer

             L    = L_NLAYRAD

             EP   = EXP(MIN((LAMDA(L)*DTAU(L)),TAUMAX)) ! CLIPPED EXPONENTIAL
             EM   = 1.0/EP
             TERM = UBARI/(1.-W0(L)*COSBAR(L))

             !C     CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
             !C     BOTTOM OF THE LAYER.  THAT IS AT DTAU  OPTICAL DEPTH

             CPMID    = B0(L)+B1(L)*DTAU(L) +B1(L)*TERM
             CMMID    = B0(L)+B1(L)*DTAU(L) -B1(L)*TERM
             FMIDP(L) = XK1(L)*EP + GAMA(L)*XK2(L)*EM + CPMID
             FMIDM(L) = XK1(L)*EP*GAMA(L) + XK2(L)*EM + CMMID

             !C     FOR FLUX WE INTEGRATE OVER THE HEMISPHERE TREATING INTENSITY CONSTANT

             FMIDP(L) = FMIDP(L)*PI
             FMIDM(L) = FMIDM(L)*PI

             !C     FLUX AT THE PTOP LEVEL

             EP   = 1.0
             EM   = 1.0
             TERM = UBARI/(1.0-W0(1)*COSBAR(1))

             !C     CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
             !C     BOTTOM OF THE LAYER.  THAT IS AT DTAU  OPTICAL DEPTH

             CPMID  = B0(1)+B1(1)*TERM
             CMMID  = B0(1)-B1(1)*TERM

             FLUXUP = XK1(1)*EP + GAMA(1)*XK2(1)*EM + CPMID
             FLUXDN = XK1(1)*EP*GAMA(1) + XK2(1)*EM + CMMID

             !C     FOR FLUX WE INTEGRATE OVER THE HEMISPHERE TREATING INTENSITY CONSTANT

             FTOPUP = (FLUXUP-FLUXDN)*PI


           END SUBROUTINE GFLUXI


           ! ----------------------------------------------------------------------

           subroutine set_planet_defaults

             use ModInputs

             DtLTERadiation = 900.0

             return

           end subroutine set_planet_defaults

           subroutine planet_limited_fluxes(iBlock)
             !! Do Nothing
           end subroutine planet_limited_fluxes

!***********************************************************************
!  nlte_tcool supporting subroutines:  November 2017
!***********************************************************************

subroutine NLTEdlvr11_ZGRID (n_gcm, &
           p_gcm, t_gcm, z_gcm, co2vmr_gcm, n2vmr_gcm, &
           covmr_gcm, o3pvmr_gcm, mmean_gcm,cpnew_gcm, &
           nl_cts_real, nzy_cts_real )

!***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     Arguments
      integer :: n_gcm             ! I
      real :: p_gcm(n_gcm), t_gcm(n_gcm) ! I
      real :: co2vmr_gcm(n_gcm), n2vmr_gcm(n_gcm) ! I
      real :: covmr_gcm(n_gcm), o3pvmr_gcm(n_gcm) ! I
      real :: z_gcm(n_gcm)         ! I
      real :: mmean_gcm(n_gcm)     ! I
      real :: cpnew_gcm(n_gcm)     ! I
      integer ::   nl_cts_real, nzy_cts_real ! O
      real :: zaux_gcm(n_gcm)

!     local variables
      integer :: i, iz
      real  :: distancia, meanm, gz, Hkm
      real  :: zmin, zmax
      real  :: mmean_nlte(n_gcm),cpnew_nlte(n_gcm)
      real  :: gg,masa,radio,kboltzman

!     functions
      external :: hrkday_convert
      real  :: 	hrkday_convert

!***********************************************************************

!     Define el working grid para MZ1D (NL, ZL, ZMIN)
!     y otro mas fino para M.Curtis (NZ, ZX, ZXMIN = ZMIN)
!     Tambien el working grid para MZESC110 (NL_cts, ZL_cts, ZMIN_cts=??)
!     Para ello hace falta una z de ref del GCM, que voy a suponer la inferior

!     Primero, construimos escala z_gcm

!      zaux_gcm(1) = z_gcm(1)             ! [km]
!      gg=6.67259e-8
!      masa=6.4163e26
!      radio=3390.
!      kboltzman=1.381e-16

!      do iz = 2, n_gcm
!         distancia = ( radio + zaux_gcm(iz-1) )*1.e5
!         gz = gg * masa / ( distancia * distancia )
!         Hkm = 0.5*( t_gcm(iz)+t_gcm(iz-1) ) /
!     $        ( mmean_gcm(iz)/6.023e23 * gz )
!         Hkm = kboltzman * Hkm *1e-5 ! [km]
!         zaux_gcm(iz) = zaux_gcm(iz-1) -
!     $        Hkm * log( p_gcm(iz)/p_gcm(iz-1) )
!      enddo


!     Segundo, definimos los límites de los 2 modelos de NLTE.
!     NLTE model completo: indices [jlowerboundary,jtopboundary]
!     NLTE CTS : indices [jbotCTS,jtopCTS] donde jbotCTS = jtopboundary-2

!!!!!!!!!Primero el NLTE completo  !!!!!!!!

                                ! Bottom boundary for NLTE model :
                                !   Pbot_atm = 2e-2 mb = 1.974e-5 atm , lnp(nb)=9.9   (see mz1d.par)

      
      jlowerboundary = 1
      do while ( p_gcm(jlowerboundary) .gt. Pbottom_atm )
         jlowerboundary = jlowerboundary + 1
         if (jlowerboundary .gt. n_gcm) then
            write (*,*) 'Error in lower boundary pressure.'
            write (*,*) ' p_gcm too low or wrong. '
	    write (*,*) ' p_gcm, Pbottom_atm =', &
                 p_gcm(n_gcm), Pbottom_atm
            stop ' Check input value "p_gcm" or modify "Pbottom_atm" '
         endif
      enddo

                                ! Top boundary for NLTE model :
                                !   Ptop_atm = 1e-9 atm                          (see mz1d.par)
      jtopboundary = jlowerboundary
      !write(*,*) "Ptop_atm:", Ptop_atm, pbottom_atm
      do while ( p_gcm(jtopboundary) .gt. Ptop_atm )
         jtopboundary = jtopboundary + 1
         if (jtopboundary .gt. n_gcm) then
            !write (*,*) '!!!!!!!! Warning in top boundary pressure. '
            !write (*,*) ' Ptop_atm too high for p_gcm. '
            !write (*,*) ' p_gcm, Ptop_atm =',  &
            !     p_gcm(n_gcm), Ptop_atm
            !write (*,*) '!!!!!!!! NLTE upper boundary modified '// &
            !     ' to match p_gcm'
            jtopboundary=n_gcm
            goto 5000
         endif
      enddo
 5000 continue

                                ! Grid steps

      zmin = z_gcm(jlowerboundary)
      zmax = z_gcm(jtopboundary)
      deltaz = (zmax-zmin) / (nl-1)
      do i=1,nl
         zl(i) = zmin + (i-1) * deltaz
      enddo


      ! Creamos el perfil del NLTE modelo completo interpolando   
      call interhunt (    pl,zl,nl,      p_gcm,z_gcm,n_gcm, 2) ! [atm]
!     call interhunt5veces &
!               ( t, co2vmr, n2vmr, covmr, o3pvmr, &
!               zl, nl, &
!               t_gcm, co2vmr_gcm, n2vmr_gcm, covmr_gcm, o3pvmr_gcm, &
!               z_gcm, n_gcm, 1)
      call interhunt5veces &
                ( tl, co2vmr, n2vmr, covmr, o3pvmr, &
                zl, nl, &
                t_gcm, co2vmr_gcm, n2vmr_gcm, covmr_gcm, o3pvmr_gcm, &
                z_gcm, n_gcm, 1)
      call interhunt ( mmean_nlte,zl,nl,mmean_gcm,z_gcm,n_gcm,1)
      call interhunt ( cpnew_nlte,zl,nl,cpnew_gcm,z_gcm,n_gcm,1)

      do i = 1, nl
!        nt(i) = 7.339e+21 * pl(i) / t(i) ! --> [cm-3]
!        nt(i) = 7.339e+21 * pl(i) / tl(i) ! --> [cm-3]
         ntl(i) = 7.339e+21 * pl(i) / tl(i) ! --> [cm-3]
         co2(i) = ntl(i) * co2vmr(i)
         n2(i) = ntl(i) * n2vmr(i)
         co(i) = ntl(i) * covmr(i)
         o3p(i) = ntl(i) * o3pvmr(i)
!     hrkday_factor(i) =  hrkday_convert( t(i),
!     $        	  co2vmr(i), o3pvmr(i), n2vmr(i), covmr(i) )
!     hrkday_factor(i) =  hrkday_convert( tl(i),
!     $        	  co2vmr(i), o3pvmr(i), n2vmr(i), covmr(i) )
         hrkday_factor(i) = hrkday_convert(mmean_nlte(i) &
              ,cpnew_nlte(i))
      enddo

                                !  Comprobar que las temps no se salen del grid del histograma

      do i=1,nl
!        if (t(i) .gt. 500.0) then
         if (tl(i) .gt. 500.0) then
            write (*,*) '!!!! WARNING    Temp higher than Histogram.'
            write (*,*) ' Histogram will be extrapolated. '
!           write (*,*) ' i, t(i), pl(i) =', i, t(i), pl(i)
            write (*,*) ' i, tl(i), pl(i) =', i, tl(i), pl(i)
         endif
!        if (t(i) .lt. 50.0) then
         if (tl(i) .lt. 50.0) then
            write (*,*) '!!!! WARNING    Temp lower than Histogram.'
            write (*,*) ' Histogram will be extrapolated. '
!           write (*,*) ' i, t(i), pl(i) =', i, t(i), pl(i)
            write (*,*) ' i, tl(i), pl(i) =', i, tl(i), pl(i)
         endif
      enddo

                                !  Fine grid for transmittance calculations

      zmin = z_gcm(jlowerboundary)
      zmax = z_gcm(jtopboundary)
      deltazy = (zmax-zmin) / (nzy-1)
      do i=1,nzy
         zy(i) = zmin + (i-1) * deltazy
      enddo
      call interhunt (    py,zy,nzy,      p_gcm,z_gcm,n_gcm, 2) ! [atm]

      call interhunt2veces ( ty,co2y, zy,nzy, &
           t_gcm,co2vmr_gcm, z_gcm,n_gcm, 1)

      do i=1,nzy
         nty(i) = 7.339e+21 * py(i) / ty(i) ! --> [cm-3]
         co2y(i) = co2y(i) * nty(i)
      enddo


!!!!!!!!!Segundo, el NLTE - CTS  !!!!!!!!

                                ! Grid steps
      deltaz_cts = deltaz
      zl_cts(1) = zl(nl-1)
      nl_cts_real = 1
      do i=2,nl_cts
         zl_cts(i) = zl_cts(1) + (i-1)*deltaz_cts
         if (zl_cts(i) .gt. z_gcm(n_gcm)) then
!            write (*,*) '!!!!!!!! Warning in top CTS layers. '
!            write (*,*) ' zl_Cts too high for z_gcm. '
!            write (*,*) ' z_gcm, zl_cts(i), i =',
!     $           z_gcm(n_gcm), zl_cts(i), i
!            write (*,*) '!!!!!!!! NLTE-CTS upper boundary modified '//
!     $           ' to match z_gcm'
            nl_cts_real=i-1
!            write (*,*) '  Original,Real NL_CTS=', nl_cts,nl_cts_real
            goto 6000
         endif
      enddo
      nl_cts_real = nl_cts
 6000 continue

                                ! Creamos perfil por interpolacion

      call interhuntlimits ( pl_cts,zl_cts,nl_cts, 1,nl_cts_real, &
                p_gcm,z_gcm,n_gcm, 2)
      call interhuntlimits5veces &
                ( t_cts, co2vmr_cts, n2vmr_cts, covmr_cts, o3pvmr_cts, &
                zl_cts, nl_cts, &
                1,nl_cts_real,  &
                t_gcm, co2vmr_gcm, n2vmr_gcm, covmr_gcm, o3pvmr_gcm, &
                z_gcm, n_gcm, 1)
      call interhuntlimits( cpnew_cts,zl_cts,nl_cts,1,nl_cts_real, &
                cpnew_gcm,z_gcm,n_gcm, 1)
      call interhuntlimits( mmean_cts,zl_cts,nl_cts,1,nl_cts_real, &
                mmean_gcm,z_gcm,n_gcm, 1)

      do i = 1, nl_cts_real
         nt_cts(i) = 7.339e+21 * pl_cts(i) / t_cts(i) ! --> [cm-3]
         co2_cts(i) = nt_cts(i) * co2vmr_cts(i)
         n2_cts(i) = nt_cts(i) * n2vmr_cts(i)
         co_cts(i) = nt_cts(i) * covmr_cts(i)
         o3p_cts(i) = nt_cts(i) * o3pvmr_cts(i)
         hrkday_factor_cts(i) =  hrkday_convert( mmean_cts(i) &
             ,cpnew_cts(i) )
      enddo

                                !  Comprobar que las temps no se salen del grid del histograma
      do i=1,nl_cts_real
         if (t_cts(i) .gt. thist_stored(1,mm_stored(1))) then
            write (*,*) '!!!! WARNING    Temp higher than Histogram.'
            write (*,*) ' ZGRID: Histogram will be extrapolated. '
!           write (*,*) ' i, t(i), pl(i) =', i, t_cts(i), pl_cts(i)
            write (*,*) ' i, tl(i), pl(i) =', i, t_cts(i), pl_cts(i)
         endif
         if (t_cts(i) .lt. 50.0) then
            write (*,*) '!!!! WARNING    Temp lower than Histogram.'
            write (*,*) ' ZGRID: Histogram will be extrapolated. '
!           write (*,*) ' i, t(i), pl(i) =', i, t_cts(i), pl_cts(i)
            write (*,*) ' i, tl(i), pl(i) =', i, t_cts(i), pl_cts(i)
         endif
      enddo

                                ! Calculo del indice maximo del GCM hasta donde llega el NLTE-CTS
      jtopCTS = jtopboundary
      do while ( p_gcm(jtopCTS) .gt. pl_cts(nl_cts_real) )
         jtopCTS = jtopCTS + 1
         if (jtopCTS .gt. n_gcm) then
            !write (*,*) '!!!!!!!! Warning in top boundary pressure. '
            !write (*,*) ' Ptop_NLTECTS too high for p_gcm. '
            !write (*,*) ' p_gcm, Ptop_NLTECTS =',  &
            !          p_gcm(n_gcm), pl_cts(nl_cts_real)
            !write (*,*) '!!!!!!!! NLTE-CTS upper boundary modified '// &
            !          ' to match p_gcm'
            jtopCTS=n_gcm
            goto 7000
         endif
      enddo
 7000 continue

                                !  Fine grid for transmittance calculations

      deltazy_cts = 0.25*deltaz_cts ! Comprobar el factor 4 en mz1d.par
      do i=1,nzy_cts
         zy_cts(i) = zl_cts(1) + (i-1) * deltazy_cts
      enddo
      nzy_cts_real = (nl_cts_real - 1)*4 + 1
      call interhuntlimits ( py_cts,zy_cts,nzy_cts, 1,nzy_cts_real, &
                p_gcm, z_gcm, n_gcm,   2) ! [atm]
      call interhuntlimits2veces &
                ( ty_cts,co2y_cts, zy_cts,nzy_cts,  1,nzy_cts_real, &
                t_gcm,co2vmr_gcm, z_gcm,n_gcm, 1)

      do i=1,nzy_cts_real
         nty_cts(i) = 7.339e+21 * py_cts(i) / ty_cts(i) ! --> [cm-3]
         co2y_cts(i) = co2y_cts(i) * nty_cts(i)
      enddo

!     write (*,*) '  NL = ', NL
!     write (*,*) '  Original,Real NL_CTS=', nl_cts,nl_cts_real
!     write (*,*) '  Original,Real NZY_CTS =', nzy_cts,nzy_cts_real

end subroutine NLTEdlvr11_ZGRID


!**********************************************************************

subroutine NLTEdlvr11_CZALU(ierr,varerr)

!***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!************common variables and constants*****************************

!     Arguments

      integer :: ierr
      real*8 :: varerr


!     local variables

!     matrixes and vectors

      real*8 :: e110(nl), e210(nl), e310(nl), e410(nl)
      real*8 :: e121(nl)
      real*8 :: f1(nl,nl)

      real*8 :: cax1(nl,nl), cax2(nl,nl), cax3(nl,nl)
      real*8 :: v1(nl), v2(nl), v3(nl)
      real*8 :: alf11(nl,nl), alf12(nl,nl)
      real*8 :: alf21(nl,nl), alf31(nl,nl), alf41(nl,nl)
      real*8 :: a11(nl), a1112(nl,nl)
      real*8 :: a1121(nl,nl), a1131(nl,nl), a1141(nl,nl)
      real*8 ::  a21(nl), a2131(nl,nl), a2141(nl,nl)
      real*8 ::	a2111(nl,nl), a2112(nl,nl)
      real*8 ::  a31(nl), a3121(nl,nl), a3141(nl,nl)
      real*8 ::	a3111(nl,nl), a3112(nl,nl)
      real*8 :: a41(nl), a4121(nl,nl), a4131(nl,nl)
      real*8 ::	a4111(nl,nl), a4112(nl,nl)
      real*8 :: a12(nl), a1211(nl,nl)
      real*8 ::	a1221(nl,nl), a1231(nl,nl), a1241(nl,nl)

      real*8 ::  aalf11(nl,nl),aalf21(nl,nl), &
           aalf31(nl,nl),aalf41(nl,nl)
      real*8 :: aa11(nl), aa1121(nl,nl), aa1131(nl,nl), aa1141(nl,nl)
      real*8 :: aa21(nl), aa2111(nl,nl), aa2131(nl,nl), aa2141(nl,nl)
      real*8 :: aa31(nl), aa3111(nl,nl), aa3121(nl,nl), aa3141(nl,nl)
      real*8 :: aa41(nl), aa4111(nl,nl), aa4121(nl,nl), aa4131(nl,nl)
      real*8 :: aa1211(nl,nl),aa1221(nl,nl), &
           aa1231(nl,nl),aa1241(nl,nl)
      real*8 ::  aa1112(nl,nl),aa2112(nl,nl), &
           aa3112(nl,nl),aa4112(nl,nl)

      real*8 :: aaalf11(nl,nl), aaalf31(nl,nl), aaalf41(nl,nl)
      real*8 :: aaa11(nl),aaa1131(nl,nl),aaa1141(nl,nl)
      real*8 :: aaa31(nl),aaa3111(nl,nl),aaa3141(nl,nl)
      real*8 :: aaa41(nl),aaa4111(nl,nl),aaa4131(nl,nl)

      real*8 :: aaaalf11(nl,nl),aaaalf41(nl,nl)
      real*8 :: aaaa11(nl),aaaa1141(nl,nl)
      real*8 ::aaaa41(nl),aaaa4111(nl,nl)


!     populations
      real*8 :: n10(nl), n11(nl), n12(nl)
      real*8 :: n20(nl), n21(nl)
      real*8 :: n30(nl), n31(nl)
      real*8 :: n40(nl), n41(nl)

!     productions and loses
      real*8 :: d19b1,d19c1
      real*8 :: d19bp1,d19cp1
      real*8 :: d19c2
      real*8 :: d19cp2
      real*8 :: d19c3
      real*8 :: d19cp3
      real*8 :: d19c4
      real*8 :: d19cp4

      real*8 :: l11, l12, l21, l31, l41
      real*8 :: p11, p12, p21, p31, p41
      real*8 :: p1112, p1211, p1221, p1231, p1241
      real*8 :: p1121, p1131, p1141
      real*8 :: p2111, p2112, p2131, p2141
      real*8 :: p3111, p3112, p3121, p3141
      real*8 :: p4111, p4112, p4121, p4131

      real*8 :: pl11, pl12, pl21, pl31, pl41

      real*8 :: minvt11, minvt21, minvt31, minvt41

!     local constants and indexes

      real*8 :: co2t, o3pdbl, codble, n2dble
      real*8 :: a12_einst(nl)
      real*8 :: a21_einst(nl), a31_einst(nl), a41_einst(nl)
      real :: tsurf

      integer :: i, isot

!     external functions and subroutines

      external :: planckdp
      real*8 ::	planckdp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!start program

      ierr = 0
      varerr = 0.d0

      call zero4v( aa11, aa21, aa31, aa41, nl)
      call zero4m( aa1121, aa1131, aa1141, aalf11, nl)
      call zero4m( aa2111, aa2131, aa2141, aalf21, nl)
      call zero4m( aa3111, aa3121, aa3141, aalf31, nl)
      call zero4m( aa4111, aa4121, aa4131, aalf41, nl)
      call zero4m( aa1112, aa2112, aa3112, aa4112, nl)
      call zero4m( aa1211, aa1221, aa1231, aa1241, nl)
      call zero3v( aaa41, aaa31, aaa11, nl )
      call zero3m( aaa4111, aaa4131, aaalf41, nl)
      call zero3m( aaa3111, aaa3141, aaalf31, nl)
      call zero3m( aaa1131, aaa1141, aaalf11, nl)
      call zero2v( aaaa11, aaaa41, nl )
      call zero2m( aaaa1141, aaaalf11, nl)
      call zero2m( aaaa4111, aaaalf41, nl)

      call zero2v (vt11,vt12,nl)
      call zero3v (vt21,vt31,vt41,nl)
      call zero2v (hr110,hr121,nl)
      call zero3v (hr210,hr310,hr410,nl)
      call zero2v (sl110,sl121,nl)
      call zero3v (sl210,sl310,sl410,nl)

      call zero4v (el11,el21,el31,el41,nl)
      call zero4v (e110,e210,e310,e410,nl)
      call zero2v (el12,e121,nl)

      call zero3m (cax1,cax2,cax3,nl)
      f1(1:nl,1:nl)=0.d0
!      call zerom (f1,nl)

      call zero3v (v1,v2,v3,nl)

      call zero4m (alf11,alf21,alf31,alf41,nl)
      alf12(1:nl,1:nl)=0.d0
!      call zerom (alf12,nl)
      call zero2v (a11,a12,nl)
      call zero3v (a21,a31,a41,nl)

      call zero3m (a1121,a1131,a1141,nl)
      a1112(1:nl,1:nl)=0.d0
!      call zerom (a1112,nl)

      call zero3m (a1221,a1231,a1241,nl)
      a1211(1:nl,1:nl)=0.d0
!      call zerom (a1211,nl)

      call zero2m (a2111,a2112,nl)
      call zero2m (a2131,a2141,nl)
      call zero2m (a3111,a3112,nl)
      call zero2m (a3121,a3141,nl)
      call zero2m (a4111,a4112,nl)
      call zero2m (a4121,a4131,nl)

      call zero2v (n11,n12,nl)
      call zero3v (n21,n31,n41,nl)

      nu11 = dble(nu(1,1))
      nu12 = dble(nu(1,2))
      nu121 =  nu12-nu11
      nu21 =  dble(nu(2,1))
      nu31 =  dble(nu(3,1))
      nu41 =  dble(nu(4,1))

!
!
      do i=1,nl
         n10(i) = dble( co2(i) * imr(1) )
         n20(i) = dble( co2(i) * imr(2) )
         n30(i) = dble( co2(i) * imr(3) )
         n40(i) = dble( co2(i) * imr(4) )
         if ( input_cza.ge.1 ) then
	    n11(i) = n10(i) *2.d0 *exp( -ee*nu11/v626t1(i) )
	    n21(i) = n20(i) *2.d0 *exp( -ee*nu21/v628t1(i) )
	    n31(i) = n30(i) *2.d0* exp( -ee*nu31/v636t1(i) )
	    n41(i) = n40(i) *2.d0* exp( -ee*nu41/v627t1(i) )
         end if
      enddo

!
!     curtis matrix calculation
!
      call zero3m (c210,c310,c410, nl)

      if ( input_cza.ge.1 ) then

         if (itt_cza.eq.15 ) then

	    call MZMC121

         elseif (itt_cza.eq.13) then

!            call zerom ( c121, nl )
            c121(1:nl,1:nl)=0.d0
            call MZESC121
            call MZTVC121!( ierr,varerr )
            !if (ierr .gt. 0) call ERRORS (ierr,varerr)

         endif

      endif

                                ! Lower Boundary
!     tsurf = t(1)
      tsurf = tl(1)
      do i=1,nl
         sl110(i) = vc110(i) * planckdp( tsurf, nu11 )
         sl210(i) = vc210(i) * planckdp( tsurf, nu21 )
         sl310(i) = vc310(i) * planckdp( tsurf, nu31 )
         sl410(i) = vc410(i) * planckdp( tsurf, nu41 )
      end do
      if (input_cza.ge.1) then
         do i=1,nl
	    sl121(i) = vc121(i) * planckdp( tsurf, nu121 )
         end do
      endif



      do 4,i=nl,1,-1            !----------------------------------------------

         co2t = dble( co2(i) *(imr(1)+imr(3)+imr(2)+imr(4)) )
         o3pdbl = dble( o3p(i) )
         n2dble = dble( n2(i) )
         codble = dble ( co(i) )

         call GETK_dlvr11 ( tl(i) )

                                ! V-T productions and losses V-T

         isot = 1
         d19b1 = k19ba(isot)*co2t + k19bb(isot)*n2dble &
              + k19bc(isot)*codble
         d19c1 = k19ca(isot)*co2t + k19cb(isot)*n2dble &
              + k19cc(isot)*codble
         d19bp1 = k19bap(isot)*co2t + k19bbp(isot)*n2dble &
              + k19bcp(isot)*codble
         d19cp1 = k19cap(isot)*co2t + k19cbp(isot)*n2dble &
              + k19ccp(isot)*codble
         isot = 2
         d19c2 = k19ca(isot)*co2t + k19cb(isot)*n2dble &
              + k19cc(isot)*codble
         d19cp2 = k19cap(isot)*co2t + k19cbp(isot)*n2dble &
              + k19ccp(isot)*codble
         isot = 3
         d19c3 = k19ca(isot)*co2t + k19cb(isot)*n2dble &
              + k19cc(isot)*codble
         d19cp3 = k19cap(isot)*co2t + k19cbp(isot)*n2dble &
              + k19ccp(isot)*codble
         isot = 4
         d19c4 = k19ca(isot)*co2t + k19cb(isot)*n2dble &
              + k19cc(isot)*codble
         d19cp4 = k19cap(isot)*co2t + k19cbp(isot)*n2dble &
              + k19ccp(isot)*codble
                                !
         l11 = d19c1 + k20c(1)*o3pdbl
         p11 = ( d19cp1 + k20cp(1)*o3pdbl ) * n10(i)
         l21 = d19c2 + k20c(2)*o3pdbl
         p21 = ( d19cp2 + k20cp(2)*o3pdbl ) *n20(i)
         l31 = d19c3 + k20c(3)*o3pdbl
         p31 = ( d19cp3 + k20cp(3)*o3pdbl ) *n30(i)
         l41 = d19c4 + k20c(4)*o3pdbl
         p41 = ( d19cp4 + k20cp(4)*o3pdbl ) *n40(i)

                                ! Addition of V-V

         l11 = l11 + k21cp(2)*n20(i) + k21cp(3)*n30(i) &
              + k21cp(4)*n40(i)
         p1121 = k21c(2) * n10(i)
         p1131 = k21c(3) * n10(i)
         p1141 = k21c(4) * n10(i)
                                !
         l21 = l21 + k21c(2)*n10(i) + k23k21c*n30(i) + k24k21c*n40(i)
         p2111 = k21cp(2) * n20(i)
         p2131 = k23k21cp * n20(i)
         p2141 = k24k21cp * n20(i)
                                !
         l31 = l31 + k21c(3)*n10(i) + k23k21cp*n20(i) + k34k21c*n40(i)
         p3111 = k21cp(3)* n30(i)
         p3121 = k23k21c * n30(i)
         p3141 = k34k21cp* n30(i)
                                !
         l41 = l41 + k21c(4)*n10(i) + k24k21cp*n20(i) + k34k21cp*n30(i)
         p4111 = k21cp(4)* n40(i)
         p4121 = k24k21c * n40(i)
         p4131 = k34k21c * n40(i)


         if ( input_cza.ge.1 ) then

	    l12 = d19b1 &
                 + k20b(1)*o3pdbl &
                 + k21b(1)*n10(i) &
                 + k33c*( n20(i) + n30(i) + n40(i) )
	    p12 = k21bp(1)*n11(i) * n11(i)
	    p1211 = d19bp1 + k20bp(1)*o3pdbl
	    p1221 = k33cp(2)*n11(i)
	    p1231 = k33cp(3)*n11(i)
	    p1241 = k33cp(4)*n11(i)

	    l11 = l11 + d19bp1 &
                 + k20bp(1)*o3pdbl &
                 + 2.d0 * k21bp(1) * n11(i) &
                 +   k33cp(2)*n21(i) + k33cp(3)*n31(i) + k33cp(4)*n41(i)
	    p1112 = d19b1  &
                 + k20b(1)*o3pdbl &
                 + 2.d0*k21b(1)*n10(i) &
                 + k33c*( n20(i) + n30(i) + n40(i) )

	    l21 = l21 + k33cp(2)*n11(i)
	    p2112 = k33c*n20(i)

	    l31 = l31 + k33cp(3)*n11(i)
	    p3112 = k33c*n30(i)

	    l41 = l41 + k33cp(4)*n11(i)
	    p4112 = k33c*n40(i)

         end if


         ! For ITT=13,15
         a21_einst(i) = a2_010_000 * 1.8d0 / 4.d0 * taustar21(i)
         a31_einst(i) = a3_010_000 * 1.8d0 / 4.d0 * taustar31(i)
         a41_einst(i) = a4_010_000 * 1.8d0 / 4.d0 * taustar41(i)

         l21 = l21 + a21_einst(i)
         l31 = l31 + a31_einst(i)
         l41 = l41 + a41_einst(i)

         ! For ITT=13
         if (input_cza.ge.1 .and. itt_cza.eq.13) then
            a12_einst(i) = a1_020_010/3.d0 * 1.8d0/4.d0 * taustar12(i)
            l12=l12+a12_einst(i)
         endif


         ! Checking for collisional severe errors
         if (l11 .le. 0.0d0) then
            ierr = 21
            varerr = l11
            return
         elseif (l21 .le. 0.0d0) then
            ierr = 22
            varerr = l21
            return
         elseif (l31 .le. 0.0d0) then
            ierr = 23
            varerr = l31
            return
         elseif (l41 .le. 0.0d0) then
            ierr = 24
            varerr = l41
            return
         endif
         if (input_cza.ge.1) then
	    if (l12 .lt. 0.0d0) then
               ierr = 25
               varerr = l12
               return
            endif
         endif
         !

         a11(i) = gamma1*nu11**3.d0 * 1.d0/2.d0 * (p11) /  &
              (n10(i)*l11)
         a1121(i,i) = (nu11/nu21)**3.d0 * n20(i)/n10(i) * p1121/l11
         a1131(i,i) = (nu11/nu31)**3.d0 * n30(i)/n10(i) * p1131/l11
         a1141(i,i) = (nu11/nu41)**3.d0 * n40(i)/n10(i) * p1141/l11
         e110(i) = 2.d0* vlight*nu11**2.d0 * 1.d0/2.d0 / &
              ( n10(i) * l11 )

         a21(i) = gamma1*nu21**3.d0 * 1.d0/2.d0 *  &
              (p21)/(n20(i)*l21)
         a2111(i,i) = (nu21/nu11)**3.d0 * n10(i)/n20(i) * p2111/l21
         a2131(i,i) = (nu21/nu31)**3.d0 * n30(i)/n20(i) * p2131/l21
         a2141(i,i) = (nu21/nu41)**3.d0 * n40(i)/n20(i) * p2141/l21
         e210(i) = 2.d0*vlight*nu21**2.d0 * 1.d0/2.d0 /  &
              ( n20(i) * l21 )

         a31(i) = gamma1*nu31**3.d0 * 1.d0/2.d0 * (p31) /   &
              (n30(i)*l31)
         a3111(i,i) = (nu31/nu11)**3.d0 * n10(i)/n30(i) * p3111/l31
         a3121(i,i) = (nu31/nu21)**3.d0 * n20(i)/n30(i) * p3121/l31
         a3141(i,i) = (nu31/nu41)**3.d0 * n40(i)/n30(i) * p3141/l31
         e310(i) = 2.d0*vlight*nu31**2.d0 * 1.d0/2.d0 /   &
              ( n30(i) * l31 )

         a41(i) = gamma1*nu41**3.d0 * 1.d0/2.d0 * (p41) /   &
              (n40(i)*l41)
         a4111(i,i) = (nu41/nu11)**3.d0 * n10(i)/n40(i) * p4111/l41
         a4121(i,i) = (nu41/nu21)**3.d0 * n20(i)/n40(i) * p4121/l41
         a4131(i,i) = (nu41/nu31)**3.d0 * n30(i)/n40(i) * p4131/l41
         e410(i) = 2.d0*vlight*nu41**2.d0 * 1.d0/2.d0 /   &
              ( n40(i) * l41 )

         if (input_cza.ge.1) then

	    a1112(i,i) = (nu11/nu121)**3.d0 * n11(i)/n10(i) *   &
                 p1112/l11
	    a2112(i,i) = (nu21/nu121)**3.d0 * n11(i)/n20(i) *   &
                p2112/l21
	    a3112(i,i) = (nu31/nu121)**3.d0 * n11(i)/n30(i) *   &
                 p3112/l31
	    a4112(i,i) = (nu41/nu121)**3.d0 * n11(i)/n40(i) *   &
                 p4112/l41
	    a12(i) = gamma1*nu121**3.d0 *2.d0/4.d0* (p12)/   &
                 (n11(i)*l12)
	    a1211(i,i) = (nu121/nu11)**3.d0 * n10(i)/n11(i) *   &
                 p1211/l12
	    a1221(i,i) = (nu121/nu21)**3.d0 * n20(i)/n11(i) *   &
                 p1221/l12
	    a1231(i,i) = (nu121/nu31)**3.d0 * n30(i)/n11(i) *   &
                 p1231/l12
	    a1241(i,i) = (nu121/nu41)**3.d0 * n40(i)/n11(i) *   &
                 p1241/l12
	    e121(i) = 2.d0*vlight*nu121**2.d0 *2.d0/4.d0 /   &
                 ( n11(i) * l12 )

         end if


 4    continue    !-------------------------------------------------------


                  !!!!!!!!!!!! Solucion del sistema

                  !! Paso 0 :  Calculo de los alphas   alf11, alf21, alf31, alf41, alf12

      call unit  ( cax2, nl )

      call diago ( cax1, e110, nl )
      call mulmmf90 ( cax3, cax1,c110, nl )
      call resmmf90 ( alf11, cax2,cax3, nl )

      call diago ( cax1, e210, nl )
      call mulmmf90 ( cax3, cax1,c210, nl )
      call resmmf90 ( alf21, cax2,cax3, nl )

      call diago ( cax1, e310, nl )
      call mulmmf90 ( cax3, cax1,c310, nl )
      call resmmf90 ( alf31, cax2,cax3, nl )

      call diago ( cax1, e410, nl )
      call mulmmf90 ( cax3, cax1,c410, nl )
      call resmmf90 ( alf41, cax2,cax3, nl )

      if (input_cza.ge.1) then
         call diago ( cax1, e121, nl )
         call mulmmf90 ( cax3, cax1,c121, nl )
         call resmmf90 ( alf12, cax2,cax3, nl )
      endif

                                !! Paso 1 :  Calculo de vectores y matrices con 1 barra (aa***)

      if (input_cza.eq.0) then  !  Skip paso 1, pues el12 no se calcula

                                ! el11
         call sypvvv( aa11, a11,e110,sl110, nl )
         call samem( aa1121, a1121, nl )
         call samem( aa1131, a1131, nl )
         call samem( aa1141, a1141, nl )
         call samem( aalf11, alf11, nl )

                                ! el21
         call sypvvv( aa21, a21,e210,sl210, nl )
         call samem( aa2111, a2111, nl )
         call samem( aa2131, a2131, nl )
         call samem( aa2141, a2141, nl )
         call samem( aalf21, alf21, nl )

                                ! el31
         call sypvvv( aa31, a31,e310,sl310, nl )
         call samem( aa3111, a3111, nl )
         call samem( aa3121, a3121, nl )
         call samem( aa3141, a3141, nl )
         call samem( aalf31, alf31, nl )

                                ! el41
         call sypvvv( aa41, a41,e410,sl410, nl )
         call samem( aa4111, a4111, nl )
         call samem( aa4121, a4121, nl )
         call samem( aa4131, a4131, nl )
         call samem( aalf41, alf41, nl )


      else                      !      (input_cza.ge.1) ,   FH !


         call sypvvv( v1, a12,e121,sl121, nl ) ! a12 + e121 * sl121

                                ! aa11
         call sypvvv( v2, a11,e110,sl110, nl )
         call trucommvv( aa11 , alf12,a1112,v2, v1, nl )

                                ! aalf11
         call invdiag( cax1, a1112, nl )
         call mulmmf90( cax2, alf12, cax1, nl ) ! alf12 * (1/a1112)
         call mulmmf90( cax3, cax2, alf11, nl )
         call resmmf90( aalf11, cax3, a1211, nl )
                                ! aa1121
         call trucodiag(aa1121, alf12,a1112,a1121, a1221, nl)
                                ! aa1131
         call trucodiag(aa1131, alf12,a1112,a1131, a1231, nl)
                                ! aa1141
         call trucodiag(aa1141, alf12,a1112,a1141, a1241, nl)


                                ! aa21
         call sypvvv( v2, a21,e210,sl210, nl )
         call trucommvv( aa21 , alf12,a2112,v2, v1, nl )

                                ! aalf21
         call invdiag( cax1, a2112, nl )
         call mulmmf90( cax2, alf12, cax1, nl ) ! alf12 * (1/a2112)
         call mulmmf90( cax3, cax2, alf21, nl )
         call resmmf90( aalf21, cax3, a1221, nl )
                                ! aa2111
         call trucodiag(aa2111, alf12,a2112,a2111, a1211, nl)
                                ! aa2131
         call trucodiag(aa2131, alf12,a2112,a2131, a1231, nl)
                                ! aa2141
         call trucodiag(aa2141, alf12,a2112,a2141, a1241, nl)


                                ! aa31
         call sypvvv ( v2, a31,e310,sl310, nl )
         call trucommvv( aa31 , alf12,a3112,v2, v1, nl )
                                ! aalf31
         call invdiag( cax1, a3112, nl )
         call mulmmf90( cax2, alf12, cax1, nl ) ! alf12 * (1/a3112)
         call mulmmf90( cax3, cax2, alf31, nl )
         call resmmf90( aalf31, cax3, a1231, nl )
                                ! aa3111
         call trucodiag(aa3111, alf12,a3112,a3111, a1211, nl)
                                ! aa3121
         call trucodiag(aa3121, alf12,a3112,a3121, a1221, nl)
                                ! aa3141
         call trucodiag(aa3141, alf12,a3112,a3141, a1241, nl)


                                ! aa41
         call sypvvv( v2, a41,e410,sl410, nl )
         call trucommvv( aa41 , alf12,a4112,v2, v1, nl )
                                ! aalf41
         call invdiag( cax1, a4112, nl )
         call mulmmf90( cax2, alf12, cax1, nl ) ! alf12 * (1/a4112)
         call mulmmf90( cax3, cax2, alf41, nl )
         call resmmf90( aalf41, cax3, a1241, nl )
                                ! aa4111
         call trucodiag(aa4111, alf12,a4112,a4111, a1211, nl)
                                ! aa4121
         call trucodiag(aa4121, alf12,a4112,a4121, a1221, nl)
                                ! aa4131
         call trucodiag(aa4131, alf12,a4112,a4131, a1231, nl)

      endif                     ! Final  caso input_cza.ge.1


                                !! Paso 2 :  Calculo de vectores y matrices con 2 barras (aaa***)

                                ! aaalf41
      call invdiag( cax1, aa4121, nl )
      call mulmmf90( cax2, aalf21, cax1, nl ) ! alf21 * (1/a4121)
      call mulmmf90( cax3, cax2, aalf41, nl )
      call resmmf90( aaalf41, cax3, aa2141, nl )
                                ! aaa41
      call trucommvv(aaa41, aalf21,aa4121,aa41, aa21, nl)
                                ! aaa4111
      call trucodiag(aaa4111, aalf21,aa4121,aa4111, aa2111, nl)
                                ! aaa4131
      call trucodiag(aaa4131, aalf21,aa4121,aa4131, aa2131, nl)

                                ! aaalf31
      call invdiag( cax1, aa3121, nl )
      call mulmmf90( cax2, aalf21, cax1, nl ) ! alf21 * (1/a3121)
      call mulmmf90( cax3, cax2, aalf31, nl )
      call resmmf90( aaalf31, cax3, aa2131, nl )
                                ! aaa31
      call trucommvv(aaa31, aalf21,aa3121,aa31, aa21, nl)
                                ! aaa3111
      call trucodiag(aaa3111, aalf21,aa3121,aa3111, aa2111, nl)
                                ! aaa3141
      call trucodiag(aaa3141, aalf21,aa3121,aa3141, aa2141, nl)

                                ! aaalf11
      call invdiag( cax1, aa1121, nl )
      call mulmmf90( cax2, aalf21, cax1, nl ) ! alf21 * (1/a1121)
      call mulmmf90( cax3, cax2, aalf11, nl )
      call resmmf90( aaalf11, cax3, aa2111, nl )
                                ! aaa11
      call trucommvv(aaa11, aalf21,aa1121,aa11, aa21, nl)
                                ! aaa1131
      call trucodiag(aaa1131, aalf21,aa1121,aa1131, aa2131, nl)
                                ! aaa1141
      call trucodiag(aaa1141, aalf21,aa1121,aa1141, aa2141, nl)


                                !! Paso 3 :  Calculo de vectores y matrices con 3 barras (aaaa***)

                                ! aaaalf41
      call invdiag( cax1, aaa4131, nl )
      call mulmmf90( cax2, aaalf31, cax1, nl ) ! aaalf31 * (1/aaa4131)
      call mulmmf90( cax3, cax2, aaalf41, nl )
      call resmmf90( aaaalf41, cax3, aaa3141, nl )
                                ! aaaa41
      call trucommvv(aaaa41, aaalf31,aaa4131,aaa41, aaa31, nl)
                                ! aaaa4111
      call trucodiag(aaaa4111, aaalf31,aaa4131,aaa4111,aaa3111, nl)

                                ! aaaalf11
      call invdiag( cax1, aaa1131, nl )
      call mulmmf90( cax2, aaalf31, cax1, nl ) ! aaalf31 * (1/aaa4131)
      call mulmmf90( cax3, cax2, aaalf11, nl )
      call resmmf90( aaaalf11, cax3, aaa3111, nl )
                                ! aaaa11
      call trucommvv(aaaa11, aaalf31,aaa1131,aaa11, aaa31, nl)
                                ! aaaa1141
      call trucodiag(aaaa1141, aaalf31,aaa1131,aaa1141,aaa3141, nl)


                                !! Paso 4 :  Calculo de vectores y matrices finales y calculo de J1

      call trucommvv(v1, aaaalf41,aaaa1141,aaaa11, aaaa41, nl)
                                !
      call invdiag( cax1, aaaa1141, nl )
      call mulmmf90( cax2, aaaalf41, cax1, nl ) ! aaaalf41 * (1/aaaa1141)
      call mulmmf90( cax3, cax2, aaaalf11, nl )
      call resmmf90( cax1, cax3, aaaa4111, nl )
                                !
      call LUdec ( el11, cax1, v1, nl, nl2 )

                                ! Solucion para el41
      call sypvmv( v1, aaaa41, aaaa4111,el11, nl )
      call LUdec ( el41, aaaalf41, v1, nl, nl2 )

                                ! Solucion para el31
      call sypvmv( v2, aaa31, aaa3111,el11, nl )
      call sypvmv( v1,    v2, aaa3141,el41, nl )
      call LUdec ( el31, aaalf31, v1, nl, nl2 )

                                ! Solucion para el21
      call sypvmv( v3, aa21, aa2111,el11, nl )
      call sypvmv( v2,   v3, aa2131,el31, nl )
      call sypvmv( v1,   v2, aa2141,el41, nl )
      call LUdec ( el21, aalf21, v1, nl, nl2 )

                                !!!
      el11(1) = planckdp( tl(1), nu11 )
      el21(1) = planckdp( tl(1), nu21 )
      el31(1) = planckdp( tl(1), nu31 )
      el41(1) = planckdp( tl(1), nu41 )
      el11(nl) = 2.d0 * el11(nl-1) - el11(nl2)
      el21(nl) = 2.d0 * el21(nl-1) - el21(nl2)
      el31(nl) = 2.d0 * el31(nl-1) - el31(nl2)
      el41(nl) = 2.d0 * el41(nl-1) - el41(nl2)

      call mulmv ( v1, c110,el11, nl )
      call sumvv ( hr110, v1,sl110, nl )

                                ! Solucion para el12
      if (input_cza.ge.1) then

         call sypvmv( v1, a12, a1211,el11, nl )
         call sypvmv( v3,  v1, a1221,el21, nl )
         call sypvmv( v2,  v3, a1231,el31, nl )
         call sypvmv( v1,  v2, a1241,el41, nl )
         call LUdec ( el12, alf12, v1, nl, nl2 )

         el12(1) = planckdp( tl(1), nu121 )
         el12(nl) = 2.d0 * el12(nl-1) - el12(nl2)

         if (itt_cza.eq.15) then
            call mulmv ( v1, c121,el12, nl )
            call sumvv ( hr121, v1,sl121, nl )
         endif

      end if



      if (input_cza.lt.1) then

         minvt11 = 1.d6
         minvt21 = 1.d6
         minvt31 = 1.d6
         minvt41 = 1.d6
         do i=1,nl
	    pl11 = el11(i)/( gamma1 * nu11**3.0d0  * 1.d0/2.d0 /n10(i) )
	    pl21 = el21(i)/( gamma1 * nu21**3.0d0  * 1.d0/2.d0 /n20(i) )
	    pl31 = el31(i)/( gamma1 * nu31**3.0d0  * 1.d0/2.d0 /n30(i) )
	    pl41 = el41(i)/( gamma1 * nu41**3.0d0  * 1.d0/2.d0 /n40(i) )
	    vt11(i) = -ee*nu11 / log( abs(pl11) / (2.0d0*n10(i)) )
	    vt21(i) = -ee*nu21 / log( abs(pl21) / (2.0d0*n20(i)) )
	    vt31(i) = -ee*nu31 / log( abs(pl31) / (2.0d0*n30(i)) )
	    vt41(i) = -ee*nu41 / log( abs(pl41) / (2.0d0*n40(i)) )
	    hr210(i) = sl210(i) -hplanck*vlight*nu21 *a21_einst(i)*pl21
	    hr310(i) = sl310(i) -hplanck*vlight*nu31 *a31_einst(i)*pl31
	    hr410(i) = sl410(i) -hplanck*vlight*nu41 *a41_einst(i)*pl41

            minvt11 = min( minvt11,vt11(i) )
	    minvt21 = min( minvt21,vt21(i) )
	    minvt31 = min( minvt31,vt31(i) )
	    minvt41 = min( minvt41,vt41(i) )
         enddo

         ! Checking for errors in Tvibs
         if (minvt11 .le. 0.d0) then
            ierr = 26
            varerr = minvt11
            return
         elseif (minvt21 .le. 0.d0) then
            ierr = 27
            varerr = minvt21
            return
         elseif (minvt31 .le. 0.d0) then
            ierr = 28
            varerr = minvt31
            return
         elseif (minvt41 .le. 0.d0) then
            ierr = 29
            varerr = minvt41
            return
         endif

         v626t1(1:nl)=vt11(1:nl)
         v628t1(1:nl)=vt21(1:nl)
         v636t1(1:nl)=vt31(1:nl)
         v627t1(1:nl)=vt41(1:nl)
!         call dinterconnection( v626t1, vt11 )
!         call dinterconnection ( v628t1, vt21 )
!         call dinterconnection ( v636t1, vt31 )
!         call dinterconnection ( v627t1, vt41 )

      else

         do i=1,nl
	    pl21 = el21(i)/( gamma1 * nu21**3.0d0 * 1.d0/2.d0 / n20(i) )
	    pl31 = el31(i)/( gamma1 * nu31**3.0d0 * 1.d0/2.d0 / n30(i) )
	    pl41 = el41(i)/( gamma1 * nu41**3.0d0 * 1.d0/2.d0 / n40(i) )
	    hr210(i) = sl210(i) -hplanck*vlight*nu21 *a21_einst(i)*pl21
	    hr310(i) = sl310(i) -hplanck*vlight*nu31 *a31_einst(i)*pl31
	    hr410(i) = sl410(i) -hplanck*vlight*nu41 *a41_einst(i)*pl41
 	    if (itt_cza.eq.13) then
               pl12 = el12(i)/( gamma1*nu121**3.0d0 * 2.d0/4.d0 /n11(i) )
               hr121(i) = - hplanck*vlight * nu121 * a12_einst(i)*pl12
               hr121(i) = hr121(i) + sl121(i)
            endif
         enddo

      endif

                                ! K/Dday
      do i=1,nl
!        hr110(i)=hr110(i)*dble( hrkday_factor(i) / nt(i) )
!        hr210(i)=hr210(i)*dble( hrkday_factor(i) / nt(i) )
!        hr310(i)=hr310(i)*dble( hrkday_factor(i) / nt(i) )
!        hr410(i)=hr410(i)*dble( hrkday_factor(i) / nt(i) )
!        hr121(i)=hr121(i)*dble( hrkday_factor(i) / nt(i) )
         hr110(i)=hr110(i)*dble( hrkday_factor(i) / ntl(i) )
         hr210(i)=hr210(i)*dble( hrkday_factor(i) / ntl(i) )
         hr310(i)=hr310(i)*dble( hrkday_factor(i) / ntl(i) )
         hr410(i)=hr410(i)*dble( hrkday_factor(i) / ntl(i) )
         hr121(i)=hr121(i)*dble( hrkday_factor(i) / ntl(i) )
      end do

end subroutine NLTEdlvr11_CZALU


!***********************************************************************

subroutine NLTEdlvr11_FB626CTS( hr110CTS, nl_cts_real )

!***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!!!!!!!!!!!!!!!!!! common variables and constants


! Arguments
      real*8 :: hr110CTS(nl_cts)   ! output
      integer  :: nl_cts_real      ! i

! local variables

      real*8 :: n11CTS(nl_cts), slopeTstar110(nl_cts)
      real*8 :: n10(nl_cts), co2t, codbl, n2dbl, o3pdbl
      real*8 :: d19c1, d19cp1, l11, p11
      real*8 :: a11_einst(nl_cts), hcv, maxslope
      integer :: i, isot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  start program

      nu11 = dble(nu(1,1))
      hcv =  hplanck*vlight*nu11

      call zero2v (hr110CTS,n11CTS,nl_cts)

      do i=1,nl_cts_real

         co2t = dble ( co2_cts(i) *(imr(1)+imr(3)+imr(2)+imr(4)) )
         n10(i) = dble( co2_cts(i) * imr(1) )
         codbl = dble(co_cts(i))
         o3pdbl = dble(o3p_cts(i))
         n2dbl = dble(n2_cts(i))

         call GETK_dlvr11 ( t_cts(i) )
         isot = 1
         d19c1 = k19ca(isot)*co2t + k19cb(isot)*n2dbl  &
                   + k19cc(isot)*codbl
         d19cp1 = k19cap(isot)*co2t + k19cbp(isot)*n2dbl  &
                   + k19ccp(isot)*codbl
         l11 = d19c1 + k20c(1)*o3pdbl
         p11 = ( d19cp1 + k20cp(1)*o3pdbl ) * n10(i)

         a11_einst(i) = a1_010_000 * 1.8d0/4.d0 * taustar11_cts(i)

         n11CTS(i) = p11 / (l11 + a11_einst(i))

         hr110CTS(i) = - n11CTS(i) * a11_einst(i) * hcv
         hr110CTS(i) = hr110CTS(i)*  &
              dble( hrkday_factor_cts(i) / nt_cts(i) ) !K/Day

      enddo


! calculo de la altura de transicion, a partir de Tstar
! y merging con el hr110(i), ya calculado con CZALU

      slopeTstar110(1) = taustar11_cts(2)-taustar11_cts(1)
      slopeTstar110(nl_cts_real) = taustar11_cts(nl_cts_real) -  &
                taustar11_cts(nl_cts_real-1)
      maxslope = max( slopeTstar110(1),slopeTstar110(nl_cts_real))
      if (nl_cts_real .gt. 2) then
         do i=2,nl_cts_real-1
            slopeTstar110(i) = ( taustar11_cts(i+1) - &
                      taustar11_cts(i-1) ) * 0.5d0
            if ( slopeTstar110(i) .gt. maxslope ) then
                !write (*,*) i, pl_cts(i), maxslope, slopeTstar110(i)
               maxslope=slopeTstar110(i)
            endif
         enddo
      endif

end subroutine NLTEdlvr11_FB626CTS

function hrkday_convert(mmean_nlte,cpmean_nlte)
!***********************************************************************
!     hrkday_convert.f   ***********paramaters instead from MGITM **********
!
!     fortran function that returns the factor for conversion from
!     hr' [erg s-1 cm-3] to hr [ k day-1 ]
!     Retain units of K/Earth_day for all heat balance terms in MGITM
!     nov 2017        swb      adapted to MGITM
!     mar 2010        fgg      adapted to GCM
!     jan 99          malv     add o2 as major component.
!     ago 98          malv     also returns cp_avg,pm_avg
!     jul 98          malv     first version.
!***********************************************************************

      use ModConstants, only :  Avogadros_Number

      implicit none

!     argumentos
      real, parameter :: n_avog = Avogadros_Number
!     real, parameter :: daysec = 88775.  !  Mars daysec
      real, parameter :: daysec = 86400.
      real  :: mmean_nlte,cpmean_nlte
      real :: hrkday_convert

!cccccccccccccccccccccccccccccccccccc

      hrkday_convert = daysec * n_avog /  &
           ( cpmean_nlte * 1.e4 * mmean_nlte )

end function hrkday_convert

subroutine ERRORS (ierr,varerr)

!***********************************************************************

      implicit none

! Arguments
      integer  :: ierr
      real*8 :: varerr

!***************

      if (ierr .eq. 15) then
         write (*,*) ' ERROR in MZESC110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   c2 < 0 after INTZHUNT_CTS'

      elseif (ierr .eq. 16) then
         write (*,*) ' ERROR in MZESC110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   p2 < 0 after INTZHUNT_CTS'

      elseif (ierr .eq. 17) then
         write (*,*) ' ERROR in MZESC110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   mr2 < 0 after INTZHUNT_CTS'

      elseif (ierr .eq. 18) then
         write (*,*) ' ERROR in MZESC110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   t2 < 0 after INTZHUNT_CTS'

      elseif (ierr .eq. 19) then
         write (*,*) ' ERROR in MZESC110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   st2 < 0 after INTZHUNT_CTS'

      elseif (ierr .eq. 33) then
         write (*,*) ' ERROR in MZESC110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   [CO2] < 0 at TOA.'

      elseif (ierr .eq. 42) then
         write (*,*) ' ERROR in MZTUD110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   Atmospheric transmittance too large. '

      elseif (ierr .eq. 43) then
         write (*,*) ' ERROR in MZTUD110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   [CO2] < 0 at  CurtisMatrix top.'

      elseif (ierr .eq. 45) then
         write (*,*) ' ERROR in MZTUD110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   c2 < 0 after INTZHUNT'

      elseif (ierr .eq. 46) then
         write (*,*) ' ERROR in MZTUD110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   p2 < 0 after INTZHUNT'

      elseif (ierr .eq. 47) then
         write (*,*) ' ERROR in MZTUD110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   mr2 < 0 after INTZHUNT'

      elseif (ierr .eq. 48) then
         write (*,*) ' ERROR in MZTUD110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   t2 < 0 after INTZHUNT'

      elseif (ierr .eq. 49) then
         write (*,*) ' ERROR in MZTUD110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   st2 < 0 after INTZHUNT'

      elseif (ierr .eq. 75) then
         write (*,*) ' ERROR in MZTUD110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   c1 < 0 after INTZHUNT'

      elseif (ierr .eq. 76) then
         write (*,*) ' ERROR in MZTUD110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   p1 < 0 after INTZHUNT'

      elseif (ierr .eq. 77) then
         write (*,*) ' ERROR in MZTUD110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   mr1 < 0 after INTZHUNT'

      elseif (ierr .eq. 78) then
         write (*,*) ' ERROR in MZTUD110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   t1 < 0 after INTZHUNT'

      elseif (ierr .eq. 79) then
         write (*,*) ' ERROR in MZTUD110.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   st1 < 0 after INTZHUNT'

      elseif (ierr .eq. 83) then
         write (*,*) ' ERROR in MZTUD121.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   [CO2] < 0 at  CurtisMatrix top.'

      elseif (ierr .eq. 85) then
         write (*,*) ' ERROR in MZTUD121.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   c1 < 0 after INTZHUNT'

      elseif (ierr .eq. 86) then
         write (*,*) ' ERROR in MZTUD121.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   p1 < 0 after INTZHUNT'

      elseif (ierr .eq. 87) then
         write (*,*) ' ERROR in MZTUD121.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   mr1 < 0 after INTZHUNT'

      elseif (ierr .eq. 88) then
         write (*,*) ' ERROR in MZTUD121.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   t1 < 0 after INTZHUNT'

      elseif (ierr .eq. 89) then
         write (*,*) ' ERROR in MZTUD121.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   st1 < 0 after INTZHUNT'

      elseif (ierr .eq. 51) then
         write (*,*) ' ERROR in MZTVC121.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   Ground transmittance vector VC < 0 '

      elseif (ierr .eq. 52) then
         write (*,*) ' ERROR in MZTVC121.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   Atmospheric transmittance too large. '

      elseif (ierr .eq. 53) then
         write (*,*) ' ERROR in MZTVC121.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   [CO2] < 0 at  CurtisMatrix top.'

      elseif (ierr .eq. 55) then
         write (*,*) ' ERROR in MZTVC121.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   c2 < 0 after INTZHUNT'

      elseif (ierr .eq. 56) then
         write (*,*) ' ERROR in MZTVC121.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   p2 < 0 after INTZHUNT'

      elseif (ierr .eq. 57) then
         write (*,*) ' ERROR in MZTVC121.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   mr2 < 0 after INTZHUNT'

      elseif (ierr .eq. 58) then
         write (*,*) ' ERROR in MZTVC121.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   t2 < 0 after INTZHUNT'

      elseif (ierr .eq. 59) then
         write (*,*) ' ERROR in MZTVC121.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   st2 < 0 after INTZHUNT'

      elseif (ierr .eq. 63) then
         write (*,*) ' ERROR in MZESC121sub.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   [CO2] < 0 at  CurtisMatrix top.'

      elseif (ierr .eq. 65) then
         write (*,*) ' ERROR in MZESC121sub.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   c2 < 0 after INTZHUNT'

      elseif (ierr .eq. 66) then
         write (*,*) ' ERROR in MZESC121sub.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   p2 < 0 after INTZHUNT'

      elseif (ierr .eq. 67) then
         write (*,*) ' ERROR in MZESC121sub.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   mr2 < 0 after INTZHUNT'

      elseif (ierr .eq. 68) then
         write (*,*) ' ERROR in MZESC121sub.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   t2 < 0 after INTZHUNT'

      elseif (ierr .eq. 69) then
         write (*,*) ' ERROR in MZESC121sub.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   st2 < 0 after INTZHUNT'

      elseif (ierr .eq. 21) then
         write (*,*) ' ERROR in CZA.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   l11 < 0 '

      elseif (ierr .eq. 22) then
         write (*,*) ' ERROR in CZA.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   l21 < 0 '

      elseif (ierr .eq. 23) then
         write (*,*) ' ERROR in CZA.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   l31 < 0 '

      elseif (ierr .eq. 24) then
         write (*,*) ' ERROR in CZA.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   l41 < 0 '

      elseif (ierr .eq. 25) then
         write (*,*) ' ERROR in CZA.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   l12 < 0 '

      elseif (ierr .eq. 26) then
         write (*,*) ' ERROR in CZA.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   Negative vibr.temp   xvt11 < 0 '

      elseif (ierr .eq. 27) then
         write (*,*) ' ERROR in CZA.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   Negative vibr.temp   xvt21 < 0 '

      elseif (ierr .eq. 28) then
         write (*,*) ' ERROR in CZA.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   Negative vibr.temp   xvt31 < 0 '

      elseif (ierr .eq. 29) then
         write (*,*) ' ERROR in CZA.    ierr=',ierr
         write (*,*) '                VAR available=', varerr
         write (*,*) '   Negative vibr.temp   xvt41 < 0 '


      endif

      stop ' Stopped in NLTE scheme due to severe error.'
      end

!end subroutine ERRORS

!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fast scheme for NLTE cooling rates at 15um by CO2 in a Martian GCM !
!                 Version dlvr11_03. 2012.                           !
! Software written and provided by IAA/CSIC, Granada, Spain,         !
! under ESA contract "Mars Climate Database and Physical Models"     !
! Person of contact: Miguel Angel Lopez Valverde  valverde@iaa.es    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     **** nlte_calc.F subroutines *****
!
!***********************************************************************
!     subroutine MZESC110 (ig,nl_cts_real, nzy_cts_real,ierr,varerr)
      subroutine MZESC110 (nl_cts_real, nzy_cts_real,ierr,varerr)
!***********************************************************************
!     Includes the following 1-d model subroutines:
!     *** Old MZESC110_dlvr11_03.f
!
!     -MZESC110_dlvr11_03.f
!     -MZTUD110_dlvr11_03.f
!     -MZMC121_dlvr11_03.f
!     -MZTUD121_dlvr11_03.f
!     -MZESC121_dlvr11_03.f
!     -MZESC121sub_dlvr11_03.f
!     -MZTVC121_dlvr11.f
!     -MZTVC121sub_dlvr11_03.f
!***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      integer ::     nl_cts_real, nzy_cts_real ! i
!     integer ::    ig

!     old arguments
      integer ::         ierr      ! o
      real*8 ::          varerr    ! o

!     local variables and constants
      integer :: 	i, iaquiHIST , iaquiZ
      integer ::   isot
      real*8 ::    argumento
      real*8 ::	   tauinf(nl_cts)
      real*8 ::	   con(nzy_cts), coninf
      real*8 ::	   c1, c2 , ccc
      real*8 ::	   t1, t2
      real*8 ::    p1, p2
      real*8 ::	   mr1, mr2
      real*8 ::	   st1, st2
      real*8 ::	   c1box(nbox_max), c2box(nbox_max)
      real*8 ::    ff      ! to avoid too small numbers
      real*8 ::	   st, beta, ts
      real*8 ::	   tyd(nzy_cts)
      real*8 ::	   correc
      real*8 ::	   deltanudbl, deltazdbl
      real*8 ::    yy

!     external function
      external  :: we_clean
      real*8    :: we_clean

!***********************************************************************
      ierr = 0
      varerr = 0.d0
!
      beta = 1.8d5
      ibcode1 = '1'
      isot = 1
      deltanudbl = dble(deltanu(1,1))
      deltazdbl = dble(deltaz_cts)
      ff=1.0d10

!
      do i=1,nzy_cts_real
         tyd(i) = dble(ty_cts(i))
         con(i) =  dble( co2y_cts(i) * imr(isot) )
         correc = 2.d0 * dexp( -ee*dble(elow(isot,2))/tyd(i) )
         con(i) = con(i) * ( 1.d0 - correc )
         mr_cts(i) = dble(co2y_cts(i)/nty_cts(i))
      end do
      if ( con(nzy_cts_real) .le. 0.0d0 ) then
         ierr = 33
         varerr = con(nzy_cts_real)
         return
      elseif ( con(nzy_cts_real-1) .le. con(nzy_cts_real) ) then
         write (*,*) ' WARNING in MZESC110 '
         write (*,*) '    [CO2] growing with altitude at TOA.'
         write (*,*) '    [CO2] @ TOA = ', con(nzy_cts_real)
         coninf = dble( con(nzy_cts_real) )
      else
         coninf = dble( con(nzy_cts_real) /  &
                log( con(nzy_cts_real-1) / con(nzy_cts_real) ) )
      endif
!
      call gethist_03 ( 1 )

!
!     tauinf
!
      call initial

      iaquiHIST = nhist/2
      iaquiZ = nzy_cts_real - 2

      do i=nl_cts_real,1,-1

         if(i.eq.nl_cts_real)then

            call intzhunt_cts (iaquiZ, zl_cts(i), nzy_cts_real,  &
                    c2,p2,mr2,t2, con)
            do kr=1,nbox
               ta(kr)=t2
	    end do
            call interstrhunt (iaquiHIST, st2,t2,ka,ta)
                                ! Check interpolation errors :
            if (c2.le.0.0d0) then
               ierr=15
               varerr=c2
               return
       	    elseif (p2.le.0.0d0) then
               ierr=16
               varerr=p2
               return
	    elseif (mr2.le.0.0d0) then
               ierr=17
               varerr=mr2
               return
	    elseif (t2.le.0.0d0) then
               ierr=18
               varerr=t2
               return
	    elseif (st2.le.0.0d0) then
               ierr=19
               varerr=st2
               return
	    endif
                                !
            aa = p2 * coninf * mr2 * (st2 * ff)
            cc = coninf * st2
            dd = t2 * coninf * st2
            do kr=1,nbox
               ccbox(kr) = coninf * ka(kr)
               ddbox(kr) = t2 * ccbox(kr)
               c2box(kr) = c2 * ka(kr) * deltazdbl
            end do
            c2 = c2 * st2 * deltazdbl

         else

            call intzhunt_cts (iaquiZ, zl_cts(i), nzy_cts_real,  &
                  c1,p1,mr1,t1, con)
            do kr=1,nbox
               ta(kr)=t1
            end do
            call interstrhunt (iaquiHIST, st1,t1,ka,ta)
            do kr=1,nbox
               c1box(kr) = c1 * ka(kr) * deltazdbl
            end do
            c1 = c1 * st1 * deltazdbl
            aa = aa + ( p1*mr1*(c1*ff) + p2*mr2*(c2*ff)) / 2.d0
            cc = cc + ( c1 + c2 ) / 2.d0
            ccc = ( c1 + c2 ) / 2.d0
            dd = dd + ( t1*c1 + t2*c2 ) / 2.d0
            do kr=1,nbox
               ccbox(kr) = ccbox(kr) +  &
                    ( c1box(kr) + c2box(kr) )/2.d0
               ddbox(kr) = ddbox(kr) +  &
                    ( t1*c1box(kr)+t2*c2box(kr) )/2.d0
            end do

            mr2 = mr1
            c2=c1
            do kr=1,nbox
               c2box(kr) = c1box(kr)
            end do
            t2=t1
            p2=p1
         end if

         pp = aa / (cc*ff)

         ts = dd/cc
         do kr=1,nbox
   	    ta(kr) = ddbox(kr) / ccbox(kr)
         end do
         call interstrhunt(iaquiHIST, st,ts,ka,ta)
         call intershphunt(iaquiHIST, alsa,alda,ta)

!
         eqw=0.0d0
         do  kr=1,nbox
            yy = ccbox(kr) * beta
            w = we_clean ( yy, pp, alsa(kr),alda(kr) )
            eqw = eqw + no(kr)*w
         end do

         argumento = eqw / deltanudbl
         tauinf(i) = dexp( - argumento )
         if (i.eq.nl_cts_real) then
            taustar11_cts(i) = 0.0d0
         else
            taustar11_cts(i) = deltanudbl * (tauinf(i+1)-tauinf(i)) &
                 / ( beta * ccc )
         endif

      end do


      call mzescape_normaliz_02 ( taustar11_cts, nl_cts_real, 2 )


      end subroutine MZESC110

!***********************************************************************
!
      subroutine MZTUD110( ierr, varerr )
!
!***********************************************************************
!     *** Old MZTUD110_dlvr11_03.f
!***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      integer  ::         ierr      ! o
      real*8   ::          varerr    ! o

!     local variables and constants
      integer  :: 	      i, in, ir, iaquiHIST , iaquiZ
      integer  ::             ib, isot
      real*8   :: 	      tau(nl,nl), argumento
      real*8   :: 	      tauinf(nl)
      real*8   :: 	      con(nzy), coninf
      real*8   :: 	      c1, c2
      real*8   :: 	      t1, t2
      real*8   :: 	      p1, p2
      real*8   ::	      mr1, mr2
      real*8   :: 	      st1, st2
      real*8   :: 	      c1box(nbox_max), c2box(nbox_max)
      real*8   ::	      ff      ! to avoid too small numbers
      real*8   ::	      tvtbs(nzy)
      real*8   :: 	      st, beta, ts
      real*8   ::  	      zld(nl), zyd(nzy), deltazdbl
      real*8   :: 	      correc
      real*8   ::	      deltanudbl
      real*8   ::             maxtau, yy

!     external function
      external  ::            we_clean
      real*8    ::            we_clean

!***********************************************************************

      ierr = 0
      varerr = 0.d0
!
      ib = 1
      beta = 1.8d5
      ibcode1 = '1'
      isot = 1
      deltanudbl = dble(deltanu(1,1))
      deltazdbl = dble(deltaz)
      ff=1.0d10

!
      do i=1,nzy
         zyd(i) = dble(zy(i))
      enddo
      do i=1,nl
         zld(i) = dble(zl(i))
      enddo
      call interhuntdp ( tvtbs,zyd,nzy, v626t1,zld,nl, 1 )
      do i=1,nzy
         con(i) =  dble( co2y(i) * imr(isot) )
         correc = 2.d0 * dexp( -ee*dble(elow(isot,2))/tvtbs(i) )
         con(i) = con(i) * ( 1.d0 - correc )
         mr(i) = dble(co2y(i)/nty(i))
      end do
      if ( con(nzy) .le. 0.0d0 ) then
         ierr = 43
         varerr = con(nzy)
         return
      elseif ( con(nzy-1) .le. con(nzy) ) then
         write (*,*) ' WARNING in MZTUD110 '
         write (*,*) '    [CO2] grows with height at CurtisMatrix top.'
         write (*,*) '    [CO2] @ top = ', con(nzy)
         coninf = dble( con(nzy) )
      else
         coninf = dble( con(nzy) / log( con(nzy-1) / con(nzy) ) )
      endif
      call mztf_correccion ( coninf, con, ib )

!
      call gethist_03 ( 1 )

!
!     tauinf
!
      call initial

      iaquiHIST = nhist/2
      iaquiZ = nzy - 2

      do i=nl,1,-1

         if(i.eq.nl)then

            call intzhunt (iaquiZ, zl(i),c2,p2,mr2,t2, con)
            do kr=1,nbox
               ta(kr)=t2
            end do
            call interstrhunt (iaquiHIST, st2,t2,ka,ta)
            ! Check interpolation errors :
            if (c2.le.0.0d0) then
               ierr=45
               varerr=c2
               return
            elseif (p2.le.0.0d0) then
               ierr=46
               varerr=p2
               return
            elseif (mr2.le.0.0d0) then
               ierr=47
               varerr=mr2
               return
            elseif (t2.le.0.0d0) then
               ierr=48
               varerr=t2
               return
            elseif (st2.le.0.0d0) then
               ierr=49
               varerr=st2
               return
            endif
                                !
            aa = p2 * coninf * mr2 * (st2 * ff)
            cc = coninf * st2
            dd = t2 * coninf * st2
            do kr=1,nbox
               ccbox(kr) = coninf * ka(kr)
               ddbox(kr) = t2 * ccbox(kr)
               c2box(kr) = c2 * ka(kr) * deltazdbl
            end do
            c2 = c2 * st2 * deltazdbl

         else

            call intzhunt (iaquiZ, zl(i),c1,p1,mr1,t1, con)
            do kr=1,nbox
               ta(kr)=t1
            end do
            call interstrhunt (iaquiHIST, st1,t1,ka,ta)
            do kr=1,nbox
               c1box(kr) = c1 * ka(kr) * deltazdbl
            end do
            ! Check interpolation errors :
            if (c1.le.0.0d0) then
               ierr=75
               varerr=c1
               return
            elseif (p1.le.0.0d0) then
               ierr=76
               varerr=p1
               return
            elseif (mr1.le.0.0d0) then
               ierr=77
               varerr=mr1
               return
            elseif (t1.le.0.0d0) then
               ierr=78
               varerr=t1
               return
            elseif (st1.le.0.0d0) then
               ierr=79
               varerr=st1
               return
            endif
	    !
            c1 = c1 * st1 * deltazdbl
            aa = aa + ( p1*mr1*(c1*ff) + p2*mr2*(c2*ff)) / 2.d0
            cc = cc + ( c1 + c2 ) / 2.d0
            dd = dd + ( t1*c1 + t2*c2 ) / 2.d0
            do kr=1,nbox
               ccbox(kr) = ccbox(kr) +  &
                   ( c1box(kr) + c2box(kr) )/2.d0
               ddbox(kr) = ddbox(kr) +  &
                   ( t1*c1box(kr)+t2*c2box(kr) )/2.d0
            end do

            mr2 = mr1
            c2=c1
            do kr=1,nbox
               c2box(kr) = c1box(kr)
            end do
            t2=t1
            p2=p1
         end if

         pp = aa / (cc*ff)

         ts = dd/cc
         do kr=1,nbox
   	    ta(kr) = ddbox(kr) / ccbox(kr)
         end do
         call interstrhunt(iaquiHIST, st,ts,ka,ta)
         call intershphunt(iaquiHIST, alsa,alda,ta)

!
         eqw=0.0d0
         do  kr=1,nbox
            yy = ccbox(kr) * beta
            w = we_clean ( yy, pp, alsa(kr),alda(kr) )
            eqw = eqw + no(kr)*w
         end do

         argumento = eqw / deltanudbl
         tauinf(i) = dexp( - argumento )

      end do

!
!     tau
!

      iaquiZ = 2
      do 1 in=1,nl-1

         call initial
         call intzhunt (iaquiZ, zl(in), c1,p1,mr1,t1, con)
         do kr=1,nbox
            ta(kr) = t1
         end do
         call interstrhunt (iaquiHIST, st1,t1,ka,ta)
         do kr=1,nbox
            c1box(kr) = c1 * ka(kr) * deltazdbl
         end do
         c1 = c1 * st1 * deltazdbl

         do 2 ir=in,nl-1

            if (ir.eq.in) then
               tau(in,ir) = 1.d0
               goto 2
            end if

            call intzhunt (iaquiZ, zl(ir), c2,p2,mr2,t2, con)
            do kr=1,nbox
               ta(kr) = t2
            end do
            call interstrhunt (iaquiHIST, st2,t2,ka,ta)
            do kr=1,nbox
               c2box(kr) = c2 * ka(kr) * deltazdbl
            end do
            c2 = c2 * st2 * deltazdbl

            aa = aa + ( p1*mr1*(c1*ff) + p2*mr2*(c2*ff)) / 2.d0
            cc = cc + ( c1 + c2 ) / 2.d0
            dd = dd + ( t1*c1 + t2*c2 ) / 2.d0
            do kr=1,nbox
               ccbox(kr) = ccbox(kr) +  &
                     ( c1box(kr) + c2box(kr) ) / 2.d0
               ddbox(kr) = ddbox(kr) +  &
                     ( t1*c1box(kr) + t2*c2box(kr) ) / 2.d0
            end do

            mr1=mr2
            t1=t2
            c1=c2
            p1=p2
            do kr=1,nbox
               c1box(kr) = c2box(kr)
            end do

            pp = aa / (cc * ff)

            ts = dd/cc
            do kr=1,nbox
               ta(kr) = ddbox(kr) / ccbox(kr)
            end do
            call interstrhunt(iaquiHIST, st,ts,ka,ta)
            call intershphunt(iaquiHIST, alsa,alda,ta)
!
            eqw=0.0d0
            do kr=1,nbox
               yy = ccbox(kr) * beta
               w = we_clean ( yy, pp, alsa(kr),alda(kr) )
               eqw = eqw + no(kr)*w
            end do

            argumento = eqw / deltanudbl
            tau(in,ir) = exp( - argumento )


 2       continue

 1    continue


!
!     tau(in,ir) for n>r
!

      in=nl

      call initial

      iaquiZ = nzy - 2
      call intzhunt (iaquiZ, zl(in), c1,p1,mr1,t1, con)
      do kr=1,nbox
         ta(kr) = t1
      end do
      call interstrhunt (iaquiHIST,st1,t1,ka,ta)
      do kr=1,nbox
         c1box(kr) = c1 * ka(kr) * deltazdbl
      end do
      c1 = c1 * st1 * deltazdbl

      do 4 ir=in-1,1,-1

         call intzhunt (iaquiZ, zl(ir), c2,p2,mr2,t2, con)
         do kr=1,nbox
            ta(kr) = t2
         end do
         call interstrhunt (iaquiHIST, st2,t2,ka,ta)
         do kr=1,nbox
            c2box(kr) = c2 * ka(kr) * deltazdbl
         end do
         c2 = c2 * st2 * deltazdbl

         aa = aa + ( p1*mr1*(c1*ff) + p2*mr2*(c2*ff)) / 2.d0
         cc = cc + ( c1 + c2 ) / 2.d0
         dd = dd + ( t1*c1 + t2*c2 ) / 2.d0
         do kr=1,nbox
            ccbox(kr) = ccbox(kr) +  &
                 ( c1box(kr) + c2box(kr) ) / 2.d0
            ddbox(kr) = ddbox(kr) +  &
                 ( t1*c1box(kr) + t2*c2box(kr) ) / 2.d0
         end do

         mr1=mr2
         c1=c2
         t1=t2
         p1=p2
         do kr=1,nbox
            c1box(kr) = c2box(kr)
         end do

         pp = aa / (cc * ff)
         ts = dd / cc
         do kr=1,nbox
            ta(kr) = ddbox(kr) / ccbox(kr)
         end do
         call interstrhunt (iaquiHIST, st,ts,ka,ta)
         call intershphunt (iaquiHIST, alsa,alda,ta)

!

         eqw=0.0d0
         do kr=1,nbox
            yy = ccbox(kr) * beta
            w = we_clean ( yy, pp, alsa(kr),alda(kr) )
            eqw = eqw + no(kr)*w
         end do

         argumento = eqw / deltanudbl
         tau(in,ir) = exp( - argumento )


 4    continue

!
      do in=nl-1,2,-1
         do ir=in-1,1,-1
            tau(in,ir) = tau(ir,in)
         end do
      end do

!
!     Tracking potential numerical errors
!
      maxtau = 0.0d0
      do in=nl-1,2,-1
         do ir=in-1,1,-1
            maxtau = max( maxtau, tau(in,ir) )
         end do
      end do
      if (maxtau .gt. 1.0d0) then
         ierr = 42
         varerr = maxtau
         return
      endif


!
      call MZCUD110 ( tauinf,tau )

      end subroutine MZTUD110


!***********************************************************************

      subroutine MZCUD110 ( tauinf,tau )

!***********************************************************************
!     *** Old file MZCUD_dlvr11.f ***
!***********************************************************************

      use ModConstants, only : pi
      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      real*8  :: 		tau(nl,nl) ! i
      real*8  ::		tauinf(nl) ! i


!     local variables
      integer :: 	i, in, ir
!     real*8  ::	a(nl,nl), cf(nl,nl), pideltanu, deltazdp, pi
      real*8  ::	a(nl,nl), cf(nl,nl), pideltanu, deltazdp

!***********************************************************************

!     pi = 3.141592
      pideltanu = pi * dble(deltanu(1,1))
      deltazdp = 2.0d5 * dble(deltaz)

      do in=1,nl
         do ir=1,nl
            cf(in,ir) = 0.0d0
            c110(in,ir) = 0.0d0
            a(in,ir) = 0.0d0
         end do
         vc110(in) = 0.0d0
      end do

!
      do in=1,nl
         do ir=1,nl

            if (ir.eq.1) then
               cf(in,ir) = tau(in,ir) - tau(in,1)
            elseif (ir.eq.nl) then
               cf(in,ir) = tauinf(in) - tau(in,ir-1)
            else
               cf(in,ir) = tau(in,ir) - tau(in,ir-1)
            end if

         end do
      end do

!
      do in=2,nl-1
         do ir=1,nl
            if (ir.eq.in+1) a(in,ir) = -1.d0
            if (ir.eq.in-1) a(in,ir) = +1.d0
            a(in,ir) = a(in,ir) / deltazdp
         end do
      end do

!
      do in=1,nl
         do ir=1,nl
	    cf(in,ir) = cf(in,ir) * pideltanu
         end do
      end do

      do in=2,nl-1
         do ir=1,nl
	    do i=1,nl
               c110(in,ir) = c110(in,ir) + a(in,i) * cf(i,ir)
	    end do
         end do
      end do

      do in=2,nl-1
         vc110(in) =  pideltanu/deltazdp *  &
              ( tau(in-1,1) - tau(in+1,1) )
      end do


!
      do in=2,nl-1
         c110(in,nl-2) = c110(in,nl-2) - c110(in,nl)
         c110(in,nl-1) = c110(in,nl-1) + 2.d0*c110(in,nl)
      end do

      end subroutine MZCUD110



!***********************************************************************

      subroutine MZMC121

!***********************************************************************
!     *** Old MZMC121_dlvr11_03.f ***
!***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

                                ! local variables

      real*8  ::  cax1(nl,nl)
      real*8  ::  v1(nl), cm_factor, vc_factor
      real    ::  nuaux1, nuaux2, nuaux3
      real*8  ::  faux2,faux3, daux2,daux3
      real*8  ::  varerr

      integer ::  i,j,ik,ib
      integer ::  ierr

!************************************************************************

      c121(1:nl,1:nl)=0.d0
!      call zerom (c121,nl)
      vc121(1:nl)=0.d0
!      call zerov (vc121,nl)

      nuaux1 = nu(1,2) - nu(1,1) ! 667.75
      nuaux2 = nu12_0200-nu(1,1) ! 618.03
      nuaux3 = nu12_1000-nu(1,1) ! 720.81
      faux2 = dble(nuaux2/nuaux1)
      faux3 = dble(nuaux3/nuaux1)
      daux2 = dble(nuaux1-nuaux2)
      daux3 = dble(nuaux1-nuaux3)

      do 11, ik=1,3

         ib=ik+1
         cax1(1:nl,1:nl)=0.d0
!         call zerom (cax1,nl)
         call MZTUD121 ( cax1,v1, ib, ierr, varerr )
         if (ierr .gt. 0) call ERRORS (ierr,varerr)

         do i=1,nl

	    if(ik.eq.1)then
               cm_factor = faux2**2.d0 * exp( daux2*ee/dble(tl(i)) )
               vc_factor = 1.d0/faux2
	    elseif(ik.eq.2)then
               cm_factor = 1.d0
               vc_factor = 1.d0
	    elseif(ik.eq.3)then
               cm_factor = faux3**2.d0 * exp( daux3*ee/dble(tl(i)) )
               vc_factor = 1.d0 / faux3
            else
               write (*,*) ' Error in 626 hot band index  ik =', ik
               stop ' ik can only be = 2,3,4.   Check needed.'
	    end if
	    do j=1,nl
               c121(i,j) = c121(i,j) + cax1(i,j) * cm_factor
	    end do

	    vc121(i) = vc121(i) + v1(i) * vc_factor

         end do

 11   continue

      end subroutine MZMC121


!***********************************************************************

      subroutine MZTUD121 ( cf,vc, ib, ierr, varerr )

!***********************************************************************
!     *** Old MZTUD121_dlvr11_03.f ***
!***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      real*8  ::  	      cf(nl,nl)	! o
      real*8  ::	      vc(nl)    ! o
      integer ::	      ib        ! i
      integer ::              ierr      ! o
      real*8  ::              varerr    ! o


!     local variables and constants
      integer :: 	      i, in, ir, iaquiHIST, iaquiZ
      integer :: 	      isot
      real*8  ::              tau(nl,nl), argumento, deltazdbl
      real*8  :: 	      tauinf(nl)
      real*8  :: 	      con(nzy), coninf
      real*8  :: 	      c1, c2
      real*8  :: 	      t1, t2
      real*8  :: 	      p1, p2
      real*8  ::	      mr1, mr2
      real*8  :: 	      st1, st2
      real*8  :: 	      c1box(nbox_max), c2box(nbox_max)
      real*8  ::	      ff      ! to avoid too small numbers
      real*8  ::	      tvtbs(nzy)
      real*8  :: 	      st, beta, ts
      real*8  ::  	      zld(nl), zyd(nzy)
      real*8  :: 	      correc
      real*8  ::  	      deltanudbl
      real*8  ::              yy

!     external function
      external ::             we_clean
      real*8   ::             we_clean


!     formats
 101  format(i1)
!***********************************************************************

      ierr = 0
      varerr = 0.d0

!     some values
      beta = 1.8d5
      isot = 1
      write (ibcode1,101) ib
      deltanudbl = dble( deltanu(isot,ib) )
      ff=1.0d10
      deltazdbl = dble(deltaz)

!!!
      do i=1,nl
         zld(i) = dble(zl(i))
      enddo
      do i=1,nzy
         zyd(i) = dble(zy(i))
      enddo

      call interhuntdp ( tvtbs,zyd,nzy, v626t1,zld,nl, 1 )

      do i=1,nzy
         con(i) =  dble( co2y(i) * imr(isot) )
         correc = 2.d0 * exp( -ee*dble(elow(isot,2))/tvtbs(i) )
         con(i) = con(i) * ( 1.d0 - correc )
         mr(i) = dble( co2y(i) / nty(i) )
      end do

      if ( con(nzy) .le. 0.0d0 ) then
         ierr = 83
         varerr = con(nzy)
         return
      elseif ( con(nzy-1) .le. con(nzy) ) then
         write (*,*) ' WARNING in MZTUD121 '
         write (*,*) '    [CO2] grows with height at CurtisMatrix top.'
         write (*,*) '    [CO2] @ top = ', con(nzy)
         coninf = dble( con(nzy) )
      else
         coninf = dble( con(nzy) / log( con(nzy-1) / con(nzy) ) )
      endif
      call mztf_correccion ( coninf, con, ib )

!!!
      call gethist_03 ( ib )


!
!     tauinf(nl)
!
      call initial

      iaquiZ = nzy - 2
      iaquiHIST = nhist / 2

      do i=nl,1,-1

         if(i.eq.nl)then

            call intzhunt ( iaquiZ, zl(i),c2,p2,mr2,t2, con)
            do kr=1,nbox
               ta(kr)=t2
            end do
            call interstrhunt (iaquiHIST, st2,t2,ka,ta)
            aa = p2 * coninf * mr2 * (st2 * ff)
            cc = coninf * st2
            dd = t2 * coninf * st2
            do kr=1,nbox
               ccbox(kr) = coninf * ka(kr)
               ddbox(kr) = t2 * ccbox(kr)
               c2box(kr) = c2 * ka(kr) * deltazdbl
            end do
            c2 = c2 * st2 * deltazdbl

         else
            call intzhunt ( iaquiZ, zl(i),c1,p1,mr1,t1, con)
            do kr=1,nbox
               ta(kr)=t1
            end do
            call interstrhunt (iaquiHIST, st1,t1,ka,ta)
            do kr=1,nbox
               c1box(kr) = c1 * ka(kr) * deltazdbl
            end do
            ! Check interpolation errors :
            if (c1.le.0.0d0) then
               ierr=85
               varerr=c1
               return
            elseif (p1.le.0.0d0) then
               ierr=86
               varerr=p1
               return
            elseif (mr1.le.0.0d0) then
               ierr=87
               varerr=mr1
               return
            elseif (t1.le.0.0d0) then
               ierr=88
               varerr=t1
               return
            elseif (st1.le.0.0d0) then
               ierr=89
               varerr=st1
               return
            endif
	    !
            c1 = c1 * st1 * deltazdbl
            aa = aa + ( p1*mr1*(c1*ff) + p2*mr2*(c2*ff)) / 2.d0
            cc = cc + ( c1 + c2 ) / 2.d0
            dd = dd + ( t1*c1 + t2*c2 ) / 2.d0
            do kr=1,nbox
               ccbox(kr) = ccbox(kr) +  &
                   ( c1box(kr) + c2box(kr) )/2.d0
               ddbox(kr) = ddbox(kr) +  &
                   ( t1*c1box(kr)+t2*c2box(kr) )/2.d0
            end do

            mr2 = mr1
            c2=c1
            do kr=1,nbox
               c2box(kr) = c1box(kr)
            end do
            t2=t1
            p2=p1
         end if

         pp = aa / (cc*ff)

         ts = dd/cc
         do kr=1,nbox
   	    ta(kr) = ddbox(kr) / ccbox(kr)
         end do
         call interstrhunt(iaquiHIST, st,ts,ka,ta)
         call intershphunt(iaquiHIST, alsa,alda,ta)

!

         eqw = 0.0d0
         do  kr=1,nbox
            yy = ccbox(kr) * beta
            w = we_clean ( yy, pp, alsa(kr),alda(kr) )
            eqw = eqw + no(kr)*w
         end do

         argumento = eqw / deltanudbl
         tauinf(i) = dexp( - argumento )


      end do                    ! i continue


!
!     tau(in,ir) for n<=r
!

      iaquiZ = 2
      do 1 in=1,nl-1

         call initial
         call intzhunt ( iaquiZ, zl(in), c1,p1,mr1,t1, con)
         do kr=1,nbox
            ta(kr) = t1
         end do
         call interstrhunt (iaquiHIST, st1,t1,ka,ta)
         do kr=1,nbox
            c1box(kr) = c1 * ka(kr) * deltazdbl
         end do
         c1 = c1 * st1 * deltazdbl

         do 2 ir=in,nl-1

            if (ir.eq.in) then
               tau(in,ir) = 1.d0
               goto 2
            end if

            call intzhunt ( iaquiZ, zl(ir), c2,p2,mr2,t2, con)
            do kr=1,nbox
               ta(kr) = t2
            end do
            call interstrhunt (iaquiHIST, st2,t2,ka,ta)
            do kr=1,nbox
               c2box(kr) = c2 * ka(kr) * deltazdbl
            end do
            c2 = c2 * st2 * deltazdbl

            aa = aa + ( p1*mr1*(c1*ff) + p2*mr2*(c2*ff)) / 2.d0
            cc = cc + ( c1 + c2 ) / 2.d0
            dd = dd + ( t1*c1 + t2*c2 ) / 2.d0
            do kr=1,nbox
               ccbox(kr) = ccbox(kr) +  &
                   ( c1box(kr) + c2box(kr) ) / 2.d0
               ddbox(kr) = ddbox(kr) +  &
                   ( t1*c1box(kr) + t2*c2box(kr) ) / 2.d0
            end do

            mr1=mr2
            t1=t2
            c1=c2
            p1=p2
            do kr=1,nbox
               c1box(kr) = c2box(kr)
            end do

            pp = aa / (cc * ff)

            ts = dd/cc
            do kr=1,nbox
               ta(kr) = ddbox(kr) / ccbox(kr)
            end do
            call interstrhunt(iaquiHIST, st,ts,ka,ta)
            call intershphunt(iaquiHIST, alsa,alda,ta)

!

            eqw = 0.0d0
            do kr=1,nbox
               yy = ccbox(kr) * beta
               w = we_clean ( yy, pp, alsa(kr),alda(kr) )
               eqw = eqw + no(kr)*w
            end do

            argumento = eqw / deltanudbl
            tau(in,ir) = dexp( - argumento )

 2       continue

 1    continue

!
!     tau(in,ir) for n>r
!

      in=nl

      call initial
      iaquiZ = nzy - 2
      call intzhunt ( iaquiZ, zl(in), c1,p1,mr1,t1, con)
      do kr=1,nbox
         ta(kr) = t1
      end do
      call interstrhunt (iaquiHIST, st1,t1,ka,ta)
      do kr=1,nbox
         c1box(kr) = c1 * ka(kr) * deltazdbl
      end do
      c1 = c1 * st1 * deltazdbl

      do 4 ir=in-1,1,-1

         call intzhunt ( iaquiZ, zl(ir), c2,p2,mr2,t2, con)
         do kr=1,nbox
            ta(kr) = t2
         end do
         call interstrhunt (iaquiHIST, st2,t2,ka,ta)
         do kr=1,nbox
            c2box(kr) = c2 * ka(kr) * deltazdbl
         end do
         c2 = c2 * st2 * deltazdbl

         aa = aa + ( p1*mr1*(c1*ff) + p2*mr2*(c2*ff)) / 2.d0
         cc = cc + ( c1 + c2 ) / 2.d0
         dd = dd + ( t1*c1 + t2*c2 ) / 2.d0
         do kr=1,nbox
            ccbox(kr) = ccbox(kr) +  &
                ( c1box(kr) + c2box(kr) ) / 2.d0
            ddbox(kr) = ddbox(kr) +  &
                ( t1*c1box(kr) + t2*c2box(kr) ) / 2.d0
         end do

         mr1=mr2
         c1=c2
         t1=t2
         p1=p2
         do kr=1,nbox
            c1box(kr) = c2box(kr)
         end do

         pp = aa / (cc * ff)
         ts = dd / cc
         do kr=1,nbox
            ta(kr) = ddbox(kr) / ccbox(kr)
         end do
         call interstrhunt (iaquiHIST, st,ts,ka,ta)
         call intershphunt (iaquiHIST, alsa,alda,ta)

!
         eqw=0.0d0
         do kr=1,nbox
            yy = ccbox(kr) * beta
            w = we_clean ( yy, pp, alsa(kr),alda(kr) )
            eqw = eqw + no(kr)*w
         end do

         argumento = eqw / deltanudbl
         tau(in,ir) = dexp( - argumento )

 4    continue

!!!

      do in=nl-1,2,-1
         do ir=in-1,1,-1
            tau(in,ir) = tau(ir,in)
         end do
      end do

!
      call MZCUD121 ( tauinf,tau, cf, vc, ib )

      end subroutine MZTUD121



!***********************************************************************

      subroutine MZCUD121 ( tauinf,tau, c,vc, ib )

!***********************************************************************
!     *** Old MZCUD121_dlvr11.f ***
!***********************************************************************

      use ModConstants, only : pi
      use ModPlanet

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      real*8 ::	   c(nl,nl) ! o
      real*8 ::	   vc(nl)  ! o
      real*8 ::	   tau(nl,nl) ! i
      real*8 ::	   tauinf(nl) ! i
      integer ::   ib      ! i

!     local variables
      integer ::   i, in, ir, isot
!     real*8  ::   a(nl,nl), cf(nl,nl), pideltanu, deltazdbl,pi
      real*8  ::   a(nl,nl), cf(nl,nl), pideltanu, deltazdbl

!***********************************************************************

!     pi=3.141592
      isot = 1
      pideltanu = pi*dble(deltanu(isot,ib))
      deltazdbl = dble(deltaz)
!
      do in=1,nl

         do ir=1,nl

            cf(in,ir) = 0.0d0
            c(in,ir) = 0.0d0
            a(in,ir) = 0.0d0

         end do

         vc(in) = 0.0d0

      end do


!
      do in=1,nl
         do ir=1,nl

            if (ir.eq.1) then
               cf(in,ir) = tau(in,ir) - tau(in,1)
            elseif (ir.eq.nl) then
               cf(in,ir) = tauinf(in) - tau(in,ir-1)
            else
               cf(in,ir) = tau(in,ir) - tau(in,ir-1)
            end if

         end do
      end do


!
      do in=2,nl-1
         do ir=1,nl
            if (ir.eq.in+1) a(in,ir) = -1.d0
            if (ir.eq.in-1) a(in,ir) = +1.d0
            a(in,ir) = a(in,ir) / ( 2.d5*deltazdbl )
         end do
      end do

!
      do in=1,nl
         do ir=1,nl
	    cf(in,ir) = cf(in,ir) * pideltanu
         end do
      end do


      do in=2,nl-1
         do ir=1,nl
	    do i=1,nl
               c(in,ir) = c(in,ir) + a(in,i) * cf(i,ir)
	    end do
         end do
         vc(in) =  pideltanu /( 2.d5*deltazdbl ) *  &
                  ( tau(in-1,1) - tau(in+1,1) )
      end do

!
      do in=2,nl-1
         c(in,nl-2) = c(in,nl-2) - c(in,nl)
         c(in,nl-1) = c(in,nl-1) + 2.d0*c(in,nl)
      end do

      end subroutine MZCUD121


!***********************************************************************

      subroutine MZESC121

!***********************************************************************
!     *** Old MZESC121_dlvr11_03.f ***
!***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     local variables
      integer ::          i,ierr
      real*8  ::          factor0200, factor0220, factor1000
      real*8  ::          aux_0200(nl), aux2_0200(nl)
      real*8  ::          aux_0220(nl), aux2_0220(nl)
      real*8  ::          aux_1000(nl), aux2_1000(nl)
      real*8  ::          varerr

!***********************************************************************

!      call zerov (taustar12, nl)
      taustar12(1:nl)=0.d0
      call zero2v(aux_0200,aux2_0200, nl)
      call zero2v(aux_0220,aux2_0220, nl)
      call zero2v(aux_1000,aux2_1000, nl)

      call MZESC121sub (aux_0200,aux2_0200, 2 , ierr, varerr)
      if (ierr .gt. 0) call ERRORS (ierr,varerr)
      call MZESC121sub (aux_0220,aux2_0220, 3 , ierr, varerr)
      if (ierr .gt. 0) call ERRORS (ierr,varerr)
      call MZESC121sub (aux_1000,aux2_1000, 4 , ierr, varerr)
      if (ierr .gt. 0) call ERRORS (ierr,varerr)

      factor0220 = 1.d0
      factor0200 = dble( (nu(1,2)-nu(1,1)) / (nu12_0200-nu(1,1)) )
      factor1000 = dble( (nu(1,2)-nu(1,1)) / (nu12_1000-nu(1,1)) )
      do i=1,nl
         taustar12(i) = taustar12(i)  &
                   + aux_0200(i) * factor0200  &
                   + aux_0220(i) * factor0220  &
                   + aux_1000(i) * factor1000
      enddo

      call mzescape_normaliz ( taustar12, 2 )

      end subroutine



!***********************************************************************

      subroutine MZESC121sub (taustar,tauinf, ib, ierr, varerr )

!***********************************************************************
!     *** Old MZESC121sub_dlvr11_03.f ***
!***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      real*8 ::               taustar(nl) ! o
      real*8 :: 	      tauinf(nl)  ! o
      integer ::	      ib          ! i
      integer ::              ierr        ! o
      real*8  ::              varerr      ! o


!     local variables and constants
      integer :: 	      i, iaquiHIST, iaquiZ, isot
      real*8  :: 	      con(nzy), coninf
      real*8  :: 	      c1, c2, ccc
      real*8  :: 	      t1, t2
      real*8  :: 	      p1, p2
      real*8  ::	      mr1, mr2
      real*8  :: 	      st1, st2
      real*8  :: 	      c1box(70), c2box(70)
      real*8  ::	      ff      ! to avoid too small numbers
      real*8  ::	      tvtbs(nzy)
      real*8  :: 	      st, beta, ts
      real*8  ::  	      zld(nl), zyd(nzy)
      real*8  :: 	      correc
      real*8  :: 	      deltanudbl, deltazdbl
      real*8  ::              yy

!     external function
      external ::             we_clean
      real*8   ::             we_clean

!     formats
 101  format(i1)

!***********************************************************************

      ierr = 0
      varerr = 0.d0
!
      beta = 1.8d5
      isot = 1
      write ( ibcode1, 101) ib
      deltanudbl = dble( deltanu(isot,ib) )
      ff=1.0d10
      deltazdbl = dble(deltaz)

!
      do i=1,nzy
         zyd(i) = dble(zy(i))
      enddo
      do i=1,nl
         zld(i) = dble(zl(i))
      enddo

      call interhuntdp ( tvtbs,zyd,nzy, v626t1,zld,nl, 1 )

      do i=1,nzy
         con(i) =  dble( co2y(i) * imr(isot) )
         correc = 2.d0 * dexp( -ee*dble(elow(isot,2))/tvtbs(i) )
         con(i) = con(i) * ( 1.d0 - correc )
         mr(i) = dble(co2y(i)/nty(i))
      end do
      if ( con(nzy) .le. 0.0d0 ) then
         ierr = 63
         varerr = con(nzy)
         return
      elseif ( con(nzy-1) .le. con(nzy) ) then
         write (*,*) ' WARNING in MZESC121sub '
         write (*,*) '    [CO2] grows with height at CurtisMatrix top.'
         write (*,*) '    [CO2] @ top = ', con(nzy)
         coninf = dble( con(nzy) )
      else
         coninf = dble( con(nzy) / log( con(nzy-1) / con(nzy) ) )
      endif
      call mztf_correccion ( coninf, con, ib )

!
      call gethist_03 ( ib )

!
!     tauinf
!
      call initial

      iaquiHIST = nhist/2
      iaquiZ = nzy - 2

      do i=nl,1,-1

         if(i.eq.nl)then

            call intzhunt (iaquiZ, zl(i),c2,p2,mr2,t2, con)
            do kr=1,nbox
               ta(kr)=t2
	    end do
            call interstrhunt (iaquiHIST, st2,t2,ka,ta)
            ! Check interpolation errors :
            if (c2.le.0.0d0) then
               ierr=65
               varerr=c2
               return
            elseif (p2.le.0.0d0) then
               ierr=66
               varerr=p2
               return
            elseif (mr2.le.0.0d0) then
               ierr=67
               varerr=mr2
               return
            elseif (t2.le.0.0d0) then
               ierr=68
               varerr=t2
               return
            elseif (st2.le.0.0d0) then
               ierr=69
               varerr=st2
               return
            endif
	    !
            aa = p2 * coninf * mr2 * (st2 * ff)
            cc = coninf * st2
            dd = t2 * coninf * st2
            do kr=1,nbox
               ccbox(kr) = coninf * ka(kr)
               ddbox(kr) = t2 * ccbox(kr)
               c2box(kr) = c2 * ka(kr) * deltazdbl
            end do
            c2 = c2 * st2 * deltazdbl

         else
            call intzhunt (iaquiZ, zl(i),c1,p1,mr1,t1, con)
            do kr=1,nbox
               ta(kr)=t1
            end do
            call interstrhunt (iaquiHIST,st1,t1,ka,ta)
            do kr=1,nbox
               c1box(kr) = c1 * ka(kr) * deltazdbl
            end do
            c1 = c1 * st1 * deltazdbl
            aa = aa + ( p1*mr1*(c1*ff) + p2*mr2*(c2*ff)) / 2.d0
            cc = cc + ( c1 + c2 ) / 2.d0
            ccc = ( c1 + c2 ) / 2.d0
            dd = dd + ( t1*c1 + t2*c2 ) / 2.d0
            do kr=1,nbox
               ccbox(kr) = ccbox(kr) +  &
                   ( c1box(kr) + c2box(kr) )/2.d0
               ddbox(kr) = ddbox(kr) +  &
                   ( t1*c1box(kr)+t2*c2box(kr) )/2.d0
            end do

            mr2 = mr1
            c2=c1
            do kr=1,nbox
               c2box(kr) = c1box(kr)
            end do
            t2=t1
            p2=p1
         end if

         pp = aa / (cc*ff)

         ts = dd/cc
         do kr=1,nbox
   	    ta(kr) = ddbox(kr) / ccbox(kr)
         end do
         call interstrhunt(iaquiHIST,st,ts,ka,ta)
         call intershphunt(iaquiHIST,alsa,alda,ta)

!
         eqw=0.0d0
         do  kr=1,nbox
            yy = ccbox(kr) * beta
            w = we_clean ( yy, pp, alsa(kr),alda(kr) )
            eqw = eqw + no(kr)*w
         end do
         tauinf(i) = dexp( - eqw / deltanudbl )
         if (tauinf(i).lt.0.d0) tauinf(i) = 0.0d0

         if (i.eq.nl) then
            taustar(i) = 0.0d0
         else
            taustar(i) = deltanudbl * (tauinf(i+1)-tauinf(i))  &
                 / ( beta * ccc  )
         endif

      end do


      end subroutine MZESC121sub


!***********************************************************************

      subroutine MZTVC121

!***********************************************************************
!     *** Old MZTVC121_dlvr11.f ***
!***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!!!!!!!!!!!!!!!!!!!!!!!
!     common variables & constants


      integer :: ierr
      real*8 :: varerr


!     local variables

      real*8 :: v1(nl), vc_factor
      integer :: i,ik,ib

!************************************************************************

!      call zerov( vc121, nl )
      vc121(1:nl)=0.d0

      do 11, ik=1,3

         ib=ik+1

         call MZTVC121sub (v1, ib, ierr,varerr )

         do i=1,nl

	    if(ik.eq.1)then
               vc_factor =  &
                   dble( (nu(1,2)-nu(1,1)) / (nu12_0200-nu(1,1)) )
	    elseif(ik.eq.2)then
               vc_factor = 1.d0
	    elseif(ik.eq.3)then
               vc_factor =  &
                  dble( (nu(1,2)-nu(1,1)) / (nu12_1000-nu(1,1)) )
	    end if

	    vc121(i) = vc121(i) + v1(i) * vc_factor

         end do

 11   continue

      end subroutine MZTVC121


!***********************************************************************

      subroutine MZTVC121sub  ( vc, ib,  ierr, varerr )

!     *** Old MZTVC121sub_dlvr11_03.f ***
!***********************************************************************

      use ModConstants, only : pi
      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none


!     arguments
      real*8 ::	            vc(nl)  ! o
      integer ::	    ib      ! i
      integer ::            ierr    ! o
      real*8 ::             varerr  ! o

!     local variables and constants
      integer :: 	    i, in, ir, iaquiHIST , iaquiZ, isot
      real*8  :: 	    tau(nl,nl), argumento
      real*8  :: 	    con(nzy), coninf
      real*8  :: 	    c1, c2
      real*8  :: 	    t1, t2
      real*8  :: 	    p1, p2
      real*8  ::	    mr1, mr2
      real*8  :: 	    st1, st2
      real*8  :: 	    c1box(70), c2box(70)
      real*8  ::	    ff      ! to avoid too small numbers
      real*8  ::	    tvtbs(nzy)
      real*8  :: 	    st, beta, ts
      real*8  ::  	    zld(nl), zyd(nzy), deltazdbl
      real*8  :: 	    correc
      real*8  :: 	    deltanudbl, pideltanu
      real*8  ::            yy
      real*8  ::            minvc, maxtau

!     external function
      external ::          we_clean
      real*8   ::          we_clean

!     formats
 101  format(i1)

!***********************************************************************

      ierr = 0
      varerr = 0.d0
!
!     pi=3.141592
      isot = 1
      beta = 1.8d5
      write (ibcode1,101) ib
      deltanudbl = dble( deltanu(isot,ib) )
      pideltanu = pi*deltanudbl
      ff=1.0d10
      deltazdbl = dble(deltaz)
!
!
!

      do i=1,nzy
         zyd(i) = dble(zy(i))
      enddo
      do i=1,nl
         zld(i) = dble(zl(i))
      enddo

      call interhuntdp ( tvtbs,zyd,nzy, v626t1,zld,nl, 1 )

      do i=1,nzy
         con(i) =  dble( co2y(i) * imr(isot) )
         correc = 2.d0 * dexp( -ee*dble(elow(isot,2))/tvtbs(i) )
         con(i) = con(i) * ( 1.d0 - correc )
         mr(i) = dble(co2y(i)/nty(i))
      end do

      if ( con(nzy) .le. 0.0d0 ) then
         ierr = 53
         varerr = con(nzy)
         return
      elseif ( con(nzy-1) .le. con(nzy) ) then
         write (*,*) ' WARNING in MZTVC121sub '
         write (*,*) '    [CO2] grows with height at CurtisMatrix top.'
         write (*,*) '    [CO2] @ top = ', con(nzy)
         coninf = dble( con(nzy) )
      else
         coninf = dble( con(nzy) / log( con(nzy-1) / con(nzy) ) )
      endif
      call mztf_correccion ( coninf, con, ib )

!!!
      call gethist_03 ( ib )

!
!     tau(1,ir)
!
      call initial

      iaquiHIST = nhist/2

      in=1

      tau(in,1) = 1.d0

      call initial
      iaquiZ = 2
      call intzhunt ( iaquiZ, zl(in), c1,p1,mr1,t1, con)
      do kr=1,nbox
         ta(kr) = t1
      end do
      call interstrhunt (iaquiHIST, st1,t1,ka,ta)
      do kr=1,nbox
         c1box(kr) = c1 * ka(kr) * deltazdbl
      end do
      c1 = c1 * st1 * deltazdbl
                                ! Check interpolation errors :
      if (c1.le.0.0d0) then
         ierr=55
         varerr=c1
         return
      elseif (p1.le.0.0d0) then
         ierr=56
         varerr=p1
         return
      elseif (mr1.le.0.0d0) then
         ierr=57
         varerr=mr1
         return
      elseif (t1.le.0.0d0) then
         ierr=58
         varerr=t1
         return
      elseif (st1.le.0.0d0) then
         ierr=59
         varerr=st1
         return
      endif
                                !

      do 2 ir=2,nl

         call intzhunt (iaquiZ, zl(ir), c2,p2,mr2,t2, con)
         do kr=1,nbox
            ta(kr) = t2
         end do
         call interstrhunt (iaquiHIST, st2,t2,ka,ta)
         do kr=1,nbox
            c2box(kr) = c2 * ka(kr) * deltazdbl
         end do
         c2 = c2 * st2 * deltazdbl

         aa = aa + ( p1*mr1*(c1*ff) + p2*mr2*(c2*ff)) / 2.d0
         cc = cc + ( c1 + c2 ) / 2.d0
         dd = dd + ( t1*c1 + t2*c2 ) / 2.d0
         do kr=1,nbox
            ccbox(kr) = ccbox(kr) + ( c1box(kr) + c2box(kr) ) /2.d0
            ddbox(kr) = ddbox(kr) +  &
                  ( t1*c1box(kr) + t2*c2box(kr) ) / 2.d0
         end do

         mr1=mr2
         t1=t2
         c1=c2
         p1=p2
         do kr=1,nbox
            c1box(kr) = c2box(kr)
         end do

         pp = aa / (cc * ff)

         ts = dd/cc
         do kr=1,nbox
   	    ta(kr) = ddbox(kr) / ccbox(kr)
         end do
         call interstrhunt(iaquiHIST, st,ts,ka,ta)
         call intershphunt(iaquiHIST, alsa,alda,ta)

         eqw=0.0d0
         do kr=1,nbox
            yy = ccbox(kr) * beta
            w = we_clean ( yy, pp, alsa(kr),alda(kr) )
            eqw = eqw + no(kr)*w
         end do

         argumento = eqw / deltanudbl
         tau(in,ir) = dexp( - argumento )

 2    continue


!
!
!
      do in=nl,2,-1
         tau(in,1) = tau(1,in)
      end do

!
      vc(1) = 0.0d0
      vc(nl) = 0.0d0
      do in=2,nl-1
         vc(in) =  pideltanu /( 2.d5*deltazdbl ) *  &
                   ( tau(in-1,1) - tau(in+1,1) )
         if (vc(in) .lt. 0.0d0) vc(in) = vc(in-1)
      end do

!
!     Tracking potential numerical errors
!
      minvc = 1.d6
      maxtau = tau(nl,1)
      do in=2,nl-1
         minvc = min( minvc, vc(in) )
         maxtau = max( maxtau, tau(in,1) )
      end do
      if (maxtau .gt. 1.0d0) then
         ierr = 52
         varerr = maxtau
         return
      else if (minvc .lt. 0.0d0) then
         ierr = 51
         varerr = minvc
         return
      endif

      end subroutine MZTVC121sub
! ***********************************************************************

      subroutine init_nlte_setup

! ***********************************************************************
!     swb     Nov 17          Adapt to MGITM on Pleiades
!     malv    Oct 09          Adapt mz1d_onlyTCR_MUCHASveces.f to "V09"
!     malv    Sep 07          Add LU deccomp & repetition option to test CPU
!     malv    Jan 07          Add new vertical fine-grid for NLTE
!     apr 06  malv            Read date,effuv from Driver. T fixed at zbott.
!     2003    fgg             Double precission in UV, Photoq, Conduct & Diff
!     oct 02  malv            V02: New scheme to allow for continuity eq.
!     dec 01  malv            See changes/progress of the code in mz1d.actual
!     nov 01  malv            adapt for parameterizations of tcr y shr
!     nov 98  malv            add chemical & photochem. processes
!     jul 98  malv            transic hiperb con zs fuera de la region
!     equil hidrostatico. smoothing en cr y sh
!     jan 98	malv		first version
!***********************************************************************
      use ModPlanet

!     include	'nlte_paramdef.h'  !  Set MGITM path to all datafiles :
!                                     /UA/DataIn/*.dat
!     include	'nlte_commons.h'

      implicit none

!***********************************************************************

!     local variables

      integer ::i, k, lun1, lun2
      real*8  ::  xx
      character	:: isotcode*2

!     formats
 132  format (i2)

!**********************************************************************

!     *** Groups old 1-d model subroutines SETTINGS and LeeESCTVCISO_dlvr11
!     *** Both were called in old NLTEdlvr11_SETUP ***

!     *** Old SETTINGS ***

      lun1 = 1
      lun2 = 2

      do k=1,nisot
         write (isotcode,132) indexisot(k)
!        open (lun1,  &
!                file=trim(datafile)//'UA/DataIn/enelow' &
!                //isotcode//'.dat',status='old')
         open (lun1, &
                 file='UA/DataIn/enelow' &
                 //isotcode//'.dat',status='old')
!        open (lun2,  &
!                file=trim(datafile)//'UA/DataIn/deltanu' &
!                //isotcode//'.dat',status='old')
         open (lun2,  &
                 file='UA/DataIn/deltanu'  &
                 //isotcode//'.dat',status='old')
         read (lun1,*)
         read (lun2,*)
         read (lun1,*) (elow(k,i), i=1,nb)
         read (lun2,*) (deltanu(k,i), i=1,nb)
         close (lun1)
         close (lun2)
      end do

      a1_010_000 = 1.3546d00
      a2_010_000 = 1.3452d00
      a3_010_000 = 1.1878d00
      a4_010_000 = 1.2455d00
      a1_020_010 = 4.35d0

!     *** Old LeeESCTVCISO_dlvr11 ***

!     open( 11, file=trim(datafile)//  &
!            'UA/DataIn/parametp_Tstar_IAA1204.dat' )
      open( 11, file=  &
             'UA/DataIn/parametp_Tstar_IAA1204.dat' )
      read (11, *)
      do i=1,nztabul
         read (11,*) lnpnbtab(i), tstar11tab(i),  &
                tstar21tab(i), tstar31tab(i), tstar41tab(i)
      enddo
      close (11)

!     open( 12, file=trim(datafile)//  &
!            'UA/DataIn/parametp_VC_IAA1204.dat' )
      open( 12, file=  &
             'UA/DataIn/parametp_VC_IAA1204.dat' )
      read (12, *)
      do i=1,nztabul
         read (12,*) xx, vc210tab(i), vc310tab(i), vc410tab(i)
      enddo
      close (12)
      xx=xx

      call LeeHISTOGRMS

      end subroutine init_nlte_setup

!***********************************************************************
      subroutine LeeHISTOGRMS
!***********************************************************************
      use ModPlanet

!     include	'nlte_paramdef.h'
!     include	'nlte_commons.h'

      implicit none

!     local variables and constants
      integer :: ihist


!***********************************************************************

                                ! Banda fundamental
                                !
!     hisfile = trim(datafile)//  &
!            'UA/DataIn/hid26-1.dat'
      hisfile = 'UA/DataIn/hid26-1.dat'
      ihist = 1
      call rhist_03 (ihist)

                                ! First Hot bands
                                !
!     hisfile = trim(datafile)//  &
!            'UA/DataIn/hid26-2.dat'
      hisfile = 'UA/DataIn/hid26-2.dat'
      ihist = 2
      call rhist_03 (ihist)

!     hisfile = trim(datafile)//  &
!            'UA/DataIn/hid26-3.dat'
      hisfile = 'UA/DataIn/hid26-3.dat'
      ihist = 3
      call rhist_03 (ihist)

!     hisfile = trim(datafile)//  &
!            'UA/DataIn/hid26-4.dat'
      hisfile = 'UA/DataIn/hid26-4.dat'
      ihist = 4
      call rhist_03 (ihist)

      end subroutine LeeHISTOGRMS


!***********************************************************************
!     *** Old GETK_dlvr11.f ***

      subroutine GETK_dlvr11(tt)

!***********************************************************************
      use ModPlanet

!     include	'nlte_paramdef.h'
!     include	'nlte_commons.h'

      implicit none

!     arguments
      real ::	tt	! i. temperature

!     ! local variables:
      real*8 :: k20x, k20xb, k20xc
      real*8 :: k19xca,k19xcb,k19xcc
      real*8 :: k19xba,k19xbb,k19xbc
      real*8 :: k21x,k21xa,k21xb,k21xc
      real*8 :: anu, factor , tdt
      integer :: 	i

!***********************************************************************

      tdt = dble(tt)

                                !! k19 & k20

      k20x = 3.d-12
      k20xc = k20x * rf20
      k20xb = 2.d0 * k20xc

      k19xca = 4.2d-12 * exp( -2988.d0/tdt + 303930.d0/tdt**2.d0 )
      if (tt.le.175.) k19xca = 3.3d-15
      k19xcb = 2.1d-12 * exp( -2659.d0/tdt + 223052.d0/tdt**2.d0 )
      if (tt.le.175.) k19xcb = 7.6d-16
      k19xca = k19xca * rf19
      k19xcb = k19xcb * rf19
      k19xcc = k19xcb

      factor = 2.5d0
      k19xba = factor * k19xca
      k19xbb = factor * k19xcb
      k19xbc = factor * k19xcc

      do i = 1, nisot

         k19ba(i) = k19xba
         k19ca(i) = k19xca
         k19bb(i) = k19xbb
         k19cb(i) = k19xcb
         k19bc(i) = k19xbc
         k19cc(i) = k19xcc

         k20b(i) = k20xb
         k20c(i) = k20xc

         anu = dble( nu(i,2)-nu(i,1) )

         k19bap(i) = k19ba(i) * 2.d0 * exp( -ee*anu/tdt )
         k19bbp(i) = k19bb(i) * 2.d0 * exp( -ee*anu/tdt )
         k19bcp(i) = k19bc(i) * 2.d0 * exp( -ee*anu/tdt )

         k20bp(i) = k20b(i)*4.d0/2.d0 * exp( -ee/tdt * anu )

         anu = dble( nu(i,1) )

         k19cap(i) = k19ca(i) * 2.d0 * exp( -ee*anu/tdt )
         k19cbp(i) = k19cb(i) * 2.d0 * exp( -ee*anu/tdt )
         k19ccp(i) = k19cc(i) * 2.d0 * exp( -ee*anu/tdt )

         k20cp(i) = k20c(i)*2.d0/1.d0 * exp( -ee/tdt * anu )

      end do


                                !! k21 &  k23k21c &  k24k21c & k34k21c

      k21x = 2.49d-11
      k21xb = k21x
      k21xa = 3.d0/2.d0 * k21xb
      k21xc = k21xb / 2.d0

      k21xa = k21xa * rf21a
      k21xb = k21xb * rf21b
      k21xc = k21xc * rf21c

      do i = 1, nisot
	 k21b(i) = k21xb
	 k21c(i) = k21xc
	 k21bp(i) = k21b(i) *  &
            exp( -ee/tdt* dble(nu(i,2)-nu(i,1)-nu(1,1)) )
	 k21cp(i) = k21c(i) *  &
            exp( -ee/tdt * dble(nu(i,1)-nu(1,1)) )
      end do

      k23k21c = k21xc
      k24k21c = k21xc
      k34k21c = k21xc
      k23k21cp = k23k21c*2.d0/2.d0 *  &
           exp( -ee/tdt* dble(nu(2,1)-nu(3,1)) )
      k24k21cp = k24k21c*2.d0/2.d0 *  &
           exp( -ee/tdt* dble(nu(2,1)-nu(4,1)) )
      k34k21cp = k34k21c*2.d0/2.d0 *  &
           exp( -ee/tdt* dble(nu(3,1)-nu(4,1)) )


                                !! k33

      k33c = k21x * rf33bc
      do i=2,nisot
	 k33cp(i) = k33c *  &
           exp( -ee/tdt * dble(nu(1,2)-nu(1,1)-nu(i,1)) )
      end do

      end subroutine GETK_dlvr11

!***********************************************************************
!  3-datasets not readin: hid27-1.dat, hid28-1.dat, hid36-1.dat
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fast scheme for NLTE cooling rates at 15um by CO2 in a Martian GCM !
!                 Version dlvr11_03. 2012.                           !
! Software written and provided by IAA/CSIC, Granada, Spain,         !
! under ESA contract "Mars Climate Database and Physical Models"     !
! Person of contact: Miguel Angel Lopez Valverde  valverde@iaa.es    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     **** nlte_aux.F subroutines *****
!
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fast scheme for NLTE cooling rates at 15um by CO2 in a Martian GCM !
!                 Version dlvr11_03. 2012.                           !
! Software written and provided by IAA/CSIC, Granada, Spain,         !
! under ESA contract "Mars Climate Database and Physical Models"     !
! Person of contact: Miguel Angel Lopez Valverde  valverde@iaa.es    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**********************************************************************

!     Includes the following old 1-D model files/subroutines

!**********************************************************************


!     *** Old MZTCRSUB_dlvr11.f ***

!************************************************************************

!***********************************************************************
      function planckdp(tp,xnu)
!***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

      real*8 :: planckdp
      real*8 :: xnu
      real :: tp

!     planckdp = gamma*xnu**3.0d0 / exp( ee*xnu/dble(tp) )
                                !erg cm-2.sr-1/cm-1.
      planckdp = gamma1*xnu**3.0d0 / exp( ee*xnu/dble(tp) )
                                !erg cm-2.sr-1/cm-1.

      end function planckdp

!***********************************************************************
      subroutine leetvt

!***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     local variables
      integer :: i
      real*8 ::	zld(nl), zyd(nzy)
      real*8  ::xvt11(nzy), xvt21(nzy), xvt31(nzy), xvt41(nzy)

!***********************************************************************

      do i=1,nzy
         zyd(i) = dble(zy(i))
         xvt11(i)= dble( ty(i) )
         xvt21(i)= dble( ty(i) )
         xvt31(i)= dble( ty(i) )
         xvt41(i)= dble( ty(i) )
      end do

      do i=1,nl
         zld(i) = dble( zl(i) )
      enddo
      call interhuntdp4veces ( v626t1,v628t1,v636t1,v627t1, zld,nl,  &
          xvt11, xvt21, xvt31, xvt41, zyd,nzy, 1 )

      end subroutine leetvt


!     ****************************************************************
!     *** MZTFSUB_dlvr11_02.f ***
!     ****************************************************************
      subroutine initial

!     ****************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     local variables
      integer :: 	i

!     ***************

      eqw = 0.0d00
      aa = 0.0d00
      cc = 0.0d00
      dd = 0.0d00

      do i=1,nbox
         ccbox(i) = 0.0d0
         ddbox(i) = 0.0d0
      end do

      end subroutine initial

!     **********************************************************************

      subroutine intershphunt (i, alsx,adx,xtemp)

!     **********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      real*8 :: alsx(nbox_max),adx(nbox_max) ! Output
      real*8 :: xtemp(nbox_max)    ! Input
      integer ::    i              ! I , O

!     local variables
      integer :: k
      real*8  :: factor
      real*8  :: temperatura     ! para evitar valores ligeramnt out of limits

!     ***********

      do 1, k=1,nbox
         temperatura = xtemp(k)
         if (abs(xtemp(k)-thist(1)).le.0.01d0) then
            temperatura=thist(1)
         elseif (abs(xtemp(k)-thist(nhist)).le.0.01d0) then
            temperatura=thist(nhist)
         elseif (xtemp(k).lt.thist(1)) then
            temperatura=thist(1)
            write (*,*) ' WARNING intershphunt/ Too low atmosph Tk:'
            write (*,*) ' WARNING      k,xtemp = ', k,xtemp(k)
            write (*,*) ' Minimum Tk in histogram used : ', thist(1)
         elseif (xtemp(k).gt.thist(nhist)) then
            temperatura=thist(nhist)
            write (*,*) ' WARNING intershphunt/ Very high atmosph Tk:'
            write (*,*) ' WARNING      k,xtemp = ', k,xtemp(k)
            write (*,*) ' Max Tk in histogram used : ', thist(nhist)
         endif
         call huntdp ( thist,nhist, temperatura, i )
         if ( i.eq.0 .or. i.eq.nhist ) then
	    write (*,*) ' HUNT/ Limits input grid:',  &
                 thist(1),thist(nhist)
	    write (*,*) ' HUNT/ location in grid:', xtemp(k)
            stop ' INTERSHP/ Interpolation error. T out of Histogram.'
         endif
         factor = 1.d0 /  (thist(i+1)-thist(i))
         alsx(k) = (( xls1(i,k)*(thist(i+1)-xtemp(k)) +  &
              xls1(i+1,k)*(xtemp(k)-thist(i)) )) * factor
         adx(k)  = (( xld1(i,k)*(thist(i+1)-xtemp(k)) +  &
              xld1(i+1,k)*(xtemp(k)-thist(i)) )) * factor
 1    continue

      end subroutine intershphunt

!     **********************************************************************

      subroutine interstrhunt (i, stx, ts, sx, xtemp )

!     **********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      real*8  ::	stx     ! output, total band strength
      real*8  ::	ts      ! input, temp for stx
      real*8  ::        sx(nbox_max) ! output, strength for each box
      real*8  ::	xtemp(nbox_max) ! input, temp for sx
      integer  :: 	i

!     local variables
      integer :: 	k
      real*8  ::        factor
      real*8  ::        temperatura

!     ***********

      do 1, k=1,nbox
         temperatura = xtemp(k)
         if (abs(xtemp(k)-thist(1)).le.0.01d0) then
            temperatura=thist(1)
         elseif (abs(xtemp(k)-thist(nhist)).le.0.01d0) then
            temperatura=thist(nhist)
         elseif (xtemp(k).lt.thist(1)) then
            temperatura=thist(1)
            write (*,*) ' WARNING interstrhunt/ Too low atmosph Tk:'
	    write (*,*) ' WARNING     k,xtemp(k) = ', k,xtemp(k)
	    write (*,*) ' Minimum Tk in histogram used : ', thist(1)
         elseif (xtemp(k).gt.thist(nhist)) then
            temperatura=thist(nhist)
            write (*,*) ' WARNING interstrhunt/ Very high atmosph Tk:'
	    write (*,*) ' WARNING     k,xtemp(k) = ', k,xtemp(k)
	    write (*,*) ' Max Tk in histogram used : ', thist(nhist)
         endif
         call huntdp ( thist,nhist, temperatura, i )
         if ( i.eq.0 .or. i.eq.nhist ) then
	    write(*,*)'HUNT/ Limits input grid:',  &
                thist(1),thist(nhist)
	    write(*,*)'HUNT/ location in grid:',xtemp(k)
            stop 'INTERSTR/1/ Interpolation error. T out of Histogram.'
         endif
         factor = 1.d0 /  (thist(i+1)-thist(i))
         sx(k) = ( sk1(i,k)   * (thist(i+1)-xtemp(k))  &
              + sk1(i+1,k) * (xtemp(k)-thist(i))  ) * factor
 1    continue


      temperatura = ts
      if (abs(ts-thist(1)).le.0.01d0) then
         temperatura=thist(1)
      elseif (abs(ts-thist(nhist)).le.0.01d0) then
         temperatura=thist(nhist)
      elseif (ts.lt.thist(1)) then
         temperatura=thist(1)
         write (*,*) ' WARNING interstrhunt/ Too low atmosph Tk:'
         write (*,*) ' WARNING            ts = ', temperatura
         write (*,*) ' Minimum Tk in histogram used : ', thist(1)
      elseif (ts.gt.thist(nhist)) then
         temperatura=thist(nhist)
         write (*,*) ' WARNING interstrhunt/ Very high atmosph Tk:'
         write (*,*) ' WARNING            ts = ', temperatura
         write (*,*) ' Max Tk in histogram used : ', thist(nhist)
      endif
      call huntdp ( thist,nhist, temperatura, i )
      if ( i.eq.0 .or. i.eq.nhist ) then
         write (*,*) ' HUNT/ Limits input grid:',  &
              thist(1),thist(nhist)
         write (*,*) ' HUNT/ location in grid:', ts
         stop ' INTERSTR/2/ Interpolat error. T out of Histogram.'
      endif
      factor = 1.d0 /  (thist(i+1)-thist(i))
      stx = 0.d0
      do k=1,nbox
         stx = stx + no(k) * ( sk1(i,k)*(thist(i+1)-ts) +  &
              sk1(i+1,k)*(ts-thist(i)) ) * factor
      end do

      end subroutine interstrhunt

!     **********************************************************************

      subroutine intzhunt (k, h, aco2,ap,amr,at, con)

!     k lleva la posicion de la ultima llamada a intz , necesario para
!     que esto represente una aceleracion real.
!     **********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      real  ::   	h       ! i
      real*8  ::	con(nzy) ! i
      real*8  ::	aco2, ap, at, amr ! o
      integer  ::	k       ! i

!     local variables
      real  ::          factor

!     ************

      call hunt ( zy,nzy, h, k )
      factor =  (h-zy(k)) /  (zy(k+1)-zy(k))
      ap = dble( exp( log(py(k)) + log(py(k+1)/py(k)) * factor ) )
      aco2 = dlog(con(k)) + dlog( con(k+1)/con(k) ) * dble(factor)
      aco2 = exp( aco2 )
      at = dble( ty(k) + (ty(k+1)-ty(k)) * factor )
      amr = dble( mr(k) + (mr(k+1)-mr(k)) * factor )

      end subroutine intzhunt

!     **********************************************************************

      subroutine intzhunt_cts (k,h,nzy_cts_real,aco2,ap,amr,at,con)

!     k lleva la posicion de la ultima llamada a intz , necesario para
!     que esto represente una aceleracion real.
!     **********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      real  ::		h       ! i
      real*8  ::	con(nzy_cts) ! i
      real*8  ::	aco2, ap, at, amr ! o
      integer  ::	k       ! i
      integer  ::       nzy_cts_real ! i

!     local variables
      real  ::          factor

!     ************

      call hunt_cts ( zy_cts,nzy_cts, nzy_cts_real, h, k )
      factor =  (h-zy_cts(k)) /  (zy_cts(k+1)-zy_cts(k))
      ap = dble( exp( log(py_cts(k)) +  &
           log(py_cts(k+1)/py_cts(k)) * factor ) )
      aco2 = dlog(con(k)) + dlog( con(k+1)/con(k) ) * dble(factor)
      aco2 = exp( aco2 )
      at = dble( ty_cts(k) + (ty_cts(k+1)-ty_cts(k)) * factor )
      amr = dble( mr_cts(k) + (mr_cts(k+1)-mr_cts(k)) * factor )


      end subroutine intzhunt_cts


!     **********************************************************************

!     real*8 function we_clean  ( y,pl, xalsa, xalda )
      real*8 function we_clean  ( y,pl2, xalsa, xalda )

!     **********************************************************************

      use ModPlanet
      use ModConstants, only : pi
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      real*8  ::	  y         ! I. path's absorber amount * strength
!     real*8  ::          pl        ! I. path's partial pressure of CO2
      real*8  ::          pl2        ! I. path's partial pressure of CO2
      real*8  ::          xalsa     ! I.  Self lorentz linewidth for 1 isot & 1 box
      real*8  ::          xalda     ! I.  Doppler linewidth        "           "

!     local variables
      integer  :: 	i
      real*8   ::	x,wl,wd,wvoigt
      real*8   ::	cn(0:7),dn(0:7)
      real*8   ::       factor, denom
!     real*8   ::       pi, pi2, sqrtpi
      real*8   ::       pi2, sqrtpi

!     data blocks
      data cn/9.99998291698d-1,-3.53508187098d-1,9.60267807976d-2,  &
           -2.04969011013d-2,3.43927368627d-3,-4.27593051557d-4,  &
           3.42209457833d-5,-1.28380804108d-6/
      data dn/1.99999898289,5.774919878d-1,-5.05367549898d-1,  &
           8.21896973657d-1,-2.5222672453,6.1007027481,  &
           -8.51001627836,4.6535116765/

!     ***********

!     pi = 3.141592
      pi2= 6.28318531
      sqrtpi = 1.77245385

!     x=y / ( pi2 * xalsa*pl )
      x=y / ( pi2 * xalsa*pl2 )


!     Lorentz
      wl=y/sqrt(1.0d0+pi*x/2.0d0)

!     Doppler
      x = y / (xalda*sqrtpi)
      if (x .lt. 5.0d0) then
         wd = cn(0)
         factor = 1.d0
         do i=1,7
            factor = factor * x
	    wd = wd + cn(i) * factor
         end do
         wd = xalda * x * sqrtpi * wd
      else
         wd = dn(0)
         factor = 1.d0 / log(x)
         denom = 1.d0
         do i=1,7
            denom = denom * factor
	    wd = wd + dn(i) * denom
         end do
         wd = xalda * sqrt(log(x)) * wd
      end if

!     Voigt
      wvoigt = wl*wl + wd*wd - (wd*wl/y)*(wd*wl/y)

      if ( wvoigt .lt. 0.0d0 ) then
       write (*,*) ' Subroutine WE/ Error in Voift EQS calculation'
       write (*,*) '  WL, WD, X, Y = ', wl, wd, x, y
       stop '  ERROR : Imaginary EQW. Revise spectral data. '
      endif

      we_clean = sqrt( wvoigt )

      end function we_clean

!     ***********************************************************************

      subroutine mztf_correccion (coninf, con, ib )

!     ***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      integer  ::	ib
      real*8  ::	con(nzy), coninf

!     local variables
      integer  :: 	i, isot
      real*8   ::	tvt0(nzy), tvtbs(nzy), zld(nl),zyd(nzy)
      real*8  ::        xqv, xes, xlower, xfactor

!     *********

      isot = 1
      nu11 = dble( nu(1,1) )

      do i=1,nzy
         zyd(i) = dble(zy(i))
      enddo
      do i=1,nl
         zld(i) = dble( zl(i) )
      end do

!     tvtbs
      call interhuntdp (tvtbs,zyd,nzy, v626t1,zld,nl, 1 )

!     tvt0
      if (ib.eq.2 .or. ib.eq.3 .or. ib.eq.4) then
         call interhuntdp (tvt0,zyd,nzy, v626t1,zld,nl, 1 )
      else
         do i=1,nzy
            tvt0(i) = dble( ty(i) )
         end do
      end if

!     factor
      do i=1,nzy

         xlower = exp( ee*dble(elow(isot,ib)) *  &
              ( 1.d0/dble(ty(i))-1.d0/tvt0(i) ) )
         xes = 1.0d0
         xqv = ( 1.d0-exp( -ee*nu11/tvtbs(i) ) ) /  &
              (1.d0-exp( -ee*nu11/dble(ty(i)) ))
         xfactor = xlower * xqv**2.d0 * xes

         con(i) = con(i) * xfactor
         if (i.eq.nzy) coninf = coninf * xfactor

      end do

      end subroutine mztf_correccion

!    ***********************************************************************

      subroutine mzescape_normaliz ( taustar, istyle )

!     ***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      real*8  :: 	taustar(nl) ! o
      integer  ::       istyle    ! i

!     local variables and constants
      integer  :: 	i, imaximum
      real*8   ::       maximum

!     ***************

      taustar(nl) = taustar(nl-1)

      if ( istyle .eq. 1 ) then
         imaximum = nl
         maximum = taustar(nl)
         do i=1,nl-1
	    if (taustar(i).gt.maximum) taustar(i) = taustar(nl)
         enddo
      elseif ( istyle .eq. 2 ) then
         imaximum = nl
         maximum = taustar(nl)
         do i=nl-1,1,-1
	    if (taustar(i).gt.maximum) then
	       maximum = taustar(i)
	       imaximum = i
	    endif
         enddo
         do i=imaximum,nl
	    if (taustar(i).lt.maximum) taustar(i) = maximum
         enddo
      endif

      do i=1,nl
         taustar(i) = taustar(i) / maximum
      enddo

      end subroutine mzescape_normaliz

!     ***********************************************************************

      subroutine mzescape_normaliz_02 ( taustar, nn, istyle )

!     ***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      real*8  :: 	taustar(nn) ! i,o
      integer ::        istyle    ! i
      integer ::        nn        ! i

!     local variables and constants
      integer  :: 	i, imaximum
      real*8   ::       maximum

!     ***************

      taustar(nn) = taustar(nn-1)

      if ( istyle .eq. 1 ) then
         imaximum = nn
         maximum = taustar(nn)
         do i=1,nn-1
	    if (taustar(i).gt.maximum) taustar(i) = taustar(nn)
         enddo
      elseif ( istyle .eq. 2 ) then
         imaximum = nn
         maximum = taustar(nn)
         do i=nn-1,1,-1
	    if (taustar(i).gt.maximum) then
	       maximum = taustar(i)
	       imaximum = i
	    endif
         enddo
         do i=imaximum,nn
	    if (taustar(i).lt.maximum) taustar(i) = maximum
         enddo
      endif

      do i=1,nn
         taustar(i) = taustar(i) / maximum
      enddo

      end subroutine mzescape_normaliz_02


!***********************************************************************
!     *** interdp_ESCTVCISO_dlvr11.f ***
!***********************************************************************

      subroutine interdp_ESCTVCISO

!***********************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     local variables
      integer  ::    i
      real*8   ::    lnpnb(nl)

!***********************************************************************

!     Use pressure in the NLTE grid but in log and in nb
      do i=1,nl
         lnpnb(i) = log( dble( pl(i) * 1013.25 * 1.e6) )
      enddo

!     Interpolations

      call interhuntdp3veces  &
           ( taustar21,taustar31,taustar41,    lnpnb, nl,  &
           tstar21tab,tstar31tab,tstar41tab, lnpnbtab, nztabul, &
           1 )

      call interhuntdp3veces ( vc210,vc310,vc410, lnpnb, nl, &
           vc210tab,vc310tab,vc410tab, lnpnbtab, nztabul, 2 )

      end subroutine interdp_ESCTVCISO


! ***********************************************************************
!     *** hunt_cts.f ***
! ***********************************************************************

      SUBROUTINE hunt_cts(xx,n,n_cts,x,jlo)
!
!     La dif con hunt es el uso de un indice superior (n_cts) mas bajito que (n)
!
!     Arguments
      INTEGER :: jlo               ! O
      INTEGER :: n                 ! I
      INTEGER :: n_cts             ! I
      REAL  :: xx(n)               ! I
      REAL  :: x                   ! I

!     Local variables
      INTEGER  :: inc,jhi,jm
      LOGICAL :: ascnd
!
      ascnd=xx(n_cts).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n_cts)then
         jlo=0
         jhi=n_cts+1
         goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
 1       jhi=jlo+inc
!     write (*,*) jlo
         if(jhi.gt.n_cts)then
            jhi=n_cts+1
!     write (*,*) jhi-1
         else if(x.ge.xx(jhi).eqv.ascnd)then
            jlo=jhi
            inc=inc+inc
!     write (*,*) jlo
            goto 1
         endif
      else
         jhi=jlo
 2       jlo=jhi-inc
!     write (*,*) jlo
         if(jlo.lt.1)then
            jlo=0
         else if(x.lt.xx(jlo).eqv.ascnd)then
            jhi=jlo
            inc=inc+inc
            goto 2
         endif
      endif
 3    if(jhi-jlo.eq.1)then
         if(x.eq.xx(n_cts))jlo=n_cts-1
         if(x.eq.xx(1))jlo=1
!     write (*,*) jlo
         return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
         jlo=jm
      else
         jhi=jm
      endif
!     write (*,*) jhi-1
      goto 3
!
      END SUBROUTINE hunt_cts


! *****************************************************
!     *** huntdp.f ***
! *****************************************************

      SUBROUTINE huntdp(xx,n,x,jlo)
!
!     Arguments
      INTEGER :: jlo  ! O
      INTEGER :: n    ! I
      REAL*8  :: xx(n)             ! I
      REAL*8  :: x                 ! I

!     Local variables
      INTEGER :: inc,jhi,jm
      LOGICAL :: ascnd
!
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
         jlo=0
         jhi=n+1
         goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
 1       jhi=jlo+inc
         if(jhi.gt.n)then
            jhi=n+1
         else if(x.ge.xx(jhi).eqv.ascnd)then
            jlo=jhi
            inc=inc+inc
            goto 1
         endif
      else
         jhi=jlo
 2       jlo=jhi-inc
         if(jlo.lt.1)then
            jlo=0
         else if(x.lt.xx(jlo).eqv.ascnd)then
            jhi=jlo
            inc=inc+inc
            goto 2
         endif
      endif
 3    if(jhi-jlo.eq.1)then
         if(x.eq.xx(n))jlo=n-1
         if(x.eq.xx(1))jlo=1
         return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
         jlo=jm
      else
         jhi=jm
      endif
      goto 3

      END SUBROUTINE huntdp

! *****************************************************
!     *** hunt.f ***
! *****************************************************

      SUBROUTINE hunt(xx,n,x,jlo)

!     Arguments
      INTEGER :: jlo               ! O
      INTEGER :: n                 ! I
      REAL  ::   xx(n)               ! I
      REAL  ::   x                   ! I

!     Local variables
      INTEGER :: inc,jhi,jm
      LOGICAL :: ascnd
!

      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
         jlo=0
         jhi=n+1
         goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
 1       jhi=jlo+inc
!     write (*,*) jlo
         if(jhi.gt.n)then
            jhi=n+1
!     write (*,*) jhi-1
         else if(x.ge.xx(jhi).eqv.ascnd)then
            jlo=jhi
            inc=inc+inc
!     write (*,*) jlo
            goto 1
         endif
      else
         jhi=jlo
 2       jlo=jhi-inc
!     write (*,*) jlo
         if(jlo.lt.1)then
            jlo=0
         else if(x.lt.xx(jlo).eqv.ascnd)then
            jhi=jlo
            inc=inc+inc
            goto 2
         endif
      endif
 3    if(jhi-jlo.eq.1)then
         if(x.eq.xx(n))jlo=n-1
         if(x.eq.xx(1))jlo=1
!     write (*,*) jlo
         return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
         jlo=jm
      else
         jhi=jm
      endif
!     write (*,*) jhi-1
      goto 3

      END SUBROUTINE hunt

! **************************************************************************
!     *** interdp_limits.f ***
! **************************************************************************

      subroutine interdp_limits ( yy, zz, m,   i1,i2,  &
           y,  z, n,   j1,j2,  opt)

! **************************************************************************
!     Interpolation soubroutine.
!     Returns values between indexes i1 & i2, donde  1 =< i1 =< i2 =< m
!     Solo usan los indices de los inputs entre j1,j2, 1 =< j1 =< j2 =< n
!     Input values: y(n) , z(n)  (solo se usarann los valores entre j1,j2)
!     zz(m) (solo se necesita entre i1,i2)
!     Output values: yy(m) (solo se calculan entre i1,i2)
!     Options:    opt=1 -> lineal ,,  opt=2 -> logarithmic
!     Difference with interdp:
!     here interpolation proceeds between indexes i1,i2 only
!     if i1=1 & i2=m, both subroutines are exactly the same
!     thus previous calls to interdp or interdp2 could be easily replaced
!     JAN 98 	MALV 		Version for mz1d
!     ***********************************************************************

      implicit none

!     Arguments
      integer 	:: n,m             ! I. Dimensions
      integer 	:: i1, i2, j1, j2, opt ! I
      real*8    :: zz(m)   ! I
      real*8    :: yy(m)   ! O
      real*8    :: z(n),y(n) ! I

!     Local variables
      integer  :: 	i,j
      real*8   ::	zmin,zzmin,zmax,zzmax

!     *******************************

!     write (*,*) ' d interpolating '
!     call mindp_limits (z,n,zmin, j1,j2)
!     call mindp_limits (zz,m,zzmin, i1,i2)
!     call maxdp_limits (z,n,zmax, j1,j2)
!     call maxdp_limits (zz,m,zzmax, i1,i2)
      zmin=minval(z(j1:j2))
      zzmin=minval(zz(i1:i2))
      zmax=maxval(z(j1:j2))
      zzmax=maxval(zz(i1:i2))

! *********************************************************************
!  Throwing error for equal quantities (to 15 significant digits)
!  (test by turning off check): S.Bougher (12/13/2017)
!     if(zzmin.lt.zmin)then
!        write (*,*) 'from d interp: new variable out of limits'
!        write (*,*) zzmin,'must be .ge. ',zmin
!        stop
! *********************************************************************
!     elseif(zzmax.gt.zmax)then
!     type *,'from interp: new variable out of limits'
!     type *,zzmax, 'must be .le. ',zmax
!     stop
!     end if
! *********************************************************************

      do 1,i=i1,i2

         do 2,j=j1,j2-1
            if(zz(i).ge.z(j).and.zz(i).lt.z(j+1)) goto 3
 2       continue

!     in this case (zz(i2).eq.z(j2)) and j leaves the loop with j=j2-1+1=j2
         if(opt.eq.1)then
            yy(i)=y(j2-1)+(y(j2)-y(j2-1))*(zz(i)-z(j2-1))/  &
                 (z(j2)-z(j2-1))
         elseif(opt.eq.2)then
            if(y(j2).eq.0.0d0.or.y(j2-1).eq.0.0d0)then
               yy(i)=0.0d0
            else
               yy(i)=exp(log(y(j2-1))+log(y(j2)/y(j2-1))*  &
                    (zz(i)-z(j2-1))/(z(j2)-z(j2-1)))
            end if
         else
            write (*,*) ' d interp : opt must be 1 or 2, opt= ',opt
         end if
         goto 1
 3       continue
         if(opt.eq.1)then
            yy(i)=y(j)+(y(j+1)-y(j))*(zz(i)-z(j))/(z(j+1)-z(j))
!     type *, ' '
!     type *, ' z(j),z(j+1) =', z(j),z(j+1)
!     type *, ' t(j),t(j+1) =', y(j),y(j+1)
!     type *, ' zz, tt =  ', zz(i), yy(i)
         elseif(opt.eq.2)then
            if(y(j+1).eq.0.0d0.or.y(j).eq.0.0d0)then
               yy(i)=0.0d0
            else
               yy(i)=exp(log(y(j))+log(y(j+1)/y(j))*  &
                    (zz(i)-z(j))/(z(j+1)-z(j)))
            end if
         else
            write (*,*) ' interp : opt must be 1 or 2, opt= ',opt
         end if
 1    continue

      end subroutine interdp_limits

! ***************************************************************************
!     *** interhunt2veces.f ***
! ***************************************************************************

      subroutine interhunt2veces ( y1,y2,  zz,m,  &
          x1,x2,  z,n,  opt)

!     interpolation soubroutine basada en Numerical Recipes HUNT.FOR
!     input values: y(n) at z(n)
!     output values: yy(m) at zz(m)
!     options: 1 -> lineal
!     2 -> logarithmic
!     ***********************************************************************

      implicit none

!     Arguments
      integer	:: n,m,opt         ! I
      real	:: zz(m),z(n)      ! I
      real      :: y1(m),y2(m)       ! O
      real      :: x1(n),x2(n)       ! I


!     Local variables
      integer :: i, j
      real    :: factor
      real    :: zaux
      
!!!!

      j = 1                     ! initial first guess (=n/2 is anothr pssblty)

      do 1,i=1,m                !

                                ! Busca indice j donde ocurre q zz(i) esta entre [z(j),z(j+1)]
         zaux = zz(i)
         if (abs(zaux-z(1)).le.0.01) then
            zaux=z(1)
         elseif (abs(zaux-z(n)).le.0.01) then
            zaux=z(n)
         endif
         call hunt ( z,n, zaux, j )
     
         if ( j.eq.0 .or. j.eq.n ) then
	    write (*,*) ' HUNT/ Limits input grid:', z(1),z(n)
	    write (*,*) ' HUNT/ location in new grid:', zz(i)
            stop ' interhunt2/ Interpolat error. zz out of limits.'
         endif

                                ! Perform interpolation
         factor = (zz(i)-z(j))/(z(j+1)-z(j))
         if (opt.eq.1) then
	    y1(i) = x1(j) + (x1(j+1)-x1(j)) * factor
	    y2(i) = x2(j) + (x2(j+1)-x2(j)) * factor
         else
	    y1(i) = exp( log(x1(j)) + log(x1(j+1)/x1(j)) * factor )
	    y2(i) = exp( log(x2(j)) + log(x2(j+1)/x2(j)) * factor )
         end if

 1    continue

      end subroutine interhunt2veces

! ***************************************************************************
!     *** interhunt5veces.f ***
! ***************************************************************************

      subroutine interhunt5veces ( y1,y2,y3,y4,y5,  zz,m,  &
           x1,x2,x3,x4,x5,  z,n,  opt)

! ***************************************************************************
!     interpolation soubroutine basada en Numerical Recipes HUNT.FOR
!     input values: y(n) at z(n)
!     output values: yy(m) at zz(m)
!     options: 1 -> lineal
!     2 -> logarithmic
!     ***********************************************************************

      implicit none
!     Arguments
      integer	:: n,m,opt         ! I
      real	:: zz(m),z(n)      ! I
      real      :: y1(m),y2(m),y3(m),y4(m),y5(m) ! O
      real      :: x1(n),x2(n),x3(n),x4(n),x5(n) ! I


!     Local variables
      integer :: i, j
      real    :: factor
      real    :: zaux

!!!!

      j = 1                     ! initial first guess (=n/2 is anothr pssblty)

      do 1,i=1,m                !

                                ! Busca indice j donde ocurre q zz(i) esta entre [z(j),z(j+1)]
         zaux = zz(i)
         if (abs(zaux-z(1)).le.0.01) then
            zaux=z(1)
         elseif (abs(zaux-z(n)).le.0.01) then
            zaux=z(n)
         endif
         call hunt ( z,n, zaux, j )
         if ( j.eq.0 .or. j.eq.n ) then
	    write (*,*) ' HUNT/ Limits input grid:', z(1),z(n)
	    write (*,*) ' HUNT/ location in new grid:', zz(i)
            stop ' interhunt5/ Interpolat error. zz out of limits.'
         endif

                                ! Perform interpolation
         factor = (zz(i)-z(j))/(z(j+1)-z(j))
         if (opt.eq.1) then
	    y1(i) = x1(j) + (x1(j+1)-x1(j)) * factor
	    y2(i) = x2(j) + (x2(j+1)-x2(j)) * factor
	    y3(i) = x3(j) + (x3(j+1)-x3(j)) * factor
	    y4(i) = x4(j) + (x4(j+1)-x4(j)) * factor
	    y5(i) = x5(j) + (x5(j+1)-x5(j)) * factor
         else
	    y1(i) = exp( log(x1(j)) + log(x1(j+1)/x1(j)) * factor )
	    y2(i) = exp( log(x2(j)) + log(x2(j+1)/x2(j)) * factor )
	    y3(i) = exp( log(x3(j)) + log(x3(j+1)/x3(j)) * factor )
	    y4(i) = exp( log(x4(j)) + log(x4(j+1)/x4(j)) * factor )
	    y5(i) = exp( log(x5(j)) + log(x5(j+1)/x5(j)) * factor )
         end if

 1    continue

      end subroutine interhunt5veces

! ***************************************************************************
!     *** interhuntdp3veces.f ***

! ***************************************************************************

      subroutine interhuntdp3veces ( y1,y2,y3, zz,m,  &
          x1,x2,x3,  z,n,  opt)

! ***************************************************************************
!     interpolation soubroutine basada en Numerical Recipes HUNT.FOR
!     input values: x(n) at z(n)
!     output values: y(m) at zz(m)
!     options: opt = 1 -> lineal
!     opt=/=1 -> logarithmic
!     ***********************************************************************
!     Arguments
      integer	:: n,m,opt         ! I
      real*8	:: zz(m),z(n)      ! I
      real*8    :: y1(m),y2(m),y3(m) ! O
      real*8    :: x1(n),x2(n),x3(n) ! I


!     Local variables
      integer ::  i, j
      real*8  ::  factor
      real*8  ::  zaux

!!!!

      j = 1                     ! initial first guess (=n/2 is anothr pssblty)

      do 1,i=1,m                !

                                ! Busca indice j donde ocurre q zz(i) esta entre [z(j),z(j+1)]
         zaux = zz(i)
         if (abs(zaux-z(1)).le.0.01d0) then
            zaux=z(1)
         elseif (abs(zaux-z(n)).le.0.01d0) then
            zaux=z(n)
         endif
         call huntdp ( z,n, zaux, j )
         if ( j.eq.0 .or. j.eq.n ) then
	    write (*,*) ' HUNT/ Limits input grid:', z(1),z(n)
	    write (*,*) ' HUNT/ location in new grid:', zz(i)
            stop ' INTERHUNTDP3/ Interpolat error. zz out of limits.'
         endif

                                ! Perform interpolation
         factor = (zz(i)-z(j))/(z(j+1)-z(j))
         if (opt.eq.1) then
	    y1(i) = x1(j) + (x1(j+1)-x1(j)) * factor
	    y2(i) = x2(j) + (x2(j+1)-x2(j)) * factor
	    y3(i) = x3(j) + (x3(j+1)-x3(j)) * factor
         else
	    y1(i) = dexp( dlog(x1(j)) + dlog(x1(j+1)/x1(j)) * factor )
	    y2(i) = dexp( dlog(x2(j)) + dlog(x2(j+1)/x2(j)) * factor )
	    y3(i) = dexp( dlog(x3(j)) + dlog(x3(j+1)/x3(j)) * factor )
         end if

 1    continue

      end subroutine interhuntdp3veces

! ***************************************************************************
!     *** interhuntdp4veces.f ***
! ***************************************************************************

      subroutine interhuntdp4veces ( y1,y2,y3,y4, zz,m,  &
           x1,x2,x3,x4,  z,n,  opt)

! ***************************************************************************
!     interpolation soubroutine basada en Numerical Recipes HUNT.FOR
!     input values: x1(n),x2(n),x3(n),x4(n) at z(n)
!     output values: y1(m),y2(m),y3(m),y4(m) at zz(m)
!     options: 1 -> lineal
!     2 -> logarithmic
!     ***********************************************************************

      implicit none

!     Arguments
      integer	:: n,m,opt         ! I
      real*8	:: zz(m),z(n)      ! I
      real*8    :: y1(m),y2(m),y3(m),y4(m) ! O
      real*8    :: x1(n),x2(n),x3(n),x4(n) ! I


!     Local variables
      integer :: i, j
      real*8  ::   factor
      real*8  ::  zaux

!!!!

      j = 1                     ! initial first guess (=n/2 is anothr pssblty)

      do 1,i=1,m                !

                                ! Caza del indice j donde ocurre que zz(i) esta entre [z(j),z(j+1)]
         zaux = zz(i)
         if (abs(zaux-z(1)).le.0.01d0) then
            zaux=z(1)
         elseif (abs(zaux-z(n)).le.0.01d0) then
            zaux=z(n)
         endif
         call huntdp ( z,n, zaux, j )
         if ( j.eq.0 .or. j.eq.n ) then
	    write (*,*) ' HUNT/ Limits input grid:', z(1),z(n)
	    write (*,*) ' HUNT/ location in new grid:', zz(i)
            stop ' INTERHUNTDP4/ Interpolat error. zz out of limits.'
         endif

                                ! Perform interpolation
         factor = (zz(i)-z(j))/(z(j+1)-z(j))
         if (opt.eq.1) then
	    y1(i) = x1(j) + (x1(j+1)-x1(j)) * factor
	    y2(i) = x2(j) + (x2(j+1)-x2(j)) * factor
	    y3(i) = x3(j) + (x3(j+1)-x3(j)) * factor
	    y4(i) = x4(j) + (x4(j+1)-x4(j)) * factor
         else
	    y1(i) = dexp( dlog(x1(j)) + dlog(x1(j+1)/x1(j)) * factor )
	    y2(i) = dexp( dlog(x2(j)) + dlog(x2(j+1)/x2(j)) * factor )
	    y3(i) = dexp( dlog(x3(j)) + dlog(x3(j+1)/x3(j)) * factor )
	    y4(i) = dexp( dlog(x4(j)) + dlog(x4(j+1)/x4(j)) * factor )
         end if

 1    continue

      end subroutine interhuntdp4veces

! **************************************************************************
!     *** interhuntdp.f ***
! **************************************************************************

      subroutine interhuntdp ( y1, zz,m, &
           x1,  z,n,  opt)

! **************************************************************************
!     interpolation soubroutine basada en Numerical Recipes HUNT.FOR
!     input values: x1(n) at z(n)
!     output values: y1(m) at zz(m)
!     options: 1 -> lineal
!     2 -> logarithmic
!     ***********************************************************************

      implicit none

!     Arguments
      integer	:: n,m,opt         ! I
      real*8	:: zz(m),z(n)      ! I
      real*8    :: y1(m)           ! O
      real*8    :: x1(n)           ! I


!     Local variables
      integer :: i, j
      real*8  :: factor
      real*8  :: zaux

!!!!

      j = 1                     ! initial first guess (=n/2 is anothr pssblty)

      do 1,i=1,m                !

                                ! Caza del indice j donde ocurre que zz(i) esta entre [z(j),z(j+1)]
         zaux = zz(i)
         if (abs(zaux-z(1)).le.0.01d0) then
            zaux=z(1)
         elseif (abs(zaux-z(n)).le.0.01d0) then
            zaux=z(n)
         endif
         call huntdp ( z,n, zaux, j )
         if ( j.eq.0 .or. j.eq.n ) then
	    write (*,*) ' HUNT/ Limits input grid:', z(1),z(n)
	    write (*,*) ' HUNT/ location in new grid:', zz(i)
            stop ' INTERHUNT/ Interpolat error. zz out of limits.'
         endif

                                ! Perform interpolation
         factor = (zz(i)-z(j))/(z(j+1)-z(j))
         if (opt.eq.1) then
	    y1(i) = x1(j) + (x1(j+1)-x1(j)) * factor
         else
	    y1(i) = dexp( dlog(x1(j)) + dlog(x1(j+1)/x1(j)) * factor )
         end if

 1    continue

      end subroutine interhuntdp


! **********************************************************************
!     *** interhunt.f ***
! **********************************************************************

      subroutine interhunt ( y1, zz,m,  &
           x1,  z,n,  opt)

!***********************************************************************
!     interpolation soubroutine basada en Numerical Recipes HUNT.FOR
!     input values: x1(n) at z(n)
!     output values: y1(m) at zz(m)
!     options: 1 -> lineal
!     2 -> logarithmic
!***********************************************************************

      implicit none

!     Arguments
      integer	:: n,m,opt         ! I
      real	:: zz(m),z(n)      ! I
      real    :: y1(m)             ! O
      real    :: x1(n)             ! I


!     Local variables
      integer :: i, j
      real    :: factor
      real    :: zaux

!!!!

      j = 1                     ! initial first guess (=n/2 is anothr pssblty)

      do 1,i=1,m                !

                                ! Busca indice j donde ocurre q zz(i) esta entre [z(j),z(j+1)]
         zaux = zz(i)
         if (abs(zaux-z(1)).le.0.01) then
            zaux=z(1)
         elseif (abs(zaux-z(n)).le.0.01) then
            zaux=z(n)
         endif
         call hunt ( z,n, zaux, j )
         if ( j.eq.0 .or. j.eq.n ) then
	    write (*,*) ' HUNT/ Limits input grid:', z(1),z(n)
	    write (*,*) ' HUNT/ location in new grid:', zz(i)
            stop ' interhunt/ Interpolat error. z out of limits.'
         endif

                                ! Perform interpolation
         factor = (zz(i)-z(j))/(z(j+1)-z(j))
         if (opt.eq.1) then
	    y1(i) = x1(j) + (x1(j+1)-x1(j)) * factor
         else
	    y1(i) = exp( log(x1(j)) + log(x1(j+1)/x1(j)) * factor )
         end if


 1    continue

      end subroutine interhunt


! **********************************************************************
!     *** interhuntlimits2veces.f ***
! **********************************************************************

      subroutine interhuntlimits2veces  &
           ( y1,y2, zz,m, limite1,limite2, &
           x1,x2,  z,n, opt)

! **********************************************************************
!     Interpolation soubroutine basada en Numerical Recipes HUNT.FOR
!     Input values:  x1,x2(n) at z(n)
!     Output values:
!     y1,y2(m) at zz(m)   pero solo entre los indices de zz
!     siguientes: [limite1,limite2]
!     Options: 1 -> linear in z and linear in x
!     2 -> linear in z and logarithmic in x
!     3 -> logarithmic in z and linear in x
!    4 -> logarithmic in z and logaritmic in x
!
!     NOTAS: Esta subrutina extiende y generaliza la usual
!     "interhunt5veces" en 2 direcciones:
!     - la condicion en los limites es que zz(limite1:limite2)
!     esté dentro de los limites de z (pero quizas no todo zz)
!     - se han añadido 3 opciones mas al caso de interpolacion
!     logaritmica, ahora se hace en log de z, de x o de ambos.
!     Notese que esta subrutina engloba a la interhunt5veces
!      ( esta es reproducible haciendo  limite1=1  y  limite2=m
!     y usando una de las 2 primeras opciones  opt=1,2 )
!***********************************************************************

      implicit none

!     Arguments
      integer	:: n,m,opt, limite1,limite2 ! I
      real	:: zz(m),z(n)      ! I
      real    :: y1(m),y2(m)       ! O
      real    :: x1(n),x2(n)       ! I


!     Local variables
      integer :: i, j
      real    :: factor
      real    :: zaux

!!!!

      j = 1                     ! initial first guess (=n/2 is anothr pssblty)

      do 1,i=limite1,limite2

                                ! Busca indice j donde ocurre q zz(i) esta entre [z(j),z(j+1)]
         zaux = zz(i)
         if (abs(zaux-z(1)).le.0.01) then
            zaux=z(1)
         elseif (abs(zaux-z(n)).le.0.01) then
            zaux=z(n)
         endif
         call hunt ( z,n, zaux, j )
         if ( j.eq.0 .or. j.eq.n ) then
	    write (*,*) ' HUNT/ Limits input grid:', z(1),z(n)
	    write (*,*) ' HUNT/ location in new grid:', zz(i)
            stop ' interhuntlimits/ Interpolat error. z out of limits.'
         endif

                                ! Perform interpolation
         if (opt.eq.1) then
	    factor = (zz(i)-z(j))/(z(j+1)-z(j))
	    y1(i) = x1(j) + (x1(j+1)-x1(j)) * factor
	    y2(i) = x2(j) + (x2(j+1)-x2(j)) * factor

         elseif (opt.eq.2) then
	    factor = (zz(i)-z(j))/(z(j+1)-z(j))
	    y1(i) = exp( log(x1(j)) + log(x1(j+1)/x1(j)) * factor )
	    y2(i) = exp( log(x2(j)) + log(x2(j+1)/x2(j)) * factor )
         elseif (opt.eq.3) then
	    factor = (log(zz(i))-log(z(j)))/(log(z(j+1))-log(z(j)))
	    y1(i) = x1(j) + (x1(j+1)-x1(j)) * factor
	    y2(i) = x2(j) + (x2(j+1)-x2(j)) * factor
         elseif (opt.eq.4) then
	    factor = (log(zz(i))-log(z(j)))/(log(z(j+1))-log(z(j)))
	    y1(i) = exp( log(x1(j)) + log(x1(j+1)/x1(j)) * factor )
	    y2(i) = exp( log(x2(j)) + log(x2(j+1)/x2(j)) * factor )
         end if


 1    continue

      end subroutine interhuntlimits2veces

! **********************************************************************
!     *** interhuntlimits5veces.f ***
! **********************************************************************

      subroutine interhuntlimits5veces  &
           ( y1,y2,y3,y4,y5, zz,m, limite1,limite2,  &
           x1,x2,x3,x4,x5,  z,n, opt)

! **********************************************************************
!     Interpolation soubroutine basada en Numerical Recipes HUNT.FOR
!     Input values:  x1,x2,..,x5(n) at z(n)
!     Output values:
!     y1,y2,...,y5(m) at zz(m)   pero solo entre los indices de zz
!     siguientes: [limite1,limite2]
!     Options: 1 -> linear in z and linear in x
!     2 -> linear in z and logarithmic in x
!     3 -> logarithmic in z and linear in x
!     4 -> logarithmic in z and logaritmic in x
!
!     NOTAS: Esta subrutina extiende y generaliza la usual
!     "interhunt5veces" en 2 direcciones:
!     - la condicion en los limites es que zz(limite1:limite2)
!     esté dentro de los limites de z (pero quizas no todo zz)
!     - se han añadido 3 opciones mas al caso de interpolacion
!     logaritmica, ahora se hace en log de z, de x o de ambos.
!     Notese que esta subrutina engloba a la interhunt5veces
!     ( esta es reproducible haciendo  limite1=1  y  limite2=m
!     y usando una de las 2 primeras opciones  opt=1,2 )
!***********************************************************************

      implicit none

!     Arguments
      integer	:: n,m,opt, limite1,limite2 ! I
      real	:: zz(m),z(n)      ! I
      real      :: y1(m),y2(m),y3(m),y4(m),y5(m) ! O
      real      :: x1(n),x2(n),x3(n),x4(n),x5(n) ! I


!     Local variables
      integer :: i, j
      real    :: factor
      real    :: zaux

!!!!

      j = 1                     ! initial first guess (=n/2 is anothr pssblty)

      do 1,i=limite1,limite2

                                ! Busca indice j donde ocurre q zz(i) esta entre [z(j),z(j+1)]
         zaux = zz(i)
         if (abs(zaux-z(1)).le.0.01) then
            zaux=z(1)
         elseif (abs(zaux-z(n)).le.0.01) then
            zaux=z(n)
         endif
         call hunt ( z,n, zaux, j )
         if ( j.eq.0 .or. j.eq.n ) then
	    write (*,*) ' HUNT/ Limits input grid:', z(1),z(n)
	    write (*,*) ' HUNT/ location in new grid:', zz(i)
            stop ' interhuntlimits/ Interpolat error. z out of limits.'
         endif

                                ! Perform interpolation
         if (opt.eq.1) then
	    factor = (zz(i)-z(j))/(z(j+1)-z(j))
	    y1(i) = x1(j) + (x1(j+1)-x1(j)) * factor
	    y2(i) = x2(j) + (x2(j+1)-x2(j)) * factor
	    y3(i) = x3(j) + (x3(j+1)-x3(j)) * factor
	    y4(i) = x4(j) + (x4(j+1)-x4(j)) * factor
	    y5(i) = x5(j) + (x5(j+1)-x5(j)) * factor
         elseif (opt.eq.2) then
	    factor = (zz(i)-z(j))/(z(j+1)-z(j))
	    y1(i) = exp( log(x1(j)) + log(x1(j+1)/x1(j)) * factor )
	    y2(i) = exp( log(x2(j)) + log(x2(j+1)/x2(j)) * factor )
	    y3(i) = exp( log(x3(j)) + log(x3(j+1)/x3(j)) * factor )
	    y4(i) = exp( log(x4(j)) + log(x4(j+1)/x4(j)) * factor )
	    y5(i) = exp( log(x5(j)) + log(x5(j+1)/x5(j)) * factor )
         elseif (opt.eq.3) then
	    factor = (log(zz(i))-log(z(j)))/(log(z(j+1))-log(z(j)))
	    y1(i) = x1(j) + (x1(j+1)-x1(j)) * factor
	    y2(i) = x2(j) + (x2(j+1)-x2(j)) * factor
	    y3(i) = x3(j) + (x3(j+1)-x3(j)) * factor
	    y4(i) = x4(j) + (x4(j+1)-x4(j)) * factor
	    y5(i) = x5(j) + (x5(j+1)-x5(j)) * factor
         elseif (opt.eq.4) then
	    factor = (log(zz(i))-log(z(j)))/(log(z(j+1))-log(z(j)))
	    y1(i) = exp( log(x1(j)) + log(x1(j+1)/x1(j)) * factor )
	    y2(i) = exp( log(x2(j)) + log(x2(j+1)/x2(j)) * factor )
	    y3(i) = exp( log(x3(j)) + log(x3(j+1)/x3(j)) * factor )
	    y4(i) = exp( log(x4(j)) + log(x4(j+1)/x4(j)) * factor )
	    y5(i) = exp( log(x5(j)) + log(x5(j+1)/x5(j)) * factor )
         end if


 1    continue

      end subroutine interhuntlimits5veces

! **********************************************************************
!     *** interhuntlimits.f ***
! **********************************************************************

      subroutine interhuntlimits ( y, zz,m, limite1,limite2,  &
           x,  z,n, opt)

! **********************************************************************
!     Interpolation soubroutine basada en Numerical Recipes HUNT.FOR
!     Input values:  x(n) at z(n)
!     Output values: y(m) at zz(m)   pero solo entre los indices de zz
!     siguientes: [limite1,limite2]
!     Options: 1 -> linear in z and linear in x
!     2 -> linear in z and logarithmic in x
!     3 -> logarithmic in z and linear in x
!     4 -> logarithmic in z and logaritmic in x
!
!     NOTAS: Esta subrutina extiende y generaliza la usual  "interhunt"
!     en 2 direcciones:
!     - la condicion en los limites es que zz(limite1:limite2)
!     esté dentro de los limites de z (pero quizas no todo zz)
!     - se han añadido 3 opciones mas al caso de interpolacion
!     logaritmica, ahora se hace en log de z, de x o de ambos.
!     Notese que esta subrutina engloba a la usual interhunt
!     ( esta es reproducible haciendo  limite1=1  y  limite2=m
!     y usando una de las 2 primeras opciones  opt=1,2 )
!***********************************************************************

      implicit none

!     Arguments
      integer	:: n,m,opt, limite1,limite2 ! I
      real	:: zz(m),z(n)      ! I
      real      :: y(m)              ! O
      real      :: x(n)              ! I


!     Local variables
      integer :: i, j
      real    :: factor
      real    :: zaux

!!!!

      j = 1                     ! initial first guess (=n/2 is anothr pssblty)

      do 1,i=limite1,limite2

                                ! Busca indice j donde ocurre q zz(i) esta entre [z(j),z(j+1)]
         zaux = zz(i)
         if (abs(zaux-z(1)).le.0.01) then
            zaux=z(1)
         elseif (abs(zaux-z(n)).le.0.01) then
            zaux=z(n)
         endif
         call hunt ( z,n, zaux, j )
         if ( j.eq.0 .or. j.eq.n ) then
	    write (*,*) ' HUNT/ Limits input grid:', z(1),z(n)
	    write (*,*) ' HUNT/ location in new grid:', zz(i)
            stop ' interhuntlimits/ Interpolat error. z out of limits.'
         endif

                                ! Perform interpolation
         if (opt.eq.1) then
	    factor = (zz(i)-z(j))/(z(j+1)-z(j))
	    y(i) = x(j) + (x(j+1)-x(j)) * factor
         elseif (opt.eq.2) then
	    factor = (zz(i)-z(j))/(z(j+1)-z(j))
	    y(i) = exp( log(x(j)) + log(x(j+1)/x(j)) * factor )
         elseif (opt.eq.3) then
	    factor = (log(zz(i))-log(z(j)))/(log(z(j+1))-log(z(j)))
	    y(i) = x(j) + (x(j+1)-x(j)) * factor
         elseif (opt.eq.4) then
	    factor = (log(zz(i))-log(z(j)))/(log(z(j+1))-log(z(j)))
	    y(i) = exp( log(x(j)) + log(x(j+1)/x(j)) * factor )
         end if


 1    continue

      end subroutine interhuntlimits

! ***************************************************************
!     *** lubksb_dp.f ***
! ***************************************************************

      subroutine lubksb_dp(a,n,np,indx,b)

      implicit none

      integer,intent(in) :: n,np
      real*8,intent(in) :: a(np,np)
      integer,intent(in) :: indx(n)
      real*8,intent(out) :: b(n)

      real*8 :: sum
      integer :: ii, ll, i, j

      ii=0
      do 12 i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0)then
            do 11 j=ii,i-1
               sum=sum-a(i,j)*b(j)
 11         continue
         else if (sum.ne.0.0) then
            ii=i
         endif
         b(i)=sum
 12   continue
      do 14 i=n,1,-1
         sum=b(i)
         if(i.lt.n)then
            do 13 j=i+1,n
               sum=sum-a(i,j)*b(j)
 13         continue
         endif
         b(i)=sum/a(i,i)
 14   continue
      end subroutine lubksb_dp

! ***************************************************
!     *** ludcmp_dp.f ***
! ***************************************************

      subroutine ludcmp_dp(a,n,np,indx,d)

      implicit none

      integer,intent(in) :: n, np
      real*8,intent(inout) :: a(np,np)
      real*8,intent(out) :: d
      integer,intent(out) :: indx(n)

      integer :: i, j, k, imax
      real*8, parameter :: tiny=1.0d-20
      integer,parameter :: nmax=100
      real*8 :: vv(nmax), aamax, sum, dum

      d=1.0d0
      do 12 i=1,n
         aamax=0.0d0
         do 11 j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
 11      continue
         if (aamax.eq.0.0) then
            write(*,*) 'ludcmp_dp: singular matrix!'
            stop
         endif
         vv(i)=1.0d0/aamax
 12   continue
      do 19 j=1,n
         if (j.gt.1) then
            do 14 i=1,j-1
               sum=a(i,j)
               if (i.gt.1)then
                  do 13 k=1,i-1
                     sum=sum-a(i,k)*a(k,j)
 13               continue
                  a(i,j)=sum
               endif
 14         continue
         endif
         aamax=0.0d0
         do 16 i=j,n
            sum=a(i,j)
            if (j.gt.1)then
               do 15 k=1,j-1
                  sum=sum-a(i,k)*a(k,j)
 15            continue
               a(i,j)=sum
            endif
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
 16      continue
         if (j.ne.imax)then
            do 17 k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
 17         continue
            d=-d
            vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(j.ne.n)then
            if(a(j,j).eq.0.0)a(j,j)=tiny
            dum=1.0d0/a(j,j)
            do 18 i=j+1,n
               a(i,j)=a(i,j)*dum
 18         continue
         endif
 19   continue
      if(a(n,n).eq.0.0)a(n,n)=tiny

      end subroutine ludcmp_dp


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     *** LUdec.f ***
!
!     Solution of linear equation without inverting matrix
!     using LU decomposition:
!     AA * xx = bb         AA, bb: known
!     xx: to be found
!     AA and bb are not modified in this subroutine
!
!     MALV , Sep 2007
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine LUdec(xx,aa,bb,m,n)

      implicit none

!     Arguments
      integer,intent(in) ::     m, n
      real*8,intent(in) ::      aa(m,m), bb(m)
      real*8,intent(out) ::     xx(m)


!     Local variables
      real*8      :: a(n,n), b(n), x(n), d
      integer     :: i, j, indx(n)


!     Subrutinas utilizadas
!     ludcmp_dp, lubksb_dp

!!!!!!!!!!!!!!!Comienza el programa !!!!!!!!!!!!!!

      do i=1,n
         b(i) = bb(i+1)
         do j=1,n
            a(i,j) = aa(i+1,j+1)
         enddo
      enddo

                                ! Descomposicion de auxm1
      call ludcmp_dp ( a, n, n, indx, d)

                                ! Sustituciones foward y backwards para hallar la solucion
      do i=1,n
         x(i) = b(i)
      enddo
      call lubksb_dp( a, n, n, indx, x )

      do i=1,n
         xx(i+1) = x(i)
      enddo

      end subroutine LUdec

! ***************************************************************************
!     *** mat_oper.f ***
! ***************************************************************************
      subroutine unit(a,n)
!     store the unit value in the diagonal of a
! ***************************************************************************
      implicit none
      real*8 :: a(n,n)
      integer :: n,i,j,k
      do 1,i=2,n-1
         do 2,j=2,n-1
	    if(i.eq.j) then
               a(i,j) = 1.d0
	    else
               a(i,j)=0.0d0
	    end if
 2       continue
 1    continue
      do k=1,n
         a(n,k) = 0.0d0
         a(1,k) = 0.0d0
         a(k,1) = 0.0d0
         a(k,n) = 0.0d0
      end do
      end subroutine unit

!     ***********************************************************************
      subroutine diago(a,v,n)
!     store the vector v in the diagonal elements of the square matrix a
!     ***********************************************************************
      implicit none

      integer :: n,i,j,k
      real*8 :: a(n,n),v(n)

      do 1,i=2,n-1
         do 2,j=2,n-1
	    if(i.eq.j) then
               a(i,j) = v(i)
	    else
               a(i,j)=0.0d0
	    end if
 2       continue
 1    continue
      do k=1,n
         a(n,k) = 0.0d0
         a(1,k) = 0.0d0
         a(k,1) = 0.0d0
         a(k,n) = 0.0d0
      end do

      end subroutine diago

!     ***********************************************************************
      subroutine invdiag(a,b,n)
!     inverse of a diagonal matrix
!     ***********************************************************************
      implicit none

      integer :: n,i,j,k
      real*8 :: a(n,n),b(n,n)

      do 1,i=2,n-1
         do 2,j=2,n-1
	    if (i.eq.j) then
               a(i,j) = 1.d0/b(i,i)
	    else
               a(i,j)=0.0d0
	    end if
 2       continue
 1    continue
      do k=1,n
         a(n,k) = 0.0d0
         a(1,k) = 0.0d0
         a(k,1) = 0.0d0
         a(k,n) = 0.0d0
      end do

      end subroutine invdiag

!     ***********************************************************************
      subroutine samem (a,m,n)

!     store the matrix m in the matrix a
!     ***********************************************************************

      implicit none

      real*8 :: a(n,n),m(n,n)
      integer :: n,i,j,k

      do 1,i=2,n-1
         do 2,j=2,n-1
            a(i,j) = m(i,j)
 2       continue
 1    continue
      do k=1,n
         a(n,k) = 0.0d0
         a(1,k) = 0.0d0
         a(k,1) = 0.0d0
         a(k,n) = 0.0d0
      end do

      end subroutine samem

!     ***********************************************************************
      subroutine mulmv(a,b,c,n)

!     do a(i)=b(i,j)*c(j). a, b, and c must be distint
!     ***********************************************************************

      implicit none

      integer :: n,i,j
      real*8 :: a(n),b(n,n),c(n),sum

      do 1,i=2,n-1
         sum=0.0d0
         do 2,j=2,n-1
	    sum = sum + b(i,j) * c(j)
 2       continue
         a(i)=sum
 1    continue
      a(1) = 0.0d0
      a(n) = 0.0d0

      end subroutine mulmv

!     ***********************************************************************
      subroutine trucodiag(a,b,c,d,e,n)

!     inputs: matrices b,c,d,e
!     output: matriz diagonal a
!     Operacion a realizar:  a = b * c^(-1) * d + e
!     La matriz c va a ser invertida
!     Todas las matrices de entrada son diagonales excepto b
!     Aprovechamos esa condicion para invertir c, acelerar el calculo, y
!     ademas, para forzar que a sea diagonal
!     ***********************************************************************

      implicit none

      real*8 :: a(n,n),b(n,n),c(n,n),d(n,n),e(n,n), sum
      integer :: n,i,j,k
      do 1,i=2,n-1
         sum=0.0d0
         do 2,j=2,n-1
	    sum=sum+ b(i,j) * d(j,j)/c(j,j)
 2       continue
         a(i,i) = sum + e(i,i)
 1    continue
      do k=1,n
         a(n,k) = 0.0d0
         a(1,k) = 0.0d0
         a(k,1) = 0.0d0
         a(k,n) = 0.0d0
      end do

      end subroutine trucodiag

!     ***********************************************************************
      subroutine trucommvv(v,b,c,u,w,n)

!     inputs: matrices b,c , vectores u,w
!     output: vector v
!     Operacion a realizar:  v = b * c^(-1) * u + w
!     La matriz c va a ser invertida
!     c es diagonal, b no
!     Aprovechamos esa condicion para invertir c, y acelerar el calculo
!     ***********************************************************************

      implicit none

      real*8 :: v(n),b(n,n),c(n,n),u(n),w(n), sum
      integer :: n,i,j
      do 1,i=2,n-1
         sum=0.0d0
         do 2,j=2,n-1
	    sum=sum+ b(i,j) * u(j)/c(j,j)
 2       continue
         v(i) = sum + w(i)
 1    continue
      v(1) = 0.d0
      v(n) = 0.d0

      end subroutine trucommvv

!     ***********************************************************************
      subroutine sypvmv(v,u,c,w,n)

!     inputs: matriz diagonal c , vectores u,w
!     output: vector v
!     Operacion a realizar:  v = u + c * w
!     ***********************************************************************

      implicit none

      real*8 :: v(n),u(n),c(n,n),w(n)
      integer :: n,i
      do 1,i=2,n-1
         v(i)= u(i) + c(i,i) * w(i)
 1    continue
      v(1) = 0.0d0
      v(n) = 0.0d0

      end subroutine sypvmv

!     ***********************************************************************
      subroutine sumvv(a,b,c,n)

!     a(i)=b(i)+c(i)
!     ***********************************************************************
      implicit none

      integer :: n,i
      real*8 :: a(n),b(n),c(n)

      do 1,i=2,n-1
         a(i)= b(i) + c(i)
 1    continue
      a(1) = 0.0d0
      a(n) = 0.0d0

      end subroutine sumvv

!     ***********************************************************************
      subroutine sypvvv(a,b,c,d,n)

!     a(i)=b(i)+c(i)*d(i)
!     ***********************************************************************
      implicit none

      real*8 :: a(n),b(n),c(n),d(n)
      integer :: n,i
      do 1,i=2,n-1
         a(i)= b(i) + c(i) * d(i)
 1    continue
      a(1) = 0.0d0
      a(n) = 0.0d0

      end subroutine sypvvv

!     ***********************************************************************
      subroutine zero4m(a,b,c,d,n)

!     a(i,j) = b(i,j) = c(i,j) = d(i,j) = 0.0
!     ***********************************************************************
      implicit none

      real*8 :: a(n,n), b(n,n), c(n,n), d(n,n)
      integer :: n
      a(1:n,1:n)=0.d0
      b(1:n,1:n)=0.d0
      c(1:n,1:n)=0.d0
      d(1:n,1:n)=0.d0
!      do 1,i=1,n
!         do 2,j=1,n
!	    a(i,j) = 0.0d0
!	    b(i,j) = 0.0d0
!	    c(i,j) = 0.0d0
!	    d(i,j) = 0.0d0
! 2       continue
! 1    continue

      end subroutine zero4m

!     ***********************************************************************
      subroutine zero3m(a,b,c,n)

!     a(i,j) = b(i,j) = c(i,j) = 0.0
!     **********************************************************************
      implicit none

      real*8 :: a(n,n), b(n,n), c(n,n)
      integer :: n
      a(1:n,1:n)=0.d0
      b(1:n,1:n)=0.d0
      c(1:n,1:n)=0.d0
!      do 1,i=1,n
!         do 2,j=1,n
!	    a(i,j) = 0.0d0
!	    b(i,j) = 0.0d0
!	    c(i,j) = 0.0d0
! 2       continue
! 1    continue

      end subroutine zero3m


!     ***********************************************************************
      subroutine zero2m(a,b,n)

!     a(i,j) = b(i,j) = 0.0
!     ***********************************************************************
      implicit none

      real*8 :: a(n,n), b(n,n)
      integer :: n
      a(1:n,1:n)=0.d0
      b(1:n,1:n)=0.d0
!      do 1,i=1,n
!         do 2,j=1,n
!	    a(i,j) = 0.0d0
!	    b(i,j) = 0.0d0
! 2       continue
! 1    continue

      end subroutine zero2m

!     ***********************************************************************
      subroutine zero4v(a,b,c,d,n)

!     a(i) = b(i) = c(i) = d(i,j) = 0.0
!     ***********************************************************************
      implicit none

      real*8 :: a(n), b(n), c(n), d(n)
      integer :: n
      a(1:n)=0.d0
      b(1:n)=0.d0
      c(1:n)=0.d0
      d(1:n)=0.d0
!      do 1,i=1,n
!         a(i) = 0.0d0
!         b(i) = 0.0d0
!         c(i) = 0.0d0
!         d(i) = 0.0d0
! 1    continue

      end subroutine zero4v


!     ***********************************************************************
      subroutine zero3v(a,b,c,n)

!     a(i) = b(i) = c(i) = 0.0
!     ***********************************************************************
      implicit none

      real*8 :: a(n), b(n), c(n)
      integer :: n
      a(1:n)=0.d0
      b(1:n)=0.d0
      c(1:n)=0.d0
!      do 1,i=1,n
!         a(i) = 0.0d0
!         b(i) = 0.0d0
!         c(i) = 0.0d0
! 1    continue

      end subroutine zero3v

!     ***********************************************************************
      subroutine zero2v(a,b,n)

!     a(i) = b(i) = 0.0
!     ***********************************************************************
      implicit none

      real*8 :: a(n), b(n)
      integer :: n
      a(1:n)=0.d0
      b(1:n)=0.d0
!      do 1,i=1,n
!         a(i) = 0.0d0
!         b(i) = 0.0d0
! 1    continue

      end subroutine zero2v


!****************************************************************************
!
!     *** suaviza.f ***
!
!*****************************************************************************
!
      subroutine suaviza ( x, n, ismooth, y )

!*****************************************************************************
!     x - input and return values
!     y - auxiliary vector
!     ismooth = 0  --> no smoothing is performed
!     ismooth = 1  --> weak smoothing (5 points, centred weighted)
!     ismooth = 2  --> normal smoothing (3 points, evenly weighted)
!     ismooth = 3  --> strong smoothing (5 points, evenly weighted)
!     august 1991
!*****************************************************************************

      implicit none

      integer	:: n, imax, imin, i, ismooth
      real*8	:: x(n), y(n)
!*****************************************************************************

      imin=1
      imax=n

      if (ismooth.eq.0) then

         return

      elseif (ismooth.eq.1) then ! 5 points, with central weighting

         do i=imin,imax
	    if(i.eq.imin)then
               y(i)=x(imin)
	    elseif(i.eq.imax)then
               y(i)=x(imax-1)+(x(imax-1)-x(imax-3))/2.d0
	    elseif(i.gt.(imin+1) .and. i.lt.(imax-1) )then
               y(i) = ( x(i+2)/4.d0 + x(i+1)/2.d0 + 2.d0*x(i)/3.d0 +  &
                    x(i-1)/2.d0 + x(i-2)/4.d0 )* 6.d0/13.d0
	    else
               y(i)=(x(i+1)/2.d0+x(i)+x(i-1)/2.d0)/2.d0
	    end if
         end do

      elseif (ismooth.eq.2) then ! 3 points, evenly spaced

         do i=imin,imax
	    if(i.eq.imin)then
               y(i)=x(imin)
	    elseif(i.eq.imax)then
               y(i)=x(imax-1)+(x(imax-1)-x(imax-3))/2.d0
	    else
               y(i) = ( x(i+1)+x(i)+x(i-1) )/3.d0
	    end if
         end do

      elseif (ismooth.eq.3) then ! 5 points, evenly spaced

         do i=imin,imax
	    if(i.eq.imin)then
               y(i) = x(imin)
	    elseif(i.eq.(imin+1) .or. i.eq.(imax-1))then
               y(i) = ( x(i+1)+x(i)+x(i-1) )/3.d0
	    elseif(i.eq.imax)then
               y(i) = ( x(imax-1) + x(imax-1) + x(imax-2) ) / 3.d0
	    else
               y(i) = ( x(i+2)+x(i+1)+x(i)+x(i-1)+x(i-2) )/5.d0
	    end if
         end do

      else

         write (*,*) ' Error in suaviza.f   Wrong ismooth value.'
         stop

      endif

!     rehago el cambio, para devolver x(i)
      do i=imin,imax
         x(i)=y(i)
      end do

      end subroutine suaviza

!     ***********************************************************************

      subroutine mulmmf90(a,b,c,n)

!     ***********************************************************************
      implicit none

      real*8 :: a(n,n), b(n,n), c(n,n)
      integer :: n

      a=matmul(b,c)
      a(1,:)=0.d0
      a(:,1)=0.d0
      a(n,:)=0.d0
      a(:,n)=0.d0

      end subroutine mulmmf90

!     ***********************************************************************

      subroutine resmmf90(a,b,c,n)

!     ***********************************************************************
      implicit none
      real*8 :: a(n,n), b(n,n), c(n,n)
      integer :: n

      a=b-c
      a(1,:)=0.d0
      a(:,1)=0.d0
      a(n,:)=0.d0
      a(:,n)=0.d0

      end subroutine resmmf90


!*******************************************************************

      subroutine gethist_03 (ihist)

!*******************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      integer  ::         ihist

!     local variables
      integer  ::         j, r

!     ***************

      nbox = nbox_stored(ihist)
      do j=1,mm_stored(ihist)
         thist(j) = thist_stored(ihist,j)
         do r=1,nbox_stored(ihist)
	    no(r) = no_stored(ihist,r)
            sk1(j,r) = sk1_stored(ihist,j,r)
            xls1(j,r) = xls1_stored(ihist,j,r)
            xld1(j,r) = xld1_stored(ihist,j,r)
         enddo
      enddo

      end subroutine gethist_03


!     *******************************************************************

      subroutine rhist_03 (ihist)

!     *******************************************************************

      use ModPlanet
!     use ModNLTE

!     include 'nlte_paramdef.h'
!     include 'nlte_commons.h'

      implicit none

!     arguments
      integer   ::         ihist

!     local variables
      integer  ::      j, r
      real*8   ::      xx

!     ***************

      open(unit=3,file=hisfile,status='old')

      read(3,*)
      read(3,*)
      read(3,*) mm_stored(ihist)
      read(3,*)
      read(3,*) nbox_stored(ihist)
      read(3,*)

      if ( nbox_stored(ihist) .gt. nbox_max ) then
         write (*,*) ' nbox too large in input file ', hisfile
         stop ' Check maximum number nbox_max in mz1d.par '
      endif

      do j=mm_stored(ihist),1,-1
         read(3,*) thist_stored(ihist,j)
         do r=1,nbox_stored(ihist)
	    read(3,*) no_stored(ihist,r),  &
                 sk1_stored(ihist,j,r),    &
                 xls1_stored(ihist,j,r),   &
                 xx, xld1_stored(ihist,j,r)
         enddo

      enddo

      close(unit=3)

      end subroutine rhist_03

! **************************************************************
