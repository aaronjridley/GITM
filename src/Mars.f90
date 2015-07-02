!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine fill_photo(photoion, photoabs, photodis)

  use ModPlanet
  use ModEUV

  implicit none

  real, intent(out) :: photoion(Num_WaveLengths_High, nIons-1)
  real, intent(out) :: photoabs(Num_WaveLengths_High, nSpeciesTotal)
  real, intent(out) :: photodis(Num_WaveLengths_High, nSpeciesTotal)


  integer :: i, iSpecies, iIon, NWH

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

  ! ---------------------------------------------------------------------
  !  Specific PhotoiDissociations (Absorption-Total_Ionization) (nSpecies)
  !  Need:  CO2, O2, and N2 dissociations (only)
  ! ---------------------------------------------------------------------

  photodis(1:NWH,iCO2_)  = PhotoAbs_CO2(1:NWH)-PhotoIon_CO2(1:NWH)
  photodis(1:NWH,iN2_)   = PhotoAbs_N2(1:NWH)-PhotoIon_N2(1:NWH)
  photodis(1:NWH,iO2_)   = PhotoAbs_O2(1:NWH)-PhotoIon_O2(1:NWH)

end subroutine fill_photo

!---------------------------------------------------------------------
! Initialize Heating Rates
!---------------------------------------------------------------------

subroutine init_heating_efficiency

  use ModEUV, only: HeatingEfficiency_CB, eHeatingEfficiency_CB

  implicit none

  HeatingEfficiency_CB  = 0.21
  eHeatingEfficiency_CB = 0.0

  call init_radcool

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
  use ModEUV, only:SunOrbitEccentricity,AveCosSza
  use ModTime

  implicit none

  integer, intent(in) :: iBlock
  integer :: nw, iLat, iLon, iAlt
  integer :: L_LAYERS, L_LEVELS, L_NLAYRAD,L_NLEVRAD

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

  if (useGravityWave)  call calc_GW(iBlock)
  
  !\ -------------------------------------------------------------------
  ! CO2 NLTE Cooling Formulation from Miguel Lopez-Valverde (2001)
  !/

  !/ Cooling ON

  call calc_radcooling(iBlock)

  !\
  ! RadCoolingRate is in K/s (from nltecool.F routine)
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
  RadCooling(1:nLons,1:nLats,1:nAlts,iBlock) = &
       RadCoolingRate(1:nLons,1:nLats,1:nAlts,iBlock)/&
       (TempUnit(1:nLons,1:nLats,1:nAlts)) 

  !/ Cooling OFF (zeroed out)
  !    RadCooling = 0.0
  !    RadCoolingRate = 0.0
  !\ -------------------------------------------------------------------


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

  ! SunOrbitEccentricity = sqrt(2.64236)
  do NW=1,L_NSPECTV
     SOL(nw) = SOLARF(NW)/(SunOrbitEccentricity**2)
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

     !	 print*, ell_s

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

           call calc_lowatmosrad(iblock,iLat,iLon,L_LAYERS,L_LEVELS,&
                L_NLAYRAD,L_NLEVRAD)

        end do !Latitude Loop
     end do    !Longitude Loop

!!$     xLowAtmosRadRate(1:nLons,1:nLats,1:nAlts,iBlock) = &
!!$          xLowAtmosRadRate(1:nLons,1:nLats,1:nAlts,iBlock)/&
!!$          (TempUnit(1:nLons,1:nLats,1:nAlts))

  endif

  !################# END MAIN COMPUTATIONAL LOOP#################

  call end_timing("calc_planet_sources")

end subroutine calc_planet_sources

!---------------------------------------------------------+
!
!---------------------------------------------------------+

!---------------------------------------------------------+
!
!---------------------------------------------------------+
!------------------------------------------------------------------
subroutine calc_radcooling(iBlock)

  !  ---------------------------------------------------------------------
  !  Purpose: to calculate NLTE CO2 15-micron cooling in Mars upper atmosphere
  !           (above 80 km, where NLTE physics becomes important)
  !  Source: nltecool.F from mtgcm16 
  !  Coders: M. Lopez-Valverde (2001)
  !          S. W. Bougher (2004-MTGCM15, 2007-Mars GITM)
  !  Inputs: from MarsGITM  (Pressure, Temperature, VMRCO2, VMRO, VMRN2CO)
  !          all on MarsGITM grid (nLons,nLats,nAlts)
  !  Inputs: from init_radcool (pnbr, ef1, ef2, co2vmr,o2pvmr,n2covmr)
  !          all on np=68 1-D NLTE code grid
  !  Output: COOLTOT (K/sec) and RadCoolingRate(1:nLons,1:nLats,1:nAlts,iBlock)
  !  =======================================================================
  !
  !     subroutine nltecool(nlayer,player,tlayer,dt)
  !
  ! This code was designed as a delivery for the "Martian Environment Models" 
  ! project ( ESA contract 11369/95/nl/jg CCN2 )
  ! Computes non-LTE heating rates from CO2 emission at 15 um
  ! in the Martian upper atmosphere.
  ! Uses a simplified model consisting of two excited levels with two
  ! emission bands, one of them stronger than the other, which correspond
  ! to the behaviours of the 626 fundamental band and the isotopic fund.bands.
  ! It uses a cool-to-space approximation with tabulated escape functions.
  ! These escape functions have been precomputed for the strong and weak bands,
  ! and are given as a function of pressure in separate files.
  ! The output values are the heating rates (actually, cooling, since they
  ! are always negative) for the two bands, i.e., the total cooling is the
  ! sum of them.
  ! Miguel A. Lopez Valverde
  ! Instituto de Astrofisica de Andalucia (CSIC), Granada, Spain
  !
  ! Version 1b.  See description above.  22-March-2000.
  ! Adapted as a subroutine for use in GCM -- PLR/SRL 6/2000
  ! Version 1c.  Inclusion of VMR in the tabulation of escape functions. 
  !              Table contains now only 1 input file -- Miguel 11/Jul/2000
  ! Nov-2001     Adapt to the MZ1D model (1D version & extended parametp_table)
  ! Nov-2005     Adapt to the MTGCM16 model (3D version & extended parametp_table)
  ! Jan-2007     Adapt to the MarsGITM model (3D version & extended parametp_table)
  !  -----------------------------------------------------------------

  use ModInputs
  use ModSources, only: RadCoolingRate
  use ModPlanet
  use ModGITM
  use ModConstants, only:  Boltzmanns_Constant, Speed_Light, &
       Planck_Constant
  use ModIndicesInterfaces

  implicit none

  ! INPUTS:----------------------------------------------------------

  integer,intent(in):: iBlock 

  !  -----------------------------------------------------------------
  ! Input variables

  integer, parameter :: nlayer = nAlts ! no. of atmospheric layers
  integer, parameter :: nlon = nLons ! no. of long gridpoints
  integer, parameter :: nlat = nLats ! no. of latitude gridpoints
  real :: player(nlayer)            ! input pressure grid [Pa]
  real :: tlayer(nlayer)            ! input temperatures [K]
  logical :: found
  
  !  -----------------------------------------------------------------

  ! Standard atmosphere variables (overlay with MarsGITM values)

  real :: nt , co2(nlayer) , o3p(nlayer) , n2co(nlayer)

  ! Vectors/indexes for the datavol2 tabulation of escape functions and VMR

  !     Number data points in Tabulation (in ModMars module)
  !     integer,parameter :: np=68 
  !     Dimensioned above from ModMars module
  !     real ::  pnbr(np), ef1(np), ef2(np), co2vmr(np), o3pvmr(np), n2covmr(np)
  !     Reference Pressure (Pascals) from 1-D NLTE code
  real ::  pnb(np) 
  !     Interpolated escape functions
  real ::  escf1(nlayer) , escf2(nlayer) 

  ! Local Constants and variables

  real  ::  n1, n2, co2t , l1, p1, p12 , l2, p2, p21
  real  ::  tt, c1, c2, ae1, ae2,  a1, a2, a12, a21
  real  ::  pl1, pl2, el1, el2, hr1, hr2, x 
  real  ::  hr(nlayer)  ,peakval

  ! Indexes

  integer  :: i,j,ii,ipeak
  integer :: iLat, iLon, iAlt

  ! Rate coefficients

  real  ::  k19xca, k19xcb, k19cap1, k19cap2, k19cbp1, k19cbp2
  real  ::  d19c, d19cp1, d19cp2, k20xc, k20cp1, k20cp2
  real  ::  k21xc, k21cp2
  !
  ! Data and settings
  !
  real,parameter :: nu1=667.38, nu2=662.3734, hplanck=6.6261e-27,  &
       gamma1=1.191e-5, vlight=3.e10, ee=1.438769
  real,parameter :: imr1=0.987, imr2=0.00408 + 0.0112,  &
       rfvt=0.1, rfvto3p=1.0, rfvv=0.1

  !  -----------------------------------------------------------------
  ! Internal fields recast from MarsGITM

  real,dimension(1:nLons,1:nLats,1:nAlts) ::    &
       T, P, vmrco2, vmro, vmrn2co, cooltot, mnd 

  !     Mars GITM real temperature (on its grid)  ----------------------

  T(1:nLons,1:nLats,1:nAlts) = &
       Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*  &
       TempUnit(1:nLons,1:nLats,1:nAlts)

  !     Mars GITM pressure in pascals  (on its grid) -------------------

  P(1:nLons,1:nLats,1:nAlts) = &
       Pressure(1:nLons,1:nLats,1:nAlts,iBlock)

  !    Volume Mixing Ratio Calculations

  mnd(1:nLons,1:nLats,1:nAlts) = &  
       NDensity(1:nLons,1:nLats,1:nAlts,iBlock)+1.0

  vmro(1:nLons,1:nLats,1:nAlts)  = & 
       NdensityS(1:nLons,1:nLats,1:nAlts,iO_,iBlock)/ &
       mnd(1:nLons,1:nLats,1:nAlts)
  vmrn2co(1:nLons,1:nLats,1:nAlts) = & 
       (NdensityS(1:nLons,1:nLats,1:nAlts,iCO_,iBlock)  &
       + NdensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock))/ &
       mnd(1:nLons,1:nLats,1:nAlts) 
  vmrco2(1:nLons,1:nLats,1:nAlts) = &  
       NdensityS(1:nLons,1:nLats,1:nAlts,iCO2_,iBlock) &
       /mnd(1:nLons,1:nLats,1:nAlts)

  !  -----------------------------------------------------------------
  !
  ! **** ARRAYS made available via init_radcool 
  ! **** Convert pnb array to pascal units here
  !
  do i=1,np
     pnb(i)=1.0e-4*exp(pnbr(i)) ! p into Pa
  enddo

  ! ****  Initialize to zero
  cooltot(1:nLons,1:nLats,1:nAlts) = 0.0 

  !-------------------------------------------------------------

  ! MAIN LOOP, for MarsGITM horizontal grid (Latitude,Longitude) 

  do iLat = 1, nlat
     do iLon = 1, nlon

        ! *****  populate 1-D arrays(nlayer) for calc_radcooling from MarsGITM 
        do  iAlt=1,nlayer
           player(iAlt)= P(iLon,iLat,iAlt)
           tlayer(iAlt)= T(iLon,iLat,iAlt)
           co2(iAlt)=vmrco2(iLon,iLat,iAlt)
           o3p(iAlt)=vmro(iLon,iLat,iAlt)
           n2co(iAlt)=vmrn2co(iLon,iLat,iAlt)
        enddo

        !        interpolate escape functions (only) to nlayer grid
        call interp1(escf2,player,nlayer,ef2,pnb,np)
        call interp1(escf1,player,nlayer,ef1,pnb,np)

        ! ALTITUDE LOOP, for each altitude (for 1-latitude/longitude profile):
        ! (from old nltecool.F routine)

        do i=1,nlayer  

           if (player(i).gt.PLONG .or. player(i).lt.4.0e-6) then 

              hr(i)=0.0
              cooltot(iLon,iLat,i)= 0

           else

              nt = player(i)/(1.381e-17*tlayer(i)) 
              co2(i)=co2(i)*nt                
              o3p(i)=o3p(i)*nt                
              n2co(i)=n2co(i)*nt
              n1 = co2(i) * imr1
              n2 = co2(i) * imr2
              co2t = n1 + n2

              tt = tlayer(i)*tlayer(i)
              k19xca = 4.2e-12 * exp( -2988.0/tlayer(i) + 303930.0/tt )
              k19xcb = 2.1e-12 * exp( -2659.0/tlayer(i) + 223052.0/tt )
              if (tlayer(i).le.175.0) then
                 k19xca = 3.3e-15
                 k19xcb = 7.6e-16
              endif
              k19xca = k19xca * rfvt
              k19xcb = k19xcb * rfvt
              k19cap1 = k19xca * 2.0 * exp( -ee*nu1/tlayer(i) )
              k19cap2 = k19xca * 2.0 * exp( -ee*nu2/tlayer(i) )
              k19cbp1 = k19xcb * 2.0 * exp( -ee*nu1/tlayer(i) )
              k19cbp2 = k19xcb * 2.0 * exp( -ee*nu2/tlayer(i) )
              d19c = k19xca*co2t + k19xcb*n2co(i)
              d19cp1 = k19cap1*co2t + k19cbp1*n2co(i)
              d19cp2 = k19cap2*co2t + k19cbp2*n2co(i)
              !
              k20xc = 3.e-12 * rfvto3p
              k20cp1 = k20xc * 2.0 * exp( -ee/tlayer(i) * nu1 )
              k20cp2 = k20xc * 2.0 * exp( -ee/tlayer(i) * nu2 )
              !
              k21xc = 2.49e-11 * 0.5 * rfvv
              k21cp2 = k21xc * exp( - ee/tlayer(i) * (nu2-nu1) )
              !
              l1 = d19c + k20xc*o3p(i) + k21cp2*n2
              p1 = ( d19cp1 + k20cp1*o3p(i) ) * n1
              p12 = k21xc*n1
              !
              l2 = d19c + k20xc*o3p(i) + k21xc*n1
              p2 = ( d19cp2 + k20cp2*o3p(i) ) * n2
              p21 = k21cp2*n2

              ae1 = 1.3546 * 1.66 / 4.0 * escf1(i)
              ae2 = ( 1.3452 + 1.1878 ) * 1.66 / 4.0 * escf2(i)
              l1 = l1 + ae1
              l2 = l2 + ae2

              c1 = gamma1*nu1**3. * 0.5
              c2 = gamma1*nu2**3. * 0.5
              a1 = c1 * p1 / (n1*l1)
              a2 = c2 * p2 / (n2*l2)
              a12 = (nu1/nu2)**3. * n2/n1 * p12/l1
              a21 = (nu2/nu1)**3. * n1/n2 * p21/l2
              el2 = (a2 + a21 * a1 ) / ( 1.0 - a21 * a12 )
              el1 = a1 + a12 * el2
              pl1 = el1 * n1 / c1
              pl2 = el2 * n2 / c2

              hr1 = - hplanck*vlight * nu1 * ae1 * pl1
              hr2 = - hplanck*vlight * nu2 * ae2 * pl2
              hr(i) = hr1 + hr2
              !        Cooling Rate (K/sec units)
              cooltot(iLon,iLat,i)= 0.1*hr(i)*tlayer(i)/(4.4*player(i))

           endif
        
        enddo  !-------- END OF MAIN ALTITUDE LOOP  
     enddo  !-------- END OF MAIN LONGITUDE LOOP  
  enddo  !-------- END OF MAIN LATITUDE LOOP  

  !-------------------------------------------------------------

  RadCoolingRate(1:nLons,1:nLats,1:nAlts,iBlock) = &
       -cooltot(1:nLons,1:nLats,1:nAlts) 

  !-------------------------------------------------------------

end subroutine calc_radcooling
!---------------------------------------------------------+

subroutine interp1(escout,p,nlayer,escin,pin,nl)
  implicit none

  ! subroutine to perform linear interpolation in pressure from 1D profile 
  ! escin(nl) sampled on pressure grid pin(nl) to profile
  ! escout(nlayer) on pressure grid p(nlayer).
  !
  ! Input args:
  !     integer,intent(in):: iBlock 

  integer,intent(in) :: nlayer,nl
  real,dimension(nl),intent(in) :: escin,pin
  real,dimension(nlayer),intent(in) :: p
  real,dimension(nlayer),intent(out) :: escout
  !
  ! Local :
  integer :: n,nm,np,n1
  real :: wm,wp
  !

  do n1=1,nlayer

     if(p(n1) .gt. 16.0 .or. p(n1) .lt. 5.0e-6) then
        escout(n1) = 0.0
     else

        do n = 1,nl-1

           if (p(n1).le.pin(n).and.p(n1).ge.pin(n+1)) then

              nm=n
              np=n+1
              wm=abs(pin(np)-p(n1))/(pin(nm)-pin(np))
              wp=1.0 - wm
           endif


        enddo


        escout(n1) = escin(nm)*wm + escin(np)*wp
     endif
  enddo

end subroutine interp1
!---------------------------------------------------------+

subroutine init_radcool
  !
  ! ---------------------------------------------------------
  ! Call to init_radcool will return corresponding set of six 1-D arrays:
  !    pnb(np),ef1(np),ef2(np),co2vmr(np),o3pvmr(np), n2covmr(np)
  !    These are 1-D NLTE cooling code inputs from Miguel Lopez-Valverde
  ! ---------------------------------------------------------

  use ModPlanet
  use ModGITM

  !\
  ! Below we prescribe the 6-arrays needed
  !/

  pnbr = (/ 12.0000, 11.0000, 10.8000, 10.6000, 10.4000, 10.2000, &
       10.0000, 9.80000, 9.60000, 9.40000, 9.20000, 9.00000, 8.80000,   & 
       8.60000, 8.40000, 8.20000, 8.00000, 7.80000, 7.60000, 7.40000,   &
       7.20000, 7.00000, 6.80000, 6.60000, 6.40000, 6.20000, 6.00000,   &
       5.80000, 5.60000, 5.40000, 5.20000, 5.00000, 4.80000, 4.60000,   &
       4.40000, 4.20000, 4.00000, 3.80000, 3.60000, 3.40000, 3.20000,   &
       3.00000, 2.80000, 2.60000, 2.40000, 2.20000, 2.00000, 1.80000,   &
       1.60000, 1.40000, 1.20000, 1.00000, 0.80000, 0.59999, 0.40000,   &
       0.20000, 0.00000,-0.20000,-0.40000,-0.60000,-0.80000,-1.00000,   &
       -1.20000,-1.40000,-1.60000,-1.80000,-2.00000,-3.00000/)

  ef1 = (/  0.000458112, 0.000458112,                              &
       0.000458112, 0.000458112, 0.000458112, 0.000458112, 0.000458112,  &
       0.000458112, 0.000458112, 0.000458112, 0.000458112, 0.000458112,  &
       0.000458112, 0.000458112, 0.000458112, 0.000458112, 0.000458112,  &
       0.000458112, 0.000458112, 0.000458112, 0.000458112, 0.000458112,  &
       0.000461707, 0.000476886, 0.000495638, 0.000520935, 0.000555511,  &
       0.000601219, 0.000663734, 0.000750691, 0.000863474, 0.001009000,  &
       0.001196420, 0.001426900, 0.001713980, 0.002066630, 0.002489740,  &
       0.003015780, 0.003643500, 0.004403230, 0.005320660, 0.006404560,  &
       0.007720690, 0.009256840, 0.011090500, 0.013237400, 0.015764300,  &
       0.018738800, 0.022207200, 0.026309900, 0.031061400, 0.036694800,  &
       0.043237300, 0.051502200, 0.062145500, 0.077721200, 0.099202700,  &
       0.131155000, 0.179470000, 0.258913000, 0.380549000, 0.530450000,  &
       0.643180000, 0.741061000, 0.826336000, 0.922787000, 0.997203000,  &
       1.000000000/)

  ef2 = (/  0.00198992,                                            &
       0.001989920, 0.001989920, 0.001989920, 0.001989920, 0.001989920,  &
       0.001989920, 0.001989920, 0.001989920, 0.001989920, 0.001989920,  &
       0.002013760, 0.002094500, 0.002229930, 0.002420560, 0.002680180,  &
       0.003043980, 0.003438960, 0.003802820, 0.004206220, 0.004761210,  &
       0.008016980, 0.011994700, 0.016914900, 0.022449700, 0.028524400,  &
       0.035481300, 0.043926400, 0.054624800, 0.067536700, 0.082993100,  &
       0.101717000, 0.123422000, 0.148468000, 0.177096000, 0.208816000,  &
       0.244003000, 0.282013000, 0.322559000, 0.365542000, 0.410518000,  &
       0.457384000, 0.505358000, 0.553627000, 0.600472000, 0.644807000,  &
       0.687185000, 0.727429000, 0.764734000, 0.798562000, 0.828699000,  &
       0.854797000, 0.877717000, 0.897874000, 0.915258000, 0.929904000,  &
       0.942381000, 0.952906000, 0.962173000, 0.970191000, 0.976437000,  &
       0.981501000, 0.985406000, 0.988560000, 0.991111000, 0.993653000,  &
       0.995561000, 1.000000000/)

  co2vmr =  (/                                                 &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.950000, 0.950000, 0.950000, 0.950000, 0.950000, 0.950000,  &
       0.949619, 0.947694, 0.945830, 0.944016, 0.940557, 0.937068,  &
       0.932366, 0.893661/)

  o3pvmr =  (/ 0.000000000, 0.000000000, 0.000000000, 0.000000000,   &
       0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000,   &
       0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000,   &
       0.000100148, 0.000118846, 0.000139681, 0.000164909, 0.000193617,   &
       0.000225161, 0.000260834, 0.000301501, 0.000344953, 0.000391011,   &
       0.000440377, 0.000490820, 0.000543200, 0.000595335, 0.000645420,   &
       0.000693166, 0.000743729, 0.000793710, 0.000844394, 0.000894318,   &
       0.000944732, 0.000994964, 0.001049010, 0.001100080, 0.001163020,   &
       0.001229890, 0.001300260, 0.001371310, 0.001455560, 0.001551860,   &
       0.001663280, 0.001778020, 0.001915460, 0.002075030, 0.002249030,   &
       0.002471170, 0.002717280, 0.002997390, 0.003335820, 0.003735070,   &
       0.004208190, 0.004768870, 0.005425580, 0.006208150, 0.007144730,   &
       0.008285450, 0.009517790, 0.010814000, 0.012235900, 0.013687000,   &
       0.015149500, 0.016719600, 0.018548500, 0.033625200/)

  n2covmr =  (/  0.027145, 0.027145,                                &
       0.0271450, 0.0271450, 0.0271450, 0.0271450, 0.0271450, 0.0271450, &
       0.0271450, 0.0271450, 0.0271450, 0.0271450, 0.0271450, 0.0271450, &
       0.0272661, 0.0272848, 0.0273054, 0.0273279, 0.0273514, 0.0273775, &
       0.0274048, 0.0274345, 0.0274672, 0.0275021, 0.0275404, 0.0275826, &
       0.0276340, 0.0277013, 0.0278220, 0.0279707, 0.0281759, 0.0284339, &
       0.0287587, 0.0291600, 0.0296561, 0.0302558, 0.0309922, 0.0318062, &
       0.0327010, 0.0335635, 0.0344388, 0.0353327, 0.0362143, 0.0370941, &
       0.0379315, 0.0387626, 0.0395524, 0.0403071, 0.0410071, 0.0416229, &
       0.0421231, 0.0425167, 0.0427964, 0.0429773, 0.0430488, 0.0429638, &
       0.0428049, 0.0426788, 0.0426822, 0.0429426, 0.0434634, 0.0442559, &
       0.0453038, 0.0465879, 0.0480262, 0.0496303, 0.0514885,            &
       0.0691651/)
  !
end subroutine init_radcool
! ----------------------------------------------------------------------

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
     KMax = 1000.0
     KMin = 100.0

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

              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) =  &
                   KMax* sqrt( NEddyMax(iLon,iLat) / NDensity(iLon,iLat,iAlt,iBlock)) 
              !
              !! \
              !! This gives an upper bound of Kmax
              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = &
                   min(KMax, KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) )
              !
              !! This gives an lower bound of Kmin
              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = &
                   max(KMin, KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) )
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

     real :: u(nlay),v(nlay),t(nlay),rhogw(nlay),rholev(nlay), &
          geot(nlay),pz(nlay),du(nlay),uav(nlay),plev(nlay), &
          tavg(nlay),theta(nlay),thtavg(nlay), &
          dtheta(nlay),stress(0:nlay),dstress(nlay), &
          dudt(nlay),cpco2

     real :: N(nlay),dh(nlay),dhsat(nlay),H(nlay), &
          dz(nlay),Ri(nlay),Rimin(nlay),eps(nlay)
     real :: Nsfc,k,C, VAR,TS,xRi,rhosfc,psfc,RCO2,gravv,gravsfc

     integer :: ilon,ilat,ialt

     PSFC = 400.0
     RCO2 = 189.0
     VAR  = 400.0**2
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
           T = Temperature(ilon,ilat,1:nAlts,iBlock)*TempUnit(ilon,ilat,1:nAlts)
           pz = rhogw * RCO2 * t
           THETA = t * (psfc/pz)**(RCO2/CpCO2)

           eps = 0.0
           dhsat = 0.0
           dh = 0.0
           Ri=0.0
           Rimin = 0.0
           stress = 0.0
           dstress = 0.0
           dudt = 0.0

           do nl=1,nlay-1
              !    Calculate the geometric distance between adjacent layer midpoints;
              !   units are  meters

              dz(nl) = geot(nl+1) - geot(nl)

              !    Calculate the average wind speed at the top of each model layer
              !   (but not the top layer), using a simple averaging of the zonal wind
              !   speeds at the two neighboring layer midpoints

              uav(nl) = 0.5 * (u(nl)+u(nl+1))

              !    Calculate the potential temperature at layer boundaries, via a simple
              !   averaging of the potential temperatures at the neighboring layer midpoints
              !    This can be imrpoved by using the layer boundary pressure (known
              !   via  psfc  and  sigma) to directly calculate the pressure term
              !   in the potential temperature calculation.. need a temperature at
              !   the layer boundary

              thtavg(nl) = 0.5 * (theta(nl)+theta(nl+1))

              !    Calculate the average temperature at layer boundaries, via a simple
              !   averaging of the temperatures at the neighboring layer midpoints

              tavg(nl) = 0.5 * (t(nl)+t(nl+1))

              !   Calculate the difference in zonal wind speed across a layer boundary

              du(nl) =   u(nl+1) - u(nl)
              !
              !   Calculate the density values at layer boundaries using the average of
              !  the pressures at the two neighboring layer midpoints and the average
              !  temperature at that same layer boundary;  this is a poor way of calculating 
              !  this average density and will be improved upon using independently
              !  determined values of pressure and temperature at layer boundaries

              plev(nl) = 0.5*(pz(nl)+pz(nl+1))/2.0
              rholev(nl) = 0.5*(pz(nl)+pz(nl+1))/(RCO2*Tavg(nl))

              !  Difference in potential temperature values across layer boundaries,
              ! calcultes by calculating the difference between the two neighboring
              ! layer midpoint potential temperatures;  this should be  (nl-1)-(nl)

              dtheta(nl) = theta(nl+1) - theta(nl)

           enddo



           !      Calculate the surface drag magnitude, using lowest layer midpoint
           !     density and wind speed values, and N (Brunt-Vaisala frequency) value
           !     at the interface between the bottom two atmosphere layers

           rhosfc = psfc / (RCO2 * T(1))

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


           do  22 nl=1,nlay-1

              gravv = -Gravity_GB(ilon,ilat,nl,iBlock)

              !    DO NOT IMPLEMENT WAVE PROPAGATION (more properly, 'BREAKING')
              !    WITHIN THE BOTTOM ONE MODEL LAYER(S)

              IF(NL.LE.1) THEN
                 STRESS(NL) = STRESS(NL-1)
                 DSTRESS(NL) = 0.0

                 !   The  ELSE below accounts for model layers above the bottom three
                 !  layers... in these upper layers wave breaking is permitted
              ELSE

                 !  Determine if a CRITICAL LEVEL situation is occurring, and if it is
                 ! deposit all of the stress within the offending layer (the layer within
                 ! which the sign of the zonal wind speed is opposite the sign of the
                 ! stress magnitude at the bottom of that layer)

                 !  DIAGNOSE CRITICAL LEVEL OCCURRENCE BASED UPON EITHER LAYER MIDPOINT
                 !  ZONAL WIND (u(nl)) OR AVERAGE ZONAL WIND AT THE TOP OF THE LAYER
                 !  (uav(nl)) BEING OPPOSITE IN SIGN TO THE STRESS VALUE AT THE BOTTOM
                 !  OF THE LAYER

                 if(u(nl)*stress(nl-1) .lt. 0.0 .or. &
                      uav(nl)*stress(nl-1) .lt. 0.0) then
                    !     stress at the layer top is set equal to zero
                    stress(nl) = 0.0
                    !     dstress at the layer top is set equal to the stress value
                    !   at the bottom of the layer

                    dstress(nl) = stress(nl-1)

                    !   for this layer within which a CRITICAL LEVEL has arisen and
                    !   wave momentum has been deposited, calculate the wind acceleration
                    !   (units of m/s/s) by dividing the  dstress value by the mass (kg) 
                    !   in the layer; a NEGATIVE sign indicates a westward acceleration,
                    !   and a positive value indicates an eastward accaleration

                    dudt(nl)=-dstress(nl)/((plev(nl)-plev(nl+1))/gravv)
                    GWAccel(ilon,ilat,nl,iEast_) = dudt(nl) * Dt

                    !       set dh(nl) and dhsat(nl) values to -9.99 in this CRITICAL LAYER
                    !     occurrence to easily help identify such layers
                    dh(nl) = -9.99
                    dhsat(nl) = -9.99
                    !     the  goto  statement below sends the code to the next layer,
                    !   but in actuality there is no further need to continue upwards 
                    !   in this column since there is no wave energy to consider...
                    !   THIS CAN BE IMPROVED AND MADE MORE EFFICIENT

                    !   WRITE to unit 17 to permit statistical look at wave breaking level

!                    write(17,*) msol,mbin,nl,dudt(nl),1.25+(nl-1)*2.5
                    !    below it would be more appropriate to  goto  statement  6  since
                    !    at this point this column is completed
                    goto 22
                 endif
                    !    the above  endif  ends the treatment of CRITICAL LEVELS

                    !  Calculate the Brunt-Vaisala frequency at the top boundary of
                    ! layer  nl

                    N(nl) = (gravv/thtavg(nl)) * (theta(nl+1)-theta(nl))/ &
                         dz(nl)

                    !  THERE IS CURRENTLY NO 'CHECK' FOR NEGATIVE N(nl) VALUES, WHICH
                    !  would arise in the presence of a superadiabatic lapse rate.. which
                    !  will not occur in teh Ames model output but could, I believe, arise
                    !  the no-hydrostatic GITM..

                    N(nl) = sqrt(N(nl))
!                    write(88,*) nl,N(nl),theta(nl+1),theta(nl),dz(nl)

                    !  Estimate the isentropic vertical displacement at the top of layer  nl
                    ! using the calculated stress value at the value at the lnect lower layer
                    ! boundary and the density, Brunt-Vaisala frequency and average zonal wind
                    ! speed at the layer boundary of interest

                    !    use ABS(stress(nl+1)) and ABS(uav(nl))

                    dh(nl) = sqrt(abs(stress(nl-1)/(rhogw(nl)*kap*N(nl)*uav(nl))))

!                    print * , '**',stress(nl-1),rhogw(nl),kap,n(nl),uav(nl),nl

                    !  Calculate the Richardson Number at the top boundary of Layer  nl

                    !   MAKE Ri  and Rimin  vectors with NLAY elements

                    Ri(nl)= (gravv/thtavg(nl))* (dtheta(nl)/dz(nl))/ &
                         ((du(nl)/dz(nl))**2)

!                    write(33,*) nl, Ri(nl),gravv,thtavg(nl),dtheta(nl), &
!                         dz(nl),du(nl),dz(nl)

                    !  Calculate the 'minimum' Richardson Number at the top of Layer  nl
                    ! using the Brunt-Vaisala frequency, estimated isentropic vertical displacement,
                    ! average zonal wind speed, and Richardson Number already calculated at that 
                    ! layer boundary of interest

                    !    implement  ABS(uav(nl))
                    Rimin(nl) = Ri(nl) * (1.0-(N(nl)*dh(nl)/ABS(uav(nl)))) / &
                         (1.0 + ((Ri(nl)**0.5)*N(nl)*dh(nl)/abs(uav(nl))))**2.0

                    !  If Rimin is greater than 0.25, than there is no wave breaking within
                    ! layer  nl  and thus the stress value at the top of Layer  nl  is the same
                    ! as the stress value at the bottom of Layre  nl  and the dstress value in
                    ! Layer  nl  is zero, and there is no wind acceleration within that layer

                    if(Rimin(nl).gt.0.25) then
                       stress(nl) = stress(nl-1)
                       dstress(nl) = 0.0
                    else

                       !   if Rimin is less than 0.25, than wave breaking is occurring and some stress
                       !  (acceleration.. really deceleration for westerly flow) is being applied to
                       !  the zonal wind in Layer  nl;

                       !  with the below inclusion of xRi I am NOT permitting
                       !  the value of  Ri  to be less than 0.25 when I calculate
                       !  the value of  epsILON... I am not sure if this is the
                       !  proper way to be proceeding here

                       xRi = max (Ri(nl),0.25000)

                       eps(nl)=(xRi**(-0.5))*(1.0+2.0*xRi**0.5)* &
                            ((2*(xRi**0.250)*(1.0+2.0*(xRi**0.50))**(-0.50))-1.0)
!                       write(15,99) nl,xRi,Ri(nl),Rimin(nl),eps(nl),2*(xRi**0.250), &
!                            (1.0+2.0*(xRi**0.50))**(-0.50),  &
!                            ((2*(xRi**0.250)*(1.0+2.0*(xRi**0.50))**(-0.50))-1.0)

                       !  eps  will not be less than zero, since if it did have a 
                       !  negative value it would result in a negative  dh  value

                       eps(nl) = max(eps(nl),0.0)

                       !  Calculate a more representative value of the isentropic displacement
                       ! at the top of Layer  nl, and then use that value to calculate a new value
                       ! of the stress at the top of Layer  nl, and then a value (dstress) that
                       ! is the change in the stress across the layer.  it is the magnitude of this
                       ! dstress value that determines the magnitude of the wind acceleration that
                       ! results

                       !  implement  ABS(uav(nl))

                       dhsat(nl)= eps(nl) / (N(nl)/ABS(uav(nl)))

                       !   Determine the stress value at a layer's top boundary if wave breaking
                       !  is diagnosed within that layer

                       stress(nl) = eps(nl) * eps(nl) * kap * rhogw(nl)*(uav(nl)**3)/ &
                            N(nl)

                       !  If the stress at the top of layer  nl  is greater than the stress value
                       ! at the bottom of Layer  nl, then the stress value at the top of Layer  nl
                       ! is set equal to zero and the  dstress value for layer  nl  is set equal
                       ! to the stress value at the bottom of Layer  nl

                       !  IS THIS CORRECT??
                       if(stress(nl).gt.stress(nl-1)) then
                          dstress(nl) = stress(nl-1)
                          stress(nl) = 0.0
                       else

                          !  If the stress value at the top of Layer  nl  is greater than
                          ! zero but less than the stress value at the bottom of layer  nl,
                          ! the stress value at the top of Layer  nl  is as calculated before 
                          ! entering this  if  statement and the dstress value is set equal to
                          ! stress(nl) - stress(nl+1)

                          dstress(nl) = stress(nl-1) - stress(nl)
                       endif
                    endif
                    !     Endif NL.LE.NLAY-3

              endif
              !    Calculate the layer midoint wind accaleration value based upon
              !   the layer's value of  dstress  and the mass within the layer

              dudt(nl) = -dstress(nl) / ((pz(nl)-pz(nl-1))/gravv)

              !    Now, end the primary DO 22 loop

              !      write out to look at distribution of breaking vs layer
              if(dudt(nl).ne.0.0) then
                 !     write(17,*) msol,mbin,nl,dudt(nl),1.25+(NLAY-nl)*2.5
!                 write(17,*) msol,mbin,nl,dudt(nl),1.25+(nl-1)*2.5
              endif

              !    Populate the GWAccel array
              !   Be sure that the GWAccel array is properly zeroed upon each entry in to
              !  subroutine GW... the  "goto 22" statement would bypass the GWAccel
              !  assignment below, as would a   "goto 6"  statement

              GWAccel(ilon,ilat,nl,iEast_) = dudt(nl) * Dt

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
!           do nl=1,nlay
!              print *, u(nl),pz(nl),t(nl),theta(nl),dz(nl), &
!                   uav(nl),thtavg(nl),N(nl),dh(nl),Ri(nl),Rimin(nl),stress(nl), &
!                   dstress(nl),ts,dhsat(nl)
!              write(13,*) nl,stress(nl),dstress(nl),dudt,dz(nl),pz(nl),rhogw(nl)
!              write(14,99) nl,stress(nl),dstress(nl),dudt(nl),dz(nl),pz(nl), &
!                   Ri(nl),Rimin(nl),u(nl),t(nl),dh(nl),dhsat(nl),uav(nl), &
!                   eps(nl),1.250+(nl-1)*2.5
!              write(24,98) nl,stress(nl),dstress(nl),dudt(nl),dz(nl),pz(nl), &
!                   Ri(nl),Rimin(nl),u(nl),t(nl),dh(nl),dhsat(nl),uav(nl), &
!                   eps(nl),1.250+(nl-1)*2.5
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
    !  Output: xLowAtmosRadRate(1:nLons,1:nLats,1:nAlts,iBlock)
    !  =======================================================================
    !
    use ModInputs
    use ModSources, only: LowAtmosRadRate
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


    REAL :: DustTime(nDustTimes),ConrathTime(nConrathTimes),invDDiff
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
    real :: FMNETV(LL_NLAYRAD), diffvt, tDiff(nDustTimes),ctDiff(nDustTimes),conrnu,tautot
    real :: fluxupv(LL_NLAYRAD), fluxdnv(LL_NLAYRAD), NFLUXTOPV

    real :: fluxvd(LL_LAYERS),HEATING(LL_LAYERS),TOTAL(LL_LAYERS),rtime,conrathrtime
    real :: fluxid(LL_LAYERS),GREENHOUSE(LL_LAYERS),XLTECORRECTION(LL_LAYERS)

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

     if (UseDustDistribution) then
       DustTime(:) = 0
       ConrathTime(:) = 0

       do N = 1, nDustTimes
          DustTime(N) = TimeDust(N)

       enddo
       do N = 1, nConrathTimes
          ConrathTime(N) = TimeConrath(N)
       enddo

       tDiff = CurrentTime - DustTime
       ctDiff = CurrentTime - ConrathTime
       
       where(tDiff .lt. 1) tDiff = 1e20
       where(ctDiff .lt. 1) ctDiff = 1e20

       iMin = minloc(tDiff)
       cMin = minloc(ctDiff)

       invDDiff = 1 / &
            (DustTime(iMin(1)+1)-DustTime(imin(1)))
       rtime = (CurrentTime - Dusttime(imin(1)))*(HorizontalDustProfile(imin(1) + 1,iLat,iBlock) - &
            HorizontalDustProfile(imin(1),iLat,iBlock)) * invDDiff

       invDDiff = 1 / &
            (ConrathTime(cMin(1)+1)-ConrathTime(cmin(1)))
       conrathrtime = (CurrentTime - Conrathtime(cmin(1)))*&
            (HorizontalConrathProfile(cmin(1) + 1,iLat,iBlock) - &
            HorizontalConrathProfile(cmin(1),iLat,iBlock)) * invDDiff


       DustDistribution(iLon,iLat,iBlock) = HorizontalDustProfile(imin(1),iLat,iBlock) + rtime
       ConrathDistribution(iLon,iLat,iBlock) = HorizontalConrathProfile(imin(1),iLat,iBlock) + conrathrtime

       tautot = DustDistribution(ilon,ilat,iblock)
       conrnu = ConrathDistribution(ilon,ilat,iblock)

    else
       tautot = tautot_temp
       conrnu = conrnu_temp
    endif

    
    CALL DUSTPROFILE(PLEV(L_LEVELS),PTROP,PLEV,TAUCUM,TAUREF,L_LEVELS,TauTot,ConrNU)


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
    end if

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

  subroutine dustprofile(PSF,PTROP,PLEV,TAUCUM,TAUREF,L_LEVELS,TauTot,ConrNU)

    !C Bob's updates 9/17/99 
    !C Reference the dust optical depth to PSF (The surface pressure, mbar),
    !C and modify the way the dust mixing ratio is calculated to more accurately 
    !C reflect the pressure-optical depth relationship.
    !C GCM2.0  Sept 2002
    !C Driver:  Jan 2003 - Modified from GCM 3-D to DRIVER 1-D

    use ModPlanet
    use ModInputs

    implicit none

    integer JSRCHGT
    external JSRCHGT

    INTEGER, PARAMETER :: NPDST = 100
    real, INTENT(in) :: tautot,conrnu
    integer :: n, k, nstar,L_LEVELS
    REAL ::  QRDST(NPDST), PRDST(NPDST), TAUDST(NPDST)
    real ::  PSF, PTROP
    real ::  TAUCUM(LL_LEVELS),TAUREF(LL_LEVELS+1),PLEV(LL_LEVELS+1)
    real ::  refpr, pave, sum, qrdst0, pstar, pstar1, pdif1, pdif2

    !C======================================================================C

    !C Calculate the Reference Pressure Grid (prdst)

    refpr    = (5.0*psf/ptrop)**(1.0/(float(npdst)-1.0))
    prdst(1) = ptrop

    do n=2,npdst
       prdst(n) = refpr*prdst(n-1)
    end do

    !C Calculate the Mixing Ratio at the Reference Pressure Level

    sum = 0.
    do n = 2,npdst
       if (prdst(n).lt.rptau) then
          pave = 0.5*(prdst(n)+prdst(n-1))
          sum  = sum + exp(conrnu*(1.-(rptau/pave)))*&
               (prdst(n)-prdst(n-1))
       end if
       if (prdst(n).ge.rptau) go to 10 
    end do

10  continue

    pave = 0.5*(rptau+prdst(n-1))
    sum  = sum + exp(conrnu*(1.-(rptau/pave)))*(rptau-prdst(n-1))

    !C  GCM1.7  6/28/01   spatially varying dust

    qrdst0 = tautot/sum

    !C Now calculate the mixing ratio at all other levels

    do n=1,npdst-1

       !C Region 1: Mixing ratio changes continuously through the layer

       if (rptau.gt.prdst(n+1)) then
          pave     = 0.5*(prdst(n+1)+prdst(n))
          qrdst(n) = qrdst0*exp(conrnu*(1.0-(rptau/pave)))
       end if

       !C Region 2: Reference pressure level within this layer. 

       if (rptau.le.prdst(n+1).and.rptau.ge.prdst(n)) then
          pave     = 0.5*(prdst(n)+rptau)
          pdif1    = rptau-prdst(n)
          pdif2    = prdst(n+1)-rptau
          qrdst(n) = qrdst0*(&
               exp(conrnu*(1.0-(rptau/pave)))*pdif1+pdif2) / &
               (prdst(n+1)-prdst(n))
       end if

       !C Region 3: Mixing ratio constant

       if (rptau.lt.prdst(n)) then
          qrdst(n) = qrdst0
       end if

    end do

    !C Now compute the optical depths (taudst).

    taudst(1) = 0.0

    do n=2,npdst
       taudst(n) = taudst(n-1) + qrdst(n-1)*(prdst(n)-prdst(n-1))
    end do

    !C  Dust optical depth at the bottom of each sub-layer.

    TAUCUM = 0.0
    TAUREF = 0.0

    DO K=4,L_LEVELS
       PSTAR     = PLEV(K)
       PSTAR1    = MAX(PSTAR,PRDST(1))
       NSTAR     = MIN0(JSRCHGT(NPDST-1,PRDST,1,PSTAR1)-1,NPDST-1)
       TAUCUM(K) = TAUDST(NSTAR)+(PSTAR1-PRDST(NSTAR))*&
            (TAUDST(NSTAR+1) - TAUDST(NSTAR))/&
            (PRDST(NSTAR+1)-PRDST(NSTAR))
       TAUREF(K) = TAUCUM(K) - TAUCUM(K-1)

    END DO


    TAUREF(L_LEVELS+1) = 0.0D0


  END SUBROUTINE DUSTPROFILE

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

  END FUNCTION JSRCHGT
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
       integer,PARAMETER :: NL=101  ! MUST BE LARGER THAN NLL

       REAL :: EM, EP
       REAL :: W0(LL_NLAYRAD), COSBAR(LL_NLAYRAD), DTAU(LL_NLAYRAD)
       REAL :: TAU(LL_NLEVRAD), WDEL(LL_NLAYRAD), CDEL(LL_NLAYRAD)
       REAL :: DTDEL(LL_NLAYRAD), TDEL(LL_NLEVRAD)
       REAL :: FMIDP(LL_NLAYRAD), FMIDM(LL_NLAYRAD)
       REAL :: LAMDA(NL), ALPHA(NL),XK1(NL),XK2(NL)
       REAL :: G1(NL), G2(NL), G3(NL), GAMA(NL), CP(NL), CM(NL), CPM1(NL)
       REAL :: CMM1(NL), E1(NL), E2(NL), E3(NL), E4(NL), EXPTRM(NL)
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

     SUBROUTINE DSOLVER(NL,GAMA,CP,CM,CPM1,CMM1,E1,E2,E3,E4,BTOP,&
          BSURF,RSF,XK1,XK2)

       !C  GCM2.0  Feb 2003
       !C
       !C DOUBLE PRECISION VERSION OF SOLVER

       implicit none

       integer,PARAMETER :: NMAX=201
       integer :: NL, L, N, LM1, LM2, I
       real :: GAMA(NL),CP(NL),CM(NL),CPM1(NL),CMM1(NL),XK1(NL)
       real :: XK2(NL),E1(NL),E2(NL),E3(NL),E4(NL),RSF,BSURF,BTOP
       real :: AF(NMAX),BF(NMAX),CF(NMAX),DF(NMAX),XK(NMAX)
!!$C*********************************************************
!!$C* THIS SUBROUTINE SOLVES FOR THE COEFFICIENTS OF THE    *
!!$C* TWO STREAM SOLUTION FOR GENERAL BOUNDARY CONDITIONS   *
!!$C* NO ASSUMPTION OF THE DEPENDENCE ON OPTICAL DEPTH OF   *
!!$C* C-PLUS OR C-MINUS HAS BEEN MADE.                      *
!!$C* NL     = NUMBER OF LAYERS IN THE MODEL                *
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

       L=2*NL

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

       AF(L) = E1(NL)-RSF*E3(NL)
       BF(L) = E2(NL)-RSF*E4(NL)
       CF(L) = 0.0
       DF(L) = BSURF-CP(NL)+RSF*CM(NL)
!if (BF(L) .ne. bf(L)) then
!   write(*,*) "BF: ",bf(L),e2(Nl),rsf,e4(NL)
!endif
       CALL DTRIDGL(L,AF,BF,CF,DF,XK)

       !C     ***UNMIX THE COEFFICIENTS****

       DO 28 N=1,NL
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
!write(*,*) taucumi(:,3,15)
!write(*,*) "dtauk: ",dtauki(:,3,15)
!write(*,*) "dtaui: ",dtaui(:,3,15)
!stop

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
!write(*,*) "before: ",taucumi(1,nw,ng),nw,ng
!write(*,*) "before: ",taucumi(:,nw,ng)
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
!if (l .eq. 45 .and. ng .eq.15 .and. nw .eq. 3) then 
!stop
!write(*,*) fmneti(l),fmupi(l),fmdi(l),dwni(nw),gweight(ng),fzeroi(nw)
!endif
                   FMNETI(L)  = FMNETI(L)+(FMUPI(L)-FMDI(L))*DWNI(NW)*&
                        GWEIGHT(NG)*(1.0-FZEROI(NW))
                   FLUXUPI(L) = FLUXUPI(L) + FMUPI(L)*DWNI(NW)*GWEIGHT(NG)*&
                        (1.0-FZEROI(NW))
                   FLUXDNI(L) = FLUXDNI(L) + FMDI(L)*DWNI(NW)*GWEIGHT(NG)*&
                        (1.0-FZEROI(NW))
!if (l .eq. 45 .and. ng .eq.15 .and. nw .eq. 3) then 
!write(*,*) fmneti(l)
!stop
!endif
!if (fmneti(l) .ne. fmneti(l)) then
!   write(*,*) l,l_nlevrad-1,ng,l_ngauss-1,nw,l_nspecti
!stop
!endif
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
             integer,PARAMETER :: NL=101 ! MUST BE LARGER THAN NLEVEL 


             INTEGER :: NLL, NLAYER, L, NW, NT, NT2
             REAL ::  TERM, CPMID, CMMID
             REAL ::  EM,EP
             REAL ::  COSBAR(LL_NLAYRAD), W0(LL_NLAYRAD), DTAU(LL_NLAYRAD)
             REAL ::  TAUCUM(LL_LEVELS), DTAUK
             REAL ::  TLEV(LL_LEVELS)
             REAL ::  WAVEN, DW, RSF
             REAL ::  BTOP, BSURF, FMIDP(LL_NLAYRAD), FMIDM(LL_NLAYRAD)
             REAL ::  B0(NL),B1(NL),ALPHA(NL),LAMDA(NL),XK1(NL),XK2(NL)
             REAL ::  GAMA(NL),CP(NL),CM(NL),CPM1(NL),CMM1(NL),E1(NL),E2(NL)
             REAL ::  E3(NL),E4(NL)
             REAL ::  TAUMAX

             REAL ::  FTOPUP, FLUXUP, FLUXDN


             DATA TAUMAX / L_TAUMAX /

             !C======================================================================C

             !C     WE GO WITH THE HEMISPHERIC CONSTANT APPROACH IN THE INFRARED
!write(*,*)"gflux :", taucum
!stop
             IF (NLL .GT. NL) STOP 'PARAMETER NL TOO SMALL IN GLUFV'

             NLAYER = L_NLAYRAD

             DO L=1,L_NLAYRAD-1
                ALPHA(L) = SQRT( (1.0-W0(L))/(1.0-W0(L)*COSBAR(L)) )
                LAMDA(L) = ALPHA(L)*(1.0-W0(L)*COSBAR(L))/UBARI

                NT2   = TLEV(2*L+2)*10.0D0-499
                NT    = TLEV(2*L)*10.0D0-499

                B1(L) = (PLANCKIR(NW,NT2)-PLANCKIR(NW,NT))/DTAU(L)
                B0(L) = PLANCKIR(NW,NT)
             END DO

             !C     Take care of special lower layer

             L        = L_NLAYRAD
             ALPHA(L) = SQRT( (1.0-W0(L))/(1.0-W0(L)*COSBAR(L)) )
             LAMDA(L) = ALPHA(L)*(1.0-W0(L)*COSBAR(L))/UBARI

             NT    = TLEV(2*L+1)*10.0D0-499
             NT2   = TLEV(2*L)*10.0D0-499
             B1(L) = (PLANCKIR(NW,NT)-PLANCKIR(NW,NT2))/DTAU(L)
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

!if (l .eq. 45) then
!   write(*,*) "in glux: ",fmidp(l),fmidm(l), taucum(2*l+1),taucum(2*l)

!endif
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
