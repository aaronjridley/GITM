!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine calc_rates(iBlock)

  use ModGITM
  use ModRates
  use ModConstants
  use ModPlanet
  use ModInputs
  use ModEUV, only : SunOrbitEccentricity
  use ModSources, only : KappaEddyDiffusion

  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iIon, iSpecies, iError, iiAlt, iLat,iLon

  real, dimension(nLons, nLats, nAlts) :: &
       Tn, Ti, TWork1, TWork2, TWork3, NO2

  real :: ScaleHeight(nLons, nLats)

  !! Titan-specific variables used to calculate
  !! ViscCoef and KappaTemp
  real, dimension(nLons,nLats,0:nAlts+1) :: &
       KappaN2,          & 
       KappaCH4,         &
       KappaH2,         &
       cpn2, cpch4,      &
       ViscN2, ViscCH4,  &
       ViscH2,           &
       TL,               &
       eta0, temp0,      &
       cn2,              &
       Gnc,Gcn,          &
       xnc,xcn,          &
       knc,kcn,          &
       EddyDiffProfile,  &
       ones

  real :: mnc,mcn 
  real :: fnc 
  real :: TimeFactor

  call report("calc_rates",2)
  call start_timing("calc_rates")

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> mean major mass", iblock

  ! We add 1 because this is in the denominator a lot, and the corners 
  ! don't have anything.

  MeanIonMass = 0.0
  MeanMajorMass = 0.0
  do iSpecies = 1, nSpecies
     MeanMajorMass = MeanMajorMass + &
          Mass(iSpecies) * &
          NDensityS(:,:,:,iSpecies,iBlock)/(NDensity(:,:,:,iBlock)+1.0)
  enddo

  ! Once again, in the corners, the meanmajormass is 0.
  where (MeanMajorMass == 0) MeanMajorMass = Mass(1)

  do iIon = 1, nIons-1
     MeanIonMass = MeanIonMass + &
          MassI(iIon) * IDensityS(:,:,:,iIon,iBlock) / &
          IDensityS(:,:,:,ie_,iBlock)
  enddo

  !TempUnit = MeanMajorMass / Boltzmanns_Constant       
  TempUnit = Mass(iN2_) / Boltzmanns_Constant       

  TL(1:nLons,1:nLats,0:nAlts+1) = &
     Temperature(1:nLons,1:nLats,0:nAlts+1,iBlock)*&
       TempUnit(1:nLons,1:nLats,0:nAlts+1)

  ones = 1.0

  fnc = 1.065/(2.0*sqrt(2.0))
  mnc = Mass(iN2_)/Mass(iCH4_)
  mcn = 1/mnc

  xnc = NDensityS(1:nLons,1:nLats,0:nAlts+1,iN2_,iBlock)/&
        NDensityS(1:nLons,1:nLats,0:nAlts+1,iCH4_,iBlock)

  xcn = NDensityS(1:nLons,1:nLats,0:nAlts+1,iCH4_,iBlock)/&
        NDensityS(1:nLons,1:nLats,0:nAlts+1,iN2_,iBlock)

!\
! ViscN2, ViscCH4 are in units of kg/m/s
! ViscCH4 (Yaws 1995)
! ViscN2 Sutherland Formula(LMNO 2003)
!/
  eta0 = 0.01781e-03 ! in kg/m/s 
  cn2 = 111.0        ! in K
  temp0 = 300.55     ! in K

  ViscN2(:,:,:) =           &
   (   eta0*(cn2 + temp0)/(TL + cn2)  )*( (TL/temp0)**(1.50) )   
      
  ViscCH4(:,:,:)  =  &
   (  3.8435*ones + (4.0112e-01)*(TL) +  &
   (1.4303e-04)*((TL)**2.0) )*(1.0e-07)

  ViscH2(:,:,:) = (1.4648e-07)*(TL**0.716)

  TimeFactor = 0.0
  ViscCoef(1:nLons,1:nLats,0:nAlts+1) = (1.0 + TimeFactor)*&
  (  Mass(iN2_)*NDensityS(1:nLons,1:nLats,0:nAlts+1,iN2_,iBlock)*ViscN2 +  &
     Mass(iCH4_)*NDensityS(1:nLons,1:nLats,0:nAlts+1,iCH4_,iBlock)*ViscCH4 + &
     Mass(iH2_)*NDensityS(1:nLons,1:nLats,0:nAlts+1,iH2_,iBlock)*ViscH2  )/ &
     Rho(1:nLons,1:nLats,0:nAlts+1,iBlock)
!\
! KappaN2, CH4 are in units of J/m/s/K (Yaws 1997)
!/
  KappaN2(:,:,:)  =  &
     0.00309*ones + (7.593e-05)*(TL) -  &
     (1.1014e-08)*((TL)**2.0)

  KappaCH4(:,:,:)  =  &
     -0.00935*ones + (1.4028e-04)*(TL) +  &
     (3.318e-08)*((TL)**2.0)

  KappaH2(:,:,:)  =  (1.262e-03)*(TL**0.876)

  knc(:,:,:) = &
     (  (   ViscN2(:,:,:)/      &
            ViscCH4(:,:,:)  )/   &       
      mcn)**(0.5)

  kcn(:,:,:) = &
     (   (   ViscCH4(:,:,:)/      &
              ViscN2(:,:,:)  )/   &       
      mnc)**(0.5)

  Gnc(:,:,:) =         &
      ( fnc/sqrt(1+mnc) ) *              &
      (   ones + knc*((mnc)**(0.25))  )**(2.0)

  Gcn(:,:,:) =         &
      ( fnc/sqrt(1+mcn) ) *              &
      (   ones + kcn*((mcn)**(0.25))  )**(2.0)
      

!  KappaTemp(:,:,0:nAlts+1,iBlock) = &
!          (KappaN2/(ones + Gnc*xcn) +   &
!          KappaCH4/(ones + Gcn*xnc)  ) 

! Use approximate form here
  KappaTemp(:,:,0:nAlts+1,iBlock) = &
   (  KappaN2(:,:,0:nAlts+1)*NDensityS(1:nLons,1:nLats,0:nAlts+1,iN2_ ,iBlock) + &
     KappaCH4(:,:,0:nAlts+1)*NDensityS(1:nLons,1:nLats,0:nAlts+1,iCH4_,iBlock) + &
      KappaH2(:,:,0:nAlts+1)*NDensityS(1:nLons,1:nLats,0:nAlts+1,iH2_ ,iBlock))/ &
      NDensity(1:nLons,1:nLats,0:nAlts+1,iBlock)

  ThermalDiffCoefS(1:nSpecies) = 0.0
  ThermalDiffCoefS(iH_)  = -0.38
  ThermalDiffCoefS(iH2_) = -0.38
  ThermalDiffCoefS(iAr_) =  0.17


  call end_timing("calc_rates")

end subroutine calc_rates

subroutine calc_collisions(iBlock)

  use ModGITM
  use ModRates
  use ModConstants
  use ModPlanet
  use ModInputs

  implicit none

  integer, intent(in) :: iBlock

  integer :: iError, iSpecies

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: &
       Ne, mnd, Te, tmp, Tn, Ti, Tr, TrAltered

  call start_timing("calc_rates")

  !\
  ! Need to get the neutral, ion, and electron temperature
  !/

  Tn = Temperature(:,:,:,iBlock)*&
       TempUnit(:,:,:)
  Ti = ITemperature(:,:,:,iBlock)

  Tr = (Tn+Ti)/2

  mnd = NDensity(:,:,:,iBlock)+1.0
  Ne  = IDensityS(:,:,:,ie_,iBlock)

  !\
  ! -----------------------------------------------------------
  ! Collision Frequencies
  ! -----------------------------------------------------------
  !/

  e_gyro = &
       Element_Charge * B0(:,:,:,iMag_,iBlock) / Mass_Electron

!
! Ion Neutral Collision Frequency (From Kelley, 1989, pp 460):
!

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> vin",iblock

  Collisions(:,:,:,iVIN_) = 2.6e-15 * (mnd + Ne)/sqrt(MeanMajorMass/AMU)

!
! From Schunk and Nagy table 4.4 & 4.5
!  TrAltered = Tr
!  where(TrAltered < 235.0) TrAltered=235.0
!  IonCollisions(:,:,:,iO_4SP_,iO_3P_) = &
!       3.67e-17 * NDensityS(:,:,:,iO_3P_,iBlock) * &
!       sqrt(TrAltered) * (1.0-0.064*log10(TrAltered))**2


!  IonCollisions(:,:,:,iO_4SP_,iO2_) =&
!       6.64e-16*NDensityS(:,:,:,iO2_,iBlock)
!  IonCollisions(:,:,:,iO_4SP_,iN2_) =&
!       6.82e-16*NDensityS(:,:,:,iN2_,iBlock)
!  IonCollisions(:,:,:,iO_4SP_,iN_4S_) =&
!       4.62e-16*NDensityS(:,:,:,iN_4S_,iBlock)
  ! This is an average of O2 and N2, since NO doesn't exist
!  IonCollisions(:,:,:,iO_4SP_,iNO_) =&
!       6.73e-16*NDensityS(:,:,:,iNO_,iBlock)

!  Collisions(:,:,:,iVIN_) = IonCollisions(:,:,:,iO_4SP_,1)
!  do iSpecies = 2, nSpecies
!     Collisions(:,:,:,iVIN_) = &
!          Collisions(:,:,:,iVIN_) + IonCollisions(:,:,:,iO_4SP_,iSpecies)
!  enddo

!
! Electron Neutral Collision Frequency
!

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> ven", iblock

  Te = eTemperature(:,:,:,iBlock)
  where(te == 0.0) te = 1000.0
  Collisions(:,:,:,iVEN_) = 5.4e-16 * (mnd)*sqrt(Te)

!
! Electron Ion Collision Frequency
!

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> vei", iblock

  tmp = (34.0 + 4.18*log((TE**3.0)/(Ne*1.0e-6)))
  Collisions(:,:,:,iVEI_) = tmp*Ne*TE**(-3.0/2.0) * 1.0e-6

  i_gyro = Element_Charge * B0(:,:,:,iMag_,iBlock) / MeanIonMass

  call end_timing("calc_rates")

end subroutine calc_collisions

subroutine calc_viscosity(iBlock)

  use ModGITM

  implicit none

  integer, intent(in) :: iBlock

  integer :: iSpecies
  real :: TimeFactor
  real, dimension(nLons,nLats,-1:nAlts+2) :: &
       ViscN2, ViscCH4, ViscH2,  &
       TL, eta0, temp0,          &
       cn2,              &
       Gnc,Gcn,          &
       xnc,xcn,          &
       knc,kcn,          &
       EddyDiffProfile,  &
       Nu_H2_H2,  &
       ones

  real :: mnc,mcn 
  real :: fnc 

  TL(1:nLons,1:nLats,-1:nAlts+2) = &
     Temperature(1:nLons,1:nLats,-1:nAlts+2,iBlock)*&
        TempUnit(1:nLons,1:nLats,-1:nAlts+2)
          
  ones = 1.0

  fnc = 1.065/(2.0*sqrt(2.0))
  mnc = Mass(iN2_)/Mass(iCH4_)
  mcn = 1/mnc

  xnc = NDensityS(1:nLons,1:nLats,-1:nAlts+2,iN2_,iBlock)/&
        NDensityS(1:nLons,1:nLats,-1:nAlts+2,iCH4_,iBlock)

  xcn = NDensityS(1:nLons,1:nLats,-1:nAlts+2,iCH4_,iBlock)/&
        NDensityS(1:nLons,1:nLats,-1:nAlts+2,iN2_,iBlock)

!\
! ViscN2, ViscCH4 are in units of kg/m/s
! ViscCH4 (Yaws 1995)
! ViscN2 Sutherland Formula(LMNO 2003)
!/
  eta0 = 0.01781e-03 ! in kg/m/s 

  cn2 = 111.0        ! in K
  temp0 = 300.55     ! in K

  ViscN2(:,:,:) =           &
      (   eta0*(cn2 + temp0)/(TL + cn2)  )*( (TL/temp0)**(1.50) )   
      
  ViscCH4(:,:,:)  =  &
      (  3.8435*ones + (4.0112e-01)*(TL) +  &
      (1.4303e-04)*((TL)**2.0) )*(1.0e-07)

!!! Shuts down the enhanced viscosity

  ViscCoef(1:nLons,1:nLats,0:nAlts+1) = &
      (Mass(iN2_)*NDensityS(1:nLons,1:nLats,0:nAlts+1,iN2_,iBlock)*ViscN2(:,:,0:nAlts+1) +  &
      Mass(iCH4_)*NDensityS(1:nLons,1:nLats,0:nAlts+1,iCH4_,iBlock)*ViscCH4(:,:,0:nAlts+1) )/&
         Rho(1:nLons,1:nLats,0:nAlts+1,iBlock)

  do iSpecies = 1, nSpecies
    ViscCoefS(1:nLons,1:nLats, 0:nAlts+1,iSpecies) = &
         ViscCoef(1:nLons,1:nLats, 0:nAlts+1)
  enddo 

  ViscCoefS(1:nLons,1:nLats, 0:nAlts+1,iN2_) = &
    ViscN2(1:nLons,1:nLats, 0:nAlts+1)

 ViscCoefS(1:nLons,1:nLats, 0:nAlts+1,iCH4_) = &
    ViscCH4(1:nLons,1:nLats, 0:nAlts+1)

     ViscCoefS(1:nLons,1:nLats,0:nAlts+1,i15N2_) = &
        ViscN2(1:nLons,1:nLats,0:nAlts+1)

     ViscCoefS(1:nLons,1:nLats, 0:nAlts+1,i13CH4_) = &
        ViscCH4(1:nLons,1:nLats, 0:nAlts+1)

!!!! Cui et al. [2008] suggested 5.5 x 10^-6 kg/m/s
!     ViscCoefS(1:nLons,1:nLats, 0:nAlts+1,iH2_) = &
!           (1.4648e-07)*(TL(1:nLons,1:nLats,0:nAlts+1)**0.716)


! Reduce Viscosity consistent with Boqueho and Blelly for minor species
!     ViscCoefS(1:nLons,1:nLats, 0:nAlts+1,iH_)  = (1.0000e-07)
     ViscCoefS(1:nLons,1:nLats, 0:nAlts+1,iH2_) = &
           (1.4648e-07)*(TL(1:nLons,1:nLats,0:nAlts+1)**0.716)

     ! Remove the atomic viscosities

!     ViscCoefS(1:nLons,1:nLats,-1:nAlts+2,iH_)  = 0.0
!     ViscCoefS(1:nLons,1:nLats,-1:nAlts+2,iN4S_)  = 0.0

     ViscCoefS(1:nLons,1:nLats,-1:nAlts+2,iH_)  = 1.0e-06

!     ViscCoefS(1:nLons,1:nLats,-1:nAlts+2,iH_)  = &
!           (2.0715e-07)*(TL**0.716)

    ! ViscCoefS(1:nLons,1:nLats,-1:nAlts+2,iH_ ) = 0.0
    ! ViscCoefS(1:nLons,1:nLats,-1:nAlts+2,iH2_) = 0.0

end subroutine calc_viscosity

