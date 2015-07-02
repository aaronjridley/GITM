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

  real :: e2

!! Titan-specific variables used to calculate
!! ViscCoef and KappaTemp

  real, dimension(nLons,nLats,0:nAlts+1) :: &
       KappaN2,          & 
       KappaCH4,         &
       KappaINGO,        & 
       cpn2, cpch4,      &
       MeanMolecMass,    & 
       ViscN2, ViscCH4,  &
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

  TempUnit = MeanMajorMass / Boltzmanns_Constant       

  !\
  ! These are needed for the Euv Heating and other thermodynamics:
  !/

  if (iDebugLevel > 4) write(*,*) "=====> cp and kappatemp", iblock

     TL(1:nLons,1:nLats,0:nAlts+1) = &
        Temperature(1:nLons,1:nLats,0:nAlts+1,iBlock)*  &
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

     ViscCoef(:,:,0:nAlts+1) = 1.0*&
     (  Mass(iN2_)*NDensityS(1:nLons,1:nLats,0:nAlts+1,iN2_,iBlock)*ViscN2 +  &
        Mass(iCH4_)*NDensityS(1:nLons,1:nLats,0:nAlts+1,iCH4_,iBlock)*ViscCH4  )/ &
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
      
     KappaTemp(:,:,0:nAlts+1,iBlock) = &
          (KappaN2/(ones + Gnc*xcn) +   &
          KappaCH4/(ones + Gcn*xnc)  ) 

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

  real, dimension(nLons, nLats, nAlts) :: Tn, Ti, e2

  integer :: iError

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: &
       Ne, mnd, Te, tmp

  call start_timing("calc_rates")

  !\
  ! Need to get the neutral, ion, and electron temperature
  !/

  Tn = Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
       TempUnit(1:nLons,1:nLats,1:nAlts)
  Ti = ITemperature(1:nLons,1:nLats,1:nAlts,iBlock)

  mnd = NDensity(:,:,:,iBlock)+1.0
  Ne  = IDensityS(:,:,:,ie_,iBlock)

  !\
  ! -----------------------------------------------------------
  ! Collision Frequencies
  ! -----------------------------------------------------------
  !/

  e_gyro = &
       Element_Charge * B0(:,:,:,iMag_,iBlock) / Mass_Electron

  e2 = Element_Charge * Element_Charge

!
! Ion Neutral Collision Frequency (From Kelley, 1989, pp 460):
!

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> vin",iblock

  Collisions(:,:,:,iVIN_) = 2.6e-15 * (mnd + Ne)/sqrt(MeanMajorMass/AMU)

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


end subroutine calc_viscosity

