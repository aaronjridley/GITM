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

  if (iDebugLevel > 4) write(*,*) "=====> kappatemp", iblock

  do iAlt = 0, nAlts+1
     
     ! JMB: 06/24/2016:  Added He Conduction contribution
     ! kT = 3.736e-03 and s = 0.648

     KappaTemp(:,:,iAlt,iBlock) = &
          (NDensityS(1:nLons,1:nLats,iAlt,iO2_,iBlock) / &
          NDensity(1:nLons,1:nLats,iAlt,iBlock) + &
          NDensityS(1:nLons,1:nLats,iAlt,iN2_,iBlock)/ &
          NDensity(1:nLons,1:nLats,iAlt,iBlock)) * ThermalConduction_AO2 * &
          (Temperature(1:nLons,1:nLats,ialt,iBlock)* &
          TempUnit(1:nLons,1:nLats,iAlt))**ThermalConduction_s + &
          (NDensityS(1:nLons,1:nLats,iAlt,iO_3P_,iBlock)/&
          NDensity(1:nLons,1:nLats,iAlt,iBlock)*ThermalConduction_AO) * &
          (Temperature(1:nLons,1:nLats,iAlt,iBlock) * &
          TempUnit(1:nLons,1:nLats,iAlt))**ThermalConduction_s + &
          (NDensityS(1:nLons,1:nLats,iAlt,iHe_,iBlock)/&
            NDensity(1:nLons,1:nLats,iAlt,iBlock)*(3.736e-03)) * &
          (Temperature(1:nLons,1:nLats,iAlt,iBlock) * &
          TempUnit(1:nLons,1:nLats,iAlt))**0.648
     
  enddo

  ! Thermal Diffusion is zero for all but the lightest species
  ! Banks and Kockarts suggest Alpha_T = -0.38 for He
  ThermalDiffCoefS(1:nSpecies) = 0.0 
  ThermalDiffCoefS(iHe_) = -0.38

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

  integer :: iError, iSpecies, iLon, iLat, iAlt

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: &
       Ne, mnd, Te, tmp, Tn, Ti, Tr, &
       NeMajor, Frac, ODensity, NitroDensity

  call start_timing("calc_rates")

  !\
  ! Need to get the neutral, ion, and electron temperature
  !/

  Tn = Temperature(:,:,:,iBlock)*TempUnit(:,:,:)
  Ti = ITemperature(:,:,:,iBlock)

  Tr = (Tn+Ti)/2
  ! Set a floor on Tr:
  where (Tr < 200) Tr = 200.0

  mnd = NDensity(:,:,:,iBlock)+1.0
  Ne  = IDensityS(:,:,:,ie_,iBlock)

  ODensity = &
       NDensityS(:,:,:,iO_3P_,iBlock) + &
       NDensityS(:,:,:,iO_1D_,iBlock)
       
  NitroDensity = &
       NDensityS(:,:,:,iN_4S_,iBlock) + &
       NDensityS(:,:,:,iN_2P_,iBlock) + &
       NDensityS(:,:,:,iN_2D_,iBlock)
       
  ! We include O+, O2+, NO+ below
  NeMajor  = 0.0
  do iSpecies = 1, nIonsAdvect
     NeMajor = NeMajor + IDensityS(:,:,:,iSpecies,iBlock)
  enddo

  !\
  ! -----------------------------------------------------------
  ! Collision Frequencies
  ! -----------------------------------------------------------
  !/

  e_gyro = &
       Element_Charge * B0(:,:,:,iMag_,iBlock) / Mass_Electron

  !
  ! Ion Neutral Collision Frequency (From Kelley, 1989, pp 460):
  ! This is the old-school method, and should not be used.
  !

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> vin",iblock

  Collisions(:,:,:,iVIN_) = 2.6e-15 * (mnd + Ne)/sqrt(MeanMajorMass/AMU)

  IonCollisions = 0.0
  
  !
  ! From Schunk and Nagy table 4.4 & 4.5
  !

  ! O+ (with O, O2, N2, N, NO)

  if (iDebugLevel > 5) write(*,*) "======> o+ ",iblock

  where (Tr > 235.0) IonCollisions(:,:,:,iO_4SP_,iO_3P_) = &
       3.67e-17 * ODensity * sqrt(Tr) * (1.0-0.064*log10(Tr))**2
  where(tr <= 235) IonCollisions(:,:,:,iO_4SP_,iO_3P_) = 8.6e-16 * ODensity
  
  IonCollisions(:,:,:,iO_4SP_,iO2_)   = 6.64e-16*NDensityS(:,:,:,iO2_,iBlock)
  IonCollisions(:,:,:,iO_4SP_,iN2_)   = 6.82e-16*NDensityS(:,:,:,iN2_,iBlock)
  IonCollisions(:,:,:,iO_4SP_,iN_4S_) = 4.62e-16*NDensityS(:,:,:,iN_4S_,iBlock)
  ! This is an average of O2 and N2, since NO doesn't exist
  IonCollisions(:,:,:,iO_4SP_,iNO_)   = 6.73e-16*NDensityS(:,:,:,iNO_,iBlock)

  IonCollisions(:,:,:,iO_2DP_,iO_3P_) = IonCollisions(:,:,:,iO_4SP_,iO_3P_)
  IonCollisions(:,:,:,iO_2DP_,iO2_)   = IonCollisions(:,:,:,iO_4SP_,iO2_)
  IonCollisions(:,:,:,iO_2DP_,iN2_)   = IonCollisions(:,:,:,iO_4SP_,iN2_)
  IonCollisions(:,:,:,iO_2DP_,iN_4S_) = IonCollisions(:,:,:,iO_4SP_,iN_4S_)
  IonCollisions(:,:,:,iO_2DP_,iNO_)   = IonCollisions(:,:,:,iO_4SP_,iNO_)
  
  IonCollisions(:,:,:,iO_2PP_,iO_3P_) = IonCollisions(:,:,:,iO_4SP_,iO_3P_)
  IonCollisions(:,:,:,iO_2PP_,iO2_)   = IonCollisions(:,:,:,iO_4SP_,iO2_)
  IonCollisions(:,:,:,iO_2PP_,iN2_)   = IonCollisions(:,:,:,iO_4SP_,iN2_)
  IonCollisions(:,:,:,iO_2PP_,iN_4S_) = IonCollisions(:,:,:,iO_4SP_,iN_4S_)
  IonCollisions(:,:,:,iO_2PP_,iNO_)   = IonCollisions(:,:,:,iO_4SP_,iNO_)
  
  ! O2+ (with O2, O, N2, N, NO)

  if (iDebugLevel > 5) write(*,*) "======> o2+ ",iblock

  where (tr > 800.) IonCollisions(:,:,:,iO2P_,iO2_) = &
       2.59e-17 * NDensityS(:,:,:,iO2_,iBlock) * tr**0.5 * (1 - 0.073*log10(tr))**2
  where (tr <= 800) IonCollisions(:,:,:,iO2P_,iO2_) = 8.2e-16 * NDensityS(:,:,:,iO2_,iBlock)

  IonCollisions(:,:,:,iO2P_,iO_3P_) = 2.31e-16 * ODensity
  IonCollisions(:,:,:,iO2P_,iN2_)   = 4.13e-16 * NDensityS(:,:,:,iN2_,iBlock)
  IonCollisions(:,:,:,iO2P_,iN_4S_) = 2.64e-16 * NDensityS(:,:,:,iN_4S_,iBlock)
  ! This is just N2, since NO doesn't exist (ave CO and N2)
  IonCollisions(:,:,:,iO2P_,iNO_)   = 4.25e-16*NDensityS(:,:,:,iNO_,iBlock)

  ! NO+ (with NO, O, N2, N, O2)

  if (iDebugLevel > 5) write(*,*) "======> no+ ",iblock
  
  ! This resonant ion-neutral is made up, since NO-NO+ doesn't exist
  where (tr > 800.) IonCollisions(:,:,:,iNOP_,iNO_) = &
       2.59e-17 * NDensityS(:,:,:,iNO_,iBlock) * tr**0.5 * (1 - 0.073*log10(tr))**2
  where (tr <= 800) IonCollisions(:,:,:,iNOP_,iNO_) = 8.2e-16 * NDensityS(:,:,:,iNO_,iBlock)

  IonCollisions(:,:,:,iNOP_,iO_3P_) = 2.44e-16 * ODensity
  IonCollisions(:,:,:,iNOP_,iN2_)   = 4.34e-16 * NDensityS(:,:,:,iN2_,iBlock)
  IonCollisions(:,:,:,iNOP_,iN_4S_) = 2.79e-16 * NDensityS(:,:,:,iN_4S_,iBlock)
  IonCollisions(:,:,:,iNOP_,iO2_)   = 4.27e-16 * NDensityS(:,:,:,iO2_,iBlock)

  ! N2+ (with N2, O2, O, N, NO)

  if (iDebugLevel > 5) write(*,*) "======> n2+ ",iblock

  IonCollisions(:,:,:,iN2P_,iN2_) = &
       5.14e-17 * NDensityS(:,:,:,iN2_,iBlock) * tr**0.5 * (1 - 0.069*log10(tr))**2

  IonCollisions(:,:,:,iN2P_,iO2_)   = 4.49e-16 * NDensityS(:,:,:,iO2_,iBlock)
  IonCollisions(:,:,:,iN2P_,iO_3P_) = 2.58e-16 * ODensity
  IonCollisions(:,:,:,iN2P_,iN_4S_) = 2.95e-16 * NDensityS(:,:,:,iN_4S_,iBlock)
  ! This is just N2, since NO doesn't exist (ave CO and O2)
  IonCollisions(:,:,:,iN2P_,iNO_)   = 4.66e-16 * NDensityS(:,:,:,iNO_,iBlock)

  ! N+ (with N, O2, N2, O, NO)

  if (iDebugLevel > 5) write(*,*) "======> n+ ",iblock

  ! where (Tr > 275.0)
  IonCollisions(:,:,:,iNP_,iN_4S_) = &
       3.83e-17 * NitroDensity * sqrt(Tr) * (1.0-0.063*log10(Tr))**2
  ! where (tr <= 275) IonCollisions(:,:,:,iO_4SP_,iO_3P_) = 8.6e-16 * ODensity

  IonCollisions(:,:,:,iNP_,iO2_)   = 7.25e-16*NDensityS(:,:,:,iO2_,iBlock)
  IonCollisions(:,:,:,iNP_,iN2_)   = 7.47e-16*NDensityS(:,:,:,iN2_,iBlock)
  IonCollisions(:,:,:,iNP_,iO_3P_) = 4.42e-16*NDensityS(:,:,:,iO_3P_,iBlock)
  ! This is an average of O2 and N2, since NO doesn't exist
  IonCollisions(:,:,:,iNP_,iNO_)   = 7.36e-16*NDensityS(:,:,:,iNO_,iBlock)

  ! I don't think that this method is right, but here we go:
  
!  Collisions(:,:,:,iVIN_) = 0.0
!
!  ! O+
!  do iSpecies = 1, nSpecies
!     Frac = IDensityS(:,:,:,iO_4SP_,iBlock) / NeMajor
!     Collisions(:,:,:,iVIN_) = &
!          Collisions(:,:,:,iVIN_) + &
!          IonCollisions(:,:,:,iO_4SP_,iSpecies) * Frac
!  enddo
!
!  ! O2+
!  do iSpecies = 1, nSpecies
!     Frac = IDensityS(:,:,:,iO2P_,iBlock) / NeMajor
!     Collisions(:,:,:,iVIN_) = &
!          Collisions(:,:,:,iVIN_) + &
!          IonCollisions(:,:,:,iO2P_,iSpecies) * Frac
!  enddo
!
!  ! NO+
!  do iSpecies = 1, nSpecies
!     Frac = IDensityS(:,:,:,iNOP_,iBlock) / NeMajor
!     Collisions(:,:,:,iVIN_) = &
!          Collisions(:,:,:,iVIN_) + &
!          IonCollisions(:,:,:,iNOP_,iSpecies) * Frac
!  enddo

!
! Electron Neutral Collision Frequency
!

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> ven", iblock

  Te = eTemperature(:,:,:,iBlock)
  where (Te <= 0.0) Te = 1000.0

  Collisions(:,:,:,iVEN_) = 5.4e-16 * mnd * sqrt(Te)
  
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
  use ModInputs, only: UseTestViscosity, TestViscosityFactor
  implicit none

  integer, intent(in) :: iBlock

  integer :: iSpecies
  ! This is Earth-based, and 

  if (UseTestViscosity) then

     ! TempUnit is mmm/boltzman
     ! visc should be sqrt(Tn * mmm / Boltz) * constant 
     ! The constant sets the viscosity to roughly 3.5 times the old
     ! viscosity for March 2-5, 2013. Not sure if this is ideal, so it
     ! should be tested a bit.

     ViscCoef(1:nLons,1:nLats,0:nAlts+1) = TestViscosityFactor * &
          0.00013 * sqrt(Temperature(1:nLons,1:nLats,0:nAlts+1,iBlock) * &
          TempUnit(1:nLons,1:nLats,0:nAlts+1) * &
          TempUnit(1:nLons,1:nLats,0:nAlts+1))

     do iSpecies = 1, nSpecies
        ViscCoefS(1:nLons,1:nLats,0:nAlts+1,iSpecies) = & 
             ViscCoef(1:nLons,1:nLats,0:nAlts+1) * &
             sqrt(Mass(iSpecies) / MeanMajorMass(1:nLons,1:nLats,0:nAlts+1))
     enddo

  else

     ViscCoef(1:nLons,1:nLats,0:nAlts+1) = 4.5e-5 * &
          (Temperature(1:nLons,1:nLats,0:nAlts+1,iBlock)*&
          TempUnit(1:nLons,1:nLats,0:nAlts+1)/ 1000.)**(-0.71)

     do iSpecies = 1, nSpecies
        ViscCoefS(1:nLons,1:nLats,0:nAlts+1,iSpecies) = & 
             ViscCoef(1:nLons,1:nLats,0:nAlts+1)
     enddo

  endif
  

end subroutine calc_viscosity

