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
  TrAltered = Tr
  where(TrAltered < 235.0) TrAltered=235.0
  IonCollisions(:,:,:,iO_4SP_,iO_3P_) = &
       3.67e-17 * NDensityS(:,:,:,iO_3P_,iBlock) * &
       sqrt(TrAltered) * (1.0-0.064*log10(TrAltered))**2

  IonCollisions(:,:,:,iO_4SP_,iO2_) =&
       6.64e-16*NDensityS(:,:,:,iO2_,iBlock)
  IonCollisions(:,:,:,iO_4SP_,iN2_) =&
       6.82e-16*NDensityS(:,:,:,iN2_,iBlock)
  IonCollisions(:,:,:,iO_4SP_,iN_4S_) =&
       4.62e-16*NDensityS(:,:,:,iN_4S_,iBlock)
  ! This is an average of O2 and N2, since NO doesn't exist
  IonCollisions(:,:,:,iO_4SP_,iNO_) =&
       6.73e-16*NDensityS(:,:,:,iNO_,iBlock)

  Collisions(:,:,:,iVIN_) = IonCollisions(:,:,:,iO_4SP_,1)
  do iSpecies = 2, nSpecies
     Collisions(:,:,:,iVIN_) = &
          Collisions(:,:,:,iVIN_) + IonCollisions(:,:,:,iO_4SP_,iSpecies)
  enddo

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

