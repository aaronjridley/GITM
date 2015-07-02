!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine fill_photo(photoion, photoabs, photodis)

  use ModPlanet
  use ModEUV

  implicit none

  real, intent(out) :: photoion(Num_WaveLengths_High, nIons-1)
  real, intent(out) :: photoabs(Num_WaveLengths_High, nSpecies)
  real, intent(out) :: photodis(Num_WaveLengths_High, nSpecies)

  integer :: iSpecies, iWave

  PhotoAbs = 0.0
  PhotoIon = 0.0
  PhotoDis = 0.0

  photoabs(:,iH2_)    = PhotoAbs_H2
  photoion(:,iH2P_)   = PhotoIon_H2

end subroutine fill_photo

subroutine calc_planet_sources(iBlock)

  use ModInputs
  use ModSources
  use ModGITM
  use ModTime
  
  implicit none

  integer, intent(in) :: iBlock

  LowAtmosRadRate = 0.0

  !\
  ! Cooling ----------------------------------------------------------
  !/


  ! None

  return

end subroutine calc_planet_sources

!---------------------------------------------------------------------
! Initialize Heating Rates
!---------------------------------------------------------------------

subroutine init_heating_efficiency

  use ModGITM, only: nLons, nLats, nAlts, nBlocks, Altitude_GB
  use ModEUV, only: HeatingEfficiency_CB, eHeatingEfficiency_CB

  implicit none

  integer :: iLon, iLat, iAlt
  !------------------------------------------------------------------

  HeatingEfficiency_CB(:,:,:,1:nBlocks) = 0.05

end subroutine init_heating_efficiency

!---------------------------------------------------------------------
! Calculate Eddy Diffusion Coefficient
!---------------------------------------------------------------------

subroutine calc_eddy_diffusion_coefficient(iBlock)

  use ModSizeGITM
  use ModGITM, only: pressure
  use ModInputs, only: EddyDiffusionPressure0,EddyDiffusionPressure1, &
       EddyDiffusionCoef
  use ModSources, only: KappaEddyDiffusion

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon

  KappaEddyDiffusion=0.

!  do iAlt = -1, nAlts+2
!
!     do iLat = 1, nLats
!        do iLon = 1, nLons
!
!           if (pressure(iLon,iLat,iAlt,iBlock) >EddyDiffusionPressure0) then
!              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = EddyDiffusionCoef
!              
!           else if (pressure(iLon,iLat,iAlt,iBlock) > &
!                EddyDiffusionPressure1) then
!
!              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = EddyDiffusionCoef * &
!                   (pressure(iLon,iLat,iAlt,iBlock) - &
!                   EddyDiffusionPressure1)/&
!                   (EddyDiffusionPressure0 - EddyDiffusionPressure1)
!
!           endif
!        enddo
!     enddo
!  enddo

end subroutine calc_eddy_diffusion_coefficient

subroutine set_planet_defaults

  use ModInputs

  return

end subroutine set_planet_defaults

subroutine planet_limited_fluxes(iBlock)

  !! Do Nothing

  integer, intent(in) :: iBlock

  return

end subroutine planet_limited_fluxes
