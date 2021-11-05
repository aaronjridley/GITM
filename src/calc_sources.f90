! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine calc_GITM_sources(iBlock)

  use ModInputs
  use ModSources
  use ModGITM
  use ModPlanet, only: IsVenus

  implicit none

  integer, intent(in) :: iBlock

  integer :: iError

  call report("calc_GITM_sources",1)

  !\
  ! ---------------------------------------------------------------
  ! calc_rate is used to determine reaction rates, heating coefficients,
  ! ion-neutral collision frequency, lambdas for ion drag, etc.
  ! ---------------------------------------------------------------
  !/

  call calc_eddy_diffusion_coefficient(iBlock)  
  call calc_rates(iBlock)
  call calc_collisions(iBlock)

  !\
  ! ---------------------------------------------------------------
  ! These terms are for Neutral Temperature
  ! ---------------------------------------------------------------
  !/

  !\
  ! Solar Heating -------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> solar heating", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)

  if (UseSolarHeating .or. UseIonChemistry) then

     ! So far, calc_physics only has stuff that is needed for solar
     ! euv, such as solar zenith angles, and local time.

     call calc_physics(iBlock)
     call euv_ionization_heat(iBlock)

  endif

  if (.not. UseSolarHeating) EuvHeating = 0.0

  !\
  ! Planet-specific Heating ---------------------------------------
  !/

  RadCooling = 0.0
  call calc_planet_sources(iBlock)

  !\
  ! Thermal Conduction --------------------------------------------
  !/

  if (UseConduction) then
     call calc_thermal_conduction(iBlock)
  else
     Conduction = 0.0
  end if

  !\
  ! IR Heating (Venus Only) ---------------------------------------
  !/

  QnirTOT(:,:,:,:) = 0.0
  if (UseIRHeating .and. isVenus) then
     call calc_ir_heating(iBlock)
  endif
 
  !\
  ! ---------------------------------------------------------------
  ! These terms are for Neutral Winds
  ! ---------------------------------------------------------------
  !/

  !\
  ! Ion Drag ----------------------------------------------------
  !/

  if (UseIonDrag) then
     call calc_ion_drag(iBlock)
  else
     IonDrag = 0.0
     VerticalIonDrag = 0.0
  endif


  !\
  ! Viscosity ----------------------------------------------------
  !/

  if (UseViscosity) then
     call calc_viscosity(iBlock)
  else
     Viscosity = 0.0
  end if

  !\
  ! ---------------------------------------------------------------
  ! These terms are for the ionosphere
  ! ---------------------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> get_potential", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  call get_potential(iBlock)

  if (iDebugLevel > 4) write(*,*) "=====> Efield", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  call calc_efield(iBlock)

  if (iDebugLevel > 4) write(*,*) "=====> Aurora", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  call aurora(iBlock)

  if (iDebugLevel > 4) write(*,*) "=====> Ion Velocity", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  call calc_ion_v(iBlock)

  !\
  ! ---------------------------------------------------------------
  ! Chemistry
  ! ---------------------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> Chemistry", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  call calc_chemistry(iBlock)

end subroutine calc_GITM_sources
