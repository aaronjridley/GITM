! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine calc_GITM_sources(iBlock)

  use ModInputs
  use ModSources
  use ModGITM
  use ModPlanet, only: IsVenus

  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iError, iDir, iLat, iLon, iSpecies, iIon
  integer :: iiAlt, iiLat, iiLon, iSza, jAlt
  real :: tmp(nLons, nLats, nAlts)
  real :: tmp3(nLons, nLats, 0:nAlts+1)
  real :: RhoI(nLons, nLats, nAlts)
  real :: Rho110(nLons, nLats, 0:nAlts+1)
  real :: Weight(nLons, nLats, nAlts)
  real :: ScaleHeight(-1:nLons+2, -1:nLats+2, -1:nAlts+2)

  real :: TmpDiff(nLons, nLats, 0:nAlts+1)
  real :: TmpMulFac(nLons, nLats, 0:nAlts+1)
  
  ! Sub-timestep used in the multi-step conduction
  real :: DtCSLocal
  ! NeuBCS used for Conduction and Bulk Winds
  logical ::NeuBCS                      ! Set to true if you want Neumann BCs

  ! Note:  NeuBCs is a setting for Neumann Boundary Conditions.
  ! If you are not specifying a Jeans Escape at the top, then
  ! NeuBCS = .true.
  ! 
  ! If you are setting a Jeans escape (i.e. for H, H2, or He)
  ! then NeuBCS MUST BE SET TO FALSE FOR THAT SPECIES
  ! Otherwise, you will not get the escape that you want!

  call report("calc_GITM_sources",1)

  ChemicalHeatingRate = 0.0
  ChemicalHeatingRateIon = 0.0
  ChemicalHeatingRateEle = 0.0

  ChemicalHeatingSpecies = 0.0

  ! calc_rate is used to determine reaction rates, heating coefficients,
  ! ion-neutral collision frequency, lambdas for ion drag, etc.

  call calc_eddy_diffusion_coefficient(iBlock)  
  call calc_rates(iBlock)
  call calc_collisions(iBlock)

  RhoI = IDensityS(1:nLons,1:nLats,1:nAlts,ie_,iBlock) * &
       MeanIonMass(1:nLons,1:nLats,1:nAlts)

  Rho110 = Rho(1:nLons, 1:nLats,0:nAlts+1, iBlock)
  
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

  ! This includes Radiative Cooling....
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
  ! IR Heating (Venus Only) --------------------------------------------
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
  ! These terms are for ionosphere
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
  if (UseAuroralHeating) then
     AuroralHeating = AuroralHeatingRate(:,:,:,iBlock) / &
          TempUnit(1:nLons,1:nLats,1:nAlts) / cp(:,:,1:nAlts,iBlock) / &
          rho(1:nLons,1:nLats,1:nAlts, iBlock)
  else
     AuroralHeating = 0.0
  endif

  if (iDebugLevel > 4) write(*,*) "=====> Ion Velocity", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  call calc_ion_v(iBlock)

  if (iDebugLevel > 4) write(*,*) "=====> Chemistry", iproc
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  ! The Emissions array was never set. Should this be here or earlier ????
  Emissions(:,:,:,:,iBlock) = 0.0
  call calc_chemistry(iBlock)

  ChemicalHeatingRate(:,:,:) = &
       ChemicalHeatingRate(:,:,:) * Element_Charge / &
       TempUnit(1:nLons,1:nLats,1:nAlts) / cp(1:nLons,1:nLats,1:nAlts,iBlock)/&
       rho(1:nLons,1:nLats,1:nAlts,iBlock)
	   
  ChemicalHeatingRateIon(:,:,:) = &
       ChemicalHeatingRateIon(:,:,:) * Element_Charge

  ChemicalHeatingSpecies = ChemicalHeatingSpecies * Element_Charge

end subroutine calc_GITM_sources
