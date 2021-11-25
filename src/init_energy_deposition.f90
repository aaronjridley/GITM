! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine init_energy_deposition

  use ModInputs, only: iDebugLevel
  use ModGITM, only: Altitude_GB
  use ModSources

  implicit none

  integer :: ierr, iAlt, i, iEnergy
  real :: a
  logical :: IsDone

  !--------------------------------------
  ! Start doing energy deposition stuff
  !--------------------------------------

  call report("init_energy_deposition",2)

  if (iDebugLevel > 2) write(*,*) "===> ED_Init"
  call ED_Init(ierr)

  if (ierr /= 0) then
     call stop_gitm("Error in initilizing the energy deposition tables")
  endif

  if (iDebugLevel > 2) write(*,*) "===> ED_Get_Grid_Size"
  call ED_Get_Grid_Size(ED_N_Alts)
  allocate(ED_grid(ED_N_Alts), stat=ierr)
  if (ierr /= 0) then
     call stop_gitm("Error allocating array ED_grid")
  endif

  allocate(ED_Ion(ED_N_Alts), stat=ierr)
  if (ierr /= 0) then
     call stop_gitm("Error allocating array ED_Ion")
  endif

  allocate(ED_Heating(ED_N_Alts), stat=ierr)
  if (ierr /= 0) then
     call stop_gitm("Error allocating array ED_Heating")
  endif

  if (iDebugLevel > 2) write(*,*) "===> ED_Get_Number_of_Energies"
  call ED_Get_Number_of_Energies(ED_N_Energies)
  allocate(ED_Energies(ED_N_Energies), stat=ierr)
  allocate(ED_Energy_edges(ED_N_Energies+1), stat=ierr)
  allocate(ED_delta_energy(ED_N_Energies), stat=ierr)
  if (ierr /= 0) then
     call stop_gitm("Error allocating array ED_Energies")
  endif

  allocate(ED_Flux(ED_N_Energies), stat=ierr)
  allocate(ED_EnergyFlux(ED_N_Energies), stat=ierr)
  allocate(ED_Ion_Flux(ED_N_Energies), stat=ierr)
  allocate(ED_Ion_EnergyFlux(ED_N_Energies), stat=ierr)
  if (ierr /= 0) then
     call stop_gitm("Error allocating array ED_Flux")
  endif

  if (iDebugLevel > 2) write(*,*) "===> ED_Get_Grid"
  call ED_Get_Grid(ED_grid, .true., ierr)

  do iAlt = 1, ED_N_Alts
     ED_grid(iAlt) = alog(ED_grid(iAlt))
  enddo

  call ED_Get_Energies(ED_Energies)

  ! Energies are from biggest to smallest.
  ! First is an appoximation:
  ED_delta_energy(1) = ED_Energies(1) - ED_Energies(2)
  ED_Energy_edges(1) = ED_Energies(1) + ED_delta_energy(1)/2.0
  iEnergy = 1
  do iEnergy = 2, ED_N_Energies
     ED_Energy_edges(iEnergy) = &
          (ED_Energies(iEnergy-1) + ED_Energies(iEnergy))/2
     ED_delta_energy(iEnergy-1) = &
          ED_Energy_edges(iEnergy-1) - ED_Energy_edges(iEnergy)
  enddo

  iEnergy = ED_N_Energies
  ED_delta_energy(iEnergy) = ED_Energies(iEnergy-1) - ED_Energies(iEnergy)

  iEnergy = ED_N_Energies+1
  ED_Energy_edges(iEnergy) = ED_Energies(iEnergy-1) - ED_delta_energy(iEnergy-1)/2
     
end subroutine init_energy_deposition
