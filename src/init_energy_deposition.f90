!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine init_energy_deposition

  use ModInputs, only: iDebugLevel
  use ModGITM, only: Altitude_GB
  use ModSources

  implicit none

  integer :: ierr, iAlt, i
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
  if (ierr /= 0) then
     call stop_gitm("Error allocating array ED_Energies")
  endif

  allocate(ED_Flux(ED_N_Energies), stat=ierr)
  if (ierr /= 0) then
     call stop_gitm("Error allocating array ED_Flux")
  endif

  if (iDebugLevel > 2) write(*,*) "===> ED_Get_Grid"
  call ED_Get_Grid(ED_grid, .true., ierr)

  do iAlt = 1, ED_N_Alts
     ED_grid(iAlt) = alog(ED_grid(iAlt))
  enddo

  call ED_Get_Energies(ED_Energies)

end subroutine init_energy_deposition
