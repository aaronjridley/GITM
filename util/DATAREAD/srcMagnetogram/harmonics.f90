!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program harmonics

  ! Transform raw magnetogram file into spherical harmonics file

  use ModMagHarmonics
  use ModMpi

  implicit none
  integer:: iError
  !----------------------------------------------------------------------------
  call MPI_init(iError)
  call read_harmonics_param
  call read_raw_magnetogram 
  call calc_harmonics
  call MPI_finalize(iError)

end program harmonics
!==============================================================================
subroutine CON_stop(TypeMessage)
  character(LEN=*),intent(in)::TypeMessage
  stop
end subroutine CON_stop
