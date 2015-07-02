!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program test_ray

  use CON_ray_trace

  implicit none

  integer :: iError

  call MPI_init(iError)
  call ray_test
  call MPI_finalize(iError)

end program test_ray

subroutine CON_stop(StringError)

  use ModMpi
  implicit none

  character (len=*), intent(in) :: StringError

  integer :: iProc,iError,nError

  !----------------------------------------------------------------------------

  write(*,'(a)')StringError
  call MPI_COMM_RANK(MPI_COMM_WORLD, iProc, iError)
  write(*,'(a,i3)')'!!! SWMF_ABORT !!! requested by processor ',iProc
  call MPI_abort(MPI_COMM_WORLD, nError, iError)
  stop

end subroutine CON_stop

