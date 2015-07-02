!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
function timing_cpu()
  implicit none
  integer, parameter     :: Real8_ = selected_real_kind(12)
  real(Real8_)           :: timing_cpu
  real(Real8_), external :: MPI_WTIME

  timing_cpu=MPI_WTIME()

end function timing_cpu

