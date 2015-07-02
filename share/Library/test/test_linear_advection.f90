!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program linear_advection_test

  use ModLinearAdvection, ONLY: test_linear_advection

  implicit none

  call test_linear_advection

end program linear_advection_test

subroutine CON_stop(String)
  implicit none
  character(len=*), intent(in):: String
  write(*,*)'ERROR: ',String
  stop
end subroutine CON_stop
