!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program test_random

  use ModRandomNumber, ONLY: test => test_random_number

  implicit none

  call test

end program test_random

subroutine CON_stop(StringError)

  implicit none

  character (len=*), intent(in) :: StringError
  !----------------------------------------------------------------------------
  write(*,'(a)')StringError
  write(*,'(a,i3)')'!!! SWMF_ABORT !!!'
  stop

end subroutine CON_stop

