!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program test_freq

  use ModFreq, ONLY: test => test_freq

  implicit none

  call test

end program test_freq

subroutine CON_stop(StringError)

  implicit none

  character (len=*), intent(in) :: StringError

  write(*,'(a)')StringError

  stop

end subroutine CON_stop
