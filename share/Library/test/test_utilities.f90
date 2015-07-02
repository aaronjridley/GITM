!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program test_utilities

use ModUtilities, ONLY: test => test_mod_utility

implicit none

call test

end program test_utilities


subroutine CON_stop(StringError)

  implicit none

  character (len=*), intent(in) :: StringError

  write(*,'(a)') 'ERROR: '//StringError
  write(*,'(a)') 'Do not actually stop in the test'

end subroutine CON_stop
