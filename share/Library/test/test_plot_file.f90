!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program test_plot_file

  use ModPlotFile, ONLY: test => test_plot_file

  implicit none

  call test

end program test_plot_file

subroutine CON_stop(StringError)

  implicit none

  character (len=*), intent(in) :: StringError

  write(*,'(a)') 'ERROR: '//StringError

  stop

end subroutine CON_stop
