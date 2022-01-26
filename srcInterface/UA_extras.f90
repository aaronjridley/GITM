!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!\
! ----------------------------------------------------------------------
!/

!! CANDIDATE FOR REMOVAL.

subroutine IE_set_inputs(Lines)

  use ModUtilities, ONLY: CON_stop
  
  implicit none

  character (len=100), dimension(100), intent(in) :: Lines

  call CON_stop('we should not use UA_extras::IE_set_inputs')
  
  return

end subroutine IE_set_inputs

!\
! ----------------------------------------------------------------------
!/

subroutine IE_Initialize(iError)

  use ModUtilities, ONLY: CON_stop
  
  implicit none

  integer, intent(out) :: iError

  iError = 0

  call CON_stop('we should not use UA_extras::IE_initialize')
  
  return

end subroutine IE_Initialize
