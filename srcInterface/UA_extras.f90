!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!\
! ----------------------------------------------------------------------
!/

subroutine IE_set_inputs(Lines)

  implicit none

  character (len=100), dimension(100), intent(in) :: Lines

  return

end subroutine IE_set_inputs

!\
! ----------------------------------------------------------------------
!/

subroutine IE_Initialize(iError)

  implicit none

  integer, intent(out) :: iError

  iError = 0

  return

end subroutine IE_Initialize
