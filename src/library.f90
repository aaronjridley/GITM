!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------

subroutine report(str, iLevel)

  use ModInputs, only : iDebugLevel
  implicit none

  character (len=*), intent(in) :: str
  integer, intent(in) :: iLevel
  character (len = 11) :: cArrow
  integer :: i

  if (iDebugLevel < iLevel .or. iDebugLevel > 10) return

  do i=1,iLevel
     cArrow(i:i) = "="
  enddo
  cArrow(iLevel+1:iLevel+1) = ">"

  write(*,*) cArrow(1:iLevel+1), " ",str

end subroutine report

!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------

subroutine stop_gitm(str)

  use ModGITM
  use ModInputs, only: IsFramework
  use ModMpi
  implicit none

  character (len=*), intent(in) :: str
  integer :: ierror, erno

  if (IsFramework) then
     call CON_stop("UA/GITM Error: "//str)
  else
     write(*,*)'Stopping execution! iProc=',iProc,' with msg=',str
     call MPI_abort(iCommGITM, erno, ierror)
     stop
  endif

end subroutine stop_gitm

!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------

subroutine i2s(iValue, cOut, iLength)

  integer, intent(in)            :: iValue, iLength
  character (len=*), intent(out) :: cOut
  character (len=10)             :: cFormat
  integer                        :: i

  if (iLength < 10) then
     write(cFormat,"('(I',I1,')')") iLength
  else
     write(cFormat,"('(I',I2,')')") iLength
  endif

  write(cOut, cFormat) iValue

  do i=1, iLength
    if (cOut(i:i) == ' ') cOut(i:i) = '0'
  enddo

end subroutine i2s
