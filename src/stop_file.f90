!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!\
! -------------------------------------------------------------------
! -------------------------------------------------------------------
!/

subroutine check_stop

  use ModGITM
  use ModTime
  use ModInputs, only: CPUTimeMax, iOutputUnit_, DoCheckStopFile
  use ModMpi

  implicit none

  real, external :: get_timing

  real*8  :: EndTimeLocal
  logical :: IsThere
  integer :: iError

  if (.not.DoCheckStopFile) return

  call report("check_stop",2)

  call start_timing("check_stop")

  EndTimeLocal = EndTime

  if (iProc == 0) then 

     inquire(file="GITM.STOP",EXIST=IsThere)
     if (IsThere) then
        write(*,*) "GITM.STOP file found. Exiting."
        EndTimeLocal = CurrentTime - 1.0
     endif

  endif

  if (get_timing("GITM") > CPUTimeMax) then
     if (iProc == 0) write(*,*) "CPUTimeMax Exceeded. Exiting."
     EndTimeLocal = CurrentTime - 1.0
     open(unit=iOutputUnit_, file="GITM.CPU", status="unknown")
     close(iOutputUnit_)
  endif

  call MPI_AllREDUCE(EndTimeLocal, EndTime,  &
       1, MPI_DOUBLE_PRECISION, MPI_MIN, iCommGITM, iError)

  call check_start

  call end_timing("check_stop")

end subroutine check_stop

!\
! -------------------------------------------------------------------
! -------------------------------------------------------------------
!/

subroutine delete_stop

  use ModGITM
  use ModInputs, only: iOutputUnit_

  implicit none

  logical :: IsThere

  call report("delete_stop",2)

  inquire(file='GITM.STOP',EXIST=IsThere)
  if (IsThere .and. iProc == 0) then
     open(iOutputUnit_, file = 'GITM.STOP', status = 'OLD')
     close(iOutputUnit_, status = 'DELETE')
  endif

  inquire(file='GITM.DONE',EXIST=IsThere)
  if (IsThere .and. iProc == 0) then
     open(iOutputUnit_, file = 'GITM.DONE', status = 'OLD')
     close(iOutputUnit_, status = 'DELETE')
  endif

  inquire(file='GITM.CPU',EXIST=IsThere)
  if (IsThere .and. iProc == 0) then
     open(iOutputUnit_, file = 'GITM.CPU', status = 'OLD')
     close(iOutputUnit_, status = 'DELETE')
  endif

end subroutine delete_stop

!\
! -------------------------------------------------------------------
! -------------------------------------------------------------------
!/

subroutine check_start

  use ModGITM
  use ModTime
  use ModInputs, only: CPUTimeMax, iOutputUnit_
  use ModMpi
  use ModUtilities, only: sleep

  implicit none

  real*8  :: EndTimeLocal
  logical :: IsThere
  integer :: iError

  if ((CurrentTime-dt) < PauseTime .and. CurrentTime > PauseTime) then

     write(*,*) "Pausing"

     IsThere = .false.

     do while (.not.IsThere) 

        inquire(file="GITM.START",EXIST=IsThere)
        if (IsThere .and. iProc == 0) then
           if (iProc == 0) then
              write(*,*) "GITM.START file found. Continuing."
              open(iOutputUnit_, file = 'GITM.START', status = 'OLD')
              close(iOutputUnit_, status = 'DELETE')
           endif
        endif

        if (.not. IsThere) call sleep(2.0)

        call MPI_BARRIER(iCommGITM,iError)

     enddo

     ! Here is where to open and read the file that will set the new pause time
     ! and update the state of GITM.

     PauseTime = PauseTime + 300.0  ! delete this when you read the new file....

  endif

end subroutine check_start

