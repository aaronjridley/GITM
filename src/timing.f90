
subroutine start_timing(cTimingNameIn)

  use ModTimingGITM
  use ModMpi
  implicit none

  character (LEN=*), intent(in) :: cTimingNameIn
  integer :: i, iTiming, iLevel
  logical :: IsFirstTime = .true.

  if (IsFirstTime) then
     IsFirstTime = .false.
     Timings = 0.0
     iTimingLevel = 0
     IsTiming = .false.
  endif

  iLevel = 0
  iTiming = -1
  do i=1,nTiming
     if (IsTiming(i)) iLevel = iLevel + 1
     if (index(cTimingNames(i),cTimingNameIn) > 0) then
        iTiming = i
     endif
  enddo

  if (iTiming < 0) then
     nTiming = nTiming + 1
     cTimingNames(nTiming) = cTimingNameIn
     iTiming = nTiming
  endif

  IsTiming(iTiming) = .true.
  if (iTimingLevel(iTiming) < iLevel) iTimingLevel(iTiming) = iLevel
  Timings(iTiming,1) = mpi_wtime()

end subroutine start_timing

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

subroutine end_timing(cTimingNameIn)

  use ModTimingGITM
  use ModMpi
  implicit none

  character (LEN=*), intent(in) :: cTimingNameIn
  integer :: i, iTiming

  iTiming = -1
  do i=1,nTiming
     if (index(cTimingNames(i),cTimingNameIn) > 0) then
        iTiming = i
     endif
  enddo

  if (iTiming < 0) then
     write(*,*) "Unknown timing : ",cTimingNameIn
     return
  endif

  IsTiming(iTiming) = .false.
  Timings(iTiming,2) = Timings(iTiming,2) + &
       (mpi_wtime() - Timings(iTiming,1))

end subroutine end_timing

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

real function get_timing(cTimingNameIn)

  use ModTimingGITM
  use ModMpi
  implicit none

  character (LEN=*), intent(in) :: cTimingNameIn
  integer :: i, iTiming

  iTiming = -1
  do i=1,nTiming
     if (index(cTimingNames(i),cTimingNameIn) > 0) then
        iTiming = i
     endif
  enddo

  if (iTiming < 0) then
     write(*,*) "Error!!"
     write(*,*) "Unknown timing : ",cTimingNameIn
     get_timing = 0.0
  else
     if (IsTiming(iTiming)) then
        get_timing = mpi_wtime() - Timings(iTiming,1)
     else
        get_timing = Timings(iTiming,2)
     endif
  endif

end function get_timing

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

subroutine report_timing(cTimingNameIn)

  use ModTimingGITM
  use ModMpi
  implicit none

  character (LEN=*), intent(in) :: cTimingNameIn
  character (len=10) :: cLevel, clevelSpace
  integer :: i, iTiming

  cLevel = "----------"
  cLevelSpace = "          "

  if (index(cTimingNameIn,'all') > 0) then
     do iTiming=1,nTiming
        write(*,"(a,a,a,a,a,f7.1,a)") &
             clevel(1:iTimingLevel(iTiming)+1),"> ",&
             cTimingNames(iTiming),&
             clevelSpace(1:5-iTimingLevel(iTiming)),&
             " took ",Timings(iTiming,2),&
             " seconds to complete"
     enddo
     return
  endif

  iTiming = -1
  do i=1,nTiming
     if (index(cTimingNames(i),cTimingNameIn) > 0) then
        iTiming = i
     endif
  enddo

  if (iTiming < 0) then
     write(*,*) "Unknown timing : ",cTimingNameIn
     return
  endif

  write(*,"(a,a,f7.1,a)") &
          cTimingNames(iTiming)," took ",Timings(iTiming,2),&
          " seconds to complete", iTimingLevel(iTiming)

end subroutine report_timing


