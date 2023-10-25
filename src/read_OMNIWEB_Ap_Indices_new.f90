! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!==============================================================================

subroutine read_OMNIWEB_Ap_Indices_new(iOutputError, StartTime, EndTime)

  use ModKind
  use ModIndices
  implicit none

  integer, intent(out)     :: iOutputError
  real(Real8_), intent(in) :: EndTime, StartTime

  integer :: ierror, iAp, iYear, iDay, iHour, ap
  logical :: done

  ! One line of input
  integer, parameter :: iCharLenGitm = 400
  character (len=iCharLenGitm) :: line

  real (Real8_) :: BufferTime

  integer, dimension(7) :: itime
  !------------------------------------------------------------------------
  iOutputError = 0

  done = .false.

  call init_mod_indices

  !write(*,*) "-->",NameOfIndexFile,"<--"
  
  open(LunIndices_, file=NameOfIndexFile, status="old", iostat = ierror)

  write(*,*) 'open : ', ierror
  if (ierror.ne.0) then
     iOutputError = 1
     return
  endif

  iAp = 1

  ! Allow a 48 hour buffer time before and after start/end time:
  BufferTime = 48.0 * 3600.0

  do while (.not.done)
     
     read(LunIndices_,*,iostat=iError) iYear, iDay, iHour, Ap

     if (ierror /= 0) then
        done = .true.
     else
     
        iTime(1) = iYear
        iTime(2) = 1
        iTime(3) = iDay
        iTime(4) = iHour
        iTime(5:7) = 0

        Indices_TV(iAp,ap_) = Ap

        call time_int_to_real(iTime,IndexTimes_TV(iAp,ap_))

        if ( IndexTimes_TV(iAp,ap_) >= StartTime-BufferTime .and. &
             IndexTimes_TV(iAp,ap_) <= EndTime+BufferTime .and. &
             iAp < MaxIndicesEntries) then

           iAp = iAp + 1

        else

           ! This means that the GITM time is all BEFORE the first 
           ! line in the file! 
           if (EndTime < IndexTimes_TV(iAp,ap_) .and. iAp == 1) then
              iAp = iAp +1
           endif

        endif

        if (iAp >= MaxIndicesEntries) done = 1

     endif
        
  enddo

  close(LunIndices_)

  nIndices_V(ap_) = iAp - 2

end subroutine read_OMNIWEB_Ap_Indices_new

