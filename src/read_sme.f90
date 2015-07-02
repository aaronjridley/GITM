!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================

subroutine read_al_onset_list(iOutputError, StartTime, EndTime)

  use ModKind
  use ModIndices
  use ModTimeConvert, ONLY: time_int_to_real
  implicit none

  integer, intent(out)     :: iOutputError
  real(Real8_), intent(in) :: EndTime, StartTime

  integer :: ierror, iAE, j, npts
  logical :: IsDone

  ! One line of input
  character (len=100) :: line

  real (Real8_) :: TimeDelay, BufferTime

  integer, dimension(7) :: itime
  !------------------------------------------------------------------------
  iOutputError = 0

  BufferTime = 24.0*3600.0 * 31.0

  IsDone = .false.

  npts = 0
  TimeDelay = 0.0

  if (NameOfSecondIndexFile == "none" .and. nIndices_V(onsetut_) == 1) return

  call init_mod_indices

  ! If we have been here before and we read the entire file, leave
  if (.not.ReReadOnsetFile .and. nIndices_V(onsetut_) > 0) return

  ! This if statement makes sure that we have been here before and 
  ! want to be here again
  if (ReReadOnsetFile) then
     ! If we still have a lot of data in memory, then don't bother
     ! reading more.
     if (StartTime + BufferTime < IndexTimes_TV(nIndices_V(onsetut_),onsetut_)) &
          return
  endif

  ! Assume that we can read the entire file
  ReReadOnsetFile = .false.

  if (nIndices_V(onsetut_) == 0) NameOfOnsetFile = NameOfSecondIndexFile

  open(LunIndices_, file=NameOfOnsetFile, status="old", iostat = ierror)

  if (ierror.ne.0) then
     iOutputError = 1
     return
  endif

  iAE = 1
  iTime = 0
  do while (.not.IsDone)
     
     read(LunIndices_,*,iostat=iError) &
          (iTime(j),j=1,5), &
          Indices_TV(iAE,onsetmlat_), &
          Indices_TV(iAE,onsetmlt_)

     iTime(6) = 0

     if (ierror /= 0) then
        IsDone = .true.

        ! This means that the GITM time is all AFTER the first 
        ! line in the file! 
        if (StartTime > IndexTimes_TV(iAE,onsetut_)) then
           iAE = iAE +1
        endif

     else

        call time_int_to_real(iTime,IndexTimes_TV(iAE,onsetut_))

        Indices_TV(iAE,onsetut_) = IndexTimes_TV(iAE,onsetut_) - StartTime

        ! This makes sure that we only store the values that we 
        ! are really interested in

        if ( IndexTimes_TV(iAE,onsetut_) >= StartTime-BufferTime .and. &
             IndexTimes_TV(iAE,onsetut_) <= EndTime+BufferTime .and. &
             iAE < MaxIndicesEntries) then

           IndexTimes_TV(iAE,onsetmlat_) = IndexTimes_TV(iAE,onsetut_)
           IndexTimes_TV(iAE,onsetmlt_) = IndexTimes_TV(iAE,onsetut_)

           iAE = iAE + 1

        else

           ! This means that the GITM time is all BEFORE the first 
           ! line in the file! 
           if (EndTime < IndexTimes_TV(iAE,onsetut_) .and. iAE == 1) then
              iAE = iAE +1
           endif

        endif

     endif

  enddo

  close(LunIndices_)

  if (iAE >= MaxIndicesEntries) ReReadOnsetFile = .true.

  nIndices_V(onsetut_)   = iAE - 1
  nIndices_V(onsetmlat_) = iAE - 1
  nIndices_V(onsetmlt_)  = iAE - 1

end subroutine read_al_onset_list

!==============================================================================

subroutine read_sme(iOutputError, StartTime, EndTime)

  use ModKind
  use ModIndices
  use ModTimeConvert, ONLY: time_int_to_real
  implicit none

  integer, intent(out)     :: iOutputError
  real(Real8_), intent(in) :: EndTime, StartTime

  integer :: ierror, iAE, j, npts
  logical :: IsDone

  ! One line of input
  character (len=100) :: line

  real (Real8_) :: TimeDelay, BufferTime = 180.0

  integer, dimension(7) :: itime
  !------------------------------------------------------------------------
  iOutputError = 0

  IsDone = .false.

  npts = 0
  TimeDelay = 0.0

  if (NameOfIndexFile == "none" .and. nIndices_V(ae_) == 1) return

  call init_mod_indices

  ! If we have been here before and we read the entire file, leave
  if (.not.ReReadSMEFile .and. nIndices_V(ae_) > 0) return

  ! This if statement makes sure that we have been here before and 
  ! want to be here again
  if (ReReadSMEFile) then
     ! If we still have a lot of data in memory, then don't bother
     ! reading more.
     if (StartTime + BufferTime < IndexTimes_TV(nIndices_V(ae_),ae_)) &
          return
  endif

  ! Assume that we can read the entire file
  ReReadSMEFile = .false.

  if (nIndices_V(ae_) == 0) NameOfSMEFile = NameOfIndexFile

  open(LunIndices_, file=NameOfSMEFile, status="old", iostat = ierror)

  if (ierror.ne.0) then
     iOutputError = 1
     return
  endif

  iAE = 1
  iTime = 0
  do while (.not.IsDone)
     
     read(LunIndices_,*,iostat=iError) &
          (iTime(j),j=1,6), &
          Indices_TV(iAE,ae_), &
          Indices_TV(iAE,al_), &
          Indices_TV(iAE,au_)

     if (ierror /= 0) then
        IsDone = .true.

        ! This means that the GITM time is all AFTER the first 
        ! line in the file! 
        if (StartTime > IndexTimes_TV(iAE,ae_)) then
           iAE = iAE +1
        endif

     else

        call time_int_to_real(iTime,IndexTimes_TV(iAE,ae_))

        ! This makes sure that we only store the values that we 
        ! are really interested in

        if ( IndexTimes_TV(iAE,ae_) >= StartTime-BufferTime .and. &
             IndexTimes_TV(iAE,ae_) <= EndTime+BufferTime .and. &
             iAE < MaxIndicesEntries) then

           iAE = iAE + 1

        else

           ! This means that the GITM time is all BEFORE the first 
           ! line in the file! 
           if (EndTime < IndexTimes_TV(iAE,ae_) .and. iAE == 1) then
              iAE = iAE +1
           endif

        endif

     endif

  enddo

  close(LunIndices_)

  if (iAE >= MaxIndicesEntries) ReReadSMEFile = .true.

  nIndices_V(ae_) = iAE - 1
  nIndices_V(au_) = iAE - 1
  nIndices_V(al_) = iAE - 1

end subroutine read_sme

