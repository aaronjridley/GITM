!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
! This routine reads two files from the ACE science center:
! 1. IMF Bx, By and Bz (date is in year DOY hour minute second)
! 2. SWEPAM H+ density and Vx (date is in year DOY hour minute second)
!

subroutine read_ACE_data(iOutputError)

  use ModKind
  use ModIndices
  use ModTimeConvert, ONLY: time_int_to_real
  implicit none

  integer, intent(out) :: iOutputError

  integer :: ierror, iIMF, iSW, j, npts
  logical :: done, done_inner

  ! One line of input
  character (len=100) :: line

  real (Real8_) :: TimeDelay

  integer, dimension(7) :: itime
  integer :: ydhm(4)
  real :: second, avevx

  !------------------------------------------------------------------------
  iOutputError = 0

  call init_mod_indices

  TimeDelay=3600.0
  iSW = 1

  write(*,*) "ACE Science Center Data Reader"

  ierror = 0

  open(LunIndices_, file=NameOfIndexFile, status="old", iostat = ierror)

  if (ierror.ne.0) then
     iOutputError = 1
     return
  endif

  done = .false.

  do while (.not.done)
     
     read(LunIndices_,'(a)', iostat = ierror ) line
     if (ierror.ne.0) done = .true.
     if (index(line,'BEGIN DATA')>0) done = .true.

  enddo

  if (ierror == 0) then

     done = .false.

     iIMF = 1

     do while (.not.done)

        read(LunIndices_,*,iostat=iError) &
             (ydhm(j),j=1,4), second, &
             Indices_TV(iIMF,imf_bx_), &
             Indices_TV(iIMF,imf_by_), &
             Indices_TV(iIMF,imf_bz_)

        itime(1) = ydhm(1)
        itime(2) = 1
        itime(3) = ydhm(2)
        itime(4) = ydhm(3)
        itime(5) = ydhm(4)
        itime(6) = second
        itime(7) = 0

        if (ierror /= 0) then
           done = .true.
        else

           if (Indices_TV(iIMF,imf_bx_) > -999.0) then 

              call time_int_to_real(iTime,IndexTimes_TV(iIMF,imf_bx_))

              IndexTimes_TV(iIMF,imf_by_) = IndexTimes_TV(iIMF,imf_bx_)
              IndexTimes_TV(iIMF,imf_bz_) = IndexTimes_TV(iIMF,imf_bx_)

              iIMF = iIMF + 1

           endif

        endif

     enddo

  endif

  close(LunIndices_)

  open(LunIndices_, file=NameOfSecondIndexFile, status="old", iostat = ierror)

  if (ierror.ne.0) then
     iOutputError = 1
     return
  endif

  done = .false.

  do while (.not.done)
     
     read(LunIndices_,'(a)', iostat = ierror ) line
     if (ierror.ne.0) done = .true.
     if(index(line,'BEGIN DATA')>0) done = .true.

  enddo

  if (ierror == 0) then

     done = .false.

     iSW = 1

     do while (.not.done)

        read(LunIndices_,*,iostat=iError) &
             (ydhm(j),j=1,4), second, &
             Indices_TV(iSW,sw_n_), &
             Indices_TV(iSW,sw_vx_)

        if (ierror /= 0) then
           done = .true.
        else

           if ((Indices_TV(iSW,sw_n_) > 0.0) .and. &
               (Indices_TV(iSW,sw_vx_) > -9999.0)) then

              Indices_TV(iSW,sw_vy_) = 0.0
              Indices_TV(iSW,sw_vz_) = 0.0

              Indices_TV(iSW,sw_v_) = sqrt( &
                   Indices_TV(iSW,sw_vx_)**2 + &
                   Indices_TV(iSW,sw_vy_)**2 + &
                   Indices_TV(iSW,sw_vz_)**2)

              itime(1) = ydhm(1)
              itime(2) = 1
              itime(3) = ydhm(2)
              itime(4) = ydhm(3)
              itime(5) = ydhm(4)
              itime(6) = second
              itime(7) = 0

              call time_int_to_real(iTime,IndexTimes_TV(iSW,sw_vx_))

              IndexTimes_TV(iSW,sw_v_)  = IndexTimes_TV(iSW,sw_v_)
              IndexTimes_TV(iSW,sw_vy_) = IndexTimes_TV(iSW,sw_vx_)
              IndexTimes_TV(iSW,sw_vz_) = IndexTimes_TV(iSW,sw_vx_)
              IndexTimes_TV(iSW,sw_n_)  = IndexTimes_TV(iSW,sw_vx_)
              
              iSW = iSW + 1

           endif

        endif

     enddo

  endif

  close(LunIndices_)

  nIndices_V(imf_bx_) = iIMF - 1
  nIndices_V(imf_by_) = iIMF - 1
  nIndices_V(imf_bz_) = iIMF - 1
  nIndices_V(sw_vx_)  = iSW - 1
  nIndices_V(sw_vy_)  = iSW - 1
  nIndices_V(sw_vz_)  = iSW - 1
  nIndices_V(sw_v_)   = iSW - 1
  nIndices_V(sw_n_)   = iSW - 1
  nIndices_V(sw_t_)   = iSW - 1

  avevx = 0.0

  do iSW = 1, nIndices_V(sw_vx_)
     avevx = avevx + Indices_TV(iSW,sw_v_)/nIndices_V(sw_vx_)
  enddo

  TimeDelay = 1.5e6 / avevx

  write(*,*) "Time Delay : ", TimeDelay/60.0, ' minutes'

  do iSW = 1, nIndices_V(sw_vx_)
     IndexTimes_TV(iSW,sw_vx_) = IndexTimes_TV(iSW,sw_vx_) &
          + TimeDelay
     IndexTimes_TV(iSW,sw_v_)  = IndexTimes_TV(iSW,sw_v_) &
          + TimeDelay
     IndexTimes_TV(iSW,sw_vy_) = IndexTimes_TV(iSW,sw_vy_) &
          + TimeDelay
     IndexTimes_TV(iSW,sw_vz_) = IndexTimes_TV(iSW,sw_vz_) &
          + TimeDelay
     IndexTimes_TV(iSW,sw_n_)  = IndexTimes_TV(iSW,sw_n_) &
          + TimeDelay
     IndexTimes_TV(iSW,sw_t_)  = IndexTimes_TV(iSW,sw_t_) &
          + TimeDelay
  enddo

  do iIMF = 1, nIndices_V(imf_bx_)
     IndexTimes_TV(iIMF,imf_bx_) = IndexTimes_TV(iIMF,imf_bx_)+&
          TimeDelay
     IndexTimes_TV(iIMF,imf_by_) = IndexTimes_TV(iIMF,imf_by_)+&
          TimeDelay
     IndexTimes_TV(iIMF,imf_bz_) = IndexTimes_TV(iIMF,imf_bz_)+&
          TimeDelay
  enddo


end subroutine read_ACE_data

