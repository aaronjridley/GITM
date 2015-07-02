!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================

subroutine read_NOAAHPI_Indices_new(iOutputError,StartTime,EndTime)

  use ModKind
  use ModIndices
  use ModTimeConvert

  implicit none

  integer, intent(out) :: iOutputError
  real(Real8_), intent(in) :: EndTime, StartTime

  integer :: ierror, i, j, npts, npts_hpi, k
  integer :: datatype
  logical :: done

  ! One line of input
  character (len=100) :: line

  real, dimension(6,MaxIndicesEntries) :: tmp

  real (Real8_), dimension(MaxIndicesEntries) :: ut_new, ut_keep, ut_tmp
  real, dimension(MaxIndicesEntries)   :: data_new, data_keep, data_tmp

  real (Real8_) :: BufferTime = 3600.0, time_now

  integer, dimension(7) :: itime

  iOutputError = 0
  !-------------------------------------------------------------------------

  ! Problem - we have to figure out whether we are being called with no IMF
  ! file at all or whether we actually have an index file to read (or reread)
  ! Can't change comparability, so we can only use certain information.

  ! We know that nIndices_V(hpi_) = 1 if there is a single index.
  ! We know that NameOfIMFFile=none if (a) we have never read the file or (b)
  !    we have no index file.

  if (nIndices_V(hpi_) == 0 .and. NameOfHPIFile == "none") &
       NameOfHPIFile = NameOfIndexFile

  if (NameOfHPIFile == "none") return
  if (nIndices_V(hpi_) == 1) return

  call init_mod_indices


  ! If we have been here before and we read the entire file, leave
  if (.not.ReReadHPIFile .and. nIndices_V(hpi_) > 0) return

  ! This if statement makes sure that we have been here before and 
  ! want to be here again
  if (ReReadHPIFile) then
     ! If we still have a lot of data in memory, then don't bother
     ! reading more.
     if (StartTime + BufferTime < IndexTimes_TV(nIndices_V(hpi_),hpi_)) return
  endif

  ! Assume that we can read the entire file
  ReReadIMFFile = .false.

  if (nIndices_V(hpi_) == 0) NameOfHPIFile = NameOfIndexFile

  open(LunIndices_, file=NameOfHPIFile, status="old", iostat = ierror)

  if (ierror.ne.0) then
     iOutputError = 1
     return
  endif

  done = .false.

  npts_hpi = 0

  do while (.not.done)

     read(LunIndices_,'(a)', iostat = ierror ) line
     if (ierror.ne.0) done = .true.

     if(index(line,'Normalizing factor')>0)then
        if(index(line,'F8.3')>0) datatype = 1
        if(index(line,'F7.2')>0) datatype = 2
        
        call read_values
        call merge_hpi_data
     endif

  enddo

  close(LunIndices_)

  if (npts_hpi > 0) then

     tmp = 0.0

     do i=1,npts_hpi
        call time_real_to_int(ut_keep(i), itime)
        tmp(1:5,i) = itime(1:5)
        tmp(6,i) = data_keep(i)
     enddo

     call Insert_into_Indices_Array(tmp, hpi_)

     nIndices_V(hpi_norm_) = nIndices_V(hpi_) 

     do i=1,nIndices_V(hpi_norm_)

        IndexTimes_TV(i,hpi_norm_) = IndexTimes_TV(i,hpi_) 

        if (Indices_TV(i,hpi_) > 0) then
           Indices_TV(i,hpi_norm_) = &
                2.09 * ALOG(Indices_TV(i,hpi_)) * 1.0475
        endif

     enddo

  endif

contains 

  subroutine read_values

    logical :: done_inner, IsFirstTime
    real :: missing
    integer :: iYear

    done_inner = .false.

    tmp = 0.0

    missing = -1.0e32

    read(LunIndices_,'(a)', iostat = ierror ) line

    if (ierror.eq.0) then

       done_inner = .false.

       i = 1 
       IsFirstTime = .true.
      
       do while (.not.done_inner)

          if(datatype.eq.1) then
             ! OLD NOAA HPI FILES
             if (IsFirstTime) then 
                read(LunIndices_,'(i4)', iostat = ierror ) iYear
                IsFirstTime = .false.
             endif
             tmp(1,i) = iYear
             tmp(2,i) = 1
             read(LunIndices_,'(a10,f3.0,f2.0,f2.0,f8.1)', iostat = ierror ) &
                  line,tmp(3:6,i)
             if (tmp(3,i) < 1) iError = 1
          endif

          if(datatype.eq.2) then
             ! NEW NOAA HPI FILES
             read(LunIndices_,'(f4.0,a1,f2.0,a1,f2.0,a1,f2.0,a1,f2.0,a1,f2.0,a15,f8.1)', &
                  iostat = ierror ) &
                  tmp(1,i),line,tmp(2,i),line,tmp(3,i),line,tmp(4,i),line,tmp(5,i),line,tmp(6,i), &
                  line,tmp(6,i)

          endif

          if (ierror /= 0) then
             done_inner = .true.
             ! This means that the last line in the file is before the 
             ! first time in GITM.
             if (StartTime > time_now) i = i + 1
          else
             itime(1:5) = tmp(1:5,i)
             itime(6:7) = 0
             call time_int_to_real(itime, time_now)
             if ( time_now >= StartTime-BufferTime .and. &
                  time_now <= EndTime+BufferTime*2 .and. &
                  i < MaxIndicesEntries) then 
                i = i + 1
             else
                ! This means that the GITM time is all BEFORE the first 
                ! line in the file! 
                if (EndTime < time_now .and. i == 1) i = i + 1
             endif

          endif

       enddo

       npts = i-1
       if (npts+1 >= MaxIndicesEntries) ReReadHPIFile = .true.

    end if

  end subroutine read_values

  subroutine merge_hpi_data

    use ModTimeConvert, ONLY: time_int_to_real

    itime = 0

    do i=1,npts 
       itime(1:5) = tmp(1:5,i)
       call time_int_to_real(itime, ut_new(i))
       data_new(i) = tmp(6,i)
    enddo

    if (npts_hpi == 0) then
       npts_hpi = npts
       ut_keep(1:npts) = ut_new(1:npts)
       data_keep(1:npts) = data_new(1:npts)
    else
              
       ut_tmp = ut_keep
       data_tmp = data_keep

       ut_keep = 0.0
       data_keep = 0.0

       j = 1
       i = 1
       k = 1

       do while (i <= npts .or. j <= npts_hpi)

          if (i > npts) then
             ut_keep(k) = ut_tmp(j) 
             data_keep(k) = data_tmp(j) 
             k = k + 1
             j = j + 1
          else
             if (j > npts_hpi) then
                ut_keep(k) = ut_new(i) 
                data_keep(k) = data_new(i) 
                k = k + 1
                i = i + 1
             else
                if (ut_tmp(j) < ut_new(i)) then
                   ut_keep(k) = ut_tmp(j) 
                   data_keep(k) = data_tmp(j) 
                   k = k + 1
                   j = j + 1
                else
                   ut_keep(k) = ut_new(i) 
                   data_keep(k) = data_new(i) 
                   k = k + 1
                   i = i + 1
                endif
             endif
          endif

       enddo

       npts_hpi = npts_hpi + npts

    endif

  end subroutine merge_hpi_data


end subroutine read_NOAAHPI_Indices_new

