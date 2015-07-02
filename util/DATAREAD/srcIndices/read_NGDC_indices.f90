!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================

subroutine read_NGDC_Indices(iOutputError)

  use ModKind
  use ModIndices
  use CON_geopack, ONLY: CON_recalc, GsmGse_DD, JulianDay
  use ModTimeConvert, ONLY: time_real_to_int
  implicit none

  integer, intent(out) :: iOutputError

  integer :: ierror, i, j, npts, npts_hpi, k
  integer :: input_coor_system, iday 
  logical :: done

  ! One line of input
  character (len=100) :: line

  real, dimension(6,MaxIndicesEntries) :: tmp

  real (Real8_), dimension(MaxIndicesEntries) :: ut_new, ut_keep, ut_tmp
  real, dimension(MaxIndicesEntries)   :: data_new, data_keep, data_tmp

  real (Real8_) :: time_now

  integer, dimension(7) :: itime
  !------------------------------------------------------------------------
  call init_mod_indices

  iOutputError = 0

  open(LunIndices_, file=NameOfIndexFile, status="old", iostat = ierror)

  if (ierror.ne.0) then
     iOutputError = 1
     return
  endif

  done = .false.

  npts_hpi = 0

  Input_Coor_System = GSM_

  do while (.not.done)

     read(LunIndices_,'(a)', iostat = ierror ) line
     if (ierror.ne.0) done = .true.

     if(index(line,'#Element: dst')>0)then
        call read_values
        call Insert_into_Indices_Array(tmp, dst_)
     endif

     if(index(line,'#Element: kp')>0)then
        call read_values
        call Insert_into_Indices_Array(tmp, kp_)
     endif

     if(index(line,'#Element: ap')>0)then
        call read_values
        call Insert_into_Indices_Array(tmp, ap_)
     endif

     if(index(line,'#Element: bx')>0)then
        read(LunIndices_,'(a)', iostat = ierror ) line
        if(index(line,'#Table: IMFMin')>0)then
           call read_values
           call Insert_into_Indices_Array(tmp, imf_bx_)
        endif
     endif

     if(index(line,'#Element: by')>0)then
        read(LunIndices_,'(a)', iostat = ierror ) line
        if(index(line,'#Table: IMFMin')>0)then
           call read_values
           call Insert_into_Indices_Array(tmp, imf_by_)
        endif
     endif

     if(index(line,'#Element: bz')>0)then
        read(LunIndices_,'(a)', iostat = ierror ) line
        if(index(line,'#Table: IMFMin')>0)then
           call read_values
           call Insert_into_Indices_Array(tmp, imf_bz_)
        endif
     endif

     if(index(line,'#Element: flow')>0)then
        call read_values
        call Insert_into_Indices_Array(tmp, sw_v_)
     endif

     if(index(line,'#Element: swVelX')>0)then
        call read_values
        call Insert_into_Indices_Array(tmp, sw_vx_)
     endif

     if(index(line,'#Element: swVelY')>0)then
        call read_values
        call Insert_into_Indices_Array(tmp, sw_vy_)
     endif

     if(index(line,'#Element: swVelZ')>0)then
        call read_values
        call Insert_into_Indices_Array(tmp, sw_vz_)
     endif

     if(index(line,'#Element: observed')>0)then
        call read_values
        call Insert_into_Indices_Array(tmp, f107_)
     endif

     if(index(line,'#Element: adjusted')>0)then
        call read_values
        call Insert_into_Indices_Array(tmp, f107_)
     endif

     if(index(line,'#Element: flux')>0 .or. index(line,'#Table: Flux')>0) then
        call read_values
        call Insert_into_Indices_Array(tmp, f107_)
     endif

     ! This has to be before the DMSP stuff because of the 
     ! extra N and S on the end of the Element line.

!     if (north) then
        if(index(line,'#Element: power')>0)then
!           read(LunIndices_,'(a)', iostat = ierror ) line
!           if(index(line,'#Table: Hpinoaa')>0)then
              call read_values
              call merge_hpi_data
!           endif
        endif
!     else
!        if(index(line,'#Element: power')>0)then
!           read(LunIndices_,'(a)', iostat = ierror ) line
!           if(index(line,'#Table: Hpinoaa')>0)then
!              call read_values
!              call merge_hpi_data
!           endif
!        endif
!     endif
!
!     if(index(line,'#Element: power')>0)then
!        read(LunIndices_,'(a)', iostat = ierror ) line
!        if(index(line,'#Table: Hpidmsp')>0)then
!           call read_values
!           call merge_hpi_data
!        endif
!     endif

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

  if (Input_Coor_System == GSE_) then

     do i=1,nIndices_V(imf_bx_)

        time_now = IndexTimes_TV(i,imf_bx_)
        call time_real_to_int(time_now, itime)

        iday = JulianDay(itime(1),itime(2),itime(3))
        call CON_recalc(itime(1),iday,itime(4),itime(5),itime(6),0)

        ! Convert from GSE to GSM
        Indices_TV(i,imf_bx_:imf_bz_) = &
             matmul(GsmGse_DD, Indices_TV(i,imf_bx_:imf_bz_))

     enddo

  endif

contains 

  subroutine read_values

    logical :: done_inner
    real :: missing, missingPlusTol, missingMinusTol

    done_inner = .false.

    tmp = 0.0

    missing = -1.0e32

    do while (.not.done_inner)

       read(LunIndices_,'(a)', iostat = ierror ) line
       if(index(line,'#yyyy-')>0) done_inner = .true.
       if (ierror.ne.0) done = .true.

       if (index(line,'#Missing value:')>0) read(line,'(15x,f7.1)') missing

       if (index(line,'GSE')>0) Input_Coor_System = GSE_

    enddo

    missingPlusTol  = missing + abs(missing)*0.01
    missingMinusTol = missing - abs(missing)*0.01

    if (ierror.eq.0) then

       done_inner = .false.

       i = 1

       do while (.not.done_inner)

          read(LunIndices_,'(f4.0,4(1x,f2.0),f7.1)', iostat = ierror ) tmp(1:6,i)

          if (ierror /= 0) then 
             done_inner = .true.
          else
             if (tmp(6,i) < missingMinusTol .or. &
                 tmp(6,i) > missingPlusTol) then
                i = i + 1
             else
                tmp(1:6,i) = 0.0
             endif
          endif

       enddo

       npts = i-1

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


end subroutine read_NGDC_Indices

