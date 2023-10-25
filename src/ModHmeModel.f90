module ModHmeModel

  use ModSizeGitm, only: nLons, nLats
  
  implicit none

  integer, parameter :: iCharLenHme_ = 10
  integer, parameter :: iCharLenFile_ = 300
  integer, parameter :: iCharLenLong_ = 1500

  character (len=iCharLenFile_) :: hmeDir = 'UA/DataIn/Hme/fitted_tide_fortran/'
  character (len=iCharLenHme_), dimension(3) :: altn

  integer, parameter :: nAltHme = 3    ! alt = 95.0, 97.5, 100 km
  integer, parameter :: nLatsHme = 59
  integer, parameter :: nHme = 49
  integer, parameter :: nDayHme = 37

  real, dimension(nHme, nLatsHme, nAltHme) :: StorageU_amp, StorageU_phase
  real, dimension(nHme, nLatsHme, nAltHme) :: StorageV_amp, StorageV_phase
  real, dimension(nHme, nLatsHme, nAltHme) :: StorageT_amp, StorageT_phase
  real, dimension(nHme, nLatsHme, nAltHme) :: StorageG_amp, StorageG_phase
  real, dimension(nHme, nLatsHme, nAltHme) :: StorageR_amp, StorageR_phase
  logical :: isReadInHme = .false.
  
  integer :: day
  integer, dimension(nDayHme) :: daysTidi = (/(day,day=1,361,10)/) 

  logical :: ReadFiles = .true.
  
  ! tidi fitted coef
  real, dimension(nHme)  :: amp_fac_arr,pha_shift_arr
  character(len=iCharLenHme_), dimension(nHme) :: hme_arr
  real, dimension(nHme,2) :: Storage_Tidi_amp
  real, dimension(nHme,2) :: Storage_Tidi_phase

  character (len=iCharLenHme_), dimension(:),allocatable :: hmes_slt
  real, dimension(:),allocatable :: amps_slt,phas_slt
  logical, dimension(nHme)  :: lpHme_slt

  ! For construct global map, lon can be change
  real :: omega = 360.0/24
  real :: pi = 3.141592
  
  ! change this based on GITM
  integer, parameter :: nLonHme = nLons + 4  ! add gost cells
  real, dimension(nLonHme) :: lon
  integer :: iday = -1
  real :: ut  ! in hours

  real, dimension(nLatsHme)  :: hmeLats
  real, dimension(nLonHme,nLatsHme,nAltHme) :: &
       u_3d, v_3d, geopt_3d, temp_3d, dr_rho_3d

  real, dimension(nAltHme)  :: hmeAlts
  
  integer :: iAlt_Hme(-1:0)
  real :: rAlt_Hme(-1:0)

  integer :: iLat_Hme(-1:nLats+2)
  real :: rLat_Hme(-1:nLats+2)
    
contains

  ! -----------------------------------------------------------------
  ! Set the longitudes to get
  ! -----------------------------------------------------------------

  subroutine set_hme_longitudes(lonsInRad)

    real, dimension(nLons+4) :: lonsInRad
    lon = lonsInRad * 180.0 / 3.141592
        
  end subroutine set_hme_longitudes
  
  ! -----------------------------------------------------------------
  ! Read inputs of HMEs (selected HMEs) 
  ! -----------------------------------------------------------------
  subroutine read_HmeInputs(NameOfFile)
   
    implicit none
    integer :: LunIndices_, ierror,iOutputError
    character (len=iCharLenFile_),intent(in) :: NameOfFile
       
    integer :: i
   
    LunIndices_ = 100
   
    open(LunIndices_, file=NameOfFile, status="old", iostat = ierror)
   
    if (ierror.ne.0) then
       iOutputError = 1
       write(*,*) 'No Hme Input file was found'
       return
    endif
   
    do i = 1,nHme
       read(LunIndices_,*) hme_arr(i),lpHme_slt(i)
    enddo
   
    close(LunIndices_)
   
  end subroutine read_HmeInputs
   
  ! -----------------------------------------------------------------
  ! Read HME file
  ! -----------------------------------------------------------------

  subroutine read_hme_file(NameOfIndexFile,u_a,u_p,v_a,v_p,&
           geopt_a,geopt_p,temp_a,temp_p,dr_rho_a,dr_rho_p)
   
    implicit none
   
    integer :: LunIndices_, ierror,iOutputError
    character (len=iCharLenFile_),intent(in) :: NameOfIndexFile
    character (len=iCharLenLong_) :: line
    
    logical :: done
    integer :: ipt
    real, dimension(11,nLatsHme) :: tmp
    real, dimension(nLatsHme)  :: u_a,u_p,v_a,v_p,geopt_a,geopt_p,&
                                 temp_a,temp_p,dr_rho_a,dr_rho_p
   
    LunIndices_ = 100
   
    open(LunIndices_, file=NameOfIndexFile, status="old", iostat = ierror)
   
    if (ierror.ne.0) then
       iOutputError = 1
       write(*,*) 'No HME file was found. Please check the directory of data'
       return
    endif
   
    done = .false.
    
    ipt = 0
    do while (.not.done)
   
       if (ipt<1) then 
           read(LunIndices_,'(a)', iostat = ierror ) line
           if (ierror.ne.0) done = .true.
       else
           
           done = .true.
       end if
   
       ipt = ipt + 1
   
    enddo
   
    read(LunIndices_,*) tmp
   
    close(LunIndices_)
    
    !['lat', 'u_a', 'u_p', 'v_a', 'v_p', 'geopt_a', 'geopt_p', 
    !  'temp_a', 'temp_p', 'dr_rho_a', 'dr_rho_p']
    hmeLats = tmp(1,:)
    u_a   = tmp(2,:)
    u_p   = tmp(3,:)
    v_a   = tmp(4,:)
    v_p   = tmp(5,:)
    geopt_a = tmp(6,:)
    geopt_p = tmp(7,:)
    temp_a  = tmp(8,:)
    temp_p  = tmp(9,:)
    dr_rho_a = tmp(10,:)
    dr_rho_p = tmp(11,:)
   
  end subroutine read_hme_file
   
  ! -----------------------------------------------------------------
  ! Read tidi fitted coefficients 
  ! -----------------------------------------------------------------
  subroutine read_coef_file(NameOfIndexFile)
   
    implicit none
    integer :: LunIndices_, ierror,iOutputError
    character (len=iCharLenFile_),intent(in) :: NameOfIndexFile
    character (len=iCharLenLong_) :: line
       
    logical :: done
    integer :: ipt
       
    integer :: i
   
    LunIndices_ = 100
   
    open(LunIndices_, file=NameOfIndexFile, status="old", iostat = ierror)
   
    if (ierror.ne.0) then
       iOutputError = 1
       write(*,*) 'No tidi coef file found. Please check the directory of data'
       return
    endif
   
    done = .false.
       
    ipt = 0
    do while (.not.done)
   
       if (ipt<1) then 
          read(LunIndices_,'(a)', iostat = ierror ) line
          if (ierror.ne.0) done = .true.
       else
              
          done = .true.
       end if
   
       ipt = ipt + 1
   
    enddo
       
    do i = 1,nHme
       read(LunIndices_,*) hme_arr(i),amp_fac_arr(i),pha_shift_arr(i)
    enddo
   
    close(LunIndices_)
   
  end subroutine read_coef_file
   
  ! -----------------------------------------------------------------
  ! Get n and s from name of HME  
  ! -----------------------------------------------------------------
  subroutine get_ns(hme_tmp,n,s) 
       
    implicit none
    integer :: n, s,tmp
    character (len=iCharLenHme_),intent(in) :: hme_tmp
       
    if (hme_tmp(1:1) .eq. 'D') then
       n = 1
    else if (hme_tmp(1:1) .eq. 'S') then
       n = 2
    else if (hme_tmp(1:1) .eq. 'T') then
       n = 3
    else
       write(*,*) 'Please check hme_tmp', hme_tmp
       call stop_gitm("I am unsure of what to do now! Stopping!")
    end if
   
    if (hme_tmp(2:2) .eq. 'W') then
       read(hme_tmp(3:3) , *) s 
   
    else if (hme_tmp(2:2) .eq. 'E') then
       read(hme_tmp(3:3) , *) tmp
       s = -tmp
   
    else if (hme_tmp(2:2) .eq. '0') then
       s = 0
    else
       write(*,*) 'Plese check hme_tmp',hme_tmp
       call stop_gitm("I am unsure of what to do now! Stopping!")
    end if
   
  end subroutine get_ns
   
  ! -----------------------------------------------------------------
  !  Take tides if input day is one of days that provide
  !  tidi fitted coefficients 
  ! -----------------------------------------------------------------
  subroutine calc_1tide_iday()
   
    use ModInputs, only : iDebugLevel

    implicit none
   
    character (len=iCharLenFile_) :: NameOfCoefFile1
    real, dimension(nHme)  :: amp_arr1,pha_arr1
       
    real, dimension(:),allocatable :: amp_fac1,pha_shift1
    real, dimension(nLonHme,nLatsHme,nAltHme) :: u1,v1,geopt1,temp1,dr_rho1
    integer :: i,j,k
    character(len=3) :: mm
       
    if (ReadFiles) then
       write(mm,'(I3.3)') iday
       NameOfCoefFile1 = trim(hmeDir) // 'tidi_coef/' // &
            'tidi_coef_2020' // mm // '.txt'
       if (iDebuglevel > 2) &
            write(*,*) 'Reading coef file for iday: ', NameOfCoefFile1
       call read_coef_file(NameOfCoefFile1)
       Storage_Tidi_amp(:,1) = amp_fac_arr
       Storage_Tidi_phase(:,1) = pha_shift_arr
    else
       ! Use the stored values:
       amp_fac_arr = Storage_Tidi_amp(:,1)
       pha_shift_arr = Storage_Tidi_phase(:,1)
    endif
    
    if (count(lpHme_slt) .eq. 0) then
       call stop_gitm('No HMEs found! Plese check selected calc_1tide_iday)')
    end if

    ! select HME component
    allocate(hmes_slt(count(lpHme_slt)))
    allocate(amps_slt(count(lpHme_slt)))
    allocate(phas_slt(count(lpHme_slt)))
    hmes_slt = pack(hme_arr, lpHme_slt)
    amps_slt = pack(amp_fac_arr, lpHme_slt)
    phas_slt = pack(pha_shift_arr, lpHme_slt)
   
    if (iDebuglevel > 2) &
         write(*,*) 'Now constructing tides for HMEs/iday: '
    call calc_1tide_1day(u1,v1,geopt1,temp1,dr_rho1) ! 1day
    deallocate(amps_slt)
    deallocate(phas_slt)
    deallocate(hmes_slt)
   
    u_3d = u1
    v_3d = v1
    geopt_3d = geopt1
    temp_3d  = temp1
    dr_rho_3d = dr_rho1
   
  end subroutine calc_1tide_iday
  
  ! -----------------------------------------------------------------
  !  Take tides of two days and intepolate them for the input day
  ! -----------------------------------------------------------------
  subroutine calc_1tide(day1,day2)

    use ModInputs, only : iDebugLevel
    
    implicit none
    integer,intent(in) :: day1,day2
   
    character (len=iCharLenFile_) :: NameOfCoefFile1,NameOfCoefFile2
    real, dimension(nLonHme,nLatsHme,nAltHme) :: u1,v1,geopt1,temp1,dr_rho1, &
         u2,v2,geopt2,temp2,dr_rho2
    integer :: i,j,k
    character(len=3) :: mm
       
    !----------------------------day1-------------------------------------
    if (ReadFiles) then
       write(mm,'(I3.3)') day1
       NameOfCoefFile1 = trim(hmeDir) // 'tidi_coef/' // &
            'tidi_coef_2020' // mm // '.txt'
       if (iDebuglevel > 2) &
            write(*,*) 'Reading coef file for day1: ', NameOfCoefFile1
       call read_coef_file(NameOfCoefFile1)
       
       if (count(lpHme_slt) .eq. 0) then
          call stop_gitm('No HMEs found! Plese check selected calc_1tide)')
       end if
       Storage_Tidi_amp(:,1) = amp_fac_arr
       Storage_Tidi_phase(:,1) = pha_shift_arr
    else
       ! Use the stored values:
       amp_fac_arr = Storage_Tidi_amp(:,1)
       pha_shift_arr = Storage_Tidi_phase(:,1)
    endif

    ! select HME component
    allocate(hmes_slt(count(lpHme_slt)))
    allocate(amps_slt(count(lpHme_slt)))
    allocate(phas_slt(count(lpHme_slt)))
    hmes_slt = pack(hme_arr,lpHme_slt)
    amps_slt = pack(amp_fac_arr,lpHme_slt)
    phas_slt = pack(pha_shift_arr,lpHme_slt)
    
    if (iDebuglevel > 2) &
         write(*,*) 'Now constructing tides for HMEs/day1: '
    call calc_1tide_1day(u1,v1,geopt1,temp1,dr_rho1) ! 1day
    deallocate(amps_slt)
    deallocate(phas_slt)
   
    !----------------------------day2-------------------------------------
    if (ReadFiles) then
       write(mm,'(I3.3)') day2
       NameOfCoefFile2 = trim(hmeDir) // 'tidi_coef/' // &
            'tidi_coef_2020' // mm // '.txt'
       if (iDebuglevel > 2) &
            write(*,*) 'Reading coef file for day2: ', NameOfCoefFile2
       call read_coef_file(NameOfCoefFile2)
       Storage_Tidi_amp(:,2) = amp_fac_arr
       Storage_Tidi_phase(:,2) = pha_shift_arr
    else
       ! Use the stored values:
       amp_fac_arr = Storage_Tidi_amp(:,2)
       pha_shift_arr = Storage_Tidi_phase(:,2)
    endif
       
    ! select HME component
    allocate(amps_slt(count(lpHme_slt)))
    allocate(phas_slt(count(lpHme_slt)))
    amps_slt = pack(amp_fac_arr,lpHme_slt)
    phas_slt = pack(pha_shift_arr,lpHme_slt)
   
    if (iDebuglevel > 2) &
         write(*,*) 'Now constructing tides for HMEs/day2: ',hmes_slt
    call calc_1tide_1day(u2,v2,geopt2,temp2,dr_rho2) ! 1day
       
    deallocate(amps_slt)
    deallocate(phas_slt)
    deallocate(hmes_slt)
   
    ! interp tides for iday   (nLonHme,nLatsHme,nAltHme)
   
    do i = 1, nLonHme
       do j = 1, nLatsHme
          do k = 1, nAltHme
             u_3d(i,j,k)     = (iday-day1)/(day2-day1) * &
                  (u2(i,j,k)-u1(i,j,k)) + u1(i,j,k)
             v_3d(i,j,k)     = (iday-day1)/(day2-day1) * &
                  (v2(i,j,k)-v1(i,j,k)) + v1(i,j,k)
             geopt_3d(i,j,k) = (iday-day1)/(day2-day1) * &
                  (geopt2(i,j,k)-geopt1(i,j,k)) + geopt1(i,j,k)
             temp_3d(i,j,k)  = (iday-day1)/(day2-day1) * &
                  (temp2(i,j,k)-temp1(i,j,k)) + temp1(i,j,k)
             dr_rho_3d(i,j,k) = (iday-day1)/(day2-day1) * &
                  (dr_rho2(i,j,k)-dr_rho1(i,j,k)) + dr_rho1(i,j,k)
          enddo
       enddo
    enddo       

  end subroutine calc_1tide
   
  ! -----------------------------------------------------------------
  !  Construct global maps of tides at three altitudes for each HME
  !  and sum them for all the selected HMEs
  ! -----------------------------------------------------------------
  subroutine calc_1tide_1day(u_all,v_all,geopt_all,temp_all,dr_rho_all)

    implicit none
       
    character (len=iCharLenHme_) :: hme_tmp
    real :: amp_fac, pha_shift
    integer :: n, s
       
    character (len=iCharLenFile_) :: NameOfHMEFile
    real, dimension(nLatsHme,nAltHme) :: u_a3,u_p3,v_a3,v_p3,geopt_a3,geopt_p3,&
         temp_a3,temp_p3,dr_rho_a3,dr_rho_p3
    real, dimension(nLonHme,nLatsHme,nAltHme) :: u_all,v_all,geopt_all,temp_all,&
         dr_rho_all,& 
         u,v,geopt,temp,dr_rho
    integer :: i,j,k,ipt,ilon
        
    ipt = 1

    u_all = 0.0
    v_all = 0.0
    geopt_all = 0.0
    temp_all = 0.0
    dr_rho_all = 0.0
    
    do i = 1,size(hmes_slt)
       hme_tmp = hmes_slt(i)
       call get_ns(hme_tmp(5:7),n,s)
         
       amp_fac = amps_slt(i)
       pha_shift = phas_slt(i)

       if (IsReadInHme) then
          u_a3 = StorageU_amp(i,:,:)
          u_p3 = StorageU_phase(i,:,:)

          v_a3 = StorageV_amp(i,:,:)
          v_p3 = StorageV_phase(i,:,:)

          temp_a3 = StorageT_amp(i,:,:)
          temp_p3 = StorageT_phase(i,:,:)

          geopt_a3 = StorageG_amp(i,:,:)
          geopt_p3 = StorageG_phase(i,:,:)

          dr_rho_a3 = StorageR_amp(i,:,:)
          dr_rho_p3 = StorageR_phase(i,:,:)
       else

          call hme_3alt(hme_tmp(5:7),u_a3,u_p3,v_a3,v_p3,geopt_a3,geopt_p3,&
               temp_a3,temp_p3,dr_rho_a3,dr_rho_p3)

          StorageU_amp(i,:,:) = u_a3
          StorageU_phase(i,:,:) = u_p3

          StorageV_amp(i,:,:) = v_a3
          StorageV_phase(i,:,:) = v_p3

          StorageT_amp(i,:,:) = temp_a3
          StorageT_phase(i,:,:) = temp_p3

          StorageG_amp(i,:,:) = geopt_a3
          StorageG_phase(i,:,:) = geopt_p3

          StorageR_amp(i,:,:) = dr_rho_a3
          StorageR_phase(i,:,:) = dr_rho_p3
          
       endif
          
       do j=1,nLatsHme
          do k = 1,nAltHme ! alt = 95.0, 97.5, 100 km
             do ilon = 1,nLonHme
   
                u(ilon,j,k) = &
                     amp_fac * u_a3(j,k) * cos((n * omega * ut + &
                     s * lon(ilon) - u_p3(j,k) * (n * omega) - &
                     pha_shift) * pi/180.0)
               
                v(ilon,j,k) = &
                     amp_fac * v_a3(j,k) * cos((n * omega * ut + &
                     s * lon(ilon) - v_p3(j,k) * (n * omega) - &
                     pha_shift) * pi/180.0)
   
                geopt(ilon,j,k) = &
                     amp_fac * geopt_a3(j,k) * &
                     cos((n * omega * ut + &
                     s * lon(ilon) - geopt_p3(j,k) * (n * omega)- &
                     pha_shift) &
                     * pi/180.0)
               
                temp(ilon,j,k) = &
                     amp_fac * temp_a3(j,k) * &
                     cos((n * omega * ut + &
                     s * lon(ilon) - temp_p3(j,k) * (n * omega) - &
                     pha_shift) * pi/180.0)
   
                dr_rho(ilon,j,k) = &
                     amp_fac * dr_rho_a3(j,k) * &
                     cos((n * omega * ut + &
                     s * lon(ilon) - dr_rho_p3(j,k) * (n * omega)- &
                     pha_shift) &
                     * pi/180.0)
             enddo
          enddo
       enddo
       if (count(isnan(geopt)) >0) then 
          write(*,*) isnan(geopt)
          write(*,*) 'i,ilon,ilat,ialt,hme: ',i,ilon,j,k,hme_tmp
          call stop_gitm('nan was found in geopt, stop')
       endif
   
       if (ipt .eq. 1) then
          u_all = u
          v_all = v
          geopt_all = geopt
          temp_all = temp
          dr_rho_all = dr_rho
       else 
          u_all = u_all + u
          v_all = v_all + v
          geopt_all = geopt_all + geopt
          temp_all = temp_all + temp
          dr_rho_all = dr_rho_all + dr_rho
       endif
   
       ipt = ipt+1
       
    enddo

    IsReadInHme = .true.
   
  end subroutine calc_1tide_1day

  ! -----------------------------------------------------------------
  !  Get HMEs at three altitudes 
  ! -----------------------------------------------------------------
  subroutine hme_3alt(thme,u_a3,u_p3,v_a3,v_p3,geopt_a3,geopt_p3,&
                       temp_a3,temp_p3,dr_rho_a3,dr_rho_p3)

    implicit none
    character (len=iCharLenHme_),intent(in) :: thme
    character (len=iCharLenFile_) :: NameOfHMEFile
    real, dimension(nLatsHme)  :: u_a,u_p,v_a,v_p,geopt_a,geopt_p,&
         temp_a,temp_p,dr_rho_a,dr_rho_p
    real, dimension(nLatsHme,nAltHme),intent(out) :: u_a3,u_p3,v_a3,v_p3, &
         geopt_a3,geopt_p3,&
         temp_a3,temp_p3,&
         dr_rho_a3,dr_rho_p3
    character (len=5) :: tmp
    integer :: i
       
    tmp = trim(thme) 
    do i = 1,nAltHme
   
       NameOfHMEFile = (trim(hmeDir) // 'HME_3alt/'// 'HME_'// trim(tmp)// &
            trim('_') // trim(altn(i)) // trim('00m-F75.txt'))
       call read_hme_file(NameOfHMEFile,u_a,u_p,v_a,v_p,&
            geopt_a,geopt_p,temp_a,temp_p,&
            dr_rho_a,dr_rho_p)
   
       u_a3(:,i) = u_a
       u_p3(:,i) = u_p
       v_a3(:,i) = v_a
       v_p3(:,i) = v_p
       geopt_a3(:,i) = geopt_a
       geopt_p3(:,i) = geopt_p
       temp_a3(:,i) = temp_a
       temp_p3(:,i) = temp_p
       dr_rho_a3(:,i) = dr_rho_a
       dr_rho_p3(:,i) = dr_rho_p
         
    enddo
          
  end subroutine hme_3alt
   
end module ModHmeModel 
