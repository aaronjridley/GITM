
module ModReadGitm3d

  integer, parameter :: Real8_ = selected_real_kind(12,100)
  integer, parameter :: i3dall_ = 1
  integer, parameter :: i3dlst_ = 2
  integer, parameter :: nGitmCharLength = 100
  integer, parameter :: nGitmVarCharLength = 40

  integer :: iRho_ = 4
  integer :: iNeutralStart_ = 5
  integer :: iTn_ = 16, iVn_ = 17
  integer :: iIonStart_ = 26
  integer :: iTe_ = 36, iVi_ = 38
  ! This is the array that is loaded with data in set_horizontal_bcs
  real, allocatable :: GitmFileData(:,:)
    
  integer, parameter :: nMaxFiles = 300, nMaxVars = 200
  integer :: nGitmFiles = 0
  character (len=nGitmCharLength) :: GitmDir
  character (len=nGitmCharLength) :: GitmFileList(nMaxFiles)
  character (len=nGitmVarCharLength) :: GitmVariables(nMaxVars)
  real(Real8_) :: GitmFileTimes(nMaxFiles)

  real, allocatable :: GitmLons(:), GitmLats(:), GitmAlts(:)
  integer :: nLonsGitm, nLatsGitm, nAltsGitm, nVarsGitm
  real :: GitmVersion
  
  integer :: iGitmFileType = i3dall_ 

  integer :: iGitmUnit = 31

  real :: dtor = 3.141592653589793 / 180.0

  ! This is for linear interpolation in time:
  integer :: iGitmIndex = 0, iGitmIndexOld = 0
  real :: rGitmFactor = 0.0
  real, allocatable :: GitmDataTwoFiles(:,:,:,:,:), GitmDataAtTime(:,:,:,:)
  real, allocatable :: GitmDataDummy(:,:,:,:)

  real, allocatable :: GitmInLons(:), GitmInLats(:), GitmInAlts(:)
  integer, allocatable :: GitmLonsIndex(:), GitmLatsIndex(:), GitmAltsIndex(:)
  real, allocatable :: GitmLonsFactor(:), GitmLatsFactor(:), GitmAltsFactor(:)
  real, allocatable :: GitmOutData(:,:)
  integer :: nPointsToGet = 0
  integer :: nPointsToGetGitm = 0
  
contains

  !---------------------------------------------------------------------------
  ! Set directory for GITM files
  !---------------------------------------------------------------------------
  
  subroutine GitmSetDir(InDir)

    implicit none

    character (len=nGitmCharLength), intent(in) :: InDir
    GitmDir = InDir
    
  end subroutine GitmSetDir

  !---------------------------------------------------------------------------
  ! Get number of variables for GITM
  !---------------------------------------------------------------------------
  
  subroutine GitmGetnVars(nVars)

    implicit none

    integer, intent(out) :: nVars
    nVars = nVarsGitm
    
  end subroutine GitmGetnVars

  !---------------------------------------------------------------------------
  ! Get variables from GITM
  !---------------------------------------------------------------------------
  
  subroutine GitmGetVars(OutVars)

    implicit none

    character (len=nGitmVarCharLength), dimension(nVarsGitm), intent(out) :: OutVars

    OutVars(1:nVarsGitm) = GitmVariables(1:nVarsGitm)
    
  end subroutine GitmGetVars

  !---------------------------------------------------------------------------
  ! Get number of variables for GITM
  !---------------------------------------------------------------------------
  
  subroutine GitmGetData(OutData)

    implicit none

    real, dimension(nPointsToGetGitm,nVarsGitm), intent(out) :: OutData

    integer :: iPoint, i, j, k
    real :: x, y, z

    do iPoint = 1, nPointsToGetGitm

       i = GitmLonsIndex(iPoint)
       j = GitmLatsIndex(iPoint)
       k = GitmAltsIndex(iPoint)

       if (i > 0 .and. j > 0 .and. k > 0) then

          x = GitmLonsFactor(iPoint)
          y = GitmLatsFactor(iPoint)
          z = GitmAltsFactor(iPoint)
       
          OutData(iPoint,:) = &
               (1-x)*(1-y)*(1-z) * GitmDataAtTime(:,i  ,j  ,k) + &
               (  x)*(1-y)*(1-z) * GitmDataAtTime(:,i-1,j  ,k) + &
               (1-x)*(  y)*(1-z) * GitmDataAtTime(:,i  ,j-1,k) + &
               (  x)*(  y)*(1-z) * GitmDataAtTime(:,i-1,j-1,k) + &
               (1-x)*(1-y)*(  z) * GitmDataAtTime(:,i  ,j  ,k-1) + &
               (  x)*(1-y)*(  z) * GitmDataAtTime(:,i-1,j  ,k-1) + &
               (1-x)*(  y)*(  z) * GitmDataAtTime(:,i  ,j-1,k-1) + &
               (  x)*(  y)*(  z) * GitmDataAtTime(:,i-1,j-1,k-1)

       else
          OutData(iPoint,:) = -1.0e32
       endif
          
    enddo
    
  end subroutine GitmGetData


  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine GitmSetGrid(InLons, InLats, InAlts)

    implicit none

    real, dimension(nPointsToGetGitm), intent(in) :: InLons
    real, dimension(nPointsToGetGitm), intent(in) :: InLats
    real, dimension(nPointsToGetGitm), intent(in) :: InAlts

    integer :: iPoint, i
    
    GitmInLons = InLons
    GitmInLats = InLats
    GitmInAlts = inAlts

    do iPoint = 1, nPointsToGetGitm

       ! Lons First
       if (InLons(iPoint) < GitmLons(1)) then
          GitmLonsIndex(iPoint) = -1
       else
          if (InLons(iPoint) > GitmLons(nLonsGitm)) then
             GitmLonsIndex(iPoint) = -1
          else
             if (InLons(iPoint) == GitmLons(nLonsGitm)) then
                i = nLonsGitm
             else
                i = 2
                do while (GitmLons(i) <= InLons(iPoint))
                   i = i + 1
                enddo
             endif
             GitmLonsIndex(iPoint) = i
             GitmLonsFactor(iPoint) = &
                  (GitmLons(i) - InLons(iPoint)) / &
                  (GitmLons(i) - GitmLons(i-1))
          endif
       endif

       ! Lats
       if (InLats(iPoint) < GitmLats(1)) then
          i = -1
       else
          if (InLats(iPoint) > GitmLats(nLatsGitm)) then
             i = -1
          else
             if (InLats(iPoint) == GitmLats(nLatsGitm)) then
                i = nLatsGitm
             else
                i = 2
                do while (GitmLats(i) <= InLats(iPoint))
                   i = i + 1
                enddo
             endif
             GitmLatsIndex(iPoint) = i
             GitmLatsFactor(iPoint) = &
                  (GitmLats(i) - InLats(iPoint)) / &
                  (GitmLats(i) - GitmLats(i-1))
          endif
       endif
       
       ! Alts
       if (InAlts(iPoint) < GitmAlts(1)) then
          i = -1
       else
          if (InAlts(iPoint) > GitmAlts(nAltsGitm)) then
             i = -1
          else
             if (InAlts(iPoint) == GitmAlts(nAltsGitm)) then
                i = nAltsGitm
             else
                i = 2
                do while (GitmAlts(i) <= InAlts(iPoint))
                   i = i + 1
                enddo
             endif
             GitmAltsIndex(iPoint) = i
             GitmAltsFactor(iPoint) = &
                  (GitmAlts(i) - InAlts(iPoint)) / &
                  (GitmAlts(i) - GitmAlts(i-1))
          endif
       endif
    enddo

  end subroutine GitmSetGrid

  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine GitmSetnPointsToGet(nPointsToGet)

    integer, intent(in) :: nPointsToGet
    
    nPointsToGetGitm = nPointsToGet

    allocate(GitmInLons(nPointsToGet))
    allocate(GitmInLats(nPointsToGet))
    allocate(GitmInAlts(nPointsToGet))

    allocate(GitmLonsIndex(nPointsToGet))
    allocate(GitmLatsIndex(nPointsToGet))
    allocate(GitmAltsIndex(nPointsToGet))

    allocate(GitmLonsFactor(nPointsToGet))
    allocate(GitmLatsFactor(nPointsToGet))
    allocate(GitmAltsFactor(nPointsToGet))

    allocate(GitmOutData(nPointsToGet, nVarsGitm))

  end subroutine GitmSetnPointsToGet

  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
    
  subroutine GitmUpdateTime(InputTime, iError)

    real(Real8_), intent(in) :: InputTime
    integer, intent(out) :: iError
    
    iError = 0

    if (GitmFileTimes(1) > InputTime) then
       iError = 1
       write(*,*) "Requested time is before first file!"
       return
    endif
    if (GitmFileTimes(nGitmFiles) < InputTime) then
       iError = 1
       write(*,*) "Requested time is after last file!"
       return
    endif

    iGitmIndex = 1
    do while (GitmFileTimes(iGitmIndex) <= InputTime .and. iGitmIndex < nGitmFiles) 
       iGitmIndex = iGitmIndex + 1
    enddo

    rGitmFactor = (GitmFileTimes(iGitmIndex) - InputTime) / &
         (GitmFileTimes(iGitmIndex) - GitmFileTimes(iGitmIndex-1))

    if (iGitmIndex /= iGitmIndexOld) then

       write(*,*) 'Updating GITM files!'

       write(*,*) 'Reading First File : ', GitmFileList(iGitmIndex-1)
       call GitmReadFile(GitmFileList(iGitmIndex-1), iError)
       GitmDataTwoFiles(1,:,:,:,:) = GitmDataDummy

       write(*,*) 'Reading Second File : ', GitmFileList(iGitmIndex)
       call GitmReadFile(GitmFileList(iGitmIndex), iError)
       GitmDataTwoFiles(2,:,:,:,:) = GitmDataDummy

       iGitmIndexOld = iGitmIndex
       
    endif
       
    GitmDataAtTime = &
         rGitmFactor * GitmDataTwoFiles(1,:,:,:,:) + &
         (1.0 - rGitmFactor) * GitmDataTwoFiles(2,:,:,:,:)

  end subroutine GitmUpdateTime

  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine GetGitmFileList(iError)

    integer, intent(out) :: iError
    character (len=nGitmCharLength) :: Dummy
    real(Real8_) :: time

    iError = 0
    
    ! First try 3DALL files:
    write(*,*) 'ls -1 '//trim(GitmDir)//'3DALL*.bin > .list_of_gitm_files 2> .gitm_err'
    call system('ls -1 '//trim(GitmDir)//'3DALL*.bin > .list_of_gitm_files 2> .gitm_err')
    nGitmFiles = 0
    open(iGitmUnit, file = '.list_of_gitm_files', status = 'old')
    iError = 0
    do while (iError == 0)
       read(iGitmUnit, '(a)', iostat = iError) Dummy
       write(*,*) 'dummy : ',dummy
       if (iError == 0) nGitmFiles = nGitmFiles + 1
    enddo    
    close(iGitmUnit)

    if (nGitmFiles == 0) then
       ! First try 3DALL files:
       call system('ls -1 '//GitmDir//'3DLST*.bin > .list_of_gitm_files 2> .gitm_err')
       nGitmFiles = 0
       open(iGitmUnit, file = '.list_of_gitm_files', status = 'old')
       iError = 0
       do while (iError == 0)
          read(iGitmUnit, '(a)', iostat = iError) Dummy
          if (iError == 0) nGitmFiles = nGitmFiles + 1
       enddo
       close(iGitmUnit)
       iGitmFileType = i3dlst_ 
    endif

    write(*,*) "nGitmFiles : " ,nGitmFiles
    write(*,*) "Filetype : ", iGitmFileType

    open(iGitmUnit, file = '.list_of_gitm_files', status = 'old')
    do iFile = 1, nGitmFiles
       read(iGitmUnit, '(a)', iostat = iError) GitmFileList(iFile)
       call GetGitmTime(GitmFileList(iFile), time)
       GitmFileTimes(iFile) = time
    enddo
    close(iGitmUnit)

    if (nGitmFiles == 0) iError = 1
    
  end subroutine GetGitmFileList

  ! ----------------------------------------------------------------
  ! File for reading GITM 
  ! ----------------------------------------------------------------
  
  subroutine GetGitmTime(file, time)

    character (len=nGitmCharLength), intent(in) :: file
    real(Real8_), intent(out) :: time
    character (len=nGitmVarCharLength) :: Dummy

    real    :: Version
    integer :: nLons, nLats, nAlts, nVars, iVar
    integer :: iYear,iMonth,iDay,iHour,iMinute,iSecond,iMilli
    integer, dimension(7) :: iTime

    open(iInputUnit_, file=file, status="old",form="unformatted")

    read(iInputUnit_) Version
    read(iInputUnit_) nLons, nLats, nAlts
    read(iInputUnit_) nVars
    do iVar=1,nVars
       read(iInputUnit_) Dummy
    enddo

    read(iInputUnit_) iYear, iMonth, iDay, iHour, iMinute, iSecond, iMilli
    iTime = (/ iYear,  iMonth, iDay, iHour, iMinute, iSecond, iMilli /)
    call time_int_to_real(itime, time)
    
  end subroutine GetGitmTime

  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine GetGitmGeneralHeaderInfo(iError)

    integer, intent(out) :: iError
    
    integer :: iVar
    integer :: iYear,iMonth,iDay,iHour,iMinute,iSecond,iMilli
    integer, dimension(7) :: iTime
    real, allocatable :: data(:,:,:)

    iError = 0
    
    open(iInputUnit_, file=GitmFileList(1), &
         status="old", form="unformatted", iostat = iError)

    if (iError == 0) then 
    
       read(iInputUnit_) GitmVersion
       read(iInputUnit_) nLonsGitm, nLatsGitm, nAltsGitm
       read(iInputUnit_) nVarsGitm
       do iVar=1,nVarsGitm
          read(iInputUnit_) GitmVariables(iVar)
       enddo

       read(iInputUnit_) iYear, iMonth, iDay, iHour, iMinute, iSecond, iMilli

       ! Only need to read in the first 3 variables, since these are (should be!)
       ! Longitude, Latitude, Altitude (in that order):
       allocate(data(nLonsGitm,nLatsGitm,nAltsGitm))
       allocate(GitmLons(nLons))
       allocate(GitmLats(nLats))
       allocate(GitmAlts(nAlts))

       ! Longitude:
       read(iInputUnit_) data(:,:,:)
       GitmLons = data(:,1,1) / dtor ! convert to degrees

       ! Latitude:
       read(iInputUnit_) data(:,:,:)
       GitmLats = data(1,:,1) / dtor ! convert to degrees

       ! Altitude:
       read(iInputUnit_) data(:,:,:)
       GitmAlts = data(1,1,:) / 1000.0  ! Convert to km

       close(iInputUnit_)
       deallocate(data)

       allocate(GitmDataTwoFiles(2,nVarsGitm,nLonsGitm,nLatsGitm,nAltsGitm))
       allocate(GitmDataAtTime(nVarsGitm,nLonsGitm,nLatsGitm,nAltsGitm))
       allocate(GitmDataDummy(nVarsGitm,nLonsGitm,nLatsGitm,nAltsGitm))
       
    endif

  end subroutine GetGitmGeneralHeaderInfo

  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine GitmReadFile(InFile, iError)

    character (len=nGitmCharLength), intent(in) :: InFile
    integer, intent(out) :: iError
    
    character (len=nGitmVarCharLength) :: Dummy

    integer :: iVar
    integer :: iYear,iMonth,iDay,iHour,iMinute,iSecond,iMilli
    integer :: nLons, nLats, nAlts, nVars
    real :: Version
    integer, dimension(7) :: iTime

    iError = 0
    
    open(iInputUnit_, file=InFile, &
         status="old", form="unformatted", iostat = iError)

    if (iError == 0) then 
    
       read(iInputUnit_) Version
       read(iInputUnit_) nLons, nLats, nAlts
       read(iInputUnit_) nVars
       do iVar=1,nVars
          read(iInputUnit_) Dummy
       enddo

       read(iInputUnit_) iYear, iMonth, iDay, iHour, iMinute, iSecond, iMilli

       do iVar = 1, nVars
          read(iInputUnit_) GitmDataDummy(iVar,:,:,:)
       enddo

       close(iInputUnit_)
       
    endif

  end subroutine GitmReadFile

  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine GitmShutDown

    deallocate(GitmDataTwoFiles, GitmDataAtTime, GitmDataDummy)
    deallocate(GitmInLons)
    deallocate(GitmInLats)
    deallocate(GitmInAlts)
    deallocate(GitmOutData)


    
  end subroutine GitmShutDown

  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine time_int_to_real(itime, timereal)

    implicit none

    integer, dimension(1:12) :: dayofmon
    integer, dimension(1:7) :: itime
    double precision :: timereal
    integer :: nyear, nleap, nmonth, nday, nhour, nmin, nsec, i

    dayofmon(1) = 31
    dayofmon(2) = 28
    dayofmon(3) = 31
    dayofmon(4) = 30
    dayofmon(5) = 31
    dayofmon(6) = 30
    dayofmon(7) = 31
    dayofmon(8) = 31
    dayofmon(9) = 30
    dayofmon(10) = 31
    dayofmon(11) = 30
    dayofmon(12) = 31

    if (mod(itime(1),4).eq.0) dayofmon(2) = 29

    timereal = 0.0

    if (itime(1) > 1900) then
       nyear = itime(1) - 1965
    else 
       if (itime(1) > 65) then
          nyear = itime(1) - 65 
       else 
          nyear = itime(1) + 100 - 65
       endif
    endif
    nleap = nyear/4

    nmonth = itime(2) - 1

    nday = 0

    do i=1, nmonth
       nday = nday + dayofmon(i)
    enddo

    nday = nday + itime(3) - 1
    nhour = itime(4)
    nmin = itime(5)
    nsec = itime(6)

    timereal = (dble(nsec) * dble(1.0)) +                  &
         (dble(nmin) * dble(60.0)) +                       &
         (dble(nhour) * dble(60.0*60.0)) +                 &
         (dble(nday) * dble(24.0*60.0*60.0)) +             &
         (dble(nleap) * dble(24.0*60.0*60.0)) +            &
         (dble(nyear) * dble(365.0*24.0*60.0*60.0)) +      &
         itime(7)/1000.0

  end subroutine time_int_to_real
  
end module ModReadGitm3d
