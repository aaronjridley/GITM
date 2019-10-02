
module ModReadSami3d

  use ModInputs, only: iCharLen_

  implicit none

  integer :: iSami_Op_ = 1
  integer :: iSami_O2p_ = 2
  integer :: iSami_NOp_ = 3
  integer :: iSami_Hp_ = 4
  integer :: iSami_Hep_ = 5
  integer :: iSami_e_ = 6
  integer :: iSami_Te_ = 7
  integer :: iSami_Ti_ = 8
  integer :: iSami_ViNV_ = 9
  integer :: iSami_ViE_ = 10
  integer :: iSami_ViB_ = 11
  
  integer, parameter :: nSamiCharLength = iCharLen_
  integer, parameter :: nSamiVarCharLength = 40

  integer, parameter :: nMaxVars = 200
  integer, parameter :: nMaxTimes = 500
  integer :: nSamiFiles = 0
  
  integer, parameter :: Real8_ = selected_real_kind(12,100)
  character (len=nSamiCharLength) :: SamiDir

  real :: dtor = 3.141592653589793 / 180.0

  integer :: iSamiUnit = 39

  character (len=nSamiCharLength) :: SamiFileList(nMaxVars)
  character (len=nSamiVarCharLength) :: SamiVariables(nMaxVars)
  real(Real8_) :: SamiTimes(nMaxTimes)

  integer :: nLonsSami, nLatsSami, nAltsSami, nTimesSami, nVarsSami
  
  real, allocatable :: SamiLons(:), SamiLats(:), SamiAlts(:)

  real(Real8_) :: SamiStartTime

  ! This is for linear interpolation in time:
  integer :: iSamiIndex = 0, iSamiIndexOld = 0
  real :: rSamiFactor = 0.0
  
  real*4, allocatable :: SamiDataDummy(:,:,:)
  real, allocatable :: SamiDataOneTime(:,:,:,:)
  real, allocatable :: SamiDataAtTime(:,:,:,:)
  real, allocatable :: SamiDataTwoTimes(:,:,:,:,:)

  real, allocatable :: SamiInLons(:), SamiInLats(:), SamiInAlts(:)
  integer, allocatable :: SamiLonsIndex(:), SamiLatsIndex(:), SamiAltsIndex(:)
  real, allocatable :: SamiLonsFactor(:), SamiLatsFactor(:), SamiAltsFactor(:)
  real, allocatable :: SamiOutData(:,:)
  integer :: nPointsToGet = 0
  integer :: nPointsToGetSami = 0

  logical :: CorotationAdded = .true.
  
contains

  !---------------------------------------------------------------------------
  ! Set directory for SAMI files
  !---------------------------------------------------------------------------

  subroutine SamiSetDir(InDir)

    implicit none

    character (len=nSamiCharLength), intent(in) :: InDir
    SamiDir = InDir
    
  end subroutine SamiSetDir

  !---------------------------------------------------------------------------
  ! Get number of variables for SAMI
  !---------------------------------------------------------------------------
  
  subroutine SamiGetnVars(nVars)

    implicit none

    integer, intent(out) :: nVars
    nVars = nVarsSami-4
    
  end subroutine SamiGetnVars

  !---------------------------------------------------------------------------
  ! Get variables from SAMI
  !---------------------------------------------------------------------------
  
  subroutine SamiGetVars(OutVars)

    implicit none

    character (len=nSamiVarCharLength), dimension(nVarsSami), intent(out) :: OutVars

    OutVars(1:nVarsSami-4) = SamiVariables(5:nVarsSami)
    
  end subroutine SamiGetVars

  !---------------------------------------------------------------------------
  ! Read input file that describes sami data
  !---------------------------------------------------------------------------
  
  subroutine SamiReadInputFile(inFile)

    implicit none

    character (len=nSamiCharLength), intent(in) :: inFile
    character (len=nSamiCharLength) :: line
    integer :: iError

    integer, dimension(7) :: iTime
    integer :: year, month, day, iVar

    open(iSamiUnit, file = infile, status = 'old')

    iError = 0

    do while (iError == 0)

       read(iSamiUnit,'(a)',iostat=iError) line

       select case (line)
       
       case ("#STARTDATE")
          read(iSamiUnit,*) year
          read(iSamiUnit,*) month
          read(iSamiUnit,*) day
          itime = 0
          itime(1) = year
          itime(2) = month
          itime(3) = day
          call time_int_to_real(iTime, SamiStartTime)
          
       case ("#COROTATIONPOTENTIALADDED")
          read(iSamiUnit,'(L)',iostat=iError) CorotationAdded

       case ("#DIR")
           read(iSamiUnit,'(a)',iostat=iError) SamiDir

       case ("#VARS")
          iVar = 1
          do while (len(trim(Line))>1)
             read(iSamiUnit,'(a)',iostat=iError) line
             SamiVariables(iVar) = trim(line)
             iVar = iVar+1
              write(*,*) 'Dir:',trim(SamiDir)
          enddo
          nVarsSami = iVar-2
          
       case ("#VARFILES")
          do iVar = 1, nVarsSami
             read(iSamiUnit,'(a)',iostat=iError) line
             SamiFileList(iVar) = line
          enddo
       endselect

    enddo
    
    close(iSamiUnit)

    do iVar = 1, nVarsSami
       if (trim(SamiVariables(iVar)) == 'time') then
          write(*,*) 'Reading Time File : ', trim(SamiFileList(iVar))
          call ReadSamiTimeFile(SamiFileList(iVar))
       endif
       if (trim(SamiVariables(iVar)) == 'lon') then
          write(*,*) 'Reading Lon File : ', trim(SamiFileList(iVar))
          call ReadSamiLonFile(SamiFileList(iVar))
       endif
       if (trim(SamiVariables(iVar)) == 'lat') then
          write(*,*) 'Reading Lat File : ', trim(SamiFileList(iVar))
          call ReadSamiLatFile(SamiFileList(iVar))
       endif
       if (trim(SamiVariables(iVar)) == 'alt') then
          write(*,*) 'Reading Alt File : ', trim(SamiFileList(iVar))
          call ReadSamiAltFile(SamiFileList(iVar))
       endif
    enddo

    allocate(SamiDataDummy(nLatsSami-2, nAltsSami, nLonsSami-2))
    allocate(SamiDataOneTime(nVarsSami-4, nLonsSami, nLatsSami, nAltsSami))
    allocate(SamiDataAtTime(nVarsSami-4, nLonsSami, nLatsSami, nAltsSami))
    allocate(SamiDataTwoTimes(2, nVarsSami-4, nLonsSami, nLatsSami, nAltsSami))
    
  end subroutine SamiReadInputFile

  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine ReadSamiTimeFile(inFile)

    implicit none

    character (len=nSamiCharLength), intent(in) :: inFile
    integer :: iError
    integer :: iT, iHour, iMinute, iSecond, n
    real :: tHours, uts

    !open(iSamiUnit, file = trim(SamiDir)//'/'//infile, status = 'old')
      open(iSamiUnit, file =infile, status = 'old')
    iError = 0
    iT = 1
    
    do while (iError == 0)

       read(iSamiUnit,*,iostat=iError) n, iHour, iMinute, iSecond, tHours
       uts = iHour * 3600.0 + iMinute*60.0 + iSecond
       SamiTimes(iT) = tHours*3600.0 + SamiStartTime
       iT = iT + 1
       
    enddo

    nTimesSami = iT-1
    
    close(iSamiUnit)
       
  end subroutine ReadSamiTimeFile

  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine ReadSamiLonFile(inFile)

    implicit none

    character (len=nSamiCharLength), intent(in) :: inFile
    integer :: iError, iLon
    real :: lon

    ! Stupid way of doing it:
      write(*,*) trim(SamiDir)
    open(iSamiUnit, file = trim(SamiDir)//'/'//infile, status = 'old')
    
    iError = 0
    iLon = 1
    do while (iError == 0)
       read(iSamiUnit,*,iostat=iError) lon
       iLon = iLon + 1
    enddo
    ! Need to add a 0 and 360, so add 2 to nLons
    nLonsSami = iLon - 2 + 2
    close(iSamiUnit)

    allocate(SamiLons(nLonsSami))
        
    SamiLons(1) = 0.0
    open(iSamiUnit, file = trim(SamiDir)//'/'//infile, status = 'old')
    iError = 0
    do iLon=2,nLonsSami-1
       read(iSamiUnit,*,iostat=iError) lon
       SamiLons(iLon) = lon
    enddo
    close(iSamiUnit)
    SamiLons(nLonsSami) = 360.0

  end subroutine ReadSamiLonFile


  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine ReadSamiLatFile(inFile)

    implicit none

    character (len=nSamiCharLength), intent(in) :: inFile
    integer :: iError, iLat
    real :: lat

    ! Stupid way of doing it:

    open(iSamiUnit, file = trim(SamiDir)//'/'//infile, status = 'old')
    iError = 0
    iLat = 1
    do while (iError == 0)
       read(iSamiUnit,*,iostat=iError) lat
       iLat = iLat + 1
    enddo
    ! Need to add a -90 and 90, so add 2 to nLats
    nLatsSami = iLat - 2 + 2
    close(iSamiUnit)

    allocate(SamiLats(nLatsSami))
        
    SamiLats(1) = -90.0
    open(iSamiUnit, file = trim(SamiDir)//'/'//infile, status = 'old')
    iError = 0
    do iLat=2,nLatsSami-1
       read(iSamiUnit,*,iostat=iError) lat
       SamiLats(iLat) = lat
    enddo
    close(iSamiUnit)
    SamiLats(nLatsSami) = 90.0

  end subroutine ReadSamiLatFile

  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine ReadSamiAltFile(inFile)

    implicit none

    character (len=nSamiCharLength), intent(in) :: inFile
    integer :: iError, iAlt
    real :: alt

    ! Stupid way of doing it:

    open(iSamiUnit, file = trim(SamiDir)//'/'//infile, status = 'old')
    iError = 0
    iAlt = 1
    do while (iError == 0)
       read(iSamiUnit,*,iostat=iError) alt
       iAlt = iAlt + 1
    enddo
    nAltsSami = iAlt-2
    close(iSamiUnit)

    allocate(SamiAlts(nAltsSami))
        
    open(iSamiUnit, file = trim(SamiDir)//'/'//infile, status = 'old')
    iError = 0
    do iAlt=1,nAltsSami
       read(iSamiUnit,*,iostat=iError) alt
       SamiAlts(iAlt) = alt
    enddo
    close(iSamiUnit)

  end subroutine ReadSamiAltFile

  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine ReadSami(inFile,iTime)

    implicit none

    character (len=nSamiCharLength), intent(in) :: inFile
    integer, intent(in) :: iTime
    integer :: iError, iT

    open(iSamiUnit, file = trim(SamiDir)//'/'//trim(infile), status = 'old', &
         form='unformatted')

    iError = 0
    do iT=1,iTime
       read(iSamiUnit,iostat=iError) SamiDataDummy
    enddo
    if (iError > 0) then
       write(*,*) 'Error Reading Sami File : ',iError
       write(*,*) 'File : ',trim(SamiDir)//'/'//trim(infile)
    endif
    close(iSamiUnit)

  end subroutine ReadSami

  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine FileSamiData(iVar)

    implicit none

    integer, intent(in) :: iVar
    integer :: iLon, iLat, iAlt

    do iLon = 2, nLonsSami-1
       do iLat = 2, nLatsSami-1
          do iAlt = 1, nAltsSami
             SamiDataOneTime(iVar, iLon, iLat, iAlt) = &
                  SamiDataDummy(iLat-1, iAlt, iLon-1)
          enddo
       enddo
    enddo

    ! South Pole - average all longitudes
    iLat = 1
    do iAlt = 1, nAltsSami
       SamiDataOneTime(iVar, :, iLat, iAlt) = &
            sum(SamiDataOneTime(iVar, 2:nLonsSami-1, iLat+1, iAlt)) / &
            (nLonsSami-2)
    enddo

    ! North Pole - average all longitudes
    iLat = nLatsSami
    do iAlt = 1, nAltsSami
       SamiDataOneTime(iVar, :, iLat, iAlt) = &
            sum(SamiDataOneTime(iVar, 2:nLonsSami-1, iLat-1, iAlt)) / &
            (nLonsSami-2)
    enddo

    do iLat = 2, nLatsSami-1
       do iAlt = 1, nAltsSami
          SamiDataOneTime(iVar, 1, iLat, iAlt) =  &
               (SamiDataOneTime(iVar, 2, iLat, iAlt) + &
               SamiDataOneTime(iVar, nLonsSami-1, iLat, iAlt))/2 
          SamiDataOneTime(iVar, nLonsSami, iLat, iAlt) =  &
               SamiDataOneTime(iVar, 1, iLat, iAlt)
       enddo
    enddo

  end subroutine FileSamiData

  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine ReadSamiOneTime(iTime)

    implicit none

    integer, intent(in) :: iTime
    integer :: iError, iVar

    do iVar = 5, nVarsSami
       call ReadSami(SamiFileList(iVar),iTime)
       call FileSamiData(iVar-4)
    enddo
    
  end subroutine ReadSamiOneTime

  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
    
  subroutine SamiUpdateTime(InputTime, iError)

    implicit none
    
    real(Real8_), intent(in) :: InputTime
    integer, intent(out) :: iError
    
    iError = 0

    if (SamiTimes(1) > InputTime) then
       iError = 1
       write(*,*) "Requested time is before first file!"
       return
    endif
    if (SamiTimes(nTimesSami) < InputTime) then
       iError = 1
       write(*,*) "Requested time is after last file!"
       return
    endif

    iSamiIndex = 1
    do while (SamiTimes(iSamiIndex) <= InputTime .and. iSamiIndex < nTimesSami) 
       iSamiIndex = iSamiIndex + 1
    enddo

    rSamiFactor = (SamiTimes(iSamiIndex) - InputTime) / &
         (SamiTimes(iSamiIndex) - SamiTimes(iSamiIndex-1))

    if (iSamiIndex /= iSamiIndexOld) then

       write(*,*) 'Updating SAMI Times!'

       write(*,*) 'First Files'
       call ReadSamiOneTime(iSamiIndex-1)
       SamiDataTwoTimes(1,:,:,:,:) = SamiDataOneTime
       write(*,*) 'Second Files'
       call ReadSamiOneTime(iSamiIndex)
       SamiDataTwoTimes(2,:,:,:,:) = SamiDataOneTime

       iSamiIndexOld = iSamiIndex
       
    endif
       
    SamiDataAtTime = &
         rSamiFactor * SamiDataTwoTimes(1,:,:,:,:) + &
         (1.0 - rSamiFactor) * SamiDataTwoTimes(2,:,:,:,:)

  end subroutine SamiUpdateTime

  !---------------------------------------------------------------------------
  ! Get number of variables for SAMI
  !---------------------------------------------------------------------------
  
  subroutine SamiGetData(OutData)

    implicit none

    real, dimension(nPointsToGetSami,nVarsSami-4), intent(out) :: OutData

    integer :: iPoint, i, j, k
    real :: x, y, z

    do iPoint = 1, nPointsToGetSami

       i = SamiLonsIndex(iPoint)
       j = SamiLatsIndex(iPoint)
       k = SamiAltsIndex(iPoint)

       if (i > 0 .and. j > 0 .and. k > 0) then

          x = SamiLonsFactor(iPoint)
          y = SamiLatsFactor(iPoint)
          z = SamiAltsFactor(iPoint)
       
          OutData(iPoint,:) = &
               (1-x)*(1-y)*(1-z) * SamiDataAtTime(:,i  ,j  ,k) + &
               (  x)*(1-y)*(1-z) * SamiDataAtTime(:,i-1,j  ,k) + &
               (1-x)*(  y)*(1-z) * SamiDataAtTime(:,i  ,j-1,k) + &
               (  x)*(  y)*(1-z) * SamiDataAtTime(:,i-1,j-1,k) + &
               (1-x)*(1-y)*(  z) * SamiDataAtTime(:,i  ,j  ,k-1) + &
               (  x)*(1-y)*(  z) * SamiDataAtTime(:,i-1,j  ,k-1) + &
               (1-x)*(  y)*(  z) * SamiDataAtTime(:,i  ,j-1,k-1) + &
               (  x)*(  y)*(  z) * SamiDataAtTime(:,i-1,j-1,k-1)

       else
          OutData(iPoint,:) = -1.0e32
       endif
          
    enddo
    
  end subroutine SamiGetData


  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine SamiSetGrid(InLons, InLats, InAlts)

    implicit none

    real, dimension(nPointsToGetSami), intent(in) :: InLons
    real, dimension(nPointsToGetSami), intent(in) :: InLats
    real, dimension(nPointsToGetSami), intent(in) :: InAlts

    integer :: iPoint, i
    
    SamiInLons = InLons
    SamiInLats = InLats
    SamiInAlts = inAlts

    ! Check to see if points are valid:
    do iPoint = 1, nPointsToGetSami
       if (SamiInLats(iPoint) < -90.0) then
          SamiInLats(iPoint) = -180.0 - SamiInLats(iPoint)
          SamiInLons(iPoint) = SamiInLons(iPoint) + 180.0
       endif
       if (SamiInLats(iPoint) > 90.0) then
          SamiInLats(iPoint) = 180.0 - SamiInLats(iPoint)
          SamiInLons(iPoint) = SamiInLons(iPoint) + 180.0
       endif
       SamiInLons(iPoint) = mod(SamiInLons(iPoint) + 360, 360.0)
    enddo

    do iPoint = 1, nPointsToGetSami
       
       ! Lons First
       if (SamiInLons(iPoint) < SamiLons(1)) then
          SamiLonsIndex(iPoint) = -1
       else
          if (SamiInLons(iPoint) > SamiLons(nLonsSami)) then
             SamiLonsIndex(iPoint) = -1
          else
             if (SamiInLons(iPoint) == SamiLons(nLonsSami)) then
                i = nLonsSami
             else
                i = 2
                do while (SamiLons(i) <= SamiInLons(iPoint))
                   i = i + 1
                enddo
             endif
             SamiLonsIndex(iPoint) = i
             SamiLonsFactor(iPoint) = &
                  (SamiLons(i) - SamiInLons(iPoint)) / &
                  (SamiLons(i) - SamiLons(i-1))
          endif
       endif

       ! Lats
       if (SamiInLats(iPoint) < SamiLats(1)) then
          i = -1
       else
          if (SamiInLats(iPoint) > SamiLats(nLatsSami)) then
             i = -1
          else
             if (SamiInLats(iPoint) == SamiLats(nLatsSami)) then
                i = nLatsSami
             else
                i = 2
                do while (SamiLats(i) <= SamiInLats(iPoint))
                   i = i + 1
                enddo
             endif
             SamiLatsIndex(iPoint) = i
             SamiLatsFactor(iPoint) = &
                  (SamiLats(i) - SamiInLats(iPoint)) / &
                  (SamiLats(i) - SamiLats(i-1))
          endif
       endif
       
       ! Alts
       if (InAlts(iPoint) < SamiAlts(1)) then
          i = -1
       else
          if (InAlts(iPoint) > SamiAlts(nAltsSami)) then
             ! This is a hack (extrapolation), which won't work some of the time
             i = nAltsSami
             SamiAltsIndex(iPoint) = i
             SamiAltsFactor(iPoint) = &
                  (SamiAlts(i) - InAlts(iPoint)) / &
                  (SamiAlts(i) - SamiAlts(i-1))
          else
             if (InAlts(iPoint) == SamiAlts(nAltsSami)) then
                i = nAltsSami
             else
                i = 2
                do while (SamiAlts(i) <= InAlts(iPoint))
                   i = i + 1
                enddo
             endif
             SamiAltsIndex(iPoint) = i
             SamiAltsFactor(iPoint) = &
                  (SamiAlts(i) - InAlts(iPoint)) / &
                  (SamiAlts(i) - SamiAlts(i-1))
          endif
       endif

       if (SamiAltsIndex(iPoint) < 1) then
          write(*,*) 'bad alt in SamiSetGrid!', &
               SamiInLons(iPoint), SamiInLats(iPoint), InAlts(iPoint)
          write(*,*) 'Sami Alts : ', SamiAlts(1), SamiAlts(nAltsSami)
          stop
       endif
       
    enddo

  end subroutine SamiSetGrid

  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  
  subroutine SamiSetnPointsToGet(nPointsToGet)

    implicit none

    integer, intent(in) :: nPointsToGet
    
    nPointsToGetSami = nPointsToGet

    allocate(SamiInLons(nPointsToGet))
    allocate(SamiInLats(nPointsToGet))
    allocate(SamiInAlts(nPointsToGet))

    allocate(SamiLonsIndex(nPointsToGet))
    allocate(SamiLatsIndex(nPointsToGet))
    allocate(SamiAltsIndex(nPointsToGet))

    allocate(SamiLonsFactor(nPointsToGet))
    allocate(SamiLatsFactor(nPointsToGet))
    allocate(SamiAltsFactor(nPointsToGet))

    allocate(SamiOutData(nPointsToGet, nVarsSami))

  end subroutine SamiSetnPointsToGet
  
  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------

  subroutine SamiShutDown

    deallocate(SamiDataDummy)
    deallocate(SamiDataOneTime)
    deallocate(SamiDataAtTime)
    deallocate(SamiDataTwoTimes)

    deallocate(SamiInLons)
    deallocate(SamiInLats)
    deallocate(SamiInAlts)
    deallocate(SamiLonsIndex)
    deallocate(SamiLatsIndex)
    deallocate(SamiAltsIndex)
    deallocate(SamiLonsFactor)
    deallocate(SamiLatsFactor)
    deallocate(SamiAltsFactor)
    deallocate(SamiOutData)

    deallocate(SamiAlts)
    deallocate(SamiLats)
    deallocate(SamiLons)

    
  end subroutine SamiShutDown
  
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

  
end module ModReadSami3d
