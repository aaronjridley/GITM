! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

program read_gitm

  use ModReadGitm3d
  
  implicit none
  
  integer, dimension(7) :: iTime
  real(Real8_) :: time
  integer :: iErr, iStep, iPoint

  integer :: nVars = 0
  integer, parameter :: nPoints = 50
  real :: lats(nPoints), lons(nPoints), alts(nPoints)
  real, allocatable :: data(:,:)

  character (len=nGitmCharLength), allocatable :: vars(:)
  
  call GetGitmFileList(iErr)
  call GetGitmGeneralHeaderInfo(iErr)
  call GitmGetnVars(nVars)

  allocate(vars(nVars))
  allocate(data(nPoints, nVars))

  call GitmGetVars(vars)

  write(*,*) vars
  
  do iPoint = 1, nPoints
     lons(iPoint) = 360 - 5.0*(iPoint-1)
     lats(iPoint) = -90.0 + 2.0*(iPoint-1)
     alts(iPoint) = 100.0 + 10.0*(iPoint-1)
  enddo

  call GitmSetnPointsToGet(nPoints)
  
  iTime = (/ 2015,  03, 14, 00, 01, 00, 00 /)
  call time_int_to_real(itime, time)

  do iStep = 1,12
     call GitmSetGrid(lons,lats,alts)
     call GitmUpdateTime(time, iErr)
     if (iErr == 0) then 
        call GitmGetData(data)
        write(*,*) data(:,1)/dtor
     endif
     time = time + 60.0
     lons = lons-1.0
  enddo

  write(*,*) data(:,2)/dtor
  write(*,*) data(:,3)/1000
  
  deallocate(data,vars)
  
  call GitmShutDown

end program read_gitm
