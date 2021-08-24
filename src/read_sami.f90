! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

program read_sami

  use ModReadSami3D

  implicit none

  character (len = nSamiCharLength) :: InFile = 'sami_var.txt'

  integer, dimension(7) :: iTime
  real(Real8_) :: time

  integer :: iStep, iErr

  real, allocatable :: data(:,:)
  character (len=nSamiCharLength), allocatable :: vars(:)
  integer :: nVars = 0, iPoint

  integer, parameter :: nPoints = 50
  real :: lats(nPoints), lons(nPoints), alts(nPoints)
  
  call SamiReadInputFile(InFile)
  call SamiGetnVars(nVars)

  write(*,*) nVars
  
  allocate(vars(nVars))
  allocate(data(nPoints,nVars))

  call SamiGetVars(vars)

  write(*,*) vars
  
  do iPoint = 1, nPoints
     lons(iPoint) = 360 - 5.0*(iPoint-1)
     lats(iPoint) = -90.0 + 2.0*(iPoint-1)
     alts(iPoint) = 100.0 + 10.0*(iPoint-1)
  enddo

  call SamiSetnPointsToGet(nPoints)
  call SamiSetGrid(lons,lats,alts)
  
  iTime = (/ 2015,  03, 16, 00, 01, 00, 00 /)
  call time_int_to_real(itime, time)

  do iStep = 1,1
     write(*,*) iStep, time
     call SamiUpdateTime(time, iErr)
     if (iErr == 0) then 
        call SamiGetData(data)
        write(*,*) data(:,1)
     endif     
     time = time + 600.0
  enddo
  
  call SamiShutDown

  stop
  
end program read_sami
