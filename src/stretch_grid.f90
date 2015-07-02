!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

! This is the function used to distribute the range of grid lines
! over the range of latitude values. The input value x will be
! in the range [0,1] and y, the output value, should also fall 
! within this range. Setting this function to y=x will 
! distribute the grid lines evenly among the latitude range, less
! linear distributions can be achived with other functions. 
! 0 is the bottom of both the latitude and grid value ranges, similarly
! 1 is the top of both. 

subroutine stretch_grid(x,y)

  use ModInputs, only : ConcentrationLatitude, &
       StretchingPercentage,StretchingFactor, &
       NewStretchedGrid, StretchWidth

  use ModConstants, only : pi

  implicit none

  real, intent(in):: x
  real, intent(out) :: y

  real :: xnorm, xnorm2, xdist, lat

  real :: latoffnorm, StretchingPercentage2, Percent
  real :: OneMSP, OneMSP2, SF, diff, FrontFactor, InvWidth

  latoffnorm = ConcentrationLatitude/90.0

  if (StretchingPercentage > 1.0) &
       StretchingPercentage = StretchingPercentage/100.0

  y = x

  xnorm = (x-0.5)*2.0

  if (xnorm > 1.0) xnorm = 2.0 - xnorm
  if (xnorm < -1.0) xnorm = -2.0 - xnorm

  lat = xnorm * 90.0

  if (NewStretchedGrid) then

     if (StretchWidth > 20.0) StretchWidth = 20.0
     if (StretchWidth < 1.0)  StretchWidth = 1.0
     Percent = 1.0 - (20.0-StretchWidth)/20.0
     InvWidth = 200.0/StretchWidth

     FrontFactor = StretchingPercentage * (0.05 + 0.3*Percent)

     diff = FrontFactor*exp(-((-1.0+latoffnorm)**2.0)*InvWidth) - &
            FrontFactor*exp(-((-1.0-latoffnorm)**2.0)*InvWidth)

     y = xnorm+FrontFactor*exp(-((xnorm+latoffnorm)**2.0)*InvWidth) - &
               FrontFactor*exp(-((xnorm-latoffnorm)**2.0)*InvWidth) + &
               diff * xnorm

     y = y/2.0 + 0.5

  else

     StretchingPercentage2 = StretchingPercentage*1.2

     SF = 1.0/StretchingFactor

     OneMSP = 1.0 - StretchingPercentage
     OneMSP2 = 1.0 - StretchingPercentage2

     if (lat > 0.0) then

        if (lat <= ConcentrationLatitude) then

           xnorm2 = lat/ConcentrationLatitude * pi/2.0 

           xdist = sin(xnorm2)**SF * latoffnorm
           xdist = StretchingPercentage*xdist + OneMSP*xnorm

        else

           xnorm2 = (1.0 - (90.0-lat)/(90.0-ConcentrationLatitude)) * pi/2.0
           xdist = (1.0 - cos(xnorm2)**SF)*(1.0-latoffnorm) + latoffnorm
           xdist = StretchingPercentage2*xdist + OneMSP2*xnorm

        endif

     else

        lat = abs(lat)

        if (lat < ConcentrationLatitude) then

           xnorm2 = lat/ConcentrationLatitude * pi/2.0 
           
           xdist = sin(xnorm2)**SF * latoffnorm
           xdist = StretchingPercentage*xdist - OneMSP*xnorm

        else

           xnorm2 = (1.0 - (90.0-lat)/(90.0-ConcentrationLatitude)) * pi/2.0
           xdist = (1.0 - cos(xnorm2)**SF)*(1.0-latoffnorm) + latoffnorm
           xdist = StretchingPercentage2*xdist - OneMSP2*xnorm

        endif

        xdist = -xdist

     endif

     y = xdist/2.0 + 0.5

  endif

  if (x > 1.0) y = 2.0 - y
  if (x < 0.0) y = -y

end subroutine stretch_grid

subroutine stretch_grid_old(x,y)

  use ModInputs, only : ConcentrationLatitude, &
       StretchingPercentage,StretchingFactor

  implicit none

  real, intent(in):: x
  real, intent(out) :: y

  real :: xnorm, xnorm2, xdist, lat
  real, parameter :: pi = 3.141592

  real :: latoffnorm, StretchingPercentage2
  real :: OneMSP, OneMSP2, SF

  latoffnorm = ConcentrationLatitude/90.0

  if (StretchingPercentage > 1.0) StretchingPercentage = StretchingPercentage/100.0

  StretchingPercentage2 = StretchingPercentage*1.2

  SF = 1.0/StretchingFactor

  OneMSP = 1.0 - StretchingPercentage
  OneMSP2 = 1.0 - StretchingPercentage2

  y = x

  xnorm = (x-0.5)*2.0

  if (xnorm > 1.0) xnorm = 2.0 - xnorm
  if (xnorm < -1.0) xnorm = -2.0 - xnorm

  lat = xnorm * 90.0

  if (lat > 0.0) then

    if (lat <= ConcentrationLatitude) then

      xnorm2 = lat/ConcentrationLatitude * pi/2.0 

      xdist = sin(xnorm2)**SF * latoffnorm
      xdist = StretchingPercentage*xdist + OneMSP*xnorm

    else

      xnorm2 = (1.0 - (90.0-lat)/(90.0-ConcentrationLatitude)) * pi/2.0
      xdist = (1.0 - cos(xnorm2)**SF)*(1.0-latoffnorm) + latoffnorm
      xdist = StretchingPercentage2*xdist + OneMSP2*xnorm

    endif

  else

    lat = abs(lat)

    if (lat < ConcentrationLatitude) then

      xnorm2 = lat/ConcentrationLatitude * pi/2.0 

      xdist = sin(xnorm2)**SF * latoffnorm
      xdist = StretchingPercentage*xdist - OneMSP*xnorm

    else

      xnorm2 = (1.0 - (90.0-lat)/(90.0-ConcentrationLatitude)) * pi/2.0
      xdist = (1.0 - cos(xnorm2)**SF)*(1.0-latoffnorm) + latoffnorm
      xdist = StretchingPercentage2*xdist - OneMSP2*xnorm

    endif

     xdist = -xdist

  endif

  y = xdist/2.0 + 0.5

  if (x > 1.0) y = 2.0 - y
  if (x < 0.0) y = -y

end subroutine stretch_grid_old


