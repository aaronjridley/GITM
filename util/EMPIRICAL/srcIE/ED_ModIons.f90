!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModIons

  integer, parameter :: DFP=SELECTED_REAL_KIND(15,307)  
  integer, parameter :: NLng=72    ! number of geographic longitude
  integer, parameter :: NLat=28    ! number of geographic latitude
  integer ::            NAlt       ! number of altitude 
  real (kind=DFP), dimension(:), allocatable :: glng, glat, alt
  real (kind=DFP), dimension(:,:,:), allocatable :: IonizationRate
  real (kind=DFP), dimension(:,:,:), allocatable :: HeatingRate

end module ModIons
