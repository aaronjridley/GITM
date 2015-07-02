!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

Module ModTides

  !============================================================================
  ! GSWM variables that are necessary for interpolation on the GITM grid. 
  ! T, u, and v are neutral temperature, zonal and meridional winds. 
  ! 
  ! EYigit: 16June09
  !============================================================================

  use ModSizeGITM, only: nLons, nLats, nBlocksMax

  implicit none

  real, allocatable :: u_gswm(:,:,:,:)
  real, allocatable :: v_gswm(:,:,:,:)
  real, allocatable :: T_gswm(:,:,:,:)
  real, allocatable :: lon_gswm(:)
  real, allocatable :: lat_gswm(:)

  integer :: nLatsGSWM, nLonsGSWM, nAltsGswm
  integer :: nLatsWaccm, nLonsWaccm, nAltsWaccm

  character(len=40), dimension(4) :: GSWM_name           ! EYigit: 18May09
  character(len=40), dimension(4) :: GSWM_file_name      ! EYigit: 18May09

  real, dimension(-1:nLons+2, -1:nLats+2, 2, nBlocksMax) :: &
       TidesNorth, TidesEast, TidesTemp, TidesOmega, &
       TidesVertical, TidesN, rAlt_waccm

  integer, dimension(-1:nLons+2, -1:nLats+2, 2, nBlocksMax) :: iAlt_waccm

  real    :: rLon_Waccm(-1:nLons+2, nBlocksMax) 
  integer :: iLon_Waccm(-1:nLons+2, nBlocksMax) 

  real    :: rLat_Waccm(-1:nLats+2, nBlocksMax) 
  integer :: iLat_Waccm(-1:nLats+2, nBlocksMax) 

  real :: dLonGswm, dLatGswm

  real, allocatable :: t_waccm(:,:,:)
  real, allocatable :: t0_waccm(:,:,:)
  real, allocatable :: ta1_waccm(:,:,:)
  real, allocatable :: ta2_waccm(:,:,:)
  real, allocatable :: tp1_waccm(:,:,:)
  real, allocatable :: tp2_waccm(:,:,:)

  real, allocatable :: u_waccm(:,:,:)
  real, allocatable :: u0_waccm(:,:,:)
  real, allocatable :: ua1_waccm(:,:,:)
  real, allocatable :: ua2_waccm(:,:,:)
  real, allocatable :: up1_waccm(:,:,:)
  real, allocatable :: up2_waccm(:,:,:)

  real, allocatable :: v_waccm(:,:,:)
  real, allocatable :: v0_waccm(:,:,:)
  real, allocatable :: va1_waccm(:,:,:)
  real, allocatable :: va2_waccm(:,:,:)
  real, allocatable :: vp1_waccm(:,:,:)
  real, allocatable :: vp2_waccm(:,:,:)

  real, allocatable :: w_waccm(:,:,:)
  real, allocatable :: omega_waccm(:,:,:)
  real, allocatable :: omega0_waccm(:,:,:)
  real, allocatable :: omegaa1_waccm(:,:,:)
  real, allocatable :: omegaa2_waccm(:,:,:)
  real, allocatable :: omegap1_waccm(:,:,:)
  real, allocatable :: omegap2_waccm(:,:,:)

  real, allocatable :: lon_waccm(:)
  real, allocatable :: lat_waccm(:)
  real, allocatable :: alt_waccm(:,:,:)
  real, allocatable :: pressure_waccm(:)
  real, allocatable :: Alt_waccm_GitmGrid(:,:,:,:)

!  Vars= ['Longitude','Latitude','Altitude', 'Pressure', $
!       'T','T_24_COS','T_24_SIN','T_12_COS','T_12_SIN', $
!       'U','U_24_COS','U_24_SIN','U_12_COS','U_12_SIN', $
!       'V','V_24_COS','V_24_SIN','V_12_COS','V_12_SIN', $
!       'OMEGA','OMEGA_24_COS','OMEGA_24_SIN','OMEGA_12_COS','OMEGA_12_SIN']

end Module ModTides


