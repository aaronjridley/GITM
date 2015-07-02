!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModMHD_Interface

  character (len=100) :: MHD_FileName
  integer :: MHD_nLats, MHD_nMlts, MHD_nTimes
  integer :: MHD_Year, MHD_Month, MHD_Day

  ! For a single file
  real*4, allocatable,dimension(:)     :: MHD_Lats, MHD_MLTs
  real*4, allocatable,dimension(:,:,:,:) :: MHD_Potential,MHD_EFlux,MHD_AveE
  real*4, allocatable,dimension(:,:,:,:) :: MHD_Value
  real*8, allocatable,dimension(:,:)     :: MHD_Time

  integer, parameter :: MHD_Closest_     = 1
  integer, parameter :: MHD_After_       = 2
  integer, parameter :: MHD_Interpolate_ = 3

  integer :: MHD_iDebugLevel = 0

  integer :: MHD_South_ = 1
  integer :: MHD_North_ = 2

end Module ModMHD_Interface
