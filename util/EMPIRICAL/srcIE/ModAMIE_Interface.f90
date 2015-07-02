!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModAMIE_Interface

  character (len=100) :: AMIE_FileName
  integer :: AMIE_nLats, AMIE_nMlts, AMIE_nTimes

  ! For a single file
  real*4, allocatable,dimension(:)     :: AMIE_Lats, AMIE_MLTs
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_Potential,AMIE_EFlux,AMIE_AveE
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_Value
  real*8, allocatable,dimension(:,:)     :: AMIE_Time

  integer, parameter :: AMIE_Closest_     = 1
  integer, parameter :: AMIE_After_       = 2
  integer, parameter :: AMIE_Interpolate_ = 3

  integer :: AMIE_iDebugLevel = 0

  integer :: AMIE_South_ = 1
  integer :: AMIE_North_ = 2

  integer, parameter :: potential_ = 1
  integer, parameter :: eflux_     = 2
  integer, parameter :: avee_      = 3

end Module ModAMIE_Interface
