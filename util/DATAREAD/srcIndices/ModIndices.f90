!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

Module ModIndices

  use ModKind
  use ModIoUnit, only : UnitTmp_

  implicit none

  integer, parameter :: MaxInputLines=10000

  integer, parameter :: GSM_ = 1
  integer, parameter :: GSE_ = 2

  integer, parameter :: f107_     =  1
  integer, parameter :: hpi_      =  2
  integer, parameter :: kp_       =  3
  integer, parameter :: au_       =  4
  integer, parameter :: al_       =  5
  integer, parameter :: ae_       =  6
  integer, parameter :: nAE_      =  7
  integer, parameter :: Dst_      =  8
  integer, parameter :: nDst_     =  9
  integer, parameter :: imf_bx_   = 10
  integer, parameter :: imf_by_   = 11
  integer, parameter :: imf_bz_   = 12
  integer, parameter :: sw_vx_    = 13
  integer, parameter :: sw_vy_    = 14
  integer, parameter :: sw_vz_    = 15
  integer, parameter :: sw_v_     = 16
  integer, parameter :: sw_n_     = 17
  integer, parameter :: sw_t_     = 18
  integer, parameter :: hpi_calc_ = 19
  integer, parameter :: jh_calc_  = 20
  integer, parameter :: hpi_norm_ = 21
  integer, parameter :: f107a_    = 22
  integer, parameter :: ap_       = 23
  integer, parameter :: onsetut_  = 24
  integer, parameter :: onsetmlat_= 25
  integer, parameter :: onsetmlt_ = 26

  integer, parameter :: MaxIndicesEntries = 60000
  integer, parameter :: nIndices  = onsetmlt_

  real, allocatable           :: Indices_TV(:,:)
  integer, dimension(nIndices)                           :: nIndices_V=0
  real (Real8_), allocatable :: IndexTimes_TV(:,:)
  real, dimension(nIndices)                              :: SavedIndices_V
  integer, dimension(nIndices)                           :: SavedErrors_V
  real (Real8_)                                          :: SavedTime

  real :: Propagation_Plane_XY, Propagation_Plane_XZ
  real :: Satellite_X_Pos, Satellite_Y_Pos, Satellite_Z_Pos

  character (len=100) :: NameOfIndexFile
  character (len=100) :: NameOfSecondIndexFile
  integer             :: LunIndices_ = UnitTmp_

  logical :: ReReadIMFFile  = .false.
  logical :: ReReadHPIFile  = .false.
  logical :: ReReadSMEFile  = .false.
  logical :: ReReadOnsetFile = .false.
  character (len=100) :: NameOfIMFFile = "none"
  character (len=100) :: NameOfHPIFile = "none"
  character (len=100) :: NameOfSMEFile = "none"
  character (len=100) :: NameOfOnsetFile = "none"

contains

  subroutine init_mod_indices

    if(allocated(Indices_TV)) return
    allocate( &
         Indices_TV(MaxIndicesEntries, nIndices), &
         IndexTimes_TV(MaxIndicesEntries, nIndices))

  end subroutine init_mod_indices

end Module ModIndices
