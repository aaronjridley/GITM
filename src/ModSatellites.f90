!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!----------------------------------------------------------------------------
! $Id: ModSatellites.f90,v 1.11 2015/07/02 13:36:25 ridley Exp $
! Author: Aaron Ridley, UMichigan
!
! Modified:
!           AGB, Oct 2013 - Added seperate routine to allocate satellite data
!                           structures if data assimilation will be performed
!           Asad, Feb 2013 - Modified to allow satellite data to be read in.
!                            Increased maximum allowable size of input file.
!
! Comments: Variables used to track satellites through the GITM grid and
!           routines used to initialize these variables
!----------------------------------------------------------------------------

module ModSatellites

  use ModInputs, only : iCharLen_
  use ModIoUnit, only : UnitTmp_

  implicit none

  integer, parameter :: nMaxSats = 115
  integer, parameter :: nMaxSatInputLines = 100000 !Asad: increased for RCAC
  integer, parameter :: nMaxSatPos = 120

  integer :: nSats = 0

  integer :: iSatUnit = UnitTmp_

  character (len=iCharLen_) :: cSatFileName(nMaxSats)

  !!! Xing
  character (len=iCharLen_) :: SatOutputType(nMaxSats)

  Integer, parameter :: dblprec = selected_real_kind(14,200)

  real (kind=dblprec), allocatable :: SatTime(:,:)
  real, allocatable                :: SatPos(:,:,:,:)
  integer, allocatable             :: nSatPos(:,:)
  real                :: SatCurrentPos(nMaxSats, 3, nMaxSatPos)
  real                :: SatDtPlot(nMaxSats)
  integer             :: nSatLines(nMaxSats)
  integer             :: iSatCurrentIndex(nMaxSats)

  real                :: CurrentSatellitePosition(3)
  character (len=50)   :: CurrentSatelliteName 

  integer :: CurrSat = 0
  integer :: nRCMRSat
  integer, dimension(nMaxSats) :: RCMRSat

  real, allocatable :: SatDat(:,:)
  real, allocatable :: SatCurrentDat(:)
  real, allocatable :: SatAltDat(:)

contains

  subroutine init_mod_satellites

    if(allocated(SatTime)) return
    allocate( &
         SatTime(nMaxSats, nMaxSatInputLines), &
         SatPos(nMaxSats, 3, nMaxSatPos, nMaxSatInputLines), &
         nSatPos(nMaxSats, nMaxSatInputLines))

  end subroutine init_mod_satellites

  subroutine init_mod_satellites_rcmr
    ! Asad: Added allocation for new variables SatDat, SatCurrentDat,
    !       and SatAltDat
    ! AGB: Defined this change as a seperate routine to simplify runs
    !      with and without RCMR

    if(allocated(SatAltDat)) return
    if(allocated(SatTime)) then
       allocate( &
            SatDat(nMaxSats, nMaxSatInputLines), &
            SatCurrentDat(nMaxSats), SatAltDat(nMaxSats))
    else
       allocate( &
            SatTime(nMaxSats, nMaxSatInputLines), &
            SatPos(nMaxSats, 3, nMaxSatPos, nMaxSatInputLines), &
            nSatPos(nMaxSats, nMaxSatInputLines), &
            SatDat(nMaxSats, nMaxSatInputLines), &
            SatCurrentDat(nMaxSats), SatAltDat(nMaxSats))
    end if

  end subroutine init_mod_satellites_rcmr

end module ModSatellites
