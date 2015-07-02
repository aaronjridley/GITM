!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!-----------------------------------------------------------------------------
! $Id: calc_tec.f90,v 1.5 2013/10/12 04:01:00 kopmanis Exp $
!
! calc_tec
!
! Author: Angeline G. Burrell (AGB), UMichigan, December 2012
!
! Modified: Alexey Morozov (AM), UMichigan, Feb 2013 - added
!                                                      calc_single_vtec_interp
!           AGB, UMichian, Feb 2013 - Removed useless lines in
!                                     calc_single_vtec_interp and added use
!                                     of BlockLocationIndex to reduce code
!                                     duplication.
!
! Comments: Routines to compute the total electron content (TEC) from
!           the GITM electron density profiles.  TEC is defined as the integral
!           of the electron density from the receiver to the satellite.
!           VTEC (vertical TEC) integrates straight up from the ground up to
!           the GPS satellite orbit height of ~20,200 km.  As GITM reaches up
!           to ~600 km, the VTEC computed here should be less than the VTEC 
!           measured by ground-based receivers.  Extending the topside 
!           ionosphere and protonosphere would make up for much of the 
!           difference.  The plasmapause (above 2,000 km) typically contributes
!           less than 5% to the VTEC.  TEC is given in TECU, where 
!           1 TECU = 10^16 m^-2
!
!           Numerical integration is performed using the Trapazoidal Rule
!
! Contains: calc_single_vtec <- computes the VTEC at a specified latitude
!                               and longitude, using the latitude, longitude
!                               and block indeces
!           calc_vtec <- computes the VTEC at all locations
!           calc_single_vtec_interp <- computes VTEC by first interpolating 
!                                      IDensityS to lonfind and latfind
!-----------------------------------------------------------------------------

subroutine calc_single_vtec(iLon, iLat, iBlock, single_vtec)

  use ModGITM

  implicit none

  integer, intent(in) :: iLon, iLat, iBlock
  real, intent(out) :: single_vtec

  integer :: m, i
  real :: height, density

  ! Perform a simple numerical integration, using the trapazoidal rule.  VTEC
  ! is in TECU while electron density is in m^-3

  m           = nAlts - 1
  single_vtec = 0.0

  do i = 1, m
     ! Calculate height incriment and average density
     height  = Altitude_GB(iLon,iLat,i+1,iBlock)-Altitude_GB(iLon,iLat,i,iBlock)
     density = 0.5 * (IDensityS(iLon,iLat,i,ie_,iBlock) &
          + IDensityS(iLon,iLat,i+1,ie_,iBlock))

     ! Sum successive incrimentations of height * density
     single_vtec = single_vtec + height * density
  enddo

  ! Convert from SI units to TEC units

  single_vtec = single_vtec * (10.0**(-16.0))
end subroutine calc_single_vtec


subroutine calc_vtec(iBlock)

  use ModGITM

  implicit none

  integer, intent(in) :: iBlock
  integer :: iLon, iLat

  interface
     subroutine calc_single_vtec(iLon, iLat, iBlock, single_vtec)
       integer, intent(in) :: iLon, iLat, iBlock
       real, intent(out) :: single_vtec
     end subroutine calc_single_vtec
  end interface

  ! Perform a simple numerical integration, summing the electron density
  ! and multiplying it by the altitude range.  VTEC is in TECU while
  ! electron density is in m^-3

  do iLon=-1,nLons+2
     do iLat=-1,nLats+2
        call calc_single_vtec(iLon, iLat, iBlock, VTEC(iLon,iLat,iBlock))
     enddo
  enddo

  return
end subroutine calc_vtec

subroutine calc_single_vtec_interp(LonFind, LatFind, single_vtec)

  use ModGITM

  implicit none

  real, intent(in) :: LonFind, LatFind
  real, intent(out) :: single_vtec

  integer :: i, m, iiLon, iiLat, iiBlock
  real :: rLon, rLat, height
  real, dimension(0:nAlts+1) :: column !what to integrate over
  real :: Tmp(0:nLons+1,0:nLats+1,0:nAlts+1)

  ! Initialize VTEC to an unpysical value as an error flag

  single_vtec = -1.0e32

  ! Identify the latitude and longitude indexes as well as their interpolated
  ! values

  call LocationIndex(LonFind, LatFind, iiBlock, iiLon, iiLat, rLon, rLat)

  if(iiLon.lt.0 .or. iiLat.lt.0 .or. iiLon.gt.nLons .or. iiLat.gt.nLats) then
     return
  end if

  ! Interpolate electron density

  Tmp = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,ie_,iiBlock)
  call inter2(Tmp,iiLon,iiLat,rlon,rlat,column)

  ! Integrate using the trapezoidal rule

  single_vtec = 0.0
  m           = nAlts - 1

  do i = 1, m
     height      = Altitude_GB(iiLon,iiLat,i+1,iiBlock) &
          - Altitude_GB(iiLon,iiLat,i,iiBlock)
     single_vtec = single_vtec + height * (column(i) + column(i+1)) / 2.0
  enddo

  !  VTEC is in TECU while electron density is in m^-3

  single_vtec = single_vtec * (10.0**(-16.0))

contains

  subroutine inter2(variable, iiLon, iiLat, rLon, rLat, out) 

    implicit none

    real, intent(in) :: variable(0:nLons+1, 0:nLats+1, 0:nAlts+1), rLon, rLat
    integer, intent(in) :: iiLon, iiLat
    real, intent(out) :: out(0:nAlts+1)

    out = &  
          (  rLon)*(  rLat)*Variable(iiLon  ,iiLat  ,:) + &
          (1-rLon)*(  rLat)*Variable(iiLon+1,iiLat  ,:) + &
          (  rLon)*(1-rLat)*Variable(iiLon  ,iiLat+1,:) + &
          (1-rLon)*(1-rLat)*Variable(iiLon+1,iiLat+1,:)

  end subroutine inter2

end subroutine calc_single_vtec_interp
