!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!------------------------------------------------------------------------------
! $Id: get_location.f90,v 1.4 2013/10/12 04:01:00 kopmanis Exp $
!
! Author: Angeline G. Burrell (AGB), UMichigan, Jan 2013
!
! Modified: AGB, UMichigan, Feb 2013 - added BlockLocationIndex
!
! LocationIndex: A routine to retireve the longitude, latitude, and block
!                indeces for a specified location.  Shamelessly stolen from
!                another place in the GITM code and put in a subroutine so that
!                it can be used in multiple places.  Exit statements were
!                added to prevent additional cycling through do-loops.
!
! Inputs: LonFind = Desired longitude
!         LatFind = Desired latitude
!
! Outputs: iiBlock = Block index containing the desired location
!          iiLon   = Longitude index for LonFind
!          iiLat   = Latitude index for LatFind
!          rLon    = Longitude interpolation scaling factor
!          rLat    = Latitude interpolation scaling factor
!
! BlockLocationIndex: A routine just like LocationIndex, but for a specified
!                     block
!
! Inputs: LonFind = Desired Longitude
!         LatFind = Desired Latitude
!         iBlock  = Block index containing the desired longitude and latitude
!
! Outputs: iiLon   = Longitude index for LonFind
!          iiLat   = Latitude index for LatFind
!          rLon    = Longitude interpolation scaling factor
!          rLat    = Latitude interpolation scaling factor
!
! BlockAltIndex: A routine similar to BlockLocationIndex, but for a specified
!                altitude
!
! Inputs: AltFind = Desired Altitude
!         iBlock  = Block index containing the desired longitude and latitude
!         iLon    = Longitude index
!         iLat    = Latitude index
!
! Outputs: iiAlt  = Altitude index for AltFind
!          rAlt   = Altitude interpolation scaling factor
!------------------------------------------------------------------------------

subroutine LocationIndex(LonFind, LatFind, iiBlock, iiLon, iiLat, rLon, rLat)

  use ModGITM

  real, intent(in) :: LonFind, LatFind
  integer, intent(out) :: iiBlock, iiLon, iiLat
  real, intent(out) :: rLon, rLat

  integer iBlock, iLon, iLat

  iiBlock = -1
  iiLon   = -1
  iiLat   = -1

  do iBlock=1,nBlocks

     if((Longitude(0,iBlock)+Longitude(1,iBlock))/2 <=LonFind .and. &
          (Longitude(nLons,iBlock)+Longitude(nLons+1,iBlock))/2 >LonFind) then

        if((Latitude(0,iBlock)+Latitude(1,iBlock))/2 <=LatFind .and. &
             (Latitude(nLats,iBlock)+Latitude(nLats+1,iBlock))/2 >LatFind) then

           iiBlock = iBlock

           do iLon = 0,nLons
              if(Longitude(iLon,iBlock) <= LonFind .and. &
                   Longitude(iLon+1,iBlock) > LonFind) then
                 iiLon = iLon
                 rLon = 1.0 - (LonFind - Longitude(iLon,iBlock)) / &
                      (Longitude(iLon+1,iBlock) - Longitude(iLon,iBlock))
                 exit
              endif
           enddo

           do iLat = 0,nLats
              if(Latitude(iLat,iBlock) <= LatFind .and. &
                   Latitude(iLat+1,iBlock) > LatFind) then
                 iiLat = iLat
                 rLat = 1.0 - (LatFind - Latitude(iLat,iBlock)) / &
                      (Latitude(iLat+1,iBlock) - Latitude(iLat,iBlock))
                 exit
              endif
           enddo

           if(iiLon >= 0 .and. iiLat >= 0) then
              exit
           end if
        end if
     end if
  end do

end subroutine LocationIndex

subroutine BlockLocationIndex(LonFind,LatFind,iBlock,iiLon,iiLat,rLon,rLat)

  use ModGITM

  real, intent(in) :: LonFind, LatFind
  integer, intent(in) :: iBlock
  integer, intent(out) :: iiLon, iiLat
  real, intent(out) :: rLon, rLat

  integer iLon, iLat

  iiLon = -1
  iiLat = -1
  rLon  = -1.0
  rLat  = -1.0

  if((Longitude(0,iBlock)+Longitude(1,iBlock))/2 <=LonFind .and. &
       (Longitude(nLons,iBlock)+Longitude(nLons+1,iBlock))/2 >LonFind) then

     if((Latitude(0,iBlock)+Latitude(1,iBlock))/2 <=LatFind .and. &
          (Latitude(nLats,iBlock)+Latitude(nLats+1,iBlock))/2 >LatFind) then

        do iLon = 0,nLons
           if(Longitude(iLon,iBlock) <= LonFind .and. &
                Longitude(iLon+1,iBlock) > LonFind) then
              iiLon = iLon
              rLon  = 1.0 - (LonFind - Longitude(iLon,iBlock)) / &
                   (Longitude(iLon+1,iBlock) - Longitude(iLon,iBlock))
              exit
           endif
        enddo

        do iLat = 0,nLats
           if(Latitude(iLat,iBlock) <= LatFind .and. &
               Latitude(iLat+1,iBlock) > LatFind) then
              iiLat = iLat
              rLat = 1.0 - (LatFind - Latitude(iLat,iBlock)) / &
                   (Latitude(iLat+1,iBlock) - Latitude(iLat,iBlock))
              exit
           endif
        enddo
     end if
  end if

end subroutine BlockLocationIndex

subroutine BlockAltIndex(AltFind,iBlock,iLon,iLat,iAlt,rAlt)

  use ModGITM

  real, intent(in)     :: AltFind
  integer, intent(in)  :: iBlock, iLon, iLat
  integer, intent(out) :: iAlt
  real, intent(out)    :: rAlt

  integer jAlt

  iAlt = -1
  rAlt = -1.0

  do jAlt = 0,nAlts
     if (Altitude_GB(iLon, iLat, jAlt, iBlock) <= AltFind .and. &
          Altitude_GB(iLon, iLat, jAlt+1,iBlock) > AltFind) then
        iAlt = jAlt
        rAlt  = 1.0 - (AltFind - Altitude_GB(iLon, iLat, iAlt, iBlock)) &
             / (Altitude_GB(iLon, iLat, iAlt+1, iBlock) &
             - Altitude_GB(iLon, iLat, iAlt, iBlock))
        exit
     endif
  enddo

end subroutine BlockAltIndex
