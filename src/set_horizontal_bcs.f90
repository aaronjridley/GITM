!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine set_horizontal_bcs(iBlock)

  use ModSizeGitm
  use ModPlanet, only : nSpecies, nIons
  use ModInputs, only : LatStart, LatEnd, LonStart, LonEnd
  use ModGITM

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iIon, iSpecies, iLat, iLon

  call report("set_horizontal_bcs",2)

  ! Western Boundary

  if (LonStart /= LonEnd .and. &
       minval(Longitude(:,iBlock)) < LonStart) then

     do iAlt = -1, nAlts+2
        do iLon = 0,-1,-1

           do iSpecies = 1, nSpecies
              VerticalVelocity(iLon,:,iAlt,iSpecies,iBlock) = &
                   VerticalVelocity(iLon+1,:,iAlt,iSpecies,iBlock)
              nDensityS(iLon,:,iAlt,iSpecies,iBlock) = &
                   nDensityS(iLon+1,:,iAlt,iSpecies,iBlock)
!              nDensityS(iLon,:,iAlt,iSpecies,iBlock) = &
!                   2*nDensityS(iLon+1,:,iAlt,iSpecies,iBlock) - &
!                   nDensityS(iLon+2,:,iAlt,iSpecies,iBlock)
           enddo

           Rho(iLon, :,iAlt,iBlock) = &
                Rho(iLon+1,:,iAlt,iBlock)
           Temperature( iLon,:,iAlt,iBlock)= &
                Temperature(iLon+1,:,iAlt,iBlock)
           iTemperature(iLon,:,iAlt,iBlock) = &
                iTemperature(iLon+1,:,iAlt,iBlock)
           eTemperature(iLon,:,iAlt,iBlock) = &
                eTemperature(iLon+1,:,iAlt,iBlock)

!           Rho(iLon, :,iAlt,iBlock) = &
!                2*Rho(iLon+1,:,iAlt,iBlock) - &
!                Rho(iLon+2,:,iAlt,iBlock)
!           Temperature( iLon,:,iAlt,iBlock)= &
!                2*Temperature(iLon+1,:,iAlt,iBlock) - &
!                Temperature(iLon+2,:,iAlt,iBlock)
!           iTemperature(iLon,:,iAlt,iBlock) = &
!                2*iTemperature(iLon+1,:,iAlt,iBlock) - &
!                iTemperature(iLon+2,:,iAlt,iBlock)
!           eTemperature(iLon,:,iAlt,iBlock) = &
!                2*eTemperature(iLon+1,:,iAlt,iBlock) - &
!                eTemperature(iLon+2,:,iAlt,iBlock)

           do iIon = 1, nIons
              IDensityS(iLon,:,iAlt,iIon,iBlock) = &
                   IDensityS(iLon+1,:,iAlt,iIon,iBlock)
!              IDensityS(iLon,:,iAlt,iIon,iBlock) = &
!                   2*IDensityS(iLon+1,:,iAlt,iIon,iBlock) - &
!                   IDensityS(iLon+2,:,iAlt,iIon,iBlock)
           enddo

           Velocity(iLon,:,iAlt,iNorth_,iBlock) = &
                Velocity(iLon+1,:,iAlt,iNorth_,iBlock)
           Velocity(iLon,:,iAlt,iEast_,iBlock) = &
                Velocity(iLon+1,:,iAlt,iEast_,iBlock)
           Velocity(iLon,:,iAlt,iUp_,iBlock) = &
                Velocity(iLon+1,:,iAlt,iUp_,iBlock)

           iVelocity(iLon,:,iAlt,iNorth_,iBlock) = &
                iVelocity(iLon+1,:,iAlt,iNorth_,iBlock)
           iVelocity(iLon,:,iAlt,iEast_,iBlock) = &
                iVelocity(iLon+1,:,iAlt,iEast_,iBlock)
           iVelocity(iLon,:,iAlt,iUp_,iBlock) = &
                iVelocity(iLon+1,:,iAlt,iUp_,iBlock)

        enddo
     enddo
  endif

  if (LonStart /= LonEnd .and. &
       maxval(Longitude(:,iBlock)) > LonEnd) then

     do iAlt = -1, nAlts+2
        do iLon = nLons+1,nLons+2

           do iSpecies = 1, nSpecies
              VerticalVelocity(iLon,:,iAlt,iSpecies,iBlock) = &
                   VerticalVelocity(iLon-1,:,iAlt,iSpecies,iBlock)
              nDensityS(iLon,:,iAlt,iSpecies,iBlock) = &
                   nDensityS(iLon-1,:,iAlt,iSpecies,iBlock)
!              nDensityS(iLon,:,iAlt,iSpecies,iBlock) = &
!                   2*nDensityS(iLon-1,:,iAlt,iSpecies,iBlock) - &
!                   nDensityS(iLon-2,:,iAlt,iSpecies,iBlock)
           enddo

           Rho(iLon, :,iAlt,iBlock) = &
                Rho(iLon-1,:,iAlt,iBlock)
           Temperature( iLon,:,iAlt,iBlock)= &
                Temperature(iLon-1,:,iAlt,iBlock)
           iTemperature(iLon,:,iAlt,iBlock) = &
                iTemperature(iLon-1,:,iAlt,iBlock)
           eTemperature(iLon,:,iAlt,iBlock) = &
                eTemperature(iLon-1,:,iAlt,iBlock)

!           Rho(iLon, :,iAlt,iBlock) = &
!                2*Rho(iLon-1,:,iAlt,iBlock) - &
!                Rho(iLon-2,:,iAlt,iBlock)
!           Temperature( iLon,:,iAlt,iBlock)= &
!                2*Temperature(iLon-1,:,iAlt,iBlock) - &
!                Temperature(iLon-2,:,iAlt,iBlock)
!           iTemperature(iLon,:,iAlt,iBlock) = &
!                2*iTemperature(iLon-1,:,iAlt,iBlock) - &
!                iTemperature(iLon-2,:,iAlt,iBlock)
!           eTemperature(iLon,:,iAlt,iBlock) = &
!                2*eTemperature(iLon-1,:,iAlt,iBlock) - &
!                eTemperature(iLon-2,:,iAlt,iBlock)

           do iIon = 1, nIons
              IDensityS(iLon,:,iAlt,iIon,iBlock) = &
                   IDensityS(iLon-1,:,iAlt,iIon,iBlock)
!              IDensityS(iLon,:,iAlt,iIon,iBlock) = &
!                   2*IDensityS(iLon-1,:,iAlt,iIon,iBlock) - &
!                   IDensityS(iLon-2,:,iAlt,iIon,iBlock)
           enddo

           Velocity(iLon,:,iAlt,iNorth_,iBlock) = &
                Velocity(iLon-1,:,iAlt,iNorth_,iBlock)
           Velocity(iLon,:,iAlt,iEast_,iBlock) = &
                Velocity(iLon-1,:,iAlt,iEast_,iBlock)
           Velocity(iLon,:,iAlt,iUp_,iBlock) = &
                Velocity(iLon-1,:,iAlt,iUp_,iBlock)

           iVelocity(iLon,:,iAlt,iNorth_,iBlock) = &
                iVelocity(iLon-1,:,iAlt,iNorth_,iBlock)
           iVelocity(iLon,:,iAlt,iEast_,iBlock) = &
                iVelocity(iLon-1,:,iAlt,iEast_,iBlock)
           iVelocity(iLon,:,iAlt,iUp_,iBlock) = &
                iVelocity(iLon-1,:,iAlt,iUp_,iBlock)

        enddo
     enddo
  endif


  ! Southern Boundary

  if ( LatStart > -pi/2.0 .and. &
       minval(Latitude(:,iBlock)) < LatStart) then

     do iAlt = -1, nAlts+2
        do iLat = 0,-1,-1

           do iSpecies = 1, nSpecies
              VerticalVelocity(:,iLat,iAlt,iSpecies,iBlock) = &
                   VerticalVelocity(:,iLat+1,iAlt,iSpecies,iBlock)
              nDensityS(:,iLat,iAlt,iSpecies,iBlock) = &
                   nDensityS(:,iLat+1,iAlt,iSpecies,iBlock)
!              nDensityS(:,iLat,iAlt,iSpecies,iBlock) = &
!                   2*nDensityS(:,iLat+1,iAlt,iSpecies,iBlock)- &
!                   nDensityS(:,iLat+2,iAlt,iSpecies,iBlock)
           enddo

           Rho(:, iLat,iAlt,iBlock) = &
                Rho(:,iLat+1,iAlt,iBlock)
           Temperature( :, iLat,iAlt,iBlock)= &
                Temperature(:,iLat+1,iAlt,iBlock)
!           Rho(:, iLat,iAlt,iBlock) = &
!                2*Rho(:,iLat+1,iAlt,iBlock) - &
!                Rho(:,iLat+2,iAlt,iBlock)
!           Temperature( :, iLat,iAlt,iBlock)= &
!                2*Temperature(:,iLat+1,iAlt,iBlock)- &
!                Temperature(:,iLat+2,iAlt,iBlock)
           iTemperature(:, iLat,iAlt,iBlock)=iTemperature(:,iLat+1,iAlt,iBlock)
           eTemperature(:, iLat,iAlt,iBlock)=eTemperature(:,iLat+1,iAlt,iBlock)

           do iIon = 1, nIons
              IDensityS(:,iLat,iAlt,iIon,iBlock) = &
                   IDensityS(:,iLat+1,iAlt,iIon,iBlock)
           enddo

           Velocity(:,iLat,iAlt,iNorth_,iBlock) = &
                Velocity(:,iLat+1,iAlt,iNorth_,iBlock)
           Velocity(:,iLat,iAlt,iEast_,iBlock) = &
                Velocity(:,iLat+1,iAlt,iEast_,iBlock)
           Velocity(:,iLat,iAlt,iUp_,iBlock) = &
                Velocity(:,iLat+1,iAlt,iUp_,iBlock)

           iVelocity(:,iLat,iAlt,iNorth_,iBlock) = &
                iVelocity(:,iLat+1,iAlt,iNorth_,iBlock)
           iVelocity(:,iLat,iAlt,iEast_,iBlock) = &
                iVelocity(:,iLat+1,iAlt,iEast_,iBlock)
           iVelocity(:,iLat,iAlt,iUp_,iBlock) = &
                iVelocity(:,iLat+1,iAlt,iUp_,iBlock)

        enddo
     enddo
  endif

  if (LatEnd < pi/2.0 .and. &
       maxval(Latitude(:,iBlocK)) > LatEnd) then

     do iAlt = -1, nAlts+2
        do iLat = nLats+1, nLats+2

           do iSpecies = 1, nSpecies
              VerticalVelocity(:,iLat,iAlt,iSpecies,iBlock) = &
                   VerticalVelocity(:,iLat-1,iAlt,iSpecies,iBlock)
              nDensityS(:,iLat,iAlt,iSpecies,iBlock) = &
                   nDensityS(:,iLat-1,iAlt,iSpecies,iBlock)
!              nDensityS(:,iLat,iAlt,iSpecies,iBlock) = &
!                   2*nDensityS(:,iLat-1,iAlt,iSpecies,iBlock) - &
!                   nDensityS(:,iLat-2,iAlt,iSpecies,iBlock)
           enddo

           Rho(:, iLat,iAlt,iBlock) = &
                Rho(:,iLat-1,iAlt,iBlock)
           Temperature( :, iLat,iAlt,iBlock)= &
                Temperature(:,iLat-1,iAlt,iBlock)
!           Rho(:, iLat,iAlt,iBlock) = &
!                2*Rho(:,iLat-1,iAlt,iBlock) - &
!                Rho(:,iLat-2,iAlt,iBlock)
!           Temperature( :, iLat,iAlt,iBlock)= &
!                2*Temperature(:,iLat-1,iAlt,iBlock) - &
!                Temperature(:,iLat-2,iAlt,iBlock)
           iTemperature(:, iLat,iAlt,iBlock)=iTemperature(:,iLat-1,iAlt,iBlock)
           eTemperature(:, iLat,iAlt,iBlock)=eTemperature(:,iLat-1,iAlt,iBlock)

           do iIon = 1, nIons
              IDensityS(:,iLat,iAlt,iIon,iBlock) = &
                   IDensityS(:,iLat-1,iAlt,iIon,iBlock)
           enddo

           Velocity(:,iLat,iAlt,iNorth_,iBlock) = &
                Velocity(:,iLat-1,iAlt,iNorth_,iBlock)
           Velocity(:,iLat,iAlt,iEast_,iBlock) = &
                Velocity(:,iLat-1,iAlt,iEast_,iBlock)
           Velocity(:,iLat,iAlt,iUp_,iBlock) = &
                Velocity(:,iLat-1,iAlt,iUp_,iBlock)

           iVelocity(:,iLat,iAlt,iNorth_,iBlock) = &
                iVelocity(:,iLat-1,iAlt,iNorth_,iBlock)
           iVelocity(:,iLat,iAlt,iEast_,iBlock) = &
                iVelocity(:,iLat-1,iAlt,iEast_,iBlock)
           iVelocity(:,iLat,iAlt,iUp_,iBlock) = &
                iVelocity(:,iLat-1,iAlt,iUp_,iBlock)

        enddo
     enddo
  endif

end subroutine set_horizontal_bcs
