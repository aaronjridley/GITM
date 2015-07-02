!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine set_bcs

  use ModGITM
  use ModInputs

  implicit none

  integer :: iBlock, iLon, iLat, iSpecies

  do iBlock = 1, nBlocks

     ! Bottom
     Velocity(:,:,-1,iUp_,iBlock) = - Velocity(:,:,2,iUp_,iBlock)
     Velocity(:,:, 0,iUp_,iBlock) = - Velocity(:,:,1,iUp_,iBlock)

     ! Fixed LogRho and Temp

     ! if we don't touch LogRho(:,:,-1:0,iBlock) it will never change

     ! Top

     do iLon = 1, nLons
        do iLat = 1, nLats

           if(Velocity(iLon,iLat,nAlts,iUp_,iBlock)>0.)then

              Velocity(iLon,iLat,nAlts+1,iUp_,iBlock)  = &
                   Velocity(iLon,iLat,nAlts,iUp_,iBlock)*0.0
              Velocity(iLon,iLat,nAlts+2,iUp_,iBlock)  = &
                   Velocity(iLon,iLat,nAlts,iUp_,iBlock)*0.0

           else

              Velocity(iLon,iLat,nAlts+1,iUp_,iBlock) = 0.0 ! -Vel(nAlts)
              Velocity(iLon,iLat,nAlts+2,iUp_,iBlock) = 0.0 ! -Vel(nAlts-1)

           endif
        enddo
     enddo

!     Temperature(:,:,nAlts+1,iBlock) = TempMax/TempUnit
!     Temperature(:,:,nAlts+2,iBlock) = TempMax/TempUnit

     Temperature(:,:,nAlts+1,iBlock) = Temperature(:,:,nAlts,iBlock)
     Temperature(:,:,nAlts+2,iBlock) = Temperature(:,:,nAlts,iBlock)


     do iSpecies=1,nSpecies
        LogNS(:,:,nAlts+1,iSpecies,iBlock) = &
             LogNS(:,:,nAlts,iSpecies,iBlock) + &
             dAlt(nAlts) * Gravity(nAlts)/Temperature(:,:,nAlts,iBlock)
        LogNS(:,:,nAlts+1,iSpecies,iBlock) = &
             LogNS(:,:,nAlts,iSpecies,iBlock) + &
             2*dAlt(nAlts) * Gravity(nAlts)/Temperature(:,:,nAlts,iBlock)
     enddo

     LogRho(:,:,nAlts+1,iBlock) = &
          LogRho(:,:,nAlts,iBlock) + &
          dAlt(nAlts) * Gravity(nAlts)/Temperature(:,:,nAlts,iBlock)
     LogRho(:,:,nAlts+1,iBlock) = &
          LogRho(:,:,nAlts,iBlock) + &
          2*dAlt(nAlts) * Gravity(nAlts)/Temperature(:,:,nAlts,iBlock)

  enddo

end subroutine set_bcs
