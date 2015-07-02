!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine init_iri

  use EUA_ModIri90, ONLY: iri90

  use ModGITM
  use ModInputs
  use ModConstants
  use ModPlanet
  use ModTime

  implicit none

  ! iri variables

  !  character, dimension(1:80) :: ccirnm, ursinm

  integer :: iBlock, iAlt, iLat, iLon, iIon

  call report("init_iri",1)

  do iBlock = 1, nBlocks
     iTemperature(:,:,:,iBlock) = 200.0
     eTemperature(:,:,:,iBlock) = 200.0
     do iAlt = -1, nAlts+2
        do iLon=-1,nLons+2
           do iLat=-1,nLats+2

              IDensityS(iLon,iLat,iAlt,:,iBlock) = 1.0

              IDensityS(iLon,iLat,iAlt,nIons,iBlock) = 0.0
              do iIon = 1, nIons-1
                 IDensityS(iLon,iLat,iAlt,nIons,iBlock) = &
                      IDensityS(iLon,iLat,iAlt,nIons,iBlock) + &
                      IDensityS(iLon,iLat,iAlt,iIon,iBlock)
              enddo

           enddo
        enddo
     enddo

  enddo

end subroutine init_iri
