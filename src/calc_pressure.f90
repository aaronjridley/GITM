!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine calc_pressure

  use ModGITM
  use ModPlanet
  use ModConstants
  use ModInputs, ONLY: iDebugLevel

  implicit none

  integer :: iSpecies, iBlock, iAlt

  call report("calc_pressure",2)

  NDensity = 0.0
  do iSpecies = 1, nSpecies
     NDensity(:,:,:,1:nBlocks) = &
          NDensity(:,:,:,1:nBlocks) + NDensityS(:,:,:,iSpecies,1:nBlocks)
  enddo

  Pressure    = Temperature * Rho

  IPressure = 0.0
  do iSpecies = 1, nIons-1
     IPressure(:,:,:,1:nBlocks) = IPressure(:,:,:,1:nBlocks) + &
          IDensityS(:,:,:,iSpecies,1:nBlocks) * &
          Boltzmanns_Constant * ITemperature(:,:,:,1:nBlocks)
  enddo

  ePressure(:,:,:,1:nBlocks) = &
       IDensityS(:,:,:,ie_,1:nBlocks) * Boltzmanns_Constant * &
       eTemperature(:,:,:,1:nBlocks)

  if (iDebugLevel > 2)write(*,*) &
       'calc_pressure iPres, ePres, iDens, iTemp, eTemp=',&
       sum(abs(IPressure(:,:,:,1:nBlocks))),     &
       sum(abs(ePressure(:,:,:,1:nBlocks))),     &
       sum(abs(IDensityS(:,:,:,ie_,1:nBlocks))), &
       sum(abs(ITemperature(:,:,:,1:nBlocks))),  &
       sum(abs(eTemperature(:,:,:,1:nBlocks)))

  !-------------------------------------------------------------------------
  !  B&K (1973) mixture method used for cp calculation
  !  GITM formulation for cp  seems to be identical to that of vtgcm2d
  !  cpmix = RGAS*AMU*(3.5*pco2/Mass(iCO2_) + &
  !          3.5*(pn2/Mass(iN2_)+pco/Mass(iCO_)  + &
  !          2.5*(po/Mass(iO_))
  !-------------------------------------------------------------------------
  !    
  ! The Vibration-2 is a convertion from cp to cv.  Cv should be
  ! used when in an altitude coordinate system.  Cp is for a 
  ! pressure based system.

  do iBlock = 1, nBlocks

     do iAlt = -1, nAlts+2
     
        cp(:,:,iAlt,iBlock) = 0.0
        Gamma(:,:,iAlt,iBlock) = 0.0  

        do iSpecies = 1, nSpecies
           
           cp(:,:,iAlt,iBlock) = cp(:,:,iAlt,iBlock) + &
                (Vibration(iSpecies)-2) * &
                NDensityS(1:nLons,1:nLats,iAlt,iSpecies,iBlock) * &
                (Boltzmanns_Constant / Mass(iSpecies))

           Gamma(1:nLons,1:nLats,iAlt,iBlock) = &
                gamma(1:nLons,1:nLats,iAlt,iBlock) + &
                NDensityS(1:nLons,1:nLats,iAlt,iSpecies,iBlock) / & 
                (Vibration(iSpecies)-2)    

        enddo

        cp(:,:,iAlt,iBlock) = cp(:,:,iAlt,iBlock) / &
             (2.0 * NDensity(1:nLons,1:nLats,iAlt,iBlock))

        gamma(1:nLons,1:nLats,iAlt,iBlock) = &
             gamma(1:nLons,1:nLats,iAlt,iBlock) *2.0/ &   
             NDensity(1:nLons,1:nLats,iAlt,iBlock) + 1 

     enddo

  enddo

end subroutine calc_pressure
