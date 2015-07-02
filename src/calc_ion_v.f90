!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine calc_ion_v(iBlock)

  use ModGITM
  use ModInputs
  use ModConstants

  implicit none

  integer, intent(in) :: iBlock
  integer :: iLon, iLat, iAlt
  integer :: imax, jmax, kmax, iError, iDir
  real    :: maxi

  real, dimension(-1:nLons+2,-1:nLats+2,-1:nAlts+2) ::           &
                  B02, ForceDotB, Nie, RhoNu, IRho, &
                  VIParallel, VNParallel, gDotB, gpDotB, UDotB

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3) ::           &
                  LocPressGrad, Force, BLocal, & ! AGB: moved PressureGradient to ModGitm
                  ForceCrossB, ForcePerp

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2):: Pressure_G
  !---------------------------------------------------------------------------

  call report("Ion Forcing Terms",1)
  call start_timing("Ion Forcing")

  IVelocity(:,:,:,:,iBlock) = 0.0

  if (iDebugLevel > 4) write(*,*) "=====> pressure gradient", iproc

  Pressure_G = IPressure(:,:,:,iBlock)+ePressure(:,:,:,iBlock)
  call UAM_Gradient_GC(Pressure_G, LocPressGrad, iBlock)

  PressureGradient(:,:,:,:,iBlock) = LocPressGrad
!  PressureGradient(:,:,nAlts+1,iUp_,iBlock) = PressureGradient(:,:,nAlts,iUp_,iBlock)
!  PressureGradient(:,:,nAlts+2,iUp_,iBlock) = PressureGradient(:,:,nAlts,iUp_,iBlock)

  if (Is1D) then
     PressureGradient(:,:,:,iEast_,iBlock) = 0.0
     PressureGradient(:,:,:,iNorth_,iBlock) = 0.0
  endif

  Force = 0.0

  IRho = IDensityS(:,:,:,ie_,iBlock) * &
       MeanIonMass(:,:,:)

  if (UseIonPressureGradient) Force = Force - PressureGradient(:,:,:,:,iBlock)

  if (UseIonGravity) then
     do iAlt = -1, nAlts+2
        Force(:,:,iAlt,iUp_) = Force(:,:,iAlt,iUp_) + &
             IRho(:,:,iAlt) * Gravity_GB(:,:,iAlt,iBlock)
     enddo
  endif

  Nie = IDensityS(:,:,:,ie_,iBlock) * Element_Charge

  BLocal = B0(:,:,:,1:3,iBlock)
  B02 = B0(:,:,:,iMag_,iBlock)**2

  if (UseExB) then
     do iDir = 1, 3
        Force(:,:,:,iDir) = Force(:,:,:,iDir) + &
             Nie * EField(:,:,:,iDir)
     enddo
  endif

  RhoNu = IRho * Collisions(:,:,:,iVIN_)

  if (UseNeutralDrag) then
     do iDir = 1, 3
        Force(:,:,:,iDir) = Force(:,:,:,iDir) + &
             RhoNu * Velocity(:,:,:,iDir,iBlock)
     enddo
  endif

  ForceDotB = sum(Force * BLocal, dim=4)

  do iDir = 1, 3
     ForcePerp(:,:,:,iDir) = Force(:,:,:,iDir) - &
          Force(:,:,:,iDir) * B0(:,:,:,iDir,iBlock) / &
          B0(:,:,:,iMag_,iBlock)
  enddo

  VIParallel = 0.0
  VNParallel = 0.0

  if (maxval(blocal) == 0) then

     IVelocity(:,:,:,iUp_,iBlock) = &
          Velocity(:,:,:,iUp_,iBlock) + &
          (Gravity_GB(:,:,:, iBlock) - &
          (PressureGradient(:,:,:,iUp_,iBlock) / IRho) / &
          Collisions(:,:,:,iVIN_))

     IVelocity(:,:,:,iEast_,iBlock) = &
          Velocity(:,:,:,iEast_,iBlock) - &
          (PressureGradient(:,:,:,iEast_,iBlock) / IRho) / &
          Collisions(:,:,:,iVIN_)

     IVelocity(:,:,:,iNorth_,iBlock) = &
          Velocity(:,:,:,iNorth_,iBlock) - &
          (PressureGradient(:,:,:,iNorth_,iBlock) / IRho) / &
          Collisions(:,:,:,iVIN_)
         
  else

  UDotB = sum(Velocity(:,:,:,:,iBlock) * BLocal, dim=4)/ &
       B0(:,:,:,iMag_,iBlock)
  gpDotB = sum(PressureGradient(:,:,:,:,iBlock) * &
       BLocal, dim=4) / B0(:,:,:,iMag_,iBlock)

  do iLon = -1,nLons+2
     do iLat = -1,nLats+2
        gDotB(iLon,iLat,:) = Gravity_GB(iLon, iLat, :, iBlock) &
             * BLocal(iLon,iLat,:,iUp_) &
             /     B0(iLon,iLat,:,iMag_,iBlock)
     enddo
  enddo

  VIParallel = UDotB + &
       ( gDotB - gpDotB / IRho) / Collisions(:,:,:,iVIN_)

  ! Let's limit the Parallel Flow to something reasonable...
! AGB: Moved MaxVParallel to an input option
!  MaxVParallel = 100.0

  if (UseNeutralDrag) then
     VIParallel = min( UDotB + MaxVParallel, VIParallel)
     VIParallel = max( UDotB - MaxVParallel, VIParallel)
  else
     VIParallel = min( MaxVParallel, VIParallel)
     VIParallel = max(-MaxVParallel, VIParallel)
  endif

  ForceCrossB(:,:,:,iEast_) = &
       Force(:,:,:,iNorth_) * BLocal(:,:,:,iUp_) - &
       Force(:,:,:,iUp_)    * BLocal(:,:,:,iNorth_)

  ForceCrossB(:,:,:,iNorth_) = &
       Force(:,:,:,iUp_)    * BLocal(:,:,:,iEast_) - &
       Force(:,:,:,iEast_)  * BLocal(:,:,:,iUp_)

  ForceCrossB(:,:,:,iUp_)    = &
       Force(:,:,:,iEast_)  * BLocal(:,:,:,iNorth_) - &
       Force(:,:,:,iNorth_) * BLocal(:,:,:,iEast_)

  do iDir = 1, 3
     IVelocity(:,:,:,iDir, iBlock) = &
          VIParallel*BLocal(:,:,:,iDir)/&
          B0(:,:,:,iMag_,iBlock) + &
          ( RhoNu * ForcePerp(:,:,:,iDir) &
          + Nie * ForceCrossB(:,:,:,iDir) &
          ) / (RhoNu**2 + Nie**2 * B02)
  enddo

endif

  IVelocity(:,:,:,:,iBlock) = min( 3000.0, IVelocity(:,:,:,:,iBlock))
  IVelocity(:,:,:,:,iBlock) = max(-3000.0, IVelocity(:,:,:,:,iBlock))

  call end_timing("Ion Forcing")

  if (iDebugLevel > 4) write(*,*) "=====> done with calc_ion_v", iproc

end subroutine calc_ion_v
