!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!\
! ------------------------------------------------------------
! advance
! ------------------------------------------------------------
!/

subroutine advance_vertical_1d_rusanov

  use ModVertical
  use ModGITM, ONLY : Dt, iCommGITM, iProc
  use ModInputs, only: UseBarriers, iDebugLevel
  implicit none
  !-----------------------------------------------------------
  integer :: iError, iAlt
  !!!!! Variables for the Runga-Kutta 4th Order Time-stepping
  real :: OrigLogNS(-1:nAlts+2,1:nSpecies)
  real :: OrigLogINS(-1:nAlts+2,1:nIons)
  real :: OrigLogRho(-1:nAlts+2)
  real :: OrigVel_GD(-1:nAlts+2,1:3)
  real :: OrigTemp(-1:nAlts+2)
  real :: OrigVS(-1:nAlts+2,1:nSpecies)

  real :: UpdatedLogNS(-1:nAlts+2,1:nSpecies)
  real :: UpdatedLogINS(-1:nAlts+2,1:nIons)
  real :: UpdatedLogRho(-1:nAlts+2)
  real :: UpdatedVel_GD(-1:nAlts+2,1:3)
  real :: UpdatedTemp(-1:nAlts+2)
  real :: UpdatedVS(-1:nAlts+2,1:nSpecies)

  real :: FinalLogNS(-1:nAlts+2,1:nSpecies)
  real :: FinalLogINS(-1:nAlts+2,1:nIons)
  real :: FinalLogRho(-1:nAlts+2)
  real :: FinalVel_GD(-1:nAlts+2,1:3)
  real :: FinalTemp(-1:nAlts+2)
  real :: FinalVS(-1:nAlts+2,1:nSpecies)

!!! RK-4 Coefficients
  real :: K1LogNS(-1:nAlts+2,1:nSpecies)
  real :: K1LogINS(-1:nAlts+2,1:nIons)
  real :: K1LogRho(-1:nAlts+2)
  real :: K1Vel_GD(-1:nAlts+2,1:3)
  real :: K1Temp(-1:nAlts+2)
  real :: K1VS(-1:nAlts+2,1:nSpecies)

  real :: K2LogNS(-1:nAlts+2,1:nSpecies)
  real :: K2LogINS(-1:nAlts+2,1:nIons)
  real :: K2LogRho(-1:nAlts+2)
  real :: K2Vel_GD(-1:nAlts+2,1:3)
  real :: K2Temp(-1:nAlts+2)
  real :: K2VS(-1:nAlts+2,1:nSpecies)

  real :: K3LogNS(-1:nAlts+2,1:nSpecies)
  real :: K3LogINS(-1:nAlts+2,1:nIons)
  real :: K3LogRho(-1:nAlts+2)
  real :: K3Vel_GD(-1:nAlts+2,1:3)
  real :: K3Temp(-1:nAlts+2)
  real :: K3VS(-1:nAlts+2,1:nSpecies)

  real :: K4LogNS(-1:nAlts+2,1:nSpecies)
  real :: K4LogINS(-1:nAlts+2,1:nIons)
  real :: K4LogRho(-1:nAlts+2)
  real :: K4Vel_GD(-1:nAlts+2,1:3)
  real :: K4Temp(-1:nAlts+2)
  real :: K4VS(-1:nAlts+2,1:nSpecies)

  real :: DtOriginal

  real :: OldBCINS(2,1:nIons)

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 6) write(*,*) "=======> vertical bcs 1", iproc

  DtOriginal = Dt  !!! Store this so that it doesn't change

!!! =================
!!! General RK4 Update:
!!! Y(n+1) = Y(n) + Dt/6*(k1 + 2k2 + 2k3 + k4)
!!! Time(n+1) = Time(n) + Dt
!!! 
!!! k1 = f(tn,yn)
!!! k2 = f(tn + Dt/2, Yn + Dt/2*k1)
!!! k3 = f(tn + Dt/2, Yn + Dt/2*k2)
!!! k4 = f(tn + Dt, Yn + Dt*k3)
!!! =================

  ! Step 1, Fill in Ghost Cells
  call set_vertical_bcs(LogRho,LogNS,Vel_GD,Temp,LogINS,IVel,VertVel)

  ! Store our original time step from GITM (CFL-limited).

!!! Set the Original State -> Orig = Y(n)
   OrigLogNS(-1:nAlts+2,1:nSpecies)    =   LogNS(-1:nAlts+2,1:nSpecies)
  OrigLogINS(-1:nAlts+2,1:nIons) =  LogINS(-1:nAlts+2,1:nIons)
  OrigLogRho(-1:nAlts+2)               =  LogRho(-1:nAlts+2)
  OrigVel_GD(-1:nAlts+2,1:3)           =  Vel_GD(-1:nAlts+2,1:3)
    OrigTemp(-1:nAlts+2)               =    Temp(-1:nAlts+2)
      OrigVS(-1:nAlts+2,1:nSpecies)    = VertVel(-1:nAlts+2,1:nSpecies)

  NewLogNS = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
  NewTemp = Temp
  NewVertVel = VertVel


!!! Now calculate, k1 = f(tn, yn)
  call advance_vertical_1stage(&
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)

!!! note that Stage 1 -> updated by a 1/2 step
!!! (NewLogNS - LogNS) = f(tn + Dt/2, yn + dt/2)

!      Dt = DtOriginal/2.0

       K1LogNS(-1:nAlts+2,1:nSpecies)      = &
     (NewLogNS(-1:nAlts+2,1:nSpecies) - LogNS(-1:nAlts+2,1:nSpecies))

      K1LogINS(-1:nAlts+2,1:nIons)  = &
    (NewLogINS(-1:nAlts+2,1:nIons) - LogINS(-1:nAlts+2,1:nIons))

      K1LogRho(-1:nAlts+2)                = &
    (NewLogRho(-1:nAlts+2) - LogRho(-1:nAlts+2))

      K1Vel_GD(-1:nAlts+2,1:3)            = & 
    (NewVel_GD(-1:nAlts+2,1:3) - Vel_GD(-1:nAlts+2,1:3))

        K1Temp(-1:nAlts+2)                  = &
      (NewTemp(-1:nAlts+2) -  Temp(-1:nAlts+2))

          K1VS(-1:nAlts+2,1:nSpecies)      = &
   (NewVertVel(-1:nAlts+2,1:nSpecies) - VertVel(-1:nAlts+2,1:nSpecies))

  !!! Now Calculate the Next Update Stage
  !!! We need Y(Updated) = Y(n) + 0.5*K1

   UpdatedVel_GD(-1:nAlts+2,1:3) = &
      OrigVel_GD(-1:nAlts+2,1:3) + &
        0.5*K1Vel_GD(-1:nAlts+2,1:3) 

    UpdatedLogNS(-1:nAlts+2,1:nSpecies) = &
       OrigLogNS(-1:nAlts+2,1:nSpecies) +  &
         0.5*K1LogNS(-1:nAlts+2,1:nSpecies)

   UpdatedLogINS(-1:nAlts+2,1:nIons) = &
      OrigLogINS(-1:nAlts+2,1:nIons) + &
        0.5*K1LogINS(-1:nAlts+2,1:nIons) 

   UpdatedLogRho(-1:nAlts+2) = &
      OrigLogRho(-1:nAlts+2) + &
        0.5*K1LogRho(-1:nAlts+2) 

     UpdatedTemp(-1:nAlts+2) = &
        OrigTemp(-1:nAlts+2) + &
         0.5*K1Temp(-1:nAlts+2) 

       UpdatedVS(-1:nAlts+2,1:nSpecies) = &
          OrigVS(-1:nAlts+2,1:nSpecies)   + &
            0.5*K1VS(-1:nAlts+2,1:nSpecies) 

!!! UpdateStage 1 Upper Boundary
  call set_vertical_bcs(UpdatedLogRho, UpdatedLogNS, UpdatedVel_GD, &
                        UpdatedTemp, UpdatedLogINS, IVel, UpdatedVS)

   LogNS  = UpdatedLogNS
   LogINS = UpdatedLogINS
   LogRho = UpdatedLogRho
   Vel_GD = UpdatedVel_GD
     Temp = UpdatedTemp
  VertVel = UpdatedVS

  NewLogNS  = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
     NewTemp = Temp
  NewVertVel = VertVel

!!!!! Calculate K2

  call advance_vertical_1stage(&
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)

!!! K2 Coefficients for RK-4
       K2LogNS(-1:nAlts+2,1:nSpecies)      = &
     (NewLogNS(-1:nAlts+2,1:nSpecies) - LogNS(-1:nAlts+2,1:nSpecies))

      K2LogINS(-1:nAlts+2,1:nIons)  = &
    (NewLogINS(-1:nAlts+2,1:nIons) - LogINS(-1:nAlts+2,1:nIons))

      K2LogRho(-1:nAlts+2)                = &
    (NewLogRho(-1:nAlts+2) - LogRho(-1:nAlts+2))

      K2Vel_GD(-1:nAlts+2,1:3)            = & 
    (NewVel_GD(-1:nAlts+2,1:3) - Vel_GD(-1:nAlts+2,1:3))

        K2Temp(-1:nAlts+2)                  = &
      (NewTemp(-1:nAlts+2) -  Temp(-1:nAlts+2))

          K2VS(-1:nAlts+2,1:nSpecies)      = &
   (NewVertVel(-1:nAlts+2,1:nSpecies) - VertVel(-1:nAlts+2,1:nSpecies))

  !!! Now we want Y(Updated) = Y(n) + 0.5*K2

   UpdatedVel_GD(-1:nAlts+2,1:3) = &
      OrigVel_GD(-1:nAlts+2,1:3) + &
    0.5*K2Vel_GD(-1:nAlts+2,1:3) 

    UpdatedLogNS(-1:nAlts+2,1:nSpecies) = &
       OrigLogNS(-1:nAlts+2,1:nSpecies) +  &
     0.5*K2LogNS(-1:nAlts+2,1:nSpecies)

   UpdatedLogINS(-1:nAlts+2,1:nIons) = &
      OrigLogINS(-1:nAlts+2,1:nIons) + &
    0.5*K2LogINS(-1:nAlts+2,1:nIons) 

   UpdatedLogRho(-1:nAlts+2) = &
      OrigLogRho(-1:nAlts+2) + &
    0.5*K2LogRho(-1:nAlts+2) 

     UpdatedTemp(-1:nAlts+2) = &
        OrigTemp(-1:nAlts+2) + &
     0.5*K2Temp(-1:nAlts+2) 

       UpdatedVS(-1:nAlts+2,1:nSpecies) = &
          OrigVS(-1:nAlts+2,1:nSpecies)   + &
        0.5*K2VS(-1:nAlts+2,1:nSpecies) 


!! Update Boundary Conditions
  call set_vertical_bcs(UpdatedLogRho, UpdatedLogNS, UpdatedVel_GD, &
                          UpdatedTemp, UpdatedLogINS, IVel, UpdatedVS)

  LogNS  = UpdatedLogNS
  LogINS = UpdatedLogINS
  LogRho = UpdatedLogRho
  Vel_GD = UpdatedVel_GD
  Temp = UpdatedTemp
  VertVel = UpdatedVS

  NewLogNS = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
  NewTemp = Temp
  NewVertVel = VertVel
!
!
!!!!!! Calculate K3
!
  call advance_vertical_1stage(&
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)
!
!!!! K3 Coefficients for RK-4
  K3LogNS(-1:nAlts+2,1:nSpecies)      = &
     (NewLogNS(-1:nAlts+2,1:nSpecies) - LogNS(-1:nAlts+2,1:nSpecies))

  K3LogINS(-1:nAlts+2,1:nIons)  = &
     (NewLogINS(-1:nAlts+2,1:nIons) - LogINS(-1:nAlts+2,1:nIons))

  K3LogRho(-1:nAlts+2)                = &
      (NewLogRho(-1:nAlts+2) - LogRho(-1:nAlts+2))

  K3Vel_GD(-1:nAlts+2,1:3)            = & 
  (NewVel_GD(-1:nAlts+2,1:3) - Vel_GD(-1:nAlts+2,1:3))

  K3Temp(-1:nAlts+2)                  =  &
    (NewTemp(-1:nAlts+2) -  Temp(-1:nAlts+2))

  K3VS(-1:nAlts+2,1:nSpecies)      = &
     (NewVertVel(-1:nAlts+2,1:nSpecies) - VertVel(-1:nAlts+2,1:nSpecies))
  !!! Now we want Y(Updated) = Y(n) + K3


   UpdatedVel_GD(-1:nAlts+2,1:3) = &
      OrigVel_GD(-1:nAlts+2,1:3) + &
        K3Vel_GD(-1:nAlts+2,1:3) 

    UpdatedLogNS(-1:nAlts+2,1:nSpecies) = &
       OrigLogNS(-1:nAlts+2,1:nSpecies) +  &
         K3LogNS(-1:nAlts+2,1:nSpecies)

   UpdatedLogINS(-1:nAlts+2,1:nIons) = &
      OrigLogINS(-1:nAlts+2,1:nIons) + &
        K3LogINS(-1:nAlts+2,1:nIons) 

   UpdatedLogRho(-1:nAlts+2) = &
      OrigLogRho(-1:nAlts+2) + &
        K3LogRho(-1:nAlts+2) 

     UpdatedTemp(-1:nAlts+2) = &
        OrigTemp(-1:nAlts+2) + &
         K3Temp(-1:nAlts+2) 

       UpdatedVS(-1:nAlts+2,1:nSpecies) = &
          OrigVS(-1:nAlts+2,1:nSpecies)   + &
            K3VS(-1:nAlts+2,1:nSpecies) 

!!!! Update Boundary Conditions
  call set_vertical_bcs(UpdatedLogRho, UpdatedLogNS, UpdatedVel_GD, &
                          UpdatedTemp, UpdatedLogINS, IVel, UpdatedVS)

  LogNS  = UpdatedLogNS
  LogINS = UpdatedLogINS
  LogRho = UpdatedLogRho
  Vel_GD = UpdatedVel_GD
  Temp = UpdatedTemp
  VertVel = UpdatedVS

  NewLogNS = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
  NewTemp = Temp
  NewVertVel = VertVel
  
!! Calculate K4 (Final Coefficient)

  call advance_vertical_1stage(&
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)
  
!!! K4 Coefficients for RK-4
  K4LogNS(-1:nAlts+2,1:nSpecies)      = &
     (NewLogNS(-1:nAlts+2,1:nSpecies) - LogNS(-1:nAlts+2,1:nSpecies))

  K4LogINS(-1:nAlts+2,1:nIons)  = &
     (NewLogINS(-1:nAlts+2,1:nIons) - LogINS(-1:nAlts+2,1:nIons))

  K4LogRho(-1:nAlts+2)                = &
      (NewLogRho(-1:nAlts+2) - LogRho(-1:nAlts+2))

  K4Vel_GD(-1:nAlts+2,1:3)            = & 
  (NewVel_GD(-1:nAlts+2,1:3) - Vel_GD(-1:nAlts+2,1:3))

  K4Temp(-1:nAlts+2)                  = &
    (NewTemp(-1:nAlts+2) -  Temp(-1:nAlts+2))

  K4VS(-1:nAlts+2,1:nSpecies)      = &
     (NewVertVel(-1:nAlts+2,1:nSpecies) - VertVel(-1:nAlts+2,1:nSpecies))
!!!!! END Update Cycle ========================================

!!! Set the Updated State:  Stage 2
  FinalLogNS  = OrigLogNS + &
     (1.0/6.0)*(K1LogNS  + 2.0*K2LogNS  + 2.0*K3LogNS + K4LogNS)

  FinalLogINS  = OrigLogINS + &
     (1.0/6.0)*(K1LogINS  + 2.0*K2LogINS  + 2.0*K3LogINS + K4LogINS)

  FinalLogRho  = OrigLogRho + &
     (1.0/6.0)*( K1LogRho  + 2.0*K2LogRho  + 2.0*K3LogRho + K4LogRho)

  FinalVel_GD  = OrigVel_GD + &
      (1.0/6.0)*(K1Vel_GD  + 2.0*K2Vel_GD  + 2.0*K3Vel_GD + K4Vel_GD)

  FinalTemp  = OrigTemp + &
    (1.0/6.0)*(K1Temp  + 2.0*K2Temp  + 2.0*K3Temp + K4Temp)

  FinalVS  = OrigVS + &
  (1.0/6.0)*(K1VS  + 2.0*K2VS  + 2.0*K3VS + K4VS)

  call set_vertical_bcs(FinalLogRho, FinalLogNS, FinalVel_GD, &
                          FinalTemp, FinalLogINS, IVel, FinalVS)

!!! Set the Updated State:  Stage 2

   LogNS = FinalLogNS
  LogINS = FinalLogINS
  LogRho = FinalLogRho
  Vel_GD = FinalVel_GD
    Temp = FinalTemp
 VertVel = FinalVS

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 7) &
       write(*,*) "========> Done with advance_vertical_1d", iproc

end subroutine advance_vertical_1d_rusanov

!=============================================================================
subroutine advance_vertical_1stage( &
     LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
     LogINS, NewLogINS, IVel, VertVel, NewVertVel)

  ! With fluxes and sources based on LogRho..Temp, update NewLogRho..NewTemp

  use ModGITM, only: &
       Dt, iEast_, iNorth_, iUp_, ThermalDiffCoefS
  use ModPlanet
  use ModSizeGitm
  use ModVertical, only : &
       EddyCoef_1D, Centrifugal, Coriolis, &
       MeanMajorMass_1d, Gamma_1d, InvRadialDistance_C, &
       KappaTemp_1d, &
       Gravity_G, Altitude_G,Cv_1D, dAlt_F
  use ModTime
  use ModInputs
  use ModConstants
  implicit none

  real, intent(in) :: LogRho(-1:nAlts+2)
  real, intent(in) :: LogNS(-1:nAlts+2,nSpecies)
  real, intent(in) :: LogINS(-1:nAlts+2,nIons)
  real, intent(in) :: Vel_GD(-1:nAlts+2,3)
  real, intent(in) :: IVel(-1:nAlts+2,3)
  real, intent(in) :: Temp(-1:nAlts+2)
  real, intent(in) :: VertVel(-1:nAlts+2,nSpecies)

  real, intent(inout) :: NewLogRho(-1:nAlts+2)
  real, intent(inout) :: NewLogNS(-1:nAlts+2,nSpecies)
  real, intent(inout) :: NewLogINS(-1:nAlts+2,nIons)
  real, intent(inout) :: NewVel_GD(-1:nAlts+2,3)
  real, intent(inout) :: NewTemp(-1:nAlts+2)
  real, intent(out) :: NewVertVel(-1:nAlts+2,nSpecies)
  real :: NS(-1:nAlts+2,nSpecies)
  real :: NS_small(1:nAlts,nSpecies)
  real :: Rho(-1:nAlts+2)
  real :: LogNum(-1:nAlts+2)

  real, dimension(1:nAlts)    :: GradLogRho, DivVel, GradTemp, GradTempKoM, &
       DiffLogRho, DiffTemp, GradTmp, DiffTmp, DiffLogNum, GradLogNum, &
       DiviVel
  real, dimension(1:nAlts,3) :: GradVel_CD, DiffVel_CD
  real, dimension(1:nAlts,3) :: GradiVel_CD, DiffiVel_CD

  real, dimension(1:nAlts,nSpecies)    :: GradLogNS, DiffLogNS, &
       GradVertVel, DiffVertVel, DivVertVel
  real, dimension(1:nAlts,nIons) :: GradLogINS, DiffLogINS
  real :: NewSumRho, NewLogSumRho

  integer :: iAlt, iSpecies, jSpecies, iDim

  real, dimension(-1:nAlts+2)    :: NT

  real :: nVel(1:nAlts,1:nSpecies)
  integer :: nFilter, iFilter
  real :: LowFilter
!\
! Parameters Used for the Sponge
! This Sponge is useful to dampen out spurious modes
! oscillating between the bottom and top of the model.
  integer :: nAltsSponge = 12
  real :: kSP, NuSP, AmpSP
  ! JMB:  Adding Eddy Diffusion Variables here
  ! Note:  These are used in the calc_neutral_friction
  !--------------------------------------------------------------------------
  !! Eddy Diffusion Variables
  real, dimension(1:nAlts,nSpecies)    :: GradLogConS
  real, dimension(-1:nAlts+2,nSpecies)    :: ConS, LogConS
  !--------------------------------------------------------------------------
  ! 4th Order Gradients on a Non-Uniform Mesh (5-point Stencil)
  ! Used for calculating the d(ln[Chi])/dr -> Log of the concentration gradient
  !--------------------------------------------------------------------------
  real :: h1, h2, h3, h4
  real :: MeshH1, MeshH2, MeshH3, MeshH4
  real :: MeshCoef0, MeshCoef1, &
          MeshCoef2, MeshCoef3, &
          MeshCoef4
  !--------------------------------------------------------------------------
  NS = exp(LogNS)
  Rho = exp(LogRho)
  LogNum = alog(sum(NS,dim=2))
  nFilter = 10
  NT(-1:nAlts+2) = exp(LogNum(-1:nAlts+2))

  call calc_rusanov_alts_rusanov(LogRho ,GradLogRho,  DiffLogRho)
  call calc_rusanov_alts_rusanov(LogNum ,GradLogNum,  DiffLogNum)
  call calc_rusanov_alts_rusanov(Temp   ,GradTemp,    DiffTemp)
  do iDim = 1, 3
     call calc_rusanov_alts_rusanov(Vel_GD(:,iDim), &
          GradVel_CD(:,iDim),DiffVel_CD(:,iDim))
     call calc_rusanov_alts_rusanov(iVel(:,iDim), &
          GradiVel_CD(:,iDim),DiffiVel_CD(:,iDim))
  enddo

  ! Add geometrical correction to gradient and obtain divergence
  DivVel = GradVel_CD(:,iUp_) + 2*Vel_GD(1:nAlts,iUp_)*InvRadialDistance_C(1:nAlts)
  DiviVel = GradiVel_CD(:,iUp_) + 2*iVel(1:nAlts,iUp_)*InvRadialDistance_C(1:nAlts)

  do iSpecies=1,nSpecies
     call calc_rusanov_alts_rusanov(LogNS(:,iSpecies),GradTmp, DiffTmp)
     GradLogNS(:,iSpecies) = GradTmp
     DiffLogNS(:,iSpecies) = DiffTmp

     call calc_rusanov_alts_rusanov(VertVel(:,iSpecies),GradTmp, DiffTmp)
     GradVertVel(:,iSpecies) = GradTmp
     DiffVertVel(:,iSpecies) = DiffTmp
     DivVertVel(:,iSpecies) = GradVertVel(:,iSpecies) + &
          2*VertVel(1:nAlts,iSpecies)*InvRadialDistance_C(1:nAlts)
  enddo

  do iSpecies=1,nIons-1
     call calc_rusanov_alts_rusanov(LogINS(:,iSpecies), GradTmp, DiffTmp)
     GradLogINS(:,iSpecies) = GradTmp
     DiffLogINS(:,iSpecies) = DiffTmp
  enddo

  ! Use Colegrove Method
  ! Step 1:  Calculate Ln(rho_{s}/Rho)
  do iSpecies = 1, nSpecies
    LogConS(-1:nAlts+2,iSpecies) = &
         alog(Mass(iSpecies)*NS(-1:nAlts+2,iSpecies)/Rho(-1:nAlts+2))
  enddo 

  do iAlt = 1, nAlts

     ! 4th Order Concentration Gradient
     ! On Non-Uniform Mesh requires 5-pt Stencil
     h1 = dAlt_F(iAlt-1)
     h2 = dAlt_F(iAlt+0)
     h3 = dAlt_F(iAlt+1)
     h4 = dAlt_F(iAlt+2)

     MeshH2 = h2 + h1
     MeshH3 = h3 + h2 + h1
     MeshH4 = h4 + h3 + h2 + h1

     MeshCoef0 = (h2*h3*(h3+h4))/(h1*MeshH2*MeshH3*MeshH4)
     MeshCoef1 = -1.0*(MeshH2*h3*(h3 + h4))/(h1*h2*(h2+h3)*(h2+h3+h4))
     MeshCoef3 = MeshH2*h2*(h4 + h3)/(MeshH3*(h2+h3)*h3*h4) 
     MeshCoef4 = -1.0*MeshH2*h2*h3/(MeshH4*(h2+h3+h4)*(h3+h4)*h4)

     MeshCoef2 = (h2*h3*(h3+h4) + &
                  MeshH2*h3*(h3+h4) - &
                  MeshH2*h2*(h3+h4) - &
                  MeshH2*h2*h3)/&
                  (MeshH2*h2*h3*(h3+h4))

     do iSpecies = 1, nSpecies
        GradLogConS(iAlt,iSpecies) =  &
          MeshCoef0*LogConS(iAlt-2,iSpecies)&
       +  MeshCoef1*LogConS(iAlt-1,iSpecies)&
       +  MeshCoef2*LogConS(iAlt  ,iSpecies)&
       +  MeshCoef3*LogConS(iAlt+1,iSpecies)&
       +  MeshCoef4*LogConS(iAlt+2,iSpecies)
     enddo  ! iSpecies Loop
  enddo ! iAlt Loop

  AmpSP = (1.0/(10.0*Dt))
  kSP = nAltsSponge + 1

  do iAlt = 1,nAlts

     do iSpecies=1,nSpecies
        NewLogNS(iAlt,iSpecies) = LogNS(iAlt,iSpecies) - Dt * &
             (DivVertVel(iAlt,iSpecies) + &
             VertVel(iAlt,iSpecies) * GradLogNS(iAlt,iSpecies) ) + &
              Dt * DiffLogNS(iAlt,iSpecies)
     enddo

     do iSpecies=1,nIonsAdvect
        if (UseImprovedIonAdvection) then
           NewLogINS(iAlt,iSpecies) = NewLogINS(iAlt,iSpecies) - Dt * &
                (DiviVel(iAlt) * LogINS(iAlt,iSpecies) + &
                IVel(iAlt,iUp_) * GradLogINS(iAlt,iSpecies) ) &
                + Dt * DiffLogINS(iAlt,iSpecies)
        else
           NewLogINS(iAlt,iSpecies) = NewLogINS(iAlt,iSpecies) - Dt * &
                (IVel(iAlt,iUp_) * GradLogINS(iAlt,iSpecies) ) &
                + Dt * DiffLogINS(iAlt,iSpecies)
        endif
     enddo


     ! Bulk Velocity defined as the mass-weighted average of the
     ! species' velocities
     NewVel_GD(iAlt,iUp_) = 0.0

     if (iAlt >= (nAlts - nAltsSponge)) then
        NuSP = AmpSP*(1.0 - cos( pi*(kSP - (nAlts - iAlt))/kSP) )
     else
        NuSP = 0.0
     endif

     if (UseDamping) then
        VertTau(iAlt) = &
             15 - (1 - exp(-1.0*altitude_G(ialt)/1000.0/40.0))*5.0
     endif

     do iSpecies=1,nSpecies
        !The tau term was added as a vertical wind damping term
        ! Version of vertical velocity with grad(p) and g here :

       
        NewVertVel(iAlt, iSpecies) = VertVel(iAlt, iSpecies) - Dt * &
             (VertVel(iAlt,iSpecies)*GradVertVel(iAlt,iSpecies) &
             - (Vel_GD(iAlt,iNorth_)**2 + Vel_GD(iAlt,iEast_)**2) &
             * InvRadialDistance_C(iAlt) + &
             Temp(iAlt)*GradLogNS(iAlt,iSpecies) * Boltzmanns_Constant / &
             Mass(iSpecies) + &
             GradTemp(iAlt) * Boltzmanns_Constant / Mass(iSpecies) &
             - Gravity_G(iAlt)) &
             + Dt * DiffVertVel(iAlt,iSpecies) - VertVel(ialt,iSpecies)/VertTau(ialt)

        if (UseCoriolis) then
           NewVertVel(iAlt,ispecies) = NewVertVel(iAlt,ispecies) + Dt * ( &
                Centrifugal / InvRadialDistance_C(iAlt) + &
                Coriolis * Vel_GD(iAlt,iEast_))
        endif

        ! Thermal Diffusion Effects (For Light Species H2, H, and He) 
        ! ThermalDiffCoefS is set in calc_rates
        ! Note:  ThermalDiffCoefS is < 0.0 for light species
        ! This forces light species into hot zones and heavy species into cold zones
           NewVertVel(iAlt,iSpecies) = NewVertVel(iAlt,iSpecies) - &
             Dt*(ThermalDiffCoefS(iSpecies)*Boltzmanns_Constant*&
                GradTemp(iAlt))/Mass(iSpecies)

     enddo

  enddo

  ! Add Neutral Friction Between Species
  ! Needed for each increment in the RK-4 Solver
!  if (UseNeutralFriction) then
!     nVel(1:nAlts,1:nSpecies) = NewVertVel(1:nAlts,1:nSpecies)
!     call calc_neutral_friction(nVel(1:nAlts,1:nSpecies), &
!                  EddyCoef_1d(1:nAlts), NT(1:nAlts), &
!                           NS(1:nAlts,1:nSpecies), &
!                  GradLogConS(1:nAlts,1:nSpecies), &
!                         Temp(1:nAlts))
!     NewVertVel(1:nAlts,1:nSpecies) = nVel(1:nAlts,1:nSpecies)
!  endif 
!

  if (UseNeutralFriction) then
     nVel(1:nAlts,1:nSpecies) = NewVertVel(1:nAlts,1:nSpecies)
     NS_small = NS(1:nAlts,1:nSpecies)
     call calc_neutral_friction(Dt,nVel(1:nAlts,1:nSpecies), &
                  EddyCoef_1d(1:nAlts), NT(1:nAlts), &
                           NS_small, &
                         GradLogConS(1:nAlts,1:nSpecies), &
                                Temp(1:nAlts))
     NewVertVel(1:nAlts,1:nSpecies) = nVel(1:nAlts,1:nSpecies)
  endif 

  do iAlt = 1, nAlts

     do iSpecies=1,nSpecies

        NewVertVel(iAlt, iSpecies) = max(-MaximumVerticalVelocity, &
             NewVertVel(iAlt, iSpecies))
        NewVertVel(iAlt, iSpecies) = min( MaximumVerticalVelocity, &
             NewVertVel(iAlt, iSpecies))

        NewVel_GD(iAlt,iUp_) = NewVel_GD(iAlt,iUp_) + &
             NewVertVel(iAlt, iSpecies) * &
             (Mass(iSpecies) * NS(iAlt,iSpecies) / Rho(iAlt))

     enddo

  enddo

  do iAlt = 1, nAlts

     ! dVphi/dt = - (V grad V)_phi
     NewVel_GD(iAlt,iEast_) = NewVel_GD(iAlt,iEast_) - Dt * &
          Vel_GD(iAlt,iUp_)*GradVel_CD(iAlt,iEast_) &
          + Dt * DiffVel_CD(iAlt,iEast_)

     ! dVtheta/dt = - (V grad V)_theta
     NewVel_GD(iAlt,iNorth_) = NewVel_GD(iAlt,iNorth_) - Dt * &
          Vel_GD(iAlt,iUp_)*GradVel_CD(iAlt,iNorth_) &
          + Dt * DiffVel_CD(iAlt,iNorth_)

     ! dT/dt = -(V.grad T + (gamma - 1) T div V +  &
     !        (gamma - 1) * g  * grad (KeH^2  * rho) /rho 
        NewTemp(iAlt)   = NewTemp(iAlt) - Dt * &
             (Vel_GD(iAlt,iUp_)*GradTemp(iAlt) + &
             (Gamma_1d(iAlt) - 1.0) * ( &
             Temp(iAlt)*DivVel(iAlt))) &
             + Dt * DiffTemp(iAlt) 

  end do
  do iAlt = 1, nAlts
     NewSumRho    = sum( Mass(1:nSpecies)*exp(NewLogNS(iAlt,1:nSpecies)) )
     NewLogRho(iAlt) = alog(NewSumRho)
  enddo

end subroutine advance_vertical_1stage

!\
! ------------------------------------------------------------
! calc_rusanov
! ------------------------------------------------------------
!/

subroutine calc_rusanov_alts_rusanov(Var, GradVar, DiffVar)

  use ModSizeGitm
  use ModVertical, only : dAlt_C, cMax
  implicit none

  real, intent(in) :: Var(-1:nAlts+2)
  real, intent(out):: GradVar(1:nAlts), DiffVar(1:nAlts)

  real, dimension(1:nAlts+1) :: VarLeft, VarRight, DiffFlux
  !------------------------------------------------------------

  call calc_facevalues_alts_rusanov(Var, VarLeft, VarRight)

  ! Gradient based on averaged Left/Right values
  GradVar = 0.5 * &
       (VarLeft(2:nAlts+1)+VarRight(2:nAlts+1) - &
       VarLeft(1:nAlts)-VarRight(1:nAlts))/dAlt_C(1:nAlts)

  ! Rusanov/Lax-Friedrichs diffusive term
  DiffFlux = 0.5 * max(cMax(0:nAlts),cMax(1:nAlts+1)) * (VarRight - VarLeft)

  DiffVar = (DiffFlux(2:nAlts+1) - DiffFlux(1:nAlts))/dAlt_C(1:nAlts)

end subroutine calc_rusanov_alts_rusanov

!\
! ------------------------------------------------------------
! calc_facevalues_alts_rusanov
! ------------------------------------------------------------
!/

subroutine calc_facevalues_alts_rusanov(Var, VarLeft, VarRight)

  use ModVertical, only: dAlt_F, InvDAlt_F
  use ModSizeGITM, only: nAlts
  use ModLimiterGitm

  implicit none
  
  real, intent(in) :: Var(-1:nAlts+2)
  real, intent(out):: VarLeft(1:nAlts+1), VarRight(1:nAlts+1)

  real :: dVarUp, dVarDown, dVarLimited(0:nAlts+1)

  real, parameter :: Factor1=0.6250000 ! 15/24
  real, parameter :: Factor2=0.0416667 !  1/24
  real :: h

  integer :: i

  do i=1,nAlts

     ! 4th order scheme for calculating face values

     h  = InvDAlt_F(i+1)*2.0
     dVarUp   = h*(Factor1*(Var(i+1)-Var(i)  ) - Factor2*(Var(i+2)-Var(i-1)))
     h  = InvDAlt_F(i)*2.0
     dVarDown = h*(Factor1*(Var(i)  -Var(i-1)) - Factor2*(Var(i+1)-Var(i-2)))

!     ! This is Gabor's scheme
!     dVarUp            = (Var(i+1) - Var(i))   * InvDAlt_F(i+1)
!     dVarDown          = (Var(i)   - Var(i-1)) * InvDAlt_F(i)

     dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)

!     write(*,*) dVarUp, dVarDown, dVarLimited(i)

  end do

  i = 0
  dVarUp            = (Var(i+1) - Var(i))   * InvDAlt_F(i+1)
  dVarDown          = (Var(i)   - Var(i-1)) * InvDAlt_F(i)
  dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)

  i = nAlts+1
  dVarUp            = (Var(i+1) - Var(i))   * InvDAlt_F(i+1)
  dVarDown          = (Var(i)   - Var(i-1)) * InvDAlt_F(i)
  dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)

  do i=1,nAlts+1
     VarLeft(i)  = Var(i-1) + 0.5*dVarLimited(i-1) * dAlt_F(i)
     VarRight(i) = Var(i)   - 0.5*dVarLimited(i)   * dAlt_F(i)
  end do

end subroutine calc_facevalues_alts_rusanov


