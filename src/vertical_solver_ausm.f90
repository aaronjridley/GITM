!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!\
! ------------------------------------------------------------
! advance
! ------------------------------------------------------------
!/

subroutine check_ion_densities(iDen)

  use ModSizeGitm
  use ModPlanet, only: nIonsAdvect, nIons
  real, intent(inout) :: iDen(-1:nAlts+2,nIons)

  do iIon = 1, nIonsAdvect
     do iAlt = 1, nAlts+2
        if (iDen(iAlt, iIon) < 0.0) then
!           write(*,*) 'Found negative ion density: ', &
!                iAlt, iIon, iDen(iAlt,iIon)
           iDen(iAlt,iIon) = max(iDen(iAlt-1,iIon)*0.99,10.0)
        endif
     enddo
  enddo

end subroutine check_ion_densities

subroutine advance_vertical_1d_ausm

  use ModVertical
  use ModPlanet, only: Mass
  use ModGITM, ONLY : Dt, iCommGITM, iProc, iUp_
  use ModInputs, only: UseBarriers, iDebugLevel, UseImprovedIonAdvection
  implicit none
  !-----------------------------------------------------------
  integer :: iError, iAlt, iSpecies, iDir
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

!!! RStage-4 Coefficients
  real :: Stage1LogNS(-1:nAlts+2,1:nSpecies)
  real :: Stage1LogINS(-1:nAlts+2,1:nIons)
  real :: Stage1LogRho(-1:nAlts+2)
  real :: Stage1Vel_GD(-1:nAlts+2,1:3)
  real :: Stage1Temp(-1:nAlts+2)
  real :: Stage1VS(-1:nAlts+2,1:nSpecies)

  real :: Stage2LogNS(-1:nAlts+2,1:nSpecies)
  real :: Stage2LogINS(-1:nAlts+2,1:nIons)
  real :: Stage2LogRho(-1:nAlts+2)
  real :: Stage2Vel_GD(-1:nAlts+2,1:3)
  real :: Stage2Temp(-1:nAlts+2)
  real :: Stage2VS(-1:nAlts+2,1:nSpecies)

  real :: Stage3LogNS(-1:nAlts+2,1:nSpecies)
  real :: Stage3LogINS(-1:nAlts+2,1:nIons)
  real :: Stage3LogRho(-1:nAlts+2)
  real :: Stage3Vel_GD(-1:nAlts+2,1:3)
  real :: Stage3Temp(-1:nAlts+2)
  real :: Stage3VS(-1:nAlts+2,1:nSpecies)

  real :: Stage4LogNS(-1:nAlts+2,1:nSpecies)
  real :: Stage4LogINS(-1:nAlts+2,1:nIons)
  real :: Stage4LogRho(-1:nAlts+2)
  real :: Stage4Vel_GD(-1:nAlts+2,1:3)
  real :: Stage4Temp(-1:nAlts+2)
  real :: Stage4VS(-1:nAlts+2,1:nSpecies)

  real :: DtOriginal
  real :: DtIn

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

   Stage1LogNS(-1:nAlts+2,1:nSpecies)    =   LogNS(-1:nAlts+2,1:nSpecies)
  Stage1LogINS(-1:nAlts+2,1:nIons) =  LogINS(-1:nAlts+2,1:nIons)
  Stage1LogRho(-1:nAlts+2)               =  LogRho(-1:nAlts+2)
  Stage1Vel_GD(-1:nAlts+2,1:3)           =  Vel_GD(-1:nAlts+2,1:3)
    Stage1Temp(-1:nAlts+2)               =    Temp(-1:nAlts+2)
      Stage1VS(-1:nAlts+2,1:nSpecies)    = VertVel(-1:nAlts+2,1:nSpecies)

  NewLogNS = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
  NewTemp = Temp
  NewVertVel = VertVel

  DtIn = 0.5*Dt  !!! Store this so that it doesn't change
  call advance_vertical_1stage_ausm(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)

!!! note that Stage 1 -> updated by a 1/2 step
!!! (NewLogNS - LogNS) = f(tn + Dt/2, yn + dt/2)
       Stage2LogNS(-1:nAlts+2,1:nSpecies)      = &
       Stage1LogNS(-1:nAlts+2,1:nSpecies)      + &
     (NewLogNS(-1:nAlts+2,1:nSpecies) - LogNS(-1:nAlts+2,1:nSpecies))

      Stage2LogINS(-1:nAlts+2,1:nIons)  = &
      Stage1LogINS(-1:nAlts+2,1:nIons)  + &
    (NewLogINS(-1:nAlts+2,1:nIons) - LogINS(-1:nAlts+2,1:nIons))

      if (UseImprovedIonAdvection) call check_ion_densities(Stage2LogINS)

      Stage2LogRho(-1:nAlts+2)                = &
      Stage1LogRho(-1:nAlts+2)                + &
    (NewLogRho(-1:nAlts+2) - LogRho(-1:nAlts+2))

      Stage2Vel_GD(-1:nAlts+2,1:3)            = & 
      Stage1Vel_GD(-1:nAlts+2,1:3)            + & 
    (NewVel_GD(-1:nAlts+2,1:3) - Vel_GD(-1:nAlts+2,1:3))

        Stage2Temp(-1:nAlts+2)                  = &
        Stage1Temp(-1:nAlts+2)                  + &
     (NewTemp(-1:nAlts+2) -  Temp(-1:nAlts+2))

          Stage2VS(-1:nAlts+2,1:nSpecies)      = &
          Stage1VS(-1:nAlts+2,1:nSpecies)      + &
   (NewVertVel(-1:nAlts+2,1:nSpecies) - VertVel(-1:nAlts+2,1:nSpecies))

  !!! Now Calculate the Next Update Stage
  !!! We need Y(Updated) = Y(n) + 0.5*Stage2
   UpdatedVel_GD(-1:nAlts+2,1:3) = &
        Stage2Vel_GD(-1:nAlts+2,1:3) 

    UpdatedLogNS(-1:nAlts+2,1:nSpecies) = &
         Stage2LogNS(-1:nAlts+2,1:nSpecies)

   UpdatedLogINS(-1:nAlts+2,1:nIons) = &
         Stage2LogINS(-1:nAlts+2,1:nIons) 

   UpdatedLogRho(-1:nAlts+2) = &
         Stage2LogRho(-1:nAlts+2) 

     UpdatedTemp(-1:nAlts+2) = &
          Stage2Temp(-1:nAlts+2) 

       UpdatedVS(-1:nAlts+2,1:nSpecies) = &
            Stage2VS(-1:nAlts+2,1:nSpecies) 

  DtIn = 0.5*Dt
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

  call advance_vertical_1stage_ausm(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)

       Stage3LogNS(-1:nAlts+2,1:nSpecies)      = &
       Stage1LogNS(-1:nAlts+2,1:nSpecies)      + &
     (NewLogNS(-1:nAlts+2,1:nSpecies) - LogNS(-1:nAlts+2,1:nSpecies))

      Stage3LogINS(-1:nAlts+2,1:nIons)  = &
      Stage1LogINS(-1:nAlts+2,1:nIons)  + &
    (NewLogINS(-1:nAlts+2,1:nIons) - LogINS(-1:nAlts+2,1:nIons))

      if (UseImprovedIonAdvection) call check_ion_densities(Stage3LogINS)

      Stage3LogRho(-1:nAlts+2)                = &
      Stage1LogRho(-1:nAlts+2)                + &
    (NewLogRho(-1:nAlts+2) - LogRho(-1:nAlts+2))

      Stage3Vel_GD(-1:nAlts+2,1:3)            = & 
      Stage1Vel_GD(-1:nAlts+2,1:3)            + & 
    (NewVel_GD(-1:nAlts+2,1:3) - Vel_GD(-1:nAlts+2,1:3))

        Stage3Temp(-1:nAlts+2)                  = &
        Stage1Temp(-1:nAlts+2)                  + &
      (NewTemp(-1:nAlts+2) -  Temp(-1:nAlts+2))

          Stage3VS(-1:nAlts+2,1:nSpecies)      = &
          Stage1VS(-1:nAlts+2,1:nSpecies)      + &
   (NewVertVel(-1:nAlts+2,1:nSpecies) - VertVel(-1:nAlts+2,1:nSpecies))

  !!! Now Calculate the Next Update Stage
  !!! We need Y(Updated) = Y(n) + 0.5*Stage3
   UpdatedVel_GD(-1:nAlts+2,1:3) = &
        Stage3Vel_GD(-1:nAlts+2,1:3) 

    UpdatedLogNS(-1:nAlts+2,1:nSpecies) = &
         Stage3LogNS(-1:nAlts+2,1:nSpecies)

   UpdatedLogINS(-1:nAlts+2,1:nIons) = &
         Stage3LogINS(-1:nAlts+2,1:nIons) 

   UpdatedLogRho(-1:nAlts+2) = &
         Stage3LogRho(-1:nAlts+2) 

     UpdatedTemp(-1:nAlts+2) = &
          Stage3Temp(-1:nAlts+2) 

       UpdatedVS(-1:nAlts+2,1:nSpecies) = &
            Stage3VS(-1:nAlts+2,1:nSpecies) 

  DtIn = Dt
!! Update Boundary Conditions

  call set_vertical_bcs(UpdatedLogRho, UpdatedLogNS, UpdatedVel_GD, &
                          UpdatedTemp, UpdatedLogINS, IVel, UpdatedVS)

!
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
  call advance_vertical_1stage_ausm(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)
!!!! K3 Coefficients for RK-4

       Stage4LogNS(-1:nAlts+2,1:nSpecies)      = &
       Stage1LogNS(-1:nAlts+2,1:nSpecies)      + &
     (NewLogNS(-1:nAlts+2,1:nSpecies) - LogNS(-1:nAlts+2,1:nSpecies))

      Stage4LogINS(-1:nAlts+2,1:nIons)  = &
      Stage1LogINS(-1:nAlts+2,1:nIons)  + &
    (NewLogINS(-1:nAlts+2,1:nIons) - LogINS(-1:nAlts+2,1:nIons))

      if (UseImprovedIonAdvection) call check_ion_densities(Stage4LogINS)

      Stage4LogRho(-1:nAlts+2)                = &
      Stage1LogRho(-1:nAlts+2)                + &
    (NewLogRho(-1:nAlts+2) - LogRho(-1:nAlts+2))

      Stage4Vel_GD(-1:nAlts+2,1:3)            = & 
      Stage1Vel_GD(-1:nAlts+2,1:3)            + & 
    (NewVel_GD(-1:nAlts+2,1:3) - Vel_GD(-1:nAlts+2,1:3))

        Stage4Temp(-1:nAlts+2)                  = &
        Stage1Temp(-1:nAlts+2)                  + &
      (NewTemp(-1:nAlts+2) -  Temp(-1:nAlts+2))

          Stage4VS(-1:nAlts+2,1:nSpecies)      = &
          Stage1VS(-1:nAlts+2,1:nSpecies)      + &
   (NewVertVel(-1:nAlts+2,1:nSpecies) - VertVel(-1:nAlts+2,1:nSpecies))

  !!! Now Calculate the Next Update Stage
  !!! We need Y(Updated) = Y(n) + 0.5*Stage4
   UpdatedVel_GD(-1:nAlts+2,1:3) = &
        Stage4Vel_GD(-1:nAlts+2,1:3) 

    UpdatedLogNS(-1:nAlts+2,1:nSpecies) = &
         Stage4LogNS(-1:nAlts+2,1:nSpecies)

   UpdatedLogINS(-1:nAlts+2,1:nIons) = &
         Stage4LogINS(-1:nAlts+2,1:nIons) 

   UpdatedLogRho(-1:nAlts+2) = &
         Stage4LogRho(-1:nAlts+2) 

     UpdatedTemp(-1:nAlts+2) = &
          Stage4Temp(-1:nAlts+2) 

       UpdatedVS(-1:nAlts+2,1:nSpecies) = &
          Stage4VS(-1:nAlts+2,1:nSpecies) 


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

  DtIn = 0.5*Dt
  call advance_vertical_1stage_ausm(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)

  ! This section ensures that our lower boundary conditions are maintained
  ! and not overwritten.
  FinalLogNS = Stage1LogNS
  FinalLogINS = Stage1LogINS
  FinalLogRho = Stage1LogRho
  FinalVel_GD = Stage1Vel_GD
  FinalTemp   = Stage1Temp
  FinalVS     = Stage1VS

!!! Set the Updated State:  Stage 2
  FinalLogNS(1:nAlts,:) = &
     (1.0/3.0)*(-1.0*Stage1LogNS(1:nAlts,:) + Stage2LogNS(1:nAlts,:)  +  &
                 2.0*Stage3LogNS(1:nAlts,:) + Stage4LogNS(1:nAlts,:)  +  &
                       (NewLogNS(1:nAlts,:) - LogNS(1:nAlts,:)) )

  FinalLogINS(1:nAlts,:) = &
     (1.0/3.0)*(-1.0*Stage1LogINS(1:nAlts,:) + Stage2LogINS(1:nAlts,:)  +  &
                 2.0*Stage3LogINS(1:nAlts,:) + Stage4LogINS(1:nAlts,:)  +  &
                       (NewLogINS(1:nAlts,:) - LogINS(1:nAlts,:)) )

  if (UseImprovedIonAdvection) call check_ion_densities(FinalLogINS)

  FinalLogRho(1:nAlts) = &
     (1.0/3.0)*(-1.0*Stage1LogRho(1:nAlts) + Stage2LogRho(1:nAlts)  +  &
                 2.0*Stage3LogRho(1:nAlts) + Stage4LogRho(1:nAlts)  +  &
                       (NewLogRho(1:nAlts) - LogRho(1:nAlts)) )

  FinalVel_GD(1:nAlts,:) = &
     (1.0/3.0)*(-1.0*Stage1Vel_GD(1:nAlts,:) + Stage2Vel_GD(1:nAlts,:)  +  &
                 2.0*Stage3Vel_GD(1:nAlts,:) + Stage4Vel_GD(1:nAlts,:)  +  &
                       (NewVel_GD(1:nAlts,:) - Vel_GD(1:nAlts,:)) )

  FinalTemp(1:nAlts) = &
     (1.0/3.0)*(-1.0*Stage1Temp(1:nAlts) + Stage2Temp(1:nAlts)  +  &
                 2.0*Stage3Temp(1:nAlts) + Stage4Temp(1:nAlts)  +  &
                       (NewTemp(1:nAlts) - Temp(1:nAlts)) )

  FinalVS(1:nAlts,:) = &
     (1.0/3.0)*(-1.0*Stage1VS(1:nAlts,:) + Stage2VS(1:nAlts,:)  +  &
                 2.0*Stage3VS(1:nAlts,:) + Stage4VS(1:nAlts,:)  +  &
                  (NewVertVel(1:nAlts,:) - VertVel(1:nAlts,:)) )

  call set_vertical_bcs(FinalLogRho, FinalLogNS, FinalVel_GD, &
                          FinalTemp, FinalLogINS, IVel, FinalVS)

   LogNS = FinalLogNS
  LogINS = FinalLogINS
  LogRho = FinalLogRho
  Vel_GD = FinalVel_GD
    Temp = FinalTemp
 VertVel = FinalVS

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 7) &
       write(*,*) "========> Done with advance_vertical_1d", iproc

end subroutine advance_vertical_1d_ausm

!=============================================================================
subroutine advance_vertical_1stage_ausm( DtIn, &
     LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
     LogINS, NewLogINS, IVel, VertVel, NewVertVel)

  ! With fluxes and sources based on LogRho..Temp, update NewLogRho..NewTemp

  use ModGITM, only: &
       Dt, iEast_, iNorth_, iUp_, ThermalDiffCoefS
  use ModPlanet
  use ModSizeGitm
  use ModVertical, only : &
       Heating, EddyCoef_1D, ViscCoef_1d, Coriolis, &
       MeanMajorMass_1d, Gamma_1d, InvRadialDistance_C, &
       KappaTemp_1d, ViscCoefS_1d, &
       !ChemSources_1d, &
       Centrifugal, &
       Gravity_G, Altitude_G,Cv_1D, dAlt_F
  use ModTime
  use ModInputs
  use ModConstants
  use ModSources, only : EddyCondAdia
  implicit none

  real, intent(in) :: DtIn
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
  real :: NewVel2_G(-1:nAlts+2)
  real, intent(inout) :: NewTemp(-1:nAlts+2)
  real, intent(out) :: NewVertVel(-1:nAlts+2,nSpecies)
  real :: NS(-1:nAlts+2,nSpecies), Pressure1D(-1:nAlts+2)
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
  real :: NewSumRho, NewLogSumRho, rat, ed

  integer :: iAlt, iSpecies, jSpecies, iDim

  real, dimension(-1:nAlts+2)    :: NT
  real, dimension(-1:nAlts+2)    :: Press, LogPress
  real, dimension(1:nAlts)    :: DiffLogPress, GradLogPress
  real, dimension(1:nAlts,nSpecies)    :: EddyDiffusionVel

  real :: nVel(1:nAlts,1:nSpecies)
  integer :: nFilter, iFilter
  real :: LowFilter

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
  real, dimension(1:nAlts,nSpecies)    :: EddyCoefRatio_1d
  !--------------------------------------------------------------------------
  ! 4th Order Gradients on a Non-Uniform Mesh (5-point Stencil)
  ! Used for calculating the d(ln[Chi])/dr -> Log of the concentration gradient
  !--------------------------------------------------------------------------
  real :: h1, h2, h3, h4
  real :: MeshH1, MeshH2, MeshH3, MeshH4
  real :: MeshCoef0, MeshCoef1, &
          MeshCoef2, MeshCoef3, &
          MeshCoef4
  ! ----------------------------------------------------
  ! JMB:  AUSM Variables
  ! ----------------------------------------------------
  real ::    RhoS(-1:nAlts+2,1:nSpecies),&
          NewRhoS(-1:nAlts+2,1:nSpecies),&
      AUSMRhoSFluxes(1:nAlts,1:nSpecies)

  real :: InvScaleHeight, MeanGravity, MeanTemp, MeanMass

  real ::   HydroNS(-1:nAlts+2,1:nSpecies),&
          HydroRhoS(-1:nAlts+2,1:nSpecies), &
            HydroNT(-1:nAlts+2),&
     HydroPressureS(-1:nAlts+2,1:nSpecies),&
      DeviationRhoS(-1:nAlts+2,1:nSpecies),&
 DeviationRhoSRatio(-1:nAlts+2,1:nSpecies),&
      HydroPressure(-1:nAlts+2), &
           HydroRho(-1:nAlts+2)

  ! Momentum Fluxes and Variables
  real ::    MomentumS(-1:nAlts+2,1:nSpecies),&
          NewMomentumS(-1:nAlts+2,1:nSpecies),&
              Momentum(-1:nAlts+2,1:3),&        ! Bulk Momentum
           NewMomentum(-1:nAlts+2,1:3),&        ! Bulk Momentum
       AUSMMomentumSFluxes(1:nAlts,1:nSpecies), &
       AUSMMomentumFluxes(1:nAlts,3)

  real :: PressureS(-1:nAlts+2,1:nSpecies), &
          NewNS(-1:nAlts+2,1:nSpecies), &
          NewNT(-1:nAlts+2)

  real :: NS_small(nAlts,nSpecies)
  
  real :: TotalEnergy(-1:nAlts+2),&
       NewTotalEnergy(-1:nAlts+2),&
       AUSMTotalEnergyFluxes(1:nAlts), &
       NewPress(-1:nAlts+2), NewRho(-1:nAlts+2)

  logical :: SubtractHydrostatic(-1:nAlts+2,1:nSpecies)   
  real :: RadialDistance_C(-1:nAlts+2)
  real :: EffectiveGravity(-1:nAlts+2)

  ! JMB:  Use these as Limiters on Winds for an initial startup
  real :: TimeFactor, Vel0, DeltaV, VelocityCap

  ! ----------------------------------------------------
  ! JMB: 06/27/2016:---- THERMAL CONDUCTION & VISCOSITY
  ! ----------------------------------------------------
  real, dimension(1:nAlts)    :: VertVisc
  real, dimension(1:nAlts)    :: ThermalCond


  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  Vel0 = 100.0 ! initial velocity in m/s
  TimeFactor = exp(-tSimulation/43200.0)
  DeltaV = MaximumVerticalVelocity - Vel0
  VelocityCap = Vel0 + DeltaV*(1.0 - TimeFactor)

  do iAlt = -1, nAlts + 2
  !   EffectiveGravity(iAlt) = &
  !      Gravity_G(iAlt) 

     EffectiveGravity(iAlt) = &
        Gravity_G(iAlt) + &
        Centrifugal / InvRadialDistance_C(iAlt) 
  enddo 

  NS = exp(LogNS)
  Rho = exp(LogRho)
  LogNum = alog(sum(NS,dim=2))
  nFilter = 10
  do iAlt = -1, nAlts + 2
     RadialDistance_C(iAlt) = 1.0/InvRadialDistance_C(iAlt)
  enddo 
  
  NT(-1:nAlts+2) = 0.0
  do iSpecies = 1, nSpecies
     NT(-1:nAlts+2) = NT(-1:nAlts+2) + &
                      NS(-1:nAlts+2,iSpecies)
  enddo 

  do iAlt = -1, nAlts + 2
    Press(iAlt) = NT(iAlt)*Boltzmanns_Constant*Temp(iAlt)
    LogPress(iAlt) = alog(Press(iAlt))
  enddo

  call calc_rusanov_alts_ausm(LogPress ,GradLogPress,  DiffLogPress)
  call calc_rusanov_alts_ausm(LogRho ,GradLogRho,  DiffLogRho)
  call calc_rusanov_alts_ausm(LogNum ,GradLogNum,  DiffLogNum)
  call calc_rusanov_alts_ausm(Temp   ,GradTemp,    DiffTemp)
  do iDim = 1, 3
     call calc_rusanov_alts_ausm(Vel_GD(:,iDim), &
          GradVel_CD(:,iDim),DiffVel_CD(:,iDim))
     call calc_rusanov_alts_ausm(iVel(:,iDim), &
          GradiVel_CD(:,iDim),DiffiVel_CD(:,iDim))
  enddo

  ! Add geometrical correction to gradient and obtain divergence
  DivVel = GradVel_CD(:,iUp_) + 2*Vel_GD(1:nAlts,iUp_)*InvRadialDistance_C(1:nAlts)
  DiviVel = GradiVel_CD(:,iUp_) + 2*iVel(1:nAlts,iUp_)*InvRadialDistance_C(1:nAlts)

  do iSpecies=1,nSpecies

     call calc_rusanov_alts_ausm(LogNS(:,iSpecies),GradTmp, DiffTmp)
     GradLogNS(:,iSpecies) = GradTmp
     DiffLogNS(:,iSpecies) = DiffTmp

     call calc_rusanov_alts_ausm(VertVel(:,iSpecies),GradTmp, DiffTmp)
     GradVertVel(:,iSpecies) = GradTmp
     DiffVertVel(:,iSpecies) = DiffTmp
     DivVertVel(:,iSpecies) = GradVertVel(:,iSpecies) + &
          2*VertVel(1:nAlts,iSpecies)*InvRadialDistance_C(1:nAlts)

  enddo

  do iSpecies=1,nIons-1
     call calc_rusanov_alts_ausm(LogINS(:,iSpecies), GradTmp, DiffTmp)
     GradLogINS(:,iSpecies) = GradTmp
     DiffLogINS(:,iSpecies) = DiffTmp
  enddo

  ! Step 1:  Calculate Ln(rho_{s}/Rho)
  do iSpecies = 1, nSpecies
    LogConS(-1:nAlts+2,iSpecies) = &
         alog(Mass(iSpecies)*NS(-1:nAlts+2,iSpecies)/Rho(-1:nAlts+2))
    !LogConS(-1:nAlts+2,iSpecies) = &
    !     alog(NS(-1:nAlts+2,iSpecies)/NT(-1:nAlts+2))
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

  !!!! JMB AUSM: BEGIN THE HYDROSTATIC BACKGROUND
  !!!! We define a hydrostatic background to subtract off
  !!!! this removes errors introduced in the Grad(P) - rho*g
  !!!! that reduces the accuracy of non-hydrostatic calculations.
  HydroNS(-1:nAlts+2,1:nSpecies) = NS(-1:nAlts+2,1:nSpecies)
  HydroNT(-1:nAlts+2)            = NT(-1:nAlts+2)
  !!! Calculate the Bulk Hydrostatic Density (NT)
  do iAlt = 1, nAlts + 2
    MeanMass = 0.5*(MeanMajorMass_1d(iAlt-1) + MeanMajorMass_1d(iAlt))
    MeanGravity = 0.5*( EffectiveGravity(iAlt) + EffectiveGravity(iAlt-1) )
    MeanTemp = 0.5*( Temp(iAlt) + Temp(iAlt-1) )
    InvScaleHeight = -1.0* MeanMass*MeanGravity/&
                     (Boltzmanns_Constant*MeanTemp)
    HydroNT(iAlt) = HydroNT(iAlt-1)*(Temp(iAlt-1)/Temp(iAlt))*&
          exp (-1.0*dAlt_F(iAlt)*InvScaleHeight)
  enddo 
  iAlt = -1
  MeanMass = 0.5*( MeanMajorMass_1d(-1) + MeanMajorMass_1d(0))
  MeanGravity = 0.5*( EffectiveGravity(iAlt+1) + EffectiveGravity(iAlt) )
  MeanTemp = 0.5*( Temp(iAlt+1) + Temp(iAlt) )
  InvScaleHeight = -1.0* MeanMass*MeanGravity/&
                  (Boltzmanns_Constant*MeanTemp)
  HydroNT(iAlt) = HydroNT(iAlt+1)*(Temp(iAlt+1)/Temp(iAlt))*&
                  exp(dAlt_F(iAlt)*InvScaleHeight)
  !!! Calculate a Background Atmosphere for each species
  !!! Note: This really only works for the major species
  !!! which is fine, since minor species will have large difference
  !!! Grad(P) - rho*g terms and are not subject to the same level
  !!! of sensitivity.
  do iSpecies =  1, nSpecies
     do iAlt = 1, nAlts + 2
     MeanMass = Mass(iSpecies)
     MeanGravity = 0.5*( EffectiveGravity(iAlt) + EffectiveGravity(iAlt-1) )
     MeanTemp = 0.5*( Temp(iAlt) + Temp(iAlt-1) )
     InvScaleHeight = -1.0* MeanMass*MeanGravity/&
                      (Boltzmanns_Constant*MeanTemp)
     HydroNS(iAlt,iSpecies) = &
           HydroNS(iAlt-1,iSpecies)*(Temp(iAlt-1)/Temp(iAlt))*&
           exp (-1.0*dAlt_F(iAlt)*InvScaleHeight)
     enddo 
  enddo 
  ! Extend the densities downward with the mean mass of the atmosphere
  do iSpecies =  1, nSpecies
     iAlt = -1
     MeanMass = 0.5*( MeanMajorMass_1d(-1) + MeanMajorMass_1d(0))
     MeanGravity = 0.5*( EffectiveGravity(iAlt+1) + EffectiveGravity(iAlt) )
     MeanTemp = 0.5*( Temp(iAlt+1) + Temp(iAlt) )
     InvScaleHeight = -1.0* MeanMass*MeanGravity/&
                      (Boltzmanns_Constant*MeanTemp)
     HydroNS(iAlt,iSpecies) = &
               HydroNS(iAlt+1,iSpecies)*(Temp(iAlt+1)/Temp(iAlt))*&
               exp(dAlt_F(iAlt)*InvScaleHeight)
  enddo 
  do iAlt = -1, nAlts + 2
     do iSpecies =  1, nSpecies
       HydroPressureS(iAlt,iSpecies) = &
           HydroNS(iAlt,iSpecies)*Boltzmanns_Constant*Temp(iAlt)
       HydroRhoS(iAlt,iSpecies) = Mass(iSpecies)*HydroNS(iAlt,iSpecies)
     enddo 
       HydroPressure(iAlt) = HydroNT(iAlt)*Boltzmanns_Constant*Temp(iAlt)
            HydroRho(iAlt) = HydroNT(iAlt)*MeanMajorMass_1d(iAlt)
  enddo 
  do iSpecies = 1, nSpecies
     RhoS(-1:nAlts+2,iSpecies) =  &
        Mass(iSpecies)*NS(-1:nAlts+2,iSpecies)
     NewRhoS(-1:nAlts+2,iSpecies) = RhoS(-1:nAlts+2,iSpecies)
     MomentumS(-1:nAlts+2,iSpecies) =  &
        Mass(iSpecies)*NS(-1:nAlts+2,iSpecies)*&
         VertVel(-1:nAlts+2,iSpecies)
     PressureS(-1:nAlts+2,iSpecies) =  &
        NS(-1:nAlts+2,iSpecies)*Temp(-1:nAlts+2)*&
        Boltzmanns_Constant
  enddo 

  do iAlt = -1, nAlts + 2
   TotalEnergy(iAlt) = &
        Press(iAlt)/(Gamma_1d(iAlt) - 1.0) + &
        0.5*Rho(iAlt)*(Vel_GD(iAlt,iUp_)**2.0 + &
                       Vel_GD(iAlt,iNorth_)**2.0 + &
                       Vel_GD(iAlt,iEast_)**2.0) 
  enddo 
  do iDim = 1, 3
     Momentum(-1:nAlts+2,iDim) = Rho(-1:nAlts+2)*Vel_GD(-1:nAlts+2,iDim)
  enddo 

  DeviationRhoSRatio(-1:nAlts+2,1:nSpecies) = & 
            abs(RhoS(-1:nAlts+2,1:nSpecies) - &
           HydroRhoS(-1:nAlts+2,1:nSpecies))/&
                RhoS(-1:nAlts+2,1:nSpecies)

  DeviationRhoS(-1:nAlts+2,1:nSpecies) = &
           RhoS(-1:nAlts+2, 1:nSpecies) - &
      HydroRhoS(-1:nAlts+2,1:nSpecies)

  do iAlt = -1, nAlts + 2
    do iSpecies = 1, nSpecies
      if ( abs(DeviationRhoSRatio(iAlt,iSpecies)) .gt. 1.0) then
          SubtractHydrostatic(iAlt,iSpecies) = .false.
      else
          SubtractHydrostatic(iAlt,iSpecies) = .true.
      endif 
    enddo 
  enddo 

  NewRho = Rho
  NewPress = Press
  NewTotalEnergy = TotalEnergy
  ! Call the AUSM Solvers
  call calc_all_fluxes_hydro(DtIn, RhoS, PressureS, HydroPressureS, HydroRhoS, &
        HydroPressure,  HydroRho, AUSMRhoSFluxes,AUSMMomentumSFluxes, &
        AUSMTotalEnergyFluxes, AUSMMomentumFluxes, RadialDistance_C,  &
        SubtractHydrostatic)

  AmpSP = (1.0/(10.0*DtIn))
  kSP = nAltsSponge + 1

  do iAlt = 1,nAlts
     NewLogRho(iAlt) = NewLogRho(iAlt) - DtIn * &
          (DivVel(iAlt) + Vel_GD(iAlt,iUp_) * GradLogRho(iAlt) ) &
          + DtIn * DiffLogRho(iAlt)

     do iSpecies=1,nSpecies
        NewRhoS(iAlt,iSpecies) = RhoS(iAlt,iSpecies) &
               - DtIn*(AUSMRhoSFluxes(iAlt,iSpecies))
        NewLogNS(iAlt,iSpecies) = alog( NewRhoS(iAlt,iSpecies)/Mass(iSpecies) )
     enddo

     do iSpecies=1,nIonsAdvect
        if (UseImprovedIonAdvection) then
           NewLogINS(iAlt,iSpecies) = NewLogINS(iAlt,iSpecies) - DtIn * &
                (DiviVel(iAlt) * LogINS(iAlt,iSpecies) + &
                IVel(iAlt,iUp_) * GradLogINS(iAlt,iSpecies) ) &
                + DtIn * DiffLogINS(iAlt,iSpecies)
        else
           NewLogINS(iAlt,iSpecies) = NewLogINS(iAlt,iSpecies) - DtIn * &
                (IVel(iAlt,iUp_) * GradLogINS(iAlt,iSpecies) ) &
                + DtIn * DiffLogINS(iAlt,iSpecies)
        endif
     enddo

  enddo !iAlt = 1,nAlts


  NewNS  = 0.0
  NewNT  = 0.0
  NewRho = 0.0

  do iAlt = -1, nAlts+2
     do iSpecies = 1, nSpecies
         NewNS(iAlt,iSpecies) = NewRhoS(iAlt,iSpecies)/Mass(iSpecies)

         NewRho(iAlt) = NewRho(iAlt) + &
                  NewRhoS(iAlt,iSpecies)

         NewNT(iAlt) = NewNT(iAlt) + &
               NewNS(iAlt,iSpecies)
     enddo 
  enddo 


  if (iAlt >= (nAlts - nAltsSponge)) then
     NuSP = AmpSP*(1.0 - cos( pi*(kSP - (nAlts - iAlt))/kSP) )
  else
     NuSP = 0.0
  endif

  if (UseDamping) then
     VertTau(iAlt) = &
          15 - (1 - exp(-1.0*altitude_G(ialt)/1000.0/40.0))*5.0
  endif

  do iAlt = 1,nAlts
     do iSpecies=1,nSpecies

        if (.not. SubtractHydrostatic(iAlt,iSpecies)) then
           NewMomentumS(iAlt,iSpecies) = MomentumS(iAlt,iSpecies) - &
                 DtIn*(AUSMMomentumSFluxes(iAlt,iSpecies)) + &
                 DtIn*RhoS(iAlt,iSpecies)*EffectiveGravity(iAlt) 
        else
           NewMomentumS(iAlt,iSpecies) = MomentumS(iAlt,iSpecies) - &
                 DtIn*(AUSMMomentumSFluxes(iAlt,iSpecies)) + &
                 DtIn*DeviationRhoS(iAlt,iSpecies)*EffectiveGravity(iAlt) 
        endif 

!        NewMomentumS(iAlt,iSpecies) = NewMomentumS(iAlt,iSpecies) + &
!              DtIn*ChemSources_1d(iAlt,iSpecies)*&
!              VertVel(iAlt,iSpecies)*Mass(iSpecies)
!------------------
        NewVertVel(iAlt,iSpecies) = &
           NewMomentumS(iAlt,iSpecies)/NewRhoS(iAlt,iSpecies) 
        ! Add Explicit Spherical Curvature Terms here
        NewVertVel(iAlt,iSpecies) = &
           NewVertVel(iAlt,iSpecies) + DtIn*&
            (Vel_GD(iAlt,iNorth_)**2 + Vel_GD(iAlt,iEast_)**2) &
            * InvRadialDistance_C(iAlt) 
        !------
        if (UseCoriolis) then
           NewVertVel(iAlt,ispecies) = NewVertVel(iAlt,ispecies) + DtIn * ( &
                Coriolis * Vel_GD(iAlt,iEast_))
        endif
        ! Thermal Diffusion Effects (For Light Species H2, H, and He) 
        ! ThermalDiffCoefS is set in calc_rates
        ! Note:  ThermalDiffCoefS is < 0.0 for light species
        ! This forces light species into hot zones and heavy species into cold zones
        NewVertVel(iAlt,iSpecies) = NewVertVel(iAlt,iSpecies) - &
          DtIn*(ThermalDiffCoefS(iSpecies)*Boltzmanns_Constant*&
             GradTemp(iAlt))/Mass(iSpecies)
     enddo ! iSpecies

  enddo ! iAlt

  if (UseNeutralFriction) then
     nVel(1:nAlts,1:nSpecies) = NewVertVel(1:nAlts,1:nSpecies)
     NS_small = NewNS(1:nAlts,1:nSpecies)
     call calc_neutral_friction(DtIn,nVel(1:nAlts,1:nSpecies), &
                         EddyCoef_1d(1:nAlts), &
                               NewNT(1:nAlts), &
                               NS_small, &
                         GradLogConS(1:nAlts,1:nSpecies), &
                                Temp(1:nAlts))
     NewVertVel(1:nAlts,1:nSpecies) = nVel(1:nAlts,1:nSpecies)
  endif 


  do iAlt = -1, nAlts+2
     do iSpecies=1,nSpecies
        NewVertVel(iAlt, iSpecies) = max(-MaximumVerticalVelocity, &
             NewVertVel(iAlt, iSpecies))
        NewVertVel(iAlt, iSpecies) = min( MaximumVerticalVelocity, &
             NewVertVel(iAlt, iSpecies))
     enddo
  enddo

  NewVel_GD(-1:nAlts+2,iUp_) = 0.0
  do iAlt = -1, nAlts+2
     do iSpecies=1,nSpecies
!        NewVertVel(iAlt, iSpecies) = max(-MaximumVerticalVelocity, &
!             NewVertVel(iAlt, iSpecies))
!        NewVertVel(iAlt, iSpecies) = min( MaximumVerticalVelocity, &
!             NewVertVel(iAlt, iSpecies))
!
!        NewVertVel(iAlt, iSpecies) = max(-VelocityCap, &
!             NewVertVel(iAlt, iSpecies))
!        NewVertVel(iAlt, iSpecies) = min( VelocityCap, &
!             NewVertVel(iAlt, iSpecies))
!
        NewVel_GD(iAlt,iUp_) = NewVel_GD(iAlt,iUp_) + &
             NewVertVel(iAlt, iSpecies)*NewRhoS(iAlt,iSpecies)/NewRho(iAlt)
     enddo
  enddo


  do iAlt = 1, nAlts

      NewMomentum(iAlt,iEast_) = Momentum(iAlt,iEast_) - &
           DtIn*(AUSMMomentumFluxes(iAlt,iEast_)) 
 
      NewVel_GD(iAlt,iEast_) = NewMomentum(iAlt,iEast_)/NewRho(iAlt)
      NewVel_GD(iAlt,iEast_) = NewVel_GD(iAlt,iEast_) &
          + 0.35*exp(-1.0*(Altitude_G(iAlt) - Altitude_G(0))/&
            (45.0e+03)) * DtIn*DiffVel_CD(iAlt,iEast_)
 
      NewMomentum(iAlt,iNorth_) = Momentum(iAlt,iNorth_) - &
           DtIn*(AUSMMomentumFluxes(iAlt,iNorth_)) 
 
      NewVel_GD(iAlt,iNorth_) = NewMomentum(iAlt,iNorth_)/NewRho(iAlt)
      NewVel_GD(iAlt,iNorth_) = NewVel_GD(iAlt,iNorth_) &
          + 0.35*exp(-1.0*(Altitude_G(iAlt) - Altitude_G(0))/&
            (45.0e+03)) * DtIn*DiffVel_CD(iAlt,iNorth_)

      ! dVphi/dt = - (V grad V)_phi
!       NewVel_GD(iAlt,iEast_) = Vel_GD(iAlt,iEast_) - DtIn * &
!            Vel_GD(iAlt,iUp_)*GradVel_CD(iAlt,iEast_) &
!            + DtIn * DiffVel_CD(iAlt,iEast_)
! 
!       ! dVtheta/dt = - (V grad V)_theta
!       NewVel_GD(iAlt,iNorth_) = Vel_GD(iAlt,iNorth_) - DtIn * &
!            Vel_GD(iAlt,iUp_)*GradVel_CD(iAlt,iNorth_) &
!            + DtIn * DiffVel_CD(iAlt,iNorth_)
 
     ! dT/dt = -(V.grad T + (gamma - 1) T div V +  &
     !        (gamma - 1) * g  * grad (KeH^2  * rho) /rho 
     ! AUSM Method
     NewTotalEnergy(iAlt)   = TotalEnergy(iAlt) - &
         DtIn*AUSMTotalEnergyFluxes(iAlt) + &
         DtIn*Rho(iAlt)*Vel_GD(iAlt,iUp_)*EffectiveGravity(iAlt) 
 
     NewPress(iAlt) = &
        (NewTotalEnergy(iAlt) - &
         0.5*NewRho(iAlt)*(NewVel_GD(iAlt,iUp_)**2.0 + &
                           NewVel_GD(iAlt,iNorth_)**2.0 + &
                           NewVel_GD(iAlt,iEast_ )**2.0 ) )*&
                           (Gamma_1d(iAlt) - 1.0)
 
     NewTemp(iAlt) = NewPress(iAlt)/(Boltzmanns_Constant*NewNT(iAlt))

!     NewTemp(iAlt) = NewTemp(iAlt) + &
!       0.25*exp(-1.0*(Altitude_G(iAlt) - Altitude_G(0))/&
!       (45.0e+03)) * DtIn*DiffTemp(iAlt)

  enddo ! iAlt

  do iAlt = 1, nAlts
     NewSumRho    = sum( Mass(1:nSpecies)*exp(NewLogNS(iAlt,1:nSpecies)) )
     NewLogRho(iAlt) = alog(NewSumRho)
  enddo

end subroutine advance_vertical_1stage_ausm

!\
! ------------------------------------------------------------
! calc_rusanov
! ------------------------------------------------------------
!/

subroutine calc_rusanov_alts_ausm(Var, GradVar, DiffVar)

  use ModSizeGitm
  use ModVertical, only : dAlt_C, cMax
  implicit none

  real, intent(in) :: Var(-1:nAlts+2)
  real, intent(out):: GradVar(1:nAlts), DiffVar(1:nAlts)

  real, dimension(1:nAlts+1) :: VarLeft, VarRight, DiffFlux
  !------------------------------------------------------------

  call calc_facevalues_alts_ausm(Var, VarLeft, VarRight)

  ! Gradient based on averaged Left/Right values
  GradVar = 0.5 * &
       (VarLeft(2:nAlts+1)+VarRight(2:nAlts+1) - &
       VarLeft(1:nAlts)-VarRight(1:nAlts))/dAlt_C(1:nAlts)

  ! Rusanov/Lax-Friedrichs diffusive term
  DiffFlux = 0.5 * max(cMax(0:nAlts),cMax(1:nAlts+1)) * (VarRight - VarLeft)

  DiffVar = (DiffFlux(2:nAlts+1) - DiffFlux(1:nAlts))/dAlt_C(1:nAlts)

end subroutine calc_rusanov_alts_ausm

!\
! ------------------------------------------------------------
! calc_facevalues_alts_ausm
! ------------------------------------------------------------
!/

subroutine calc_facevalues_alts_ausm(Var, VarLeft, VarRight)

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

end subroutine calc_facevalues_alts_ausm


subroutine calc_all_fluxes_hydro(DtIn, RhoS, PressureS, HydroPressureS, HydroRhoS, &
            HydroPressure, HydroRho, RhoSFlux, MomentumSFlux, &
            EnergyFlux, MomentumFlux, &
            RadDist, SubtractHydrostatic)

  use ModSizeGitm
  use ModVertical, only : dAlt_C, cMax, VertVel, Gamma_1D, &
                          LogRho, Vel_GD, MeanMajorMass_1d, &
                          Temp, Altitude_G


  use ModPlanet, only : nSpecies, nIonsAdvect, Mass, RBody, &
                        iN2_, cSpecies
  use ModGITM, only : Dt,iUp_, iEast_, iNorth_
  use ModConstants, only : Boltzmanns_Constant

  implicit none

! Passed from Vertical_Solver In
  real, intent(in) :: DtIn
  real, intent(in) :: RhoS(-1:nAlts+2, 1:nSpecies)
  real, intent(in) :: PressureS(-1:nAlts+2,1:nSpecies)
  real, intent(in) :: HydroRhoS(-1:nAlts+2, 1:nSpecies)
  real, intent(in) :: HydroPressureS(-1:nAlts+2,1:nSpecies)
  real, intent(in) :: HydroRho(-1:nAlts+2) 
  real, intent(in) :: HydroPressure(-1:nAlts+2)
  real, intent(out):: RhoSFlux(1:nAlts,1:nSpecies)
  real, intent(out):: MomentumSFlux(1:nAlts,1:nSpecies)
  real, intent(out):: EnergyFlux(1:nAlts)
  real, intent(out):: MomentumFlux(1:nAlts,1:3)
  real, intent(in) :: RadDist(-1:nAlts+2)
  logical, intent(in) :: SubtractHydrostatic(-1:nAlts+2,1:nSpecies)

  real, dimension(1:nAlts,1:nSpecies) :: RhoSLeft_M12, RhoSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: RhoSLeft_P12, RhoSRight_P12

  real, dimension(-1:nAlts+2, 1:nSpecies) :: LogRhoS
  real, dimension(-1:nAlts+2, 1:nSpecies) :: LogPS

  real, dimension(1:nAlts,1:nSpecies) :: LogRhoSLeft_M12, LogRhoSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: LogRhoSLeft_P12, LogRhoSRight_P12

  real, dimension( 1:nAlts  ,1:nIonsAdvect) :: RhoILeft_M12, RhoIRight_M12
  real, dimension( 1:nAlts  ,1:nIonsAdvect) :: RhoILeft_P12, RhoIRight_P12
  real, dimension(-1:nAlts+2,1:nIonsAdvect) :: LogRhoI
  real, dimension( 1:nAlts  ,1:nIonsAdvect) :: LogRhoILeft_M12, LogRhoIRight_M12
  real, dimension( 1:nAlts  ,1:nIonsAdvect) :: LogRhoILeft_P12, LogRhoIRight_P12


  real, dimension(1:nAlts,1:nSpecies) :: PressureSLeft_M12, PressureSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: PressureSLeft_P12, PressureSRight_P12

  real, dimension(1:nAlts,1:nSpecies) :: LogPressureSLeft_M12, LogPressureSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: LogPressureSLeft_P12, LogPressureSRight_P12

  real, dimension(1:nAlts,1:nSpecies) :: VelLeft_M12, VelRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: VelLeft_P12, VelRight_P12

  real, dimension(1:nAlts,3) :: VelGDLeft_M12, VelGDRight_M12
  real, dimension(1:nAlts,3) :: VelGDLeft_P12, VelGDRight_P12

  real, dimension(1:nAlts,3) :: IVelLeft_M12, IVelRight_M12
  real, dimension(1:nAlts,3) :: IVelLeft_P12, IVelRight_P12

  real, dimension(1:nAlts) :: PLeft_M12, PRight_M12
  real, dimension(1:nAlts) :: PLeft_P12, PRight_P12

  real, dimension(1:nAlts) :: GammaLeft_M12, GammaRight_M12
  real, dimension(1:nAlts) :: GammaLeft_P12, GammaRight_P12

  real, dimension(1:nAlts) :: ELeft_M12, ERight_M12
  real, dimension(1:nAlts) :: ELeft_P12, ERight_P12

  real, dimension(1:nAlts) :: RhoLeft_M12, RhoRight_M12
  real, dimension(1:nAlts) :: RhoLeft_P12, RhoRight_P12

  real, dimension(1:nAlts) :: CSLeft_M12, CSRight_M12
  real, dimension(1:nAlts) :: CSLeft_P12, CSRight_P12

  real, dimension(1:nAlts,1:nSpecies) :: RhoSFlux_M12, RhoSFlux_P12
  real, dimension(1:nAlts,1:nSpecies) :: MomentumSFlux_M12, MomentumSFlux_P12
  real, dimension(1:nAlts,3) :: Momentum_M12, Momentum_P12
  real, dimension(1:nAlts,1:nIonsAdvect) :: RhoIFlux_M12, RhoIFlux_P12

  real, dimension(1:nAlts) :: EnergyFlux_M12, EnergyFlux_P12

!!! Hydrostatic Variables
  real, dimension(-1:nAlts+2,1:nSpecies) :: LogHydroPressureS
  real, dimension(1:nAlts,1:nSpecies) :: LogHydroPressureSLeft_M12, &
                                         LogHydroPressureSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: LogHydroPressureSLeft_P12, &
                                         LogHydroPressureSRight_P12

  real, dimension(1:nAlts,1:nSpecies) :: HydroPressureSLeft_M12, &
                                        HydroPressureSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: HydroPressureSLeft_P12, &
                                         HydroPressureSRight_P12

  real, dimension(1:nAlts,1:nSpecies) :: MeanHydroPressureS_M12
  real, dimension(1:nAlts,1:nSpecies) :: MeanHydroPressureS_P12

  real :: SubCs
  integer :: iSpecies, iAlt, iDim
  !------------------------------------------------------------

  ! ==================== AUSM Flux Variables

  real, dimension( 1:nAlts,1:nSpecies) :: NumericalVelocity_P12, &
                                          NumericalVelocity_M12
  real, dimension( 1:nAlts,1:nSpecies) :: NumericalPressure_P12, &
                                          NumericalPressure_M12   

  real, dimension( 1:nAlts) :: BulkNumericalVelocity_P12, &
                               BulkNumericalVelocity_M12    
  real, dimension( 1:nAlts) :: BulkIVel_P12, BulkIVel_M12    
  real, dimension( 1:nAlts) :: MeanCS_P12, MeanCS_M12    

  real, dimension( 1:nAlts, 1:nSpecies) :: MeanPressureS_P12,&
                                           MeanPressureS_M12    
  real, dimension( 1:nAlts, 1:nSpecies) :: MeanRhoS_P12, &
                                           MeanRhoS_M12    

  real :: Kp(1:nSpecies), Ku(1:nAlts,1:nSpecies)
  real :: LiouKp, LiouKu

  real :: LiouKpS(1:nAlts,1:nSpecies), LiouKuS(1:nAlts,1:nSpecies)
  real :: MaxKpS(1:nAlts,1:nSpecies), MaxKuS(1:nAlts,1:nSpecies)
  real :: MinKpS(1:nAlts,1:nSpecies), MinKuS(1:nAlts,1:nSpecies)
  real :: KpWidth, KpAltMidPoint
  integer :: AltIndex

  real, dimension( 1:nAlts) :: LeftRadius, RightRadius    
  real, dimension( 1:nAlts) :: AreaFunction_P12, AreaFunction_M12    
  real, dimension( 1:nAlts) :: LocalCellVolume

!!!!!! ====================  JMB:  Liou Stuff ===============================
  real, dimension(1:nAlts) :: LiouCSLeft_M12, LiouCSRight_M12
  real, dimension(1:nAlts) :: LiouCSLeft_P12, LiouCSRight_P12

  real, dimension(1:nAlts) :: LiouEnthalpyLeft_M12, LiouEnthalpyRight_M12
  real, dimension(1:nAlts) :: LiouEnthalpyLeft_P12, LiouEnthalpyRight_P12

  real, dimension(1:nAlts) :: InterfaceCS_M12, InterfaceCS_P12

  real, dimension(1:nAlts,1:nSpecies) :: MLeft_M12, MRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: MLeft_P12, MRight_P12

  real, dimension(1:nAlts,1:nSpecies) :: M2Bar_M12, M2Bar_P12
  real, dimension(1:nAlts,1:nSpecies) :: M2Zero_M12, M2Zero_P12

  real, dimension(1:nAlts,1:nSpecies) :: MZero_M12, MZero_P12
  real:: MInf, LiouBeta
  real, dimension(1:nAlts,1:nSpecies):: ModifiedZeta

  real, dimension(1:nAlts,1:nSpecies) :: FA_M12, FA_P12

!!! Polynomial Mach Functions

!!! First Order Polynomial
  real, dimension(1:nAlts,1:nSpecies) :: MF1P_Left_M12, MF1N_Left_M12
  real, dimension(1:nAlts,1:nSpecies) :: MF1P_Right_M12, MF1N_Right_M12

  real, dimension(1:nAlts,1:nSpecies) :: MF1P_Left_P12, MF1N_Left_P12
  real, dimension(1:nAlts,1:nSpecies) :: MF1P_Right_P12, MF1N_Right_P12

!!!! 2nd Order Polynomial
  real, dimension(1:nAlts,1:nSpecies) :: MF2P_Left_M12, MF2N_Left_M12
  real, dimension(1:nAlts,1:nSpecies) :: MF2P_Right_M12, MF2N_Right_M12

  real, dimension(1:nAlts,1:nSpecies) :: MF2P_Left_P12, MF2N_Left_P12
  real, dimension(1:nAlts,1:nSpecies) :: MF2P_Right_P12, MF2N_Right_P12

!!!! 4th Order Polynomial
  real, dimension(1:nAlts,1:nSpecies) :: MF4P_Left_M12, MF4N_Left_M12
  real, dimension(1:nAlts,1:nSpecies) :: MF4P_Right_M12, MF4N_Right_M12

  real, dimension(1:nAlts,1:nSpecies) :: MF4P_Left_P12, MF4N_Left_P12
  real, dimension(1:nAlts,1:nSpecies) :: MF4P_Right_P12, MF4N_Right_P12

!!!! Pressure Mach Number
  real, dimension(1:nAlts,1:nSpecies) :: MPress_M12, MPress_P12

!!!! Interface Mach Number
  real, dimension(1:nAlts,1:nSpecies) :: InterfaceMach_M12, InterfaceMach_P12

  real, dimension(1:nAlts,1:nSpecies) :: LiouNumericalVelocity_M12, &
                                         LiouNumericalVelocity_P12

!!! Extract Mean Atmosphere Values

  real :: LogKpS(1:nAlts,1:nSpecies), LogKuS(1:nAlts,1:nSpecies)
  real :: LogMinKpS(1:nAlts,1:nSpecies), LogMinKuS(1:nAlts,1:nSpecies)
  real :: LogMaxKpS(1:nAlts,1:nSpecies), LogMaxKuS(1:nAlts,1:nSpecies)

!  write(*,*) 'Key Variables =================='
!  do iAlt = -1, nAlts + 2
!       write(*,*) 'iAlt, RadialDist, Gamma_1d =', iAlt, &
!             RadDist(iAlt), Gamma_1d(iAlt)
!  enddo 

!  do iSpecies = 1, nSpecies
!     write(*,*) 'SPECIES ==================', cSpecies(iSpecies)
!     do iAlt = -1, nAlts + 2
!          write(*,*) 'iAlt, Subtract Hydro =', iAlt, &
!                SubtractHydrostatic(iAlt,iSpecies)
!     enddo 
!  enddo 
!
!  do iSpecies = 1, nSpecies
!     write(*,*) 'SPECIES ==================', cSpecies(iSpecies)
!     do iAlt = -1, nAlts + 2
!          write(*,*) 'iAlt, RhoS, VertVelS =', iAlt, &
!                RhoS(iAlt,iSpecies), VertVel(iAlt,iSpecies)
!     enddo 
!  enddo 
!
!
!stop

  MInf = 1.0e-19
  LiouBeta = 1.0/8.0

  LogRhoS(-1:nAlts+2,1:nSpecies) = alog(RhoS(-1:nAlts+2,1:nSpecies))
  LogPS(-1:nAlts+2,1:nSpecies) = alog(PressureS(-1:nAlts+2,1:nSpecies))
  LogHydroPressureS(-1:nAlts+2,1:nSpecies) = &
                            alog(HydroPressureS(-1:nAlts+2,1:nSpecies))

  Kp(1:nSpecies) = 0.10             !! Ullrich et al. [2011]
  Ku(1:nAlts,1:nSpecies) = 2.0             !! Ullrich et al. [2011]

!!! 1-D Settings
  MinKuS(1:nAlts,1:nSpecies) = 0.25
  MaxKuS(1:nAlts,1:nSpecies) = 0.50

  MinKpS(1:nAlts,1:nSpecies) = 1.0e-19
  MaxKpS(1:nAlts,1:nSpecies) = 1.0e-09

  LogMaxKuS = alog(MaxKuS)
  LogMinKuS = alog(MinKuS)

  LogMaxKpS = alog(MaxKpS)
  LogMinKpS = alog(MinKpS)

   KpWidth = 15.0e+03
   KpAltMidPoint = 160.0e+03
 
  do iAlt = 1, nAlts
     do iSpecies = 1, nSpecies
         LogKpS(iAlt,iSpecies) = LogMinKpS(iAlt,iSpecies) +  &
               0.5*(LogMaxKpS(iAlt,iSpecies) - LogMinKpS(iAlt,iSpecies))*&
               ( 1.0 + tanh(  (Altitude_G(iAlt) - KpAltMidPoint)/KpWidth ) )

         LogKuS(iAlt,iSpecies) = LogMinKuS(iAlt,iSpecies) +  &
               0.5*(LogMaxKuS(iAlt,iSpecies) - LogMinKuS(iAlt,iSpecies))*&
               ( 1.0 + tanh(  (Altitude_G(iAlt) - KpAltMidPoint)/KpWidth ) )

        LiouKpS(iAlt,iSpecies) = exp(LogKpS(iAlt,iSpecies))
             Ku(iAlt,iSpecies) = exp(LogKuS(iAlt,iSpecies))
     enddo 
  enddo 

!!!! Grab the left and right states of the Variables
!!!!  on boh Interfaces (P12 = +1/2 and M12 = -1/2)
    do iSpecies = 1, nSpecies
          !! Calculate the Left and Right Faces of the RhoS 
           call calc_kt_facevalues(LogRhoS(-1:nAlts+2,iSpecies), &
                           LogRhoSLeft_M12( 1:nAlts  ,iSpecies), &
                          LogRhoSRight_M12( 1:nAlts  ,iSpecies), &
                           LogRhoSLeft_P12( 1:nAlts  ,iSpecies), &
                          LogRhoSRight_P12( 1:nAlts  ,iSpecies) )

           RhoSLeft_M12(:,iSpecies) = exp( LogRhoSLeft_M12(:,iSpecies)) 
          RhoSRight_M12(:,iSpecies) = exp(LogRhoSRight_M12(:,iSpecies)) 

           RhoSLeft_P12(:,iSpecies) = exp( LogRhoSLeft_P12(:,iSpecies)) 
          RhoSRight_P12(:,iSpecies) = exp(LogRhoSRight_P12(:,iSpecies)) 

    enddo 

!    write(*,*) '==================== CALC_HYDRO_FLUXES ================'
!    write(*,*) '======= M12 FACEVALUES ========'
!    do iSpecies = 1, nSpecies
!       write(*,*) '  SPECIES TYPE :    ',iSpecies, cSpecies(iSpecies)
!    do iAlt = 1, nAlts
!        write(*,*) 'iAlt, RhoS(i-1), Left_M12, Right_M12, RhoS(i) = ', &
!                 iAlt, RhoS(iAlt-1,iSpecies), RhoSLeft_M12(iAlt,iSpecies), &
!          RhoSRight_M12(iAlt,iSpecies), RhoS(iAlt,iSpecies)  
!    enddo 
!    enddo 
!
!    write(*,*) '======= P12 FACEVALUES ========'
!    do iSpecies = 1, nSpecies
!       write(*,*) '  SPECIES TYPE :    ',iSpecies, cSpecies(iSpecies)
!    do iAlt = 1, nAlts
!        write(*,*) 'iAlt, RhoS(i-1), Left_P12, Right_P12, RhoS(i) = ', &
!                 iAlt, RhoS(iAlt,iSpecies), RhoSLeft_P12(iAlt,iSpecies), &
!          RhoSRight_P12(iAlt,iSpecies), RhoS(iAlt+1,iSpecies)  
!    enddo 
!    enddo 
!
!stop


    do iSpecies = 1, nSpecies
          !! Calculate the Left and Right Faces of the RhoS 
           call calc_kt_facevalues(LogHydroPressureS(-1:nAlts+2,iSpecies), &
                              LogHydroPressureSLeft_M12(1:nAlts,iSpecies), &
                              LogHydroPressureSRight_M12(1:nAlts,iSpecies), &
                               LogHydroPressureSLeft_P12(1:nAlts,iSpecies), &
                              LogHydroPressureSRight_P12(1:nAlts,iSpecies) )

           HydroPressureSLeft_M12(:,iSpecies) = &
                 exp( LogHydroPressureSLeft_M12(:,iSpecies)) 
          HydroPressureSRight_M12(:,iSpecies) = &
                 exp(LogHydroPressureSRight_M12(:,iSpecies)) 

           HydroPressureSLeft_P12(:,iSpecies) = &
                 exp( LogHydroPressureSLeft_P12(:,iSpecies)) 
          HydroPressureSRight_P12(:,iSpecies) = &
                 exp(LogHydroPressureSRight_P12(:,iSpecies)) 

    enddo 


   do iSpecies = 1, nSpecies
         !! Calculate the Left and Right Faces of the PressureS
      call calc_kt_facevalues(LogPS(:,iSpecies), &
          LogPressureSLeft_M12(:,iSpecies), LogPressureSRight_M12(:,iSpecies), &
          LogPressureSLeft_P12(:,iSpecies), LogPressureSRight_P12(:,iSpecies) )

      PressureSLeft_M12(:,iSpecies) = exp( LogPressureSLeft_M12(:,iSpecies)) 
      PressureSRight_M12(:,iSpecies) = exp(LogPressureSRight_M12(:,iSpecies)) 

      PressureSLeft_P12(:,iSpecies) = exp( LogPressureSLeft_P12(:,iSpecies)) 
      PressureSRight_P12(:,iSpecies) = exp(LogPressureSRight_P12(:,iSpecies)) 

   enddo 

   do iSpecies = 1, nSpecies
         !! Calculate the Left and Right Faces of the Var (Rho) 
     call calc_kt_facevalues(VertVel(:,iSpecies), &
                         VelLeft_M12(:,iSpecies), VelRight_M12(:,iSpecies), &
                         VelLeft_P12(:,iSpecies), VelRight_P12(:,iSpecies) )
   enddo 

!    write(*,*) '======= M12 VELOCITY FACEVALUES ========'
!    do iSpecies = 1, nSpecies
!       write(*,*) '  SPECIES TYPE :    ',iSpecies, cSpecies(iSpecies)
!    do iAlt = 1, nAlts
!        write(*,*) 'iAlt, Vel(i-1), VelLeft_M12, VelRight_M12, Vel(i) = ', &
!                 iAlt, VertVel(iAlt-1,iSpecies), VelLeft_M12(iAlt,iSpecies), &
!          VelRight_M12(iAlt,iSpecies), VertVel(iAlt,iSpecies)  
!    enddo 
!    enddo 
!
!    write(*,*) '======= P12 VELOCITY FACEVALUES ========'
!    do iSpecies = 1, nSpecies
!       write(*,*) '  SPECIES TYPE :    ',iSpecies, cSpecies(iSpecies)
!    do iAlt = 1, nAlts
!        write(*,*) 'iAlt, Vel(i-1), VelLeft_P12, VelRight_P12, Vel(i) = ', &
!                 iAlt, VertVel(iAlt-1,iSpecies), VelLeft_P12(iAlt,iSpecies), &
!          VelRight_P12(iAlt,iSpecies), VertVel(iAlt,iSpecies)  
!    enddo 
!    enddo 


   RhoLeft_M12(:) = 0.0
   RhoRight_M12(:) = 0.0

   RhoLeft_P12(:) = 0.0
   RhoRight_P12(:) = 0.0

   do iSpecies = 1, nSpecies
      RhoLeft_M12(1:nAlts) = RhoLeft_M12(1:nAlts) + &
                            RhoSLeft_M12(1:nAlts,iSpecies)
      RhoRight_M12(1:nAlts) = RhoRight_M12(1:nAlts) + &
                             RhoSRight_M12(1:nAlts,iSpecies)

      RhoLeft_P12(1:nAlts) = RhoLeft_P12(1:nAlts) + &
                            RhoSLeft_P12(1:nAlts,iSpecies)
      RhoRight_P12(1:nAlts) = RhoRight_P12(1:nAlts) + &
                             RhoSRight_P12(1:nAlts,iSpecies)
   enddo 

   PLeft_M12(:) = 0.0
   PRight_M12(:) = 0.0
 
   PLeft_P12(:) = 0.0
   PRight_P12(:) = 0.0

   do iSpecies = 1, nSpecies
      PLeft_M12(1:nAlts) = PLeft_M12(1:nAlts) + &
                   PressureSLeft_M12(1:nAlts,iSpecies)
      PRight_M12(1:nAlts) = PRight_M12(1:nAlts) + &
                    PressureSRight_M12(1:nAlts,iSpecies)

      PLeft_P12(1:nAlts) = PLeft_P12(1:nAlts) + &
                   PressureSLeft_P12(1:nAlts,iSpecies)
      PRight_P12(1:nAlts) = PRight_P12(1:nAlts) + &
                    PressureSRight_P12(1:nAlts,iSpecies)
   enddo 


   do iDim = 1, 3
      !! Calculate the Left and Right Faces of the Var (Rho) 
       call calc_kt_facevalues(Vel_GD(:,iDim), &
                        VelGDLeft_M12(:,iDim), &
                       VelGDRight_M12(:,iDim), &
                     VelGDLeft_P12(:,iDim), &
                       VelGDRight_P12(:,iDim) )
   enddo 

   VelGDLeft_M12(:,iUp_) = 0.0
   VelGDRight_M12(:,iUp_) = 0.0

   VelGDLeft_P12(:,iUp_) = 0.0
   VelGDRight_P12(:,iUp_) = 0.0

   do iSpecies = 1, nSpecies

        VelGDLeft_M12(1:nAlts,iUp_) = &
                 VelGDLeft_M12(1:nAlts,iUp_) + &
                 RhoSLeft_M12(1:nAlts,iSpecies)*&
                 VelLeft_M12(1:nAlts,iSpecies)/RhoLeft_M12(1:nAlts)

        VelGDRight_M12(1:nAlts,iUp_) = &
                 VelGDRight_M12(1:nAlts,iUp_) + &
                 RhoSRight_M12(1:nAlts,iSpecies)*&
                 VelRight_M12(1:nAlts,iSpecies)/RhoRight_M12(1:nAlts)

        VelGDLeft_P12(1:nAlts,iUp_) = &
                 VelGDLeft_P12(1:nAlts,iUp_) + &
                 RhoSLeft_P12(1:nAlts,iSpecies)*&
                 VelLeft_P12(1:nAlts,iSpecies)/RhoLeft_P12(1:nAlts)

        VelGDRight_P12(1:nAlts,iUp_) = &
                 VelGDRight_P12(1:nAlts,iUp_) + &
                 RhoSRight_P12(1:nAlts,iSpecies)*&
                 VelRight_P12(1:nAlts,iSpecies)/RhoRight_P12(1:nAlts)
   enddo 

        !! Calculate the Left and Right Faces of the Pressure 
   call calc_kt_facevalues(Gamma_1d(:), GammaLeft_M12(:), GammaRight_M12(:), &
                                        GammaLeft_P12(:), GammaRight_P12(:) )

    do iAlt = 1, nAlts
!!!! ============= Bulk Values
        ELeft_M12(iAlt) = &
            ( 1.0/(GammaLeft_M12(iAlt) - 1.0))*PLeft_M12(iAlt) + &
              0.5*RhoLeft_M12(iAlt)*&
             (VelGDLeft_M12(iAlt,iUp_)**2.0 + VelGDLeft_M12(iAlt,iEast_)**2.0 + &
              VelGDLeft_M12(iAlt,iNorth_)**2.0)

        ERight_M12(iAlt) = &
            ( 1.0/(GammaRight_M12(iAlt) - 1.0))*PRight_M12(iAlt) + &
              0.5*RhoRight_M12(iAlt)*&
             (VelGDRight_M12(iAlt,iUp_)**2.0 + VelGDRight_M12(iAlt,iEast_)**2.0 + &
              VelGDRight_M12(iAlt,iNorth_)**2.0)

        ELeft_P12(iAlt) = &
            ( 1.0/(GammaLeft_P12(iAlt) - 1.0))*PLeft_P12(iAlt) + &
              0.5*RhoLeft_P12(iAlt)* &
             (VelGDLeft_P12(iAlt,iUp_)**2.0 + VelGDLeft_P12(iAlt,iEast_)**2.0 + &
              VelGDLeft_P12(iAlt,iNorth_)**2.0)

        ERight_P12(iAlt) = &
            ( 1.0/(GammaRight_P12(iAlt) - 1.0))*PRight_P12(iAlt) + &
              0.5*RhoRight_P12(iAlt)* &
             (VelGDRight_P12(iAlt,iUp_)**2.0 + VelGDRight_P12(iAlt,iEast_)**2.0 + &
              VelGDRight_P12(iAlt,iNorth_)**2.0)

    enddo 



!!!! Liou et al. [2006] suggest using the Enthalpy for the 
!!!! numerical speed of sound.
!!!! Calculate the Enthalpy at the cell faces here.
!    write(*,*) 'ENTHALPY CALC:=============='
    do iAlt = 1, nAlts 

!      write(*,*) 'Gamma(i-1), GammaLeft_M12, GammaRight_M12, Gamma(i) =', &
!            Gamma_1d(iAlt-1), GammaLeft_M12(iAlt), GammaRight_M12(iAlt), &
!            Gamma_1d(iAlt)

       LiouEnthalpyLeft_M12(iAlt) = &
         0.5*(VelGDLeft_M12(iAlt,iUp_)**2.0 + &
              VelGDLeft_M12(iAlt,iEast_)**2.0 + &
              VelGDLeft_M12(iAlt,iNorth_)**2.0) + &
            (GammaLeft_M12(iAlt)/(GammaLeft_M12(iAlt) - 1.0))*&
            PLeft_M12(iAlt)/RhoLeft_M12(iAlt)

       LiouEnthalpyRight_M12(iAlt) = &
            0.5*(VelGDRight_M12(iAlt,iUp_)**2.0 + &
                 VelGDRight_M12(iAlt,iEast_)**2.0 + &
                 VelGDRight_M12(iAlt,iNorth_)**2.0) + &
            (GammaRight_M12(iAlt)/(GammaRight_M12(iAlt) - 1.0))*&
            PRight_M12(iAlt)/RhoRight_M12(iAlt)


       LiouEnthalpyLeft_P12(iAlt) = &
            0.5*(VelGDLeft_P12(iAlt,iUp_)**2.0 + &
                 VelGDLeft_P12(iAlt,iEast_)**2.0 + &
                 VelGDLeft_P12(iAlt,iNorth_)**2.0) + &
            (GammaLeft_P12(iAlt)/(GammaLeft_P12(iAlt) - 1.0))*&
            PLeft_P12(iAlt)/RhoLeft_P12(iAlt)


       LiouEnthalpyRight_P12(iAlt) = &
            0.5*(VelGDRight_P12(iAlt,iUp_)**2.0 + &
                 VelGDRight_P12(iAlt,iEast_)**2.0 + &
                 VelGDRight_P12(iAlt,iNorth_)**2.0) + &
            (GammaRight_P12(iAlt)/(GammaRight_P12(iAlt) - 1.0))*&
            PRight_P12(iAlt)/RhoRight_P12(iAlt)

!      write(*,*) 'Gamma(i), GammaLeft_P12, GammaRight_P12, Gamma(i+1) =', &
!            Gamma_1d(iAlt), GammaLeft_P12(iAlt), GammaRight_P12(iAlt), &
!            Gamma_1d(iAlt+1)

!      write(*,*) ' PLeft_M12(iAlt), PRight_M12, PLeft_P12, PRight_P12 =', &
!            PLeft_M12(iAlt), PRight_M12(iAlt), PLeft_P12(iAlt), &
!            PRight_P12(iAlt)
    enddo 

!!! Liou Methodology

   do iAlt = 1, nAlts

       SubCs = & 
         sqrt(2.0*( (GammaLeft_M12(iAlt) - 1.0 )/(GammaLeft_M12(iAlt) + 1.0)) *&
         LiouEnthalpyLeft_M12(iAlt) )

         LiouCSLeft_M12(iAlt) = & 
             (SubCs**2.0)/max(SubCs, VelGDLeft_M12(iAlt,iUp_))

       SubCs = &
       sqrt(2.0*( (GammaRight_M12(iAlt) - 1.0 )/(GammaRight_M12(iAlt) + 1.0)) *&
        LiouEnthalpyRight_M12(iAlt) )

       LiouCSRight_M12(iAlt) = &
           (SubCs**2.0)/max(SubCs, -1.0*VelGDRight_M12(iAlt,iUp_))

       InterfaceCS_M12(iAlt) = min(LiouCSLeft_M12(iAlt), LiouCSRight_M12(iAlt))

      SubCs = &
         sqrt(2.0*( (GammaLeft_P12(iAlt) - 1.0 )/(GammaLeft_P12(iAlt) + 1.0)) *&
            LiouEnthalpyLeft_P12(iAlt) )
      LiouCSLeft_P12(iAlt) = &
            (SubCs**2.0)/max(SubCs, VelGDLeft_P12(iAlt,iUp_))

      SubCs = &
      sqrt(2.0*( (GammaRight_P12(iAlt) - 1.0 )/(GammaRight_P12(iAlt) + 1.0)) *&
         LiouEnthalpyRight_P12(iAlt) )

         LiouCSRight_P12(iAlt) = &
              (SubCs**2.0)/max(SubCs, -1.0*VelGDRight_P12(iAlt,iUp_))

      InterfaceCS_P12(iAlt) = min(LiouCSLeft_P12(iAlt), LiouCSRight_P12(iAlt))


!      write(*,*) 'iAlt, InterfaceCS_M12(iAlt), InterfaceCS_P12(iAlt) =',&
!                  iAlt, InterfaceCS_M12(iAlt), InterfaceCS_P12(iAlt)
   enddo 


!stop
  
  MeanPressureS_P12(1:nAlts,1:nSpecies) = &
        0.5*(PressureSLeft_P12(1:nAlts,1:nSpecies) + &
            PressureSRight_P12(1:nAlts,1:nSpecies))

  MeanPressureS_M12(1:nAlts,1:nSpecies) = &
        0.5*(PressureSLeft_M12(1:nAlts,1:nSpecies) + &
            PressureSRight_M12(1:nAlts,1:nSpecies))

  MeanHydroPressureS_P12(1:nAlts,1:nSpecies) = &
        0.5*(HydroPressureSLeft_P12(1:nAlts,1:nSpecies) + &
            HydroPressureSRight_P12(1:nAlts,1:nSpecies))

  MeanHydroPressureS_M12(1:nAlts,1:nSpecies) = &
        0.5*(HydroPressureSLeft_M12(1:nAlts,1:nSpecies) + &
            HydroPressureSRight_M12(1:nAlts,1:nSpecies))

  MeanRhoS_P12(1:nAlts,1:nSpecies) = &
       0.5*(RhoSLeft_P12(1:nAlts,1:nSpecies) + RhoSRight_P12(1:nAlts,1:nSpecies))

  MeanRhoS_M12(1:nAlts,1:nSpecies) = &
       0.5*(RhoSLeft_M12(1:nAlts,1:nSpecies) + RhoSRight_M12(1:nAlts,1:nSpecies))

  do iAlt = 1, nAlts 
     MeanCS_P12(iAlt) = InterfaceCS_P12(iAlt)
     MeanCS_M12(iAlt) = InterfaceCS_M12(iAlt)
  enddo 

!!! Next, define local mach numbers at the interfaces

   do iAlt = 1, nAlts 
     do iSpecies = 1, nSpecies

        MLeft_M12(iAlt,iSpecies) = &
           VelLeft_M12(iAlt,iSpecies)/MeanCS_M12(iAlt)

        MRight_M12(iAlt,iSpecies) = &
           VelRight_M12(iAlt,iSpecies)/MeanCS_M12(iAlt)

        MLeft_P12(iAlt,iSpecies) = &
           VelLeft_P12(iAlt,iSpecies)/MeanCS_P12(iAlt)

        MRight_P12(iAlt,iSpecies) = &
           VelRight_P12(iAlt,iSpecies)/MeanCS_P12(iAlt)
    
     enddo 
   enddo 

     do iSpecies = 1, nSpecies
       M2Bar_M12(1:nAlts,iSpecies) = &
        0.5*(MLeft_M12(1:nAlts,iSpecies)**2.0 + &
             MRight_M12(1:nAlts,iSpecies)**2.0 )

        M2Bar_P12(1:nAlts,iSpecies) = &
           0.5*(MLeft_P12(1:nAlts,iSpecies)**2.0 + &
               MRight_P12(1:nAlts,iSpecies)**2.0 )
     enddo 

     do iSpecies = 1, nSpecies
      do iAlt = 1, nAlts 

        M2Zero_M12(iAlt,iSpecies) = &
              min(1.0, max(M2Bar_M12(iAlt,iSpecies), MInf)) 
        M2Zero_P12(iAlt,iSpecies) = &
              min(1.0, max(M2Bar_P12(iAlt,iSpecies), MInf)) 

        MZero_M12(iAlt,iSpecies) = sqrt(M2Zero_M12(iAlt,iSpecies))
        MZero_P12(iAlt,iSpecies) = sqrt(M2Zero_P12(iAlt,iSpecies))

      enddo 
     enddo 

     do iSpecies = 1, nSpecies
      do iAlt = 1, nAlts 

        FA_M12(iAlt,iSpecies) = &
               MZero_M12(iAlt,iSpecies)*(2.0 - MZero_M12(iAlt,iSpecies))
        FA_P12(iAlt,iSpecies) = &
               MZero_P12(iAlt,iSpecies)*(2.0 - MZero_P12(iAlt,iSpecies))

      enddo 
     enddo 

   do iSpecies = 1, nSpecies
     do iAlt = 1, nAlts

       MF1P_Left_M12(iAlt,iSpecies) = &
             0.5*(MLeft_M12(iAlt,iSpecies) + abs(MLeft_M12(iAlt,iSpecies)) )
       MF1N_Left_M12(iAlt,iSpecies) = &
             0.5*(MLeft_M12(iAlt,iSpecies) - abs(MLeft_M12(iAlt,iSpecies)) )

       MF1P_Right_M12(iAlt,iSpecies) = &
             0.5*(MRight_M12(iAlt,iSpecies) + abs(MRight_M12(iAlt,iSpecies)) )
       MF1N_Right_M12(iAlt,iSpecies) = &
             0.5*(MRight_M12(iAlt,iSpecies) - abs(MRight_M12(iAlt,iSpecies)) )

       MF1P_Left_P12(iAlt,iSpecies) = &
             0.5*(MLeft_P12(iAlt,iSpecies) + abs(MLeft_P12(iAlt,iSpecies)) )
       MF1N_Left_P12(iAlt,iSpecies) = &
             0.5*(MLeft_P12(iAlt,iSpecies) - abs(MLeft_P12(iAlt,iSpecies)) )

       MF1P_Right_P12(iAlt,iSpecies) = &
             0.5*(MRight_P12(iAlt,iSpecies) + abs(MRight_P12(iAlt,iSpecies)) )
       MF1N_Right_P12(iAlt,iSpecies) = &
             0.5*(MRight_P12(iAlt,iSpecies) - abs(MRight_P12(iAlt,iSpecies)) )
     enddo 
  enddo 

   do iSpecies = 1, nSpecies
     do iAlt = 1, nAlts

      MF2P_Left_M12(iAlt,iSpecies) =  0.25*(MLeft_M12(iAlt,iSpecies) + 1.0)**2.0
      MF2N_Left_M12(iAlt,iSpecies) = -0.25*(MLeft_M12(iAlt,iSpecies) - 1.0)**2.0

      MF2P_Right_M12(iAlt,iSpecies) =  0.25*(MRight_M12(iAlt,iSpecies) + 1.0)**2.0
      MF2N_Right_M12(iAlt,iSpecies) = -0.25*(MRight_M12(iAlt,iSpecies) - 1.0)**2.0

      MF2P_Left_P12(iAlt,iSpecies) =  0.25*(MLeft_P12(iAlt,iSpecies) + 1.0)**2.0
      MF2N_Left_P12(iAlt,iSpecies) = -0.25*(MLeft_P12(iAlt,iSpecies) - 1.0)**2.0

      MF2P_Right_P12(iAlt,iSpecies) =  0.25*(MRight_P12(iAlt,iSpecies) + 1.0)**2.0
      MF2N_Right_P12(iAlt,iSpecies) = -0.25*(MRight_P12(iAlt,iSpecies) - 1.0)**2.0

     enddo 
  enddo 


!!! Begin 4th Order Mach Number Polynomials
  do iSpecies = 1, nSpecies
    do iAlt = 1, nAlts 
      if ( abs(MLeft_M12(iAlt,iSpecies)) .ge. 1.0) then 
         MF4P_Left_M12(iAlt,iSpecies) = MF1P_Left_M12(iAlt,iSpecies)
      else
         MF4P_Left_M12(iAlt,iSpecies) = MF2P_Left_M12(iAlt,iSpecies)*&
                   (1.0 - 16.0*LiouBeta*MF2N_Left_M12(iAlt,iSpecies))
      endif 

      if ( abs(MRight_M12(iAlt,iSpecies)) .ge. 1.0) then 
         MF4N_Right_M12(iAlt,iSpecies) = MF1N_Right_M12(iAlt,iSpecies)
      else
         MF4N_Right_M12(iAlt,iSpecies) = MF2N_Right_M12(iAlt,iSpecies)*&
                    (1.0 + 16.0*LiouBeta*MF2P_Right_M12(iAlt,iSpecies))
      endif 

      if ( abs(MLeft_P12(iAlt,iSpecies)) .ge. 1.0) then 
         MF4P_Left_P12(iAlt,iSpecies) = MF1P_Left_P12(iAlt,iSpecies)
      else
         MF4P_Left_P12(iAlt,iSpecies) = MF2P_Left_P12(iAlt,iSpecies)*&
                   (1.0 - 16.0*LiouBeta*MF2N_Left_P12(iAlt,iSpecies))
      endif 

      if ( abs(MRight_P12(iAlt,iSpecies)) .ge. 1.0) then 
         MF4N_Right_P12(iAlt,iSpecies) = MF1N_Right_P12(iAlt,iSpecies)
      else
         MF4N_Right_P12(iAlt,iSpecies) = MF2N_Right_P12(iAlt,iSpecies)*&
                    (1.0 + 16.0*LiouBeta*MF2P_Right_P12(iAlt,iSpecies))
      endif 


    enddo 
  enddo 

!!! Begin 2nd Order Mach Number Polynomials
  do iSpecies = 1, nSpecies
    do iAlt = 1, nAlts 

       ModifiedZeta(iAlt,iSpecies) = 1.0

!!!! Hydrostatic Species Use This one
       MPress_M12(iAlt,iSpecies) = &
         LiouKpS(iAlt,iSpecies)*max( (1.0 - M2Bar_M12(iAlt,iSpecies)), 0.0)*&
         ((PressureSRight_M12(iAlt, iSpecies) - &
         HydroPressureSRight_M12(iAlt,iSpecies) ) - &
         (PressureSLeft_M12(iAlt, iSpecies) - &
          HydroPressureSLeft_M12(iAlt,iSpecies) ) )/ &
         ( MeanRhoS_M12(iAlt,iSpecies)*MeanCS_M12(iAlt)*&
         (FA_M12(iAlt,iSpecies)*MeanCS_M12(iAlt) + &
          ModifiedZeta(iAlt,iSpecies)*dAlt_C(iAlt)/DtIn)  ) 

       MPress_P12(iAlt,iSpecies) = &
          LiouKpS(iAlt,iSpecies)*max( (1.0 - M2Bar_P12(iAlt,iSpecies)), 0.0)*&
          (  (PressureSRight_P12(iAlt, iSpecies) - &
          HydroPressureSRight_P12(iAlt,iSpecies) ) - &
          (PressureSLeft_P12(iAlt, iSpecies) - &
          HydroPressureSLeft_P12(iAlt,iSpecies) ) )/ &
          ( MeanRhoS_P12(iAlt,iSpecies)*MeanCS_P12(iAlt)*&
          (FA_P12(iAlt,iSpecies)*MeanCS_P12(iAlt) + &
           ModifiedZeta(iAlt,iSpecies)*dAlt_C(iAlt)/DtIn)  ) 

!!!!! Non-Hydrostatic Species Should use this version

       if (.not. SubtractHydrostatic(iAlt,iSpecies)) then

       MPress_P12(iAlt,iSpecies) = &
          LiouKpS(iAlt,iSpecies)*max( (1.0 - M2Bar_P12(iAlt,iSpecies)), 0.0)*&
          (PressureSRight_P12(iAlt, iSpecies) - PressureSLeft_P12(iAlt,iSpecies) )/&
          ( MeanRhoS_P12(iAlt,iSpecies)*MeanCS_P12(iAlt)*&
          (FA_P12(iAlt,iSpecies)*MeanCS_P12(iAlt) + &
          ModifiedZeta(iAlt,iSpecies)*dAlt_C(iAlt)/DtIn)  ) 

       MPress_M12(iAlt,iSpecies) = &
         LiouKpS(iAlt,iSpecies)*max( (1.0 - M2Bar_M12(iAlt,iSpecies)), 0.0)*&
         (PressureSRight_M12(iAlt, iSpecies) - PressureSLeft_M12(iAlt,iSpecies) )/&
         ( MeanRhoS_M12(iAlt,iSpecies)*MeanCS_M12(iAlt)*&
         (FA_M12(iAlt,iSpecies)*MeanCS_M12(iAlt) + &
          ModifiedZeta(iAlt,iSpecies)*dAlt_C(iAlt)/DtIn)  ) 


        endif 
     enddo ! iSpecies
   enddo ! iAlt


   do iAlt = 1, nAlts 
      do iSpecies = 1, nSpecies

          InterfaceMach_M12(iAlt,iSpecies) =  &
                 MF4P_Left_M12(iAlt,iSpecies) + MF4N_Right_M12(iAlt,iSpecies) &
                    - MPress_M12(iAlt,iSpecies)

          LiouNumericalVelocity_M12(iAlt,iSpecies) = &
                   MeanCS_M12(iAlt)*InterfaceMach_M12(iAlt,iSpecies) 

          InterfaceMach_P12(iAlt,iSpecies) =  &
                 MF4P_Left_P12(iAlt,iSpecies) + MF4N_Right_P12(iAlt,iSpecies) &
                    - MPress_P12(iAlt,iSpecies)
 
          LiouNumericalVelocity_P12(iAlt,iSpecies) = &
                   MeanCS_P12(iAlt)*InterfaceMach_P12(iAlt,iSpecies) 
      enddo ! iSpecies 
   enddo  ! iAlt

      NumericalVelocity_M12(1:nAlts,1:nSpecies) = &
  LiouNumericalVelocity_M12(1:nAlts,1:nSpecies)

      NumericalVelocity_P12(1:nAlts,1:nSpecies) = &
  LiouNumericalVelocity_P12(1:nAlts,1:nSpecies)


!! NUMERICAL PRESSURE
  do iSpecies = 1, nSpecies
       do iAlt = 1, nAlts

             NumericalPressure_P12(iAlt,iSpecies) = &
                  (MeanPressureS_P12(iAlt,iSpecies) - &
              MeanHydroPressureS_P12(iAlt,iSpecies) )&
                - 0.5*Ku(iAlt,iSpecies)*MeanCS_P12(iAlt)*&
               ( (RhoSRight_P12(iAlt,iSpecies) )*VelRight_P12(iAlt,iSpecies) - &
                  (RhoSLeft_P12(iAlt,iSpecies) )* VelLeft_P12(iAlt,iSpecies) )

             NumericalPressure_M12(iAlt,iSpecies) = &
               (MeanPressureS_M12(iAlt,iSpecies) - &
           MeanHydroPressureS_M12(iAlt,iSpecies) )&
                - 0.5*Ku(iAlt,iSpecies)*MeanCS_M12(iAlt)*&
               ( (RhoSRight_M12(iAlt,iSpecies) )*VelRight_M12(iAlt,iSpecies) - &
                  (RhoSLeft_M12(iAlt,iSpecies) )* VelLeft_M12(iAlt,iSpecies) )

      
       if (.not. SubtractHydrostatic(iAlt,iSpecies)) then

             NumericalPressure_P12(iAlt,iSpecies) = &
               (MeanPressureS_P12(iAlt,iSpecies) )&
                - 0.5*Ku(iAlt,iSpecies)*MeanCS_P12(iAlt)*&
               ( (RhoSRight_P12(iAlt,iSpecies) )*VelRight_P12(iAlt,iSpecies) - &
                  (RhoSLeft_P12(iAlt,iSpecies) )* VelLeft_P12(iAlt,iSpecies) )

             NumericalPressure_M12(iAlt,iSpecies) = &
               (MeanPressureS_M12(iAlt,iSpecies) )&
                - 0.5*Ku(iAlt,iSpecies)*MeanCS_M12(iAlt)*&
               ( (RhoSRight_M12(iAlt,iSpecies) )*VelRight_M12(iAlt,iSpecies) - &
                  (RhoSLeft_M12(iAlt,iSpecies) )* VelLeft_M12(iAlt,iSpecies) )

       endif 
           
       enddo !!iAlt = 1, nAlts
  enddo !!! nSpecies


  do iSpecies = 1, nSpecies
     do iAlt = 1, nAlts 
        if (NumericalVelocity_P12(iAlt,iSpecies) .ge. 0.0) then 
           RhoSFlux_P12(iAlt,iSpecies) = &
               (RhoSLeft_P12(iAlt,iSpecies)*NumericalVelocity_P12(iAlt,iSpecies)) 

           MomentumSFlux_P12(iAlt,iSpecies) = &
               (RhoSLeft_P12(iAlt,iSpecies)*VelLeft_P12(iAlt,iSpecies)*&
               NumericalVelocity_P12(iAlt,iSpecies)) 
        else
           RhoSFlux_P12(iAlt,iSpecies) = &
               (RhoSRight_P12(iAlt,iSpecies)*NumericalVelocity_P12(iAlt,iSpecies)) 

           MomentumSFlux_P12(iAlt,iSpecies) = &
               (RhoSRight_P12(iAlt,iSpecies)*VelRight_P12(iAlt,iSpecies)*&
               NumericalVelocity_P12(iAlt,iSpecies)) 
        endif
        if (NumericalVelocity_M12(iAlt,iSpecies) .ge. 0.0) then 
           RhoSFlux_M12(iAlt,iSpecies) = &
               (RhoSLeft_M12(iAlt,iSpecies)*NumericalVelocity_M12(iAlt,iSpecies)) 

           MomentumSFlux_M12(iAlt,iSpecies) = &
               (RhoSLeft_M12(iAlt,iSpecies)*VelLeft_M12(iAlt,iSpecies)*&
               NumericalVelocity_M12(iAlt,iSpecies)) 
        else
           RhoSFlux_M12(iAlt,iSpecies) = &
               (RhoSRight_M12(iAlt,iSpecies)*NumericalVelocity_M12(iAlt,iSpecies)) 

           MomentumSFlux_M12(iAlt,iSpecies) = &
               (RhoSRight_M12(iAlt,iSpecies)*VelRight_M12(iAlt,iSpecies)*&
               NumericalVelocity_M12(iAlt,iSpecies)) 
        endif
     enddo 
  enddo 


    BulkNumericalVelocity_P12(1:nAlts) = 0.0
    BulkNumericalVelocity_M12(1:nAlts) = 0.0

    do iSpecies = 1, nSpecies
        BulkNumericalVelocity_P12(1:nAlts) = &
        BulkNumericalVelocity_P12(1:nAlts) + &
           ( RhoSLeft_P12(1:nAlts,iSpecies) + RhoSRight_P12(1:nAlts,iSpecies) )*&
             NumericalVelocity_P12(1:nAlts,iSpecies)/&
             (RhoLeft_P12(1:nAlts) + RhoRight_P12(1:nAlts) )

        BulkNumericalVelocity_M12(1:nAlts) = &
        BulkNumericalVelocity_M12(1:nAlts) + &
           ( RhoSLeft_M12(1:nAlts,iSpecies) + RhoSRight_M12(1:nAlts,iSpecies) )*&
             NumericalVelocity_M12(1:nAlts,iSpecies)/&
             (RhoLeft_M12(1:nAlts) + RhoRight_M12(1:nAlts) )
    enddo 


    do iAlt = 1, nAlts 
       if (BulkNumericalVelocity_P12(iAlt) .ge. 0.0) then 
           EnergyFlux_P12(iAlt) = &
             ( ELeft_P12(iAlt) + PLeft_P12(iAlt) ) &
               *BulkNumericalVelocity_P12(iAlt) 

           Momentum_P12(iAlt,1:3) = &
                RhoLeft_P12(iAlt)*VelGDLeft_P12(iAlt,1:3)*&
                BulkNumericalVelocity_P12(iAlt) 
       else
           EnergyFlux_P12(iAlt) = &
             ( ERight_P12(iAlt) + PRight_P12(iAlt) ) &
               *BulkNumericalVelocity_P12(iAlt) 

           Momentum_P12(iAlt,1:3) = &
                RhoRight_P12(iAlt)*VelGDRight_P12(iAlt,1:3)*&
                BulkNumericalVelocity_P12(iAlt) 
       endif
       if (BulkNumericalVelocity_M12(iAlt) .ge. 0.0) then 
           EnergyFlux_M12(iAlt) = &
             ( ELeft_M12(iAlt) + PLeft_M12(iAlt) ) &
               *BulkNumericalVelocity_M12(iAlt) 

           Momentum_M12(iAlt,1:3) = &
                RhoLeft_M12(iAlt)*VelGDLeft_M12(iAlt,1:3)*&
                BulkNumericalVelocity_M12(iAlt) 
       else
           EnergyFlux_M12(iAlt) = &
             ( ERight_M12(iAlt) + PRight_M12(iAlt) ) &
               *BulkNumericalVelocity_M12(iAlt) 

           Momentum_M12(iAlt,1:3) = &
                RhoRight_M12(iAlt)*VelGDRight_M12(iAlt,1:3)*&
                BulkNumericalVelocity_M12(iAlt) 
       endif
    enddo !! End iAlt Loop 


    do iAlt = 1, nAlts
          LeftRadius(iAlt) = 0.5*(RadDist(iAlt) + RadDist(iAlt-1))
          RightRadius(iAlt) = 0.5*(RadDist(iAlt) + RadDist(iAlt+1))
    enddo 

    do iAlt = 1, nAlts
          AreaFunction_P12(iAlt) = RightRadius(iAlt)**2.0
          AreaFunction_M12(iAlt) = LeftRadius(iAlt)**2.0
          LocalCellVolume(iAlt) = &
             (1.0/3.0)*(RightRadius(iAlt)**3.0 - LeftRadius(iAlt)**3.0)
    enddo 

    do iSpecies = 1, nSpecies
       RhoSFlux(1:nAlts,iSpecies) = &
            ( (AreaFunction_P12(1:nAlts))*RhoSFlux_P12(1:nAlts,iSpecies) - &
              (AreaFunction_M12(1:nAlts))*RhoSFlux_M12(1:nAlts,iSpecies) )/&
              LocalCellVolume(1:nAlts)

       MomentumSFlux(1:nAlts,iSpecies) = &
            ( (AreaFunction_P12(1:nAlts))*MomentumSFlux_P12(1:nAlts,iSpecies) - &
              (AreaFunction_M12(1:nAlts))*MomentumSFlux_M12(1:nAlts,iSpecies) )/&
              LocalCellVolume(1:nAlts)

       MomentumSFlux(1:nAlts,iSpecies) = &
       MomentumSFlux(1:nAlts,iSpecies) + &
           (NumericalPressure_P12(1:nAlts,iSpecies) - &
            NumericalPressure_M12(1:nAlts,iSpecies))/dAlt_C(1:nAlts)
    enddo 

    EnergyFlux(1:nAlts) = &
         ( (AreaFunction_P12(1:nAlts))*EnergyFlux_P12(1:nAlts) - &
           (AreaFunction_M12(1:nAlts))*EnergyFlux_M12(1:nAlts))/&
           LocalCellVolume(1:nAlts)

     do iDim = 1, 3
       MomentumFlux(1:nAlts,iDim) = &
            ( AreaFunction_P12(1:nAlts)*Momentum_P12(1:nAlts,iDim) - &
              AreaFunction_M12(1:nAlts)*Momentum_M12(1:nAlts,iDim))/&
              LocalCellVolume(1:nAlts)
     enddo 

end subroutine calc_all_fluxes_hydro


 subroutine calc_kt_facevalues(Var, VarLeft_M12, VarRight_M12, &
                                    VarLeft_P12, VarRight_P12)
 
 
   use ModVertical, only: dAlt_F, InvDAlt_F
   use ModSizeGITM, only: nAlts
   use ModLimiterGitm
 
   implicit none
   
   real, intent(in) :: Var(-1:nAlts+2)
   real, intent(out):: VarLeft_P12(1:nAlts), VarRight_P12(1:nAlts)
   real, intent(out):: VarLeft_M12(1:nAlts), VarRight_M12(1:nAlts)
 
   real :: dVarUp, dVarDown, dVarLimited(0:nAlts+1)
 
   integer :: i

!!!!!! TREAT THE LOWER BOUNDARY SEPARATELY
!        write(*,*) 'IN CALC_FACEVALUES'

        do i=0,nAlts+1
             dVarUp            = (Var(i+1) - Var(i  ))   * InvDAlt_F(i+1)
             dVarDown          = (Var(i  ) - Var(i-1)) * InvDAlt_F(i)
             dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)
        end do

        do i=1,nAlts
            VarLeft_M12(i) = Var(i-1) + 0.5*dVarLimited(i-1) * dAlt_F(i-1)
           VarRight_M12(i) = Var(i  ) - 0.5*dVarLimited(i  ) * dAlt_F(i  )
      
            VarLeft_P12(i) = Var(i  ) + 0.5*dVarLimited(i  ) * dAlt_F(i  )
           VarRight_P12(i) = Var(i+1) - 0.5*dVarLimited(i+1) * dAlt_F(i+1)
        end do
 
 
!        do i=1,nAlts
!            VarLeft_M12(i) = Var(i-1) + 0.5*dVarLimited(i-1) * dAlt_F(i  )
!           VarRight_M12(i) = Var(i  ) - 0.5*dVarLimited(i  ) * dAlt_F(i  )
!      
!            VarLeft_P12(i) = Var(i  ) + 0.5*dVarLimited(i  ) * dAlt_F(i+1)
!           VarRight_P12(i) = Var(i+1) - 0.5*dVarLimited(i+1) * dAlt_F(i+1)
!        end do

 end subroutine calc_kt_facevalues


