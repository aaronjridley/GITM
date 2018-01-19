!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine advance_horizontal(iBlock)

  use ModConstants, only : pi
  use ModSizeGitm
  use ModPlanet, only : nSpecies, nIonsAdvect, OmegaBody
  use ModGITM
  use ModInputs
  use ModSources, only : HorizontalTempSource
  
  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iIon, iLon, iLat, iSpecies
  logical :: IsFound

  real :: MaxDiff

  real :: cp_C(1:nLons,1:nLats)
  real :: Rho_C(-1:nLons+2,-1:nLats+2)
  real :: Temp_C(-1:nLons+2,-1:nLats+2)
  real :: gamma_C(-1:nLons+2,-1:nLats+2)
  real :: Vel_CD(-1:nLons+2,-1:nLats+2,3)
  real :: IVel_CD(-1:nLons+2,-1:nLats+2,3)
  real :: Num_CV(-1:nLons+2,-1:nLats+2,nSpecies)
  real :: VertVel_CV(-1:nLons+2,-1:nLats+2,nSpecies)
  real :: INum_CV(-1:nLons+2,-1:nLats+2,nIonsAdvect)

  real :: NewRho_C(nLons,nLats)
  real :: NewTemp_C(nLons,nLats)
  real :: NewVel_CD(nLons,nLats, 3)
  real :: NewNum_CV(nLons,nLats, nSpecies)
  real :: NewVertVel_CV(nLons,nLats, nSpecies)
  real :: NewINum_CV(nLons,nLats, nIonsAdvect)

  ! Vertical derivative of current variable (needed for topography only)
  real :: dVarDAlt_C(nLons,nLats) 

  ! JMB: 06/2016.  
  ! The following Variables enable us to use RK-4 time-stepping
  real :: OrigRho_C(nLons,nLats)
  real :: OrigTemp_C(nLons,nLats)
  real :: OrigVel_CD(nLons,nLats, 3)
  real :: OrigNum_CV(nLons,nLats, nSpecies)
  real :: OrigVertVel_CV(nLons,nLats, nSpecies)
  real :: OrigINum_CV(nLons,nLats, nIonsAdvect)

  real :: UpdatedRho_C(nLons,nLats)
  real :: UpdatedTemp_C(nLons,nLats)
  real :: UpdatedVel_CD(nLons,nLats, 3)
  real :: UpdatedNum_CV(nLons,nLats, nSpecies)
  real :: UpdatedVertVel_CV(nLons,nLats, nSpecies)
  real :: UpdatedINum_CV(nLons,nLats, nIonsAdvect)

  real :: K1Rho_C(nLons,nLats)
  real :: K1Temp_C(nLons,nLats)
  real :: K1Vel_CD(nLons,nLats, 3)
  real :: K1Num_CV(nLons,nLats, nSpecies)
  real :: K1VertVel_CV(nLons,nLats, nSpecies)
  real :: K1INum_CV(nLons,nLats, nIonsAdvect)

  real :: K2Rho_C(nLons,nLats)
  real :: K2Temp_C(nLons,nLats)
  real :: K2Vel_CD(nLons,nLats, 3)
  real :: K2Num_CV(nLons,nLats, nSpecies)
  real :: K2VertVel_CV(nLons,nLats, nSpecies)
  real :: K2INum_CV(nLons,nLats, nIonsAdvect)

  real :: K3Rho_C(nLons,nLats)
  real :: K3Temp_C(nLons,nLats)
  real :: K3Vel_CD(nLons,nLats, 3)
  real :: K3Num_CV(nLons,nLats, nSpecies)
  real :: K3VertVel_CV(nLons,nLats, nSpecies)
  real :: K3INum_CV(nLons,nLats, nIonsAdvect)

  real :: K4Rho_C(nLons,nLats)
  real :: K4Temp_C(nLons,nLats)
  real :: K4Vel_CD(nLons,nLats, 3)
  real :: K4Num_CV(nLons,nLats, nSpecies)
  real :: K4VertVel_CV(nLons,nLats, nSpecies)
  real :: K4INum_CV(nLons,nLats, nIonsAdvect)
  !----------------------------------------------------------------------------
  MaxDiff = 0.0

  call report("advance_horizontal",2)

  do iAlt=1,nAlts

     cp_c       = cp(:,:,iAlt,iBlock)
     Gamma_c    = gamma(:,:,iAlt,iBlock) 
     Rho_C      = Rho(:,:,iAlt,iBlock)
     Vel_CD     = Velocity(:,:,iAlt,:,iBlock)
     VertVel_CV = VerticalVelocity(:,:,iAlt,1:nSpecies,iBlock)
     Num_CV     = nDensityS(:,:,iAlt,1:nSpecies,iBlock)
     Temp_C     = Temperature(:,:,iAlt,iBlock)

     IVel_CD = IVelocity(:,:,iAlt,:,iBlock)
     INum_CV = IDensityS(:,:,iAlt,1:nIonsAdvect,iBlock)
     
     !------------
     ! BEGIN RK-4
     NewRho_C      = Rho_C(1:nLons,1:nLats)
     NewVel_CD     = Vel_CD(1:nLons,1:nLats,:)
     NewNum_CV     = Num_CV(1:nLons,1:nLats,:)
     NewVertVel_CV = VertVel_CV(1:nLons,1:nLats,:)
     NewTemp_C     = Temp_C(1:nLons,1:nLats)
     NewINum_CV    = INum_CV(1:nLons,1:nLats,:)

     OrigRho_C      = Rho_C(1:nLons,1:nLats)
     OrigVel_CD     = Vel_CD(1:nLons,1:nLats,:)
     OrigNum_CV     = Num_CV(1:nLons,1:nLats,:)
     OrigVertVel_CV = VertVel_CV(1:nLons,1:nLats,:)
     OrigTemp_C     = Temp_C(1:nLons,1:nLats)
     OrigINum_CV    = INum_CV(1:nLons,1:nLats,:)

!!! 1st Update Step
     call horizontal_solver
     
     K1Rho_C      = NewRho_C(1:nLons,1:nLats) - Rho_C(1:nLons,1:nLats)
     K1Vel_CD     = NewVel_CD(1:nLons,1:nLats,1:3) - &
                   Vel_CD(1:nLons,1:nLats,1:3)
     K1Num_CV     = NewNum_CV(1:nLons,1:nLats,1:nSpecies) - &
                   Num_CV(1:nLons,1:nLats,1:nSpecies)
     K1VertVel_CV = NewVertVel_CV(1:nLons,1:nLats,1:nSpecies) - &
                   VertVel_CV(1:nLons,1:nLats,1:nSpecies)
     K1Temp_C     = NewTemp_C(1:nLons,1:nLats) - Temp_C(1:nLons,1:nLats)
     K1INum_CV     = NewINum_CV(1:nLons,1:nLats,1:nIonsAdvect) - &
                    INum_CV(1:nLons,1:nLats,1:nIonsAdvect)

     UpdatedRho_C   = 0.5*K1Rho_C(1:nLons,1:nLats) + &
                        OrigRho_C(1:nLons,1:nLats)
     UpdatedVel_CD  = 0.5*K1Vel_CD(1:nLons,1:nLats,1:3) + &
                        OrigVel_CD(1:nLons,1:nLats,1:3)
     UpdatedNum_CV  = 0.5*K1Num_CV(1:nLons,1:nLats,1:nSpecies) + &
                        OrigNum_CV(1:nLons,1:nLats,1:nSpecies)
     UpdatedVertVel_CV = 0.5*K1VertVel_CV(1:nLons,1:nLats,1:nSpecies) + &
                           OrigVertVel_CV(1:nLons,1:nLats,1:nSpecies)
     UpdatedTemp_C  = 0.5*K1Temp_C(1:nLons,1:nLats) + &
                        OrigTemp_C(1:nLons,1:nLats)
     UpdatedINum_CV = 0.5*K1INum_CV(1:nLons,1:nLats,1:nIonsAdvect) + &
                        OrigINum_CV(1:nLons,1:nLats,1:nIonsAdvect)

     Rho_C(1:nLons,1:nLats)      = UpdatedRho_C
     Vel_CD(1:nLons,1:nLats,:)     = UpdatedVel_CD
     Num_CV(1:nLons,1:nLats,:)     = UpdatedNum_CV
     VertVel_CV(1:nLons,1:nLats,:) = UpdatedVertVel_CV
     Temp_C(1:nLons,1:nLats)     = UpdatedTemp_C
     INum_CV(1:nLons,1:nLats,:)    = UpdatedINum_CV

     NewRho_C      = UpdatedRho_C
     NewVel_CD     = UpdatedVel_CD
     NewNum_CV     = UpdatedNum_CV
     NewVertVel_CV = UpdatedVertVel_CV
     NewTemp_C     = UpdatedTemp_C
     NewINum_CV    = UpdatedINum_CV


!!! 2nd Update Step
     call horizontal_solver
     
     K2Rho_C      = NewRho_C(1:nLons,1:nLats) - Rho_C(1:nLons,1:nLats)
     K2Vel_CD     = NewVel_CD(1:nLons,1:nLats,1:3) - &
                   Vel_CD(1:nLons,1:nLats,1:3)
     K2Num_CV     = NewNum_CV(1:nLons,1:nLats,1:nSpecies) - &
                   Num_CV(1:nLons,1:nLats,1:nSpecies)
     K2VertVel_CV = NewVertVel_CV(1:nLons,1:nLats,1:nSpecies) - &
                   VertVel_CV(1:nLons,1:nLats,1:nSpecies)
     K2Temp_C     = NewTemp_C(1:nLons,1:nLats) - Temp_C(1:nLons,1:nLats)
     K2INum_CV     = NewINum_CV(1:nLons,1:nLats,1:nIonsAdvect) - &
                    INum_CV(1:nLons,1:nLats,1:nIonsAdvect)


     UpdatedRho_C   = 0.5*K2Rho_C(1:nLons,1:nLats) + &
                        OrigRho_C(1:nLons,1:nLats)
     UpdatedVel_CD  = 0.5*K2Vel_CD(1:nLons,1:nLats,1:3) + &
                        OrigVel_CD(1:nLons,1:nLats,1:3)
     UpdatedNum_CV  = 0.5*K2Num_CV(1:nLons,1:nLats,1:nSpecies) + &
                        OrigNum_CV(1:nLons,1:nLats,1:nSpecies)
     UpdatedVertVel_CV = 0.5*K2VertVel_CV(1:nLons,1:nLats,1:nSpecies) + &
                           OrigVertVel_CV(1:nLons,1:nLats,1:nSpecies)
     UpdatedTemp_C  = 0.5*K2Temp_C(1:nLons,1:nLats) + &
                        OrigTemp_C(1:nLons,1:nLats)
     UpdatedINum_CV = 0.5*K2INum_CV(1:nLons,1:nLats,1:nIonsAdvect) + &
                     OrigINum_CV(1:nLons,1:nLats,1:nIonsAdvect)

     Rho_C(1:nLons,1:nLats)      = UpdatedRho_C
     Vel_CD(1:nLons,1:nLats,:)     = UpdatedVel_CD
     Num_CV(1:nLons,1:nLats,:)     = UpdatedNum_CV
     VertVel_CV(1:nLons,1:nLats,:) = UpdatedVertVel_CV
     Temp_C(1:nLons,1:nLats)     = UpdatedTemp_C
     INum_CV(1:nLons,1:nLats,:)    = UpdatedINum_CV


     NewRho_C      = UpdatedRho_C
     NewVel_CD     = UpdatedVel_CD
     NewNum_CV     = UpdatedNum_CV
     NewVertVel_CV = UpdatedVertVel_CV
     NewTemp_C     = UpdatedTemp_C
     NewINum_CV    = UpdatedINum_CV

!!! 3rd Update Step
     call horizontal_solver
     
     K3Rho_C      = NewRho_C(1:nLons,1:nLats) - Rho_C(1:nLons,1:nLats)
     K3Vel_CD     = NewVel_CD(1:nLons,1:nLats,1:3) - &
                   Vel_CD(1:nLons,1:nLats,1:3)
     K3Num_CV     = NewNum_CV(1:nLons,1:nLats,1:nSpecies) - &
                   Num_CV(1:nLons,1:nLats,1:nSpecies)
     K3VertVel_CV = NewVertVel_CV(1:nLons,1:nLats,1:nSpecies) - &
                   VertVel_CV(1:nLons,1:nLats,1:nSpecies)
     K3Temp_C     = NewTemp_C(1:nLons,1:nLats) - Temp_C(1:nLons,1:nLats)
     K3INum_CV     = NewINum_CV(1:nLons,1:nLats,1:nIonsAdvect) - &
                    INum_CV(1:nLons,1:nLats,1:nIonsAdvect)

     UpdatedRho_C   = K3Rho_C(1:nLons,1:nLats) + &
                        OrigRho_C(1:nLons,1:nLats)
     UpdatedVel_CD  = K3Vel_CD(1:nLons,1:nLats,1:3) + &
                        OrigVel_CD(1:nLons,1:nLats,1:3)
     UpdatedNum_CV  = K3Num_CV(1:nLons,1:nLats,1:nSpecies) + &
                    OrigNum_CV(1:nLons,1:nLats,1:nSpecies)
     UpdatedVertVel_CV = K3VertVel_CV(1:nLons,1:nLats,1:nSpecies) + &
                           OrigVertVel_CV(1:nLons,1:nLats,1:nSpecies)
     UpdatedTemp_C  = K3Temp_C(1:nLons,1:nLats) + &
                        OrigTemp_C(1:nLons,1:nLats)
     UpdatedINum_CV = K3INum_CV(1:nLons,1:nLats,1:nIonsAdvect) + &
                     OrigINum_CV(1:nLons,1:nLats,1:nIonsAdvect)

     Rho_C(1:nLons,1:nLats)      = UpdatedRho_C
     Vel_CD(1:nLons,1:nLats,:)     = UpdatedVel_CD
     Num_CV(1:nLons,1:nLats,:)     = UpdatedNum_CV
     VertVel_CV(1:nLons,1:nLats,:) = UpdatedVertVel_CV
     Temp_C(1:nLons,1:nLats)     = UpdatedTemp_C
     INum_CV(1:nLons,1:nLats,:)    = UpdatedINum_CV


     NewRho_C      = UpdatedRho_C
     NewVel_CD     = UpdatedVel_CD
     NewNum_CV     = UpdatedNum_CV
     NewVertVel_CV = UpdatedVertVel_CV
     NewTemp_C     = UpdatedTemp_C
     NewINum_CV    = UpdatedINum_CV

!!! 4th (Final) Update Step
     call horizontal_solver

     K4Rho_C      = NewRho_C(1:nLons,1:nLats) - Rho_C(1:nLons,1:nLats)
     K4Vel_CD     = NewVel_CD(1:nLons,1:nLats,1:3) - &
                   Vel_CD(1:nLons,1:nLats,1:3)
     K4Num_CV     = NewNum_CV(1:nLons,1:nLats,1:nSpecies) - &
                   Num_CV(1:nLons,1:nLats,1:nSpecies)
     K4VertVel_CV = NewVertVel_CV(1:nLons,1:nLats,1:nSpecies) - &
                   VertVel_CV(1:nLons,1:nLats,1:nSpecies)
     K4Temp_C     = NewTemp_C(1:nLons,1:nLats) - Temp_C(1:nLons,1:nLats)
     K4INum_CV     = NewINum_CV(1:nLons,1:nLats,1:nIonsAdvect) - &
                    INum_CV(1:nLons,1:nLats,1:nIonsAdvect)

     UpdatedRho_C   = OrigRho_C + (1.0/6.0)*&
           (K1Rho_C + 2.0*K2Rho_C + 2.0*K3Rho_C + K4Rho_C) 

     UpdatedVel_CD   = OrigVel_CD + (1.0/6.0)*&
           (K1Vel_CD + 2.0*K2Vel_CD + 2.0*K3Vel_CD + K4Vel_CD) 

     UpdatedNum_CV   = OrigNum_CV + (1.0/6.0)*&
           (K1Num_CV + 2.0*K2Num_CV + 2.0*K3Num_CV + K4Num_CV) 

     UpdatedVertVel_CV   = OrigVertVel_CV + (1.0/6.0)*&
           (K1VertVel_CV + 2.0*K2VertVel_CV + 2.0*K3VertVel_CV + K4VertVel_CV) 

     UpdatedTemp_C   = OrigTemp_C + (1.0/6.0)*&
           (K1Temp_C + 2.0*K2Temp_C + 2.0*K3Temp_C + K4Temp_C) 

     UpdatedINum_CV   = OrigINum_CV + (1.0/6.0)*&
           (K1INum_CV + 2.0*K2INum_CV + 2.0*K3INum_CV + K4INum_CV) 

     NewRho_C      = UpdatedRho_C
     NewVel_CD     = UpdatedVel_CD
     NewNum_CV     = UpdatedNum_CV
     NewVertVel_CV = UpdatedVertVel_CV
     NewTemp_C     = UpdatedTemp_C
     NewINum_CV    = UpdatedINum_CV

     ! END RK-4
     ! ---------
     ! JMB:  06/2016
     ! After the RK-4 Update, we update the state variables
     Rho(1:nLons,1:nLats,iAlt,iBlock)                     = NewRho_C
     Velocity(1:nLons,1:nLats,iAlt,:,iBlock)              = NewVel_CD

     HorizontalTempSource(1:nLons,1:nLats,iAlt)           = &
          NewTemp_C-Temperature(1:nLons,1:nLats,iAlt,iBlock)

     Temperature(1:nLons,1:nLats,iAlt,iBlock)             = NewTemp_C
     VerticalVelocity(1:nLons,1:nLats,iAlt,1:nSpecies,iBlock)         = &
         NewVertVel_CV

     if (minval(NewNum_CV) < 0.0) then
        write(*,*) "Negative Density after horizontal advection!!"
        write(*,*) "Correcting...."
        do iLon = 1, nLons
           do iLat = 1, nLats
              IsFound = .false.
              do iSpecies = 1, nSpecies
                 if (NewNum_CV(iLon, iLat, iSpecies) < 0.0) then
                    write(*,*) "Species : ", iSpecies, iLon, iLat, iBlock
                    stop
                    NewNum_CV(iLon, iLat, iSpecies) = 1.0
                    IsFound=.true.
                 endif
              enddo
              If (IsFound) then
                 Rho(iLon,iLat,iAlt,iBlock) = 0.0
                 do iSpecies = 1, nSpecies
                    Rho(iLon,iLat,iAlt,iBlock) = &
                         Rho(iLon,iLat,iAlt,iBlock) + &
                         NewNum_CV(iLon, iLat, iSpecies)*Mass(iSpecies)
                 enddo
              endif
           enddo
        enddo
     endif
 
     NDensityS(1:nLons,1:nLats,iAlt,1:nSpecies,iBlock)    = NewNum_CV

     if (UseIonAdvection) then

        if (minval(NewINum_CV) < 1.0e2) then
!           write(*,*) "Negative Ion Density after horizontal advection!!"
!           write(*,*) "Correcting...."
           do iLon = 1, nLons
              do iLat = 1, nLats
                 do iIon = 1, nIonsAdvect
                    if (NewINum_CV(iLon, iLat, iIon) < 1.0e2) then
!                       write(*,*) "Location : ", &
!                            Longitude(iLon,iBlock)*180/pi, &
!                            Latitude(iLon,iBlock)*180/pi, &
!                            Altitude(iAlt)/1000.0, iIon
                       NewINum_CV(iLon, iLat, iIon) = 1.0e2
                    endif
                 enddo
              enddo
           enddo
        endif

        IDensityS(1:nLons,1:nLats,iAlt,1:nIonsAdvect,iBlock) = NewINum_CV
        !\
        ! New Electron Density
        !/
        IDensityS(1:nLons,1:nLats,iAlt,ie_,iBlock) = 0.0
        do iIon = 1, nIons-1
           IDensityS(1:nLons,1:nLats,iAlt,ie_,iBlock) = &
                IDensityS(1:nLons,1:nLats,iAlt,ie_,iBlock) + &
                IDensityS(1:nLons,1:nLats,iAlt,iIon,iBlock)
        enddo
     endif

  end do

  if (iDebugLevel > 2) &
       write(*,*) "===> advance_horizontal, MaxDiff, IVel, IDens, Dt : ", &
       MaxDiff, maxval(abs(IVel_CD)), &
       maxval(IDensityS(1:nLons,1:nLats,1:nAlts,1:nIonsAdvect,iBlock)), dt

contains

  subroutine horizontal_solver

    use ModInputs, only: UseCoriolis

    ! Solve horizontal equations for a single block

    integer :: iDim, iSpc

    real, dimension(-1:nLons+2,-1:nLats+2):: &
         Quantity 
    real, dimension(nLons,nLats)          :: &
         GradLonQuantity, GradLatQuantity, DiffQuantity
    real, dimension(nLons,nLats)          :: &
         GradLonRho_C, GradLatRho_C,   &
         DiffLonRho_C, DiffLatRho_C
    real, dimension(nLons,nLats)          :: &
         GradLonTemp_C, GradLatTemp_C, &
         DiffLonTemp_C, DiffLatTemp_C, DivVel_C
    real, dimension(nLons,nLats,3)        :: &
         GradLonVel_CD, GradLatVel_CD, &
         DiffLonVel_CD, DiffLatVel_CD
    real, dimension(nLons,nLats,nSpecies) :: &
         GradLonVertVel_CV, GradLatVertVel_CV, &
         DiffLonVertVel_CV, DiffLatVertVel_CV
    real, dimension(nLons,nLats,nSpecies) :: &
         GradLonNum_CV, GradLatNum_CV, &
         DiffLonNum_CV, DiffLatNum_CV
    real, dimension(nLons,nLats,nIonsAdvect) :: &
         GradLonINum_CV, GradLatINum_CV, &
         DiffLonINum_CV, DiffLatINum_CV

    real :: CoriolisSin, CoriolisCos, CentrifugalParameter

    real :: RhoTest, CosLat(nLats), SinLat(nLats)

    real :: HalfInvDAlt_C(nLons, nLats)

    integer :: iLon, iLat
    !--------------------------------------------------------------------------
    if(UseTopography)then
       HalfInvDAlt_C = 0.5/dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)

       dVarDAlt_C = HalfInvDAlt_C * &
            ( Rho(1:nLons,1:nLats,iAlt+1,iBlock) &
            - Rho(1:nLons,1:nLats,iAlt-1,iBlock) )
    end if
    call calc_rusanov_lons( Rho_C,  GradLonRho_C,  DiffLonRho_C)
    call calc_rusanov_lats( Rho_C,  GradLatRho_C,  DiffLatRho_C)

    if(UseTopography) dVarDAlt_C = HalfInvDAlt_C * &
         ( Temperature(1:nLons,1:nLats,iAlt+1,iBlock) &
         - Temperature(1:nLons,1:nLats,iAlt-1,iBlock) )
    call calc_rusanov_lons( Temp_C, GradLonTemp_C, DiffLonTemp_C)
    call calc_rusanov_lats( Temp_C, GradLatTemp_C, DiffLatTemp_C)

    do iDim = 1,3
       if(UseTopography) dVarDAlt_C = HalfInvDAlt_C * &
            ( Velocity(1:nLons,1:nLats,iAlt+1,iDim,iBlock) &
            - Velocity(1:nLons,1:nLats,iAlt-1,iDim,iBlock) )
       call calc_rusanov_lons(Vel_CD(:,:,iDim), &
            GradLonVel_CD(:,:,iDim), DiffLonVel_CD(:,:,iDim))
       call calc_rusanov_lats(Vel_CD(:,:,iDim), &
            GradLatVel_CD(:,:,iDim), DiffLatVel_CD(:,:,iDim))
    end do

    do iLat = 1, nLats
       DivVel_C(:,iLat) = &
            GradLatVel_CD(:,iLat,iNorth_) + GradLonVel_CD(:,iLat,iEast_) &
            - TanLatitude(iLat,iBlock) * Vel_CD(1:nLons,iLat,iNorth_) * &
            InvRadialDistance_GB(1:nLons,iLat,iAlt,iBlock)
    end do

    do iSpc = 1,nSpecies
       if(UseTopography) dVarDAlt_C = HalfInvDAlt_C * &
            ( nDensityS(1:nLons,1:nLats,iAlt+1,iSpc,iBlock) &
            - nDensityS(1:nLons,1:nLats,iAlt-1,iSpc,iBlock) )
       call calc_rusanov_lons(Num_CV(:,:,iSpc), &
            GradLonNum_CV(:,:,iSpc), DiffLonNum_CV(:,:,iSpc))
       call calc_rusanov_lats(Num_CV(:,:,iSpc), &
            GradLatNum_CV(:,:,iSpc), DiffLatNum_CV(:,:,iSpc))

       if(UseTopography) dVarDAlt_C = HalfInvDAlt_C * &
            ( VerticalVelocity(1:nLons,1:nLats,iAlt+1,iSpc,iBlock) &
            - VerticalVelocity(1:nLons,1:nLats,iAlt-1,iSpc,iBlock) )
       call calc_rusanov_lons( VertVel_CV(:,:,iSpc), &
            GradLonVertVel_CV(:,:,iSpc), DiffLonVertVel_CV(:,:,iSpc))
       call calc_rusanov_lats( VertVel_CV(:,:,iSpc), &
            GradLatVertVel_CV(:,:,iSpc), DiffLatVertVel_CV(:,:,iSpc))
       
    end do

    do iSpc = 1,nIonsAdvect
       if(UseTopography) dVarDAlt_C = HalfInvDAlt_C * &
            ( IDensityS(1:nLons,1:nLats,iAlt+1,iSpc,iBlock) &
            - IDensityS(1:nLons,1:nLats,iAlt-1,iSpc,iBlock) )
       call calc_rusanov_lons(INum_CV(:,:,iSpc), &
            GradLonINum_CV(:,:,iSpc), DiffLonINum_CV(:,:,iSpc))
       call calc_rusanov_lats(INum_CV(:,:,iSpc), &
            GradLatINum_CV(:,:,iSpc), DiffLatINum_CV(:,:,iSpc))
    end do

    SinLat = sin(Latitude(1:nLats,iBlock))
    CosLat = CosLatitude(1:nLats,iBlock)

    do iLat=1,nLats

       CoriolisSin = SinLat(iLat) * 2 * OmegaBody
       CoriolisCos = CosLat(iLat) * 2 * OmegaBody

       CentrifugalParameter = OmegaBody**2 * cosLat(iLat) * &
            sinLat(iLat)

       do iLon=1,nLons

          ! drho/dt = -(rho*div V + grad rho * V)
          !         = -[rho*(dv/dLat - v*tan(Lat) + (du/dLon)/cos(Lat))
          !             drho/dLat * vLat + drho/dLon * vLon/cos(Lat)]
          !            / r

          NewRho_C(iLon,iLat) = NewRho_C(iLon,iLat) - Dt * ( &
               Rho_C(iLon,iLat) * DivVel_C(iLon,iLat) &
               + GradLatRho_C(iLon,iLat)*Vel_CD(iLon,iLat,iNorth_) & 
               + GradLonRho_C(iLon,iLat)*Vel_CD(iLon,iLat,iEast_)) & 
               + Dt * (DiffLonRho_C(iLon,iLat)+DiffLatRho_C(iLon,iLat))

          do iSpc = 1, nSpecies
             NewNum_CV(iLon,iLat,iSpc) = NewNum_CV(iLon,iLat,iSpc) - Dt * ( &
                  Num_CV(iLon,iLat,iSpc) * DivVel_C(iLon,iLat) &
                  + GradLatNum_CV(iLon,iLat,iSpc)*Vel_CD(iLon,iLat,iNorth_) & 
                  + GradLonNum_CV(iLon,iLat,iSpc)*Vel_CD(iLon,iLat,iEast_)) & 
                  + Dt * (&
                  DiffLonNum_CV(iLon,iLat,iSpc)+DiffLatNum_CV(iLon,iLat,iSpc))
          enddo

          do iSpc = 1, nIonsAdvect
             NewINum_CV(iLon,iLat,iSpc) = NewINum_CV(iLon,iLat,iSpc) - Dt * ( &
!                  INum_CV(iLon,iLat,iSpc) * DivVel_C(iLon,iLat) &
                  + GradLatINum_CV(iLon,iLat,iSpc)*IVel_CD(iLon,iLat,iNorth_) &
                  + GradLonINum_CV(iLon,iLat,iSpc)*IVel_CD(iLon,iLat,iEast_)) &
                  + Dt * (&
                  DiffLonINum_CV(iLon,iLat,iSpc)+&
                  DiffLatINum_CV(iLon,iLat,iSpc))
          enddo

          RhoTest = sum(Mass(1:nSpecies) * NewNum_CV(iLon,iLat,1:nSpecies))

!          if (abs(RhoTest - NewRho_C(iLon,iLat))/RhoTest > 0.1) then
!             write(*,*) "Problem!! ", RhoTest, NewRho_C(iLon,iLat)
!!             call stop_gitm("Have to stop")
!          endif

          ! dv_phi/dt = -(V grad V + (1/rho) grad P)_phi
          ! (1/rho) grad p = grad T + T/rho grad rho 

          NewVel_CD(iLon,iLat,iEast_) = NewVel_CD(iLon,iLat,iEast_) - Dt * ( &
               Vel_CD(iLon,iLat,iNorth_)*GradLatVel_CD(iLon,iLat,iEast_) + &
               Vel_CD(iLon,iLat,iEast_ )*GradLonVel_CD(iLon,iLat,iEast_) + &
               Vel_CD(iLon,iLat,iEast_)*(Vel_CD(iLon,iLat,iUp_) &
               - TanLatitude(iLat,iBlock)*Vel_CD(iLon,iLat,iNorth_)) &
               * InvRadialDistance_GB(iLon,iLat,iAlt,iBlock) + &
               GradLonTemp_C(iLon,iLat) + & 
               GradLonRho_C(iLon,iLat)*Temp_C(iLon,iLat)/Rho_C(iLon,iLat)) &
               + Dt * (&
               DiffLonVel_CD(iLon,iLat,iEast_)+DiffLatVel_CD(iLon,iLat,iEast_))


          ! dv_theta/dt = -(V grad V + (1/rho) grad P)_theta
          ! (1/rho) grad p = grad T + T/rho grad rho 

          NewVel_CD(iLon,iLat,iNorth_) = NewVel_CD(iLon,iLat,iNorth_) &
               - Dt * ( &
               Vel_CD(iLon,iLat,iNorth_)*GradLatVel_CD(iLon,iLat,iNorth_) + &
               Vel_CD(iLon,iLat,iEast_ )*GradLonVel_CD(iLon,iLat,iNorth_) + &
               (Vel_CD(iLon,iLat,iNorth_)*Vel_CD(iLon,iLat,iUp_) &
               + TanLatitude(iLat,iBlock)*Vel_CD(iLon,iLat,iEast_)**2 &
               ) * InvRadialDistance_GB(iLon,iLat,iAlt,iBlock) + &
               GradLatTemp_C(iLon,iLat) + & 
               GradLatRho_C(iLon,iLat)*Temp_C(iLon,iLat)/Rho_C(iLon,iLat)) &
               + Dt * ( &
               DiffLonVel_CD(iLon,iLat,iNorth_)+ &
               DiffLatVel_CD(iLon,iLat,iNorth_))

!          if (iLon == 1) then
!             write(*,*) "vel before cor : ",NewVel_CD(iLon,iLat,iNorth_) 
!          endif

          ! dv_r/dt = -(V grad V)_r

          NewVel_CD(iLon,iLat,iUp_) = NewVel_CD(iLon,iLat,iUp_) &
               - Dt * ( &
               Vel_CD(iLon,iLat,iNorth_)*GradLatVel_CD(iLon,iLat,iUp_) + &
               Vel_CD(iLon,iLat,iEast_ )*GradLonVel_CD(iLon,iLat,iUp_)) &
               + Dt * (&
               DiffLonVel_CD(iLon,iLat,iUp_)+DiffLatVel_CD(iLon,iLat,iUp_))

          ! Same as bulk vertical velocity (above) but for each species

          do iSpc = 1, nSpecies
             NewVertVel_CV(iLon,iLat,iSpc) = NewVertVel_CV(iLon,iLat,iSpc) &
                  - Dt * ( &
                  Vel_CD(iLon,iLat,iNorth_)*GradLatVertVel_CV(iLon,iLat,iSpc)+&
                  Vel_CD(iLon,iLat,iEast_ )*GradLonVertVel_CV(iLon,iLat,iSpc))&
                  + Dt * ( &
                  DiffLatVertVel_CV(iLon,iLat,iSpc) + &
                  DiffLonVertVel_CV(iLon,iLat,iSpc))
          enddo

          if (UseCoriolis) then

             NewVel_CD(iLon,iLat,iEast_) = NewVel_CD(iLon,iLat,iEast_) + Dt*( &
                  + CoriolisSin * Vel_CD(iLon,iLat,iNorth_) &
                  - CoriolisCos * Vel_CD(iLon,iLat,iUp_))

             NewVel_CD(iLon,iLat,iNorth_) = NewVel_CD(iLon,iLat,iNorth_) - &
                  Dt*( &
                  CentrifugalParameter &
                  * RadialDistance_GB(iLon,iLat,iAlt,iBlock) &
                  + CoriolisSin * Vel_CD(iLon,iLat,iEast_))
          endif

          ! dT/dt = -(V.grad T + (gamma - 1) T div V

          NewTemp_C(iLon,iLat) = NewTemp_C(iLon,iLat) - Dt * ( &
               (gamma_c(iLon,iLat)-1) * Temp_C(iLon,iLat) &
               * DivVel_C(iLon,iLat) &
               + GradLatTemp_C(iLon,iLat)*Vel_CD(iLon,iLat,iNorth_) & 
               + GradLonTemp_C(iLon,iLat)*Vel_CD(iLon,iLat,iEast_)) & 
               + Dt * (DiffLonTemp_C(iLon,iLat)+DiffLatTemp_C(iLon,iLat))

       end do
    end do

  end subroutine horizontal_solver

  !=====================================================================
  subroutine calc_rusanov_lats(Var, GradVar, DiffVar)

    use ModSizeGitm
    use ModGITM

    implicit none
  
    real, intent(in)  :: Var(-1:nLons+2, -1:nLats+2)
    real, intent(out) :: GradVar(nLons, nLats)
    real, intent(inout) :: DiffVar(nLons, nLats)

    real :: DiffLocP(-1:nLats+1), InvDiffLocP(-1:nLats+1)

    real, dimension(1:nLats+1) :: VarNorth, VarSouth, DiffFlux
    real :: InvdLat(nLats), InvdLon
    real :: TempVar(-1:nLats+2)

    integer :: iLon

    ! Calculate gradient and diffusive flux with respect to latitude

    do iLon = 1, nLons
       TempVar(-1:nLats+2) = Var(iLon,-1:nLats+2)
       InvdLat = InvDLatDist_GB(iLon, 1:nLats,iAlt,iBlock)

       call calc_facevalues_lats(iLon, iAlt, iBlock, TempVar, &
            VarSouth, VarNorth)

!       call calc_facevalues_lats(iLon, iAlt, iBlock, Var(iLon,:), &
!            VarSouth, VarNorth)

       ! Gradient based on averaged Left/Right values

       GradVar(iLon,:) = 0.5 * &
            (VarSouth(2:nLats+1)+VarNorth(2:nLats+1) - &
            VarSouth(1:nLats)-VarNorth(1:nLats)) * InvdLat

       ! Rusanov/Lax-Friedrichs diffusive term
       DiffFlux = 0.5 * max(&
            cMax_GDB(iLon,0:nLats  ,iAlt,iNorth_,iBlock),&
            cMax_GDB(iLon,1:nLats+1,iAlt,iNorth_,iBlock)) &
            * (VarNorth - VarSouth)

       DiffVar(iLon,:) = &
            (DiffFlux(2:nLats+1) - DiffFlux(1:nLats)) * &
            InvdLat

    end do

    if(UseTopography) &
         GradVar = GradVar - dVarDAlt_C * dAltDLat_CB(:,:,iAlt,iBlock)

  end subroutine calc_rusanov_lats

  !===========================================================================
  subroutine calc_rusanov_lons(Var, GradVar, DiffVar)

    use ModSizeGitm
    use ModGITM

    implicit none

    real, intent(in)  :: Var(-1:nLons+2, -1:nLats+2)
    real, intent(out) :: GradVar(nLons, nLats)
    real, intent(out) :: DiffVar(nLons, nLats)

    real, dimension(1:nLons+1) :: VarEast,  VarWest,  DiffFlux

    integer :: iLat

    real :: InvdLat(nLats), InvdLon(nLons)
    !--------------------------------------------------------------------------

    ! Calculate gradient and diffusive flux with respect to longitude
    do iLat = 1, nLats

       InvdLon = InvDLonDist_GB(1:nLons,iLat,iAlt,iBlock)

       call calc_facevalues_lons(iLat, iAlt, iBlock, &
            Var(:,iLat), VarWest, VarEast)

       ! Gradient based on averaged West/East values

       GradVar(:,iLat) = 0.5 * &
            (VarWest(2:nLons+1) + VarEast(2:nLons+1) - &
            VarWest(1:nLons)    - VarEast(1:nLons)) * InvdLon

       ! Rusanov/Lax-Friedrichs diffusive term
       DiffFlux = 0.5 * max(&
            cMax_GDB(0:nLons,   iLat, iAlt, iEast_, iBlock), &
            cMax_GDB(1:nLons+1, iLat, iAlt, iEast_, iBlock)) &
            * (VarEast - VarWest)

       DiffVar(:,iLat) = &
            (DiffFlux(2:nLons+1) - DiffFlux(1:nLons)) * InvdLon

    end do

    if(UseTopography) &
         GradVar = GradVar - dVarDAlt_C * dAltDLon_CB(:,:,iAlt,iBlock)

  end subroutine calc_rusanov_lons

end subroutine advance_horizontal

!=============================================================================

subroutine calc_facevalues_lats(iLon, iAlt, iBlock, Var, VarLeft, VarRight)

  use ModSizeGITM, only: nLats
  use ModGITM, only: dLatDist_FB, InvDLatDist_FB
  use ModLimiterGitm

  implicit none
  
  integer, intent(in) :: iLon, iAlt, iBlock

  real, intent(in)    :: Var(-1:nLats+2)
  real, intent(out)   :: VarLeft(1:nLats+1), VarRight(1:nLats+1)

  real :: dVarUp, dVarDown, dVarLimited(0:nLats+1)

  real, parameter :: Factor1=0.6250000 ! 15/24
  real, parameter :: Factor2=0.0416667 !  1/24
  real :: h

  integer :: i
  !---------------------------------------------------------------------------

  i = 0

  h = InvDLatDist_FB(iLon,i+1,iAlt,iBlock)*2
  dVarUp   = h*(Factor1*(Var(i+1)-Var(i)  ) - Factor2*(Var(i+2)-Var(i-1)))
  dVarDown = (Var(i)   - Var(i-1)) * InvDLatDist_FB(iLon,i  ,iAlt,iBlock)
  dVarLimited(i)= Limiter_mc(dVarUp, dVarDown)

  do i=1,nLats

     h = InvDLatDist_FB(iLon,i+1,iAlt,iBlock)*2
     dVarUp   = h*(Factor1*(Var(i+1)-Var(i)  ) - Factor2*(Var(i+2)-Var(i-1)))
     h = InvDLatDist_FB(iLon,i  ,iAlt,iBlock)*2
     dVarDown = h*(Factor1*(Var(i)  -Var(i-1)) - Factor2*(Var(i+1)-Var(i-2)))

!     dVarUp   = (Var(i+1) - Var(i))   * InvDLatDist_FB(iLon,i+1,iAlt,iBlock)
!     dVarDown = (Var(i)   - Var(i-1)) * InvDLatDist_FB(iLon,i  ,iAlt,iBlock)

     dVarLimited(i)= Limiter_mc(dVarUp, dVarDown)
  end do

  i = nLats+1
  dVarUp   = (Var(i+1) - Var(i))   * InvDLatDist_FB(iLon,i+1,iAlt,iBlock)
  h = InvDLatDist_FB(iLon,i  ,iAlt,iBlock)*2
  dVarDown = h*(Factor1*(Var(i)  -Var(i-1)) - Factor2*(Var(i+1)-Var(i-2)))

  dVarLimited(i)= Limiter_mc(dVarUp, dVarDown)

  do i=1,nLats+1
     VarLeft(i) =Var(i-1)+0.5*dVarLimited(i-1)*dLatDist_FB(iLon,i,iAlt,iBlock)
     VarRight(i)=Var(i)  -0.5*dVarLimited(i)  *dLatDist_FB(iLon,i,iAlt,iBlock)
  end do

end subroutine calc_facevalues_lats

!=============================================================================

subroutine calc_facevalues_lons(iLat, iAlt, iBlock, Var, VarLeft, VarRight)

  use ModSizeGITM, only: nLons
  use ModGITM, only: dLonDist_FB, InvDLonDist_FB
  use ModLimiterGitm

  implicit none
  
  real, intent(in)    :: Var(-1:nLons+2)
  integer, intent(in) :: iLat,iAlt,iBlock
  real, intent(out)   :: VarLeft(1:nLons+1), VarRight(1:nLons+1)

  real :: dVarUp, dVarDown, dVarLimited(0:nLons+1)
  real, parameter :: Factor1=0.6250000 ! 15/24
  real, parameter :: Factor2=0.0416667 !  1/24
  real :: h

  integer :: i

  i = 0

  h  = InvDLonDist_FB(i+1,iLat,iAlt,iBlock)*2
  dVarUp   = h*(Factor1*(Var(i+1)-Var(i)  ) - Factor2*(Var(i+2)-Var(i-1)))
  dVarDown = (Var(i)   - Var(i-1))*InvDLonDist_FB(i  ,iLat,iAlt,iBlock)
  dVarLimited(i)= Limiter_mc(dVarUp, dVarDown)

  do i=1,nLons
     h  = InvDLonDist_FB(i+1,iLat,iAlt,iBlock)*2
     dVarUp   = h*(Factor1*(Var(i+1)-Var(i)  ) - Factor2*(Var(i+2)-Var(i-1)))
     h  = InvDLonDist_FB(i  ,iLat,iAlt,iBlock)*2
     dVarDown = h*(Factor1*(Var(i)  -Var(i-1)) - Factor2*(Var(i+1)-Var(i-2)))
!     dVarUp   = (Var(i+1) - Var(i))  *InvDLonDist_FB(i+1,iLat,iAlt,iBlock)
!     dVarDown = (Var(i)   - Var(i-1))*InvDLonDist_FB(i  ,iLat,iAlt,iBlock)
     dVarLimited(i)= Limiter_mc(dVarUp, dVarDown)
  end do

  i = nLons+1

  dVarUp   = (Var(i+1) - Var(i))  *InvDLonDist_FB(i+1,iLat,iAlt,iBlock)
  h  = InvDLonDist_FB(i  ,iLat,iAlt,iBlock)*2
  dVarDown = h*(Factor1*(Var(i)  -Var(i-1)) - Factor2*(Var(i+1)-Var(i-2)))
  dVarLimited(i)= Limiter_mc(dVarUp, dVarDown)

  do i=1,nLons+1
     VarLeft(i) =Var(i-1)+0.5*dVarLimited(i-1)*dLonDist_FB(i,iLat,iAlt,iBlock)
     VarRight(i)=Var(i)  -0.5*dVarLimited(i)  *dLonDist_FB(i,iLat,iAlt,iBlock)
  end do

end subroutine calc_facevalues_lons


