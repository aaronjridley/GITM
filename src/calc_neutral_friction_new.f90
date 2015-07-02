!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine calc_neutral_friction(oVel, EddyCoef_1d, NDensity_1d, NDensityS_1d, &
                                 GradLogCon, EddyCoefRatio_1d, Temp, Gravity_1d )

  use ModGITM
  use ModSources
  use ModPlanet, only: Diff0, DiffExp, IsEarth
  use ModInputs, only: UseNeutralFriction, UseBoquehoAndBlelly, UseEddyCorrection

  implicit none

  real,intent(inout) :: oVel(1:nAlts,1:nSpecies)
  real,intent(in) :: EddyCoef_1d(1:nAlts)
  real,intent(in) :: NDensity_1d(1:nAlts)
  real,intent(in) :: NDensityS_1d(1:nAlts,1:nSpecies)
  real,intent(in) :: GradLogCon(1:nAlts,1:nSpecies)
  real,intent(inout) :: EddyCoefRatio_1d(1:nAlts,1:nSpecies)
  real,intent(in) :: Temp(1:nAlts)
  real,intent(in) :: Gravity_1d(1:nAlts)

  integer :: iSpecies, jSpecies
  real :: CoefMatrix(nSpecies, nSpecies), kTOverM
  real :: Matrix(nSpecies, nSpecies)
  real :: Vel(nSpecies), Parity, Fraction=1.0, N0 = 1.28E+19
  integer :: iPivot(nSpecies)

! Added by Jared 11-29-2007
  real :: TempDij
  real :: InvDij(nSpecies)
  real :: Dij(nSpecies)
  real :: denscale
  real :: mscale
  real :: mms
  real :: mmwos(nSpecies)

  real :: EddyDiffCorrection(nSpecies), BulkWind, MassDensity, r

  integer :: iAlt

  call report("calc_neutral_friction",4)

  if (.not.UseNeutralFriction) return

    EddyCoefRatio_1d(1:nAlts,1:nSpecies) = 0.0

!    write(*,*) 'Now in Calc_Neutral_Friction'
  do iAlt = 1, nAlts

!    write(*,*) 'iAlt = ', iAlt
     Vel = oVel(iAlt,1:nSpecies)
     CoefMatrix = 0.0

     mms = 0.0
     mmwos = 0.0
     InvDij = 0.0

     EddyDiffCorrection(1:nSpecies) = 0.0

     do iSpecies = 1, nSpecies

        InvDij(iSpecies) = 0.0
        kTOverM = Boltzmanns_Constant * Temp(iAlt) / Mass(iSpecies)
        denscale = 1.0/NDensity_1d(iAlt) 

        do jSpecies = 1, nSpecies

           if (jSpecies == iSpecies) cycle
           
           ! \
           ! Please note that TempDij is the Dij binary coefficients
           ! Based upon the formulation by Banks and Kokarts.
           ! These coefficients demand that 
           ! (1) NDensity be in cm^-3 (hence the 1.0e-06) factor below
           ! (2) Additionally, the Dij's are in cm^2/s, thus the 1.0e-04 factor
           !         TempDij = (1.0e-04)*&   ! Scales the Dij from cm^2/s -> m^2/s
           !           (   Diff0(iSpecies,jSpecies)*
           !           ( Temp(iAlt)**DiffExp(iSpecies,jSpecies) )   ) / &
           !           (   NDensity_1d(iAlt)*(1.0e-06) )     ! Converts to #/cm^-3

           if (IsEarth) then 

              ! TempDij is the old Dij coefficient according to Ridley [2006]

              TempDij = Diff0(iSpecies, jSpecies)*&
                   (N0/NDensity_1d(iAlt))*&
                   (Temp(iAlt) / Temp(1))** DiffExp(iSpecies, jSpecies) 

              if (UseBoquehoAndBlelly) then

                 CoefMatrix(iSpecies, jSpecies) = &
                      kTOverM * denscale * NDensityS_1d(iAlt, jSpecies) / &
                      ( TempDij + EddyCoef_1d(iAlt)  )

                 InvDij(iSpecies) = InvDij(iSpecies) + &
                      denscale * NDensityS_1d(iAlt, jSpecies) / &
                      ( TempDij + EddyCoef_1d(iAlt)  )
              else

                 CoefMatrix(iSpecies, jSpecies) = &
                      kTOverM * denscale * NDensityS_1d(iAlt, jSpecies) / &
                      ( TempDij )

                 InvDij(iSpecies) = InvDij(iSpecies) + &
                      denscale * NDensityS_1d(iAlt, jSpecies) / &
                      ( TempDij )

              endif   !! UseBoqueho And Blelly

           else   !! Not Earth

              TempDij = (1.0e-04) * &       ! Scales the Dij from cm^2/s -> m^2/s
                   (Diff0(iSpecies,jSpecies)*&
                   (Temp(iAlt)**DiffExp(iSpecies,jSpecies) )   ) / &
                   ( NDensity_1d(iAlt)*(1.0e-06) )     ! Converts to #/cm^-3

              if ( UseBoquehoAndBlelly) then
                 CoefMatrix(iSpecies, jSpecies) = &
                      kTOverM * denscale * &
                      NDensityS_1d(iAlt, jSpecies) / &
                      ( TempDij + EddyCoef_1d(iAlt))
                 InvDij(iSpecies) = InvDij(iSpecies) + &
                      denscale*NDensityS_1d(iAlt, jSpecies)/ &
                      ( TempDij + EddyCoef_1d(iAlt) )
              else
                 CoefMatrix(iSpecies, jSpecies) = &
                      kTOverM * denscale * NDensityS_1d(iAlt, jSpecies) / TempDij
                 InvDij(iSpecies) = InvDij(iSpecies) + &
                      denscale*NDensityS_1d(iAlt, jSpecies)/ &
                      ( TempDij )
              endif

              if (UseBoquehoAndBlelly) then
                 EddyDiffCorrection(iSpecies) = &
                      EddyDiffCorrection(iSpecies) + &
                      Dt*InvDij(jSpecies)*EddyCoef_1d(iAlt)*&
                      GradLogCon(iAlt,jSpecies)
              else
                 EddyDiffCorrection(iSpecies) = &
                      EddyDiffCorrection(iSpecies) + &
                      Dt*CoefMatrix(iSpecies,jSpecies)* &
                      EddyCoef_1d(iAlt)*GradLogCon(iAlt,jSpecies)
              endif
              
           endif   ! Check if IsEarth
           
        enddo  ! End DO over jSpecies

        if (UseBoquehoAndBlelly) then

           EddyCoefRatio_1d(iAlt,iSpecies) =  &
                Dt*(EddyCoef_1d(iAlt)*InvDij(iSpecies))*GradLogCon(iAlt,iSpecies)

           if(UseEddyCorrection) then
              EddyCoefRatio_1d(iAlt,iSpecies) =  &
                   EddyCoefRatio_1d(iAlt,iSpecies) - EddyDiffCorrection(iSpecies)
           endif

        else

           EddyCoefRatio_1d(iAlt,iSpecies) =  &
                -1.0*Dt*(Boltzmanns_Constant*Temp(iAlt)/Mass(iSpecies) ) * &
                EddyCoef_1d(iAlt)*InvDij(iSpecies)*GradLogCon(iAlt,iSpecies) 

           if(UseEddyCorrection) then
              EddyCoefRatio_1d(iAlt,iSpecies) =  &
                   EddyCoefRatio_1d(iAlt,iSpecies) + EddyDiffCorrection(iSpecies)
           endif

        endif

     enddo  !End DO Over iSpecies

     Vel(1:nSpecies) = Vel(1:nSpecies) + EddyCoefRatio_1d(iAlt,1:nSpecies)

     Fraction = 1.4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CoefMatrix = Fraction * CoefMatrix

     Matrix = -dt*CoefMatrix

     do iSpecies = 1, nSpecies
        Matrix(iSpecies,iSpecies) = &
             1.0 + dt*(sum(CoefMatrix(iSpecies,:)))
     enddo

     BulkWind = 0.0
     MassDensity = 0.0
     do iSpecies = 1, nSpecies
        BulkWind = BulkWind + &
             vel(iSpecies) * mass(iSpecies) * NDensityS_1d(iAlt,iSpecies)
        MassDensity = MassDensity + mass(iSpecies)*NDensityS_1d(iAlt,iSpecies)
     enddo
     BulkWind = BulkWind/MassDensity

     do iSpecies = 1, nSpecies
        r = Matrix(iSpecies,iSpecies)/sum(Matrix(iSpecies,:))-1.0
!        write(*,*) iAlt, iSpecies, vel(iSpecies), r, BulkWind
        if (r > 1.0) vel(iSpecies) = BulkWind
     enddo

     call ludcmp(Matrix, nSpecies, nSpecies, iPivot, Parity)
     call lubksb(Matrix, nSpecies, nSpecies, iPivot, Vel)

     oVel(iAlt, 1:nSpecies) = Vel(1:nSpecies) !+ EddyCoefRatio_1d(iAlt,1:nSpecies) 

  enddo

end subroutine calc_neutral_friction
