subroutine get_msis_temperature(lon, lat, alt, t, h)

  use ModTime
  use ModInputs
  use ModPlanet
  use ModGITM

  implicit none

  real, intent(in) :: lon, lat, alt
  real, intent(out) :: t, h

!! At Titan the only major species are CH4, N2, and H2 for the background (H2 is still pretty minor)
  real :: nN2, nCH4, nH2, m, r, g
  integer :: i
  real :: TempThermo, TempBase
  real :: Width, AltMidPoint

  i = 1
  do while ( (alt >= newalt(i)) .and.  (i <= nAlts+2) ) 
     i = i + 1
  enddo
  i = i - 1

  TempBase = 175.0
  TempThermo = 150.0

  AltMidPoint = 900.0e+03
  Width = 150.0e+03
  t = TempBase - 0.5*(TempBase - TempThermo)*&
      (1.0 + TANH( (alt - AltMidPoint)/Width) )

  m = Mass(iN2_)
  r = RBody + alt
  g = Gravitational_Constant * (RBody/r) ** 2
  h = Boltzmanns_Constant * t / (m*g)

end subroutine get_msis_temperature


subroutine init_msis

  use ModPlanet
  use ModGITM
  use ModEUV
  use ModInputs, only: iDebugLevel, iInputUnit_, EddyDiffusionCoef

  implicit none

  integer :: iBlock

  integer :: iiLon,iiLat,iiAlt
  integer :: iLon,iLat,iAlt, iSpecies, iIon, jSpecies

! Variables for the Asymmetric Wind Profile at the
! Lower Boundary

  real :: PeakLat
  real :: MinVelocity
  real :: RiseBeginLat, RiseEndLat
  real :: StandardT, RequiredT
  real :: RatioT
  real :: DelLat
  real :: NewLat
  real :: Freq

!! CH4 Mixing Ratio Fixes
  real :: DesiredCH4
  real :: TestCH4
  real :: CH4Scaling

!! CH4 Mixing Ratio Fixes
  real :: DesiredH2
  real :: TestH2
  real :: H2Scaling

!! CH4 Mixing Ratio Fixes
  real :: DesiredAr
  real :: TestAr
  real :: ArScaling

!! HCN Mixing Ratio Fixes
  real :: DesiredHCN
  real :: TestHCN
  real :: HCNScaling

!! CH4 Mixing Ratio Fixes
  real :: DesiredH
  real :: TestH
  real :: HScaling

  real :: GlobalScaling 
  real :: CIRSDensityValue

  real :: MinorMixing(1:nSpeciesTotal) 
  real :: TotalMixing


!!! Setting up a hyberbolic Tangent temperature field
  real :: TempBase
  real :: TempThermo
  real :: AltMidPoint, Width

!!! Used to create our initial density profiles
!!! Note that this method is not perfect
!!! Note also that it doesn't need to be.

  real :: InvHs, MeanTemp, MeanGravity, MeanMass
  real :: InvHa, MeanLambda(1:nSpecies), MeanDS(1:nSpecies)

  real :: TempDij, InvDij, MeanKE

  real :: Centrifugal, CentrifugalAccel(-1:nAlts+2)
  real :: EffectiveGravity(-1:nAlts+2)

  real :: HorizArg

  real :: DensityScaleFactor(-1:nLons+2,-1:nLats+2)
  real :: MeanMass0

!!! This is our Anchor Density Value

  !CIRSDensityValue = 1.20e+20
  CIRSDensityValue = 9.00e+19

  MinorMixing(1:nSpeciesTotal) = 1.0e-20

  MinorMixing(iAr_) = 3.30e-05
  MinorMixing(iCH4_) = 0.0140
  MinorMixing(iHCN_) = 1.0e-06

  MinorMixing(iH2_) = 0.0020
  MinorMixing(iH_) = 5.0e-07

  MinorMixing(iN4S_) = 1.0e-10
  MinorMixing(iC2H4_) = 1.0e-10

  MinorMixing(i13CH4_) = (1.0/91.0)*MinorMixing(iCH4_)


!! Mixing Ratio of the Major Background Species

  if (nSpecies .gt. 1) then
    MinorMixing(iN2_) = 1.0 
    do iSpecies = 2, nSpecies
      MinorMixing(iN2_) = MinorMixing(iN2_) - MinorMixing(iSpecies)
    enddo 
  else
      MinorMixing(1) = 1.0
  endif
  MinorMixing(i15N2_) = (2.0/167.7)*MinorMixing(iN2_)


  ! Just do this once, since nBlocks = 1 most of the time
  iBlock = 1
  MeanMass0 = 0.0
  do iSpecies = 1, nSpecies
      MeanMass0 = MeanMass0 + &
          Mass(iSpecies)*MinorMixing(iSpecies)
  enddo 

  TempBase = 175.0

!  call calc_densityscaling(iBlock, TempBase, MeanMass0, DensityScaleFactor)
!!! First Stipulate a Background Temperature
DensityScaleFactor(-1:nLons+2,-1:nLats+2) = 1.0

TempBase = 175.0
TempThermo = 150.0

!TempThermo = 152.0


!write(*,*) '=========================================================='
!write(*,*) ' NOW IN INIT_MSIS ===='

AltMidPoint = 900.0e+03
Width = 150.0e+03
do iBlock = 1, nBlocks
   do iLat = -1, nLats + 2
      do iLon = -1, nLons + 2
        do iAlt = -1, nAlts + 2
          Temperature(iLon,iLat,iAlt,iBlock) = &
               TempBase - 0.5*(TempBase - TempThermo)*&
                (1.0 + TANH( (Altitude_GB(iLon,iLat,iAlt,iBlock) - AltMidPoint)/Width) )
          InTemp(iAlt) = Temperature(iLon,iLat,iAlt,iBlock)
          !write(*,*) 'iAlt, Alt, RadDist, Temperature = ',&
          !            iAlt, Altitude_GB(iLon,iLat,iAlt,iBlock), &
          !            RadialDistance_GB(iLon,iLat,iAlt,iBlock), &
          !                  Temperature(iLon,iLat,iAlt,iBlock) 
        enddo 


       ! do iAlt = -1, nAlts + 2
       !   write(*,*) 'iAlt, Gravity = ',iAlt, Gravity_GB(iLon,iLat,iAlt,iBlock) 
       ! enddo 

          Temperature(iLon,iLat,0,iBlock) = Temperature(iLon,iLat,1,iBlock)
         Temperature(iLon,iLat,-1,iBlock) = Temperature(iLon,iLat,1,iBlock)


           eTemperature(iLon,iLat,-1:nAlts+2,iBlock) =  &
                Temperature(iLon,iLat,-1:nAlts+2,iBlock)

!           eTemperature(iLon,iLat,-1:nAlts+2,iBlock) =  &
!                IneTemp(-1:nAlts+2)

           ITemperature(iLon,iLat,-1:nAlts+2,iBlock) =  &
                Temperature(iLon,iLat,-1:nAlts+2,iBlock)

        enddo !iLon = -1, nLons + 2
     enddo !iLat = -1, nLats + 2
  enddo !iBlock = 1, nBlocks


!!! Next, Integrate upward from the -1 Cell
  do iBlock = 1, nBlocks
   do iLat = -1, nLats + 2
     do iLon = -1, nLons + 2

       !!!!! Separately Calculate the -1 Cell Values
       iAlt = -1

       MeanMajorMass(iLon,iLat,iAlt) = 0.0
       do iSpecies = 1, nSpecies
          MeanMajorMass(iLon,iLat,iAlt) = &  
          MeanMajorMass(iLon,iLat,iAlt) + &  
            Mass(iSpecies)*MinorMixing(iSpecies)
       enddo 
       !write(*,*) 'iAlt, MeanMass = ', MeanMajorMass(iLat,iLon,iAlt)/AMU

       NDensity(iLon,iLat,iAlt,iBlock) = DensityScaleFactor(iLon,iLat)*CIRSDensityValue
       NDensityS(iLon,iLat,iAlt,1:nSpeciesTotal,iBlock) = &
          DensityScaleFactor(iLon,iLat)*CIRSDensityValue*MinorMixing(1:nSpeciesTotal)

       MeanKE = EddyDiffusionCoef*&
          ( NDensity(iLon,iLat,-1,iBlock)/NDensity(iLon,iLat,iAlt,iBlock) )**0.50

      !!!! Next Calculate the Molecular Diffusion Coefficient
       do iSpecies = 1, nSpecies
          InvDij = 0.0
          do jSpecies = 1, nSpecies
              if (iSpecies .eq. jSpecies) cycle
                    TempDij = (1.0e+02)*&                        
                     (   Diff0(iSpecies,jSpecies)*&
                        ( Temperature(iLon,iLat,iAlt,iBlock)**DiffExp(iSpecies,jSpecies) )   ) / &
                         NDensity(iLon,iLat,iAlt,iBlock)    

                    InvDij = InvDij + &
                        NDensityS(iLon,iLat,iAlt,jSpecies,iBlock)/&
                       ( NDensity(iLon,iLat,iAlt,iBlock)*TempDij)
          enddo !jSpecies = 1, nSpecies
          MeanDs(iSpecies) = 1.0/InvDij
          MeanLambda(iSpecies) = MeanKE/MeanDs(iSpecies)
       enddo !iSpecies = 1, nSpecies

!!! Add the Centrifugal Forces to the Gravity

   do iAlt = -1, nAlts + 2
      EffectiveGravity(iAlt) = Gravity_GB(iLon,iLat,iAlt,iBlock) 
   enddo 

   ! write(*,*) 'INIT_MSIS:  NDENSITY LOOP:'
    do iAlt = 0, nAlts + 2

      MeanGravity = -0.5*(EffectiveGravity(iAlt) + EffectiveGravity(iAlt-1))
         MeanTemp = 0.5*(Temperature(iLon,iLat,iAlt,iBlock) + Temperature(iLon,iLat,iAlt-1,iBlock))

      MeanMass = 0.0
      do iSpecies = 1, nSpecies
           MeanMass = MeanMass + &
             NDensityS(iLon,iLat,iAlt-1,iSpecies,iBlock)*Mass(iSpecies)/&
              NDensity(iLon,iLat,iAlt-1,iBlock)
      enddo 

      InvHa = MeanMass*MeanGravity/(MeanTemp*Boltzmanns_Constant)

!!! MeanMassCalculation
      NDensity(iLon,iLat,iAlt,iBlock) = 0.0

      do iSpecies = 1, nSpecies
            NDensityS(iLon,iLat,iAlt,iSpecies,iBlock) = &
               NDensityS(iLon,iLat,iAlt-1,iSpecies,iBlock)*&
              (Temperature(iLon,iLat,iAlt-1,iBlock)/Temperature(iLon,iLat,iAlt,iBlock))*&
               exp( -1.0*(Altitude_GB(iLon,iLat,iAlt  ,iBlock) - &
                          Altitude_GB(iLon,iLat,iAlt-1,iBlock))*InvHa)
 
           NDensity(iLon,iLat,iAlt,iBlock) = &
           NDensity(iLon,iLat,iAlt,iBlock) + NDensityS(iLon,iLat,iAlt,iSpecies,iBlock)
      enddo    !iSpecies

      MeanKE = EddyDiffusionCoef*&
         ( NDensity(iLon,iLat,-1,iBlock)/NDensity(iLon,iLat,iAlt,iBlock) )**0.50

      do iSpecies = 1, nSpecies
         InvDij = 0.0

         do jSpecies = 1, nSpecies
             if (iSpecies .eq. jSpecies) cycle
              
              TempDij = (1.0e+02)*&                        
                    (   Diff0(iSpecies,jSpecies)*&
                       ( Temperature(iLon,iLat,iAlt,iBlock)**DiffExp(iSpecies,jSpecies) )   ) / &
                        NDensity(iLon,iLat,iAlt,iBlock)    

               InvDij = InvDij + &
                        NDensityS(iLon,iLat,iAlt,jSpecies,iBlock)/&
                       ( NDensity(iLon,iLat,iAlt,iBlock)*TempDij)

         enddo !jSpecies = 1, nSpecies

         MeanDs(iSpecies) = 1.0/InvDij

         MeanLambda(iSpecies) = MeanKE/MeanDs(iSpecies)

      enddo !iSpecies = 1, nSpecies


               MeanMajorMass(iLon,iLat,iAlt) = 0.0
           do iSpecies = 1, nSpecies
               MeanMajorMass(iLon,iLat,iAlt) = &
               MeanMajorMass(iLon,iLat,iAlt) + &
                 NDensityS(iLon,iLat,iAlt,iSpecies,iBlock)*Mass(iSpecies)/&
                 NDensity(iLon,iLat,iAlt,iBlock)
           enddo 



       !write(*,*) 'iAlt, Altitude_GB',iAlt, Altitude_GB(iLon,iLat,iAlt,iBlock)
     enddo  !iAlt
   enddo   ! Lats
 enddo    ! Lons
enddo    ! Blocks

!!! Set all of the very minor species to fixed mixing ratios (initially)
do iBlock = 1, nBlocks 
 do iLon = -1, nLons + 2
   do iLat = -1, nLats + 2
     do iAlt = -1, nAlts + 2
       do iSpecies = nSpecies+1, nSpeciesTotal
              NDensityS(iLon,iLat,iAlt,iSpecies,iBlock) = &
               max(1.0e+01, NDensity(iLon,iLat,iAlt,iBlock) * MinorMixing(iSpecies))
       enddo !iSpecies 
     enddo  !iAlt
   enddo   ! Lats
 enddo    ! Lons
enddo    ! Blocks


  do iBlock = 1, nBlocks
     Rho(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock) = &
          MeanMajorMass(-1:nLons+2,-1:nLats+2,-1:nAlts+2)* &
          NDensity(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock)

  enddo

!!!!!! Input a Frozen HCN as our initial Settings
!
!
do iBlock = 1, nBlocks
    do iLat = -1, nLats + 2
      do iLon = -1, nLons + 2
         do iAlt = -1, nAlts + 2

           !TempUnit(iLon,iLat,iAlt) = MeanMajorMass(iLon,iLat,iAlt)/Boltzmanns_Constant
           TempUnit(iLon,iLat,iAlt) = Mass(iN2_)/Boltzmanns_Constant
!!! Scale Temps to the GITM Format
           Temperature(iLon,iLat,iAlt,iBlock) = &
           Temperature(iLon,iLat,iAlt,iBlock)/TempUnit(iLon,iLat,iAlt)
!!! Reduce HCN by factor of 10 from VDH's values


            do iIon = 1, nIons-1 
            IDensityS(iLon,iLat,iAlt,iIon,iBlock) = &
                     1.0
            enddo 
            IDensityS(iLon,iLat,iAlt,ie_,iBlock) = &
                 nIons-1
           
   
         enddo 
      enddo 
    enddo 
enddo 

!iBlock = 1
!        do iAlt = -1, nAlts + 2 
!          write(*,*) 'iAlt, Alt, RadDist, Temperature = ',&
!                      iAlt, Altitude_GB(1,1,iAlt,iBlock), &
!                      RadialDistance_GB(1,1,iAlt,iBlock), &
!                            Temperature(1,1,iAlt,iBlock) 
!        enddo 
!write(*,*) '---------- END OF INIT_MSIS CHECK '

end subroutine init_msis


subroutine calc_densityscaling(iBlock,T0, MeanMass, DenScaling)
   use ModSizeGITM, only : nLats, nLons
   use ModGITM, only : Latitude, Longitude, RadialDistance_GB
   use ModInputs, only : EquatorialRadius, PolarRadius
   use ModConstants
   use ModPlanet, only : OMEGABody
 
   implicit none

   ! For now, we assume that the Mean Mass is uniform at the lower boundary
   ! Otherwise, this approximation does not work and we must numericall integrate
   integer, intent(in) :: iBlock
   real, intent(out) :: DenScaling(-1:nLons+2,-1:nLats+2)
   real, intent(in) :: MeanMass
   real, intent(in) :: T0

   ! Major Orbital Parameters for the Titan-Saturn Tidal Equilibrium 
   real :: a_titan, e_titan, Mass_Saturn 
   real :: omega_titan
   real :: GravConstant
   real :: ExponentArg
   integer :: iAlt, iLat, iLon, iDir

   a_titan = 1.22183e+09   ! Titan orbital semi-major axis (in m) 
   e_titan = 0.0292        ! Titan orbital eccentricity (e)
   omega_titan = OMEGABody ! Titan orbital angular frequency (rads/s)
   Mass_Saturn = 5.685e+26 ! Mass of Saturn
   GravConstant = 6.67259e-11  ! Graviational constant in SI units (N m^2/kg^2)
   !Note:  Titan is phase-locked so that it's rotation rate and orbit rate are the same


   do iLat = -1, nLats+2
     do iLon = -1, nLons+2
         ExponentArg = &
            1.5*(GravConstant*Mass_Saturn/a_titan)*&
            (   (RadialDistance_GB(iLon,iLat,0,iBlock)/a_titan)**2.0  )*&
            ( MeanMass/(Boltzmanns_Constant*T0) )

         DenScaling(iLon,iLat) = &
             exp(-1.0*ExponentArg*( &
                    sin(Longitude(iLon,iBlock))**2.0 + &
                    (cos(Longitude(iLon,iBlock))**2.0)*&
                    (sin(Latitude(iLat,iBlock))**2.0) ) ) 
         
     enddo !iLon = -1, nLons+2
   enddo !iLat = -1, nLats+2
 

end subroutine calc_densityscaling

subroutine msis_bcs(iJulianDay,UTime,Alt,Lat,Lon,Lst, &
             F107A,F107,AP,LogNS, Temp, LogRho)

  write(*,*) "You can not use MSIS with any planet except Earth!!!"
  write(*,*) "If you ARE running Earth, then make the code again, using"
  write(*,*) "configure Earth ; make"
  call stop_gitm("I can not continue...")

end subroutine msis_bcs
