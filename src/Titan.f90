!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine fill_photo(photoion, photoabs, photodis)

  use ModPlanet
  use ModEUV
  use ModInputs

  implicit none

  real, intent(out) :: photoion(Num_WaveLengths_High, nIons-1)
  real, intent(out) :: photoabs(Num_WaveLengths_High, nSpecies)
  real, intent(out) :: photodis(Num_WaveLengths_High, nSpeciesTotal)

  integer :: iSpecies, iWave, NWH, i, iIon

  NWH = Num_WaveLengths_High

  photoabs(:,:) = 0.0
  photodis(:,:) = 0.0
  photoion(:,:) = 0.0

  PhotoAbs = 0.0
  PhotoIon = 0.0
  PhotoDis = 0.0

  photoabs(1:NWH,iN2_) = PhotoAbs_N2(1:NWH)
  photoabs(1:NWH,iCH4_) = PhotoAbs_CH4(1:NWH)

!  photoabs(1:NWH,iH2_) = PhotoAbs_H2(1:NWH)

!! Photo-dissociation cross sections (N2 and CH4 Only for Now)

   photodis(1:NWH,iN2_) = QuantumYield_N2_N4S(1:NWH) * &
                          PhotoAbs_N2(1:NWH)

  photodis(1:NWH,iCH4_)=  (   QuantumYield_CH4_CH3(1:NWH)  + &
                              QuantumYield_CH4_3CH2(1:NWH) + &
                              QuantumYield_CH4_1CH2(1:NWH) + &
                              QuantumYield_CH4_CH(1:NWH)  ) *  &
                              PhotoAbs_CH4(1:NWH)

  photodis(1:NWH,iH2_)= PhotoAbs_H2(1:NWH)

!\
! Photoproduction of N(4S) from N2:---------------------------------
!/
  photodis(1:NWH,iN4S_)= QuantumYield_N2_N4S(1:NWH) * &
                         PhotoAbs_N2(1:NWH)

!/
!\
! Photoproduction of CH3, 3CH2, 1CH2, and CH from CH4:--------------
!
  
  photodis(1:NWH,iCH3_)=  QuantumYield_CH4_CH3(1:NWH) * &
                          PhotoAbs_CH4(1:NWH)

  photodis(1:NWH,i3CH2_)= QuantumYield_CH4_3CH2(1:NWH) * &
                          PhotoAbs_CH4(1:NWH)

  photodis(1:NWH,i1CH2_)= QuantumYield_CH4_1CH2(1:NWH) * &
                          PhotoAbs_CH4(1:NWH)

  photodis(1:NWH,iCH_)= QuantumYield_CH4_CH(1:NWH) * &
                          PhotoAbs_CH4(1:NWH)
!
  photoion(1:NWH,iNP_)=  QuantumYield_N2_NPlus(1:NWH)  * &
                         PhotoAbs_N2(1:NWH) 

  photoion(1:NWH,iN2P_)= QuantumYield_N2_N2Plus(1:NWH) * &
                         PhotoAbs_N2(1:NWH) 

  photoion(1:NWH,iCH3P_)=   QuantumYield_CH4_CH3Plus(1:NWH)*  &
                            PhotoAbs_CH4(1:NWH)   

end subroutine fill_photo

subroutine calc_planet_sources(iBlock)

  use ModInputs
  use ModSources
  use ModGITM
  use ModTime
  use ModPlanet
  
  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iError, iDir, iLat, iLon

  LowAtmosRadRate = 0.0

  call calc_radcooling(iBlock)
  call calc_aerosols(iBlock)


  RadCooling(1:nLons,1:nLats,1:nAlts,iBlock) = &
    RadCoolingRate(1:nLons,1:nLats,1:nAlts,iBlock)/&
    (  TempUnit(1:nLons,1:nLats,1:nAlts) * &
       cp(1:nLons,1:nLats,1:nAlts,iBlock)*&
       rho(1:nLons,1:nLats,1:nAlts,iBlock)  )


  NOCooling = 0.0
  OCooling = 0.0

end subroutine calc_planet_sources

subroutine planet_limited_fluxes(iBlock)

  use ModInputs
  use ModSources
  use ModGITM
  use ModTime
  use ModPlanet

  implicit none
  
  integer, intent(in) :: iBlock

  integer :: iAlt,iSpecies, iLon,iLat


    do iLat = 1, nLats
      do iLon = 1, nLons
  
  do iAlt = -1,4

     VerticalVelocity(iLon,iLat,iAlt,iN2_,iBlock) = ((2575.0/3075.0)**2.0)*(2.5e+11)/NDensityS(iLon,iLat,iAlt,iN2_,iBlock)

     VerticalVelocity(iLon,iLat,iAlt,iCH4_,iBlock) = ((2575.0/3075.0)**2.0)*(2.5e+13)/NDensityS(iLon,iLat,iAlt,iCH4_,iBlock)

     VerticalVelocity(iLon,iLat,iAlt,iH2_,iBlock) = ((2575.0/3075.0)**2.0)*(3.0e+13)/NDensityS(iLon,iLat,iAlt,iH2_,iBlock)

     VerticalVelocity(iLon,iLat,iAlt,iHCN_,iBlock) = -1.0*((2575.0/3075.0)**2.0)*(2.0e+12)/NDensityS(iLon,iLat,iAlt,iHCN_,iBlock)

  enddo

      enddo 
    enddo 

end subroutine planet_limited_fluxes

!---------------------------------------------------------------------
! Initialize Heating Rates
!---------------------------------------------------------------------

subroutine init_heating_efficiency

  use ModPlanet
  use ModGITM, only: nLons, nLats, nAlts, nBlocks, Altitude_GB
  use ModEUV, only: HeatingEfficiency_CB, eHeatingEfficiency_CB
  use ModIoUnit, only: UnitTmp_

  implicit none

  integer :: iLon, iLat, iAlt
  real, Dimension(-1:nAlts+2) :: Heffalt
  real, Dimension(-1:nAlts+2) :: Heffin
  real, Dimension(-1:nAlts+2) :: HeffPH
  !------------------------------------------------------------------

  HeatingEfficiency_CB(:,:,:,1:nBlocks) = 0.05
  eHeatingEfficiency_CB(:,:,:,1:nBlocks) = 0.05

  open(UNIT = UnitTmp_,FILE = 'DataIn/Titan_Heff_500km_10km.txt',STATUS='OLD',ACTION='READ')
135 FORMAT(F6.1,1X,ES10.3,1X,ES10.3)

  do iAlt = -1,nAlts+2
     read(UnitTmp_,135) &
     Heffalt(iAlt),&
     Heffin(iAlt), &
     HeffPH(iAlt)
  enddo

  close(Unit = UnitTmp_)

  do iLon = 1, nLons
    do iLat = 1, nLats
      do iAlt = 1, nAlts

    HeatingEfficiency_CB(iLon,iLat,iAlt,1:nBlocks) = Heffin(iAlt)
    eHeatingEfficiency_CB(iLon,iLat,iAlt,1:nBlocks) = Heffin(iAlt)

      enddo 
    enddo 
  enddo


end subroutine init_heating_efficiency

!---------------------------------------------------------------------
! Calculate Eddy Diffusion Coefficient
!---------------------------------------------------------------------

subroutine calc_eddy_diffusion_coefficient(iBlock)

  use ModSizeGITM
  use ModGITM, only: NDensity
  use ModInputs, only: EddyDiffusionPressure0,EddyDiffusionPressure1, &
       EddyDiffusionCoef
  use ModSources, only: KappaEddyDiffusion

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon

  KappaEddyDiffusion=0.

! 10.0e+03 seemed too high

  do iAlt = -1, nAlts+2
     do iLat = 1, nLats
        do iLon = 1, nLons

          KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = EddyDiffusionCoef*&
                 ( NDensity(iLon,iLat,1,iBlock)/NDensity(iLon,iLat,iAlt,iBlock) )**0.65

        enddo
     enddo
  enddo

end subroutine calc_eddy_diffusion_coefficient

!! Aerosol Calculations BEGIN

subroutine init_aerosol

  use ModPlanet
  use ModConstants
  use ModIoUnit, only:  UnitTmp_

  real , Dimension(1:nAlts) :: AerosolDensity
  real , Dimension(1:nAlts) :: ClassA
  real , Dimension(1:nAlts) :: ClassB
  real , Dimension(1:nAlts) :: ClassC
  real , Dimension(1:nAlts) :: AeroAlt
  integer :: iLon,iLat,iAlt,iSpecies,iIon

  85 FORMAT(F7.2, 1X, ES10.3, 1X, ES10.3, 1X, ES10.3)
  open(UNIT = UnitTmp_,FILE = 'DataIn/AerosolDensity_Akiva_10km.txt',STATUS='OLD',ACTION='READ')

  do iAlt = 1,nAlts
     read(UnitTmp_,85) &
     AeroAlt(iAlt),&
     ClassA(iAlt), &
     ClassB(iAlt), &
     ClassC(iAlt)
  enddo

  close(UNIT = UnitTmp_)

  do iBlock = 1,nBlocks
   do iLon = 1,nLons
    do iLat = 1,nLats
     do iAlt = 1,nAlts

        AerosolClassA(iLon,iLat,iAlt,iBlock) = &
           ClassA(iAlt)

        AerosolClassB(iLon,iLat,iAlt,iBlock) = &
           ClassB(iAlt)

        AerosolClassC(iLon,iLat,iAlt,iBlock) = &
           ClassC(iAlt)

     enddo
    enddo
   enddo
  enddo 

!!! VERIFIED THAT THE AEROSOL VALUES ARE BEING READ CORRECLTY!!! JMB 9-17-2008

  end subroutine init_aerosol


  subroutine calc_aerosols(iBlock)


  use ModPlanet
  use ModConstants
  use ModGITM

  implicit none

  integer,intent(in):: iBlock 

  integer :: iLon,iLat,iAlt,iSpecies,iIon

  real :: TrappingRate           !  The cross-section Frequency for CH4
  real :: HeterogeneousRate       !  The cross-section Frequency for CH4

  real :: SigmaH                 !  The cross-section for Hydrogen Abstraction
  real :: NHaze                 !  The # of Haze Abstraction Sites available

  real , Dimension(1:nLons,1:nLats,1:nAlts) :: Vth
  real , Dimension(1:nLons,1:nLats,1:nAlts) :: TempLocal



    do iLon = 1, nLons
      do iLat = 1, nLats
        do iAlt = 1, nAlts 

            TrappingRate = 0.0

! This is based on Gas-trapping by Aerosols (Akiva Bar-Nun et al (2008))
!
! First AerosolClassA, ClassB, and ClassC are the # of particles/m^3 of each
! type of Aerosol Particle Group
!
! Class A:  Each as ~10-15 Polymeric Moecules
!      Polymeric Moecule => C8H8
!
! Class B:  Each as ~500 Embryos (ClassA)
!      Each ClassA Again consists of ~10 Polymeric Molecules 
!      Polymeric Moecule => C8H8  => 8 Carbon Atoms
!!
! Class C:  Each consists of ~1,100 Class B Particles
!           Each B has ~500 Embryos
!           Each Embryo ~ 10 Polymeric Molecules
!           Each Polymeric Molecule => 8 C atoms


! Efficiency of trapping CH4 is roughly 7.0e-02 (same as Xe)

! According to BAR-Nun
! Efficiencies Based upon Amorphous Ice Trapping

! Class A Loss Rates

                        TrappingRate =   8.0 * &                       ! # C Atoms/Polymeric Molecule
                                        20.0 * &                       ! # Polymeric Molecules/Embryo     
      AerosolClassA(iLon,iLat,iAlt,iBlock)                           ! # of Embryos/m^3

! Class B Loss Rates
                TrappingRate =  TrappingRate + &
                                         8.0 * &                       ! # C Atoms/Polymeric Molecule
                                        20.0 * &                       ! # Polymeric Molecules/Embryo     
                                       600.0 * &                       ! # Embryos /ClassB     
      AerosolClassB(iLon,iLat,iAlt,iBlock)                           ! # of ClassB/m^3

! Class C Loss Rates
                TrappingRate =  TrappingRate + &
                                         8.0 * &                       ! # C Atoms/Polymeric Molecule
                                        20.0 * &                       ! # Polymeric Molecules/Embryo     
                                       600.0 * &                       ! # Embryos /ClassB     
                                      1800.0 * &                       ! # ClassB /ClassC     
      AerosolClassC(iLon,iLat,iAlt,iBlock)                           ! # of ClassC / m^3

          AerosolTrappingLoss(iLon,iLat,iAlt,iCH4_,iBlock) = &
               TrappingRate*(7.0e-02)                                 ! CH4 Trapping Efficiency at Xenons (7.0e-02)

          AerosolTrappingLoss(iLon,iLat,iAlt,i13CH4_,iBlock) = &
               TrappingRate*(7.0e-02)*sqrt(Mass(i13CH4_)/Mass(iCH4_)) ! 13CH4 Trapping Efficiency at Xenons (7.0e-02) Scaled by Mass

          AerosolTrappingLoss(iLon,iLat,iAlt,iN2_,iBlock) = &
               TrappingRate*((1.0e-04)/70.0)                          ! N2 Trapping Efficiency at (1/70) Argons 

          AerosolTrappingLoss(iLon,iLat,iAlt,i15N2_,iBlock) = &
               TrappingRate*((1.0e-04)/70.0)*sqrt(Mass(i15N2_)/Mass(iN2_))   ! 15N2 Trapping Efficiency at Xenons (7.0e-02) Scaled by Mass

          AerosolTrappingLoss(iLon,iLat,iAlt,iAr_,iBlock) = &
               TrappingRate*(1.0e-04)                                 ! Ar Trapping Efficiency 

          AerosolTrappingLoss(iLon,iLat,iAlt,iH2_,iBlock) =  0.0
          AerosolTrappingLoss(iLon,iLat,iAlt,iHCN_,iBlock) =  0.0

          AerosolTrappingLoss(iLon,iLat,iAlt,1:nSpecies,iBlock) =  0.0
        enddo
      enddo
    enddo


!! Heterogeneous Processes from the Study by Lebonnois et al (2003)

Vth = 0.0
TempLocal = 0.0
SigmaH = 0.0
NHaze = 0.0

    do iLon = 1, nLons
      do iLat = 1, nLats
        do iAlt = 1, nAlts 

            HeterogeneousRate = 0.0

         TempLocal(iLon,iLat,iAlt) = &
             Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt)

         SigmaH = (1.0e-15)*(1.0e-04)*  &           ! (1.0e-04) converts from cm^2 -> m^2
                    exp(-1700.0/TempLocal(iLon,iLat,iAlt))                  

         Vth(iLon,iLat,iAlt) = &
               sqrt( (3.0*Boltzmanns_Constant/Mass(iH_) ) * TempLocal(iLon,iLat,iAlt) ) 

!! Class A Species
                NHaze = 4.0*PI*( (1.2e-09)**2)/(2.5e-20)     ! Bar-Nun (2008):  Radius of Particles A ~1.2 x 10^-9 m
                                                             ! Lebonnios (2003) uses a surface area of 2.5e-16 cm^2
            HeterogeneousRate = HeterogeneousRate + &
                (0.5)*NDensityS(iLon,iLat,iAlt,iH_,iBlock)*&
                Vth(iLon,iLat,iAlt)* &
                AerosolClassA(iLon,iLat,iAlt,iBlock)*&
                SigmaH * NHaze

! Class B Species (R ~ 1.5e-8 m):  Bar-Nun (2008)
                NHaze = 4.0*PI*( (1.5e-08)**2)/(2.5e-20)   

            HeterogeneousRate = HeterogeneousRate + &
                (0.5)*NDensityS(iLon,iLat,iAlt,iH_,iBlock)*&
                Vth(iLon,iLat,iAlt)* &
                AerosolClassB(iLon,iLat,iAlt,iBlock)*&
                SigmaH * NHaze

! Class C Species (R ~ 1.5e-8 m):  Bar-Nun (2008)
                NHaze = 4.0*PI*( (2.5e-06)**2)/(2.5e-20)   

            HeterogeneousRate = HeterogeneousRate + &
                (0.5)*NDensityS(iLon,iLat,iAlt,iH_,iBlock)*&
                Vth(iLon,iLat,iAlt)* &
                AerosolClassC(iLon,iLat,iAlt,iBlock)*&
                SigmaH * NHaze

!! Total Heterogeneous Production Rate of H2

           H2AerosolProduction(iLon,iLat,iAlt,iBlock) = HeterogeneousRate


        enddo
      enddo
    enddo

  end subroutine calc_aerosols

subroutine init_isochem

  use ModPlanet
  use ModGITM
  use ModInputs
  use ModIoUnit, only : UnitTmp_

  real , Dimension(1:nAlts) :: tempalt, IsoChem
  integer :: iiAlt, iiBlock, iiLon, iiLat

  open(UNIT = UnitTmp_, FILE = 'DataIn/IsoChem_10km.txt', STATUS='OLD',ACTION = 'READ')
   165 FORMAT(F8.3,1X,F6.2)

     do iiAlt = 1, nAlts
          read(UnitTmp_,165)   tempalt(iiAlt), &
                               IsoChem(iiAlt)
      enddo

    close(Unit = UnitTmp_)
!
    do iiBlock = 1, nBlocks
       do iiAlt = 1, nAlts
         do iiLon = 1, nLons
           do iiLat = 1, nLats

               IsotopeScaling(iiLon,iiLat,iiAlt,iiBlock) = IsoChem(iiAlt)

           enddo 
         enddo 
       enddo
    enddo


end subroutine init_isochem

subroutine init_magheat

  use ModPlanet
  use ModGITM
  use ModEUV
  use ModIoUnit, only : UnitTmp_

  implicit none

  integer :: iiLon,iiLat,iiAlt
  integer :: iLon,iLat,iAlt,iSpecies,iBlock,iIon

  real , Dimension(-1:nAlts + 2) :: eheatalt, eheat
  real , Dimension(-1:nAlts + 2,nSpeciesTotal) :: InMagDiss 
  real , Dimension(-1:nAlts + 2,nIons - 1) :: InMagIon


  eheatalt = 0.0
  eheat = 0.0
  InMagDiss = 0.0
  InMagIon = 0.0

   open(UNIT = UnitTmp_, FILE = 'DataIn/MagForcing_500km_10km.txt', STATUS='OLD',ACTION = 'READ')
   155 FORMAT(     F7.2, &
                   1X,ES10.3, 1X,ES10.3, 1X,ES10.3, 1X,ES10.3, 1X,ES10.3, &
                   1X,ES10.3, 1X,ES10.3, 1X,ES10.3, 1X,ES10.3 )

     do iiAlt = 1, nAlts
          read(UnitTmp_,155)   eheatalt(iiAlt), &
                        InMagDiss(iiAlt,iN4S_), &
                        InMagDiss(iiAlt,iCH_), &
                        InMagDiss(iiAlt,i3CH2_), &
                        InMagDiss(iiAlt,i1CH2_), &
                        InMagDiss(iiAlt,iCH3_), &
                         InMagIon(iiAlt,iNP_), &
                         InMagIon(iiAlt,iN2P_), &
                         InMagIon(iiAlt,iCH3P_), &
                         eheat(iiAlt)
      enddo
    close(Unit = UnitTmp_)

! Below, the 0.05 is the 5% of the equinox SLT = 12.00 solar EUV
! value. 
! Tom E. Cravens (2006, private communication).

   do iBlock = 1,nBlocks 

  MagDissRateS(1:nLons,1:nLats,1:nAlts,1:nSpeciesTotal,iBlock) = 0.0
  MagIonRateS(1:nLons,1:nLats,1:nAlts,1:nIons - 1,iBlock) = 0.0

    do iAlt = 1,nAlts 
     do iLat = 1,nLats 
      do iLon = 1,nLons 

       EHeatingRate(iLon,iLat,iAlt,iBlock) = &
           eheat(iAlt)*(0.05)*HeatingEfficiency_CB(iLon,iLat,iAlt,iBlock)

         do iSpecies = 1,nSpeciesTotal
           MagDissRateS(iLon,iLat,iAlt,iSpecies,iBlock) = &
                     InMagDiss(iAlt,iSpecies)*(0.05)
         enddo

         do iIon = 1,nIons - 1
           MagIonRateS(iLon,iLat,iAlt,iIon, iBlock) = &
                     InMagIon(iAlt,iIon)*(0.05)
         enddo


      enddo
     enddo
    enddo
   enddo

end subroutine init_magheat


subroutine set_planet_defaults

  use ModInputs

  return

end subroutine set_planet_defaults


! \
! BEGIN The HCN RADIATIVE COOLING CODE

      subroutine calc_radcooling(iBlock)
!  ---------------------------------------------------------------------
!  Purpose: 
!  =======
!  This is the Cooling Routine to be used for the Titan TGCM
!  This version incorporates all heating and cooling terms into
!  the full line by line radiative transfer calculations.
!     Date            Programmer                   Version
!  ==========     ==================            ==============
!  09-08-04           Jared Bell                 Original(1.0) 
!  12-02-04           Jared Bell                 Update (1.1)   
!  12-07-04           Jared Bell                 Update (1.2)   
!  01-27-05           Jared Bell                 Consistency Check   
!  01-27-05           Jared Bell                 Added Heating 1 
!  01-27-05           Jared Bell                 Added Heating 2 
!  01-28-05           Jared Bell                 Added Dynamic Memory 
!  09-07-05           Jared Bell                 Modify to use TAU 
!  09-12-05           Jared Bell                 Utilized Linear Interp. 
!  09-16-05           Jared Bell                 Expint more accurate 
!                                                and more efficient 
!  03-02-06           Jared Bell                 GITM                 
!  0l-07-07           Jared Bell                 Computational Efficiency
!                                                Update 
!  -----------------------------------------------------------------

      use ModPlanet
      use ModGITM
      use ModConstants, only:  Boltzmanns_Constant, Speed_Light, &
                               Planck_Constant
      use ModInputs, only:  iDebugLevel
      use ModSources, only: RadCoolingRate
      use ModTime, only: tSimulation, CurrentTime
      use ModIndicesInterfaces

      implicit none

! INPUTS:----------------------------------------------------------

      integer,intent(in):: iBlock 

! LOCAL VARS:------------------------------------------------------

      integer :: i,j,l,m,k,n,ii,tj 
      integer :: iLat, iLon, iAlt, iLine, iFreq

      integer, parameter :: NLVLS = nAlts
      integer, parameter :: NLON = nLons
      integer, parameter :: NLINES = rotlines 
      integer, parameter :: NFREQ = rotfreqs 
      integer, parameter :: NPTS = rotpts 
      integer, parameter :: blon = 1 
      integer, parameter :: elon = nLons 
      integer, parameter :: blat = 1 
      integer, parameter :: elat = nLats 

      real,dimension(1:nLons,1:nAlts) :: &
       TCOOLa

      real ::   &
       FTH

      real,dimension(1:nLons,1:nLats,1:nAlts) ::    &
       T,       & 
       NHCN,    &
       fcp,     &
       TCOOL2,  & 
       TCOOL2a, & 
       COOL,    &
       COOL2,   &
       HEAT1,   &
       HEAT2

      real,dimension(1:nLons,1:nLats,1:NLINES,1:nAlts) ::  &
       gcool,   &
       gheat1,  &
       gheat2,  &
       INTEN,   &
       ALPHAD,  &
       DELTAV,  &
       VULIM,   &
       VLLIM,   &
       VJM,     &
       TM,      &
       ZM,      &
       VTH,     &
       NHCNM

      real,dimension(1:nLons,1:NLINES,1:nAlts) ::  &
       NEW, &
       OLD, &
       Y1,  &
       Y2,  &
       Y

      real,dimension(1:nLons,1:NLINES,1:nAlts) ::  &
       PRODCOOL,   &
       PRODHEAT1,  &
       PRODHEAT2

      real,dimension(1:nLons,1:nAlts) ::    &
       SMCOOL,   &
       SMHEAT1,  &
       SMHEAT2
!
!
      real,dimension(1:nLons,1:nAlts) ::  &
       NHCNT 
!   
! --------------BEGIN TAU VARS--------------------------------------- 
! The Following Variables are used solely within the Tau Integration
! 
      real :: DFACT
!
      real,dimension(1:NFREQ,1:nLons,1:nLats,1:NLINES,1:nAlts) ::  &
       TAU,  &
       DOP,  &  ! the DOPPLER Lineshape
       X,  &
       V, &
       DopFactor

!\
! Used to calculate the Planck Blackbody function (Bv) at each
! Gaussian Frequency Point
!/
      real,dimension(1:NFREQ,1:nLons,1:nLats,1:NLINES,1:nAlts) ::  &
       Bv,  &
       Balpha

      real :: &
       Bfactor
!
!      VOI,  &  ! The VOIGT Lineshape
!      Y,  &   ! Only for Vogit Profiles
!
! --------------END TAU VARS--------------------------------------- 
! ----------------------------------------------------------------- 
!\ ----------------------------------------------------- 
!  Begin Substituting in for Local variables
!/ ----------------------------------------------------- 

      
          if(floor((tSimulation - dT)/dTCooling) == floor(tSimulation/dTCooling)) return
!\
! This will shut down RadCooling For Testing
!/
!      RadCoolingRate(:,:,:,iBLock) = 0.0
!      return
!\
! -------------------------------------------
!/
  call start_timing("calc_radcooling")

      HEAT1 = 0.0D0
      HEAT2 = 0.0D0

      T(1:nLons,1:nLats,1:nAlts) = &
       Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*  &
       TempUnit(1:nLons,1:nLats,1:nAlts)

      NHCN(1:nLons,1:nLats,1:nAlts) = &
       NDensityS(1:nLons,1:nLats,1:nAlts,iHCN_,iBlock)*1.0E-06 !m^3 to cm^3


! ------------------------------------------------------------------------
! INTERPOLATING INTENSITY DATA TO TEMPERATURE GRID -----------------------
! ------------------------------------------------------------------------

! Define s(1:NLINES,1:nAlts), and inten(nlon,nlats,nlines,nlev)

!      if(iDebugLevel > 4) then
!      write(*,*) 'Begin Intensity Interoplations' 
!      endif

      do iAlt = 1, nAlts
       do iLine = 1, NLINES
        do iLat = 1, nLats 
         do iLon = 1, nLons
!\
! The factor of 100.0 converts speed of light in m/s to cm/s
! this gives INTEN in SI units
!/

         INTEN(iLon,iLat,iLine,iAlt) =             &
           ( Speed_Light*100.0)*(1.0e-4)*          &
           ( pmat(iLine,9)                         &
           + pmat(iLine,8)*T(iLon,iLat,iAlt)       &
           + pmat(iLine,7)*(T(iLon,iLat,iAlt)**2)  &
           + pmat(iLine,6)*(T(iLon,iLat,iAlt)**3)  &
           + pmat(iLine,5)*(T(iLon,iLat,iAlt)**4)  &
           + pmat(iLine,4)*(T(iLon,iLat,iAlt)**5)  &
           + pmat(iLine,3)*(T(iLon,iLat,iAlt)**6)  &
           + pmat(iLine,2)*(T(iLon,iLat,iAlt)**7)  & 
           + pmat(iLine,1)*(T(iLon,iLat,iAlt)**8)  )   


          enddo
         enddo
        enddo
       enddo
!\
!  INTENSITY is the Intensity (1:nlon,1:nlats,1:nlines,1:nAlts)
!/

     if(iDebugLevel > 4) then
       write(*,*) '==> Finished with Intensity Interpolation'
     endif

!     if(iDebugLevel > 4) then
!
!       do k = 1,nAlts
!       write(*,*) '==> Max(Inten(:,:,:',k,') = ', &
!                 Maxval(Inten(1:nLons,1:nLats,1:NLINES,k)) 
!       enddo
!
!     endif

! ---------------
! FOR LOOP BEGIN-
! ---------------
!      FACTOR1 = ATMTOPAS*KB*(Speed_Light*100.0)*(1.0D6)
!      DO i = 1,NLVLS
!         DO j = 1,NLINES
!      SUBALPHL = FACTOR1*GAMMA(j)*T(i)*N(i)
!      ALPHAL(j,i)= SUBALPHL*(1.0 - NHCN(i)/N(i))*SQRT(300./T(i))
!         END DO
!      END DO
!-------------------I NEED TOTAL DENSITY SOMEHOW FOR THE VOIGT !!!!------
! ---------------
! END FOR LOOP --
! ---------------

! ------------------------------------------------------------------------
! BUILDING FULL 4-D MATRICES TM,VJM,ZM,NHCNM (1:nLons,1:nLats,1:nlines,1:nAlts)
! ------------------------------------------------------------------------
! ---------------
! FOR LOOP BEGIN-
! ---------------
!
      do iAlt = 1, nAlts 
       do iLine = 1, NLINES  
        do iLat = 1, nLats 
         do iLon = 1, nLons 

         VJM(iLon,iLat,iLine,iAlt) = freqhz(iLine)        ! Lines in Hz
          ZM(iLon,iLat,iLine,iAlt) = Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0   ! km -> m
          TM(iLon,iLat,iLine,iAlt) = T(iLon,iLat,iAlt)       ! Temp in  K
       NHCNM(iLon,iLat,iLine,iAlt) = &
             NDensityS(iLon,iLat,iAlt,iHCN_,iBlock)*1.0E-06     ! m^-3 -> cm^-3

         enddo
        enddo
       enddo
      enddo
!
! ---------------
! END FOR LOOP --
! ---------------
! ------------------------------------------------------------------------
! END BUILDING 4-D MATRICES ZM,TM,VJM,NHCNM 
! ------------------------------------------------------------------------
!
!
! ------------------------------------------------------------------------
! DOPPLER HALFWIDTHS AT GIVEN TEMPERATURES -------------------------------
! ------------------------------------------------------------------------

        FTH = SQRT(2.0*(Boltzmanns_Constant/Mass(iHCN_)))
        VTH = FTH*SQRT(TM)
        ALPHAD = (VJM*(VTH/(Speed_Light)))
        DELTAV = 5.0*ALPHAD
        VULIM = VJM + DELTAV
        VLLIM = VJM - DELTAV

! ------------------------------------------------------------------------
! END DOPPLER HALFWIDTHS AT GIVEN TEMPERATURES ---------------------------
! ------------------------------------------------------------------------
!
! ------------------------------------------------------------------------
! BEGIN CALCULATING TAU 
! ------------------------------------------------------------------------
!

 
!----------------
!WHILE LOOP BEGIN
!----------------
      m = 1
      do 
      if (2**m .EQ. NFREQ) exit
         m = m + 1
      enddo
!----------------
!END WHILE LOOP 
!---------------

! ---------------
! FOR LOOP BEGIN-
! ---------------
      do k = 1 , NFREQ/2

      V(2*k - 1,:,:,:,:)=(Qf(k,m)*(VULIM - VLLIM) + VULIM + VLLIM)*.5
      V(2*k,:,:,:,:)=(-1.0D0*Qf(k,m)*(VULIM-VLLIM) + VULIM+ VLLIM)*.5

      enddo

! ---------------
! END FOR LOOP --
! ---------------

      DFACT = SQRT(log(2.0)/PI)

! Now Begin The Line Shape Calculations:---------------------------------
!      DO k = 1, NFREQ
!         X(k,:,:,:,:) = SQRT(log(2.0))*((V(k,:,:,:,:) - VJM)/ALPHAD) 
!         DOP(k,:,:,:,:) = DFACT*(EXP(-( X(k,:,:,:,:)**2 ))/ALPHAD)
!         DopFactor(k,:,:,:,:) = &
!             1.0/(sqrt(PI)*ALPHAD(:,:,:,:) )
!      END DO
! END Line Shape Calculations:---------------------------------


   do iAlt = 1, nAlts
      do iLine = 1, NLINES
         do iLat = 1, nLats
            do iLon = 1, nLons
               do iFreq = 1, NFREQ

         X(iFreq,iLon,iLat,iLine,iAlt) =  &
             SQRT( log(2.0) )*  &
         (    &
           ( V(iFreq,iLon,iLat,iLine,iAlt)  - VJM(iLon,iLat,iLine,iAlt) )/ & 
              ALPHAD(iLon,iLat,iLine,iAlt)    &
          ) 

         DOP(iFreq,iLon,iLat,iLine,iAlt) =  &
               (DFACT/ALPHAD(iLon,iLat,iLine,iAlt))*    & 
               EXP( -( X(iFreq,iLon,iLat,iLine,iAlt)**2) ) 
                   

         DopFactor(iFreq,iLon,iLat,iLine,iAlt) =  1.0/ &
             (  sqrt(PI) * &
               ALPHAD(iLon,iLat,iLine,iAlt) &
              )

       Balpha(iFreq,iLon,iLat,iLine,iAlt)  =             &
               (  Planck_Constant/Boltzmanns_Constant)*  &
               (  V(ifreq,iLon,iLat,iLine,iAlt) /        &
                       TM(iLon,iLat,iLine,iAlt)  )

       Bfactor = 1.0/ &
                 ( EXP(Balpha(iFreq,iLon,iLat,iLine,iAlt) ) - 1.0 ) 

       Bv(iFreq,iLon,iLat,iLine,iAlt) =                   &
              ( ( 2.0*Planck_Constant ) / 	          & 
              ( (Speed_Light)**2 )  ) *                   &
              ( V(iFreq,iLon,iLat,iLine,iAlt)**3 )*Bfactor

               enddo
             enddo
          enddo
       enddo
    enddo

   do iLat = 1,nLats   ! Major Outer Latitude Loop

!\
! Calculate Optical Depth Profiles at each Frequency:--------------------
! We send in INTEN:  Svj
!            NHCNM:  nHCN
!              DOP:  Phi (line shape)
! We get out Optical Depth:  Tau (freqency gausspts, lons, lats, lines, alts)
! /

      CALL fgausstau(INTEN(1:nLons,iLat,1:NLINES,1:nAlts), &
                     NHCNM(1:nLons,iLat,1:NLINES,1:nAlts), &
                     ZM(1:nLons,iLat,1,1:nAlts), &
                     DOP(1:NFREQ,1:nLons,iLat,1:NLINES,1:nAlts), &
                     TAU(1:NFREQ,1:nLons,iLat,1:NLINES,1:nAlts) )

     if(iDebugLevel > 4) then
       write(*,*) '==> Finished Calculating TAU !'
     endif


! END Calculate Optical Depth Profiles at each Frequency:----------------

! ------------------------------------------------------------------------
! Cooling CALCULATIONS
! ------------------------------------------------------------------------

   !   CALL fgausscool(VLLIM(1:nLons,iLat,1:NLINES,1:nAlts), &
   !                   VULIM(1:nLons,iLat,1:NLINES,1:nAlts), &
   !                   Bv(1:NFREQ,1:nLons,iLat,1:NLINES,1:nAlts), &
   !                   DOP(1:NFREQ,1:nLons,iLat,1:NLINES,1:nAlts), &
   !                    gcool(1:nLons,iLat,1:NLINES,1:nAlts) )
! \
! ------------------------------------------------------------------------
! Cooling CALCULATIONS
! Cool-To-Space Portion: 
! ------------------------------------------------------------------------
!/

      NEW = 0.0
      OLD = 0.0

      j = 1

      DO 
        IF (2**j .EQ. NFREQ) EXIT
        j = j + 1
      END DO

      do iFreq = 1 , NFREQ/2
        Y1 = Bv(2*iFreq - 1,:,iLat,:,:)*DOP(2*iFreq - 1,:,iLat,:,:)
        Y2 = Bv(2*iFreq,:,iLat,:,:)*DOP(2*iFreq,:,iLat,:,:) 
        Y = Y1 + Y2
        OLD = NEW 
        NEW = OLD + Qd(iFreq,j)*Y
      enddo

      gcool(:,iLat,:,:) = (VULIM(:,iLat,:,:) - VLLIM(:,iLat,:,:) )*NEW*0.5


!\
! Multiply by the Intensity of the lines and a geometric factor of 4pi
!/
      PRODCOOL(1:nLons,1:NLINES,1:nAlts) =  &
                       gcool(1:nLons,iLat,1:NLINES,1:nAlts)* &
                       INTEN(1:nLons,iLat,1:NLINES,1:nAlts)*(4.0D0*PI)
!\
! Sum over the lines
!/
      SMCOOL(1:nLons,1:nAlts) = &
                SUM(PRODCOOL(1:nLons,1:NLINES,1:nAlts),2)
!\
! Multiply by HCN density 
!/

      COOL(1:nLons,iLat,1:nAlts)  =  &
                      SMCOOL(1:nLons,1:nAlts)* &
                        NHCN(1:nLons,iLat,1:nAlts)*1.0e+06   ! J/(m^3 s) 

!      COOL(1:nLons,iLat,1:nAlts)  =  &
!                      SMCOOL(1:nLons,1:nAlts)* &
!                        NHCN(1:nLons,iLat,1:nAlts)      ! J/(cm^3 s) 
!
!      COOL(1:nLons,iLat,1:nAlts)  =  &
!                      SMCOOL(1:nLons,1:nAlts)* &
!                        NHCN(1:nLons,iLat,1:nAlts)*1.0e7   ! ergs/(cm^3 s) 

! \
! ------------------------------------------------------------------------
! END The Cool-To-Space Calculation: 
! ------------------------------------------------------------------------
!/
      CALL heat1gauss(Bv(1:NFREQ,1:nLons,iLat,1:NLINES,1:nAlts), &
                  INTEN(1:nLons,iLat,1:NLINES,1:nAlts), &
                  NHCNM(1:nLons,iLat,1:NLINES,1:nAlts), &
                 gheat1(1:nLons,iLat,1:NLINES,1:nAlts), &
            TAU(1:NFREQ,1:nLons,iLat,1:NLINES,1:nAlts), &
            DOP(1:NFREQ,1:nLons,iLat,1:NLINES,1:nAlts), &
                  VULIM(1:nLons,iLat,1:NLINES,1:nAlts), &
                  VLLIM(1:nLons,iLat,1:NLINES,1:nAlts) )

      PRODHEAT1(1:nLons,1:NLINES,1:nAlts) =  &
                       gheat1(1:nLons,iLat,1:NLINES,1:nAlts)* &
                       INTEN(1:nLons,iLat,1:NLINES,1:nAlts)*(2.0D0*PI)

      SMHEAT1(1:nLons,1:nAlts) = &
                SUM(PRODHEAT1(1:nLons,1:NLINES,1:nAlts),2)

      HEAT1(1:nLons,iLat,1:nAlts)  =  &
                      SMHEAT1(1:nLons,1:nAlts)* &
                         NHCN(1:nLons,iLat,1:nAlts)*1.0e+06  ! J/(m^3 s) 

!      HEAT1(1:nLons,iLat,1:nAlts)  =  &
!                      SMHEAT1(1:nLons,1:nAlts)* &
!                         NHCN(1:nLons,iLat,1:nAlts)      ! J/(cm^3 s) 
!
!
!      HEAT1(1:nLons,iLat,1:nAlts)  =  &
!                      SMHEAT1(1:nLons,1:nAlts)* &
!                         NHCN(1:nLons,iLat,1:nAlts)*1.0e7   ! ergs/(cm^3 s) 
!
      CALL heat2gauss( &
              V(1:NFREQ,1:nLons,iLat,1:NLINES,1:nAlts), &
                  INTEN(1:nLons,iLat,1:NLINES,1:nAlts), &
                 ALPHAD(1:nLons,iLat,1:NLINES,1:nAlts), &
                  NHCNM(1:nLons,iLat,1:NLINES,1:nAlts), &
                 gheat2(1:nLons,iLat,1:NLINES,1:nAlts), &
                     TM(1:nLons,iLat,1:NLINES,1:nAlts), &
                     ZM(1:nLons,iLat,1:NLINES,1:nAlts), &
                    VJM(1:nLons,iLat,1:NLINES,1:nAlts), &
            TAU(1:NFREQ,1:nLons,iLat,1:NLINES,1:nAlts), &
            DOP(1:NFREQ,1:nLons,iLat,1:NLINES,1:nAlts), &
                  VULIM(1:nLons,iLat,1:NLINES,1:nAlts), &
                  VLLIM(1:nLons,iLat,1:NLINES,1:nAlts), &
                  NLON,NLVLS,                           &
                  Bv(1:NFREQ,1:nLons,iLat,1:NLINES,1:nAlts) )

      PRODHEAT2(1:nLons,1:NLINES,1:nAlts) =  &
                       gheat2(1:nLons,iLat,1:NLINES,1:nAlts)* &
                       INTEN(1:nLons,iLat,1:NLINES,1:nAlts)*(2.0D0*PI)
      SMHEAT2(1:nLons,1:nAlts) = &
                SUM(PRODHEAT2(1:nLons,1:NLINES,1:nAlts),2)

      HEAT2(1:nLons,iLat,1:nAlts)  =  &
                      SMHEAT2(1:nLons,1:nAlts)* &
                         NHCN(1:nLons,iLat,1:nAlts)*1.0E+06      ! J/(m^3 s) 
!
!
!      HEAT2(1:nLons,iLat,1:nAlts)  =  &
!                      SMHEAT2(1:nLons,1:nAlts)* &
!                         NHCN(1:nLons,iLat,1:nAlts)      ! J/(cm^3 s) 
!      HEAT2(1:nLons,iLat,1:nAlts)  =  &
!                      SMHEAT2(1:nLons,1:nAlts)* &
!                         NHCN(1:nLons,iLat,1:nAlts)*1.0E+07      ! Ergs/(cm^3 s) 
!

    RadCoolingRate(1:nLons,iLat,1:nAlts,iBlock) = &
               COOL(1:nLons,iLat,1:nAlts) - &
               HEAT1(1:nLons,iLat,1:nAlts) - &
               HEAT2(1:nLons,iLat,1:nAlts) 

   enddo  ! End Outer Latitude-Loop

!open(Unit = 16, file = 'radcool.txt',status = 'new')
!      do i = 1,nAlts
!       write(16,*) RadCoolingRate(1,1,i,iBlock)*10.0  ! in Ergs/cm^3/s
!      enddo
!close(Unit = 16)

   if(iDebugLevel > 4) then

       	  do k = 1,nAlts
              write(*,*) '==> Max(TAU(:,:,:,:',k,'), Min = ', &
                   Maxval(TAU(1:NFREQ,1:nLons,1:nLats,1:NLINES,k)), & 
                   Minval(TAU(1:NFREQ,1:nLons,1:nLats,1:NLINES,k))
          enddo

      do i = 1,nAlts
       write(*,*) 'Maxval(COOL(1:nLons,1:nLats,i)) =', &
                   Maxval(COOL(1:nLons,1:nLats,i)), 'J/m^3/s'
      enddo

      do i = 1,nAlts
       write(*,*) 'Maxval(HEAT1(1:nLons,1:nLats,i)) =', &
                   Maxval(HEAT1(1:nLons,1:nLats,i)), 'J/m^3/s'
      enddo

      do i = 1,nAlts
       write(*,*) 'Maxval(HEAT2(1:nLons,1:nLats,i)) =', &
                   Maxval(HEAT2(1:nLons,1:nLats,i)), 'J/m^3/s'
      enddo

    endif

!      do i = 1,nAlts
!       write(*,*) 'Heat1 =', &
!       Heat1(1,1,i), 'J/m^3/s'
!      enddo
!      do i = 1,nAlts
!       write(*,*) 'Heat2 =', &
!       Heat2(1,1,i), 'J/m^3/s'
!      enddo
!      do i = 1,nAlts
!       write(*,*) 'RadCoolingRate =', &
!       RadCoolingRate(1,1,i,iBlock), 'J/m^3/s'
!      enddo
!      do i = 1,nAlts
!       write(*,*) 'Log10(RadCoolingRate) =', &
!       Log10(10.0*RadCoolingRate(1,1,i,iBlock)), 'Log(Ergs/cm^3/s)'
!      enddo
!    do i = 1,nLons
!     do j = 1,nLats
!      do k = 1,nAlts
!             write(*,*) 'Cool, Heat1, Heat2, TOTAL =', &
!             Cool(i,j,k),HEAT1(i,j,k),HEAT2(i,j,k), &
!             Cool(i,j,k) - HEAT1(i,j,k) - HEAT2(i,j,k)
!        write(*,*) 'Log10(TOTAL) =', &
!        LOG10(Cool(i,j,k) - HEAT1(i,j,k) - HEAT2(i,j,k))
!
!             write(*,*) 'RadCoolingRate =', &
!                  RadCoolingRate(1,1,k,iBlock), 'J/m^3/s'
!
!      enddo
!     enddo
!    enddo



  call end_timing("calc_radcooling")

      contains


!---------------------------------------------------------------------------
! HEAT1GAUSS DECLARATION
!---------------------------------------------------------------------------
!
      subroutine heat1gauss(Planck,SJ,NHCNM,gheat1,TAU2,PHI,B,A)

      implicit none
!
! INPUT Variables:------------------------------------------
!
      REAL, INTENT(IN), DIMENSION(nLons,NLINES,nAlts) ::   &
       SJ,      &
       NHCNM,   &
       B,       &
       A  
!
      REAL,INTENT(IN),DIMENSION(NFREQ,nLons,NLINES,nAlts) ::  &
       Planck,  &
       TAU2,  &
       PHI
!
      REAL, INTENT(OUT), DIMENSION(nLons,NLINES,nAlts) ::   &
       gheat1
!
! LOCAL Variables:------------------------------------------
!
      INTEGER :: i,j,k,m
      INTEGER :: iAlt,iLine,iLon,iFreq
      INTEGER :: iiAlt
!
      REAL, DIMENSION(NFREQ,nLons,NLINES,nAlts) ::   &
       TOTAL

      REAL ::   &
       TOL
!
      REAL, DIMENSION(nLons,NLINES,nAlts) ::   &
       TOTAL0,  &
       FACTOR20
!
      REAL, DIMENSION(nLons,NLINES,nAlts - 1) ::   &
       TAUUP,  &
       TAUDWN,  &
       DTAU,  &
       EXPTAU

      REAL, DIMENSION(nLons,NLINES,nAlts) ::  &
        S,  &
        OLD 

        TOTAL0 = 0.0
        TOTAL = 0.0
!
! heat1gauss begin:------------------------------------------
      i = 1
      j = 1
      k = 1
      m = 1
      iiAlt = 1
      TOL = 1.0e-5

      FACTOR20 = 0.0D0
      FACTOR20(:,:,1) = 1.0D0

      do iFreq = 1, NFREQ   ! Outer Most Loop over Frequencies

          do iAlt = 1, nAlts
           do iLine = 1, NLINES 
            do iLon = 1,nLons 

              TOTAL0(iLon,iLine,iAlt) =  &
                 Planck(iFreq,iLon,iLine,iiAlt)* &
                    PHI(iFreq,iLon,iLine,iAlt)

             enddo ! iLon
            enddo ! iLine
           enddo ! iAlt
!
           do iAlt = 1, nAlts - 1
             do iLine = 1, NLINES 
               do iLon = 1,nLons 

            TAUUP(iLon,iLine,iAlt) = TAU2(iFreq,iLon,iLine,1)

            TAUDWN(iLon,iLine,iAlt) = TAU2(iFreq,iLon,iLine,iAlt + 1)

!            DTAU(iLon,iLine,iAlt) =  ABS(TAUUP(iLon,iLine,iAlt) - TAUDWN(iLon,iLine,iAlt) )

            DTAU(iLon,iLine,iAlt) =   &
                  max( TOL, TAUUP(iLon,iLine,iAlt) - TAUDWN(iLon,iLine,iAlt) )

               enddo
             enddo
           enddo 

      CALL expint1(DTAU,EXPTAU,nLons,nAlts - 1)

           do iAlt = 1, nAlts - 1
             do iLine = 1, NLINES 
               do iLon = 1,nLons 

             FACTOR20(iLon,iLine,iAlt + 1) =    &
            EXP(-DTAU(iLon,iLine,iAlt) ) -  &
                 DTAU(iLon,iLine,iAlt)*EXPTAU(iLon,iLine,iAlt) 


               enddo
             enddo
           enddo 

           do iAlt = 1, nAlts - 1
             do iLine = 1, NLINES 
               do iLon = 1,nLons 

                  TOTAL(iFreq,iLon,iLine,iAlt) =  &
                     FACTOR20(iLon,iLine,iAlt)*  &
                       TOTAL0(iLon,iLine,iAlt)

               enddo
             enddo
           enddo 
!
      enddo ! Outer iFreq loop

      S = 0.0D0
      OLD = 0.0D0
      m = 1

      DO 
        IF (2**m .EQ. NFREQ) EXIT
        m = m + 1
      END DO

      DO k = 1 , NFREQ/2
        S = OLD + Qd(k,m)*(TOTAL(2*k,:,:,:) + TOTAL(2*k - 1,:,:,:))
        OLD = S
      END DO

      gheat1 = (0.5)*(B-A)*S

      end subroutine heat1gauss
! --------------------------------------------------------------------------
! FGAUSSTAU DECLARATION
!---------------------------------------------------------------------------
      subroutine fgausstau(SJ,NHCNM,Z,SHAPE,ITAU)

      IMPLICIT NONE
! INPUTS ----------------------------------------------
!
      REAL,INTENT(IN),DIMENSION(nLons,nAlts) ::  &
       Z
!
      REAL,INTENT(IN),DIMENSION(nLons,NLINES,nAlts) ::  &
       SJ,  &
       NHCNM
!
      REAL,INTENT(IN),DIMENSION(NFREQ,nLons,NLINES,nAlts) ::  &
       SHAPE
!
      REAL,INTENT(OUT),DIMENSION(NFREQ,nLons,NLINES,nAlts)::  &
       ITAU 
!
! END INPUTS ----------------------------------------------
!
! Local Vars ----------------------------------------------
      INTEGER ::   &
       j,  &
       k,  &
       m,  &
       n,  &
       i,  &
       loc

      INTEGER ::   &
       iAlt, &
       iLon, &
       iLine

      INTEGER ::   &
       start,  &
       end,  &
       rate

!
      REAL,DIMENSION(NFREQ,nLons,NLINES,nAlts) ::  &
       S,  &
       OLD,  &
       BF,  &
       AF
!
      REAL,DIMENSION(nLons,NLINES,nAlts) ::  &
       TOTAL
!
      REAL,DIMENSION(NLINES,nAlts) ::  &
       ITERP2 
!
      REAL,DIMENSION(NPTS,NFREQ,nLons,NLINES,nAlts) ::  &
       DEPTH
!
      REAL,DIMENSION(NPTS,nLons,nAlts) ::  &
        x
!
      REAL,DIMENSION(nLons,nAlts) ::   &
       A,  &
       B
!
!----------END Local Vars-----------------------

!----------------
!WHILE LOOP BEGIN
!----------------
      m = 1
      DO 
      IF (2**m .EQ. NPTS) EXIT
         m = m + 1
      END DO
!----------------
!END WHILE LOOP 
!---------------

      A = Z
      do iAlt = 1, nAlts
         B(:,iAlt) = Z(:,nAlts) 
      enddo
!

      DO k = 1,NFREQ
        DO j = 1, NLINES
         AF(k,:,j,:) = A
         BF(k,:,j,:) = B
        END DO
      END DO
!
! Now Our Endpoints are B, A which are NLON,NLVLS in Size
!
      DO k = 1 , NPTS/2

      x(2*k - 1,:,:)= 0.5*(Qf(k,m)*(B - A) + B + A)
      x(2*k,:,:)= 0.5*(-Qf(k,m)*(B - A) + B + A)

      END DO

      DO k = 1, NFREQ/2
           TOTAL = (SHAPE(2*k - 1,:,:,:)*SJ*NHCNM)*1.0D9
         DO i = 1,nLons
            DO j = 1,nAlts

! x(:,i,j) is size NLVLS
! Z(i,:) is size NLVLS
! TOTAL(i,:,:) is size NLINES, NLVLS

! Linear Interpolation of Variables to the Gaussian Quadrature points:

              DO n = 1,NPTS

               CALL locate(Z(i,:),x(n,i,j),loc,NLON,NLVLS)
      DEPTH(n,2*k - 1,i,:,j) = TOTAL(i,:,loc) +  &
       (x(n,i,j) - Z(i,loc))* &
       ((TOTAL(i,:,loc+1) - TOTAL(i,:,loc))/(Z(i,loc+1) - Z(i,loc)))  

      DEPTH(n,2*k,i,:,j) = DEPTH(n,2*k - 1,i,:,j) 

              END DO ! n loop
            END DO ! j loop
         END DO ! i loop
      END DO !k loop


! END Linear Interpolation 

      OLD = 0.0D0
      DO k = 1, NPTS/2

        S = OLD + Qd(k,m)*(DEPTH(2*k - 1,:,:,:,:) +  &
        DEPTH(2*k,:,:,:,:))

        OLD = S
      END DO

      ITAU = 0.5*(BF-AF)*S


      end subroutine fgausstau
!---------------------------------------------------------------------------
! FGAUSSHEAT2 DECLARATION
!---------------------------------------------------------------------------
!
      subroutine heat2gauss(V,SJ,ALPHAD,NHCNM,gheat2,TM,  &
                           ZM,VJM,TAU,PHI,VULIM,VLLIM,NLON,NLVLS, &
                           Planck)

      IMPLICIT NONE
!
! INPUT Variables:------------------------------------------
!
      INTEGER, INTENT(IN) :: NLVLS 
      INTEGER, INTENT(IN) :: NLON 
!
      REAL, INTENT(IN), DIMENSION(NLON,NLINES,NLVLS) ::  &
       ALPHAD,  &
       SJ,  &
       NHCNM,  &
       TM,  &
       ZM,  &
       VJM,  &
       VULIM,  &
       VLLIM  
!
      REAL,INTENT(IN),DIMENSION(NFREQ,NLON,NLINES,NLVLS) ::  &
       V,  &
       TAU,  &
       Planck,  &
       PHI
!
      REAL, INTENT(OUT), DIMENSION(NLON,NLINES,NLVLS) ::  &
       gheat2
!
! LOCAL Variables:------------------------------------------

      INTEGER :: i,j,k,m,p,t,loc,ilon,iline

      INTEGER ::  & 
       start,  &
       end,  &
       rate 

      REAL,DIMENSION(NLON,NLINES,NLVLS) ::  &
       ATAU,  &
       BTAU,  &
       A1TAU,  &
       B1TAU

      REAL,DIMENSION(NLON,NLINES,NLVLS) ::   &
       ALPHA,   &
       B,   &
       FACTOR,   &
       BH,   &
       BL,   &
       SH,   &
       SL,   &
       OLDH,   &
       OLDL,   &
       S,   &
       OLD
!
      REAL,DIMENSION(NPTS,NLON,NLINES,NLVLS) ::  &
       XH,  & ! Higher Integral (above the value of Z)
       XL,  &  ! Lower Integral (below the value of Z)
       TOTH,  & ! Higher Integral (above the value of Z)
       TOTL,  &  ! Lower Integral (below the value of Z)
       TESTH,  &
       TESTL,  &
       TAUV,  &
       EXPH,  &
       EXPL
!
      REAL,DIMENSION(NFREQ,NLON,NLINES,NLVLS) ::  &
       QTOT 

     ilon = 1
     iline = 1
     i = 1

      m = 1
!----------------
!WHILE LOOP BEGIN
!----------------
      DO 
        IF (2**m .EQ. NPTS) EXIT
        m = m + 1
      END DO
!----------------
!END WHILE LOOP 
!---------------
!
! ---------------
! FOR LOOP BEGIN-
! ---------------
! BEGIN OUTER FREQUENCY LOOPING----------------------------------
! CAVEAT:  Note Problems associateed with moving m, loop above
! outside of k-loop below!!  
!
!       ALPHA = (Planck_Constant/Boltzmanns_Constant)*&
!               (V(1,:,:,:)/TM)
!       FACTOR = 1.0D0/(EXP(ALPHA) - 1.0D0) 
!       B = ((2.0D0*Planck_Constant)/((Speed_Light)**2))*&
!            (V(1,:,:,:)**3.0D0)*FACTOR

      B = Planck(1,:,:,:)

      DO k = 1,NFREQ/2

       TOTH = 0.0D0
       TOTL = 0.0D0


       BTAU = TAU(k,:,:,:)
       A1TAU = TAU(k,:,:,:)

       DO ilon = 1,NLON
          DO iline = 1,NLINES
             DO i = 1,NLVLS
        ATAU(ilon,iline,i) = TAU(k,ilon,iline,1)
        B1TAU(ilon,iline,i) = TAU(k,ilon,iline,NLVLS)
             END DO
          END DO
       END DO

! Gives the Gaussain Points in terms of TAU Values

      DO p = 1 , NPTS/2
      XH(2*p - 1,:,:,:) = (Qf(p,m)*(A1TAU - B1TAU) + A1TAU + B1TAU)/2.0D0 
      XH(2*p,:,:,:) = (-Qf(p,m)*(A1TAU - B1TAU) + A1TAU + B1TAU) &
          /2.0D0
      XL(2*p - 1,:,:,:) = (Qf(p,m)*(ATAU - BTAU) + ATAU + BTAU)/2.0D0 
      XL(2*p,:,:,:) = (-Qf(p,m)*(ATAU - BTAU) + ATAU + BTAU) &
       /2.0D0
      END DO

! XH, XL CHECK OUT !


! Linear Interpolation of the Black Body Function:---------------------
! N.B. Below, we do not explicitly call or loop over NLINES variable for
! the location variable, loc.  It was found that (to a very good approximation)
! that all key variables depended on the variable NLINES in exactly the same
! fashion.  Thus, a loop over NLINES was redundant and only cost us
! CPU time.  However, this trend is violated near the top of the model. 
! However, the heating function is less important up at these altitudes and
! so represents a minor error, if any at all.
! 
      DO p = 1 , NPTS
         DO i = 1,NLON 
             DO t = 1,NLVLS

               CALL locate(TAU(k,i,1,:),XH(p,i,1,t),loc,NLON,NLVLS)
      BH(i,:,t) = B(i,:,loc) +  &
       (XH(p,i,:,t) - TAU(k,i,:,loc))* &
       ((B(i,:,loc+1) - B(i,:,loc))/(TAU(k,i,:,loc+1) - TAU(k,i,:,loc)))  

               CALL locate(TAU(k,i,1,:),XL(p,i,1,t),loc,NLON,NLVLS)
      BL(i,:,t) = B(i,:,loc) +  &
       (XL(p,i,:,t) - TAU(k,i,:,loc))* &
       ((B(i,:,loc+1) - B(i,:,loc))/(TAU(k,i,:,loc+1) - TAU(k,i,:,loc)))  

                END DO ! t loop
            END DO ! i loop
         END DO ! p loop


! END Linear Interpolation of the Black Body Function:---------------------

      EXPH = 0.0D0 
      EXPL = 0.0D0 

      DO p = 1 , NPTS
      TAUV(p,:,:,:) = TAU(k,:,:,:)
      END DO ! p loop

      TESTH = ABS(TAUV - XH)  
      TESTL = ABS(TAUV - XL)  

      CALL expint2(TESTL(:,:,:,2:NLVLS),EXPL(:,:,:,2:NLVLS), &
                   NLON,NLVLS - 1)
      CALL expint2(TESTH(:,:,:,1:NLVLS-1),EXPH(:,:,:,1:NLVLS-1), &
                   NLON,NLVLS - 1)
  
      DO p = 1 , NPTS
      TOTH(p,:,:,:) = BH*EXPH(p,:,:,:)
      TOTL(p,:,:,:) = BL*EXPL(p,:,:,:)
      END DO ! p loop


! Testing the New Expint Function----------------------------------


! Summing over Gauss Altitude Points-----------------------
!
      OLDH = 0.0D0
      OLDL = 0.0D0
      DO p = 1 , NPTS/2
        SH = OLDH + Qd(p,m)*(TOTH(2*p,:,:,:) + TOTH(2*p - 1,:,:,:)) 
        SL = OLDL + Qd(p,m)*(TOTL(2*p,:,:,:) + TOTL(2*p - 1,:,:,:)) 
        OLDH = SH
        OLDL = SL
      END DO ! p loop

      QTOT(2*k - 1,:,:,:) = (ABS(ATAU - BTAU))*(SL)/2.0D0 + &
                      (ABS(A1TAU - B1TAU))*(SH)/2.0D0

      QTOT(2*k,:,:,:) = QTOT(2*k - 1,:,:,:) 

      END DO ! k loop End

! END Altitude Integrations----------------------------------- 

! Next, we need to sum over the FREQUENCY Gauss Points:-------
      OLD = 0.0D0
      DO k = 1 , NFREQ/2
       S = OLD +  &
           Qd(k,m)*QTOT(2*k,:,:,:)*PHI(2*k,:,:,:) +   &
           Qd(k,m)*QTOT(2*k - 1,:,:,:)*PHI(2*k - 1,:,:,:) 
       OLD = S
      END DO ! k loop End


      gheat2 = 0.5D0*(VULIM-VLLIM)*S
!
      end subroutine heat2gauss
! ---------------------------------------------------------------------
      subroutine expint2(x,ans,NLON,NLVLS)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NLON 
      INTEGER, INTENT(IN) :: NLVLS 
!
! Input Arguments:--------------------------------------
!
      REAL, INTENT(IN), DIMENSION(NPTS,NLON,NLINES,NLVLS) ::  & 
        x
      REAL, INTENT(OUT), DIMENSION(NPTS,NLON,NLINES,NLVLS) ::  &
        ans

! LOCAL Variables:---------------------------------------
      
!      INTEGER, PARAMETER :: MAXIT = 100 
      INTEGER, PARAMETER :: MAXIT = 10
!      INTEGER, PARAMETER :: MAXIT = 5 
!      INTEGER, PARAMETER :: MAXIT = 2 
!      INTEGER, PARAMETER :: MAXIT = 1 
      INTEGER :: i,j,k 

      REAL, PARAMETER :: EULER = .577215664901532860D0 
      REAL, PARAMETER :: EPS = epsilon(x(1,1,1,1)) 
      REAL, PARAMETER :: BIG = huge(x(1,1,1,1))*EPS 

      REAL, DIMENSION(NPTS,NLON,NLINES,NLVLS) ::  & 
       ID

      REAL, DIMENSION(NPTS,NLON,NLINES,NLVLS) ::   &
       a,  &
       b,  &
       c,  &
       d,  &
       h,  &
       del,  &
       fact,  &
       EXPH,  &
       EXPL

       ID = 1.0D0

! ---------------------------------------------------------
! BEGIN CALCULATION:---------------------------------------
! FORM FOR X >= 1.0 
! ---------------------------------------------------------
       b = x + ID
       c = BIG*ID 
       d = ID/b 
       h = d 
       DO i = 1,MAXIT
          a = -1.0D0*i*i*ID
          b = b + 2.0D0*ID
          d = ID/(a*d + b) 
          c = b + a/c
          del = c*d
          h = h*del
       ENDDO

       EXPH = h*exp(-x)

! ---------------------------------------------------------
! BEGIN CALCULATION:---------------------------------------
! FORM FOR X < 1.0 
! ---------------------------------------------------------

       EXPL = -log(x) - EULER 

       fact = 1.0D0 
       del = 0.0D0 
       DO i = 1, MAXIT
          fact = -(fact*x)/i 
          del = -fact/i
          EXPL = EXPL + del
       ENDDO

       WHERE(x < 1.0D0)
         ans = EXPL
       ELSEWHERE
         ans = EXPH
       END WHERE

      end subroutine expint2
! ---------------------------------------------------------------------
      subroutine expint1(x,ans,NLON,NLVLS)

      implicit none

      INTEGER, INTENT(IN) :: NLON 
      INTEGER, INTENT(IN) :: NLVLS 
!
! Input Arguments:--------------------------------------
!
      REAL, INTENT(IN), DIMENSION(NLON,NLINES,NLVLS) ::  & 
        x
      REAL, INTENT(OUT), DIMENSION(NLON,NLINES,NLVLS) ::  &
        ans

! LOCAL Variables:---------------------------------------
      
!      INTEGER, PARAMETER :: MAXIT = 100 
      INTEGER, PARAMETER :: MAXIT = 10
!      INTEGER, PARAMETER :: MAXIT = 5 
!      INTEGER, PARAMETER :: MAXIT = 2 
!      INTEGER, PARAMETER :: MAXIT = 1 
      INTEGER :: i,j,k 

      REAL, PARAMETER :: EULER = .577215664901532860D0 
      REAL, PARAMETER :: EPS = epsilon(x(1,1,1)) 
      REAL, PARAMETER :: BIG = huge(x(1,1,1))*EPS 

      REAL, DIMENSION(NLON,NLINES,NLVLS) ::  & 
       ID

      REAL, DIMENSION(NLON,NLINES,NLVLS) ::   &
       a,  &
       b,  &
       c,  &
       d,  &
       h,  &
       del,  &
       fact,  &
       EXPH,  &
       EXPL

       ID = 1.0D0

! ---------------------------------------------------------
! BEGIN CALCULATION:---------------------------------------
! FORM FOR X >= 1.0 
! ---------------------------------------------------------
       b = x + ID
       c = BIG*ID 
       d = ID/b 
       h = d 
       DO i = 1,MAXIT
          a = -1.0D0*i*i*ID
          b = b + 2.0D0*ID
          d = ID/(a*d + b) 
          c = b + a/c
          del = c*d
          h = h*del
          IF ( MAXVAL(ABS(del) - 1.0) <= EPS*1.0D13) EXIT
       ENDDO

       EXPH = h*exp(-x)

! ---------------------------------------------------------
! BEGIN CALCULATION:---------------------------------------
! FORM FOR X < 1.0 
! ---------------------------------------------------------

       EXPL = -log(x) - EULER 

       fact = 1.0D0 
       del = 0.0D0 
       DO i = 1, MAXIT
          fact = -(fact*x)/i 
          del = -fact/i
          EXPL = EXPL + del
          IF ( MAXVAL(ABS(del)) <= MAXVAL(ABS(EXPL))*EPS) EXIT
       ENDDO

       WHERE(x < 1.0D0)
         ans = EXPL
       ELSEWHERE
         ans = EXPH
       END WHERE

      END subroutine expint1
! ---------------------------------------------------------------------
      subroutine locate(xx,x,loc,NLON,NLVLS)
! ---------------------------------------------------------------------
!  PURPOSE:
!  =======
!     This is a simple binary search algorithm to be used with splint 
!     Code Adapted from the Press et al Text on Numerical Recipes.
!
!     Given array xx(1:n) and given a value, x, this routine returns the
!     value j such that x lies between xx(j) and xx(j+1).  xx must be 
!     monotonic. j = 0, N mean that x is out of range! 
! ---------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NLVLS 
      INTEGER, INTENT(IN) :: NLON 
      REAL, DIMENSION(NLVLS), INTENT(IN) :: xx
      REAL, INTENT(IN) :: x
      INTEGER, INTENT(OUT) :: loc 

      INTEGER :: n, jl, jm, ju 

      LOGICAL :: ascnd
! ---------------------------------------------------------------------
      n = size(xx)
      ascnd = (xx(n) >= xx(1))
      jl = 0
      ju = n+1
      DO
         IF(ju - jl <= 1) EXIT
         jm = (ju + jl)/2
         IF (ascnd .EQV. (x >= xx(jm))) THEN
               jl = jm
         ELSE 
               ju = jm
         END IF
      END DO

      IF ( x == xx(1)) THEN
          loc = 1
      ELSE IF ( x == xx(n)) THEN
          loc = n - 1
      ELSE
          loc = jl
      END IF

      end subroutine locate

!----------------------------------------------------------------------
!----------------------------------------------------------------------

     end subroutine calc_radcooling 




