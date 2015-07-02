!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CRASH_ModMultiGroup

  use CRASH_ModIonMix
  use CRASH_ModOpacityVoigt,   ONLY:  &
       line_width, voigt_profile, UseVoigt
  use CRASH_ModAtomicDataMix,  ONLY:  &
       nMix, nZ_I, nMixMax,           &
       Concentration_I,               & ! (1:nMixMax)
       IonizPotential_II                ! (1:nZMax,1:nMixMax)
  use CRASH_ModAtomicMass, ONLY: &
       cAtomicMass_I                    ! (1:nZMax)
  use CRASH_ModPartition,      ONLY:  &
       Population_II                    ! (0:nZMax,1:nMixMax)
  use CRASH_ModExcitationData, ONLY:  &
       UseDeltaNEq0Transition,        &
       UseCoreElectron,               &
       n_ground,  n_screened            ! (iZ,nZ)
  use CRASH_ModExcitation, ONLY:      &
       Partition_III,                 & ! (nExcitation,0:nZMax,nMixMax)
       ExcitationEnergy_III,          & ! (nExcitation,0:nZMax-1,nMixMax)
       IonizationPotentialLowering_I, & ! (0:nZMax)
       nExcitation_II                   ! (0:nZMax,nMixMax)
  use CRASH_ModPartition, ONLY: &
       Na, Te, zAv, iZMin_I, iZMax_I    ! (1:nMixMax)
  use CRASH_ModFermiGas, ONLY: &
       LogGe
  use ModConst

  implicit none

  SAVE

  private !Except

  !Public members
  public :: meshhv !Creates the grid of photon energies
  public :: abscon !Calculates the absorption, emission, and scattering 
  public :: opacys !Calculates opacities
  public :: nGroup, OpacityPlanck_I, OpacityRosseland_I
  public :: set_multigroup
  logical, public:: IsLogMultiGroupGrid = .true.

  !For test:
  public :: PhotonEnergy_I, AbsorptionCoefficient_I, nPhoton, EnergyGroup_I
  public :: OpacityPlanckTotal, OpacityRosselandTotal, ScatteringCoefficient_I

  !       nPhotonMax  - photon energy mesh points
  !       nGroupMax  - opacity groups     
  !       nfrqbb  - number of photon energy points near a line center   
  !       at which absorption coefficients will be computed  

  integer,parameter:: nPhotonMax = 10000, nGroupMax = 100, nfrqbb = 9

  !Radiation frequency groups:
  integer :: nGroup = nGroupMax 
  real :: EnergyGroup_I(0:nGroupMax)
  real :: DeltaLogFrequency

  public:: read_opacity_parameters
  public:: get_energy_g_from_temperature, get_temperature_from_energy_g
  public:: get_planck_g_from_temperature 
  public:: planck_opacity_integral, ross_opacity_integral

  integer,parameter :: nptspg = 30

  ! Meshs to evaluate the absorption coefficients
  integer::nPhoton
  real ::PhotonEnergy_I(nPhotonMax)

  ! The absorption coefficients
  !       AbsorptionCoefficient_I  -  array of absorption coefficients (cm**-1)            
  !       emscfs  -  array of emission coefficients (cm**-1)         
  !       ScatteringCoefficient_I  -  array of scattering coefficients (cm**-1) 


  real,dimension(nPhotonMax) ::AbsorptionCoefficient_I,ScatteringCoefficient_I

  !\
  !The opacities averaged over photon groups
  !/
  real :: OpacityPlanck_I(nGroupMax),OpacityRosseland_I(nGroupMax) 

  !\
  !the same, averaged over the whole frequency range
  !/
  real :: OpacityPlanckTotal, OpacityRosselandTotal  


  !\
  ! LOGICALS
  !/
 
  
  !Determine if we account for the corrections used in HYADES for
  !the oscillator strength

  logical :: UseHYADESCorrection4Strength = .true.

  !Switch to add or not to add the controbutions from the line core if
  !it is added somewhere else
  logical,public :: DoNotAddLineCore = .true.

  logical,public :: UseBremsstrahlung = .true.
  logical,public :: UsePhotoionization = .true.
  

  !If the logical below is set to .true. then the
  !photoionization from EXCITED states accounts for both the
  !process in which the EXCITED electron is detached and that in
  !which the VALENCE electron is
  !
  logical,public :: DoEnhancePhotoIonization =.true.
  logical,public :: UseScattering      = .true.
  logical,public :: UseAveragedRosselandOpacity = .true.

  !Correction factor for the Photoionization from the ground state:
  real,public    :: BoundFreeCorrection = 2.0
  real,public    :: TRadMin = 500   !K 
  real,public    :: ERadMin = CRadiation * 500.0**4
  real           :: TgMin_W( nGroupMax )
contains
  !============================================================================
  subroutine read_opacity_parameters

    ! Usage (with recommended values shown):
    !
    ! #OPACITY
    ! T                     UseExcitation
    ! T                     UseCoulombCorrection
    ! T                     DoStateElimination

    use CRASH_ModAtomicDataMix, ONLY: UseExcitation
    use CRASH_ModPartition,     ONLY: UseCoulombCorrection
    use CRASH_ModExcitation,    ONLY: DoStateElimination
    use ModReadParam,           ONLY: read_var

    !--------------------------------------------------------------------------
    call read_var('UseExcitation', UseExcitation)
    call read_var('UseCoulombCorrection', UseCoulombCorrection)
    call read_var('DoStateElimination',DoStateElimination)

  end subroutine read_opacity_parameters
  !============================================================================
  real function planck_opacity_integral(TeInK, OpacityPlanck_I)
    real, intent(in) :: TeInK    !Temperature in K
    real, intent(in) :: OpacityPlanck_I(nGroup)
    
    real :: Weight_I(nGroup)
    integer:: iGroup
    !------------------------------------------------------------------------
    do iGroup = 1, nGroup
       call get_planck_g_from_temperature(iGroup=iGroup, TeIn=TeInK, &
            EgSI=Weight_I(iGroup))
    end do
    planck_opacity_integral = sum(OpacityPlanck_I*Weight_I)/sum(Weight_I)
  end function planck_opacity_integral
  !===========================================================================
  real function ross_opacity_integral(TeInK, OpacityRoss_I)
    real, intent(in) :: TeInK    !Temperature in K
    real, intent(in) :: OpacityRoss_I(nGroup)
    
    real :: Weight_I(nGroup)
    integer:: iGroup
    !-------------------------------------------------------------------------
    do iGroup = 1, nGroup
       call get_planck_g_from_temperature(iGroup=iGroup, TeIn=TeInK, &
            CgSI=Weight_I(iGroup))
    end do
    ross_opacity_integral =sum(Weight_I)/sum(Weight_I/OpacityRoss_I)
  end function ross_opacity_integral
  !============================================================================
  subroutine get_planck_g_from_temperature(iGroup, TeIn, EgSI, CgSI)
    !\
    !Input parameters
    !/
    integer,intent(in):: iGroup
    real,   intent(in):: TeIn    ! electron temperature [K]

    !\
    !Output parameters
    !/
    real, optional, intent(out) :: EgSI !Radiation energy per group, J/m3
    real, optional, intent(out) :: CgSI !Radiation spec. heat/group, J/(K.m3)
    real :: xMin, xMax
    !-------------------------------------------------------------------------

    if(nGroup == 1)then
       if(present(EgSI))EgSI = cRadiation*TeIn**4
       if(present(CgSI))CgSI = 4*cRadiation*TeIn**3
       RETURN
    end if

    xMin = EnergyGroup_I(iGroup - 1)/(TeIn * cKToEV)
    xMax = EnergyGroup_I(iGroup    )/(TeIn * cKToEV)

    if(present(EgSI))EgSI = cNormG5*gint(5,xMin,xMax)*(  cRadiation*TeIn**4)
    if(present(CgSI))CgSI = cNormG6*gint(6,xMin,xMax)*(4*cRadiation*TeIn**3)

  end subroutine get_planck_g_from_temperature
  !============================================================================
  subroutine get_energy_g_from_temperature(iGroup, TgSIIn, EgSI, CgSI)
    !\
    !Input parameters
    !/
    integer,intent(in):: iGroup
    real,   intent(in):: TgSIIn    !Group temperature [K]

    !\
    !Output parameters
    !/
    real, optional, intent(out):: EgSI !Radiation energy per group, J/m3
    real, optional, intent(out):: CgSI !Radiation spec. heat/group, J/(K.m3)
    real :: xMin, xMax, TgSI
    !--------------------------------------------------------------------------

    TgSI = max(TGSIIn , TgMin_W(iGroup))

    xMin = EnergyGroup_I(iGroup - 1)/(TgSI * cKToEV)
    xMax = EnergyGroup_I(iGroup    )/(TgSI * cKToEV)

    if(present(EgSI))EgSI = cNormG5*gint(5,xMin,xMax)*(  cRadiation * TgSI**4)
    if(present(CgSI))CgSI = cNormG6*gint(6,xMin,xMax)*(4*cRadiation * TgSI**3)

  end subroutine get_energy_g_from_temperature
  !============================================================================
  subroutine get_temperature_from_energy_g(iGroup, EgSIIn, TgSIOut, CgSIOut)
    !\
    !Input parameters
    !/
    integer,intent(in):: iGroup
    real,   intent(in):: EgSIIn    !Group temperature [K]

    !\
    !Output parameters
    !/
    real, optional, intent(out):: TgSIOut !Radiation energy per group, J/m3
    real, optional, intent(out):: CgSIOut !Radiation spec. heat/group, J/(K.m3)

    real :: xMin, xMax, FreqMin, FreqMax
    real :: TgSI, CgSI, ToleranceEg, DeltaEg, EgSI

    real, parameter:: cTolerance = 1.0E-3

    integer, parameter :: nIter = 10
    integer :: iIter
    !--------------------------------------------------------------------------
    EgSI = max(EgSiIn, ERadMin)

    FreqMin = EnergyGroup_I(iGroup - 1) * cEVToK 
    FreqMax = EnergyGroup_I(iGroup    ) * cEVToK

    !\
    !Approximation to start:
    !/
    TgSi = sqrt(FreqMin*FreqMax)&
         /log(1.0 + cRadiation * FreqMin**2 * FreqMax**2 * cNormG5 * &
         log(FreqMax/FreqMin) / EgSI)

    ToleranceEg = cTolerance * EgSI

    iIter = 0
    DeltaEg = 2.0 * ToleranceEg !To start Newton-Rapson iterations
    do while (abs(DeltaEg) > ToleranceEg.and.iIter < nIter)
       xMin = FreqMin /TgSI
       xMax = FreqMax /TgSI

       DeltaEg = EgSI - cNormG5 * gint(5, xMin, xMax)*   cRadiation*TgSI**4
       CgSI    =        cNormG6 * gint(6, xMin, xMax)*(4*cRadiation*TgSI**3) 
       TgSI = TgSI + DeltaEg/CgSI

       iIter = iIter + 1
    end do


    if(present(TgSIOut))TgSIOut = TgSI
    if(present(CgSIOut))CgSIOut = CgSI

  end subroutine get_temperature_from_energy_g

  !===========================================================================

  subroutine set_multigroup(nGroupIn, FreqMinSI, FreqMaxSI, EnergyEv_I)

    ! Set the values of photon energy grid
    ! Either provide minimum and maximum frequencies,
    ! or minimum and maximum energies, 
    ! or an array of nGroup+1 energies 

    integer, optional, intent(in) :: nGroupIn ! Number of energy groups

    real, optional, intent(in):: FreqMinSI ! Minimim frequency [Hz]
    real, optional, intent(in):: FreqMaxSI ! Maximum frequency [Hz]

    real, optional, intent(in):: EnergyEv_I(:) ! Energy limits [eV]

    real:: eLogMin, eLogMax, eLog
    integer:: iGroup

    character(len=*), parameter:: NameSub = 'set_multigroup'
    !-------------------------------------------------------------------------

    if(present(nGroupIn)) nGroup = nGroupIn

    if(present(FreqMinSi)) EnergyGroup_I(0)      = FreqMinSI*cHPlanckEv 
    if(present(FreqMaxSi)) EnergyGroup_I(nGroup) = FreqMaxSi*cHPlanckEv

    IsLogMultiGroupGrid = .true.
    if(present(EnergyEv_I))then
       if(size(EnergyEv_I) == 2)then
          EnergyGroup_I(0)      = EnergyEv_I(1)
          EnergyGroup_I(nGroup) = EnergyEv_I(2)
       else if(size(EnergyEv_I) == nGroup + 1)then
          EnergyGroup_I(0:nGroup) = EnergyEv_I
          IsLogMultiGroupGrid = .false.
       else
          write(*,*) NameSub//' ERROR: size(EnergyEv_I), nGroup=', &
               size(EnergyEv_I), nGroup
          call CON_stop(NameSub// &
               ': size of EnergyEv_I should be 2 or nGroup+1')
       end if
    end if

    if(IsLogMultiGroupGrid .and. nGroup > 1)then

       ! Create a logarithmic energy grid
       !eRatio = (EnergyGroup_I(nGroup)/EnergyGroup_I(0))**(1.0/nGroup)
       !DeltaLogFrequency = log(eRatio)
       !do iGroup = 1, nGroup-1
       !   EnergyGroup_I(iGroup) = EnergyGroup_I(iGroup-1)*eRatio
       !end do

       eLogMin = log( EnergyGroup_I(0) )                            
       eLogMax = log( EnergyGroup_I(nGroup) )                    
       DeltaLogFrequency = (eLogMax - eLogMin)/nGroup                   
       eLog = eLogMin
       do iGroup=1, nGroup-1
          eLog = eLog + DeltaLogFrequency
          EnergyGroup_I(iGroup) = exp(eLog)
       end do

    end if

    ! Set EradMin and TgMin_W
    EradMin = cRadiation * TradMin**4 / nGroup
    do iGroup = 1, nGroup
       call get_temperature_from_energy_g(iGroup, ERadMin, TgMin_W(iGroup))
    end do

  end subroutine set_multigroup
  !======================================================================
  real function oscillator_strength(nI,nF)
    integer,intent(in):: nI,nF

    ! ... oscillator strength taken from Zeldovich      
    !     & Raizer for the upward transition from "ni" to "nf")             
    !      oscstr(ni,nf) = 1.96 / ni**5 / nf**3 / (1./ni**2-1./nf**2)**3
    !Here 1.96 = 32/(3\pi\sqrt{3})

    !Corrections used in HYADES 
    !f(1,2) = 0.4246 
    !f(1,3) = 0.0808 
    !f(1,4) = 0.0296 
    !f(1,5) = 0.1418 
    !f(2,3) = 0.6500 
    !f(2,4) = 0.1214 
    !f(2,5) = 0.0452 
    !f(3,4) = 0.8580 
    !f(3,5) = 0.1530 
    !f(4,5) = 1.058 
    real,parameter,dimension(1:4,2:5):: HYADESCorrection_II = reshape( (/&
         0.4246 , 0.0    , 0.0    , 0.0,   &  !nF=2
         0.0808 , 0.6500 , 0.0    , 0.0,   &  !nF=3
         0.0296 , 0.1214 , 0.8580 , 0.0,   &  !nF=4
         0.0142 , 0.0452 , 0.1530 , 1.058  &  !nF=5
         /),(/4,4/))
    !    nI =1  ! nI=2   ! nI=3   ! nI=4   !

    real :: rNI, rNF, rNI2, rNF2
    !---------------------------
    if(UseHYADESCorrection4Strength &
         .and. nI<=4 .and. nF<=5) then
       oscillator_strength = HYADESCorrection_II(nI,nF)
    else
       rNI = real(nI)
       rNF = real(nF)
       rNI2 = rNI * rNI
       rNF2 = rNF * rNF

       oscillator_strength = 1.96 * rNI * (rNF * rNF2)/ (rNF2 - rNI2)**3 
    end if
  end function oscillator_strength



  !==========================================
  subroutine meshhv  
    !                                                                       
    ! ... set a mesh of photon energies at which 
    !     we would like to evaluate  
    !     absorption coefficients                                           
    !   


    ! ... input variables:                                                  
    !       Te      -  plasma temperature (eV)                              
    !       densnn  -  number density of all nuclei (cm**-3)                
    !       densne  -  electron density (cm**-3)                            
    !       EnergyGroup_I  -  photon energy group boundaries (nGroup+1) (eV)       
    !       nGroup  -  number of photon energy groups                       
    !       nptspg  -  minimum number of mesh points per energy group       
    !

    real :: DensNN, DensNE

    real,parameter:: dPlus =  1.001, dMinus = 0.999
    !\
    !Loop variables
    !/
    integer :: iGroup    !Enumerates energy groups
    integer :: iSubGrid  !Photon energy grid for each energy group
    integer :: nMax      !

    integer :: iMix      !Over the mixture components
    integer :: iZ        !Ion charge 
    integer :: iN        !Principal quantum number
    integer :: iNUpper   !Principal quantum number for the excited state

    integer :: iPointPerLine

    !Enumerate the energies of the photoionization transitions
    !from inner shells
    integer :: nCore 

    !\
    !Grid parameters, for each group
    !/
    real :: HvMin, HvMax, dLogHv, LogHv

    !CutOff energy (for h\nu=\hbar\omega_{pe}
    real :: hvCut


    ! Averaged energy for transitions between the states with no
    ! change in the principal quantum number

    real :: ennp

    real::    Edge !Cutoff energy for photoionization from the ground state
    integer:: nBound !The number of bound electrons

    !A gap between the line centered frequency and two neighboring 
    !frequency points used to represent its profile.
    real:: dHNu
    !A width of the line
    real :: Gamma

    !The principal quantum number of the ion in the
    !ground state
    integer :: nGround 

    !Misc
    real::Dum1,Dum2,avoigt,dnudop !non-used output parameters
    real :: KK
    !------------------------------------                                                                 


    ! ... initialize variables     

    nPhoton = 0                                                        

    PhotonEnergy_I( 1:nPhotonMax ) = 0.0

    DensNN = Na * 1.0e-6
    DensNe = DensNN * zAv


    ! ... set up initial grid within each photon energy group
    do iGroup = 0,nGroup-1                                                

       hvmin = EnergyGroup_I(iGroup)                                             
       hvmax = EnergyGroup_I(iGroup+1)                                           
       dloghv = log( hvmax/hvmin ) / nptspg                           
       LogHv = log( hvmin ) - dloghv                                  
       if ( iGroup.lt.nGroup - 1 ) then                                       
          nmax = nptspg                                               
       else                                                           
          nmax = nptspg + 1                                           
       endif

       do  iSubGrid = 1, nmax                                                
          LogHv = LogHv + dloghv                                      
          PhotonEnergy_I(nPhoton+1) = exp( LogHv )                              
          nPhoton = nPhoton + 1                                           
       end do

    end do

    ! ... add 2 points near the plasma cutoff frequency                     

    hvcut = sqrt( densne / 7.25e20 )                                  
    if ( hvcut.gt.EnergyGroup_I(0) .and. hvcut.lt.EnergyGroup_I(nGroup) ) then    
       PhotonEnergy_I(nPhoton+1) = hvcut * dminus                               
       PhotonEnergy_I(nPhoton+2) = hvcut * dplus                                
       nPhoton = nPhoton + 2                                              
    endif

    ! ... add points near line centers and ionization edges                 

    do iMix=1,nMix                                              
       if ( Concentration_I(iMix) .lt. con(2) ) CYCLE                     

       do  iZ=iZMin_I(iMix), min(iZMax_I(iMix), nZ_I(iMix) - 1)                                           
          if ( Concentration_I(iMix)*Population_II(iZ,iMix) .lt. con(3) )CYCLE

          ! ...       find the principal quantum number 
          !           of the valence electrons  
          !           in their ground state                                       
          nbound = nZ_I(iMix) - iZ                               
          nGround = n_ground(iZ,nZ_I(iMix))                                     

          do iN = nGround, nExcitation_II(iZ,iMix)-1     

             if (Concentration_I(iMix) * &
                  Population_II(iZ,iMix) * &
                  Partition_III(iN, iZ, iMix)& 
                  < con(4) )CYCLE
             if(.not.DoNotAddLineCore)then        
                do  iNUpper = iN + 1, nExcitation_II(iZ,iMix)                              

                   ennp = ExcitationEnergy_III(iNUpper, iZ,iMix)- &
                        ExcitationEnergy_III(iN, iZ,iMix)
                   if(ennp<=0.0)then
                      write(*,*)'ennp<0',ennp
                      write(*,*)'iZ,iMix',iZ,nZ_I(iMix)
                      write(*,*)'Upper:',inUpper, ExcitationEnergy_III(iNUpper, iZ,iMix)
                      write(*,*)'Low:',iN, ExcitationEnergy_III(iN, iZ,iMix)
                      call CON_stop('')
                   end if
                   if ( ennp .gt. EnergyGroup_I(0) .and. &                     
                        ennp .lt. EnergyGroup_I(nGroup) ) then              

                      call line_width ( Te,densnn,ennp,cAtomicMass_I(nZ_I(iMix)), &      
                           gamma,avoigt,dnudop ) 

                      ! ...set points at and near the line center           
                      dhnu = con(9) * 4.14e-15 * gamma / 12.57         
                      PhotonEnergy_I(nPhoton+1) = ennp                           
                      nPhoton = nPhoton + 1                                

                      kk = 1.0
                      do  iPointPerLine = 1, nfrqbb / 2 - 1                                                                  
                         PhotonEnergy_I(nPhoton+1) = ennp - kk * dhnu             
                         PhotonEnergy_I(nPhoton+2) = ennp + kk * dhnu             
                         nPhoton = nPhoton + 2                              
                         kk = kk * 2.0
                      end do

                      ! ... stop if too many mesh pts are requested          
                      if ( nPhoton .gt. nPhotonMax - 5 - nfrqbb )&           
                           stop ' too many pts in -meshhv-'              
                   end if
                end do !iNUpper

             end if
             ! ... add 2 points near the ionization edge                    
             edge = IonizPotential_II(iZ+1,iMix) - ExcitationEnergy_III(iN,iZ,iMix)   

             if ( edge.gt.EnergyGroup_I(0) .and. edge.lt.EnergyGroup_I(nGroup) )&  
                  then                                                  
                PhotonEnergy_I(nPhoton+1) = edge * dminus                       
                PhotonEnergy_I(nPhoton+2) = edge * dplus                        
                nPhoton = nPhoton + 2                                     
             endif

          end do   !Over iN 

       end do    !Over iZ

       ! ... add 2 points near photoionization edges of core electrons     
       if(.not.UseCoreElectron)CYCLE
       do ncore=1,nZ_I(iMix)-1                                   
          edge = IonizPotential_II(nZ_I(iMix)+1-ncore,iMix)                        
          if ( edge.gt.EnergyGroup_I(0) .and. &                              
               edge.lt.EnergyGroup_I(nGroup) ) then                        
             PhotonEnergy_I(nPhoton+1) = edge * dminus                          
             PhotonEnergy_I(nPhoton+2) = edge * dplus                           
             nPhoton = nPhoton + 2                                        
          endif
       end do
    end do !Over iMix

    ! ... now, arrange the photon energies in monotonically increasing      
    !     order                                                             
    call sort 
  contains
    !==============
    subroutine sort 
      integer:: nCycle, j, i, iMinus1
      real :: rSave
      logical:: DoRepeat    
      !-----------------
      ncycle = 0                                                        

      ! ... first, throw out those values < or = zero                         

      j = 0                                                             
      do  i=1,nPhoton                                                    
         if ( PhotonEnergy_I(i)<= 0.0 ) CYCLE                                     
         j = j + 1                                                   
         PhotonEnergy_I(j) = PhotonEnergy_I(i)                                             
      end do
      nPhoton = j                                                          
      DoRepeat = .true.                                                                
      do while(DoRepeat)
         DoRepeat = .false.
         do i = 2,nPhoton                                                 
            iMinus1 = i - 1                                                 
            if ( PhotonEnergy_I(iMinus1) .gt. PhotonEnergy_I(i) ) then                            
               rSave = PhotonEnergy_I(iMinus1)                                          
               PhotonEnergy_I(iMinus1) = PhotonEnergy_I(i)                                        
               PhotonEnergy_I(i) = rSave                                            
               DoRepeat = .true.                                           
            endif
         end do
         nCycle = nCycle + 1                                           

         ! ...    check for error                                                

         if ( ncycle .gt. nPhoton ) call CON_stop( ' error in -sort-')                                                                                                  
         ! ...    if any of the elements have been swapped, try again
      end do
    end subroutine sort
    !===============
  end subroutine meshhv
  !====================

  subroutine abscon  
    use CRASH_ModPartition, ONLY: Z2_I

    ! ... this routine calculates the absorption, emission, and scattering 
    !     coefficients for an array of photon energies                      
    !                                                                      
    ! ... input variables:                                                  
    !       Te      -  plasma temperature (eV)                              
    !       densnn  -  number density of all nuclei (cm**-3)                
    !       densne  -  electron density (cm**-3)                            
    !       PhotonEnergy_I  -  array of photon energies at which the absorption
    !                  coefficient will be evaluated (eV)                   
    !       nPhoton   -  number of elements in the "PhotonEnergy_I" array   
    !                                                                      
    ! ... output variables:                                                 
    !       AbsorptionCoefficient_I  -  array of absorption coefficients (cm**-1)            
    !       emscfs  -  array of emission coefficients (cm**-1)         
    !       ScatteringCoefficient_I  -  array of scattering coefficients (cm**-1)          
    !                                                                      

    real:: DensNN, DensNe    ![cm -3]

    !\
    !Arrays to collect contributions from different effects for a given component
    real,dimension(nMixMax) :: BremsStrahlung_I = 0.0
    real,dimension(nMixMax) :: PhotoIonizationAbs_I = 0.0


    !\
    ! Loop variables
    !/
    integer :: iMix, & !runs over the mixture component
         iZ,   & !runs over the charge number
         iN,   & ! runs over the quntum principal number for the lower level
         nBound, &!Number of bound electrons
         nGround,&!For a given iZ, the principal number of the ground state.
         iPhoton

    real :: eTransition, ETransitionPerTe

    !\
    !Partial sums
    !/
    real:: SumOverZ, SumOverN 

    !\
    ! Bound-bound contributions
    !/
    real :: abslns

    !\
    ! Cutoff energy for the photon with the frequency correspondent to
    ! the plasma electron frequency
    !/
    real :: HNuCut

    !Misc: functions of the photon energy
    real :: ExpMinusHNuPerT,  PhotonEnergyCubed

    !Misc: Density times Z2
    real :: DensityZ2

    !Misc: coefficients used in free-free gaunt factor
    real :: Log10OfGamma2, GauntFactorBrems

    !Number of screened and valence electrons, as well as iZEff
    integer :: nValence, iQ

    ! densities to calculate photoionization
    real :: AbsorberDensity 


    !Misc: Scattering parameters
    real :: ScatNe, tScatt, pScatt

    !Misc: Photoionization parameters
    real :: PhotoIonizationConst, BoundFreeGauntFactor

    !---------------------------------------------------

    !Initialization

    AbsorptionCoefficient_I  = 0.0  !-  array of absorption coefficients (cm**-1)                    
    ScatteringCoefficient_I  = 0.0  !-  array of scattering coefficients (cm**-1)     

    DensNN = Na * 1.0e-6
    DensNe = DensNN * zAv
   

    ! ... loop over photon energies
    do iPhoton = 1,nPhoton                                              
       !\
       !The dependence of the cross-sections \propto (Photon Energy)^{-3} is typical
       !Therefore, introduce:
       !/                           
       PhotonEnergyCubed = PhotonEnergy_I(iPhoton)**3

       !\
       !The factor needed to account for the stimulated emission:
       !/
       ExpMinusHNuPerT = exp( -PhotonEnergy_I(iPhoton) / Te )                                    

       ! ...    loop over gas species                                          

       IMIXLOOP: do iMix = 1,nMix                                           

          BremsStrahlung_I(iMix) = 0.0                                            


          ! ...       Bremsstrahlung                                              
          !           --------------                                              

          if(Z2_I(iMix) > 1.0e-30)then

             DensityZ2 = Z2_I(iMix) * Concentration_I(iMix) * densnn * 1.e-16               

             ! ...       the free-free gaunt factor is a simple fit to the results   
             !           of Karzas and Latter (Ap. J. Suppl., 6, 167 (1961))         
             
             !
             Log10OfGamma2 = log10( 13.6 * Z2_I(iMix) / Te )                           
             GauntFactorBrems  = 1. + 0.44 * exp( -0.25*(Log10OfGamma2+0.25)**2 )          
             
             BremsStrahlung_I(iMix) = 2.4e-21 * DensityZ2 * GauntFactorBrems * densne * &          
                  (1.0 - ExpMinusHNuPerT) / ( sqrt( Te ) * PhotonEnergyCubed )
             !!NOTE: there is a more accurate calculation of the Gaunt factor in 
             !ggff.f file from HYADES
          end if

          ! ...       photoionization                                             
          !           ----------------                                            

          ! ...       sum over ionization levels for the transition from          
          !           "iZ" to "iZ+1"                                      

          PhotoIonizationAbs_I(iMix) = 0.0                                            
          SumOverZ = 0.0                                                                        

          IZLOOP: do iZ = iZMin_I(iMix), min( iZMax_I(iMix), nZ_I(iMix) - 1)                               


             ! ...          find the principal quantum number of the valence electron
             !              in the ground state for the ion before ("nprin0");       
             !              "nbound" is the number of electrons bound to the ion     
             !              after another is captured (or, before ionization).       

             nBound = nZ_I(iMix) - iZ                            
             nGround = n_ground( iZ, nZ_I(iMix) )                                  

             ! ...          first, consider the contibution from valence shell       
             !              electrons                                                

             SumOverN = 0.0                                               

             if ( Concentration_I(iMix) * Population_II(iZ,iMix) < con(3) ) &
                  CYCLE IZLOOP     

             ! ...            sum over quantum states                                
             do  iN = nGround, nExcitation_II(iZ,iMix)  

                AbsorberDensity = Population_II(iZ,iMix) * Partition_III(iN, iZ, iMix)                      

                !  calculate the energy to excite the outermost electron 
                !  into  the continuum 

                eTransition = IonizPotential_II(iZ+1, iMix) - &
                     ExcitationEnergy_III(iN, iZ, iMix)         

                ! ...  the photon energy must exceed the binding energy 
                !Acc ount for the ionization potential lowering while 
                !treating the photoionization     

                if ( PhotonEnergy_I(iPhoton) >=  eTransition &
                     .or. (iN==nGround .and.                 &
                     PhotonEnergy_I(iPhoton) >=  eTransition &
                      - IonizationPotentialLowering_I(iZ) )) then  


                   ! ... find the number of "screening" electrons     
                   !     and the number of electrons in the outermost shell   
                   !                 
                   
                   if ( iN == nGround ) then                          
                      ! ...ground state ion                                                              
                      nValence = nBound - n_screened(iZ, nZ_I(iMix))
                      BoundFreeGauntFactor = BoundFreeCorrection
                   else                                                 
                      ! ...ion is excited                                                                    
                      nValence = 1    
                      BoundFreeGauntFactor = 1.0                                     
                   endif
                   
                   ETransitionPerTe = eTransition / Te                                   

                   !\
                   ! The version more close to that implemented in
                   ! Emilio Minguez et al. With this approach we substitute
                   ! eTransition for the combination Ry * Z_eff^3/iN^2, which
                   ! is implied in the ionmix.f.
                   ! By this account, our factor PhotoIonizationConst
                   ! differs from the version in the ionmix by a foctor of
                   ! (1/13.6)**3 where 13.6 = cRyEV.
                   ! The last multiplier accounts for effects of the Fermi
                   ! statistics in an electron gas. The formula is available in
                   ! PhotoIonization.pdf document. Note, that our LogGe = -\mu/T
                   !/         
                   SumOverN = SumOverN + nValence * iN/(real(iZ+1)**2) * AbsorberDensity&
                        * BoundFreeGauntFactor &
                        * (1.-ExpMinusHNuPerT) * (eTransition**3)/&
                        (1.0 + ExpMinusHNuPerT*exp(ETransitionPerTe - LogGe))
                end if
                !  calculate the energy to excite the valence electron 
                !  into  the continuum 
                if(DoEnhancePhotoIonization .and. &
                     iZ <= nZ_I(iMix) - 2 .and. &
                     iN > nGround) then
                   eTransition = eTransition + &
                        ExcitationEnergy_III(iN, iZ+1, iMix)         

                   ! ...  the photon energy must exceed the binding energy     

                   if ( PhotonEnergy_I(iPhoton) <  eTransition ) CYCLE  


                   ! ... find the number of "screening" electrons     
                   !     and the number of electrons in the outermost shell   
                   !                 
                   
                                             
                   nValence = nBound - 1 - n_screened(iZ+1, nZ_I(iMix))                         
                  
                   
                   ETransitionPerTe = eTransition / Te                                   


                   !\
                   ! The version more close to that implemented in
                   ! Emilio Minguez et al. With this approach we substitute
                   ! eTransition for the combination Ry * Z_eff^3/iN^2, which
                   ! is implied in the ionmix.f.
                   ! By this account, our factor PhotoIonizationConst
                   ! differs from the version in the ionmix by a foctor of
                   ! (1/13.6)**3 where 13.6 = cRyEV.
                   ! The last multiplier accounts for effects of the Fermi
                   ! statistics in an electron gas. The formula is available in
                   ! PhotoIonization.pdf document. Note, that our LogGe = -\mu/T
                   !/         
                   SumOverN = SumOverN + nValence * n_ground(iZ+1,nZ_I(iMix))&
                        /(real(iZ+2)**2) * AbsorberDensity&
                        * BoundFreeCorrection &
                        * (1.-ExpMinusHNuPerT) * (eTransition**3)/&
                        (1.0 + ExpMinusHNuPerT*exp(ETransitionPerTe - LogGe))
                end if


             end do


             ! ...          now, add core electron photoionization cross-sections    
             !              to the absorption term                                   

             ! ...          loop over inner shells (K,L,M,...); each inner shell     
             !              is assume to be full                                     


             if ( & !nshels > 0 .and. 
                  UseCoreElectron ) then
                call CON_stop('UseCoreElectron should be set to .false.')
                !  do ishell=1,nshels                                 

                ! ... determine the photoionization cutoff energy (eV); use
                !     the ionization potential of the outermost bound      
                !     electron; "nocc" is the number of electrons occupying
                !     shell "ishell", "nscren" is the number of electrons  
                !     screening shell "ishell"                             

                !    nocc = 2*ishell*ishell                               
                !    nscren = nscrsh( ishell )                            
                !    izeff = izgas(iMix) - nscren                         
                !    enpi = pot(iMix,izgas(iMix)+1-nscrsh(ishell+1))      
                !    if ( PhotonEnergy_I(iPhoton) .ge. enpi ) then                  
                !       sum1ac = sum1ac + nocc * izeff**4 / ishell**5      
                !    endif

                ! end do

                ! sum1a = sum1a + sum1ac * Population_II(iZ,iMix)             

             endif

             SumOverZ = SumOverZ + SumOverN                                                                       

          end do IZLOOP
          !\
          ! Commented out version from the ionmix:
          ! const = (1.99e-14*densnn) * Concentration_I(iMix) / PhotonEnergyCubed
          ! Insetad we use here the expression from Zel'dovich and Raizer, Eq.5.34,
          ! which differs by a factor of (1/13.6)**3, es we explained above.
          ! The exact espression for the photoionization cross-section is:
          ! 7.9e-18 cm^2= \frac{64\pi}{3\sqrt{3}}\alpha a_0^2, where \alpha is the fine
          ! structure constant and a_0 is the Bohr radius
          !/
          PhotoIonizationConst = 7.9e-18 * Densnn * Concentration_I(iMix) / PhotonEnergyCubed   
          PhotoIonizationAbs_I(iMix) = PhotoIonizationConst * SumOverZ

       end do IMIXLOOP   !Over iMix                                                                    


       ! ...    Scattering contributions                                       
       !        ------------------------                                       

       ! ...    Thomson scattering contribution                                
       !        first, find the "effective" electron density for scattering;   
       !        i.e., if the photon energy is greater than the binding energy  
       !        of bound electrons, include bound electrons to the density.    
       scatne = 0.0                                                  
       do  iMix=1,nMix                                           
          iq = 0                                                      
          do  iZ = 0,nZ_I(iMix) -1                              

             if ( PhotonEnergy_I(iPhoton) > IonizPotential_II(iZ+1, iMix)) then            
                iq = iq + 1                                           
                scatne = scatne + Concentration_I(iMix)
             else
                EXIT
             endif
          end do
          !iQ is the last state at which the ionization potential 
          !is less than the photon energy
          !All the other state give only the free electron contribution to
          !

          do iZ= max(iZMin_I(iMix),iQ+1), iZMax_I(iMix)
             scatne = scatne + (iZ-iq) * Population_II(iZ,iMix)*&                    
                  Concentration_I(iMix) 
          end do
       end do
       scatne = scatne * densnn                                       
       tscatt = 6.66e-25 * scatne !Thomson cros-section                                     

       ! ...    contribution from plasma oscillations                          
       hnucut = sqrt( densne / 7.25e20 )      

       if ( PhotonEnergy_I(iPhoton) .lt. hnucut ) then                          
          pscatt = 5.05e4 * sqrt( hnucut**2 - PhotonEnergy_I(iPhoton)**2 )      
       else                                                           
          pscatt = 0.                                                 
       endif

       if ( UseScattering) ScatteringCoefficient_I(iPhoton) = tscatt + pscatt            

       !\
       !Sum up al the contributions                                                               
       AbsorptionCoefficient_I(iPhoton) = 0.0                                                                                      
       if(UseBremsstrahlung)then                                         
          AbsorptionCoefficient_I(iPhoton) = AbsorptionCoefficient_I(iPhoton)+sum(BremsStrahlung_I(1:nMix)) 
       end if
       if(UsePhotoionization)then

          AbsorptionCoefficient_I(iPhoton) = AbsorptionCoefficient_I(iPhoton)+sum(PhotoIonizationAbs_I(1:nMix))  
       end if



       ! ...    line contributions                                             
       !        ------------------                                             

       ! ...    add in the contribution from bound-bound transitions           

       call lines ( PhotonEnergy_I(iPhoton), abslns)      

       AbsorptionCoefficient_I(iPhoton) = AbsorptionCoefficient_I(iPhoton) + abslns                                                            


    end do  !over iPhoton
  contains
    !======================================================================
    subroutine lines (PhotonEnergy,AbsorptionBB)       

      !                                                                       
      ! ... compute the contribution to the absorption coefficient from       
      !     all lines (bound-bound transitions)                               
      !                                                                       
      ! ... input variables:                                                 
      !       Te      -  plasma temperature (eV)                             
      !       densnn  -  number density of all nuclei (cm**-3)               
      !       densne  -  electron density (cm**-3)                            
      !       PhotonEnergy   -  photon energy (eV)                                   
      real,intent(in) :: PhotonEnergy
      real :: densnn,densne
      !                                                                       
      ! ... output variables:                                                 
      !       AbsorptionBB  -  absorption coefficient due to lines (cm**-1)                    
      !   

      real,intent(out) :: AbsorptionBB                            

      integer:: iSav !Integer to count the total number of lines involved
      integer,parameter:: nSavMax = 200 !The upper bound for iSav

      !\
      ! Loop variables
      !/
      integer :: iMix, & !runs over the mixture component
           iZ,   & !runs over the charge number
           iN,   & ! runs over the quntum principal number for the lower level
           iNUpper,& !the same, for upper level 
           nBound, &!Number of bound electrons
           nGround  !For a given iZ, the principal number of the ground state.

      real:: denlq  !Density of ions of a given sort, [cm-3]
      real:: denlqn !The same, for the ion in the lower state, for a given transition


      !The controbutions to the total absorption,
      !calculated by the 'abslin' subroutine
      real:: abscof

      real:: TransitionEnergy
      real:: gamma,avoigt,dnudop,vvoigt  !Line width parameters

      !Difference between the line center and the photon energy, 
      !for which to calculate the absorption

      real:: DeltaNu 

      !The line height at a given frequency (the line shape)
      real:: Shape  
      !-------------------


      iSav = 0                                                          
      AbsorptionBB = 0.                                                       


      DensNN = Na * 1.0e-6  !To convert to cm-3
      DensNe = DensNN * zAv !Electron density, in cm-3

      ! ... loop over gas species                                             
      IMIXLOOP: do iMix =1,nMix                                              
         if (Concentration_I(iMix) .lt. con(2) ) CYCLE IMIXLOOP

         ! ...    loop over ionization states                         
         IZLOOP: do  iZ=iZMin_I(iMix), min(iZMax_I(iMix), nZ_I(iMix) - 1)

            if ( Concentration_I(iMix)*Population_II(iZ,iMix) .lt. con(3))&
                 CYCLE IZLOOP

            ! ... find the principal quantum number of the valence electrons  
            !           in their ground state                                       

            nGround = n_ground(iZ,nZ_I(iMix))                                   

            denlq  = densnn * Concentration_I(iMix) * Population_II(iZ,iMix)          

            ! ...  loop over initial quantum states                            
            INLOOP: do iN=nGround,nExcitation_II(iZ,iMix)-1                              
               if ( Concentration_I(iMix) * &
                    Population_II(iZ,iMix)* &
                    Partition_III(iN,iZ,iMix)  &
                    <  con(4) ) CYCLE  INLOOP                           

               ! ... compute the density of ions with an electron in          
               !             in state "n"                                             
               denlqn = denlq * Partition_III(iN,iZ,iMix)                  

               !         loop over final quantum states      

               do iNUpper=iN+1,nExcitation_II(iZ,iMix)                                 

                  TransitionEnergy = ExcitationEnergy_III(iNUpper,iZ,iMix) - ExcitationEnergy_III(iN,iZ,iMix)                  


                  ! ... compute the line widths for natural, Doppler, and pressure        
                  !     broadening to be used with Lorentzian line shape 
                  if(TransitionEnergy<=0.0)then
                      write(*,*)'TransitionEnergy<0'
                      call CON_stop('')
                   end if
                  call line_width ( Te, densnn, TransitionEnergy, cAtomicMass_I(nZ_I(iMix)),  &                            
                       gamma, avoigt, dnudop )        

                  ! ... compute this contribution if the photon energy is not far         
                  !     from the line center   

                  DeltaNu = ( TransitionEnergy - PhotonEnergy ) / cHPlanckEV                                  
                  if ( abs( deltaNu ) > con(5)*gamma ) then 
                     !Line is too far, ignore it
                     abscof = 0.                                                                                                       
                     CYCLE
                  end if

                  ! ...    if the contribution from line "cores" will be added            
                  !        elsewhere (-opacbb-), then use the value at the                
                  !        line core boundary.                                            

                  if (DoNotAddLineCore .and. abs(DeltaNu).lt.con(6)*gamma ) then     
                     DeltaNu = con(6)*gamma                                       
                  endif

                  if(UseVoigt)then
                     ! ...       use Voigt profile                                           
                     vvoigt = deltaNu / dnudop                                    
                     shape = voigt_profile ( avoigt,vvoigt ) / dnudop / 1.7725           
                  else                                                           
                     ! ...       use Lorentzian profile                                      
                     shape = (gamma/39.48) / (DeltaNu**2 + (gamma/12.57)**2)      
                  endif



                  AbsorptionBB = AbsorptionBB + &
                       
                       2.65e-2 * oscillator_strength ( iN, iNUpper) * shape * denlqn    




               end do
            end do INLOOP
         end do IZLOOP
      end do IMIXLOOP

      !Correct for stimulated emission)
      AbsorptionBB = AbsorptionBB * (1.0 - exp( -PhotonEnergy/Te ))

      !====================================
    end subroutine lines
  end subroutine abscon
  !====================



  subroutine opacys ( TRadIn )                    

    ! ... this routine computes the Planck and Rosseland opacities for      
    !     the photon energy group "EnergyGroup_I".  "opacpm" and "OpacityRosseland_Im" are the  
    !     Planck and Rosseland mean opacities (integrated over frequency).  


    real, optional, intent(in) :: TRadIn

    ! ... input variables:                                                  
    !       Te     - plasma temperature (eV)                                
    !       densnn - number density of all nuclei (cm**-3)                 
    !       densne - electron density (cm**-3)                              
    !       photen - photon energies (eV)                                   
    !       nphot  - number of elements in "photen" array                   
    !       AbsorptionCoefficient_I - array of absorption coefficients (cm**-1)              
    !       emscfs - array of emission coefficients (cm**-1)                
    !       ScatteringCoefficient_I - array of scattering coefficients (cm**-1)              
    !       EnergyGroup_I - photon energy group boundaries (eV) (nGroup+1)  
    !       nGroup - number of energy bins                                  
    !       ntrad  - number of radiation temperatures                       
    !       trad   - array of radiation temperatures (eV)                   


    ! ... output variables:                                                 
    !       OpacityPlanck_I - Planck group opacities for absorption (cm-1)        
    !       OpacityRosseland_I  - Rosseland group opacity (cm-1)                      
    !       OpacityPlanckTotal - Planck mean opacity for absorption (cm-1)                        
    !       OpacityRosselandTotal  - Rosseland mean opacity (cm-1)                       
    !       culrat - plasma cooling rate (erg cm**3/sec)                    
    !       op2tp  - Planck 2-temperature opacities (cm**2/g)               
    !       op2tr  - Rosseland 2-temperature opacities (cm**2/g)            
    !                                                                       

    real :: DensNE, DensNN
    real :: TRad
    !\
    !Loop Variables
    !/
    integer:: iGroup

    real :: XGroupMin, XGroupMax, LineCoreOpacity
    !---------------------
    DensNN = Na * 1.0e-6  !To convert to cm-3
    DensNe = DensNN * zAv !Electron density, in cm-3
    if(present(TRadIn))then
       TRad = TRadIn
    else
       TRad = Te
    end if

    ! ... compute the Planck and Rosseland opacities for each photon        
    !     energy group                                                      

    OpacityPlanckTotal    = 0.0
    OpacityRosselandTotal = 0.0
    do iGroup = 1, nGroup                                                 
       XGroupMin = EnergyGroup_I(iGroup-1) / TRad
       XGroupMax = EnergyGroup_I(iGroup  ) / TRad

       call opacgp (OpacityPlanck_I(iGroup),OpacityRosseland_I(iGroup) )  
       if(UseAveragedRosselandOpacity)&
            OpacityRosseland_I(iGroup) = sqrt(&
            OpacityPlanck_I(iGroup)*&
            OpacityRosseland_I(iGroup)       )
       if ( DoNotAddLineCore ) then                                     
          ! ...       use analytic solution to bound-bound opacities 
          call opacbb (LineCoreOpacity)
          if(LineCoreOpacity<0.0)then
             write(*,*)'iGroup =',iGroup
             write(*,*)'LineCoreOpacity =', LineCoreOpacity
             call CON_stop('Negative contribution in opacbb')
          end if
          OpacityPlanck_I(iGroup) = OpacityPlanck_I(iGroup) + LineCoreOpacity
       endif
       OpacityPlanckTotal = OpacityPlanckTotal + gint(5,XGroupMin,XGroupMax) &
            * OpacityPlanck_I(iGroup)
       OpacityRosselandTotal  = OpacityRosselandTotal  &
            + gint(6,XGroupMin,XGroupMax) / OpacityRosseland_I(iGroup)
    end do
    ! ... compute the total opacities based on the group opacities         
    OpacityPlanckTotal = OpacityPlanckTotal * cNormG5
    OpacityRosselandTotal = 1.0/ (OpacityRosselandTotal * cNormG6)

  contains

    subroutine opacgp ( OpacityPlanck, OpacityRosseland ) 
      !... this routine computes the Planck and Rosseland opacities for
      !     the photon energy group from "XGroupMin" to "XGroupMax".
      !
      ! ... input variables:
      !       Te     - plasma temperature (eV)                         
      !       rho    - mass density (g/cm**3)                                 
      !       trad   - radiation temperature (eV)                             
      !       PhotonEnergy_I - photon energies at which the absorption coefficients   
      !                were calculated                                        
      !       nphot  - number of photon energy mesh points                    
      !       AbsorptionCoefficient_I - array of absorption coefficients (cm**-1) (w/o scatt)  
      !       emscfs - array of emission coefficients (cm**-1) (w/o scatt)    
      !       ScatteringCoefficient_I - array of scattering coefficients (cm**-1)             
      !       XGroupMin  - minimum photon energy (in units of kT)                
      !       XGroupMax  - maximum photon energy (in units of kT)   



      !                                                                      
      ! ... output variables:                                                 
      !       OpacityPlanck - Planck opacity for absorption (cm**-1)                                
      !       OpacityRosseland  - Rosseland opacity (cm**-1)                             
      !                                                                       

      real, intent(out) :: OpacityPlanck, OpacityRosseland
      !\
      !Loop variables
      !/
      integer :: iPhot, iPhot0

      !\
      !Partial sums
      !/
      real :: dsumpa, dsmpa0, dsumr, dsumr0

      !Misc:
      real :: xg1,xg2, emx, fpa2,fpa20,fr2, fr20,sumpa,sumpa0,sumr,sumr0
      real :: fpa1, fpa10,fr1, fr10, dxg, g5, g6
      !-----------------------------                                     


      ! ... find the index where the photon energy equals the lower group     
      !     boundary (note: mesh pts should be placed at each boundary)      
      do iPhot=1,nPhoton                                                   
         if ( PhotonEnergy_I(iPhot) >= XGroupMin * Trad * 0.99999 ) then
            iPhot0 = iPhot
            EXIT
         end if
      end do

      ! ... integrate to get the group opacities                              
      !     ------------------------------------                              

      ! ... initialize values for the first point                             

      xg2 = PhotonEnergy_I(iphot0) / trad                                       
      emx = exp( -xg2 )                                                 

      fpa2 = xg2**3 * emx * AbsorptionCoefficient_I(iphot0) / ( 1.-emx )                 
      fpa20 = AbsorptionCoefficient_I(iphot0)

      fr2 = xg2**4 * emx / (1.-emx)**2 / &
           ( AbsorptionCoefficient_I(iphot0) + ScatteringCoefficient_I(iphot0) )                         
      fr20 = 1. / ( AbsorptionCoefficient_I(iphot0) + ScatteringCoefficient_I(iphot0) )                   

      ! ... loop over photon energy mesh points                              

      sumpa  = 0.0                                                      
      sumpa0 = 0.0                                                      

      sumr   = 0.0                                                      
      sumr0  = 0.0                                                      

      do iPhot = iPhot0 + 1, nPhoton

         if ( PhotonEnergy_I(iPhot) > XGroupMax * trad  * 1.00001 ) EXIT

         xg1 = xg2                                                      
         xg2 = PhotonEnergy_I(iPhot) / TRad                                     

         fpa1  = fpa2                                                   
         fpa10 = fpa20                                                  

         fr1   = fr2                                                    
         fr10  = fr20                                                   

         if ( xg2 .lt. 1.e20 ) then                                     

            emx = exp( -xg2 )                                           

            fpa2  = xg2**3 * emx * AbsorptionCoefficient_I(iphot) / ( 1.-emx )           
            fpa20 = AbsorptionCoefficient_I(iphot)                                      

            fr2   = xg2**4 * emx / (1.-emx)**2 / &                       
                 ( AbsorptionCoefficient_I(iPhot) + ScatteringCoefficient_I(iPhot) )                   
            fr20  = 1. / ( AbsorptionCoefficient_I(iphot) + ScatteringCoefficient_I(iphot) )              

            dxg = xg2 - xg1                                             

            ! ...       integrate using a logarithmic interpolation scheme   

            ! ...       Planck mean absorption opacities                     
            if ( abs( fpa1-fpa2 ) .gt. 1.e-3*fpa2 ) then                
               if ( fpa2 .ne. 0. ) then                                  
                  dsumpa = dxg * ( fpa2-fpa1 ) / log( fpa2/fpa1 )         
               else                                                      
                  dsumpa = 0.                                             
               endif
            else                                                       
               dsumpa = dxg * fpa1                                       
            endif

            ! ...       mean (non-weighted) absorption opacities                 
            if ( abs( fpa10-fpa20 ) .gt. 1.e-3*fpa20 ) then             
               if ( fpa20 .ne. 0. ) then                                 
                  dsmpa0 = dxg * ( fpa20-fpa10 ) / log( fpa20/fpa10 )     
               else                                                      
                  dsmpa0 = 0.                                             
               endif
            else                                                        
               dsmpa0 = dxg * fpa10                                      
            endif


            ! ...       Rosseland mean opacities                                 
            if ( abs( fr1-fr2 ) .gt. 1.e-3*fr2 ) then                   
               if ( fr2 .ne. 0. ) then                                   
                  dsumr = dxg * ( fr2-fr1 ) / log( fr2/fr1 )              
               else                                                      
                  dsumr = 0.                                              
               endif
            else                                                        
               dsumr = dxg * fr1                                         
            endif

            ! ...       mean (non-weighted) "transport" opacities                   
            if ( abs( fr10-fr20 ) .gt. 1.e-3*fr20 ) then                
               if ( fr20 .ne. 0. ) then                                  
                  dsumr0 = dxg * ( fr20-fr10 ) / log( fr20/fr10 )        
               else                                                      
                  dsumr0 = 0.                                             
               endif
            else                                                        
               dsumr0 = dxg * fr10                                       
            endif

            sumpa  = sumpa + dsumpa                                    
            sumpa0 = sumpa0 + dsmpa0                                    

            sumr   = sumr + dsumr                                       
            sumr0  = sumr0 + dsumr0                                     

         endif
      end do
      ! ... normalize to get the group opacities                  

      g5 = gint(5,XGroupMin,XGroupMax)                                          
      g6 = gint(6,XGroupMin,XGroupMax)                                          
      if ( XGroupMax.ge.1./con(7) .and. XGroupMin.le.con(7) .and. &              
           g5.gt.0. ) then                                              
         ! ...    weight absorption coef. by Planck function          
         OpacityPlanck = sumpa  / g5                                                                            
      else                                                             
         ! ...    use straight average of absorption coef.                       
         OpacityPlanck = sumpa0  / ( XGroupMax-XGroupMin )                        
      endif

      if ( XGroupMax.ge.1./con(7) .and. XGroupMin.le.con(7) .and. &              
           sumr.gt.0. .and. g6.gt.0. ) then                             
         ! ...    weight absorption and scattering coefs. by Planck function 
         !        derivative wrt temperature                                     
         OpacityRosseland = g6  / sumr                                        
      else if ( sumr0.gt.0. ) then                                      
         ! ...    use straight average of absorption plus scattering coefs.      
         OpacityRosseland = ( XGroupMax-XGroupMin )  / sumr0                          
      else                                                              
         ! ...    in the very high PhotonEnergy_I energy limit, this just becomes equal 
         !        to Thomson scattering contribution (assuming the absorption  
         !        is much less)                                            
         OpacityRosseland = ( ScatteringCoefficient_I(nPhoton) + AbsorptionCoefficient_I(nPhoton) )                 
      endif


    end subroutine opacgp
    !=============================
    subroutine opacbb (OpacityPlanck)  
      !                                                                       
      ! ... compute the contribution to the group opacity from                
      !     all lines (bound-bound transitions)                               
      !                                                                       
      ! ... input variables:                                             
      !       Te      -  plasma temperature (eV)                    
      !       densnn  -  number density of all nuclei (cm**-3)       
      !       densne  -  electron density (cm**-3)                
      !       trad    -  radiation temperature (eV)                  
      !       XGroupMin   -  minimum photon energy of group (in units of kT) 
      !       XGroupMax   -  maximum photon energy of group (in units of kT) 



      !                                                                
      ! ... output variables:                                  
      !       OpacityPlanck   -  Planck absorption due to lines (cm^-1)                   
      !       OpacityRosseland  -  Rosseland opacity due to lines (cm^)            
      !---------------------------------------------                                                                

      real, intent(out) :: OpacityPlanck           


      !\
      ! Loop variables
      !/
      integer :: iMix, & !runs over the mixture component
           iZ,   & !runs over the charge number
           iN,   & ! runs over the quntum principal number for the lower level
           iNUpper,& !the same, for upper level 
           nBound, &!Number of bound electrons
           nGround  !For a given iZ, the principal number of the ground state.

      real :: TransitionEnergy, Energy2TRadRatio 
      real :: ExpOfEnergy2TRadRatio
      !Misc:
      real :: g5, const, const1, opacij, fnnp


      !----------------------------


      OpacityPlanck = 0.0                                                                                                           

      g5  = gint(5,XGroupMin, XGroupMax )                                      

      if ( g5 < 1e-30 ) return

      const = 1.10e-16 * DensNN / ( Te * g5 )                     

      ! ... loop over gas species                                             
      IMIXLOOP: do iMix = 1,nMix                                              
         if ( Concentration_I(iMix) .lt. con(2) ) CYCLE IMIXLOOP

         ! ...    loop over ionization states                                  
         IZLOOP:   do iZ = iZMin_I(iMix), min( iZMax_I(iMix), nZ_I(iMix) - 1)                                                                              
            if ( Concentration_I(iMix) * Population_II(iZ,iMix) .lt. con(3)) &
                 CYCLE IZLOOP  

            ! ...       find the principal quantum number of the valence electrons  
            !           in their ground state                                                                            

            const1 = const * Concentration_I(iMix) * Population_II(iZ,iMix)       


            ! ...       loop over initial quantum states 

            INITIAL: do iN = n_ground(iZ, nZ_I(iMix)) , nExcitation_II(iZ,iMix) -1
               if ( Concentration_I(iMix) * Population_II(iZ,iMix) * Partition_III(iN, iZ, iMix)&
                    .lt. con(4) ) CYCLE INITIAL                             

               ! ... calculate opacity due to "delta n" = 0 transitions       

               !call bbneq0 ( Te,izgas(lgas),nbound,&
               !     enn,fnn,gnn )        
               !if ( isw(11).ne.0 ) enn = 0.                             
               !Energy2TRadRatio = enn / trad                                          
               !if ( Energy2TRadRatio.gt.XGroupMin .and. Energy2TRadRatio.le.XGroupMax ) then                
               !   ExpOfEnergy2TRadRatio = exp( -Energy2TRadRatio )                                      
               !   opacnn = const1 * fnn * fraclv(lgas,izp1,n) *&
               !        ExpOfEnergy2TRadRatio * Energy2TRadRatio**3
               ! ...  correct for the "effective" stimulated emission to    
               !      give the proper form for the cooling rate             
               !xnjoni = gnn * densne * ExpOfEnergy2TRadRatio / &
               !   ( gnn * densne + enn**3 * sqrt(Te) * 2.74e12 ) 
               !   corsee = xnjoni / ExpOfEnergy2TRadRatio                                 
               !   corsea = ( 1.-xnjoni ) / ( 1.-ExpOfEnergy2TRadRatio )                   
               !   OpacityPlanck = OpacityPlanck + opacnn * corsea                       
               !   opems = opems + opacnn * corsee                      
               !endif
               ! ...          loop over final quantum states

               FINAL: do iNUpper = iN+1, nExcitation_II(iZ,iMix)                               

                  ! ...            compute the transition energy     
                  TransitionEnergy = ExcitationEnergy_III(iNUpper, iZ,iMix)- &
                       ExcitationEnergy_III(iN, iZ,iMix)
                  Energy2TRadRatio = TransitionEnergy / TRad                                       
                  if ( Energy2TRadRatio .le. XGroupMin .or. Energy2TRadRatio .gt. XGroupMax ) CYCLE FINAL    

                  ExpOfEnergy2TRadRatio = exp( -Energy2TRadRatio )                                      
                  fnnp = oscillator_strength(iN, iNUpper)

                  opacij = const1 * fnnp * Partition_III(iN, iZ, iMix)  * &
                       ExpOfEnergy2TRadRatio * Energy2TRadRatio**3                                  

                  !\
                  !Account for the stimulated emission
                  !/
                  ! ...                             


                  OpacityPlanck = OpacityPlanck  +  &
                       opacij * ( 1.-ExpOfEnergy2TRadRatio )

               end do FINAL
            end do INITIAL
         end do IZLOOP
      end do IMIXLOOP
    end subroutine opacbb
  end subroutine opacys
  !====================
end module CRASH_ModMultiGroup

