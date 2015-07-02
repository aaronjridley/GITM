!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module CRASH_ModPartition
  use CRASH_ModAtomicDataMix
  use CRASH_ModExcitation,ONLY : &
       get_excitation_levels,&
       get_excitation_levels_zero,&
       nMixMax, LogGi_II,&
       ExtraEnergyAv_II, VirialCoeffAv_II, Cov2ExtraEnergy_II,&
       CovExtraEnergyVirial_II, Cov2VirialCoeff_II, SecondVirialCoeffAv_II
  use ModConst
  implicit none
  SAVE
  PRIVATE !Except
  logical :: DoInit = .true.             !Set to false after the initialization

 

  !\
  !For each component in the mixture:
  !/

  integer,dimension(nMixMax) :: iZMin_I  ! Numbers of the ionization states, such that the population
  integer,dimension(nMixMax) :: iZMax_I  ! of ion states with iZ<iZMin or iZ>iZMax is negligible.
  public:: iZMin_I,iZMax_I   !Open an access to these data, to optimize the spectroscopic calculations

  ! Array of energies needed to create i-level ion from a neutral atom
  real,dimension(0:nZMax,nMixMax) :: IonizEnergyNeutral_II 

  !\
  !The distribution over the ion states, for each of the elements:
  !/
  real,dimension(0:nZMax,nMixMax) :: Population_II ! Array of the populations of ions


  !The combination of fundamental constants, used to calculate the electron
  !statistical weight: electron statistical weight at the temperature of 1eV
  !and the conscentration of 1 particle per m^3
  ! eWeight1eV1m3 = 2/(Lambda^3)
  ! where 1/Lambda = sqrt(cTwoPi * (cElectronMass/cPlanckH) * (cEV/cPlanckH)),
  ! so that
  ! 1/Lambda * sqrt(cBoltzmann/cEV * T)/Ne**(1/3) would be the inverse of
  ! the DeBroigle wavelength         
  !\
  ! Numerical value: eWeight = 6.0372842638557644E+27 [1/(m^3eV^{3/2})]   
  !                1/eWeight = 1.6563738864953516E-28 [m^3eV^{3/2}]
  !/

  real :: eWeight1eV1m3  

  !Input parameters (note though that the temperature may be not
  !directly assigned and should be found from the internal energy 
  !or pressure in this case):
  real :: Te = 1.0        ! the electron temperature [eV] = ([cKToeV] * [Te in K])
  real :: Na =-1.0        ! The density of heavy particles in the plasma [#/m^3]


  !Averages:

  real :: zAv,   &  ! The average charge per ion - <Z> (elementary charge units)
          EAv,   &  ! The average ionization energy level of ions
          Z2PerA,&  ! The average of Z^2/A
          Z2,    &  ! The average of Z^2
          Z2_I(nMixMax), &  ! The average of Z^2 for each component      
          VirialCoeffAv, &  ! The value of <-\frac{V}{T} \frac{\partial E}{\partial V}>
          SecondVirialAv    ! The value of \frac{V^2}{T} < \frac{\partial^2 E}{\partial V^2} >

  !Covariances between \delta i), delta E and \delta (Virial coeff):
  real :: DeltaZ2Av      ! The value of <delta^2 i>=<i^2>-<i>^2
  real :: DeltaETeInv2Av ! The value of <(delta E/Te)^2>
  real :: ETeInvDeltaZAv ! The value of <delta E/Te delta i>
  real :: CovEnergyVirial! The value of -\frac{V}{T^2} < \delta \frac{\partial E}{\partial V} \delta E >
  real :: Cov2Virial     ! The value of \frac{V^2}{T^2} < \delta^2 \frac{\partial E}{\partial V} >
  real :: CovVirialZ     ! The value of \frac{V}{T} < \delta \frac{\partial E}{\partial V} \delta i >



  integer :: iIter   !To provide the output for the convergence efficiency, if needed



  !Auxiliary arrys of numerical constants:
  real,dimension(0:nZMax) :: N_I ! array of consecutive integers (with type real)
  real,dimension(nZMax)   :: LogN_I ! array of natural logarithms of consecutive integers

  public:: Te,&             ! The electron temperature [eV] = ([cKToeV] * [Te in K])
           Na,&             ! The density of heavy particles in the plasma [#/m^3]
           zAv,&            ! The average charge per ion - <Z> (elementary charge units)
           EAv,&            ! The average ionization energy level of ions
           Z2PerA,&         ! The average of Z^2/A
           Z2,    &         ! The average of Z^2
           Z2_I, &          ! The average of Z^2 for each component
           VirialCoeffAv,&  ! The value of <-\frac{V}{T} \frac{\partial E}{\partial V}>
           SecondVirialAv,& ! The value of \frac{V^2}{T} < \frac{\partial^2 E}{\partial V^2} >          
           DeltaZ2Av,&      ! The value of <(delta i)^2>=<i^2>-<i>^2
           DeltaETeInv2Av,& ! The value of <(delta E/Te)^2>
           ETeInvDeltaZAv,& ! The value of <Delta E/Te Delta i>
           CovEnergyVirial,&! The value of -\frac{V}{T^2} < \delta \frac{\partial E}{\partial V} \delta E >
           Cov2Virial,&     ! The value of \frac{V^2}{T^2} < \delta^2 \frac{\partial E}{\partial V} >
           CovVirialZ,&     ! The value of \frac{V}{T} < \delta \frac{\partial E}{\partial V} \delta i >
           iIter,&          ! To provide the output for the convergence efficiency, if needed
           Population_II,&  ! To calculate the spectroscopic data
           nMix, IonizPotential_II
  real,parameter::StrangeTemperature = -777.0

  ! Virial coefficients:
  !\
  ! The Madelung energy: the electrostatic energy of Z electrons
  ! uniformly filling in the iono-sphere diveded by Te

  real,public :: eMadelungPerTe  

  ! Correlation energy in the Debye-Huekel thoery of 
  ! correlations in weakly non-ideal plasma 

  real,public :: eDebyeHuekel  !In eV, per atomic cell

  
  !\
  ! The logical to handle the Coulomb corrections
  ! should be declared in this module, because
  ! it must be applied as the (1) correction to the
  ! energy levels or (2) the shift in the ionization
  ! equilibrium, rather than in the pressure/energy only.
  !/
  logical,public :: UseCoulombCorrection = .false.

  !The "Coulomb logariphm". May be defined quantitatively
  !only in a conjunction with the choice for the model
  !of Coulomb interactions (DH or iono-sphere or 
  !"average atom")
  real, public :: CoulombLog = 2.0

  real    :: IonizationEnergyLowering_I(0:nZMax) = 0.0

  ! Inverse of the iono-sphere radius,
  ! measured in Bohr radius

  real, public :: rIonoSphereInv

  !Public methods:

  !Set the elements and their Ionization Potentials
  public:: set_mixture      

  ! Find the final values of zAv and the ion populations 
  ! from temperature and heavy particle density.

  public:: init_madelung
  public:: set_ionization_equilibrium

  public:: set_zero_ionization
  public:: check_applicability, get_virial_coef

  !Tolerance parameters: should be reset if the tolerance in
  !CRASH_ModStatSum is reset
  real, public :: ToleranceZ = 0.001!  !Accuracy of Z needed
  real, public :: StatSumToleranceLog = 7.0
Contains
  !========================================================================!
  !Several initialazations of constants, such as                           !
  !calculating the natural logarithms of the first nZMax integers          !
  !========================================================================!
  subroutine mod_init
    use CRASH_ModFermiGas,ONLY: UseFermiGas, init_Fermi_function
    use CRASH_ModExcitationData,ONLY: init_excitation
    use CRASH_ModIonization,ONLY: init_ioniz_potential
    integer:: iZ  !Used for loops
    real   :: DeBroglieInv
    !-----------------
    LogN_I = (/(log(real(iZ)), iZ = 1,nZMax)/)
    N_I    = (/(real(iZ), iZ = 0,nZMax)/)
    IonizationEnergyLowering_I(0:nZMax) = 1.80 * cRyToEV  * N_I**2

    DeBroglieInv = sqrt(cTwoPi * (cElectronMass/cPlanckH) * (cEV/cPlanckH)) 
    !*sqrt(cBoltzmann/cEV * T) - temperature in eV

    eWeight1eV1m3 = 2*DeBroglieInv**3 ! 2/(Lambda^3)
   
    if(UseFermiGas) call init_Fermi_function
    call init_ioniz_potential
    call init_excitation
    DoInit=.false.

  end subroutine mod_init
  !=========================================================================!
  subroutine check_applicability(iErrorOut)
    integer,optional,intent(out) :: iErrorOut
    integer::iError
    !-------------------------------------
    !if(.not.present(iErrorOut))return
    iError = 0
    !Check relativistic temperature
    if (Te > 50000.0) iError = 2

    if (Te <= 0.0) then
       eMadelungPerTe = 0.0
       eDebyeHuekel = 0.0
       if (iError /= 2) iError = 1
    else
       call get_virial_coef
    end if

    if(present(iErrorOut))iErrorOut = iError
  end subroutine check_applicability
  !==========================================================================!
  subroutine get_virial_coef
    !Calculate the variables dependent on the iono-sphere radius

    ! Inverse of the Debye radius,
    ! measured in Bohr radius
    real:: rDebyeInv  
    !------------------------------

    eMadelungPerTe = 0.60 * &
         cRyToEv * rIonoSphereInv / Te

    rDebyeInv = sqrt(6.0* Z2 * cRyToEv/Te*rIonoSphereInv**3)

    eDebyeHuekel = Z2 * &
         cRyToEv * rDebyeInv

  end subroutine get_virial_coef
  !==========================================================================
  !Set the elements and their Ionization Potentials
  !==========================================================================
  subroutine set_mixture(nMixIn, nZIn_I, ConcentrationIn_I)
    integer,intent(in)::nMixIn
    integer,dimension(nMixIn),intent(in) :: nZIn_I
    real,dimension(nMixIn),intent(in) :: ConcentrationIn_I

    integer            :: iZ, iMix  ! for loops
    !--------------------------!

    if(DoInit)call mod_init
    call  set_data_for_mixture(nMixIn, nZIn_I, ConcentrationIn_I)

    IonizEnergyNeutral_II(0,1:nMix) = 0.0
    do iMix=1,nMix
       do iZ = 1,nZ_I(iMix)
          IonizEnergyNeutral_II(iZ,iMix) = &
               IonizEnergyNeutral_II(iZ-1,iMix) + IonizPotential_II(iZ,iMix)
       end do
    end do
    
    !call set_zero_ionization
  end subroutine set_mixture

  !=========================================================================!
  subroutine set_zero_ionization
    use CRASH_ModFermiGas,ONLY: RMinus, RPlus, LogGe

    integer :: iMix
    !--------------

    Population_II=0.0; Population_II(0,1:nMix) = 1.0

    !Do not overwrite LogGe if calculated but initialize it 
    !otherwise (at low Te)

    if(Te>0.0)LogGe = 99.0

    iZMin_I(1:nMix) = 0
    iZMax_I(1:nMix) = 0
    call get_excitation_levels_zero(rIonoSphereInv)


    zAv = 0.0
    EAv = 0.0
    Z2  = 0.0
    Z2PerA = 0.0
    Z2_I   = 0.0
    VirialCoeffAv    = 0.0
    SecondVirialAv   = 0.0

    DeltaZ2Av        = 0.0
    DeltaETeInv2Av   = 0.0
    ETeInvDeltaZAv   = 0.0
    CovEnergyVirial  = 0.0
    Cov2Virial       = 0.0
    CovVirialZ       = 0.0


    RMinus = 1.0; RPlus = 1.0


  end subroutine set_zero_ionization

  !========================================================================!
  ! Initialize data required for calculation of Coulomb correction
  ! according to the Madelung approximation
  subroutine init_madelung(NaIn)
    use CRASH_ModExcitation,ONLY:IonizationPotentialLowering_I

    real,intent(in) :: NaIn
    !----------------------
    if(Na/=NaIn) rIonoSphereInv = &
         ( (4.0 * cPi/3.0) * cBohrRadius**3 * NaIn)**(1.0/3)
    Na = NaIn

    if(UseCoulombCorrection)then
       IonizationPotentialLowering_I = 1.80 * rIonoSphereInv *&
         cRyToEv * (2.0 * N_I + 1.0)
    else
       IonizationPotentialLowering_I = 0.0
    end if

  end subroutine init_madelung

  !========================================================================!  
  ! Find the final values of zAv and the ion populations from temperature and 
  ! heavy particle density

  subroutine set_ionization_equilibrium(TeIn, NaIn, iError )

    use CRASH_ModFermiGas,ONLY: UseFermiGas, LogGe, &
         LogGeMinBoltzmann,LogGeMinFermi
    use CRASH_ModMolecularData, ONLY: UseDiatomicMolecules,&
         get_chemical_equilibrium

    ! Concentration of heavy particles (atoms+ions) in the plasma 
    ! (# of particles per m^3):
    real, intent(in)             ::  NaIn,& ![1/m^3]
                                     TeIn   !electron temperature [eV]
    integer,optional,intent(out) :: iError


    !Internal variables
    real :: lnC1  ! natural log C1 
    real :: TeInv

    real :: DiffZ

    !---------------------------------------------------------

    Te = TeIn
    call init_madelung(NaIn)

    if( Te <= 0.02 * minval( IonizPotential_II( 1,1:nMix) ) )then
       call set_zero_ionization
       call check_applicability(iError)
       return
    end if

    TeInv = cOne / TeIn        ! 1/kT; units: [1/eV]
    lnC1  = log( eWeight1eV1m3  * sqrt(TeIn)*TeIn / Na)

    call set_Z
    call set_averages_and_deviators(DoZOnly=.false.)
    call check_applicability(iError)
    
    !The chemical equilibrium is calculated only for weakly ionized plasma
    !However, the fractions of neutral atoms are used which are not exactly
    !equal to 0ne:

    if(UseDiatomicMolecules.and.zAv < 0.1)&
         call get_chemical_equilibrium( Te, Na, Population_II(0,1:nMix))
    !Calculate Coulomb logarithm
    CoulombLog = max(2.0, 24.0-log(sqrt(zAv*Na*1.0e-6)/Te))
  contains

    ! Calculating Z averaged iteratively
    subroutine set_Z

      use CRASH_ModFermiGas,ONLY: iterate_ge, RMinus, RPlus

      !Internal variables
      real    :: ZTrial, ZNew            ! The trial values of Z for iterations
      real,dimension(nMixMax) :: InitZ_I ! The initial approximation of Z
      real :: ZJ                         ! Average Z for each component
      integer :: iMix,iZ(1)              ! Misc
      !----------------------------------------------------------
      ! First approximate the value of Z by finding for which i=Z 
      ! the derivative of the populations sequence~0 (population is maximum):
      initial: do iMix = 1,nMix
         iZ = minloc( abs(lnC1 - LogN_I(1:nZ_I(iMix)) - &
              IonizPotential_II(1:nZ_I(iMix),iMix)*TeInv) )


         if(iZ(1)==1)then
            !Find zAv in the case when Z<=1
            InitZ_I(iMix)  = min( 1.0, exp(cHalf*( lnC1 -IonizPotential_II(1,iMix)*TeInv) ) )
         else
            !Apply the above estimate
            InitZ_I(iMix)  = real(iZ( 1)) -cHalf

         end if

      end do initial

      !Average the initial approximation across the mixture compenents:
      ZNew = sum(Concentration_I(1:nMix) * InitZ_I(1:nMix))


      ! Use Newton's method to iteratively get a better approximation of Z:
      ! Organize the iteration loop:
      iIter  =  0
      DiffZ = 2.0*ToleranceZ

      if(UseFermiGas)then
         !Iterate LogGe
         LogGe = lnC1 - log(ZNew)
         do while (abs(DiffZ) >= ToleranceZ .and. iIter < 10)

            !Set the ion population with the electron stqtistical weight which may expressed
            !in terms of ZTrial, for UseFermiGas==.false. , or iterated directly at UseFermiGas=.true.

            call set_populations(LogGe)

            !Take an average over the ion populations:
            call set_averages_and_deviators(DoZOnly = .true.)

            call iterate_ge(zAv, DeltaZ2Av, lnC1, DiffZ)

            iIter = iIter+1
         end do
      else
         !Iterate Z
         do while (abs(DiffZ) >= ToleranceZ .and. iIter < 10)

            ZTrial = ZNew
            LogGe = lnC1 - log(ZNew)


            !Set the ion population with the electron stqtistical weight which may expressed
            !in terms of ZTrial, for UseFermiGas==.false. , or iterated directly at UseFermiGas=.true.

            call set_populations(LogGe)

            !Take an average over the ion populations:
            call set_averages_and_deviators(DoZOnly = .true.)


            DiffZ = (ZTrial - zAv)/(cOne + DeltaZ2Av/ZTrial)
            ZNew = ZTrial - DiffZ


            iIter = iIter+1
         end do
      end if

    end subroutine set_Z

    !==============================================
    ! Finding the populations of the ion states
    subroutine set_populations(GeLogIn)

      real, intent(in) :: GeLogIn ! Natural logarithm of the electron stat weight:
      !  log(1/(Ne*lambda^3)) 

      real :: StatSumTermMax,StatSumTermMin

      real,dimension(0:nZMax) :: StatSumTermLog_I

      integer,dimension(1) :: iZDominant    !Most populated ion state

      ! ln(1.0e-3), to truncate terms of the statistical sum, which a factor of 
      ! 1e-3 less than the maximal one:

      integer :: iZ, iMix  !Loop variables
      real    :: PITotal   !Normalization constant, to set the populations total to equal 1.
      real    :: GeLog
      !--------------------------------------!

      if(UseFermiGas)then
         GeLog = GeLogIn
         if(LogGeMinFermi>-3.0)GeLog = max(LogGeMinFermi, GeLogIn)
      else
         ! The previous version does not stop simulation
         ! when the Fermi degeneracy occurs (the electron 
         ! statictical weight is not large), but it sets the
         ! low bound for the weight to be e^2\approx 8

         GeLog = max(LogGeMinBoltzmann, GeLogIn )
      end if

      mixture: do iMix=1,nMix

         ! First, make the sequence of ln(StatSumTerm) values; let ln(P0)=0 )
         StatSumTermLog_I(0) = 0.0

         do iZ = 1, nZ_I(iMix)       !Fill up the sequence using the following equation:

            StatSumTermLog_I(iZ) = StatSumTermLog_I(iZ-1) &
                 - IonizPotential_II(iZ,iMix)*TeInv + GeLog
         end do

         !Account for the ionization potential lowering resulting from the
         !electrostatic interaction:

         if(UseCoulombCorrection)StatSumTermLog_I(0:nZ_I(iMix)) = &
              StatSumTermLog_I(0:nZ_I(iMix)) + TeInv * rIonoSphereInv *&
              IonizationEnergyLowering_I(0:nZ_I(iMix))

         ! Find the location of that maximum value
         iZDominant = maxloc(StatSumTermLog_I(0:nZ_I(iMix)))-1 

         StatSumTermMax = StatSumTermLog_I(iZDominant(1))

         StatSumTermMin = StatSumTermMax - StatSumToleranceLog 


         ! Find the lower boundary of the array 
         ! below which the values of Pi can be neglected

         iZMin_I(iMix) = count( StatSumTermLog_I(0:iZDominant(1)) < StatSumTermMin )


         !Find the similar upper boundary
         iZMax_I(iMix) = max(nZ_I(iMix) - &
              count(StatSumTermLog_I( iZDominant(1):nZ_I(iMix) ) < StatSumTermMin),1)


         call get_excitation_levels(iMix, iZMin_I(iMix), min(iZMax_I(iMix), nZ_I(iMix)-1), &
              nZ_I(iMix), TeInv, rIonoSphereInv)


         StatSumTermLog_I(0:nZ_I(iMix)) = &
              StatSumTermLog_I(0:nZ_I(iMix)) + LogGi_II(0:nZ_I(iMix),iMix)

         !Initialize the population array to zeros
         Population_II(0:nZ_I(iMix), iMix) = 0.0


         !Convert the array into the Pi values from ln(Pi)
         Population_II(iZMin_I(iMix):iZMax_I(iMix), iMix) = &
              exp(StatSumTermLog_I(iZMin_I(iMix):iZMax_I(iMix)) - StatSumTermMax)

         !Normalize the Pi-s so that their sum =1

         !Add up all the values of Pi found so far
         PITotal = sum(Population_II(iZMin_I(iMix):iZMax_I(iMix),iMix))

         Population_II(iZMin_I(iMix):iZMax_I(iMix),iMix) = &
              Population_II(iZMin_I(iMix):iZMax_I(iMix),iMix)/PITotal

      end do mixture
    end subroutine set_populations
    !-----------------------------
    !=========================================================================!  
    !For known ion polulations find average Z, E, and bi-linear deviators
    subroutine set_averages_and_deviators(DoZOnly)
      use CRASH_ModAtomicMass,ONLY: cAtomicMass_I
      logical,intent(in)::DoZOnly
      
      ! < Ei/Te>_J (Ei - energy levels, Te - electron temperature [eV])
      real :: ETeInvAv_I(nMixMax)   ,& ! averaged i^2 for current iMix
              Z_I(nMixMax)          ,& ! averaged Z for current iMix                 
              VirialCoeff_I(nMixMax),& ! The value of  -\frac{V}{T} \frac{\partial E}{\partial V}
              DeltaZ2I              ,& ! \delta^2 iaveraged for current iMix 
              InvTR                    ! 1/(Te*R_{iono})

      ! Array of ionization energy levels of ions divided by the temperature in eV
      real,dimension(0:nZMax) :: ETeInv_I
      real,dimension(0:nZMax) :: VirialCoeffFull_I, SecondVirialCoeff_I
      integer :: iMix
      real,parameter:: cThird = 1.0/3.0, cFourNineth = 4.0/9.0
      !-----------------

      zAv = 0.0
      EAv = 0.0
      Z2  = 0.0
      Z2_I= 0.0
      Z2PerA = 0.0
      VirialCoeffAv    = 0.0
      SecondVirialAv   = 0.0

      DeltaZ2Av        = 0.0
      DeltaETeInv2Av   = 0.0
      ETeInvDeltaZAv   = 0.0
      CovEnergyVirial  = 0.0
      Cov2Virial       = 0.0
      CovVirialZ       = 0.0

      VirialCoeffFull_I=0.0
      SecondVirialCoeff_I =0.0

      InvTR            = TeInv * rIonoSphereInv

      do iMix = 1, nMix


         !Calculate average vaues of Z, for this component:
         Z_I(iMix) = sum(Population_II(iZMin_I(iMix):iZMax_I(iMix),iMix )* &
              N_I(iZMin_I(iMix):iZMax_I(iMix)))
         zAv       = zAv + Concentration_I(iMix) * Z_I(iMix)

         !Calculate <\delta i>^2>
         DeltaZ2I = sum( Population_II(iZMin_I(iMix):iZMax_I(iMix), iMix) *  &
              (N_I(iZMin_I(iMix):iZMax_I(iMix))-Z_I(iMix))**2 )

         !Calculate a << (\delta i)^2 >>
         DeltaZ2Av = DeltaZ2Av + Concentration_I(iMix) * &
              DeltaZ2I
        
         if(DoZOnly) CYCLE

         !Calculate <Z^2>
         Z2_I(iMix) = DeltaZ2I + Z_I(iMix)**2
         Z2 = Z2 +  Concentration_I(iMix) * Z2_I(iMix)
              

         !Calculate <Z^2/A>
         Z2PerA = Z2PerA +  Concentration_I(iMix) * &
              Z2_I(iMix)/cAtomicMass_I( nZ_I(iMix) )

         !Fill in the array of average energies (per Te)
         !for all iMix and i taken into account
         ETeInv_I(iZMin_I(iMix):iZMax_I(iMix)) = &
              (IonizEnergyNeutral_II( iZMin_I(iMix):iZMax_I(iMix),iMix ) + &
              ExtraEnergyAv_II( iZMin_I(iMix):iZMax_I(iMix),iMix )) * TeInv

         !Fill in the array with values of - \frac{V}{T} \frac{\partial E}{\partial V} 
         !for all iMix and i taken into account, where the mean value is taken over n and l.
         VirialCoeffFull_I(iZMin_I(iMix):iZMax_I(iMix) ) = &
              VirialCoeffAv_II(iZMin_I(iMix):iZMax_I(iMix), iMix) * TeInv

         !The same for the second virial coeff
         SecondVirialCoeff_I(iZMin_I(iMix):iZMax_I(iMix) ) = &
              TeInv * SecondVirialCoeffAv_II(iZMin_I(iMix):iZMax_I(iMix), iMix)

 

         if(UseCoulombCorrection)then
            ETeInv_I( iZMin_I(iMix):iZMax_I(iMix) ) = &
                 ETeInv_I( iZMin_I(iMix):iZMax_I(iMix) ) -InvTR  * &
                 IonizationEnergyLowering_I( iZMin_I(iMix):iZMax_I(iMix) )


            VirialCoeffFull_I(iZMin_I(iMix):iZMax_I(iMix)) = &
                 VirialCoeffFull_I(iZMin_I(iMix):iZMax_I(iMix)) - &
                 cThird * InvTR  * &
                 IonizationEnergyLowering_I( iZMin_I(iMix):iZMax_I(iMix) )
            SecondVirialCoeff_I(iZMin_I(iMix):iZMax_I(iMix)) = &
                 SecondVirialCoeff_I(iZMin_I(iMix):iZMax_I(iMix)) - &
                 cFourNineth * InvTR  * &
                 IonizationEnergyLowering_I( iZMin_I(iMix):iZMax_I(iMix) )
         end if
         
      
         !Calculate average value of E
        
         ETeInvAv_I(iMix) = sum(Population_II(iZMin_I(iMix):iZMax_I(iMix),iMix )* &
              ETeInv_I(iZMin_I(iMix):iZMax_I(iMix))) 

         EAv = EAv + Concentration_I(iMix) * ETeInvAv_I(iMix)               
  
         !Calculate <-\frac{V}{T} \frac{\partial E}{\partial V}>
         VirialCoeff_I(iMix) = sum( Population_II(iZMin_I(iMix):iZMax_I(iMix), iMix) *&
              VirialCoeffFull_I(iZMin_I(iMix):iZMax_I(iMix)))

         VirialCoeffAv = VirialCoeffAv + Concentration_I(iMix) *&
              VirialCoeff_I(iMix)


         !Calculate \frac{V^2}{T} < \frac{\partial^2 E}{\partial V^2} >
         SecondVirialAv = SecondVirialAv + Concentration_I(iMix) *&
              sum( Population_II(iZMin_I(iMix):iZMax_I(iMix), iMix) *&
              SecondVirialCoeff_I(iZMin_I(iMix):iZMax_I(iMix)) )

         !Covariances

         !Calculate <(delta E/Te)^2>
         DeltaETeInv2Av   = DeltaETeInv2Av + Concentration_I(iMix)*&
              sum( Population_II(iZMin_I(iMix):iZMax_I(iMix), iMix) * &
              ((ETeInv_I(iZMin_I(iMix):iZMax_I(iMix)) -ETeInvAv_I(iMix))**2 + &
              Cov2ExtraEnergy_II(iZMin_I(iMix):iZMax_I(iMix),iMix) * TeInv * TeInv)) !correction to account for
                                                                                     !energy levels divergence
         !Calculate <delta E/Te delta i>
         ETeInvDeltaZAv   = ETeInvDeltaZAv + Concentration_I(iMix)*&
              sum( Population_II(iZMin_I(iMix):iZMax_I(iMix), iMix) *&
              ETeInv_I(iZMin_I(iMix):iZMax_I(iMix)) * &
              (N_I(iZMin_I(iMix):iZMax_I(iMix)) - Z_I(iMix)) )


         !Calculate -\frac{V}{T^2} < \delta \frac{\partial E}{\partial V} \delta E >
         CovEnergyVirial = CovEnergyVirial + Concentration_I(iMix)*&
              sum( Population_II(iZMin_I(iMix):iZMax_I(iMix), iMix) *&
              ( TeInv * TeInv * CovExtraEnergyVirial_II(iZMin_I(iMix):iZMax_I(iMix), iMix) + &
              (VirialCoeffFull_I(iZMin_I(iMix):iZMax_I(iMix)) - VirialCoeff_I(iMix)) * &
              ETeInv_I(iZMin_I(iMix):iZMax_I(iMix)) ) )

         !Calculate \frac{V^2}{T^2} < \delta^2 \frac{\partial E}{\partial V} >
         Cov2Virial = Cov2Virial + Concentration_I(iMix) *&
              sum( Population_II(iZMin_I(iMix):iZMax_I(iMix), iMix) *&
              (TeInv*TeInv * Cov2VirialCoeff_II(iZMin_I(iMix):iZMax_I(iMix), iMix) + &
              (VirialCoeffFull_I(iZMin_I(iMix):iZMax_I(iMix)) - VirialCoeff_I(iMix))**2) )

         !Calculate \frac{V}{T} < \delta \frac{\partial E}{\partial V} \delta i >
         CovVirialZ = CovVirialZ - Concentration_I(iMix) *&
              sum( Population_II(iZMin_I(iMix):iZMax_I(iMix), iMix) *&
              (N_I(iZMin_I(iMix):iZMax_I(iMix)) - Z_I(iMix)) * &
              VirialCoeffFull_I(iZMin_I(iMix):iZMax_I(iMix)))
   

      end do
      
      if(nMix > 1 &
!           .and..false.& !Temporarily disabled
           )then
         do iMix=1,nMix
            DeltaZ2Av = DeltaZ2Av + Concentration_I(iMix)*&
                 (Z_I(iMix) - zAv) **2
            if(DoZOnly)CYCLE
            DeltaETeInv2Av   = DeltaETeInv2Av + Concentration_I(iMix)*&
                 (ETeInvAv_I(iMix) - EAv)**2

            ETeInvDeltaZAv   = ETeInvDeltaZAv + Concentration_I(iMix)*&
                 (ETeInvAv_I(iMix) - EAv) * (Z_I(iMix) - zAv)

            CovEnergyVirial = CovEnergyVirial + Concentration_I(iMix)*&
                 (ETeInvAv_I(iMix) - EAv) * (VirialCoeff_I(iMix) - VirialCoeffAv)
            
            Cov2Virial = Cov2Virial + Concentration_I(iMix)*&
                 (VirialCoeff_I(iMix) - VirialCoeffAv)**2
            
            CovVirialZ = CovVirialZ - Concentration_I(iMix)*&
                 (Z_I(iMix) - zAv) * (VirialCoeff_I(iMix) - VirialCoeffAv)
         end do
      end if

      EAv = EAv * Te

    end subroutine set_averages_and_deviators
  end subroutine set_ionization_equilibrium
  !========================================!

end module CRASH_ModPartition
