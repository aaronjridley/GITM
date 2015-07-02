!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CRASH_ModEos

  ! Equation Of State (EOS) for ionized plasma
  !
  ! Thermodynamic variables and other notations
  !
  !    Rho - the mass density
  !    E - internal energy of the unit of mass
  !    e, i - electron, ion
  !    V, vol - volume or volumetric
  !    \left(\frac{\partial...}{\partial...}\right)_V - thermodynamical
  !         derivative at constant volume
  !    Te, Ti - electron and ion temperature
  !    iMaterial - integer index of the material:
  !    iMaterial=0 - xenon
  !    iMaterial=1 - beryllium   
  !    iMaterial=2 - plastic (or polyimide, if more than one plastic is used).
  !    iMaterial=3 - gold
  !    iMaterial=4 - Acrylic (acronym is Ay_, since Ac and Ar are both in use.
  !    iMaterial=90 - plasma with eos E=aT^4/4; p=aT^4/12; C_V=aT^3 
  !
  ! In the initial CRASH treatment of materials,
  !
  ! 1. We can treat any given material as having a single average Z.
  !
  ! 2. Mixtures can be treated either as an average material or as mixtures.
  !
  ! 3. If mixtures are treated as mixtures, then collisional rates should be 
  !    calculated using the "effective Z", which is the average of Z squared 
  !    divided by the average Z.
  !
  ! 4. The average Z can be determined from the Saha equation, but must not 
  !    exceed the nuclear charge of the material in question.
  !
  ! 5. We do not need to account for electron degeneracy in the initial model.
  !
  ! 6. In our regime of interest, the electrons behave as an ideal gas in an 
  !    ion-sphere environment within which Coulomb interactions do affect the 
  !    electron pressure and internal energy. The electron pressure and 
  !    internal energy are best calculated using equations 3.47 through 3.50 
  !    in R. P. Drake, High Energy Density Physics
  !
  ! 7. The ion pressure is the ideal gas pressure. The ion internal energy 
  !    includes the particle energy of random motion and the energy of 
  !    ionization. The model of eqs. 3.74 through 3.76 in the mentioned book 
  !    is acceptable. Alternatively, a more complex model using actual 
  !    ionization energies would be acceptable.
  !
  !    !!! WARNING!!! Correction in this item !!!
  !
  !     Since the ionization partition function is controlled by the 
  !     electron temperature, the ionization energy as well as not mentioned
  !     excitation energy are both included to the ELECTRON ENERGY DENSITY. 
  !     See detail in util/CRASH/doc/HEDP.pdf. 
  !     To make this document go to util/CRASH/doc/Tex directory and make PDF
  !
  ! 8. The materials that matter are
  !    - Beryllium
  !    - Xenon
  !    - Polyimide (C_22 H_10 N_2 O_5)
  !
  ! Error flag (iError) values:
  ! 0 - OK 
  ! 2 - relativistic temperature
  !\
  !!   WARNING !!!
  !You cannot use total pressure and total energy density as input or output
  !parameters, if the electron temperature is not equal to ion temperature.
  !In this case ONLY electron energy density and electron pressure may be 
  !used.
  !/
  use CRASH_ModPolyimide,ONLY:cAPolyimide
  use CRASH_ModAcrylic, ONLY:cAAcrylic
  use CRASH_ModAtomicDataMix, ONLY: nMixMax
  use CRASH_ModStatSum
  use CRASH_ModAtomicMass
  use CRASH_ModPowerLawEos
  use CRASH_ModFermiGas, ONLY: UseFermiGas, LogGeMinBoltzmann, LogGeMinFermi
  use CRASH_ModMultiGroup, ONLY: meshhv, abscon, nGroup, &
       OpacityPlanck_I, OpacityRosseland_I, opacys

  implicit none

  private !Except

  public:: read_eos_parameters   ! Read parameters from input parameter file

  integer, public, parameter:: Xe_=0      ! Xenon
  integer, public, parameter:: Be_=1      ! Beryllium
  integer, public, parameter:: Plastic_=2 ! Polyimide (C_22 H_10 N_2 O_5)
  integer, public, parameter:: Au_=3      ! Gold
  integer, public, parameter:: Ay_=4      ! Acrylic

  public:: cAtomicMass_I, cAPolyimide ! inherited from ModPolyimide
  public:: cAAcrylic

  public:: UsePreviousTe ! inherited from CRASH_ModStatSum

  !This the main eos function, which may be implemented both via 
  !internal or external  eos tables and via the built-in EOS model 
  public:: eos

  interface eos
     module procedure eos_material
     module procedure eos_mixture
  end interface

  ! Local variables

  ! test material with the EOS e \propto T^4, p \propto T^4
  integer, parameter:: Test_ = 90 

  ! Maximum and actual number of materials
  integer, public, parameter:: nMaterialMax = 5
  integer, public           :: nMaterialEos = -1

  ! array of atomic numbers
  integer:: nZ_I(0:nMaterialMax-1)=&
       (/54,&  !Xenon 
       4,   &  !Beryllium
       -4,  &  !Minus number of elements in polyimide
       79,  &  !Gold
       -3/)    !Minus number of elements in acrylic

  real, dimension(0:nMaterialMax-1), public :: cAtomicMassCRASH_I=&
       (/cAtomicMass_I(54),                                            &!  Xe
       cAtomicMass_I( 4),                                              &!  Be
       (cAtomicMass_I(6)*22 + cAtomicMass_I(1)*10 + cAtomicMass_I(7)*2 &
       + cAtomicMass_I(8)*5) / (22 + 10 + 2 + 5),                      &!  Pl
       cAtomicMass_I(79),                                              &!  Au
       (cAtomicMass_I(6)*5 + cAtomicMass_I(8)*2 + cAtomicMass_I(1)*8)/15/)!Ay
 
  character(LEN=2), public ::&
       NameMaterial_I(Xe_:Ay_) = (/'Xe','Be','Pl','Au','Ay'/)

  !Arrays for mixtures
  integer, dimension(1:nMixMax, 0:nMaterialMax-1):: nZMix_II=reshape((/&
                                   54, 0, 0, 0, 0, 0, &
                                   4, 0, 0, 0, 0, 0, &
                                   6, 1, 7, 8, 0, 0, &
                                   79, 0, 0, 0, 0, 0, &
                                   6, 8, 1, 0, 0, 0/),(/6,5/))

  real, dimension(1:nMixMax, 0:nMaterialMax-1) :: cMix_II=reshape((/&
                                   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                   22.0/39, 10.0/39, 2.0/39, 5.0/39,0.0,0.0,&
                                   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                   5.0/15, 2.0/15, 8.0/15,0.0,0.0,0.0/),(/6,5/))

  public :: nZMix_II, cMix_II   !To calculate the average atomic number


  ! The logicals determine if we use or not tabulated EOS
  ! and opacities  
  
  logical, public, dimension(0:nMaterialMax-1) :: &
       UseEosTable_I = .false., UseOpacityTable_I = .false.
 
  !The columns in the EOS table
  integer,public :: P_      =  1, &
             E_      =  2, &
             Pe_     =  3, &
             Ee_     =  4, &
             Cv_     =  5, &
             Cve_    =  6, &
             Gamma_  =  7, &
             GammaE_ =  8, &
             TeTi_   =  9, &
             Cond_   = 10, &
             Z_      = 11, &
             Z2_     = 12, &
             DPOverDRho_ = 13, &
             DPOverDT_   = 14, &
             DPEOverDRho_= 15, &
             DPEOverDT_  = 16

  !The number of columns in the EOS table
  integer,public :: nVarEos =16   
  
  ! Arrays which relate the iTable for the EOS table with 
  ! the material number:
  ! integer:: iMaterial4EosTable_I(MaxTable) = -1
  integer,public:: iTableEos4Material_I(0:nMaterialMax-1) = -1

  ! Arrays which relate the iTable for the opacity table with 
  ! the material number:
  ! integer:: iMaterial4OpacTable_I(MaxTable) = -1
  integer,public:: iTableOpac4Material_I(0:nMaterialMax-1) = -1 

  !\
  ! Miscellaneous subroutines (probably, redundant)
  !/
  public:: fix_hyades_state

contains
  !============================================================================
  subroutine read_eos_parameters(nMaterial, NameCommand)

    use ModReadParam, ONLY: read_var

    integer,          intent(in):: nMaterial
    character(len=*), intent(in):: NameCommand

    logical:: UseTable

    character(len=*), parameter:: NameSub = 'read_eos_parameters'
    !--------------------------------------------------------------------------
    nMaterialEos = nMaterial

    select case(NameCommand)
    case("#EOS")
       call read_var('UseFermiGas',       UseFermiGas      )
       call read_var('LogGeMinBoltzmann', LogGeMinBoltzmann)
       call read_var('LogGeMinFermi',     LogGeMinFermi    )
    case("#EOSTABLE", "#USEEOSTABLE")
       call read_var('UseEosTable', UseTable)
       UseEosTable_I(0:nMaterial-1) = UseTable
    case("#OPACITYTABLE", "#USEOPACTABLE")
       call read_var('UseOpacityTable', UseTable)
       UseOpacityTable_I(0:nMaterial-1) = UseTable
    case("#MATERIAL")
       call read_name_material
    case default
       call CON_stop(NameSub//' unknown command='//NameCommand)
    end select

  end subroutine read_eos_parameters
  !============================================================================
  subroutine read_name_material

    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: split_string
    use CRASH_ModAtomicNotation, ONLY: MaterialMin_, i_material,&
         nZMixStored_II, cMixStored_II

    ! #MATERIAL
    ! Ar            ! 2-symbol name of chemical element
    ! H_            ! 2-symbol name of chemical element 
    ! N             ! Will be extended to N_
    ! Ay            ! The chemical formula for acrylic is stored
    ! Wa:H 2 O 1    ! The chemical formula will be used to make and save
    !    or         ! the table Wa_eos as Wa_eos_CRASH.dat. Next time
    ! Wa            ! the use of Wa in the PARAM.in file is allowed, the 
    !               ! atomic mass will be read from the table

    integer:: iMaterial, iMix
    character(LEN=30):: Name

    integer, parameter:: MaxString = 200
    character(LEN=20):: NameVar_I(MaxString)
    !------------------------------------------------------------------------
    cAtomicMassCRASH_I = -1.0
    nZ_I = 0
    cMix_II  = 0.0
    nZMix_II = 0
    do iMaterial = 0, nMaterialEos-1
       call read_var('Name material',Name)
       if(index(Name,':')>0)then
          call CON_stop('The use of new mixtures with CRASH-created '//&
               'tables is to be tested')
          NameMaterial_I(iMaterial) = Name(1:2)
          !Truncate symbols till the colon
          call split_string( Name(1:len_trim(Name)), &
               MaxString,  NameVar_I, nZ_I(iMaterial))
          
          !Re-assign nZ_I(iMaterial) to -nMix
          nZ_I(iMaterial) = - nZ_I(iMaterial)/2
          
          
          
          !Calculate atomic mass
          cAtomicMassCRASH_I(iMaterial) = 0.0
          
          do iMix = 1,- nZ_I(iMaterial)
             cAtomicMassCRASH_I(iMaterial) = &
                  cAtomicMassCRASH_I(iMaterial) +&
                  cAtomicMass_I( nZMix_II(iMix,iMaterial) ) * &
                  cMix_II(iMix,iMaterial)
          end do
       else
          !Extend the name if needed
 
          if(len_trim(Name)==1)then
             NameMaterial_I(iMaterial)=trim(Name)//'_'
          else
             NameMaterial_I(iMaterial)=trim(Name)
          end if
          !Check if this is a chemical element or a known mixture

          nZ_I(iMaterial) = i_material(NameMaterial_I(iMaterial))
          if(nZ_I(iMaterial) < MaterialMin_)then
             !No internal information
             nZ_I(iMaterial)=0

          elseif(nZ_I(iMaterial) <=0)then
             !Now nZ_I(iMaterial) is the reference number in the
             !list of the materials with stored data
             !Apply storde data:

             nZMix_II(:,iMaterial) = nZMixStored_II(:,nZ_I(iMaterial))
             cMix_II( :,iMaterial) = cMixStored_II( :,nZ_I(iMaterial))
             
             !Re-assign nZ_I(iMaterial) to -nMix
             nZ_I(iMaterial) = -count( nZMix_II(:,iMaterial)/=0)

             !Calculate atomic mass
             cAtomicMassCRASH_I(iMaterial) = 0.0

             do iMix = 1,- nZ_I(iMaterial)
                cAtomicMassCRASH_I(iMaterial) = &
                     cAtomicMassCRASH_I(iMaterial) +&
                     cAtomicMass_I( nZMix_II(iMix,iMaterial) ) * &
                     cMix_II(iMix,iMaterial)
             end do
          else
             !This is a chemical element, nZ_I is properly set
             !Find the atomic mass:
             cAtomicMassCRASH_I(iMaterial) = &
                  cAtomicMass_I(nZ_I(iMaterial))
             nZMix_II(1,iMaterial) = nZ_I(iMaterial)
             cMix_II( 1,iMaterial) = 1.0
          end if
       end if
    end do
  end subroutine read_name_material

  !============================================================================
  ! Beginning of EOS functions
  !
  !For two different kinds of input parameters for
  !to characterize the material (either the material number,
  !or the vector of their contents in a mixture) we apply
  !eos_material or eos_mixture. Both subroutine call eos_generic
  !and have very similar set of optional input and output parameters.
  !============================================================================

  subroutine eos_material(iMaterial,Rho,&
       TeIn, eTotalIn, pTotalIn, eElectronIn, pElectronIn,   &
       TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
       eElectronOut, pElectronOut, GammaEOut, CvElectronOut, &
       OpacityPlanckOut_I, OpacityRosselandOut_I,            &
       HeatCond, TeTiRelax, Ne, zAverageOut, z2AverageOut,   &
       DPOverDRho, DPOverDT, DPEOverDRho, DPEOverDT, iError)

    use ModLookupTable, ONLY: interpolate_lookup_table
    ! Eos function for single material

    integer, intent(in):: iMaterial     ! index of material
    real,    intent(in):: Rho           ! mass density [kg/m^3]
    !\
    !!   WARNING !!!
    !You cannot use total pressure and total energy density as input or output
    !parameters, if the electron temperature is not equal to ion temperature.
    !In this case ONLY electron energy density and electron pressure may be 
    !used.
    !/

    ! One of the following five energetic input parameters must be present
    real,    optional, intent(in)  :: TeIn         ! temperature SI[K]
    real,    optional, intent(in)  :: eTotalIn     ! internal energy density
    real,    optional, intent(in)  :: pTotalIn     ! pressure
    real,    optional, intent(in)  :: eElectronIn  ! internal energu density of electrons
    real,    optional, intent(in)  :: pElectronIn  ! pressure of electrons
    
    ! One or more of the output parameters can be present
    real,    optional, intent(out) :: TeOut        ! temperature
    real,    optional, intent(out) :: pTotalOut    ! pressure
    real,    optional, intent(out) :: eTotalOut    ! internal energy density
    real,    optional, intent(out) :: GammaOut     ! polytropic index
    real,    optional, intent(out) :: CvTotalOut   ! specific heat / unit volume
                                                   ! Electrons !!!!!!   
    real,    optional, intent(out) :: pElectronOut ! pressure
    real,    optional, intent(out) :: eElectronOut ! internal energy density
    real,    optional, intent(out) :: GammaEOut    ! polytropic index
    real,    optional, intent(out) :: CvElectronOut! specific heat / unit volume
    real,    optional, intent(out) :: Ne           ! electron concentration [m-3]
    real,    optional, intent(out) :: zAverageOut  ! <z>
    real,    optional, intent(out) :: z2AverageOut ! <z^2>

    real,    optional, intent(out), &              ! Opacities
                   dimension(nGroup) :: OpacityPlanckOut_I, OpacityRosselandOut_I

    real,    optional, intent(out) :: HeatCond     ! electron heat conductivity (SI)
    real,    optional, intent(out) :: TeTiRelax    ! electron-ion interaction rate (SI)
    real,    optional, intent(out) :: DPOverDRho   ! (\rho/P)(\partial P/\partial \rho)_T 
    real,    optional, intent(out) :: DPOverDT     ! (1/N_a k_B)(\partial P/\partial T)_\rho
    real,    optional, intent(out) :: DPEOverDRho  ! (\rho/P_e)(\partial P_e/\partial \rho)_T 
    real,    optional, intent(out) :: DPEOverDT    ! (1/N_a k_B)(\partial P_e/\partial T)_\rho
    integer, optional, intent(out) :: iError       ! error flag

    real   :: Natomic
    real   :: TeEV, Value_V(1:nVarEos), Opacity_V(1:2*nGroup)
    integer:: iTable
    character(LEN=*), parameter:: NameSub = 'eos_material'
    !-------------------------------------------------------------------------
    if(iMaterial == Test_)then
       call eos_esimt4(TeIn, eTotalIn, pTotalIn, &
            TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut)
       if(present(iError))iError = 0
       RETURN
    end if

    !Get the atomic concentration
    Natomic = Rho /  ( cAtomicMass * cAtomicMassCRASH_I(iMaterial) )

    if(UseEosTable_I(iMaterial))then
       iTable = iTableEos4Material_I(iMaterial)
       
       if(present(TeIn))then

          TeEV = TeIn * cKToEV
          call interpolate_lookup_table(iTable, TeEV, Natomic, Value_V, DoExtrapolate=.false.)

       elseif(present(eTotalIn))then

          ! Get an energy per the atomic cell, express in eV
          ! Find temperature from dentity and internal energy
          call interpolate_lookup_table(iTable, E_, eTotalIn/ (cEV * Natomic), Natomic, &
               Value_V, Arg1Out = TeEV, DoExtrapolate=.false.)
       

       elseif(present(pTotalIn))then
          ! Divide pressure by Na , express in eV
          !Find temperature from dentity and pressure
          call interpolate_lookup_table(iTable, P_,  pTotalIn / (cEV * Natomic), Natomic, &
               Value_V, Arg1Out = TeEV, DoExtrapolate=.false.)
          
     
       elseif(present(eElectronIn))then
          ! Get an energy per the atomic cell, express in eV
          ! Find temperature from dentity and internal energy
          call interpolate_lookup_table(iTable, Ee_, eElectronIn/ (cEV * Natomic), Natomic, &
               Value_V, Arg1Out = TeEV, DoExtrapolate=.false.)


       elseif(present(pElectronIn))then

          ! Divide pressure by Na , express in eV
          !Find temperature from dentity and pressure
          call interpolate_lookup_table(iTable, Pe_,  pElectronIn / (cEV * Natomic), Natomic, &
               Value_V, Arg1Out = TeEV, DoExtrapolate=.false.)
          
      
       else
          call CON_stop(NameSub// &
               ': none of Te, eTotal, or pTotal is among the input parameters')
       end if

       if(present(TeOut))      TeOut     = TeEV*cEVToK
       if(present(eTotalOut))  eTotalOut = Natomic*cEV*Value_V(E_)
       if(present(pTotalOut))  pTotalOut = Natomic*cEV*Value_V(P_)
       if(present(eElectronOut)) eElectronOut = Natomic*cEV*Value_V(Ee_)
       if(present(pElectronOut)) pElectronOut = Natomic*cEV*Value_V(Pe_)
       if(present(GammaEOut))  GammaEOut = Value_V(GammaE_)
       if(present(GammaOut))   GammaOut  = Value_V(Gamma_)
       if(present(CvTotalOut)) CvTotalOut = (Natomic*cBoltzmann)*Value_V(Cv_)
       if(present(CvElectronOut)) CvElectronOut = (Natomic*cBoltzmann)*Value_V(Cve_)
   
       if(present(HeatCond))   HeatCond  = Value_V(Cond_)
       if(present(TeTiRelax))  TeTiRelax = Value_V(TeTi_)
       if(present(Ne))         Ne = Value_V(Z_) * NAtomic
       if(present(zAverageOut))zAverageOut = Value_V(Z_)
       if(present(z2AverageOut))z2AverageOut = Value_V(Z2_)
       if(present(DPOverDRho))DPOverDRho = Value_V(DPOverDRho_)
       if(present(DPOverDT))DPOverDT = Value_V(DPOverDT_)
       if(present(DPEOverDRho))DPEOverDRho = Value_V(DPEOverDRho_)
       if(present(DPEOverDT))DPEOverDT = Value_V(DPEOverDT_)

       if(present(iError))then
          iError=0
          if(Z_>=1)then
             if(Value_V(Z_)<=cZMin)iError=4
          end if
       end if
       if(present(OpacityPlanckOut_I).or.&
            present(OpacityRosselandOut_I))then
          if(UseOpacityTable_I(iMaterial))then
             iTable = iTableOpac4Material_I(iMaterial)
             call interpolate_lookup_table(iTable, Rho, TeEV, Opacity_V, DoExtrapolate=.false.)
             if(present(OpacityPlanckOut_I))OpacityPlanckOut_I=Opacity_V(1:nGroup) * Rho
             if(present(OpacityRosselandOut_I))&
                  OpacityRosselandOut_I=Opacity_V(nGroup+1:2*nGroup) * Rho
             return
          end if
          !Else: we need to calculate opacities only
       else
          !Opacities are not needed, all the other parameters are calculated
          return
       end if
    end if

    if(nZ_I(iMaterial)<0)then
       call set_mixture(-nZ_I(iMaterial),&
            nZMix_II(1:-nZ_I(iMaterial),iMaterial),&
            cMix_II(1:-nZ_I(iMaterial), iMaterial) )
    else
       call set_element(nZ_I(iMaterial))
       
    end if
    if(UseEosTable_I(iMaterial))then
       !We need only opacities
       call eos_generic(Natomic, &
                  TeIn = TeEV * cEvToK, &
                  OpacityPlanckOut_I = OpacityPlanckOut_I, &
                  OpacityRosselandOut_I = OpacityRosselandOut_I , &
                  iError= iError)
    else
       call eos_generic(Natomic, &
            TeIn, eTotalIn, pTotalIn, eElectronIn, pElectronIn,   &
            TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
            eElectronOut, pElectronOut, GammaEOut, CvElectronOut, & 
            OpacityPlanckOut_I, OpacityRosselandOut_I,            &
            HeatCond, TeTiRelax, Ne, zAverageOut, z2AverageOut,   &
            DPOverDRho, DPOverDT, DPEOverDRho, DPEOverDT, iError)
    end if
  end subroutine eos_material

  !============================================================================
  !Cannot be used for mixed-cell simulations if gold and/or Acrylic is used
  subroutine eos_mixture(RhoToARatio_I,&
       TeIn, eTotalIn, pTotalIn, eElectronIn, pElectronIn,   &
       TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
       eElectronOut, pElectronOut, GammaEOut, CvElectronOut, &
       OpacityPlanckOut_I, OpacityRosselandOut_I,            & 
       HeatCond, TeTiRelax, Ne, zAverageOut, z2AverageOut,   &
       DPOverDRho, DPOverDT, DPEOverDRho, DPEOverDT, iError)
    !\
    !!   WARNING !!!
    !You cannot use total pressure and total energy density as input or output
    !parameters, if the electron temperature is not equal to ion temperature.
    !In this case ONLY electron energy density and electron pressure may be 
    !used.
    !/
    use CRASH_ModPolyimide
    ! Eos function for mixed material
    real, intent(in) :: RhoToARatio_I(Xe_:Plastic_) ! Mass densities/A

    ! One of the following five energetic input parameters must be present
    real,    optional, intent(in)  :: TeIn         ! temperature
    real,    optional, intent(in)  :: eTotalIn     ! internal energy density
    real,    optional, intent(in)  :: pTotalIn     ! pressure
    real,    optional, intent(in)  :: eElectronIn  ! internal energu density of electrons
    real,    optional, intent(in)  :: pElectronIn  ! pressure of electrons

    ! One or more of the output parameters can be present
    real,    optional, intent(out) :: TeOut        ! temperature
    real,    optional, intent(out) :: pTotalOut    ! pressure
    real,    optional, intent(out) :: eTotalOut    ! internal energy density 
    real,    optional, intent(out) :: GammaOut     ! polytropic index
    real,    optional, intent(out) :: CvTotalOut   ! specific heat per volume
    ! Electrons !!!!!!   
    real,    optional, intent(out) :: pElectronOut ! pressure
    real,    optional, intent(out) :: eElectronOut ! internal energy density
    real,    optional, intent(out) :: GammaEOut    ! polytropic index
    real,    optional, intent(out) :: CvElectronOut! specific heat / unit volume
    real,    optional, intent(out) :: Ne           ! electron concentration, [m-3]
    real,    optional, intent(out) :: zAverageOut  ! <z>
    real,    optional, intent(out) :: z2AverageOut ! <z^2>
    real,    optional, intent(out) :: DPOverDRho   ! (\rho/P)(\partial P/\partial \rho)_T 
    real,    optional, intent(out) :: DPOverDT     ! (1/N_a k_B)(\partial P/\partial T)_\rho
    real,    optional, intent(out) :: DPEOverDRho  ! (\rho/P_e)(\partial P_e/\partial \rho)_T 
    real,    optional, intent(out) :: DPEOverDT    ! (1/N_a k_B)(\partial P_e/\partial T)_\rho

    real,    optional, intent(out), &              !Opacities m^-1
                   dimension(nGroup) :: OpacityPlanckOut_I, OpacityRosselandOut_I


    real,    optional, intent(out) :: HeatCond     ! electron heat conductivity (SI)
    real,    optional, intent(out) :: TeTiRelax    ! electron-ion interaction rate (SI)
    integer, optional, intent(out) :: iError       ! error flag

    real :: RhoToATotal, Natomic

    integer, parameter :: nAll = 1 + 1 + nPolyimide   

    integer, parameter :: nZAll_I(nAll) = (/54 , &  !Xe
         4 , &  !Be
         6 , &  !C
         1 , &  !H
         7 , &  !N
         8 /)   !O
    real :: ConcentrationAll_I(nAll)
    !-------------------------------------------------------------------------
    RhoToATotal = sum( RhoToARatio_I ) 

    !Relative atomic concentrations of Xe, Be and polyimide:
    ConcentrationAll_I( 1:3 ) = RhoToARatio_I / RhoToATotal

    !Specify concentrations for C, H, N, O

    ConcentrationAll_I( 3:6 ) = ConcentrationAll_I( 3 ) * CPolyimide_I

    call set_mixture(nAll, nZAll_I, ConcentrationAll_I)

    !Get the atomic concentration
    Natomic = RhoToATotal / cAtomicMass 

    call eos_generic(Natomic, &
         TeIn, eTotalIn, pTotalIn, eElectronIn, pElectronIn,   &
         TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
         eElectronOut, pElectronOut, GammaEOut, CvElectronOut, &
         OpacityPlanckOut_I, OpacityRosselandOut_I,            & 
         HeatCond, TeTiRelax, Ne, zAverageOut, z2AverageOut,   &
         DPOverDRho, DPOverDT, DPEOverDRho, DPEOverDT,  iError)

  end subroutine eos_mixture

  !============================================================================

  subroutine eos_generic(Natomic, &
       TeIn, eTotalIn, pTotalIn, eElectronIn, pElectronIn,   &
       TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
       eElectronOut, pElectronOut, GammaEOut, CvElectronOut, & 
       OpacityPlanckOut_I, OpacityRosselandOut_I,            &
       HeatCond, TeTiRelax, Ne, zAverageOut, z2AverageOut,   &
       DPOverDRho, DPOverDT, DPEOverDRho, DPEOverDT, iError)
    use CRASH_ModTransport, ONLY: electron_heat_conductivity, te_ti_relaxation
    use CRASH_ModPartition, ONLY: zAv, Z2
    !\
    !!   WARNING !!!
    !You cannot use total pressure and total energy density as input or output
    !parameters, if the electron temperature is not equal to ion temperature.
    !In this case ONLY electron energy density and electron pressure may be 
    !used.
    !/

    real,              intent(in)  :: Natomic      ! Atomic concentration
    real,    optional, intent(in)  :: TeIn         ! temperature
    real,    optional, intent(in)  :: eTotalIn     ! internal energy density
    real,    optional, intent(in)  :: PtotalIn     ! pressure
    real,    optional, intent(in)  :: eElectronIn  ! internal energu density of electrons
    real,    optional, intent(in)  :: pElectronIn  ! pressure of electrons


    real,    optional, intent(out) :: TeOut        ! temperature
    real,    optional, intent(out) :: pTotalOut    ! pressure
    real,    optional, intent(out) :: eTotalOut    ! internal energy density 
    real,    optional, intent(out) :: GammaOut     ! polytropic index
    real,    optional, intent(out) :: CvTotalOut   ! Specific heat per volume
    ! Electrons !!!!!!   
    real,    optional, intent(out) :: pElectronOut ! pressure
    real,    optional, intent(out) :: eElectronOut ! internal energy density
    real,    optional, intent(out) :: GammaEOut    ! polytropic index
    real,    optional, intent(out) :: CvElectronOut! specific heat / unit volume
    real,    optional, intent(out) :: Ne           ! electron concentration, [m-3]
    real,    optional, intent(out) :: zAverageOut  ! <z>
    real,    optional, intent(out) :: z2AverageOut ! <z^2>

    real,    optional, intent(out), &              ! Opacities m^-1
                   dimension(nGroup) :: OpacityPlanckOut_I, OpacityRosselandOut_I


    real,    optional, intent(out) :: HeatCond     ! electron heat conductivity (SI)
    real,    optional, intent(out) :: TeTiRelax    ! electron-ion interaction rate (SI)
    real,    optional, intent(out) :: DPOverDRho   ! (\rho/P)(\partial P/\partial \rho)_T 
    real,    optional, intent(out) :: DPOverDT     ! (1/N_a k_B)(\partial P/\partial T)_\rho
    real,    optional, intent(out) :: DPEOverDRho  ! (\rho/P_e)(\partial P_e/\partial \rho)_T 
    real,    optional, intent(out) :: DPEOverDT    ! (1/N_a k_B)(\partial P_e/\partial T)_\rho
    integer, optional, intent(out) :: iError       ! error flag

    real :: ePerAtom, pPerAtom,TeInEV  !All in eV

    character (len=*), parameter:: NameSub='CRASH_ModEos::eos'
    !----------------------------------------------------------------------!
    if(present(TeIn))then

       TeInEV = TeIn * cKToEV
       call set_ionization_equilibrium(TeInEV, Natomic, iError)

    elseif(present(eTotalIn))then

       ! Get an energy per the atomic cell, express in eV
       ePerAtom = eTotalIn/ (cEV * Natomic)

       ! Find temperature from dentity and internal energy

       call set_temperature(ePerAtom, Natomic, iError)

    elseif(present(pTotalIn))then
       ! Divide pressure by Na , express in eV
       pPerAtom = pTotalIn / (cEV * Natomic)

       !Find temperature from dentity and pressure
       call pressure_to_temperature(pPerAtom, Natomic, iError)
    elseif(present(eElectronIn))then

       ! Get an energy per the atomic cell, express in eV
       ePerAtom = eElectronIn/ (cEV * Natomic)

       ! Find temperature from dentity and internal energy

       call u_e_to_temperature(ePerAtom, Natomic, iError)

    elseif(present(pElectronIn))then
       ! Divide pressure by Na , express in eV
       pPerAtom = pElectronIn / (cEV * Natomic)

       !Find temperature from dentity and pressure
       call pressure_e_to_temperature(pPerAtom, Natomic, iError)
    else
       call CON_stop(NameSub// &
            ': none of Te, eTotal, or pTotal is among the input parameters')
    end if

    if(present(TeOut))      TeOut     = Te*cEVToK
    if(present(eTotalOut))  eTotalOut = Natomic*cEV*internal_energy()
    if(present(pTotalOut))  pTotalOut = pressure()
    if(present(GammaOut))   call get_gamma(GammaSOut=GammaOut)
    if(present(eElectronOut)) eElectronOut = Natomic*cEV*internal_energy_e()
    if(present(pElectronOut))  pElectronOut = pressure_e()
    if(present(GammaEOut))   call get_gamma(GammaSeOut=GammaEOut)
    if(present(OpacityPlanckOut_I).or.present(OpacityRosselandOut_I))then
       call meshhv
       call abscon
       call opacys(TRadIn = Te)
       if(present(OpacityPlanckOut_I)) &
            OpacityPlanckOut_I = OpacityPlanck_I(1:nGroup)*100.0       ![m^-1]
       if(present(OpacityRosselandOut_I)) &
            OpacityRosselandOut_I = OpacityRosseland_I(1:nGroup)*100.0 ![m^-1]
    end if

    if(present(HeatCond))   HeatCond = electron_heat_conductivity()
    if(present(TeTiRelax))  TeTiRelax = te_ti_relaxation()
    if(present(CvTotalOut)) CvTotalOut = (Natomic*cBoltzmann)*heat_capacity()
    if(present(CvElectronOut)) CvElectronOut = (Natomic*cBoltzmann)*heat_capacity_e()
    if(present(Ne)) Ne = NAtomic * zAv
    if(present(zAverageOut)) zAverageOut = zAv
    if(present(z2AverageOut))z2AverageOut = Z2
    if(present(DPOverDRho))DPOverDRho=compressibility_at_const_te()
    if(present(DPEOverDRho))DPEOverDRho=compressibility_at_const_te_e()
    if(present(DPOverDT))DPOverDT=d_pressure_over_d_te()
    if(present(DPEOverDT))DPEOverDT=d_pressure_e_over_d_te()

    if(present(iError).and.zAv <= cZMin)iError=4
  end subroutine eos_generic
  !===========================================================================
  !End of EOS functions
  !===========================================================================
  !============================================================================
  subroutine fix_hyades_state(iMaterial, StateCgs_V, PMinSi)
    use ModConst
    integer,intent(in)   :: iMaterial
    real   ,intent(inout):: StateCgs_V(4) !Rho[Cgs], P[Cgs], Te[KeV], Ti[Kev]

    real, optional, intent(in) :: PMinSi

    real:: DensitySi, NAtomicSi, PressureSi, TeSi, TiSi
    !--------------------------------------------------------------------------
    DensitySi  = 1.0e3 * StateCgs_V(1)
    NAtomicSi  = DensitySi/(cAtomicMassCRASH_I(iMaterial)*cAtomicMass)
    TeSi       = 1.0e3 * StateCgs_V(3) * ceVToK
    TiSi       = 1.0e3 * StateCgs_V(4) * ceVToK
    call eos(iMaterial, DensitySi, TeIn = TeSi, pTotalOut = PressureSi)
    PressureSi = PressureSi + (TiSi - TeSi) * nAtomicSi * cBoltzmann
    StateCgs_V(2) = 10.0 * PressureSi
    if(present(PMinSi))then
       if(PressureSi < PMinSi)then
          PressureSi = PMinSi
          StateCgs_V(2) = 10.0 * PressureSi
          call eos(iMaterial, DensitySi, pTotalIn = PressureSi, TeOut = TeSi)
          StateCgs_V(3:4) = TeSi * cKToeV * 1.0e-3
       end if
    end if
  end subroutine fix_hyades_state

end module CRASH_ModEos
