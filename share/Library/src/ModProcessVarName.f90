!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! ^CFG COPYRIGHT UM
! ======================================================
module ModProcessVarName

  implicit none

  private
  public:: process_var_name
  public:: nVarMax

  interface process_var_name
     module procedure process_var_list, process_var_string
  end interface
  

  integer,parameter :: nVarMax = 100   ! maximum number of state variables
  integer,parameter :: nSubstance = 32 ! number of distinct fluids/species

  ! Number of state variables associated with each substance to be standarized
  integer,parameter :: nVarPerSubstance = 7

  ! Number of allowed alternative names for each variable
  integer  :: nSynonym = 3

  ! State variables not associated with a specific fluid/ specie
  integer,parameter  :: nVarExtra = 11

  ! Named indices for all substances (species or fluids)
  integer, parameter :: &
       H_    = 1,  &
       Hp_   = 2,  &
       HpSw_ = 3,  &
       H2p_  = 4,  &
       O_    = 5,  &
       Op_   = 6,  &
       O2p_  = 7,  & 
       He_   = 8,  &
       He2p_ = 9,  &
       OHp_  = 10, &
       N_    = 11, &
       COp_  = 12, &
       CO2p_ = 13, &
       H2O_  = 14, &
       H2Op_ = 15, &
       H3Op_ = 16, &
       Mp_   = 17, &
       Lp_   = 18, &
       MHCp_ = 19, &
       HHCp_ = 20, &
       HNIp_ = 21, &
       Sw_   = 22, &
       Iono_ = 23, &
       Neu1_ = 24, &
       Neu2_ = 25, &
       Neu3_ = 26, &
       Neu4_ = 27, &
       Pui1_ = 28, &
       Pui2_ = 29, &
       Pui3_ = 30, &
       Pui4_ = 31, &
       Main_ = 32 ! main component, MHD/HD

  ! String array storing the standard names of all substances
  character(len = 6) :: NameSubstance_I(nSubstance) = (/ &
       'H   ',  &
       'Hp  ',  &
       'HpSw',  &
       'H2p ',  &
       'O   ',  &
       'Op  ',  &
       'O2p ',  & 
       'He  ',  &
       'He2p',  &
       'OHp ',  &
       'N   ',  &
       'COp ',  &
       'CO2p',  &
       'H2O ',  &
       'H2Op',  &
       'H3Op',  &
       'Mp  ',  &
       'Lp  ',  &
       'MHCp',  &
       'HHCp',  &
       'HNIp',  &
       'Sw  ',  &
       'Iono',  &
       'Neu1',  &
       'Neu2',  &
       'Neu3',  &
       'Neu4',  &
       'Pui1',  &
       'Pui2',  &
       'Pui3',  &
       'Pui4',  &
       '    '  /) ! main component, MHD / HD 
          
  ! named indices for basic state variables associated with a substance
  integer,parameter :: &
       Rho_   = 1, &
       RhoUx_ = 2, &
       RhoUy_ = 3, &
       RhoUz_ = 4, &
       p_     = 5, &
       Ppar_  = 6, &
       Energy_= 7

  ! string array containing basic state variables associated with a substance
  character(len = 6) :: NameSubstanceVar_I(nVarPerSubstance) = (/ &
       'Rho   ', &
       'Mx    ', &
       'My    ', &
       'Mz    ', &
       'P     ', &
       'Ppar  ', &
       'Energy' /)

  ! string arrays containing variables not associated with any substance
  character(len=5) :: NameVarExtra_I(nVarExtra) = (/ &
       'bx   ', &
       'by   ', &
       'bz   ', &
       'pe   ', &
       'te0  ', &
       'ew   ', &
       'eint ', &
       'ehot ', &
       'hyp  ', &
       'sign ', &
       'lperp' /)

 character(len=5) :: NameVarExtraStandardized_I(nVarExtra) = (/ &
       'Bx   ', &
       'By   ', &
       'Bz   ', &
       'Pe   ', &
       'Te0  ', &
       'Ew   ', &
       'Eint ', &
       'Ehot ', &
       'Hyp  ', &
       'Sign ', &
       'Lperp' /)

  ! Array storing standarized variable names for all species / fluids
  character(len = 20),allocatable :: SubstanceStandardName_II(:,:)

  ! Array storing all possible names 
  character(len = 20),allocatable :: Dictionary_III(:, :, :)
  ! -------------------------------------------------------------------------
contains

  subroutine process_var_string(NameVar,  &
       nDensity, nSpeed, nP, nPpar, nWave, nMaterial)

    use ModUtilities,  ONLY: split_string, join_string

    character(len=*), intent(inout) :: NameVar
    integer,intent(out)             :: nDensity, nSpeed, nP, nPpar
    integer,intent(out)             :: nWave, nMaterial

    integer :: nVarName
    integer, parameter:: MaxNameVar = 100
    character(len=20):: NameVar_V(MaxNameVar)
    !-----------------------------------------------------------------------

    call split_string(NameVar, MaxNameVar, NameVar_V, nVarName)

    call process_var_list(nVarName, NameVar_V,  &
         nDensity, nSpeed, nP, nPpar, nWave, nMaterial)

    call join_string(nVarName, NameVar_V(1:nVarName), NameVar)

  end subroutine process_var_string
  !==========================================================================
  subroutine process_var_list(nVarName, NameVar_V,  &
       nDensity, nSpeed, nP, nPpar, nWave, nMaterial)

    use ModUtilities,  ONLY: lower_case

    integer,intent(in)                :: nVarName
    character(len=*), intent(inout)   :: NameVar_V(nVarName)
    integer,intent(out)               :: nDensity, nSpeed, nP, nPpar
    integer,intent(out)               :: nWave, nMaterial

    ! DESCRIPTION:
    ! ------------
    ! 1. Creates standard names and a dictionary for each standard name.
    !    The dictionary only contains the basic hydro quantities for 
    !    different substances. Other quantities, e.g. magnetic field, 
    !    that are not associated with a specific substance are
    !    stored separately. This allows us to avoid searching the complete 
    !    dictionary when it is not needed.

    ! The dictionary is a string array:
    !    Dictionary_III(nSubstance, nVarPerSubstance, nSynonym)
    !
    !    where:
    !    - nSubstance is the number of possible species/ fluids
    !    - nVarPerSubstance enumarates the variables associated with each substance.
    !    - nSynonym is the number of alternative names representing the same
    !              physical quantity, used by different ModEquation files.
    !
    ! 2. Look up the elements of NameVar_V and replace them with standard names
    !    The look up procedure in the dictionary is done by 
    !    call find_substance_replace_name
    !    Once a specific NameVarIn is found to be identical to a dictionary item:
    !    it is replaced with SubstanceStandardName_II(iSubstance,iVarPerSubstance)
    !
    ! 3. The number of fluids and species found are returned by 
    !    nDensity and nSpeed.

    integer                   :: nDistinctSubstanceVar_I(nVarPerSubstance)
    character(len=15)                 :: NameVarIn
    integer                           :: iName, iVar, iSubstanceFound = 0
    logical                           :: IsFoundVar 

    character(len=*), parameter:: NameSub = 'process_var_name'
    ! ------------------------------------------------------------------------
    nDistinctSubstanceVar_I(:) = 0
    nWave = 0
    nMaterial = 0

    ! create standard names and dictionary arrays
    allocate(SubstanceStandardName_II(nSubstance, nVarPerSubstance))
    allocate(Dictionary_III(nSubstance, nVarPerSubstance, nSynonym))
    call create_dictionary

    ! Look up each var name
    NAMELOOP: do iName = 1, nVarName
       ! init search
       IsFoundVar = .false.
       NameVarIn = NameVar_V(iName)
       call lower_case(NameVarIn)

       ! Don't look up in dictionary for: bx, by, bz, EInt, ew, pe, hyp
       do iVar = 1, nVarExtra
          if(NameVarIn == NameVarExtra_I(iVar)) then
             NameVar_V(iName) = NameVarExtraStandardized_I(iVar)
             IsFoundVar = .true.
             CYCLE NAMELOOP
          end if
       end do

       ! check dictionary ( loop over density, momentum. pressure, energy)
       do iVar = 1, nVarPerSubstance 
          call find_substance_replace_name
          if(IsFoundVar) then
             ! Count how many distinct substance variables are present
             nDistinctSubstanceVar_I(iVar) = &
                  nDistinctSubstanceVar_I(iVar) +1 
             CYCLE NAMELOOP
          end if
       end do

       ! variable name may correspond to numbered wave/material
       ! These names are created  in BATSRUS:MH_set_parameters 
       ! and need not be changed
       if (lge(NameVarIn, 'i01') .and. lle(NameVarIn, 'i99')) then
          nWave = nWave + 1
          IsFoundVar = .true.
          CYCLE NAMELOOP
       end if

       if (lge(NameVarIn, 'm1') .and. lle(NameVarIn, 'm9')) then
          nMaterial = nMaterial + 1
          IsFoundVar = .true.
          CYCLE NAMELOOP
       end if

       if(.not. IsFoundVar) then 
          write(*,*) 'ERROR: Var name not in dictionary: ',NameVarIn
          write(*,*) 'Please use standard variable names in ModEquation '// &
               'and recompile:'
          !write(*,*) SubstanceStandardName_II
          write(*,*) ''
          call CON_stop(NameSub//': unknown variable '//NameVarIn)
       end if

    end do NAMELOOP

    nDensity = nDistinctSubstanceVar_I(Rho_)
    nSpeed   = nDistinctSubstanceVar_I(RhoUx_)
    nP       = nDistinctSubstanceVar_I(P_)
    nPpar    = nDistinctSubstanceVar_I(Ppar_)

    deallocate(Dictionary_III)
    deallocate(SubstanceStandardName_II)

  contains
    !==========================================================================
    subroutine find_substance_replace_name

      ! lookup var name in dictionary, replace with standard name

      use ModUtilities,  ONLY: lower_case

      implicit none

      integer             :: iSubstance, iSynonym
      character(len=15)   :: DictionaryItem
      !----------------------------------------------------------------
      do iSubstance = 1, nSubstance 
         do iSynonym = 1, nSynonym
            DictionaryItem = Dictionary_III(iSubstance, iVar, iSynonym)
            if(len_trim(DictionaryItem) > 0) then
               call lower_case(DictionaryItem)
               if(NameVarIn ==  DictionaryItem) then
                  iSubstanceFound = iSubstance
                  IsFoundVar = .true.
                  NameVar_V(iName) = &
                       SubstanceStandardName_II(iSubstanceFound, iVar)
                  RETURN
               end if
            end if
         end do
      end do
    end subroutine find_substance_replace_name

  end subroutine process_var_list
  ! =========================================================================  
  subroutine create_standard_name

    implicit none

    integer   :: iVar, iSubstance
    ! ---------------------------------------------------------------------
 
    ! loop over all possible species/fluids to fill in Name arrays
    do iSubstance = 1, nSubstance
       do iVar = 1, nVarPerSubstance
          SubstanceStandardName_II(iSubstance,iVar) = &
              ''//trim(NameSubstance_I(iSubstance))//NameSubstanceVar_I(iVar)
       
       end do
    end do
    
  end subroutine create_standard_name
  ! =========================================================================
  subroutine create_dictionary

    implicit none

    integer  :: iSubstance
    ! --------------------------------------------------------------------
    Dictionary_III(:,:,:) = ''

    ! first page in dictionary is a 2 by 2 array of standard names
    call create_standard_name
    Dictionary_III(:,:,1) = SubstanceStandardName_II

    !\
    ! fill in alternative names
    !/
    ! The names below are alternative names to the standard names, as
    ! used by existing ModEquation files.
    ! The use of standard names in equation files is encouraged.

    ! Alternative names for energy for all substances
    do iSubstance = 1, nSubstance
       Dictionary_III(iSubstance, Energy_, 2) = &
            ''//trim(NameSubstance_I(iSubstance))//'e'
    end do

    ! main plasma fluid
    Dictionary_III(Main_, RhoUx_,    2) = 'rhoux'
    Dictionary_III(Main_, RhoUy_,    2) = 'rhouy'
    Dictionary_III(Main_, RhoUz_,    2) = 'rhouz'
    
    ! H atoms
    Dictionary_III(H_, Rho_,   2) = 'rhoh'

    ! H+ ions
    Dictionary_III(Hp_, Rho_,   2) = 'h1p'
    Dictionary_III(Hp_, Rho_,   3) = 'hp'
    Dictionary_III(Hp_, RhoUx_, 2) = 'hpux'
    Dictionary_III(Hp_, RhoUy_, 2) = 'hpuy'
    Dictionary_III(Hp_, RhoUz_, 2) = 'hpuz'

    ! H2+ ions
    Dictionary_III(H2p_, Rho_,    2) = 'h2p'

    ! He atoms
    Dictionary_III(He_, Rho_,     2) = 'rhohe'

    ! O atoms
    Dictionary_III(O_, Rho_,      2) = 'rhoo'

    ! O+ ions
    Dictionary_III(Op_, Rho_,   2) = 'op'

    ! O2+ ions
    Dictionary_III(O2p_, Rho_,   2) = 'o2p'
   
    ! CO+ ions
    Dictionary_III(COp_, Rho_,   2) = 'cop'
   
    ! CO2+ ions
    Dictionary_III(CO2p_, Rho_,   2) = 'co2p'

    ! H2O molecules
    Dictionary_III(H2O_, Rho_,    2) = 'rhoh2o'
    
    ! H2O+ ions
    Dictionary_III(H2Op_, Rho_,   2) = 'h2op'
    Dictionary_III(H2Op_, Rho_,   3) = 'rhoh2op'

    ! H3O+ ions
    Dictionary_III(H3Op_, Rho_,   2) = 'h3op'

    ! OH+ ions
    Dictionary_III(OHp_, Rho_,   2) = 'ohp'
   
    ! Saturn fluids
    Dictionary_III(N_,    Rho_,   2) = 'rhon'

    ! Titan ions
    Dictionary_III(Mp_,   Rho_,   2) = 'mp'
    Dictionary_III(Lp_,   Rho_,   2) = 'lp'
    Dictionary_III(MHCp_, Rho_,   2) = 'mhcp'
    Dictionary_III(HHCp_, Rho_,   2) = 'hhcp'
    Dictionary_III(HNIp_, Rho_,   2) = 'hnip'

    ! solar wind
    Dictionary_III(Sw_, Rho_,   2) = 'swhrho'
    Dictionary_III(Sw_, RhoUx_, 2) = 'swhmx'
    Dictionary_III(Sw_, RhoUy_, 2) = 'swhmy'
    Dictionary_III(Sw_, RhoUz_, 2) = 'swhmz'
    Dictionary_III(Sw_, p_,     2) = 'swhp'
    Dictionary_III(Sw_, Energy_,2) = 'swhe'

    Dictionary_III(Sw_, Rho_,   3) = 'rhosw'
    Dictionary_III(Sw_, Energy_,3) = 'swe'

    ! ionosphere
    Dictionary_III(Iono_, Rho_,   2) = 'rhoion'

    ! Outer Heliosphere Pop1 / arbitrary neutral
    Dictionary_III(Neu1_, Rho_,   2) = 'neurho'
    Dictionary_III(Neu1_, RhoUx_, 2) = 'neumx'
    Dictionary_III(Neu1_, RhoUy_, 2) = 'neumy'
    Dictionary_III(Neu1_, RhoUz_, 2) = 'neumz'
    Dictionary_III(Neu1_, p_,     2) = 'neup'
    Dictionary_III(Neu1_, Energy_,2) = 'neue'

    ! Create Alternate names for arbitrary neutral
    Dictionary_III(Neu1_, Rho_,   3) = 'Neu1Rho'
    Dictionary_III(Neu1_, RhoUx_, 3) = 'Neu1Mx'
    Dictionary_III(Neu1_, RhoUy_, 3) = 'Neu1My'
    Dictionary_III(Neu1_, RhoUz_, 3) = 'Neu1Mz'
    Dictionary_III(Neu1_, p_,     3) = 'Neu1P'
    Dictionary_III(Neu1_, Energy_,3) = 'Neu1E'

    ! Outer Heliosphere Pop2 / arbitrary neutral
    Dictionary_III(Neu2_, Rho_,   2) = 'ne2rho'
    Dictionary_III(Neu2_, RhoUx_, 2) = 'ne2mx'
    Dictionary_III(Neu2_, RhoUy_, 2) = 'ne2my'
    Dictionary_III(Neu2_, RhoUz_, 2) = 'ne2mz'
    Dictionary_III(Neu2_, p_,     2) = 'ne2p'
    Dictionary_III(Neu2_, Energy_,2) = 'ne2e'

    ! Outer Heliosphere Pop3 / arbitrary neutral
    Dictionary_III(Neu3_, Rho_,   2) = 'ne3rho'
    Dictionary_III(Neu3_, RhoUx_, 2) = 'ne3mx'
    Dictionary_III(Neu3_, RhoUy_, 2) = 'ne3my'
    Dictionary_III(Neu3_, RhoUz_, 2) = 'ne3mz'
    Dictionary_III(Neu3_, p_,     2) = 'ne3p'
    Dictionary_III(Neu3_, Energy_,2) = 'ne3e'

    ! Outer Heliosphere Pop4 / arbitrary neutral
    Dictionary_III(Neu4_, Rho_,   2) = 'ne4rho'
    Dictionary_III(Neu4_, RhoUx_, 2) = 'ne4mx'
    Dictionary_III(Neu4_, RhoUy_, 2) = 'ne4my'
    Dictionary_III(Neu4_, RhoUz_, 2) = 'ne4mz'
    Dictionary_III(Neu4_, p_,     2) = 'ne4p'
    Dictionary_III(Neu4_, Energy_,2) = 'ne4e'
    
    ! Outer Heliosphere Pop1 / arbitrary pick up ion
    Dictionary_III(Pui1_, Rho_,   2) = 'pu1rho'
    Dictionary_III(Pui1_, RhoUx_, 2) = 'pu1mx'
    Dictionary_III(Pui1_, RhoUy_, 2) = 'pu1my'
    Dictionary_III(Pui1_, RhoUz_, 2) = 'pu1mz'
    Dictionary_III(Pui1_, p_,     2) = 'pu1p'
    Dictionary_III(Pui1_, Energy_,2) = 'pu1e'
   
    ! Outer Heliosphere Pop2 / arbitrary pick up ion
    Dictionary_III(Pui2_, Rho_,   2) = 'pu2rho'
    Dictionary_III(Pui2_, RhoUx_, 2) = 'pu2mx'
    Dictionary_III(Pui2_, RhoUy_, 2) = 'pu2my'
    Dictionary_III(Pui2_, RhoUz_, 2) = 'pu2mz'
    Dictionary_III(Pui2_, p_,     2) = 'pu2p'
    Dictionary_III(Pui2_, Energy_,2) = 'pu2e'

    ! Outer Heliosphere Pop3 / arbitrary pick up ion
    Dictionary_III(Pui3_, Rho_,   2) = 'pu3rho'
    Dictionary_III(Pui3_, RhoUx_, 2) = 'pu3mx'
    Dictionary_III(Pui3_, RhoUy_, 2) = 'pu3my'
    Dictionary_III(Pui3_, RhoUz_, 2) = 'pu3mz'
    Dictionary_III(Pui3_, p_,     2) = 'pu3p'
    Dictionary_III(Pui3_, Energy_,2) = 'pu3e'

    ! Outer Heliosphere Pop4 / arbitrary pick up ion
    Dictionary_III(Pui4_, Rho_,   2) = 'pu4rho'
    Dictionary_III(Pui4_, RhoUx_, 2) = 'pu4mx'
    Dictionary_III(Pui4_, RhoUy_, 2) = 'pu4my'
    Dictionary_III(Pui4_, RhoUz_, 2) = 'pu4mz'
    Dictionary_III(Pui4_, p_,     2) = 'pu4p'
    Dictionary_III(Pui4_, Energy_,2) = 'pu4e'
    
  end subroutine create_dictionary
  
  ! =========================================================================
 
end module ModProcessVarName
