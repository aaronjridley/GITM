!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CRASH_ModEosTable

  use CRASH_ModEos
  use ModLookupTable, ONLY: MaxTable
  use CRASH_ModMultiGroup, ONLY: nGroup
  use CRASH_ModInterfaceNlte
  use ModConst

  implicit none

  PRIVATE
  SAVE

  !The following subroutine:
  !1. Checks if the tables are available for these 
  !   materials for which UseEosTable_I is True
  !2. If the table is described in the PARAM.in file and should be
  !   made, it is made  
  !3. If the table is not described in the PARAM.in file form it with
  !   the default ranges and fill in using the built-in EOS model

  public:: check_eos_table

  !The following subroutine:
  !1. Checks if the opacity tables are available for the
  !   materials for which UseOpacityTable_I is True
  !2. If the table is described in the PARAM.in file and should be
  !   made, it is made  
  !3. If the table is not described in the PARAM.in file form it with
  !   the default ranges and fill in using the built-in opacity model

  public:: check_opac_table

  ! Arrays which relate the iTable for the EOS table with 
  ! the material number:
  integer:: iMaterial4EosTable_I(MaxTable) = -1

  ! Arrays which relate the iTable for the opacity table with 
  ! the material number:
  integer:: iMaterial4OpacTable_I(MaxTable) = -1

  !The variable names in the EOS table
  character(LEN=100):: NameVarEos = &
       'P E Pe Ee Cv Cve Gamma GammaE TeTi Cond Z Z2 '//&
       'DPOverDRho DPOverDT DPEOverDRho DPEOverDT'

  !The following subroutine may be used to reset the
  !list of variables tabulated in the EOS tables.
  !Accordingly the named indexes are resen as well as
  !nVarEos. The list may be shortened but not extended.
  !The named indices for the excluded variables are set to
  !-1.

  public:: read_name_var_eos 

  !The number of columns in the Opac table
  integer :: nVarOpac =2   

  !The variable names in the opacity table
  character(LEN=100):: NameVarOpac = &
       'Opac(1) Ross(1)'

  !\
  ! Defaultparameters for EOS tables:
  !/
  integer, parameter :: Min_=1, Max_=2
  integer, public    :: IndexDefaultEos_I(2) = (/201, 201/)
  integer, public    :: IndexDefaultOpac_I(2)= (/201, 201/)
  real,dimension(Min_:Max_,0:nMaterialMax-1), parameter::&
       TeDefaultEos_II = reshape(&     ! original minimum
       (/1.0e-2, 1.0e+3, & !Xe_     1e-2
       1.0e-3, 4.0e+3, & !Be_     1e-3
       1.0e-3, 1.0e+2, & !Plastic 1e-3
       1.0e-3, 1.0e+2, & !Au_     1e-3
       1.0e-3, 1.0e+2  & !Ay_     1e-3
       /), (/2,5/)),    &
       TeDefaultOpac_II = reshape(&     ! original minimum
       (/3.0e-2, 1.0e+3, & !Xe_     1e-2
       3.0e-2, 4.0e+3, & !Be_     1e-3
       5.0e-2, 1.0e+2, & !Plastic 1e-3
       2.0e-1, 1.0e+2, & !Au_     1e-3
       5.0e-2, 1.0e+2  & !Ay_     1e-3
       /), (/2,5/)),    &
       NaDefault_II = reshape(&
       (/1.0e+24, 1.0e+29, & !Xe_
       1.0e+23, 2.0e+29, & !Be_
       1.0e+24, 1.5e+29, & !Plastic
       1.0e+24, 1.2e+29, & !Au_
       1.0e+24, 1.5e+29  & !Ay_
       /), (/2,5/))

  !Note that at 1 Atm and at the room temperature
  !In gas: N~3.10^{25} m-3
  !In solids: N<10^{29} m-3
  ! electron temperature of 1e-3 eV is approximately 11.6 K

contains
  !============================================================================
  ! This subroutine may be used to exclude undesired columns 
  ! from the eos lookup tables
  ! ==========================================================================
  subroutine read_name_var_eos

    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: split_string, lower_case

    integer, parameter:: MaxString = 200
    character(LEN=20):: NameVar_I(MaxString)
    integer:: iVar
    !--------------------------------------------------------------------------

    call read_var('NameVarEos', NameVarEos)
    call split_string(NameVarEos, MaxString,  NameVar_I, nVarEos)

    !Reset named indices
    P_ = -1; E_ = -1; Pe_ = -1; Ee_ = -1
    Cv_ = -1; Cve_ = -1; Gamma_ = -1; GammaE_ = -1
    Cond_ = -1; TeTi_ = -1; Z_ = -1; Z2_ = -1

    do iVar = 1, nVarEos
       call lower_case(NameVar_I(iVar))
       select case(trim(NameVar_I(iVar)))
       case('p')
          P_ = iVar
       case('e')
          E_ = iVar
       case('pe')
          Pe_ = iVar
       case('ee')
          Ee_ = iVar
       case('cv')
          Cv_ = iVar
       case('cve')
          Cve_= iVar
       case('gamma')
          Gamma_= iVar
       case('gammae')
          GammaE_ = iVar
       case('cond')
          Cond_ = iVar
       case('teti', 'tite')
          TeTi_ = iVar
       case('z')
          Z_ = iVar
       case('z2')
          Z2_= iVar
       case('dpdrho')
          DPOverDRho_ = iVar
       case('dpedrho')
          DPEOverDRho_ = iVar
      case('dpdt')
          DPOverDT_ = iVar
       case('dpedt')
          DPEOverDT_ = iVar
       case default
          call CON_stop(NameVar_I(iVar)// &
               ' is not allowed in the eos lookup tables')
       end select
    end do

  end subroutine read_name_var_eos
  !===========================================================================
  subroutine check_eos_table(iComm,TypeFileIn)

    use ModLookupTable, ONLY: Table_I, TableType, &
         i_lookup_table, init_lookup_table, make_lookup_table

    integer, optional, intent(in) :: iComm
    character(LEN=*),optional,intent(in)::TypeFileIn

    integer:: iMaterial, iTable, i
    character(len=2):: NameMaterial
    type(TableType), pointer:: Ptr
    character(len=5)::TypeFile

    character(len=*), parameter:: NameSub = 'check_eos_table'
    !------------------------------------------------------------------------
    if(present(TypeFileIn))then
       TypeFile = TypeFileIn
    else
       TypeFile = 'real8'
    end if
    do iMaterial = 0, nMaterialEos-1

       if(.not.UseEosTable_I(iMaterial))CYCLE
       iTable =  i_lookup_table(NameMaterial_I(iMaterial)//'_eos')

       if(iTable < 0)then
          NameMaterial = NameMaterial_I(iMaterial)

          ! initialize the EOS table with the default parameters
          call init_lookup_table(                                      &
               NameTable = NameMaterial//'_eos',                       &
               NameCommand = 'save',                                   &
               NameVar = 'logTe logNa '//trim(NameVarEos)//' Mass',    &
               nIndex_I = IndexDefaultEos_I,                           &
               IndexMin_I = (/TeDefaultEos_II(Min_, iMaterial),        &
               NaDefault_II(Min_, iMaterial)/),                        &
               IndexMax_I = (/TeDefaultEos_II(Max_, iMaterial),        &
               NaDefault_II(max_, iMaterial)/),                        &
               NameFile = NameMaterial//'_eos_CRASH.dat',              &
               TypeFile = TypeFile,                                    &
               StringDescription = 'CRASH EOS for '//NameMaterial,     &
               nParam = 1,                                             &
               Param_I = (/ cAtomicMassCrash_I(iMaterial) /)   )

          iTable =  i_lookup_table(NameMaterial//'_eos')
       else
          ! Obtain the atomic mass from the lookup table
          Ptr => Table_I(iTable)
          cAtomicMassCrash_I(iMaterial) = 0.0
          if(Ptr%nParam == 1)then
             ! Read the atomic mass from the table parameter
             cAtomicMassCrash_I(iMaterial) = Ptr%Param_I(1)
          end if
          ! For backwards compatibility try reading the description string
          if(cAtomicMassCrash_I(iMaterial) == 0.0)then
             ! Read the atomic mass from the description string
             i = index(Ptr%StringDescription,'Mass =')
             if(i > 0) read(Ptr%StringDescription(i+6:100),*) &
                  cAtomicMassCrash_I(iMaterial)
          end if
          ! Check if the mass was read
          if(cAtomicMassCrash_I(iMaterial) == 0.0) call CON_stop(NameSub// &
               ': atomic mass is missing from table '// &
               trim(Ptr%NameTable))
       end if

       ! The table is now initialized. Set indexes:
       iTableEos4Material_I(iMaterial) = iTable
       iMaterial4EosTable_I(iTable)    = iMaterial

       !Temporary disable the use of table in EOS
       UseEosTable_I(iMaterial) = .false.

       !Fill in the table
       call make_lookup_table(iTable, calc_eos_table, iComm)

       !Recover the true value for UseEos:
       UseEosTable_I(iMaterial) = .true.
    end do

  end subroutine check_eos_table
  !===========================================================================
  subroutine calc_eos_table(iTable, Arg1, Arg2, Value_V)
    integer, intent(in):: iTable
    real, intent(in)   :: Arg1, Arg2
    real, intent(out)  :: Value_V(:)

    real:: ValueTmp_V(-1:nVarEos)
    real:: Rho, Te
    integer::iMaterial
    !--------------------------------------------------------------------------
    iMaterial = iMaterial4EosTable_I(iTable)

    Rho = Arg2 * cAtomicMass * cAtomicMassCRASH_I(iMaterial)
    Te  = Arg1 * cEvToK

    ValueTmp_V = 0.0

    if(.not.UseNLTE)then
       call eos(iMaterial,Rho,&
       TeIn = Te, &
       eTotalOut = ValueTmp_V(E_)    ,&
       pTotalOut = ValueTmp_V(P_)    ,&
       GammaOut  = ValueTmp_V(Gamma_),&
       CvTotalOut= ValueTmp_V(Cv_)   ,&
       eElectronOut = ValueTmp_V(Ee_),&
       pElectronOut = ValueTmp_V(Pe_),&
       GammaEOut    = ValueTmp_V(GammaE_),&
       CvElectronOut = ValueTmp_V(Cve_)  ,& 
       HeatCond      = ValueTmp_V(Cond_) ,&
       TeTiRelax    = ValueTmp_V(TeTi_)  ,&
       zAverageOut   = ValueTmp_V(Z_)    ,&
       z2AverageOut  = ValueTmp_V(Z2_)   ,&
       DPOverDRho    = ValueTmp_V(DPOverDRho_ ), &
       DPOverDT      = ValueTmp_V(DPOverDT_   ), &
       DPEOverDRho   = ValueTmp_V(DPEOverDRho_), &
       DPEOverDT     = ValueTmp_V(DPEOverDT_  ))
       
    else
       call nlte_eos(iMaterial,Rho,&
       TeIn = Te, &
       eTotalOut = ValueTmp_V(E_)    ,&
       pTotalOut = ValueTmp_V(P_)    ,&
       GammaOut  = ValueTmp_V(Gamma_),&
       CvTotalOut= ValueTmp_V(Cv_)   ,& 
       HeatCond      = ValueTmp_V(Cond_) ,&
       TeTiRelax    = ValueTmp_V(TeTi_)  ,&
       zAverageOut   = ValueTmp_V(Z_)    ,&
       z2AverageOut  = ValueTmp_V(Z2_)   )
       call nlte_eos(iMaterial,Rho,&
       TeIn = Te, &
       eElectronOut = ValueTmp_V(Ee_),&
       pElectronOut = ValueTmp_V(Pe_),&
       GammaEOut    = ValueTmp_V(GammaE_),&
       CvElectronOut = ValueTmp_V(Cve_))
    end if
    ValueTmp_V(-1:0) =  0.0

    ValueTmp_V(E_)   =  ValueTmp_V(E_)/(cEV * Arg2)
    ValueTmp_V(P_)   =  ValueTmp_V(P_)/(cEV * Arg2)
    ValueTmp_V(Ee_)  =  ValueTmp_V(Ee_)/(cEV * Arg2)
    ValueTmp_V(Pe_)  =  ValueTmp_V(Pe_)/(cEV * Arg2)
    ValueTmp_V(Cv_)  =  ValueTmp_V(Cv_)/(cBoltzmann * Arg2)
    ValueTmp_V(Cve_) =  ValueTmp_V(Cve_)/(cBoltzmann * Arg2)

    Value_V(1:nVarEos) = ValueTmp_V(1:nVarEos)
  end subroutine calc_eos_table

  !=========================================================================

  subroutine check_opac_table(FreqMinSi, FreqMaxSi, iComm, &
       EGroupIn_I, TypeFileIn)

    use ModConst,       ONLY: cHPlanckEv
    use ModLookupTable, ONLY: Table_I, TableType, &
         i_lookup_table, init_lookup_table, make_lookup_table
    use CRASH_ModMultiGroup,  ONLY: set_multigroup
    use ModMpi,         ONLY: MPI_COMM_RANK
    use ModIoUnit,      ONLY: UnitTmp_

    !Input minimum and maximum frequency (Hz) of photons
    real,    optional, intent(in):: FreqMinSi, FreqMaxSi

    !Input communicator for parallel computations
    integer, optional, intent(in):: iComm

    !Input energy grid, in eV. In case there are only two
    !values they are minimum and maximum energy in electron-volts 
    real,    optional, intent(in):: EGroupIn_I(:)

    !For test the capability to save ascii file should be kept
    character(LEN=*), optional, intent(in):: TypeFileIn

    ! Minimum and maximum group energies in electron volts
    real:: EvMin, EvMax
    real, allocatable:: EGroup_I(:) !Energy Grid [eV]

    integer:: iMaterial, iTable, iProc, iError
    character(len=2):: NameMaterial
    type(TableType), pointer:: Ptr

    character(LEN=5)::TypeFile='real8'

    !Used to split NameVar
    integer:: iPosition
    character(LEN=100)Name,Name1,Name2

    logical:: UseLogarithmicGrid, DoSave
    character:: String1, String2

    character(len=*), parameter:: NameSub = 'check_opac_table'
    !----------------------------------------------------------
   
    if(.not.any(UseOpacityTable_I(0:nMaterialEos-1)))return

    ! Reset if there in an input energy array (0:nGroup)
    UseLogarithmicGrid = .true. 

    ! set iProc for less verbose error messages
    if(present(iComm))then
       call MPI_COMM_RANK(iComm, iProc, iError)
    else
       iProc = 0
    end if

    EvMin = 0.1
    if(present(FreqMinSi)) EvMin = FreqMinSi * cHPlanckEV
    EvMax = 20000.0

    if(present(FreqMaxSi)) EvMax = FreqMaxSi * cHPlanckEV

    if(present(EGroupIn_I))then
       allocate(EGroup_I(size(EGroupIn_I)))
       EGroup_I = EGroupIn_I
       UseLogarithmicGrid = size(Egroup_I) /= (ngroup + 1)
    else
       allocate(EGroup_I(2))
       EGroup_I = (/EvMin, EvMax/)
    end if
       
    if(present(TypeFileIn))TypeFile = TypeFileIn

    ! Construct NameVarOpac
    nVarOpac = 2*nGroup
    NameVarOpac = ''

    if(nGroup == 1)then
       NameVarOpac = 'Planck Ross'
    else
       ! Number of digits needed for nGroup
       write(String1,'(i1)') 1 + int(alog10(real(nGroup)))
       if(UseLogarithmicGrid)then
          write(NameVarOpac,'(a,i'//String1//',a,i'//String1//',a)')&
               'Planck(', nGroup, ') Ross(',nGroup,') EvMin EvMax'
       else
          ! Number of digits needed for nGroup+1
          write(String2,'(i1)') 1 + int(alog10(real(nGroup+1)))
          write(NameVarOpac, &
               '(a,i'//String1//',a,i'//String1//',a,i'//String2//',a)')&
               'Planck(', nGroup, ') Ross(',nGroup,') Ev(',nGroup+1,')'
       end if
    end if

    do iMaterial = 0, nMaterialEos-1
       if(.not.UseOpacityTable_I(iMaterial))CYCLE

       NameMaterial = NameMaterial_I(iMaterial)
       iTable =  i_lookup_table(NameMaterial//'_opac')

       ! Check if table is already set by #LOOKUPTABLE command
       if(iTable < 0)then
          ! initialize the opacity table with the default parameters

          call init_lookup_table(                                       &
               NameTable = NameMaterial//'_opac',                       &
               NameCommand = 'save',                                    &
               NameVar = 'logRho logTe '//NameVarOpac,                  &
               nIndex_I = IndexDefaultOpac_I,                           &
               IndexMin_I = (/NaDefault_II(Min_, iMaterial)*            &
               cAtomicMass * cAtomicMassCRASH_I(iMaterial),             &
               TeDefaultOpac_II(Min_, iMaterial)/),                     &
               IndexMax_I = (/NaDefault_II(Max_, iMaterial)*            &
               cAtomicMass * cAtomicMassCRASH_I(iMaterial),             &
               TeDefaultOpac_II(Max_, iMaterial)/),                     &
               NameFile = NameMaterial//'_opac_CRASH.dat',              &
               TypeFile = TypeFile,                                     &
               StringDescription = 'CRASH Opacity for '//NameMaterial,  &
               nParam = size(EGroup_I),                                 &
               Param_I =  EGroup_I                           )

          ! Get the table index and the pointer to the table
          iTable =  i_lookup_table(NameMaterial//'_opac')
          Ptr => Table_I(iTable)

       else
          Ptr => Table_I(iTable)
          ! Check table for number of opacity values
          if(Ptr%nValue /= 2*nGroup .and. iProc==0)then
             write(*,*)NameSub,' ERROR in table ' // trim(Ptr%NameTable)
             write(*,*)NameSub,': nValue=', Ptr%nValue, &
                  ' should be equal to 2*nGroup=', 2*nGroup
             call CON_stop(NameSub//' change number of groups or table')
          end if

          ! Check that the table parameters contain the energy group boundaries
          if(Ptr%nParam < 2)then
             ! For sake of backward compatibility add energy limits as
             ! the two table parameters here
             Ptr%nParam = 2
             if(allocated(Ptr%Param_I))deallocate(Ptr%Param_I)
             allocate(Ptr%Param_I(2))
             Ptr%Param_I = (/EvMin, EvMax/)
          end if

          if(Ptr%nParam /= 2 .and. Ptr%nParam /= nGroup+1 .and. iProc==0)then
             write(*,*)NameSub,' ERROR in table ' // trim(Ptr%NameTable)
             write(*,*)NameSub,': nParam=', Ptr%nParam, &
                  ' should be equal to 2 or nGroup+1=', nGroup+1
             call CON_stop(NameSub//' change number of groups or table')
          end if
       end if

       ! The table is now initialized. Set indexes:
       iTableOpac4Material_I(iMaterial) = iTable
       iMaterial4OpacTable_I(iTable)    = iMaterial

       ! Set up the energy grid
       if(iMaterial == 0) call set_multigroup(EnergyEv_I=Ptr%Param_I)

       ! Fill in the table. Note: nothing is done if table is loaded from file
       ! The UseOpacityTable_I has to be switched off
       DoSave = Ptr%NameCommand=='save'
       UseOpacityTable_I(iMaterial) = .false.
       call make_lookup_table(iTable, calc_opac_table, iComm)
       UseOpacityTable_I(iMaterial) = .true.

       if(.not.DoSave .or. iProc/=0)CYCLE

       ! Save an include parameter file that documents what the table contains
       ! and allows the user to recreate the table if desired

       open(UnitTmp_,file=NameMaterial//'_opac_CRASH.head')
       write(UnitTmp_,'(a)')'#LOOKUPTABLE'
       write(UnitTmp_,'(a)')NameMaterial//&
            '_opac                                 NameTable'
       write(UnitTmp_,'(a)') &
            'use param                               NameCommand'
       write(UnitTmp_,'(a)')'Tables/'//NameMaterial//'_opac_CRASH.dat'// &
            '                NameFile'
       write(UnitTmp_,'(a)') &
            'real8                                   TypeFile'
       Name      = Ptr%NameVar
       iPosition = index(Name,'Ross')
       iPosition = iPosition + index(Name(iPosition+1:len_trim(Name)),')')
       Name1 = Name(1:iPosition)
       Name2 = Name(iPosition+1:len_trim(Name))
       write(UnitTmp_,'(a)') trim(Name2)//&
            '                                 NameTableParam'
       write(UnitTmp_,&
            '(es13.7,"                           TableParam [eV]")') &
            Ptr%Param_I
       write(UnitTmp_,'(a)')'Opacity(rho,Te) for '//NameMaterial//&
            '                  StringDescription'
       write(UnitTmp_,'(a)') trim(Name1)//'         NameVar'
       write(UnitTmp_,'(i4,a)')Ptr%nIndex_I(1),&
            '                                    nIndex1'
       write(UnitTmp_,'(es13.7,a)')10.0**(Ptr%IndexMin_I(1)),&
            '                           Index1Min [kg/m3]'
       write(UnitTmp_,'(es13.7,a)')10.0**(Ptr%IndexMax_I(1)),& 
            '                           Index1Max [kg/m3]'
       write(UnitTmp_,'(i4,a)')Ptr%nIndex_I(2),&
            '                                    nIndex2'
       write(UnitTmp_,'(es13.7,a)')10.0**(Ptr%IndexMin_I(2)),&
            '                           Index2Min [eV]'
       write(UnitTmp_,'(es13.7,a)')10.0**(Ptr%IndexMax_I(2)),&
            '                           Index2Max [eV]'
       close(UnitTmp_)
    end do

  end subroutine check_opac_table
  !============================================================================
  subroutine calc_opac_table(iTable, Arg1, Arg2, Value_V)
    integer, intent(in):: iTable
    real, intent(in)   :: Arg1, Arg2
    real, intent(out)  :: Value_V(:)

    real:: Rho, Te, PlanckTmp_I(nGroup), RosselandTmp_I(nGroup)
    integer::iMaterial
    !--------------------------------------------------------------------------
    iMaterial = iMaterial4OpacTable_I(iTable)

    Rho = Arg1 
    Te  = Arg2 * cEvToK

    Value_V = 0.0

    if(.not.UseNLTE)then
       call eos(iMaterial,Rho,&
       TeIn = Te, &
       OpacityPlanckOut_I = PlanckTmp_I, &
       OpacityRosselandOut_I = RosselandTmp_I )
    else
       call eos(iMaterial,Rho,&
       TeIn = Te, &
       OpacityPlanckOut_I = PlanckTmp_I, &
       OpacityRosselandOut_I = RosselandTmp_I )
    end if


    Value_V(1:nGroup) = PlanckTmp_I/Arg1
    Value_V(1+nGroup:2*nGroup) = RosselandTmp_I/Arg1

  end subroutine calc_opac_table
  !============================================================================

end module CRASH_ModEosTable
