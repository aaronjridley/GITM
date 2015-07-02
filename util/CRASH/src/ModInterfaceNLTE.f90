!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CRASH_ModInterfaceNLTE
  use CRASH_ModMultiGroup, ONLY:nGroup, EnergyGroup_I
  use CRASH_M_EOS,   ONLY: UseNLTE=>UseCrashEos !eos_material.f90
  use CRASH_M_RADIOM,      ONLY: printversion
  implicit none
  PRIVATE !Except
  public:: UseNLTE  !Logical describing if the nonLTE is used 
  public:: read_nlte,check_nlte !reads UseNlte, makes a check if UseNlte==.true.
  public:: printversion
  public:: NLTE_EOS !Full list of the eos function parameters (no pIn)
contains
  subroutine read_nlte
    !Use in PARAM.in:
    !#NLTE
    !T/F               UseNLTE
    !
    use ModReadParam,  ONLY: read_var
    !----------------------
    call read_var('UseNLTE',UseNLTE)
  end subroutine read_nlte
  !====================
  subroutine check_nlte
    use CRASH_M_EOS,   ONLY: SetOptions
    use CRASH_M_expTab,ONLY: exp_tab8
    use ModConst,            ONLY: cHPlanckEV
    use CRASH_M_NLTE,only : ng_rad
    use CRASH_M_RADIOM, only : prep_projE, prepCorrUbar
    logical,save:: DoInit = .true.
    !---------------------
    if(.not.DoInit)return
    DoInit = .false.
    !Initialize NLTE calculations
    call exp_tab8()    !ModExpTable.f90
    call setoptions(brent=.false., EELog=.false.) !eos_material.f90
       
    !\
    ! Coefficients for transforming from the user defined grid to
    ! the refined logrithmic-uniform internal fixed grid
    !/ 
    call prep_projE(EnergyGroup_I(0:nGroup),nGroup) !radiom.f90

    !\
    ! Initialize and calculate some internal arrays
    !/
    call prepCorrUbar()                             !
   
    ng_rad=nGroup

  end subroutine check_nlte
  !==========================
  subroutine NLTE_EOS(& !Full list of the eos function parameters (no pIn)
       iMaterialIn,Rho,&
       TeIn, eTotalIn, eElectronIn,   &
       EoBIn_I,                                              &
       TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
       eElectronOut, pElectronOut, GammaEOut, CvElectronOut, &
       OpacityPlanckOut_I, OpacityRosselandOut_I,            &
       HeatCond, TeTiRelax, Ne, zAverageOut, z2AverageOut,UseERadInput)

    use CRASH_M_EOS,   ONLY: iMaterial, set_kbr !eos_material.f90
    use CRASH_M_NLTE,only : ng_rad,EoBIn, NLTE=>NLTE_EOS, setErad 
    use CRASH_ModEos,ONLY: eos, cAtomicMassCRASH_I, &
                           nZMix_II, cMix_II
    use CRASH_M_localProperties,only : atoNum,atoMass !ModLocalProperties.f90
    use CRASH_M_Radiom,ONLY: UseERadInCalTz=>UseERadInput
    use ModConst
    ! Eos function for single material

    integer, intent(in):: iMaterialIn     ! index of material
    real,    intent(in):: Rho             ! mass density [kg/m^3]
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
    real,    optional, intent(in)  :: eElectronIn  ! internal energu density of electrons

    ! E-over-B ratio for all known groups
    real,    optional, intent(in)  :: EoBIn_I(1:ng_rad)

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
    !For indirect EOS either ERad or ERad/B(Te) is used as inputs
    logical, optional, intent(in) :: UseERadInput
    
    real:: Tz, NAtomic, Te, EIn, TzSi, ZBar    !in eV, cm-3, eV, erg/cm3 
    real:: pNlte, pENlte, pLte, pELte, CvTotal, CvElectron
    real:: DPOverDRho, DPEOverDRho, DPOverDT, DPEOverDT, DTzOverDRho
    !---------------
    !Set iMaterial and dependent variables
    !Initialize
    pLte=0.0; pELte=0.0; peNlte = 0.0; pNlte=0.0
    iMaterial = iMaterialIn

    if(present(UseERadInput))then
       UseERadInCalTz = UseERadInput
    else
       UseERadInCalTz = .false.
    end if

    !Calculate atomic density
    NAtomic = Rho/( cAtomicMassCRASH_I(iMaterial)*cAtomicMass ) & !In 1/m3
         * 1.0e-6                                   ![cm-3] !Convert units

    atomass = cAtomicMassCRASH_I(iMaterial)
   
    atonum  = sum(nZMix_II(:,iMaterial)*&
                      cMix_II(:,iMaterial))
    
    if(present(EoBIn_I))then
       EoBIn(1:nGroup) = EoBIn_I
    else
       EoBIn(1:nGroup)=0.0  !Zero radiation energy
    end if

    call set_kbr(NAtom=NAtomic)
    call setErad(eg_o_bg= EoBIn(1:nGroup),&
                 hnug=EnergyGroup_I(0:nGroup),&          
                 ng=nGroup)

    if(present(TeIn))then
       !Convert temperature to eV
       Te=TeIn * cKToeV
       if(present(EElectronOut).or.present(PElectronOut))then
          !get Tz
          call NLTE(Natom=NAtomic,&
               Te_in=Te,             &
               Zbar_out=zBar, &
               Tz_out=Tz,            &
               Ee_out=EElectronOut,  &
               Pe_out=pENlte,        &
               Cv_out=CvElectronOut)
          if(present(pElectronOut))pElectronOut=pENlte
          !erg/cm3=0.1 J/m3
          pENlte = pENlte * 0.10
       else
          call NLTE(Natom=NAtomic,   &
               Te_in=Te,             &
               Zbar_out=zBar, &
               Tz_out=Tz,            &
               Et_out=ETotalOut,     &
               Pt_out=pNlte,         &
               Cv_out=CvTotalOut)
          if(present(pTotalOut))pTotalOut=pNlte
          !erg/cm3=0.1 J/m3
          pNlte = pNlte * 0.10
       end if
    elseif(present(EElectronIn))then
       !Convert J/m3 = 10^7erg/10^6cm3=10 erg/cm3
       EIn = EElectronIn * 10.0

       !Get Tz
       call NLTE(Natom=NAtomic,&
         Ee_in=EIn,           &
         Zbar_out=zBar,&
         Tz_out=Tz,           &
         Te_out=Te,           &
         Pe_out=pENlte,       &
         Cv_out=CvElectronOut)
       if(present(pElectronOut))pElectronOut=pENlte
       !erg/cm3=0.1 J/m3
       pENlte = pENlte * 0.10
    elseif(present(ETotalIn))then
       !Convert J/m3 = 10^7erg/10^6cm3=10 erg/cm3
        EIn = ETotalIn * 10.0

       !Get Tz

       call NLTE(Natom=NAtomic,&
         Et_in=EIn,           &
         Zbar_out=zBar,       &
         Tz_out=Tz,           &
         Te_out=Te,           &
         Pt_out=pNlte,        &
         Cv_out=CvTotalOut)
       if(present(pTotalOut))pTotalOut=pNlte
       !erg/cm3=0.1 J/m3
          pNlte = pNlte * 0.10
    else
       call CON_stop(&
            'Stop NLTE_eos: TeIn or EElectronIn or ETotalIn must be present')
    end if
    !\
    ! CONVERT
    !/

    !eV to K
    TzSi = Tz*ceVToK
    
    if(present(TeOut))then
       TeOut = Te*ceVToK
    end if
    !erg/cm3=0.1 J/m3
    if(present(EElectronOut))&
       EElectronOut = EElectronOut*0.10
    
    if(present(ETotalOut   ))&
       ETotalOut    = ETotalOut   *0.10

    if(present(PElectronOut))&
       PElectronOut = PElectronOut*0.10

    if(present(PTotalOut   ))&
       PTotalOut    = PTotalOut   *0.10
    if(present(CvElectronOut))&
       CvElectronOut = CvElectronOut*0.10& !erg/cm3=0.1 J/m3
                       *cKToeV             !1/eV=1/K*(K/eV)
    if(present(CvTotalOut))&
       CvTotalOut = CvTotalOut*0.10& !erg/cm3=0.1 J/m3
                       *cKToeV             !1/eV=1/K*(K/eV)
    if(present(zAverageOut))zAverageOut=zBar
    if(&
         present(GammaOut).or.      &
         present(GammaEOut).or.     &
         present(OpacityPlanckOut_I).or. &
         present(OpacityRosselandOut_I).or. &
         present(HeatCond).or.      &
         present(TeTiRelax).or.     &
         present(Ne).or.            &
         present(z2AverageOut) )    &
         call eos(&
         iMaterial=iMaterialIn,     &
         Rho=Rho,                   &
         TeIn=TzSi,                 &
         pTotalOut=pLte,            &
         pElectronOut=pELte,        &
         OpacityPlanckOut_I=OpacityPlanckOut_I,       &
         OpacityRosselandOut_I=OpacityRosselandOut_I, &
         HeatCond=HeatCond,           &
         TeTiRelax=TeTiRelax,         &
         Ne=Ne,                       &
         z2AverageOut=z2AverageOut,   &
         DPOverDRho=DPOverDRho,       &
         DPOverDT=DPOverDT,           &
         DPEOverDRho=DPEOverDRho,     &
         DPEOverDT=DPEOverDT          )
    
    if(present(GammaOut))then
       !Calculate DTzOverDRho
       
       GammaOut = (pLte/pNlte)*DPOverDRho +(NAtomic*1.0e6)*Te*cEV/pNlte*&
            ((DPOverDT-(zBar+1))*DTzOverDRho/Te + (zBar+1)*(1-Tz/Te)+&
            (DPOverDT*Tz/Te + (zBar+1)*(1-Tz/Te))**2/&
            (CvTotal/(NAtomic*1.0e6*cBoltzmann)+1.50*(zBar+1)*(1-Tz/Te)))
    end if
    if(present(GammaEOut))then
       !Calculate DTzOverDRho
       DTzOverDRho = 0.0
       
       GammaOut = (pELte/pENlte)*DPEOverDRho +(NAtomic*1.0e6)*Te*cEV/pENlte*&
            ((DPEOverDT-zBar)*DTzOverDRho/Te + zBar*(1-Tz/Te)+&
            (DPEOverDT*Tz/Te + zBar*(1-Tz/Te) )**2/&
            (CvElectron/(NAtomic*1.0e6*cBoltzmann)+1.50*zBar*(1-Tz/Te)))
    end if
    !Check positivity
    if(present(pElectronOut))then
       if(pElectronOut<0.0)then
          write(*,*)'Negative electron pressure=',pElectronOut,' n/m2'
          call print_input_and_stop
       end if
    end if
    if(present(pTotalOut))then
       if(pTotalOut<=0.0)then
          write(*,*)'Non-positive total pressure=',pTotalOut,' n/m2'
          call print_input_and_stop
       end if
    end if
    contains
      subroutine print_input_and_stop
        if(present(TeIn))then
           write(*,*)'TeIn=',TeIn,' K'
        else
           write(*,*)'TeIn is not present'
        end if
        if(present(eTotalIn))then
           write(*,*)'eTotalIn=',eTotalIn,' J/m3'
        else
           write(*,*)'eTotalIn is not present'
        end if
        if(present(eElectronIn))then
           write(*,*)'eElectronIn=',eElectronIn,' J/m3'
        else
           write(*,*)'eElectronIn is not present'
        end if
        if(present(EoBIn_I))then
           write(*,*)'EoBIn_I=',EoBIn_I
        else
           write(*,*)'EoBIn_I is not present'
        end if
        write(*,*)'NAtomic=', nAtomic*1e6,' 1/m3'
        write(*,*)'Te=',Te,' eV'
        write(*,*)'Tz=',Tz,' eV'
        write(*,*)'zBar=',zBar
        write(*,*)'pNLte  pressure=',pNLte,' n/m2'
        write(*,*)'peNLte pressure=',peNlte,' n/m2'
        write(*,*)'pLte   pressure=',pLte,' n/m2'
        write(*,*)'peLte  pressure=',peLte,' n/m2' 
        write(*,*)'iMaterial=',iMaterialIn
        call CON_stop('Unphysical output parameter(s) from NLTE EOS')
      end subroutine print_input_and_stop
    !TBD  Correct heat conduction and TeTiRelax with (Te/Tz)factors
  
    end subroutine NLTE_EOS
end module CRASH_ModInterfaceNLTE
