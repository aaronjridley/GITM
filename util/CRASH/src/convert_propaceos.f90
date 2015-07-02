!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program PROPACEOS
  use  CRASH_ModAtomicNotation
  use CRASH_ModAtomicMass,ONLY: cAtomicMass_I
  use ModConst
  use ModPlotFile
  use CRASH_ModPartition,ONLY: &
       Z2PerA, Z2, ZAv, CoulombLog, Na, Te
  use CRASH_ModAtomicDataMix,ONLY: nMix,Concentration_I, nZ_I
  use CRASH_ModPartition,ONLY: Population_II,iZMin_I
  use CRASH_ModTransport,ONLY: electron_heat_conductivity, te_ti_relaxation
  use CRASH_ModStatSum, ONLY: cZMin
  use ModUtilities, ONLY: split_string
  implicit none

  integer,parameter::iUnit = 11, nDensity=201, nTemperature = 201
  integer,parameter::nFrequency=30
  real,  &
       dimension(nFrequency+1) :: hNu_I

  real,  &
       dimension(nTemperature,nDensity) :: &
       zAvr_II, ETotal_II, dETotalOverDT_II, dETotalOverDRho_II, EIon_II, EElectron_II, &
       dEIonOverDT_II, dEElectronOverDT_II, pTotal_II,&
       pIon_II, pElectron_II, dPIonOverDT_II, dPElectronOverDT_II,&
       dPTotalOverRho_II, dPElectronOverRho_II
  
  real,  &
       dimension(nTemperature,nDensity) :: &
       RossOpacity_II, &
       PlanckOpacity_II, PlanckOpacEms_II
  
  real,  &
       dimension(nTemperature) :: &
       Temperature_I
  real,  &
       dimension(nDensity) :: &
       Density_I, Rho_I     
  
  real,  &
       dimension(nFrequency,nTemperature,nDensity) :: &
       PlanckOpacity_III, RossOpacity_III, PlanckOpacEms_III

  real, allocatable:: Value_VII(:,:,:) 

  character(len=80)  :: StringHeader  
  integer::iString, iRho, iTe
  character(LEN=13)::NameFile
  character(LEN=2)::NameMaterial,NameMixture
  integer:: iMix
  character(LEN=30)::NameToRead

  integer, parameter:: MaxString = 200
  character(LEN=20):: NameVar_I(MaxString)
  integer:: iMaterial
  real::AtomicMass


  !The columns in the EOS table
  integer :: P_      =  1, &
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
             Z2_     = 12

  !The number of columns in the EOS table
  integer :: nVarEos =12
  character(LEN=100):: NameVarEos = &
       'P E Pe Ee Cv Cve Gamma GammaE TeTi Cond Z Z2'
  
  character(LEN=100)::NameDescription
  !Commented out:
  ! logical,parameter:: DoHeaderOnly=.false., DoCut=.false.
  ! integer,parameter:: iStringStart = -1, iStringLast = -1
 
  logical, parameter :: UseLogInterpolation = .true.
  !--------------

  Population_II= 0.0; Concentration_I = 0.0
  nZ_I = 0

  write(*,*)'Enter the name of element (e.g. Ar ) - note than Ar.prp should be present'
  write(*,*)'or the chemical formula (e.g. Ay: C_ 1 H_ 1 ) - note than C_1H_1.prp should be present'

  read(*,'(a)')NameToRead
  if(len_trim(NameToRead)<=2)then
     NameMaterial=trim(NameToRead)
     if(len_trim(NameMaterial)==1)NameMaterial=trim(NameMaterial)//'_'
     iMaterial = i_material(NameMaterial)
     if(iMaterial>0)then
        write(*,*)'iMaterial =', iMaterial
        nZ_I(1) = iMaterial
        AtomicMass = cAtomicMass_I(iMaterial)
        write(*,*)'Atomic weight =',AtomicMass
        NameFile = NameMaterial//'.prp'
        nMix = 1
        iZMin_I = 0
        Concentration_I = 0.0; Concentration_I(1) = 1.0
     else
        call CON_stop('The electron heat conductivity is not supported for mixtures')
     end if
  else
     !\
     ! MIXTURE
     !/
     call split_string( trim(NameToRead), &
               MaxString,  NameVar_I, nMix)
     NameToRead = NameVar_I(1)
     NameMixture = NameToRead(1:2)
     nMix = nMix-1
     do iMix = 1,nMix; NameVar_I(iMix)=NameVar_I(iMix+1);end do
     if(2*(nMix/2)/=nMix)call CON_stop(&
          'For mixtures the formula should be enterd similarly to Wa: H_ 2 O_ 1') 
     nMix = nMix/2
     NameFile='.prp'
     iZMin_I = 0
     Concentration_I = 0.0
     AtomicMass = 0.0
     do iMix=nMix,1,-1
        read(NameVar_I(2*iMix),*)Concentration_I(iMix)
        NameMaterial=trim(NameVar_I(2*iMix-1))
        if(len_trim(NameMaterial)==1)NameMaterial=trim(NameMaterial)//'_'
        iMaterial = i_material(NameMaterial)
        if(iMaterial>0)then
           write(*,*)'iMaterial =', iMaterial
           nZ_I(iMix) = iMaterial
           AtomicMass = AtomicMass + &
                cAtomicMass_I(iMaterial)*&
                Concentration_I(iMix)
           NameToRead = trim(NameFile)
           write(NameFile,'(a,i1,a)') NameMaterial,&
                int(Concentration_I(iMix)),&
                trim(NameToRead)
 
        else
           call CON_stop('Incorrect input for mixtures')
        end if

     end do
     AtomicMass = AtomicMass/sum(Concentration_I)
     Concentration_I = Concentration_I/sum(Concentration_I)
     NameMaterial = NameMixture
     write(*,*)'Atomic weight =',AtomicMass, 'NameMixture=',NameMaterial     
  end if
  write(*,*)'NameFile=',trim(NameFile)
  open(11,file=NameFile,status='old')
  
  !Only prints out the header lines with the numbers. 
  !May be used to find unneeded pieces
  !Commented out!!!!! 
  !if(DoHeaderOnly)then
  !   iString = 0
  !   do 
  !      read(11,'(a)',err=1111,end=1111)StringHeader
  !      iString = iString + 1
  !      if(index(StringHeader,'*')>0&
  !           .and..not.index(StringHeader,'Mean')>0)&
  !           write(*,'(i6,a)')iString,':'//StringHeader
  !
  !   end do
  ! 1111 continue
  !   close (11)
  !   stop
  !end if
  !May be used to delete unneeded piece

  !if(DoHeaderOnly)then
  !   open(12, file = NameMaterial//'.new',status='replace') 
  !   iString = 0
  !   do 
  !      read(11,'(a)',err=2222,end=2222)StringHeader
  !      iString = iString + 1
  !      if(iString<iStringStart.or.iString>iStringLast)&
  !           write(12,'(a)')trim(StringHeader)
  !
  !   end do
!!!!!2222 continue
  !   close (11)
  !   close (12)
  !   stop
  !end if


  do iString =1,24
     read(11,'(a)')StringHeader
     write(*,'(a)')StringHeader
  end do

     call read_eosopa_file_main ( iUnit, &  
          nDensity, nTemperature, &
          nFrequency, &                                        
          hNu_I, zAvr_II, &
          ETotal_II, dETotalOverDT_II, dETotalOverDRho_II, EIon_II, &
          EElectron_II, dEIonOverDT_II, dEElectronOverDT_II, &
          pIon_II, pElectron_II, dPIonOverDT_II, dPElectronOverDT_II, &
          RossOpacity_II, &
          PlanckOpacity_II, &
          PlanckOpacEms_II, &
          Temperature_I, Density_I, &
          PlanckOpacity_III, &
          RossOpacity_III, PlanckOpacEms_III)

  close(11)
  !Rescale Density array

  Density_I = 1.0e6 * Density_I  ! 1/cm3 = 10+6 1/m3
  Rho_I     = cAtomicMass * AtomicMass * Density_I  !Comes in kg/m3 
  
  !============================Opacities===============================
  !Save opacities
  !Swap indexes (Rho<=>Te), merge to a sigle table
  !Convert units: cm2/g=0.10 m2/kg
  
  allocate(Value_VII(2*nFrequency,nDensity,nTemperature))
  do iTe = 1, nTemperature
     do iRho = 1, nDensity
        Value_VII(1:nFrequency, iRho, iTe) = PlanckOpacity_III(:,iTe,iRho) * 0.10
        Value_VII(1+nFrequency:2*nFrequency, iRho, iTe) = RossOpacity_III(:,iTe,iRho) * 0.10
     end do
  end do

  call save_plot_file( &
         NameMaterial//'_opac_PRISM.dat',                                &
         TypeFileIn     = 'real8',                     &
         StringHeaderIn = 'PROPACEOS Opacity for'//NameMaterial, &
         NameVarIn      = 'logRho logTe Planck(30) Ross(30) EvMin EvMax',&
         ParamIn_I      = (/0.1, 2.0e4/), &
         CoordMinIn_D   = (/log10(Rho_I(1)), log10(Temperature_I(1))/),             &                             
         CoordMaxIn_D   = (/log10(Rho_I(nDensity)), log10(Temperature_I(nTemperature))/),             &
         VarIn_VII      = Value_VII)
  deallocate(Value_VII)
  
  !=========================EOS====================
  !Convert energy in J/g to J/kg, 1/g=10+3 1/kg
  
  ETotal_II = 1.0e3* ETotal_II
  EElectron_II = 1.0e3* EElectron_II

  !Convert from J/kg to J/m3
  do iRho = 1, nDensity
     ETotal_II(:,iRho) = ETotal_II(:,iRho) * Rho_I(iRHo)
     EElectron_II(:,iRho) = EElectron_II(:,iRho) * Rho_I(iRHo)
  end do

  !Convert pressure from dyne/cm2 to Pa: dyne/cm2=0.10 * Pa:
  pElectron_II = pElectron_II *0.10
  pIon_II      = pIon_II      *0.10

  !Check positivity!
  write(*,*) 'Total energy is not positive in ',count(ETotal_II<=0.0), ' points'
  write(*,*) 'Electron energy is not positive in ',count(EElectron_II<=0.0), ' points'
  write(*,*) 'Ion pressure is not positive in ',count(pIon_II<=0.0), ' points'
  write(*,*) 'Electron pressure is not positive in ',count(pElectron_II<=0.0), ' points'
  

  pTotal_II =  pIon_II +  pElectron_II
  write(*,*) 'Total pressure is not positive in ',count(pTotal_II<=0.0), ' points'

  if(.not.UseLogInterpolation)then
     !Convert specific heat from J/g/eV to J/m^3/K
     !Convert energy in J/g to J/kg, 1/g=10+3 1/kg
     !1/eV = 1/ceVToK 
     dETotalOverDT_II    = 1.0e3 * dETotalOverDT_II   /cEVToK
     dEElectronOverDT_II = 1.0e3 * dEElectronOverDT_II/cEVToK
     
     !Convert from J/kg to J/m3
     do iRho = 1, nDensity
        dETotalOverDT_II(:,iRho) = dETotalOverDT_II(:,iRho) * Rho_I(iRHo)
        dEElectronOverDT_II(:,iRho) = dEElectronOverDT_II(:,iRho) * Rho_I(iRHo)
     end do
     
     !Pressure derivatives
     dPTotalOverRho_II(:,1) = Rho_I(1)/ pTotal_II(:,1) * &
          (pTotal_II(:,2) - pTotal_II(:,1))/(Rho_I(2)-Rho_I(1))
     
     dPElectronOverRho_II(:,1) = Rho_I(1)/ pElectron_II(:,1) * &
          (pElectron_II(:,2) - pElectron_II(:,1))/(Rho_I(2)-Rho_I(1))
     
     do iRho = 2, nDensity -1
        dPTotalOverRho_II(:,iRho) = Rho_I(iRho)/ pTotal_II(:,iRho) * &
             (pTotal_II(:,iRho+1) - pTotal_II(:,iRho-1))/&
             (Rho_I(iRho+1)-Rho_I(iRho -1))
        
        dPElectronOverRho_II(:,iRho) = Rho_I(iRho)/ pElectron_II(:,iRho) * &
             (pElectron_II(:,iRho+1) - pElectron_II(:,iRho-1))/&
             (Rho_I(iRho+1)-Rho_I(iRho -1))
     end do
     
     dPTotalOverRho_II(:,nDensity) = Rho_I(nDensity)/ pTotal_II(:,nDensity) * &
          (pTotal_II(:,nDensity) - pTotal_II(:,nDensity-1))/&
          (Rho_I(nDensity)-Rho_I(nDensity-1))
     
     dPElectronOverRho_II(:,nDensity) = Rho_I(nDensity)/ pElectron_II(:,nDensity) * &
          (pElectron_II(:,nDensity) - pElectron_II(:,nDensity-1))/&
          (Rho_I(nDensity)-Rho_I(nDensity-1))
     
     !dyne/cm2 /eV = 0.1/ceVToK 1/K
     dPIonOverDT_II = dPIonOverDT_II * 0.10/cEVToK
     dPElectronOverDT_II = dPElectronOverDT_II * 0.10/cEVToK
  else
     call d_over_dt(In_II=ETotal_II,Out_II=dETotalOverDT_II,IsElectronV=.false.)
     call d_over_dt(In_II=EElectron_II,Out_II=dEElectronOverDT_II,IsElectronV=.true.)

     call d_over_dn(In_II=PTotal_II, Out_II=dPTotalOverRho_II,IsElectronV=.false.)
     call d_over_dn(In_II=PElectron_II, Out_II=dPElectronOverRho_II,IsElectronV=.true.)

     call d_over_dt(In_II=PIon_II,Out_II=dPIonOverDT_II,IsElectronV=.false.)
     call d_over_dt(In_II=PElectron_II,Out_II=dPElectronOverDT_II,IsElectronV=.true.)
  end if
  
    !Check positivity!
  write(*,*) 'Total specific heat is not positive in ',&
       count(dETotalOverDT_II<=0.0.and.zAvr_II>cZMin), ' points'
  write(*,*) 'Electron specific heat is not positive in ',&
       count(dEElectronOverDT_II<=0.0.and.zAvr_II>cZMin), ' points'
  !write(*,*) 'Ion pressure is not positive in ',count(pIon_II<=0.0), ' points'
  !write(*,*) 'Electron pressure is not positive in ',count(pElectron_II<=0.0), ' points'

  allocate( Value_VII(nVarEos, nTemperature, nDensity) ) 
  Value_VII = 0.0
  do iRho = 1, nDensity
     do iTe = 1, nTemperature
        Value_VII(E_, iTe, iRho) = ETotal_II( iTe, iRho)/(cEV * Density_I(iRho))

        Value_VII(P_, iTe, iRho) = pTotal_II( iTe, iRho)/(cEV * Density_I(iRho))

        Value_VII(Ee_, iTe, iRho) = ETotal_II( iTe, iRho)/(cEV * Density_I(iRho))

        Value_VII(Pe_, iTe, iRho) = pElectron_II( iTe, iRho)/(cEV * Density_I(iRho))
        
        Value_VII(Cv_, iTe, iRho) = dETotalOverDT_II(iTe, iRho)/(cBoltzmann * Density_I(iRho))
        Value_VII(Cve_, iTe, iRho) = dEElectronOverDT_II(iTe, iRho)/(cBoltzmann * Density_I(iRho))

        Value_VII(Z_,  iTe, iRho) = zAvr_II(iTe, iRho)
        Value_VII(Z2_,  iTe, iRho) = zAvr_II(iTe, iRho)**2
        ZAv =  Value_VII(Z_,  iTe, iRho)
        Z2  =  ZAv **2
        Z2PerA = Z2/AtomicMass
        Na     = Density_I(iRho)
        Te     = Temperature_I(iTe)
 
        !Concentration of neutrals is switched off at 2 eV 
        Population_II(0,1:nMix) = max(0.0, min(1.0,2.0 - Te))
        Value_VII(Cond_,  iTe, iRho) = electron_heat_conductivity()
        Value_VII(TeTi_,  iTe, iRho) = te_ti_relaxation()

        !Use formula (24) from HEDP.pdf:
        !write(*,*)ZAv,dPElectronOverRho_II(iTe, iRho),&
        !     dPElectronOverDT_II( iTe, iRho),pElectron_II(iTe, iRho),&
        !     dEElectronOverDT_II(iTe, iRho),iTe,iRho, EElectron_II(iTe,iRho), EElectron_II(iTe+1,iRho)
        if(zAv > cZMin) then
           Value_VII(GammaE_,  iTe, iRho) = &
             dPElectronOverRho_II(iTe, iRho) + &
             dPElectronOverDT_II( iTe, iRho)**2&
             *Temperature_I(iTe) * cEVToK /(pElectron_II(iTe, iRho)*&
             dEElectronOverDT_II(iTe, iRho))
        else
           Value_VII(GammaE_,  iTe, iRho) = 5.0/3.0
        end if
           
       ! write(*,*)iRho,iTe,pTotal_II(iTe, iRho),&
       !      dETotalOverDT_II(iTe, iRho),dEElectronOverDT_II(iTe, iRho),dEIonOverDT_II(iTe,iRho)
        Value_VII(Gamma_,  iTe, iRho) = &
             dPTotalOverRho_II(iTe, iRho) + &
             (dPElectronOverDT_II( iTe, iRho) +&
             dPIonOverDT_II(iTe, iRho))**2&
             *Temperature_I(iTe) * cEVToK /(pTotal_II(iTe, iRho)*&
             dETotalOverDT_II(iTe, iRho))

        ! GammaSe = compressibility_at_const_te_e() + &
        !    d_pressure_e_over_d_te()**2 * (Na * Te * cEV) /(heat_capacity_e() * pressure_e())
        ! GammaS =  compressibility_at_const_te  () + &
        !    d_pressure_over_d_te  ()**2 * (Na * Te * cEV) /(heat_capacity  () * pressure  ())
        ! Where:
        !  compressibility_at_const_t = 
        !(Na /P) \left(\partial P / \partial Na \right)_{T=const}
        !   the temperuture derivative of the specific pressure (P/Na)
        ! (1/Na)*(\partial P/\partial T)_{\rho=const} 
        ! Heat capacity - dimensionless, per atom
     end do
  end do

  write(NameDescription,'(a,e13.7)')&
       'PROPACEOS EOS for '//NameMaterial//&
       'Atomic Mass = ',AtomicMass
  call save_plot_file( &
       NameMaterial//'_eos_PRISM.dat',                                &
       TypeFileIn     = 'real8',                     &
       StringHeaderIn = trim(NameDEscription), &
       NameVarIn      = 'logTe logNa '//NameVarEos, &
       CoordMinIn_D   = (/log10(Temperature_I(1)),log10(Density_I(1))/),             &                             
       CoordMaxIn_D   = (/log10(Temperature_I(nTemperature)),log10(Density_I(nDensity))/),&
       VarIn_VII      = Value_VII)
  open(11,file='PARAM.'//NameMaterial//'.PRISM',status='replace')
  write(11,'(a)')'----------------EOS TABLE-------------'
  write(11,'(a)')'    '
  write(11,'(a)')'#LOOKUPTABLE'
  write(11,'(a)')NameMaterial//'_eos                     NameTable'
  write(11,'(a)')'load                  NameCommand'
  write(11,'(a)')'Tables/'//NameMaterial//'_eos_PRISM.dat'
  write(11,'(a)')'real8                 TypeFile'
  write(11,'(a)')'EOS(Te, Na) for '//NameMaterial
  write(11,'(a)')'logTe logNa P E Pe Ee Cv Cve Gamma GammaE TeTi Cond Z Z2'
  write(11,'(a)')'201                   nIndex1'
  write(11,'(e13.7,a)')Temperature_I(1),'                       Index1Min (eV)'
  write(11,'(e13.7,a)')Temperature_I(nTemperature),&
                                        '                       Index1Max (eV)'
  write(11,'(a)')'201                   nIndex2'
  write(11,'(e13.7,a)')Density_I(1),    '                       Index2Min (m-3)'
  write(11,'(e13.7,a)')Density_I(nDensity),&
                                        '                       Index2Max (m-3)'
  write(11,'(a)')'  '
  write(11,'(a)')'----------------OPACITY TABLE----------'
  write(11,'(a)')'  '
  write(11,'(a)')'#LOOKUPTABLE'
  write(11,'(a)')NameMaterial//'_opac                   NameTable'
  write(11,'(a)')'load                  NameCommand'
  write(11,'(a)')'Tables/'//NameMaterial//'_opac_PRISM.dat'
  write(11,'(a)')'real8                 TypeFile'
  write(11,'(a)')'Opacity(rho,Te) for '//NameMaterial
  write(11,'(a)')'logrho logT Planck(30) Ross(30)'
  write(11,'(a)')'201                   nIndex1'
  write(11,'(e13.7,a)')Rho_I(1),        '                       Index1Min (kg/m3)'
  write(11,'(e13.7,a)')Rho_I(nDensity), '                       Index1Max (kg/m3)'
  write(11,'(a)')'201                   nIndex2'
  write(11,'(e13.7,a)')Temperature_I(1),'                       Index2Min (eV)'
  write(11,'(e13.7,a)')Temperature_I(nTemperature),&
                                        '                       Index2Max (eV)'
  write(11,'(a)')'  '
  write(11,'(a)')'#END'
  write(11,'(a)')'  '
  close(11)


contains
  subroutine d_over_dt(In_II,Out_II,IsElectronV)
    real, intent(in)  :: In_II(  nTemperature, nDensity)
    real, intent(out) :: Out_II( nTemperature, nDensity)
    logical, intent(in)::IsElectronV
    integer:: iTe,iRho
    !-----------------
    do iRho = 1, nDensity
       !Check if the input values can be used 
       if (zAvr_II(1, iRho)<=cZMin .and. IsElectronV) then
          !The point itself is fake
          Out_II(1 , iRho) = 0.0
       elseif(zAvr_II(2, iRho)<=cZMin.and.IsElectronV)then
          Out_II(1 , iRho) = In_II(1, iRho) /&
                             Temperature_I(1)
       else
          Out_II(1 , iRho) = ( In_II(2, iRho)   + In_II(1, iRho) )/&
                          ( Temperature_I(2) + Temperature_I(1))*&
                          ( log(In_II(2, iRho)) - log(In_II(1, iRho)))/&
                          (log(Temperature_I(2)) - log(Temperature_I(1)))
       end if
          
          
          
       do iTe = 2, nTemperature-1
          !Check if the input values can be used 

          if (zAvr_II(iTe, iRho)<=cZMin.and.IsElectronV) then
             !The point itself is fake
             Out_II(iTe , iRho) = 0.0
          elseif(&
               (zAvr_II(iTe-1, iRho) > cZMin .and. zAvr_II(iTe+1, iRho) > cZMin) &
               .or.(.not.IsElectronV) )then

             !Both neighboring points used to calculate the derivatives are good

             Out_II(iTe , iRho) = In_II(iTe, iRho)/&
                          Temperature_I(iTe)*&
               ( log(In_II(iTe+1, iRho)) - log(In_II(iTe-1, iRho)))/&
               (log(Temperature_I(iTe+1)) - log(Temperature_I(iTe-1)))
          elseif(zAvr_II(iTe-1, iRho) > cZMin)then

             !Only iTe+1 point is bad, do not use it

             Out_II(iTe , iRho) = In_II(iTe, iRho)/&
                          Temperature_I(iTe)*&
               ( log(In_II(iTe, iRho)) - log(In_II(iTe-1, iRho)))/&
               (log(Temperature_I(iTe)) - log(Temperature_I(iTe-1)))
          elseif(zAvr_II(iTe+1, iRho) > cZMin)then

             !Only iTe-1 point is bad, do not use it

             Out_II(iTe , iRho) = In_II(iTe, iRho)/&
                          Temperature_I(iTe)*&
               ( log(In_II(iTe+1, iRho)) - log(In_II(iTe, iRho)))/&
               (log(Temperature_I(iTe+1)) - log(Temperature_I(iTe)))
          else
             Out_II(iTe , iRho) = In_II(iTe, iRho) /&
                             Temperature_I(iTe)
          end if
       end do
       !Check if the input values can be used 
       if (zAvr_II(nTemperature, iRho)<=cZMin .and. IsElectronV) then
          !The point itself is fake
          Out_II(nTemperature , iRho) = 0.0
       elseif(zAvr_II(nTemperature-1, iRho)<=cZMin.and.IsElectronV)then
          Out_II(nTemperature , iRho) = In_II(nTemperature, iRho) /&
                             Temperature_I(nTemperature)
       else
          Out_II(nTemperature , iRho) = ( In_II(nTemperature, iRho)   + &
                                     In_II(nTemperature-1, iRho) )/&
                          ( Temperature_I(nTemperature) + &
                                       Temperature_I(nTemperature-1))*&
                          ( log(In_II(nTemperature, iRho)) -&
                                  log(In_II(nTemperature-1, iRho)))/&
                          (log(Temperature_I(nTemperature)) - &
                                      log(Temperature_I(nTemperature-1)))
       end if
    end do
    Out_II=Out_II/cEVToK
  end subroutine d_over_dt
  
  !-----------------------!
  
  subroutine d_over_dn(In_II, Out_II, IsElectronV)
    real, intent(in)  :: In_II(  nTemperature, nDensity)
    real, intent(out) :: Out_II( nTemperature, nDensity)
    logical,intent(in) :: IsElectronV
    integer:: iRho,iTe
    do iTe=1, nTemperature
       !Check if the input values can be used 

       if (zAvr_II(iTe, 1)<=cZMin .and. IsElectronV) then
          !The point itself is fake
          Out_II(iTe , 1) = 1.0
       elseif(zAvr_II(iTe, 2)<=cZMin.and.IsElectronV)then
          Out_II(iTe , 1) = 1.0
       else
          Out_II(iTe,1)=( log(In_II(iTe,2))  - log(In_II(iTe,1)) )/&
                ( log(Density_I(2)) - Log(Density_I(1)) )
       end if
    end do
    do iRho = 2, nDensity-1
       do iTe=1, nTemperature

          if (zAvr_II(iTe, iRho)<=cZMin.and.IsElectronV) then
             !The point itself is fake
             Out_II(iTe , iRho) = 1.0
          elseif(&
               (zAvr_II(iTe, iRho-1) > cZMin .and. zAvr_II(iTe, iRho+1) > cZMin) &
               .or.(.not.IsElectronV) )then

             !Both neighboring points used to calculate the derivatives are good
          
             Out_II(iTe,iRho)=( log(In_II(iTe,iRho+1))  - log(In_II(iTe,iRho-1)) )/&
                  ( log(Density_I(iRho+1)) - Log(Density_I(iRho-1)))
          elseif(zAvr_II(iTe, iRho-1) > cZMin )then
             
             !Only iRho+1 point is bad, ignore it
          
             Out_II(iTe,iRho)=( log(In_II(iTe,iRho))  - log(In_II(iTe,iRho-1)) )/&
                  ( log(Density_I(iRho)) - Log(Density_I(iRho-1)))
          elseif(zAvr_II(iTe, iRho+1) > cZMin )then
             
             !Only iRho-1 point is bad, ignore it
          
             Out_II(iTe,iRho)=( log(In_II(iTe,iRho+1))  - log(In_II(iTe,iRho)) )/&
                  ( log(Density_I(iRho+1)) - Log(Density_I(iRho)))
          else
             Out_II(iTe , iRho) = 1.0
          end if
       end do
    end do
    do iTe=1, nTemperature
       if (zAvr_II(iTe, nDensity)<=cZMin .and. IsElectronV) then
          !The point itself is fake
          Out_II(iTe , nDensity) = 1.0
       elseif(zAvr_II(iTe, nDensity-1)<=cZMin.and.IsElectronV)then
          Out_II(iTe , nDensity) = 1.0
       else
          Out_II(iTe,nDensity)=( log(In_II(iTe,nDensity))  - log(In_II(iTe,nDensity-1)) )/&
                ( log(Density_I(nDensity)) - Log(Density_I(nDensity-1)))
       end if
    end do
  end subroutine d_over_dn
end program PROPACEOS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                                  +
! + Copyright (c) Prism Computational Sciences, Inc., 1998-2000.     +
! +                                                                  +
! + File:         read.f90                          +
! +                                                                  +
! + Written by:   J. J. MacFarlane                                   +
! +                                                                  +
! + Created:      7/1/99                                             +
! + Modified:                                                        +
! + Modified:                                                        +
! +                                                                  +
! + Description:  Reads in equation of state and opacity data.       +
! +               The opacity data file contains:                    +
! +               electron and ion internal energies and derivatives,+
! +               electron and ion pressure and derivatives,         +
! +               mean charge states, and                            +
! +               multigroup opacities.                              +
! +                                                                  +
! + Note: All real variables and arrays in output other than hnu_in  +
! +       are small_real_kind                                        +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine read_eosopa_file_main ( iUnit, &  
     nDensity, nTemperature, &
     nFrequency, &                                        
     hNu_I, zAvr_II, &
     ETotal_II, dETotalOverDT_II, dETotalOverDRho_II, EIon_II, &
     EElectron_II, dEIonOverDT_II, dEElectronOverDT_II, &
     pIon_II, pElectron_II, dPIonOverDT_II, dPElectronOverDT_II, &
     RossOpacity_II, &
     PlanckOpacity_II, &
     PlanckOpacEms_II, &
     Temperature_I, Density_I, &
     PlanckOpacity_III, &
     RossOpacity_III, PlanckOpacEms_III)


  ! ... use statements


  ! ... implicit statements
  implicit none

  !
  ! ... input variables:
  integer, &
       intent(in) :: iUnit
  integer, &
       intent(in) :: nDensity
  integer, &
       intent(in) :: nTemperature
  integer, &
       intent(in) :: nFrequency     
 

  !
  ! ... output variables

  real, intent(out), &
       dimension(nFrequency+1) :: hNu_I

  real, intent(out), &
       dimension(nTemperature,nDensity) :: &
       zAvr_II, ETotal_II, dETotalOverDT_II, dETotalOverDRho_II, EIon_II, EElectron_II, &
       dEIonOverDT_II, dEElectronOverDT_II, &
       pIon_II, pElectron_II, dPIonOverDT_II, dPElectronOverDT_II

  real, intent(out), &
       dimension(nTemperature,nDensity) :: &
       RossOpacity_II, &
       PlanckOpacity_II, PlanckOpacEms_II

  real, intent(out), &
       dimension(nTemperature) :: &
       Temperature_I
  real, intent(out), &
       dimension(nDensity) :: &
       Density_I     

  real, intent(out), &
       dimension(nFrequency,nTemperature,nDensity) :: &
       PlanckOpacity_III, RossOpacity_III, PlanckOpacEms_III


  ! ... type declaration statements:

  integer :: num_freq_grps_local
  integer,parameter::TableInID_=3

  character(len=80)  :: StringHeader


  integer :: i, nfgp1, l, it, id, ig, m, &
       n_temp_opc_tbl, n_dens_opc_tbl, &
       n_temp_eos_tbl, n_dens_eos_tbl

  real :: rho0es
  
  logical,parameter:: DoReadOpacity = .true.

  ! *******************************************************************
  !
  !                          begin execution
  !
  ! *******************************************************************

  ! ... initializations

  hNu_I(:) = 0.0

  zAvr_II(:,:) = 0.0
  ETotal_II(:,:) = 0.0
  dETotalOverDT_II(:,:) = 0.0
  dETotalOverDRho_II(:,:) = 0.0
  EIon_II(:,:) = 0.0
  EElectron_II(:,:) = 0.0
  dEIonOverDT_II(:,:) = 0.0
  dEElectronOverDT_II(:,:) = 0.0
  pIon_II(:,:) = 0.0
  pElectron_II(:,:) = 0.0
  dPIonOverDT_II(:,:) = 0.0
  dPElectronOverDT_II(:,:) = 0.0
  RossOpacity_II(:,:) = 0.0
  PlanckOpacity_II(:,:) = 0.0
  PlanckOpacEms_II(:,:) = 0.0
  Temperature_I(:) = 0.0
  Density_I(:) = 0.0
  

  if ( ( TableInID_ .eq. 1 ) .or. & 
       ( TableInID_ .eq. 2 ) .or. & 
       ( TableInID_ .eq. 3 ) .or. &
       ( TableInID_ .eq. 4 ) .or. &
       ( TableInID_ .eq. 5 ) .or. &
       ( TableInID_ .eq. 6 ) .or. &
       ( TableInID_ .eq. 7 ) ) then

     ! ...    format is already checked in get_opacity_table_grid


     ! ...    read temperature and density grid parameters,
     !        and number of frequency groups
     !        ---------------------------------------------

     ! ...    first, the EOS grid parameters



     read (iUnit,*) n_temp_eos_tbl
     write(*,*)'n_temp_eos_tbl=',n_temp_eos_tbl
     
     read (iUnit,*) (Temperature_I(it),it=1,n_temp_eos_tbl)
     write(*,*)'Temperature range:',Temperature_I(1),Temperature_I(n_temp_eos_tbl)
     read (iUnit,*) n_dens_eos_tbl
     write(*,*)'n_dens_eos_tbl=',n_dens_eos_tbl
     read (iUnit,*) (Density_I(it),it=1,n_dens_eos_tbl)
      write(*,*)'Density range:',Density_I(1),Density_I(n_dens_eos_tbl)
     read (iUnit,*) rho0es
     !write(*,*)rho0es
     



     if ( DoReadOpacity ) then
        do i=1,4
           read (iUnit,802) StringHeader
           write(*,*)StringHeader
        enddo
     end if



     ! ...    next, the opacity grid parameters

     if ( DoReadOpacity ) then


        read (iUnit,*) n_temp_opc_tbl
        write(*,*) 'n_temp_opc_tbl=',n_temp_opc_tbl
        read (iUnit,*) (Temperature_I(it),it=1,n_temp_opc_tbl)
        write(*,*)'Temperature range:',Temperature_I(1),Temperature_I(n_temp_opc_tbl)
        read (iUnit,*) n_dens_opc_tbl
        write(*,*)'n_dens_opc_tbl=',n_dens_opc_tbl
        read (iUnit,*) (Density_I(it),it=1,n_dens_opc_tbl)
        write(*,*)'Density range:',Density_I(1),Density_I(n_dens_opc_tbl)
        read (iUnit,*) num_freq_grps_local
        write(*,*)'nFrequency=',num_freq_grps_local
        
     end if



     ! ...    radiation transport is assumed to be a multi-group model;
     !        read in the frequency group boundaries (eV)
     !        ---------------------------------------------------------

     if ( DoReadOpacity ) then

        nfgp1 = num_freq_grps_local + 1


        read (iUnit,802) StringHeader
        write(*,*)StringHeader
        read (iUnit,*) (hNu_I(l),l=1,nfgp1)
        write(*,*)'Frequency range:', hNu_I(1),hNu_I(nfgp1)
       
     end if

  else


     return

  endif


  ! ... read in the charge states, ion and electron specific energies,
  !     ion and electron pressures, and their derivatives wrt temperature;
  !     store in log form
  !     ------------------------------------------------------------------

 

  ! Zbar (esu)
  !write(*,*)'Start to read zAvr'
  read (iUnit,802) StringHeader
  read (iUnit,*) ((zAvr_II(l,m),  l=1,n_temp_eos_tbl), &
       m=1,n_dens_eos_tbl)
  write(*,*)StringHeader
 

  if ( DoReadOpacity ) then
     ! ... added frequency integrated mean opacities. prp format ID = 3 (November 2003 PRW)
     if ( TableInID_ .gt. 2 ) then
        ! Frequency integrated Rosseland group opacity (cm2/g)
   !     write(*,*)'Start to read integral Rosseland opac'
        read (iUnit,802) StringHeader
        read (iUnit,*) ((RossOpacity_II(l,m), l=1,n_temp_opc_tbl), &
             m=1,n_dens_opc_tbl)
        write(*,*)StringHeader


        ! Frequency integrated Planck group absorption opacity (cm2/g)
        write(*,*)'Start to read integral Planck opac'
        read (iUnit,802) StringHeader
        read (iUnit,*) ((PlanckOpacity_II(l,m), l=1,n_temp_opc_tbl), &
             m=1,n_dens_opc_tbl)
        write(*,*)StringHeader

        ! Frequency integrated Planck group emission opacity (cm2/g)
        write(*,*)'Start to read integral Planck emissivity'
        read (iUnit,802) StringHeader
        read (iUnit,*) ((PlanckOpacEms_II(l,m), l=1,n_temp_opc_tbl), &
             m=1,n_dens_opc_tbl)      
        write(*,*)StringHeader
         
     end if
     ! ... end add

  end if


  ! Total internal energy (J/g)
  !write(*,*)'Start to read internal energy density'
  read (iUnit,802) StringHeader
  read (iUnit,*) ((ETotal_II(l,m), l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
  write(*,*)StringHeader

  if ( TableInID_ .le. 5 ) then
     ! write(*,*)'Start to read internal energy density T-derivative'
     ! Total energy T-derivative (J/g/eV)         
     read (iUnit,802) StringHeader
     read (iUnit,*) ((dETotalOverDT_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

    ! write(*,*)'Start to read internal energy density scaled rho-derivative'

     ! Scaled total energy rho-derivative (1/eV)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((dETotalOverDRho_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

     !write(*,*)'Start to read ion energy density'
     ! Ion internal energy (J/g)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((EIon_II(l,m), l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

    ! write(*,*)'Start to read electron energy density'
     ! Electron internal energy (J/g)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((EElectron_II(l,m), l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

     ! write(*,*)'Start to read ion energy density T-derivative'
     ! Ion energy T-derivative (J/g/eV)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((dEIonOverDT_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

     ! write(*,*)'Start to read electron energy density T-derivative'
     ! Electron energy T-derivative (J/g/eV)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((dEElectronOverDT_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

     !write(*,*)'Start to read ion pressure'
     ! Ion pressure (erg/cm**3)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((pIon_II(l,m), l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

    ! write(*,*)'Start to read electron pressure'
     ! Electron pressure (erg/cm**3)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((pElectron_II(l,m), l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

     !write(*,*)'Start to read ion pressure T-derivative'
     ! Ion pressure T-derivative (erg/cm**3/eV)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((dPIonOverDT_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

     !write(*,*)'Start to read electron pressure T-derivative'
     ! Electron pressure T-derivative (erg/cm**3/eV)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((dPElectronOverDT_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

  endif



  !     ----------------------------
  ! ... now the multigroup opacities
  !     ----------------------------


  RossOpacity_III       = 0.
  PlanckOpacity_III = 0.
  PlanckOpacEms_III= 0.

  ! ... read in the rosseland, planck emission, and planck absorption tables;
  !     values are in (cm**2/g)

  if ( DoReadOpacity ) then



     do id=1, n_dens_opc_tbl
        do it=1, n_temp_opc_tbl
           read (iUnit,802) StringHeader
!           write(*,*)'Read Ross Opac:'//StringHeader
           read (iUnit,*) (RossOpacity_III(ig,it,id),ig=1,num_freq_grps_local)
           read (iUnit,802) StringHeader
!           write(*,*)'Read Planck Opac Emis:'//StringHeader
           read (iUnit,*) (PlanckOpacEms_III(ig,it,id),ig=1,num_freq_grps_local)
           read (iUnit,802) StringHeader
!           write(*,*)'Read Planck Opac Abs'//StringHeader
           read (iUnit,*) (PlanckOpacity_III(ig,it,id),ig=1,num_freq_grps_local)

        enddo
     enddo
!     read (iUnit,802) StringHeader
!     write(*,*)StringHeader
     
  end if
  return

  ! ... format statement to read the header
802 format (a80)

end subroutine read_eosopa_file_main
!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError
  write(*,*)StringError
  stop
end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
end subroutine CON_set_do_test
