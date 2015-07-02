!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program abs
  use CRASH_ModStatSum
  use CRASH_ModPartition
  use CRASH_ModPolyimide
  use CRASH_ModFermiGas
  use CRASH_ModMultiGroup
  use CRASH_ModExcitation
  use CRASH_ModAtomicDataMix
  use CRASH_ModExcitationData,ONLY : n_ground, cExcitationN_III
  use CRASH_ModIonMix
  use ModConst
  implicit NONE

  integer, parameter :: unit = 24

  real:: vTe = 10.0 !eV
  real::NaTrial = 1.0e22
  integer:: iPlot,iError
  integer :: iL, iN, iZ, iMix, i
  integer :: nGround, iGroup
  real :: TSI, Eg_W(100), CFromT_W(100),cFromE,TOutSI
  real,dimension(1:nZMax) :: IonizPotential_I
  real::RhoAl = 0.1 ![g/cm3]
  !---------------

  open(unit,file='../doc/excited_levels.tex',status='replace')
  write(unit,'(a)')'\newcolumntype{x}[1]{>{\centering\hspace{0pt}}p{#1}}'
  write(unit,'(a)')'\begin{tabular}'//&
       '{|x{1cm}|x{1cm}||x{2cm}||x{2cm}|x{2cm}|x{2cm}|x{2cm}|x{2cm}|}'
  write(unit,'(a)')'\hline'
  write(unit,'(a)')'\multirow{2}{*}{i} & \multirow{2}{*}{n} & '//&
       '\multirow{2}{*}{formula} & '//&
       '\multicolumn{5}{c|}{database, for different values of $l$} \tabularnewline'
  write(unit,'(a)')'\cline{4-8}'
  write(unit,'(a)')' & & & s & p & d & f & g \tabularnewline'
  write(unit,'(a)')'\hline'
  write(unit,'(a)')'\hline'

  call get_ioniz_potential(7, IonizPotential_I(1:7))

  do iZ = 0, 6
     write(unit,'(a,i3,a)') '\multirow{4}{*}{', iZ, '} '

     do iN = 2, 5

        nGround = n_ground(iZ, 7)

        write(unit,'(a,i3,a,f6.1,5(a,f6.1),a)') ' & ', iN, ' & ', &
             max(IonizPotential_I(iZ+1) - &
             cRyToEV * (real(iZ+1)/iN)**2, 0.0), &
             (' & ', cExcitationN_III(iL, iN, iZ), iL = 0, 4), &
             ' \tabularnewline'

     end do

     write(unit,'(a)')'\hline'

  end do

  write(unit,'(a)')'\end{tabular}'
  close(unit)

  open(unit,file='../doc/normalized_planckian.dat',status='replace')
  write(unit,*)'Notmalized integral Planckian, \int_0^x{Planckian(y)dy}'
  do i=1,300
     write(unit,*)0.1*i,gint5(0.1*i)*cNormG5
  end do
  close(unit)
  open(unit,file='../doc/normalized_planckian_derivative.dat',status='replace')
  write(unit,*)'Notmalized integral of the Planckian derivative, \int_0^x{(d Planckian(y)/dT)dy}'
  do i=1,300
     write(unit,*)0.1*i,gint6(0.1*i)*cNormG6
  end do
  close(unit)

  UseExcitation = .true.
  call set_mixture(nPolyimide, nZPolyimide_I, CPolyimide_I)
  UsePreviousTe = .false.

  DoNotAddLineCore = .false.
  UseBremsstrahlung = .true.
  UsePhotoionization = .true.
 
  UseCoulombCorrection = .true.
  call set_ionization_equilibrium(vTe,NaTrial*1000000.0,iError)
  do iZ = 0, 0
     write(*,'(a,i3,a)') '\multirow{4}{*}{', iZ, '} '

     do iN = 1, 10

        nGround = n_ground(iZ, 1)

        write(*,'(a,i3,a,f6.1,5(a,f6.1),a)') ' & ', iN, ' & ', &
             max(IonizPotential_II(iZ+1,2) - &
             cRyToEV * (real(iZ+1)/iN)**2, 0.0), &
             ' & ', ExcitationEnergy_III(iN, iZ,2), &
             ' \tabularnewline'

     end do

     write(*,'(a)')'\hline'

  end do
  call save_report('report.txt')

  call set_multigroup(100, 1.0/cHPlanckEV, 10000.0/cHPlanckEV)
  call meshhv
  call abscon
  
  open(unit,file='../doc/polyimide_absorption.dat')
  write(unit,'(a,i6,a)') &
       'Photon energy [eV]  Absorbtion Coeff cm-1, in ',nPhoton,' points' 
  do iPlot = 1, nPhoton
     if(PhotonEnergy_I(iPlot)< 0.1*Te.or.PhotonEnergy_I(iPlot)>1000.0*Te)&
          CYCLE
     write(unit,*)log10(PhotonEnergy_I(iPlot)),&
          log10(AbsorptionCoefficient_I(iPlot))
  end do
  close(unit)
  !==================  Xeonon test===============
  DoStateElimination = .true.
  call set_element(54)
  UsePreviousTe = .false.

  DoNotAddLineCore = .false.
  UseBremsstrahlung = .true.
  UsePhotoionization = .true.
 
  UseCoulombCorrection = .true.
  call set_ionization_equilibrium(vTe,3.0e25,iError)
  call save_report('report_xe.txt')
 
  call set_multigroup(100, 1.0/cHPlanckEV, 10000.0/cHPlanckEV)
  UseScattering = .true.
  DoNotAddLineCore =.false.
  call meshhv
  call abscon
  
  call save_absorption_coef('../doc/xenon_absorption.dat')

  DoNotAddLineCore =.true.
  call meshhv
  call abscon
  
  call opacys(TRadIn = Te)
  open(unit,file='../doc/xenon_opacities.dat')
  do iGroup = 1, nGroup 
     write(unit,*)0.5*log10(EnergyGroup_I(iGroup))&
     +0.5*log10(EnergyGroup_I(iGroup-1)),&
          log10(OpacityPlanck_I(iGroup)),&
          log10(OpacityRosseland_I(iGroup))
  end do
  close(unit)
 
  write(*,*)'OpacityPlanckTotal=',OpacityPlanckTotal,&
      '    OpacityRosselandTotal=',OpacityRosselandTotal, '   Z=',ZAv
  TSI = 1.0e5   !\approx 9 eV
  open(unit,file='conversion_report.txt') 
  write(unit,*)'TSI=', TSI, '  Total Radiation=',cRadiation * TSI**4, '  Total specific heat=', 4*cRadiation * TSI**3
  write(unit,*)'log10(hv) [eV]      log10(Eg) [J/m^3]  Tg  CFromT CFromE'
  do iGroup = 1,nGroup

     call get_energy_g_from_temperature(iGroup,TSI,Eg_W(iGroup),CFromT_W(iGroup))
     call get_temperature_from_energy_g(iGroup, Eg_W(iGroup), TOutSI,CFromE)
     write(unit,'(5e10.3)') sqrt(EnergyGroup_I(iGroup - 1)*&
          EnergyGroup_I(iGroup)),Eg_w(iGroup),TOutSI,CFromT_W(iGroup),CFromE
  end do
  write(unit,*)'Energy sum up:',sum(Eg_W(1:nGroup)), ' Total specific heat:=', &
       sum(CFromT_W(1:nGroup))
  close(unit)
  !====================ALUMINIUM TEST====================================!

  call set_element(13)

  NaTrial = RhoAl/(1.0e3*cAtomicMass * cAtomicMass_I(13))
  DoStateElimination = .true.
  
  vTe = 100.0 ![eV]=0.1 keV
  UseCoulombCorrection = .false.
  call set_ionization_equilibrium(vTe,NaTrial*1000000.0,iError)
  
  call save_report('report_al.txt')
  DoNotAddLineCore =.false.
  DoEnhancePhotoIonization = .true.
  
  call meshhv
  call abscon
  call save_absorption_coef('../doc/al_absoprtion.dat')
  call opacys(TRadIn = Te)
  call save_opacity_sesame('../doc/al_opacities.dat')
 

  !====================BE TEST====================================!

  call set_element(4)

  NaTrial = RhoAl/(1.0e3*cAtomicMass * cAtomicMass_I(4))
  
  vTe = 100.0 ![eV]=0.1 keV
  UseCoulombCorrection = .true.
  UseExcitation = .true.
  UseGroundStatWeight = .false.
  BoundFreeCorrection = 1.0
  call set_ionization_equilibrium(vTe,NaTrial*1000000.0,iError)
 
  call save_report('report_be.txt')
 
  DoNotAddLineCore =.false.
  DoEnhancePhotoIonization = .true.
  
  call meshhv
  call abscon
  call save_absorption_coef('../doc/Be_absoprtion.dat')
  call opacys(TRadIn = Te)
  call save_opacity_sesame('../doc/Be_opacities.dat')
  !=================ZERO temperature====================!
  UseExcitation = .true.
  call set_mixture(nPolyimide, nZPolyimide_I, CPolyimide_I)
  UsePreviousTe = .false.

  DoNotAddLineCore = .false.
  UseBremsstrahlung = .true.
  UsePhotoionization = .true.
 
  UseCoulombCorrection = .true.
  call set_ionization_equilibrium(0.1,NaTrial*1000000.0,iError)
  call save_report('report_zero.txt')
  do iMix = 1, 0!nMix
     write(*,*)'iMix,nZ)I(iMix)',iMix,nZ_I(iMix)
  do iZ = iZMin_I(iMix),iZMax_I(iMix) 
     write(*,'(a,i3,a)') '\multirow{4}{*}{', iZ, '} '

     do iN = 1, 10

        nGround = n_ground(iZ, nZ_I(iMix))

        write(*,'(a,i3,a,f6.1,5(a,f6.1),a)') ' & ', iN, ' & ', &
             max(IonizPotential_II(iZ+1,iMix) - &
             cRyToEV * (real(iZ+1)/iN)**2, 0.0), &
             ' & ', ExcitationEnergy_III(iN, iZ,iMix), &
             ' \tabularnewline'

     end do

     write(*,'(a)')'\hline'

  end do
end do
  call set_multigroup(100, 1.0/cHPlanckEV, 10000.0/cHPlanckEV)
  call meshhv
  call abscon
  
  open(unit,file='../doc/polyimide_absorption_zero.dat')
  write(unit,'(a,i6,a)') &
       'Photon energy [eV]  Absorbtion Coeff cm-1, in ',nPhoton,' points' 
  do iPlot = 1, nPhoton
     if(PhotonEnergy_I(iPlot)< 0.1*Te.or.PhotonEnergy_I(iPlot)>1000.0*Te)&
          CYCLE
     write(unit,*)log10(max(PhotonEnergy_I(iPlot),1e-10)),&
          log10(max(AbsorptionCoefficient_I(iPlot),1e-10))
  end do
 
contains
  !========================
  subroutine save_report(NameFile)
    character(LEN=*),intent(in)::NameFile
    !------------------------------------
    open(unit,file=NameFile)
    do iMix = 1,nMix
       write(unit,*)'iMix,nZ_I(iMix),iZMin_I(iMix), iZMax_I(iMix)',&
            iMix,nZ_I(iMix),iZMin_I(iMix), iZMax_I(iMix)
       do iZ = iZMin_I(iMix), iZMax_I(iMix)
          write(unit,*)' iZ = ', iZ
          write(unit,*)Partition_III(:,iZ,iMix)*Population_II(iZ,iMix)
       end do
    end do
    close(unit)
  end subroutine save_report
  !=========================
  subroutine save_opacity_sesame(NameFile)
    character(LEN=*),intent(in)::NameFile
    !------------------------------------
    write(*,*)'OpacityPlanckTotal=',OpacityPlanckTotal/RhoAl,&
         '    OpacityRosselandTotal=',OpacityRosselandTotal/RhoAl, '   Z=',ZAv
    open(unit,file=NameFile)
    do iGroup = 1, nGroup 
       write(unit,*)1.0e-3 * EnergyGroup_I(iGroup-1),&
            OpacityRosseland_I(iGroup)/RhoAl,&
            OpacityPlanck_I(iGroup)/RhoAl
    end do
    close(unit)
  end subroutine save_opacity_sesame
  !=========================
  subroutine  save_absorption_coef(NameFile)
    character(LEN=*),intent(in)::NameFile
    !------------------------------------
    open(unit,file=NameFile)
    write(unit,'(a,i6,a)') &
         'Photon energy [eV]  Absorbtion Coeff cm-1, in ',nPhoton,' points' 
    do iPlot = 1, nPhoton
       if(PhotonEnergy_I(iPlot)< 0.1*Te.or.PhotonEnergy_I(iPlot)>1000.0*Te)&
            CYCLE
       write(unit,*)log10(PhotonEnergy_I(iPlot)),&
            log10(AbsorptionCoefficient_I(iPlot)),log10(ScatteringCoefficient_I(iPlot))
    end do
    close(unit)
  end subroutine save_absorption_coef
end program abs
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

