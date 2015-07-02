!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program abs
  use CRASH_ModStatSum
  use CRASH_ModPartition
  use CRASH_ModPolyimide
  use CRASH_ModFermiGas
  use CRASH_ModMultiGroup,ONLY : set_multigroup, EnergyGroup_I
  use CRASH_ModMultiGroup,ONLY :  UseBremsstrahlung, UsePhotoionization ,&
        UseScattering
  use CRASH_ModExcitation
  use CRASH_ModAtomicDataMix
  use CRASH_ModEos
  use ModIoUnit, ONLY: io_unit_new
  use ModConst
  implicit NONE


  !Table dimensions
  integer,parameter :: nTemperature=20
  real,parameter    :: dTe = 10 ![eV]
  integer,parameter :: nDensity = 5
  real,parameter    :: DensityStart =1.0e19   ![cm3]
  real,parameter    :: DensityRatio = 10.0
  integer, parameter:: nGroup = 100
  !Description of CRASH meterials:

  integer :: iGroup,iMaterial,iT, iNa
  real :: NaTable,Rho,T
  real::OpacityPlanck_I(nGroup),OpacityRosseland_I(nGroup)
  character(LEN=24)::NameFile 
  
  !---------------

  call set_multigroup(nGroup, 1.0/cHPlanckEV, 10000.0/cHPlanckEV)
  UseExcitation = .true.
  UsePreviousTe = .false.

  UseBremsstrahlung = .true.
  UsePhotoionization = .true.
 
  UseCoulombCorrection = .true.
  UseScattering = .true.
  DoStateElimination = .true.



  NaTable = DensityStart /DensityRatio
  do iNa = 1, nDensity
     NaTable = NaTable * DensityRatio
     do iMaterial = Xe_,Plastic_
        !Get Mass Density
        Rho = NaTable * 1.0e6 * &! [1/m3]
             cAtomicMassCRASH_I(iMaterial) * cAtomicMass
        write(*,*)'Density of '//NameMaterial_I(iMaterial)//' is ',Rho,' kg/m3'
        do iT = 1, nTemperature
           T= iT * dTe        ![eV]

           call eos(iMaterial = iMaterial, &!
                    Rho  = Rho           , &! [kg/m3]
                    TeIn = T * cEVToK    , &! Convert to K
                    OpacityPlanckOut_I=    &
                    OpacityPlanck_I      , &! [1/m]
                    OpacityRosselandOut_I= &
                    OpacityRosseland_I     )! [1/m]
           
           write(NameFile,'(a,e8.2,a,e8.2,a)')NameMaterial_I(iMaterial)//'_',NaTable,'_',T,'.dat'
           call save_opacity_sesame
        end do
     end do
  end do
  
contains
  subroutine save_opacity_sesame
    integer:: unit
    !------------------------------------
    unit = io_unit_new()
    open(unit,file=NameFile)
    write(unit,'(a)')' Opacity data for '//NameMaterial_I(iMaterial)
    write(unit,'(a,e8.2,a,e8.2,a,e8.2,a)')&
         ' atomic concentration = ', NaTable,' [1/cm3], mass density = ', Rho,&
         ' [kg/m3], temperature = ', T, ' [eV].' 
    write(unit,'(a)') ' Energy group [keV]  Rosseland [cm2/g]  Planck [cm2/g]' 
    do iGroup = 1, nGroup 
       write(unit,'(3e12.2)')1.0e-3 * EnergyGroup_I(iGroup-1),&! [keV]
            OpacityRosseland_I(iGroup)*10.0/Rho,             &! [cm2/g]
            OpacityPlanck_I(iGroup)*10.0/Rho                  ! [cm2/g]
    end do
    close(unit)
  end subroutine save_opacity_sesame
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

