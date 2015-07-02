!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModReadArtep
  use ModIoUnit, ONLY: io_unit_new
  implicit none

  integer:: nRho = 201, nTe = 201, nGroup = 30

  SAVE
  real, allocatable:: Value_VII(:,:,:) 
  real, allocatable:: EGroup_I(:)
  integer::line
  real:: TMin, TMax, RhoMin, RhoMax
  character(LEN=15)NameEmpty
contains
  !============================================================================
  subroutine read_artep
    integer:: iFile, iRho, iTe
    real:: Aux_V(7), Aux
    !-----------------------------------------
    iFile = io_unit_new()
    open(iFile, file='header.x4',status = 'old')
    do line = 1,13
       read(iFile,*)
    end do
    read(iFile,*)NameEmpty, TMin, TMax, nTe
    write(*,*)'NameEmpty, TMin, TMax, nTe=',NameEmpty, TMin, TMax, nTe
    read(iFile,*)NameEmpty, RhoMin, RhoMax, nRho
    write(*,*)'NameEmpty, RhoMin, RhoMax, nRho=', &
         NameEmpty, RhoMin, RhoMax, nRho
    read(iFile,*)
    read(iFile,*)
    read(iFile,*)
    read(iFile,*)NameEmpty,NameEmpty,NameEmpty,nGroup
    write(*,*)'nGroup=', nGroup
    allocate(Value_VII(2*nGroup,nRho,nTe), EGroup_I(1+nGroup))
   
    read(iFile,*)EGroup_I
    close(iFile)
    write(*,*)EGroup_I

    iFile = io_unit_new()
    open(iFile, file='x4_opac.artep.rev', form='formatted',status='old')
    line=0
    do iTe = 1, nTe
       do iRho = 1, nRho
          line = line+1;write(*,*)line
          read(iFile,*) Aux_V, Value_VII(1:nGroup,iRho,iTe), &
               Aux,            Value_VII(nGroup+1:2*nGroup,iRho,iTe)
       end do
    end do
    close(iFile)

  end subroutine read_artep

end Module ModReadArtep
!==============================================================================
program save_eos_table

  use CRASH_ModStatSum
  use CRASH_ModPartition, ONLY: UseCoulombCorrection
  use CRASH_ModPolyimide
  use CRASH_ModFermiGas
  use CRASH_ModMultiGroup
  use CRASH_ModExcitation
  use CRASH_ModAtomicDataMix
  use CRASH_ModExcitationData,ONLY : n_ground, cExcitationN_III
  use CRASH_ModIonMix
  use ModReadArtep
  use ModPlotFile
  use ModMpi
  use CRASH_ModAtomicMass,ONLY: cAtomicMass_I
  use ModConst
  use ModIoUnit
  implicit none

  integer:: iError, iMaterial, iFile
  real:: AtomicMass
  !----------------------------------------------------------------------------------

  call read_artep
  Value_VII = 0.10 * Value_VII  !cm^2/g = 1e-4 m^2/(1e-3 kg) = 0.1 m2/kg
  RhoMin = 1.0e6 * RhoMin * cAtomicMass * cAtomicMass_I(54) !1/cm3=>kg/m3
  RhoMax = 1.0e6 * RhoMax * cAtomicMass * cAtomicMass_I(54) !1/cm3=>kg/m3
  call save_plot_file( &
         'Xe_opac.dat',                                                      &
         TypeFileIn     = 'real8',                                           &
         StringHeaderIn = 'ARTEP Opacity for Xe',                            &
         NameVarIn      = 'logRho logTe Planck(44) Ross(44) EGropu00 EGroup(44)',&
         CoordMinIn_D   = (/log10(RhoMin), log10(TMin)/),        &
         CoordMaxIn_D   = (/log10(RhoMax), log10(TMax)/),        &
         ParamIn_I      = EGroup_I,                                     &
         VarIn_VII      = Value_VII)
  iFile = io_unit_new()
  open(iFile,file='Xe_opac_ARTEP44.head')
       write(iFile,'(a)')'  '
       write(iFile,'(a)')'----------------OPACITY TABLE----------'
       write(iFile,'(a)')'  '
       write(iFile,'(a)')'#LOOKUPTABLE'
       write(iFile,'(a)')'Xe_opac                 NameTable'
       write(iFile,'(a)')'use param               NameCommand'
       write(iFile,'(a)')'Tables/Xe_opac_ARTEP44.dat'
       write(iFile,'(a)')'real8                   TypeFile'
       write(iFile,'(a)')'EGroup00 EGroup(44)'
       write(iFile,'(e13.7)')EGroup_I
       write(iFile,'(a)')'Opacity(rho,Te) for Xe'
       write(iFile,'(a)')'logRho LogT Planck(44) Ross(44)'
       write(iFile,'(i4,a)')nRho,'                    nIndex1'
       write(iFile,'(e13.7,a)')RhoMin,&
            '           Index1Min (kg/m3)'
       write(iFile,'(e13.7,a)')RhoMax,& 
            '           Index1Max (kg/m3)'
       write(iFile,'(i4,a)')nTe,'                    nIndex2'
       write(iFile,'(e13.7,a)')TMin,&
            '           Index2Min (eV)'
       write(iFile,'(e13.7,a)')TMax,&
            '           Index2Max (eV)'
       write(iFile,'(a)')'  '
       write(iFile,'(a)')'#END'
       write(iFile,'(a)')'  '
       close(iFile)

end program save_eos_table
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
