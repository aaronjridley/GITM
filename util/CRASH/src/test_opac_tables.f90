!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

program save_opac_table
  use CRASH_ModStatSum
  use CRASH_ModPartition, ONLY: UseCoulombCorrection
  use CRASH_ModPolyimide
  use CRASH_ModFermiGas
  use CRASH_ModMultiGroup
  use CRASH_ModExcitation
  use CRASH_ModAtomicDataMix
  use CRASH_ModExcitationData,ONLY : n_ground, cExcitationN_III
  use CRASH_ModIonMix
  use CRASH_ModEos, ONLY:UseEosTable_I, UseOpacityTable_I, nMaterialEos
  use CRASH_ModEosTable, ONLY:&
                         check_opac_table,IndexDefaultOpac_I
  use ModMpi
  implicit none
  integer:: iError, iMaterial
  !----------------
  UseCoulombCorrection = .false.
  UseExcitation = .true.

  UseEosTable_I = .false.
  UseOpacityTable_I = .true.

  !Test: small grid, ascii output
  !\
  !Set reduced table sizes
  !/
  !  IndexDefaultOpac_I = (/9,9/)
  !  nGroup = 10
  !  nMaterialEos = 1
  !  call check_opac_table(iComm=MPI_COMM_WORLD,EGroupIn_I=(/0.10,20000.0/),&
  !       TypeFileIn='ascii')

  !\
  ! Production run
  !/
  !
  !\
  ! As in CRASH_data/LookupTable/README
  !/

  DoStateElimination = .true.

  !Five materials, 30 groups regular grid

  nMaterialEos = 5

  !nGroup = 30
  !call check_opac_table(iComm=MPI_COMM_WORLD,EGroupIn_I=(/0.10,20000.0/))

  !\
  ! Process the input file with the energy grid
  !/
  nGroup = -1
  open(11, file='EGroup30',status='old')
  do
     read(11,*,err=100,end=100) EnergyGroup_I(nGroup+1)
     nGroup = nGroup +1
  end do
100 continue
  
  call check_opac_table(iComm=MPI_COMM_WORLD,EGroupIn_I=EnergyGroup_I(0:nGroup))
  !  ,TypeFileIn='ascii')

  call MPI_Finalize(iError)
end program save_opac_table
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
