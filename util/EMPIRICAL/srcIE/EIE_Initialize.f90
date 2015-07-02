!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine EIE_Initialize(iOutputError)

  use ModErrors
  use ModEIE_Interface
  use ModEIEFiles
  use ModIoUnit, only : UnitTmp_

  use EIE_ModWeimer, only: readcoef96, readcoef01
  use read_data,only: read_potential,read_schatable,read_bndy

  implicit none

  integer, intent(out) :: iOutputError
  character (len=100)  :: inFileName
  integer              :: iError

  integer, parameter  :: South_ = 1
  integer, parameter  :: North_ = 2

  logical :: IsFound_EFieldModel
  integer :: nAmieLats, nAmieMlts, nAmieBlocks

  iError = 0
  iOutputError = 0

  IsFound_EFieldModel = .false.

  call set_error_codes

  !\
  ! --------------------------------------------------------------------
  ! Electric Field Models
  ! --------------------------------------------------------------------
  !/

  if (iDebugLevel > 1) &
       write(*,*) "==> Efield Model : ",EIE_NameOfEFieldModel

  if (iDebugLevel > 1) &
       write(*,*) "==> Model Directory : ",EIE_NameOfModelDir

  if (index(EIE_NameOfEFieldModel,'zero') > 0) then
     IsFound_EFieldModel = .true.
  endif

  if (index(EIE_NameOfEFieldModel,'weimer96') > 0) then
     IsFound_EFieldModel = .true.
     call merge_str(EIE_NameOfModelDir, weimer96_file)
     open(UnitTmp_,file=weimer96_file,status='old', iostat = iError)
     if (iError /= 0) then
        write(6,*) 'Error opening file :',weimer96_file
        iOutputError = ecFileNotFound_
     endif
     call ReadCoef96(UnitTmp_)
     close(UnitTmp_)
  endif

  if (index(EIE_NameOfEFieldModel,'weimer01') > 0) then
     IsFound_EFieldModel = .true.
     call merge_str(EIE_NameOfModelDir, weimer01_file)
     open(UnitTmp_,file=weimer01_file,status='old',&
          form='unformatted', iostat = iError)
     if (iError /= 0) then
        write(6,*) 'Error opening file :',weimer01_file
        iOutputError = ecFileNotFound_
     endif
     call ReadCoef01(UnitTmp_)
     close(UnitTmp_)
  endif

  if (index(EIE_NameOfEFieldModel,'weimer05') > 0) then
     IsFound_EFieldModel = .true.
     inFileName = 'W05scEpot.dat'
     call merge_str(EIE_NameOfModelDir, inFileName)
     call read_potential(inFileName)
     inFileName = 'SCHAtable.dat'
     call merge_str(EIE_NameOfModelDir, inFileName)
     call read_schatable(inFileName)
     inFileName = 'W05scBndy.dat'
     call merge_str(EIE_NameOfModelDir, inFileName)
     call read_bndy(inFileName)
  endif

!  if (index(EIE_NameOfEFieldModel,'samie') > 0) then
!     IsFound_EFieldModel = .true.
!     call merge_str(EIE_NameOfModelDir, stat_amie_file)
!     open(UnitTmp_,file=stat_amie_file,status='old', iostat = iError)
!     if (iError /= 0) then
!        write(6,*) 'Error opening file :',stat_amie_file
!        iOutputError = ecFileNotFound_
!     endif
!     call read_amies(UnitTmp_)
!     close(UnitTmp_)
!  endif

  if (index(EIE_NameOfEFieldModel,'millstone_hpi') > 0) then
     IsFound_EFieldModel = .true.
     call merge_str(EIE_NameOfModelDir, millstone_hill_i_file)
     open(UnitTmp_,file=millstone_hill_i_file,status='old', iostat = iError)
     if (iError /= 0) then
        write(6,*) 'Error opening file :',millstone_hill_i_file
        iOutputError = ecFileNotFound_
     endif
     call mhinit(1, UnitTmp_, 1, iDebugLevel)
     close(UnitTmp_)
  endif

  if (index(EIE_NameOfEFieldModel,'millstone_imf') > 0) then
     IsFound_EFieldModel = .true.
     call merge_str(EIE_NameOfModelDir, millstone_hill_s_file)
     open(UnitTmp_,file=millstone_hill_s_file,status='old', iostat = iError)
     if (iError /= 0) then
        write(6,*) 'Error opening file :',millstone_hill_s_file
        iOutputError = ecFileNotFound_
     endif
     call mhinit(2, UnitTmp_, 1, iDebugLevel)
     close(UnitTmp_)
  endif

  if (index(EIE_NameOfEFieldModel,'hmr89') > 0) then
     IsFound_EFieldModel = .true.
     call merge_str(EIE_NameOfModelDir, hepner_maynard_file)
     open(UnitTmp_,file=hepner_maynard_file,status='old', iostat = iError)
     if (iError /= 0) then
        write(6,*) 'Error opening file :',hepner_maynard_file
        iOutputError = ecFileNotFound_
     endif
     call gethmr(UnitTmp_)
     close(UnitTmp_)
  endif

  if (index(EIE_NameOfEFieldModel,'izmem') > 0) then
     IsFound_EFieldModel = .true.
     call merge_str(EIE_NameOfModelDir, izmem_file)
     open(UnitTmp_,file=izmem_file,status='old', iostat = iError)
     if (iError /= 0) then
        write(6,*) 'Error opening file :',izmem_file
        iOutputError = ecFileNotFound_
     endif
     call izinit(UnitTmp_)
     close(UnitTmp_)
  endif

  !\
  ! --------------------------------------------------------------------
  ! Conductance Models
  ! --------------------------------------------------------------------
  !/

  if (iDebugLevel > 1) &
       write(*,*) "==> Conductance Model : ",EIE_NameOfAuroralModel

  if (index(EIE_NameOfAuroralModel,'ihp') > 0)  &
       call read_conductance_model(iError)
  if (index(EIE_NameOfAuroralModel,'hpi') > 0)  &
       call read_conductance_model(iError)
  if (index(EIE_NameOfAuroralModel,'pem') > 0)  &
       call read_conductance_model(iError)

  if (iDebugLevel > 4) write(*,*) "=====> Back from read conductance"

  if (index(EIE_NameOfEFieldModel,'amie') > 0) then

     IsFound_EFieldModel = .true.
     UseGridBasedEIE = .true.
     UAl_UseGridBasedEIE = .true.

     call AMIE_SetFileName(AMIEFileNorth)
     call readAMIEOutput(North_, .false., iError)

     if (index(AMIEFileSouth,'mirror') > 0) then
        call AMIE_SetFileName(AMIEFileNorth)
        call readAMIEOutput(South_, .true., iError)
     else
        call AMIE_SetFileName(AMIEFileSouth)
        call readAMIEOutput(South_, .false., iError)
     endif

     call AMIE_GetnLats(nAmieLats)
     call AMIE_GetnMLTs(nAmieMlts)
     nAmieBlocks = 2

     call EIE_InitGrid(nAmieLats, nAmieMlts, nAmieBlocks, iOutputError)

     call AMIE_GetLats(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs,&
          EIEr3_HaveLats,iError)

     call AMIE_GetMLTs(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs,&
          EIEr3_HaveMLTs,iError)

     return

  endif

  if (.not.IsFound_EFieldModel) then
     iOutputError = ecEFieldModelNotFound_
  endif

  if (iDebugLevel > 3) write(*,*) "====> Done with EIE_Initialize"

end subroutine EIE_Initialize

!------------------------------------------------------------------------
!------------------------------------------------------------------------

subroutine EIE_InitGrid(nLats, nMlts, nBlocks, iOutputError)

  use ModErrors
  use ModEIE_Interface
  use ModEIEFiles

  implicit none

  integer, intent(in)  :: nLats, nMlts, nBlocks
  integer, intent(out) :: iOutputError

  integer :: iError

  iError = 0

  EIEi_HavenLats = nLats
  EIEi_HavenMLTs = nMlts
  EIEi_HavenBLKs = nBlocks

  if (iDebugLevel > 1) then
     write(*,*) "=> EIEi_HavenBLKs : ", EIEi_HavenBLKs
     write(*,*) "=> EIEi_HavenLats : ", EIEi_HavenLats
     write(*,*) "=> EIEi_HavenMLTs : ", EIEi_HavenMLTs
  endif

  allocate(EIEr3_HaveLats(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs), &
       stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array EIEr3_HaveLats in Interface"
     stop
  endif

  allocate(EIEr3_HaveMlts(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs), &
       stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array EIEr3_HaveMlts in Interface"
     stop
  endif

  allocate(EIEr3_HavePotential(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs), &
       stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array EIEr3_HavePotential in Interface"
     stop
  endif

  allocate(EIEr3_HaveEFlux(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs), &
       stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array EIEr3_HaveEFlux in Interface"
     stop
  endif

  allocate(EIEr3_HaveAveE(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs), &
       stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array EIEr3_HaveAveE in Interface"
     stop
  endif

  iOutputError = iError

end subroutine EIE_InitGrid

!------------------------------------------------------------------------
!------------------------------------------------------------------------

subroutine EIE_FillLats(Lats, iOutputError)

  use ModErrors
  use ModEIE_Interface

  implicit none

  real, intent(in)  :: Lats(EIEi_HavenLats)
  integer, intent(out) :: iOutputError

  integer :: iMlt

  iOutputError = 0

  do iMlt=1,EIEi_HavenMlts
     EIEr3_HaveLats(iMlt,:,1) = Lats
  enddo

end subroutine EIE_FillLats

!------------------------------------------------------------------------
!------------------------------------------------------------------------

subroutine EIE_FillMltsOffset(Mlts, iOutputError)

  use ModErrors
  use ModEIE_Interface

  implicit none

  real, intent(in)  :: Mlts(EIEi_HavenMlts-1)
  integer, intent(out) :: iOutputError

  integer :: iMlt,iLat

  iOutputError = 0

  do iLat=1,EIEi_HavenLats
     EIEr3_HaveMlts(1,iLat,1) = Mlts(EIEi_HavenMlts-2)-360.0
     EIEr3_HaveMlts(2:EIEi_HavenMlts,iLat,1) = Mlts
  enddo

end subroutine EIE_FillMltsOffset

