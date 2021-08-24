! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine init_mpi_coup

  use ModGITM
  use ModTime
  use ModMpi
  use ModInputs, only: cInputFile, iInputUnit_, iCharLen_
  use ModCoupSAMI3, only: DtCouple, iProcGlobal

  !use parameter_mod,only: numwork
  implicit none
  
  integer :: global_grp,local_grp1,local_grp2

  integer, allocatable:: ranks1(:),ranks2(:)

  integer :: iError
  integer  :: i

  logical :: IsDone = .false., IsThere
  logical :: HaveGrid = .false., HaveSami = .false.

  character (len=iCharLen_) :: cLine
  integer :: nBlksLons, nBlksLats

  numsami = 13
  numgitm = 16

  ! setup mpi
  call MPI_INIT(iError)

  iCommGlobal = MPI_COMM_WORLD
   
  ! create global group
  call mpi_comm_group(iCommGlobal, global_grp, iError)
  call MPI_COMM_RANK(iCommGlobal, iProcGlobal, iError)

  ! ---------------------------------------------------------------------------
  ! This is a total hack, and I am ashamed.

  inquire(file=cInputFile,EXIST=IsThere)
  if (.not.IsThere) &
       call stop_gitm(cInputFile//" cannot be found by read_inputs")

  open(iInputUnit_,file=cInputFile,status="old")

  do while (.not.IsDone)

     read(iInputUnit_,'(a)',iostat=iError) cLine

     if (iError /= 0) then
        IsDone = .true.
     else

        if (cLine(1:5) == "#GRID") then
           read(iInputUnit_,*) nBlksLons
           read(iInputUnit_,*) nBlksLats
           numgitm = nBlksLons*nBlksLats
           HaveGrid = .true.
        endif

        if (cLine(1:5) == "#SAMI") then
           read(iInputUnit_,*) numsami
           read(iInputUnit_,*) DtCouple
           HaveSami = .true.
        endif

        if (HaveGrid.and.HaveSami) IsDone = .true.

     endif

  enddo

  close(iInputUnit_)

  if (.not.HaveGrid) then
     if (iProcGlobal == 0) then
        write(*,*) "Can't seem to find #GRID in UAM.in file...??? How???"
     endif
     stop
  endif

  if (.not.HaveSami) then
     if (iProcGlobal == 0) then
        write(*,*) "Need to add lines in UAM.in:"
        write(*,*) "#SAMI"
        write(*,*) "13        numworker + 1 from SAMI - We need to fix this.."
        write(*,*) "300.0     dtcouple"
     endif
     stop
  endif

  ! ---------------------------------------------------------------------------

  allocate(ranks1(numgitm))
  allocate(ranks2(numsami))

  ranks1=(/(i,i = 0,numgitm-1)/)
  ranks2=(/(i,i = numgitm,numgitm+numsami-1)/)

  SamiMaster = numgitm
   
  ! extract 4 ranks from global group to make a local group
  call mpi_group_incl(global_grp, numgitm, ranks1, local_grp1, iError)
  call mpi_group_incl(global_grp, numsami, ranks2, local_grp2, iError)

  ! make new communicator based on local group
  call mpi_comm_create(iCommGlobal, local_grp1, iCommGITM, iError)
  call mpi_comm_create(iCommGlobal, local_grp2, iCommSAMI0, iError)

end subroutine init_mpi_coup

subroutine init_mpi_gitm

  use ModGITM
  use ModTime
  use ModMpi

  implicit none

  integer :: iError

  ! setup mpi
  !call MPI_INIT(iError)
  !iCommGITM = MPI_COMM_WORLD
  call MPI_COMM_RANK(iCommGITM, iProc, iError)
  call MPI_COMM_SIZE(iCommGITM, nProcs, iError)

end subroutine init_mpi_gitm

