! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine init_mpi

  use ModGITM
  use ModTime
  use ModMpi

  implicit none

  integer :: iError

  ! setup mpi
  call MPI_INIT(iError)
  iCommGITM = MPI_COMM_WORLD
  call MPI_COMM_RANK(iCommGITM, iProc, iError)
  call MPI_COMM_SIZE(iCommGITM, nProcs, iError)

end subroutine init_mpi




