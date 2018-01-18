!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP
!MODULE: ModMpiOrig and ModMpi - the MPI variables and functions
!DESCRIPTION:
! In Fortran 90 it is customary to use a module instead of including files.
! The ModMpiOrig and ModMpi modules provide interfaces to mpif.h.
!
! In ModMpi the MPI\_REAL parameter is set to the value of the 
! MPI\_DOUBLE\_PRECISION parameter if the code is compiled with double
! precision accuracy (Config.pl -double). In this case the iRealPrec
! parameter is set to 1.
!
! If ModMpi is used, it should be compiled with the same precision as the F90
! code using it.
!
! This module also provides simple interfaces to 
! MPI\_reduce with MPI\_IN\_PLACE option for real and integer 
! scalars and arrays.
!
!REVISION HISTORY:
! 07/02/2003 G. Toth <gtoth@umich.edu> - initial version of ModMpi
! 07/20/2003 G. Toth - change the MPI_REAL definition
! 07/30/2004 G. Toth - updated the description for the modified mpif90.h files.
! 05/13/2011 G. Toth - modified to use original mpif.h file
! 03/10/2016 G. Toth - added MPI_reduce_* routines for in-place use.
!INTERFACE:
!BOP
!INTERFACE:
module ModMpi
  !EOP
  !BOC
  use ModMpiInterfaces

  use ModMpiOrig, &
       MPI_REAL_ORIG => MPI_REAL, MPI_COMPLEX_ORIG => MPI_COMPLEX

  implicit none

  ! iRealPrec = 1 if the code is compiled with 8 byte reals and 0 otherwise
  integer, parameter :: iRealPrec = (1.00000000011 - 1.0)*10000000000.0

  integer, parameter :: MPI_REAL = &
       iRealPrec*MPI_DOUBLE_PRECISION + (1-iRealPrec)*MPI_REAL_ORIG

  integer, parameter :: MPI_COMPLEX = &
       iRealPrec*MPI_DOUBLE_COMPLEX + (1-iRealPrec)*MPI_COMPLEX_ORIG

  !EOC

contains
  !============================================================================
  subroutine mpi_reduce_real_array(Buffer_I, nSize, iOp, iRoot, iComm, iError)

    real, intent(inout):: Buffer_I(*)
    integer, intent(in):: nSize
    integer, intent(in):: iOp
    integer, intent(in):: iRoot
    integer, intent(in):: iComm
    integer, intent(out):: iError

    external MPI_reduce
    integer:: iRank
    real:: Recv_I(1)
    !------------------------------------------------------------------------
    call MPI_comm_rank(iComm, iRank, iError)

    if(iRoot == iRank)then
       call MPI_reduce(MPI_IN_PLACE, Buffer_I, nSize, MPI_REAL, iOp, &
            iRoot, iComm, iError)
    else
       call MPI_reduce(Buffer_I, Recv_I, nSize, MPI_REAL, iOp, &
            iRoot, iComm, iError)
    end if

  end subroutine mpi_reduce_real_array
  !============================================================================
  subroutine mpi_reduce_real_scalar(Value, iOp, iRoot, iComm, iError)

    real, intent(inout):: Value
    integer, intent(in):: iOp
    integer, intent(in):: iRoot
    integer, intent(in):: iComm
    integer, intent(out):: iError

    external MPI_reduce
    integer:: iRank
    real :: Recv
    !------------------------------------------------------------------------
    call MPI_comm_rank(iComm, iRank, iError)

    if(iRoot == iRank)then
       call MPI_reduce(MPI_IN_PLACE, Value, 1, MPI_REAL, iOp, &
            iRoot, iComm, iError)
    else
       call MPI_reduce(Value, Recv, 1, MPI_REAL, iOp, &
            iRoot, iComm, iError)
    end if
    
  end subroutine mpi_reduce_real_scalar
  !============================================================================
  subroutine mpi_reduce_integer_array(&
       iBuffer_I, nSize, iOp, iRoot, iComm, iError)

    integer, intent(inout):: iBuffer_I(*)
    integer, intent(in):: nSize
    integer, intent(in):: iOp
    integer, intent(in):: iRoot
    integer, intent(in):: iComm
    integer, intent(out):: iError

    external MPI_reduce
    integer:: iRank, iRecv_I(1)
    !------------------------------------------------------------------------
    call MPI_comm_rank(iComm, iRank, iError)

    if(iRoot == iRank)then
       call MPI_reduce(MPI_IN_PLACE, iBuffer_I, nSize, MPI_INTEGER, &
            iOp, iRoot, iComm, iError)
    else
       call MPI_reduce(iBuffer_I, iRecv_I, nSize, MPI_INTEGER, &
            iOp, iRoot, iComm, iError)
    end if

  end subroutine mpi_reduce_integer_array
  !============================================================================
  subroutine mpi_reduce_integer_scalar(iValue, iOp, iRoot, iComm, iError)

    integer, intent(inout):: iValue
    integer, intent(in):: iOp
    integer, intent(in):: iRoot
    integer, intent(in):: iComm
    integer, intent(out):: iError

    external MPI_reduce
    integer:: iRank, iRecv
    !------------------------------------------------------------------------
    call MPI_comm_rank(iComm, iRank, iError)

    if(iRoot == iRank)then
       call MPI_reduce(MPI_IN_PLACE, iValue, 1, MPI_INTEGER, &
            iOp, iRoot, iComm, iError)
    else
       call MPI_reduce(iValue, iRecv, 1, MPI_INTEGER, &
            iOp, iRoot, iComm, iError)
    end if

  end subroutine mpi_reduce_integer_scalar

end module ModMpi
