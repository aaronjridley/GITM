!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!!QUOTE: \clearpage
!
!BOP
!
!QUOTE: \section{util/NOMPI: Library for Sequential Execution}
!
!MODULE: NOMPI - a library replacing MPI for sequential execution
!
!DESCRIPTION:
!
! The NOMPI library can be used for debugging on a single processor,
! and possibly running small problems on a machine with no MPI library.
!
! The NOMPI library contains all the MPI subroutines and functions used 
! in SWMF and its parallel components. This library can be used instead 
! of the real MPI library for sequential execution. You have to compile it
! with 'make NOMPI' select MPILIB = ... -lNOMPI in Makefile.conf
! and relink the executable.
!
! Most of the subroutines and functions in this library do nothing,
! or almost nothing, or stop, as it is appropriate, and they are
! needed by the compiler.
!
!REVISION HISTORY:
! 05/11/2001 G. Toth <gtoth@umich.edu> - initial version of NOMPI
!
!EOP
!=============================================================================
! core subroutines and functions used by the real MPI subroutines/functions
!=============================================================================
!BOP
!
!ROUTINE: MPI_TYPE_SIZE - return size of datatype in bytes
!
!INTERFACE:
!
integer function MPI_TYPE_SIZE(datatype)

  !USES:
  use ModMpiOrig, ONLY: MPI_REAL, MPI_DOUBLE_PRECISION, &
       MPI_INTEGER, MPI_LOGICAL, MPI_CHARACTER
  implicit none

  !INPUT ARGUMENTS:
  integer, intent(in) :: datatype
  !EOP

  integer, save :: byte_c=0, byte_l, byte_i, byte_r, byte_d
  !---------------------------------------------------------------------------

  if(byte_c==0)then
     inquire(iolength=byte_c)' '
     inquire(iolength=byte_l).true.
     inquire(iolength=byte_i)1
     inquire(iolength=byte_r)1.0
     inquire(iolength=byte_d)1.0D0
     if(byte_c/=1)then
        write(*,*)'Strange byte size for character variables in MPI_TYPE_SIZE'
        write(*,*)'byte_c, byte_l, byte_i, byte_r, byte_d=',&
             byte_c, byte_l, byte_i, byte_r, byte_d
        write(*,*)'Normalizing with byte_c'
        byte_l=byte_l/byte_c
        byte_i=byte_i/byte_c
        byte_r=byte_r/byte_c
        byte_d=byte_d/byte_c
     end if
  end if

  ! Some datatype values may coincide, so we use IF ELSEIF instead of CASE

  if(datatype==MPI_LOGICAL)then
     MPI_TYPE_SIZE=byte_l
  elseif(datatype==MPI_CHARACTER)then
     MPI_TYPE_SIZE=1
  elseif(datatype==MPI_INTEGER)then
     MPI_TYPE_SIZE=byte_i
  elseif(datatype==MPI_REAL)then
     MPI_TYPE_SIZE=byte_r
  elseif(datatype==MPI_DOUBLE_PRECISION)then
     MPI_TYPE_SIZE=byte_d
  else
     write(*,*)'Error in MPI_TYPE_SIZE: unknown datatype=',datatype
     MPI_TYPE_SIZE=-1
  end if

end function MPI_TYPE_SIZE

!BOP =========================================================================
!ROUTINE: MPI_SIMPLE_COPY - copy data between two buffers
!INTERFACE:
subroutine MPI_SIMPLE_COPY(caller,sendbuf,sendcount,sendtype,&
                                  recvbuf,recvcount,recvtype)

  !DESCRIPTION:
  ! This subroutine is used by all the collective communication routines.
  ! Although the type of sendbuf and recvbuf are character (LEN=*),
  ! they are cast into character from arbitrary data types described by
  ! sendtype and recvtype. The number of bytes copied from sendbuf to
  ! recvbuf is determined by sendcount and sendtype. The number of bytes
  ! corresponding to a given data type is calculated by the function
  ! MPI\_TYPE\_SIZE. recvcount and recvtype are used to check if the recieve
  ! buffer is large enough.

  implicit none

  !INPUT ARGUMENTS:
  character (LEN=*), intent(in) :: caller
  integer, intent(in) :: sendcount, sendtype, recvcount, recvtype

  !INPUT/OUTPUT ARGUMENTS:
  character (LEN=*)   :: sendbuf
  character (LEN=*)   :: recvbuf
  !EOP

  integer :: sendbyte,recvbyte

  integer, external :: MPI_TYPE_SIZE
  !--------------------------------------------------------------------------

  sendbyte=MPI_TYPE_SIZE(sendtype)*sendcount
  recvbyte=MPI_TYPE_SIZE(recvtype)*recvcount

  if(recvbyte<0 .or. sendbyte<0)then
     write(*,*)'sendtype,recvtype,sendbyte,recvbyte=',&
          sendtype,recvtype,sendbyte,recvbyte

     write(*,*)'Error in MPI_SIMPLE_COPY called from '//caller//&
          'unknown data type(s)'
     stop
  end if

  if(recvbyte<sendbyte)then
     write(*,*)'sendcount,sendtype,recvcount,recvtype,sendbyte,recvbyte=',&
          sendcount,sendtype,recvcount,recvtype,sendbyte,recvbyte
     write(*,*)'Error in MPI_SIMPLE_COPY called from '//caller//&
          ': recvbyte<sendbyte'
     stop
  end if

  recvbuf(1:sendbyte)=sendbuf(1:sendbyte)

end subroutine MPI_SIMPLE_COPY

!BOP ==========================================================================
!ROUTINE: MPI_LOCAL_MSG - used by MPI send and receive subroutines
!INTERFACE:
subroutine MPI_LOCAL_MSG(caller,buf,count,datatype,rank,tag)

  !DESCRIPTION:
  ! This subroutine is used by all MPI SEND and RECV subroutines.
  ! We read and write messages of size count and type datatype from/into buf.
  ! The messages are identified by tag.
  ! Two dynamically resized buffers are used for storage. 
  ! Stop if RECV is issued before the corresponding SEND.
  
  !USES:
  use ModMpiOrig, ONLY: MPI_ANY_TAG
  implicit none

  !INPUT ARGUMENTS:
  character (LEN=*), intent(in)   :: caller
  integer, intent(in)             :: count,datatype,rank,tag

  !INPUT/OUTPUT ARGUMENTS:
  character (LEN=*)               :: buf
  !EOP

  character, dimension(:,:), allocatable, save :: buffer,  buffer2
  integer,   dimension(:),   allocatable, save :: msg_tag, msg_tag2,&
                                                  msg_len, msg_len2

  integer :: max_msg_length=16, max_msg_number=16, i, j

  integer :: msg_length, msg_number=0

  integer, external :: MPI_TYPE_SIZE

  !---------------------------------------------------------------------------

  !write(*,*)'LOCAL_MSG:,caller,count,datatype,rank,tag=',&
  !     caller,count,datatype,rank,tag
 
  if(rank/=0)then
     write(*,*)'count,datatype,tag,rank=',count,datatype,tag,rank
     write(*,*)'Error in MPI_LOCAL_MSG called by '//caller// &
       ': source/dest rank is not 0'
     stop
  end if

  if(tag<0.and.tag/=MPI_ANY_TAG)then
     write(*,*)'count,datatype,tag=',count,datatype,tag
     write(*,*)'Error in MPI_LOCAL_MSG called by '//caller//': invalid tag'
     stop
  end if

  msg_length=MPI_TYPE_SIZE(datatype)*count

  if(msg_length<1)then
     write(*,*)'count,datatype,tag,msg_length=',count,datatype,tag,msg_length
     write(*,*)'Error in MPI_LOCAL_MSG called by '//caller//': msg_length<1'
     stop
  end if

  if(.not.allocated(buffer))then
     !write(*,*)'initial allocation: max_msg_length,max_msg_number=',&
     !     max_msg_length,max_msg_number
     allocate(&
          buffer(max_msg_length,max_msg_number),&
          msg_tag(max_msg_number),&
          msg_len(max_msg_number))
     msg_tag=MPI_ANY_TAG-1
     msg_len=0
  end if

  if(index(caller,'RECV')>0)then
     !write(*,*)'Receiving'

     msg_number=msg_number-1

     do i=1,max_msg_number

        if(msg_tag(i)==tag.or.(tag==MPI_ANY_TAG.and.msg_len(i)>0))then

           !write(*,*)'found corresponding tag at i=',i

           if(msg_length<msg_len(i))then
              write(*,*)'tag,count,datatype,recvlength,sendlength=',&
                   tag,count,datatype,msg_length,msg_len(i)
              write(*,*)'Error in MPI_LOCAL_MSG called by '//caller// &
                   ': receive buffer is shorter than message'
              stop
           end if
           do j=1,msg_len(i)
              buf(j:j)=buffer(j,i)
           end do
           msg_tag(i)=MPI_ANY_TAG-1
           msg_len(i)=0
           RETURN
        end if
     end do
     ! No matching message tag was found
     write(*,*)'tag,count,datatype=',tag,count,datatype
     write(*,*)caller,' was issued before corresponding SEND'
     if(caller=='MPI_IRECV')then
        stop 'Error in MPI_LOCAL_COPY: change order of IRECV/SEND'
     else
        stop 'Error in MPI_LOCAL_COPY: incorrect order of RECV/SEND'
     end if
  end if ! RECV

  ! SEND: store message into buffer

  msg_number=msg_number+1

  !write(*,*)'Sending message, msg_number=',msg_number

  ! Adjust buffer size if necessary
  if(msg_length>max_msg_length)then
     ! Increase first dimension of buffer
     allocate(buffer2(max_msg_length,max_msg_number))
     buffer2=buffer
     deallocate(buffer)
     allocate(buffer(2*msg_length,max_msg_number))
     buffer(1:max_msg_length,:)=buffer2
     deallocate(buffer2)

     !write(*,*)'Increased max_msg_length from,to=',max_msg_length,2*msg_length
     max_msg_length=2*msg_length
  end if

  if(msg_number>max_msg_number)then
     ! Increase second dimension of buffer, and length of msg_tag, msg_len
     allocate(&
          buffer2(max_msg_length,max_msg_number),&
          msg_tag2(max_msg_number),&
          msg_len2(max_msg_number))
     buffer2=buffer
     msg_tag2=msg_tag
     msg_len2=msg_len
     deallocate(buffer,msg_tag,msg_len)
     allocate(&
          buffer(msg_length,2*max_msg_number),&
          msg_tag(2*max_msg_number),&
          msg_len(2*max_msg_number))
     buffer(:,1:max_msg_number)=buffer2
     msg_tag=MPI_ANY_TAG-1
     msg_tag(1:max_msg_number)=msg_tag2
     msg_len=0
     msg_len(1:max_msg_number)=msg_len2
     deallocate(buffer2,msg_tag2,msg_len2)

     !write(*,*)'Increased max_msg_number from,to=',&
     !     max_msg_number,2*max_msg_number
     max_msg_number=2*max_msg_number
  end if

  ! Store message into buffer
  do i=1,msg_number
     if(msg_len(i)==0)then
        ! write(*,*)'Storing message into line i=',i
        do j=1,msg_length
           buffer(j,i)=buf(j:j)
        end do
        msg_tag(i)=tag
        msg_len(i)=msg_length
        RETURN
     end if
  end do

  stop 'Impossible to get here in MPI_LOCAL_MSG'

end subroutine MPI_LOCAL_MSG

!=============================================================================
! subroutines that copy between send and recv buffers
!=============================================================================

subroutine MPI_ALLGATHER(sendbuf,sendcount,sendtype,&
                         recvbuf,recvcount,recvtype,comm,ierror)

  integer, intent(in) :: sendcount,sendtype,recvcount,recvtype,comm
  integer, intent(out):: ierror
  character (LEN=*)   :: sendbuf
  character (LEN=*)   :: recvbuf

  call MPI_SIMPLE_COPY('MPI_ALLGATHER',sendbuf,sendcount,sendtype,&
                                       recvbuf,recvcount,recvtype)

  ierror=0
end subroutine MPI_ALLGATHER

subroutine MPI_ALLGATHERV(sendbuf,sendcount,sendtype,&
     recvbuf,recvcount,recvdisp,recvtype,comm,ierror)

  integer, intent(in) :: sendcount,sendtype,recvcount,recvdisp,recvtype,comm
  integer, intent(out):: ierror
  character (LEN=*)   :: sendbuf
  character (LEN=*)   :: recvbuf

  call MPI_SIMPLE_COPY('MPI_ALLGATHERV',sendbuf,sendcount,sendtype,&
       recvbuf,recvcount,recvtype)

  ierror=0
end subroutine MPI_ALLGATHERV

subroutine MPI_ALLREDUCE(sendbuf,recvbuf,count,datatype,op,comm,ierror)

  integer, intent(in) :: count,datatype,op,comm
  integer, intent(out):: ierror
  character (LEN=*)   :: sendbuf
  character (LEN=*)   :: recvbuf

  call MPI_SIMPLE_COPY('MPI_ALLREDUCE',sendbuf,count,datatype,&
                                       recvbuf,count,datatype)
  ierror=0
end subroutine MPI_ALLREDUCE

subroutine MPI_GATHER(sendbuf,sendcount,sendtype,&
                      recvbuf,recvcount,recvtype,root,comm,ierror)

  integer, intent(in) :: sendcount,sendtype,recvcount,recvtype,root,comm
  integer, intent(out):: ierror
  character (LEN=*)   :: sendbuf
  character (LEN=*)   :: recvbuf

  call MPI_SIMPLE_COPY('MPI_GATHER',sendbuf,sendcount,sendtype,&
                                    recvbuf,recvcount,recvtype)

  ierror=0
end subroutine MPI_GATHER

subroutine MPI_GATHERV(sendbuf,sendcount,sendtype,&
     recvbuf,recvcount,displs,recvtype,root,comm,ierror)

  integer, intent(in) :: sendcount,sendtype,recvcount,displs,recvtype,root,comm
  integer, intent(out):: ierror
  character (LEN=*)   :: sendbuf
  character (LEN=*)   :: recvbuf

  call MPI_SIMPLE_COPY('MPI_GATHERV', &
       sendbuf,sendcount,sendtype,    &
       recvbuf,recvcount,recvtype)

  ierror=0
end subroutine MPI_GATHERV

subroutine MPI_REDUCE(sendbuf,recvbuf,count,datatype,op,root,comm,ierror)

  integer, intent(in) :: count,datatype,op,root,comm
  integer, intent(out):: ierror
  character (LEN=*)   :: sendbuf
  character (LEN=*)   :: recvbuf

  call MPI_SIMPLE_COPY('MPI_REDUCE',sendbuf,count,datatype,&
                                    recvbuf,count,datatype)
  ierror=0
end subroutine MPI_REDUCE


!=============================================================================
! Subroutines that send and receive messages
!=============================================================================

subroutine MPI_IRECV(buf,count,datatype,source,tag,comm,request,ierror)
  integer, intent(in) :: count,datatype,source,tag,comm
  character (LEN=*)   :: buf
  integer, intent(out):: request,ierror

  call MPI_LOCAL_MSG('MPI_IRECV',buf,count,datatype,source,tag)
  ierror=0

end subroutine MPI_IRECV

subroutine MPI_ISEND(buf,count,datatype,dest,tag,comm,request,ierror)
  integer, intent(in) :: count,datatype,dest,tag,comm
  character (len=*)   :: buf
  integer, intent(out):: request,ierror

  call MPI_LOCAL_MSG('MPI_ISEND',buf,count,datatype,dest,tag)
  ierror=0

end subroutine MPI_ISEND

subroutine MPI_IBSEND(buf,count,datatype,dest,tag,comm,request,ierror)
  integer, intent(in) :: count,datatype,dest,tag,comm
  character (len=*)   :: buf
  integer, intent(out):: request,ierror

  call MPI_LOCAL_MSG('MPI_IBSEND',buf,count,datatype,dest,tag)
  ierror=0

end subroutine MPI_IBSEND

subroutine MPI_IRSEND(buf,count,datatype,dest,tag,comm,request,ierror)
  integer, intent(in) :: count,datatype,dest,tag,comm
  character (len=*)   :: buf
  integer, intent(out):: request,ierror

  call MPI_LOCAL_MSG('MPI_IRSEND',buf,count,datatype,dest,tag)
  ierror=0

end subroutine MPI_IRSEND

subroutine MPI_ISSEND(buf,count,datatype,dest,tag,comm,request,ierror)
  integer, intent(in) :: count,datatype,dest,tag,comm
  character (len=*)   :: buf
  integer, intent(out):: request,ierror

  call MPI_LOCAL_MSG('MPI_ISSEND',buf,count,datatype,dest,tag)
  ierror=0

end subroutine MPI_ISSEND

subroutine MPI_RECV(buf,count,datatype,source,tag,comm,status,ierror)
  integer, intent(in) :: count,datatype,source,tag,comm
  character (LEN=*)   :: buf
  integer, intent(out):: status(*),ierror


  call MPI_LOCAL_MSG('MPI_RECV',buf,count,datatype,source,tag)
  ierror=0

end subroutine MPI_RECV

subroutine MPI_SEND(buf,count,datatype,dest,tag,comm,ierror)
  integer, intent(in) :: count,datatype,dest,tag,comm
  character (len=*)   :: buf
  integer, intent(out):: ierror

  call MPI_LOCAL_MSG('MPI_SEND',buf,count,datatype,dest,tag)
  ierror=0

end subroutine MPI_SEND

subroutine MPI_BSEND(buf,count,datatype,dest,tag,comm,ierror)
  integer, intent(in) :: count,datatype,dest,tag,comm
  character (len=*)   :: buf
  integer, intent(out):: ierror

  call MPI_LOCAL_MSG('MPI_BSEND',buf,count,datatype,dest,tag)
  ierror=0

end subroutine MPI_BSEND

subroutine MPI_RSEND(buf,count,datatype,dest,tag,comm,ierror)
  integer, intent(in) :: count,datatype,dest,tag,comm
  character (len=*)   :: buf
  integer, intent(out):: ierror

  call MPI_LOCAL_MSG('MPI_RSEND',buf,count,datatype,dest,tag)
  ierror=0

end subroutine MPI_RSEND

subroutine MPI_SSEND(buf,count,datatype,dest,tag,comm,ierror)
  integer, intent(in) :: count,datatype,dest,tag,comm
  character (len=*)   :: buf
  integer, intent(out):: ierror

  call MPI_LOCAL_MSG('MPI_SSEND',buf,count,datatype,dest,tag)
  ierror=0

end subroutine MPI_SSEND

!=============================================================================
! Subroutines that set some output variables and return
!=============================================================================

subroutine MPI_COMM_RANK(comm,rank,ierror)
  integer, intent(in) :: comm
  integer, intent(out):: rank,ierror
  rank=0
  ierror=0
end subroutine MPI_COMM_RANK

subroutine MPI_COMM_SIZE(comm,size,ierror)
  integer, intent(in) :: comm
  integer, intent(out):: size, ierror
  size=1
  ierror=0
end subroutine MPI_COMM_SIZE

!=============================================================================
! Subroutines that set ierror=0 only and simply return
!=============================================================================

subroutine MPI_BARRIER(comm, ierror)
  integer, intent(in) :: comm
  integer, intent(out):: ierror
  ierror=0
end subroutine MPI_BARRIER

subroutine MPI_BCAST(buffer,count,datatype,root,comm,ierror)
  integer, intent(in) :: count,datatype,root,comm
  character (LEN=*)   :: buffer
  integer, intent(out):: ierror
  ierror=0
end subroutine MPI_BCAST

subroutine MPI_CART_CREATE(comm_old,ndims,dims,periods,reorder,comm_cart,&
     ierror)
  integer, intent(in) :: comm_old,ndims
  integer, intent(in) :: dims(ndims)
  logical, intent(in) :: periods(ndims),reorder
  integer, intent(out):: comm_cart,ierror
  comm_cart=comm_old
  ierror=0
end subroutine MPI_CART_CREATE

subroutine MPI_CART_GET(comm, maxdims, dims, periods, coords, ierror)
  integer, intent(in) :: comm, maxdims
  integer, intent(out):: dims(maxdims),coords(maxdims),ierror
  logical, intent(out):: periods(maxdims)
  ierror=0
end subroutine MPI_CART_GET

subroutine MPI_CART_RANK(comm, coords, rank, ierror)
  integer, intent(in) :: comm, coords(*)
  integer, intent(out):: rank, ierror
  rank=0
  ierror=0
end subroutine MPI_CART_RANK

subroutine MPI_CART_SHIFT(comm,direction, disp, rank_source, rank_dest, ierror)
  integer, intent(in) :: comm,direction, disp
  integer, intent(out):: rank_source, rank_dest, ierror
  rank_source=0
  rank_dest=0
  ierror=0
end subroutine MPI_CART_SHIFT

subroutine MPI_COMM_CREATE(comm,group,newcomm,ierror)
  integer, intent(in) :: comm,group
  integer, intent(out):: newcomm,ierror
  newcomm=0
  ierror=0
end subroutine MPI_COMM_CREATE

subroutine MPI_COMM_FREE(comm,ierror)
  integer, intent(inout):: comm
  integer, intent(out)  :: ierror
  ierror=0 
end subroutine MPI_COMM_FREE

subroutine MPI_COMM_GROUP(comm,group,ierror)
  integer, intent(in) :: comm
  integer, intent(out):: group,ierror
  group=0
  ierror=0
end subroutine MPI_COMM_GROUP

subroutine MPI_ERRHANDLER_SET(comm,errhandler,ierror)
  integer, intent(in) :: comm, errhandler
  integer, intent(out):: ierror
  ierror=0
end subroutine MPI_ERRHANDLER_SET

subroutine MPI_ERROR_CLASS(errorcode,errorclass,ierror)
  integer, intent(in) :: errorcode
  integer, intent(out):: errorclass, ierror
  errorclass=0
  ierror=0
end subroutine MPI_ERROR_CLASS

subroutine MPI_FINALIZE(ierror)
  integer, intent(out):: ierror
  ierror=0
end subroutine MPI_FINALIZE

subroutine MPI_GROUP_DIFFERENCE(group1,group2,newgroup,ierror)
  integer, intent(in) :: group1, group2
  integer, intent(out):: newgroup, ierror
  newgroup=0
  ierror=0
end subroutine MPI_GROUP_DIFFERENCE

subroutine MPI_GROUP_FREE(group,ierror)
  integer, intent(inout):: group
  integer, intent(out)  :: ierror
  group=0
  ierror=0
end subroutine MPI_GROUP_FREE

subroutine MPI_GROUP_INCL(group,n,ranks,newgroup,ierror)
  integer, intent(in) :: group, n, ranks
  integer, intent(out):: newgroup, ierror
  newgroup=0
  ierror=0
end subroutine MPI_GROUP_INCL

subroutine MPI_GROUP_RANGE_INCL(group,n,ranges,newgroup,ierror)
  integer, intent(in) :: group, n
  integer, intent(in) :: ranges(3,n)
  integer, intent(out):: newgroup, ierror
  newgroup=0
  ierror=0
end subroutine MPI_GROUP_RANGE_INCL

subroutine MPI_GROUP_RANK(group,rank,ierror)
  integer, intent(in) :: group
  integer, intent(out):: rank, ierror
  rank=0
  ierror=0
end subroutine MPI_GROUP_RANK

subroutine MPI_GROUP_SIZE(group, size, ierror)

  integer, intent(in) :: group
  integer, intent(out):: size, ierror

  size   = 1
  ierror = 0

end subroutine MPI_GROUP_SIZE

subroutine MPI_GROUP_TRANSLATE_RANKS(group1,n,ranks1,group2,ranks2,ierror)

  integer, intent(in) :: group1, n
  integer, intent(in) :: ranks1(n)
  integer, intent(in) :: group2
  integer, intent(out):: ranks2(n)
  integer, intent(out):: ierror

  ranks2 = ranks1
  ierror = 0

end  subroutine MPI_GROUP_TRANSLATE_RANKS

subroutine MPI_GROUP_UNION(group1, group2, newgroup, ierror)

  integer, intent(in) :: group1, group2
  integer, intent(out):: newgroup, ierror

  newgroup = group1
  ierror   = 0

end subroutine MPI_GROUP_UNION

subroutine MPI_INIT(ierror)
  integer, intent(out):: ierror
  ierror=0
end subroutine MPI_INIT

subroutine MPI_WAIT(request,status,ierror)
  integer, intent(inout) :: request
  integer, intent(out)   :: status(*), ierror
  ierror=0
end subroutine MPI_WAIT

subroutine MPI_WAITALL(count,array_of_requests,array_of_statuses,ierror)
  integer, intent(in)   :: count
  integer, intent(inout):: array_of_requests(*)
  integer, intent(out)  :: array_of_statuses(*), ierror
  ierror=0
end subroutine MPI_WAITALL

!=============================================================================
! Subroutines that stop
!=============================================================================

subroutine MPI_ABORT(errorcode, comm, ierror)
  integer, intent(in) :: errorcode, comm
  integer, intent(out):: ierror

  write(*,*)'errcode=',errcode
  stop
end subroutine MPI_ABORT

!=============================================================================
! Functions 
!=============================================================================

real*8 function MPI_WTIME()
  integer:: clock,clockrate,count_max
  call system_clock(clock,clockrate,count_max)
  MPI_WTIME=dble(clock)/clockrate
end function MPI_WTIME

!=============================================================================
! MPI_FILE*
!=============================================================================

subroutine MPI_FILE_WRITE

  write(*,*)'MPI_FILE_WRITE is not implemented in NOMPI library'

end subroutine MPI_FILE_WRITE

subroutine MPI_FILE_READ

  write(*,*)'MPI_FILE_WRITE is not implemented in NOMPI library'

end subroutine MPI_FILE_READ
