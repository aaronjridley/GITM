!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP
!MODULE: CON_ray_trace - trace field and stream lines in parallel
!INTERFACE:
module CON_ray_trace

  !DESCRIPTION:
  ! Provides the infrastructure for parallel tracing of a vector field:
  ! for velocity field stream lines, for magnetic field the field lines,
  ! in general rays.
  !
  ! Each processor is working on the ray segments inside their subdomain.
  ! The processors periodically exchange ray information with the other 
  ! processors. The ray information contains the starting position,
  ! the rank of the starting processor, the current position, 
  ! the direction of the ray relative to the vector field (parallel or
  ! antiparallel), the status of the ray tracing (done or in progress).

  !USES:
  use ModMpi

  implicit none

  save    ! save all variables

  private ! except

  !PUBLIC MEMBER FUNCTIONS:
  public :: ray_init          ! Initialize module
  public :: ray_clean         ! Clean up storage
  public :: ray_put           ! Put information about a ray to be continued
  public :: ray_get           ! Get information about a ray to be continued
  public :: ray_exchange      ! Exchange ray information with other processors
  public :: ray_test          ! Unit tester

  !REVISION HISTORY:
  ! 31Jan04 - Gabor Toth <gtoth@umich.edu> - initial prototype/prolog/code
  ! 26Mar04 - Gabor Toth added the passing of nValue real values
  ! 31Mar04 - Gabor Toth removed the passing of nValue real values
  !                      changed XyzStart_D(3) into iStart_D(4)
  !EOP

  ! Private constants
  character(len=*),  parameter :: NameMod='CON_ray_trace'

  ! Named indexes for ray position
  integer, parameter :: &
       RayStartI_=1, RayStartJ_=2, RayStartK_=3, &
       RayStartBlock_=4, RayStartProc_=5, &
       RayEndX_=6, RayEndY_=7, RayEndZ_=8, RayLength_=9, RayDir_=10, &
       RayDone_=11

  ! Minimum dimensionality of ray info
  integer, parameter :: nRayInfo = 11

  ! Private type
  ! Contains nRay rays with full information
  type RayPtr
     integer       :: nRay, MaxRay
     real, pointer :: Ray_VI(:,:)
  end type RayPtr

  ! Private variables

  ! Ray buffers
  type(RayPtr), pointer :: Send_P(:)     ! Rays to send to other PE-s
  type(RayPtr)          :: Recv          ! Rays received from all other PE-s
  integer, pointer      :: nRayRecv_P(:) ! Number of rays recv from other PE-s

  ! MPI variables
  integer               :: iComm=MPI_COMM_NULL, nProc, iProc ! MPI group info
  integer               :: nRequest
  integer, allocatable  :: iRequest_I(:), iStatus_II(:,:)

  integer               :: iError              ! generic MPI error value

contains

  !BOP =======================================================================
  !IROUTINE: ray_init - (re)initialize CON_ray_trace
  !INTERFACE:
  subroutine ray_init(iCommIn)

    !INPUT ARGUMENTS:
    integer, intent(in) :: iCommIn  ! MPI communicator for the processors

    !DESCRIPTION:
    ! Initialize the ray tracing for the MPI communicator iCommIn.
    ! If called multiple times, it checks if iCommIn is different from the
    ! current communicator. If the same, nothing is done, if different,
    ! the current storage is removed and new storage is allocated.
    !EOP

    character (len=*), parameter :: NameSub = NameMod//'::ray_init'
    integer :: jProc
    !------------------------------------------------------------------------

    if(iComm == iCommIn) RETURN               ! nothing to do if the same comm

    if(iComm /= MPI_COMM_NULL) call ray_clean ! clean up previous allocation

    ! Store MPI information
    iComm = iCommIn
    call MPI_COMM_SIZE (iComm, nProc,  iError)
    call MPI_COMM_RANK (iComm, iProc,  iError)

    ! Allocate MPI variables for non-blocking exchange
    ! At most nProc messages are sent
    allocate(iRequest_I(nProc), iStatus_II(MPI_STATUS_SIZE, nProc))

    ! Initialize send buffers
    allocate(Send_P(0:nProc-1))
    do jProc=0, nProc-1
       Send_P(jProc) % nRay   = 0
       Send_P(jProc) % MaxRay = 0
       nullify(Send_P(jProc) % Ray_VI)
    end do

    ! Initialize receive buffer
    allocate(nRayRecv_P(0:nProc-1))
    Recv % nRay   = 0
    Recv % MaxRay = 0
    nullify(Recv % Ray_VI)

  end subroutine ray_init

  !BOP ========================================================================
  !IROUTINE: ray_clean - clean up storage for CON_ray_trace
  !INTERFACE:
  subroutine ray_clean

    !DESCRIPTION:
    ! Remove storage allocated in CON\_ray\_trace.
    !EOP

    character (len=*), parameter :: NameSub = NameMod//'::ray_clean'
    integer :: jProc
    !------------------------------------------------------------------------

    ! Deallocate MPI variables
    deallocate(iRequest_I, iStatus_II)

    ! Deallocate send buffers
    if(associated(Send_P))then
       do jProc = 0, nProc-1
          if(associated(Send_P(jProc) % Ray_VI)) &
               deallocate(Send_P(jProc) % Ray_VI)
       end do
    end if
    deallocate(Send_P)

    ! Deallocate recv buffer
    if(associated(nRayRecv_P))deallocate(nRayRecv_P)
    if(associated(Recv % Ray_VI))deallocate(Recv % Ray_VI)
    Recv % nRay   = 0
    Recv % MaxRay = 0

    iComm = MPI_COMM_NULL

  end subroutine ray_clean

  !BOP =======================================================================
  !IROUTINE: ray_put - put ray information into send buffer
  !INTERFACE:
  subroutine ray_put(&
       iProcStart,iStart_D,iProcEnd,XyzEnd_D,Length,IsParallel,DoneRay)

    !INPUT ARGUMENTS:
    integer, intent(in) :: iProcStart,iProcEnd ! PE-s for start and end pos.
    integer, intent(in) :: iStart_D(4)         ! Indexes i,j,k,iBlock for start
    real,    intent(in) :: XyzEnd_D(3)         ! End posistion
    real,    intent(in) :: Length              ! Length of the ray so far
    logical, intent(in) :: IsParallel,DoneRay  ! Direction and status of trace

    !DESCRIPTION:
    ! Put ray information into send buffer. If DoneRay is true, the
    ! information will be sent ot iProcStart, otherwise it will be
    ! sent to iProcEnd. 
    !EOP

    character (len=*), parameter :: NameSub = NameMod//'::ray_put'

    integer :: iProcTo
    !-------------------------------------------------------------------------
    ! Where should we send the ray
    if(DoneRay)then
       iProcTo = iProcStart  ! Send back result to the PE that started tracing
    else
       iProcTo = iProcEnd    ! Send to PE which can continue the tracing
    end if

    if(iProcTo<0)&
         call CON_stop(NameSub//' SWMF_error: PE lookup to be implemented')

    call append_ray(Send_P(iProcTo))

  contains

    !========================================================================
    subroutine append_ray(Send)

      ! Append a new element to the Send buffer

      type(RayPtr) :: Send
      integer      :: iRay
      !---------------------------------------------------------------------
      iRay = Send % nRay + 1
      if(iRay > Send % MaxRay) call extend_buffer(Send, iRay+100)

      Send % Ray_VI(RayStartI_:RayStartBlock_,iRay) = real(iStart_D)
      Send % Ray_VI(RayStartProc_            ,iRay) = iProcStart
      Send % Ray_VI(RayEndX_:RayEndZ_        ,iRay) = XyzEnd_D
      Send % Ray_VI(RayLength_               ,iRay) = Length

      if(IsParallel)then
         Send % Ray_VI(RayDir_,iRay)                =  1
      else
         Send % Ray_VI(RayDir_,iRay)                = -1
      end if

      if(DoneRay)then
         Send % Ray_VI(RayDone_,iRay)               =  1
      else
         Send % Ray_VI(RayDone_,iRay)               =  0
      end if

      Send % nRay = iRay

    end subroutine append_ray

  end subroutine ray_put

  !BOP =======================================================================
  !IROUTINE: ray_get - get ray information from the receive buffer
  !INTERFACE:
  subroutine ray_get(&
       IsFound,iProcStart,iStart_D,XyzEnd_D,Length,IsParallel,DoneRay)

    !OUTPUT ARGUMENTS:
    logical, intent(out) :: IsFound            ! true if there are still rays
    integer, intent(out) :: iProcStart         ! PE-s for start and end pos.
    integer, intent(out) :: iStart_D(4)        ! Indexes i,j,k,iBlock for start
    real,    intent(out) :: XyzEnd_D(3)        ! End position
    real,    intent(out) :: Length             ! Length of the current ray
    logical, intent(out) :: IsParallel,DoneRay ! Direction and status of trace

    !DESCRIPTION:
    ! Provide the last ray for the component to store or to work on.
    ! If no ray is found in the receive buffer, IsFound=.false. is returned.
    !EOP

    character (len=*), parameter :: NameSub = NameMod//'::ray_get'

    integer :: iRay
    !-------------------------------------------------------------------------

    iRay    = Recv%nRay 
    IsFound = iRay > 0

    if(.not.IsFound) RETURN  ! No more rays in the buffer

    ! Copy last ray into output arguments
    iStart_D     = nint(Recv % Ray_VI(RayStartI_:RayStartBlock_,iRay))
    iProcStart   = nint(Recv % Ray_VI(RayStartProc_,iRay))
    XyzEnd_D     = Recv % Ray_VI(RayEndX_:RayEndZ_,iRay)
    Length       = Recv % Ray_VI(RayLength_, iRay)
    IsParallel   = Recv % Ray_VI(RayDir_,iRay)  > 0.0
    DoneRay      = Recv % Ray_VI(RayDone_,iRay) > 0.5

    ! Remove ray from buffer
    Recv % nRay = iRay - 1

  end subroutine ray_get

  !BOP =======================================================================
  !IROUTINE: ray_exchange - send information from send to receive buffers
  !INTERFACE:
  subroutine ray_exchange(DoneMe,DoneAll)

    !INPUT ARGUMENTS:
    logical, intent(in) :: DoneMe

    !OUTPUT ARGUMENTS:
    logical, intent(out):: DoneAll

    !DESCRIPTION:
    ! Send the Send\_P buffers to Recv buffers, empty the Send\_P buffers.
    ! Also check if there is more work to do. If the input argument DoneMe 
    ! is false on any PE-s (i.e. it has more rays to do), 
    ! or if there are any rays passed, the output argument DoneAll is 
    ! set to .false. for all PE-s. 
    ! If all PE-s have DoneMe true and all send buffers are
    ! empty then DoneAll is set to .true.
    !EOP

    logical :: DoneLocal

    integer, parameter :: iTag = 1
    integer :: jProc, iRay, nRayRecv

    character (len=*), parameter :: NameSub = NameMod//'::ray_exchange'

    !-------------------------------------------------------------------------

    ! Exchange number of rays in the send buffer

    ! Local copy (in case ray remains on the same PE)
    nRayRecv_P(iProc)=Send_P(iProc) % nRay

    nRequest = 0
    iRequest_I = MPI_REQUEST_NULL
    do jProc = 0, nProc-1
       if(jProc==iProc)CYCLE
       nRequest = nRequest + 1
       call MPI_irecv(nRayRecv_P(jProc),1,MPI_INTEGER,jProc,&
            iTag,iComm,iRequest_I(nRequest),iError)
    end do

    ! Wait for all receive commands to be posted for all processors
    call MPI_barrier(iComm,iError)

    ! Use ready-send
    do jProc = 0, nProc-1
       if(jProc==iProc)CYCLE
       call MPI_rsend(Send_P(jProc) % nRay,1,MPI_INTEGER,jProc,&
            iTag,iComm,iError)
    end do

    ! Wait for all messages to be received
    if(nRequest > 0)call MPI_waitall(nRequest,iRequest_I,iStatus_II,iError)

    nRayRecv = Recv % nRay + sum(nRayRecv_P)

    ! Extend receive buffer as needed
    if(nRayRecv > Recv % MaxRay) call extend_buffer(Recv,nRayRecv+100)

    ! Exchange ray information
    iRay = Recv % nRay + 1

    ! Local copy if any
    if(nRayRecv_P(iProc) > 0)then
       Recv % Ray_VI(:,iRay:iRay+nRayRecv_P(iProc)-1) = &
            Send_P(iProc) % Ray_VI(:,1:Send_P(iProc) % nRay)
       iRay = iRay + nRayRecv_P(iProc)
    end if

    nRequest   = 0
    iRequest_I = MPI_REQUEST_NULL
    do jProc = 0, nProc-1
       if(jProc==iProc)CYCLE
       if(nRayRecv_P(jProc)==0)CYCLE
       nRequest = nRequest + 1

       call MPI_irecv(Recv % Ray_VI(1,iRay),nRayRecv_P(jProc)*nRayInfo,&
               MPI_REAL,jProc,iTag,iComm,iRequest_I(nRequest),iError)
       iRay = iRay + nRayRecv_P(jProc)
    end do

    ! Wait for all receive commands to be posted for all processors
    call MPI_barrier(iComm, iError)

    do jProc = 0, nProc-1
       if(jProc==iProc)CYCLE
       if(Send_P(jProc) % nRay == 0) CYCLE

       call MPI_rsend(Send_P(jProc) % Ray_VI(1,1),&
            Send_P(jProc) % nRay*nRayInfo,MPI_REAL,jProc,iTag,iComm,iError)
    enddo

    ! Wait for all messages to be received
    if(nRequest > 0)call MPI_waitall(nRequest,iRequest_I,iStatus_II,iError)

    ! Update number of received rays
    Recv % nRay = nRayRecv

    ! Reset send buffers
    do jProc = 0, nProc-1
       Send_P(jProc) % nRay = 0
    end do

    ! Check if there is more work to do on this PE
    DoneLocal = DoneMe .and. (Recv % nRay == 0)

    ! Check if all PE-s are done
    call MPI_allreduce(DoneLocal, DoneAll, 1, MPI_LOGICAL, MPI_LAND, &
         iComm, iError)

  end subroutine ray_exchange

  !===========================================================================

  subroutine extend_buffer(Buffer, nRayNew)

    ! Extend buffer size to nRayNew (or more)

    type(RayPtr)        :: Buffer
    integer, intent(in) :: nRayNew

    real, pointer :: OldRay_VI(:,:)
    !------------------------------------------------------------------------
    if(.not.associated(Buffer % Ray_VI))then
       allocate(Buffer % Ray_VI(nRayInfo,nRayNew))    ! allocatenew buffer
       Buffer % nRay   = 0                            ! buffer is empty
       Buffer % MaxRay = nRayNew                      ! set buffer size
    else
       OldRay_VI => Buffer % Ray_VI                   ! store old values
       allocate(Buffer % Ray_VI(nRayInfo,nRayNew))    ! allocate new buffer
       Buffer % Ray_VI(:,1:Buffer % nRay) = &
            OldRay_VI(:,1:Buffer % nRay)              ! copy old values
       deallocate(OldRay_VI)                          ! free old storage
       Buffer % MaxRay = nRayNew                      ! change buffer size
    end if

  end subroutine extend_buffer

  !BOP ========================================================================
  !IROUTINE: ray_test - unit tester
  !INTERFACE:
  subroutine ray_test

    !DESCRIPTION:
    ! Test the CON\_ray\_trace module. This subroutine should be called from
    ! a small stand alone program.
    !EOP

    logical :: IsFound
    integer :: iProcStart
    integer :: iStart_D(4)
    real    :: XyzEnd_D(3), Length
    logical :: IsParallel,DoneRay, DoneAll
    integer :: jProc
    !------------------------------------------------------------------------

    call ray_init(MPI_COMM_WORLD)

    if(iProc==0) write(*,'(a,i4,i4)')'ray_init done, iProc,nProc=',iProc,nProc

    write(*,"(a,i2,a,4i4,a,i2,a,3f5.0,a,f5.0,a,2l2)") &
         " Sending from iProc=",iProc,&
         " iStart=",(/110+iProc,120+iProc,130+iProc,140+iProc/),&
         " to jProc=",mod(iProc+1,nProc),&
         " XyzEnd=",(/210.+iProc,220.+iProc,230.+iProc/), &
         " Length=",10.0*iProc, &
         " IsParallel, DoneRay=",.true.,.false.

    call ray_put(iProc, (/110+iProc,120+iProc,130+iProc,140+iProc/), &
         mod(iProc+1,nProc), (/210.+iProc,220.+iProc,230.+iProc/), &
         10.0*iProc, .true.,.false.)

    if(iProc==0) write(*,"(a,i2,a,4i4,a,i2,a,3f5.0,a,f5.0,a,2l2)") &
         " Sending from iProc=",iProc,&
         " iStart=",(/110+iProc,120+iProc,130+iProc,140+iProc/),&
         " to jProc=",mod(nProc+iProc-1,nProc),&
         " XyzEnd=",(/210.+iProc,220.+iProc,230.+iProc/), &
         " Length=",10.0*iProc+1.0, &
         " IsParallel, DoneRay=",.false.,.false.

    call ray_put(iProc, (/110+iProc,120+iProc,130+iProc,140+iProc/), &
         mod(nProc+iProc-1,nProc), (/210.+iProc,220.+iProc,230.+iProc/), &
         10.0*iProc+1.0,.false.,.false.)

    do jProc = 0, nProc-1
       write(*,"(a,i2,i2,i4,i4,100f5.0)")'iProc,jProc,Send_P(jProc)=',&
            iProc,jProc,Send_P(jProc) % MaxRay,&
            Send_P(jProc) % nRay, &
            Send_P(jProc) % Ray_VI(:,1:Send_P(jProc) % nRay)
    end do

    if(iProc==0) write(*,'(a)')'ray_put done'

    call ray_exchange(.true., DoneAll)

    write(*,*)'ray_exchange done, DoneAll=',DoneAll

    do
       call ray_get(IsFound,iProcStart,iStart_D,XyzEnd_D,Length, &
            IsParallel,DoneRay)
       if(.not.IsFound) EXIT
       write(*,"(a,i2,a,4i4,a,i2,a,3f5.0,a,f5.0,a,2l2)")&
            'iProc ',iProc,' received iStart=',iStart_D,&
            ' iProcStart=',iProcStart,' XyzEnd=',XyzEnd_D,' Length=',Length,&
            ' Isparallel,DoneRay=',IsParallel,DoneRay
    end do

    if(iProc==0) write(*,'(a)')'ray_get done'

    call ray_exchange(.true., DoneAll)

    if(iProc==0) write(*,'(a,l1)')'ray_exchange repeated, DoneAll=',DoneAll

    call ray_clean

    if(iProc==0) write(*,'(a)')'ray_clean done'

  end subroutine ray_test

  !===========================================================================

end module CON_ray_trace
