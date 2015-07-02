!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine exchange_messages_sphere

  use ModSphereInterface

  implicit none

  logical :: isOk
  integer :: iBlock

  call report("Exchange Messages",2)
  call start_timing("Message Pass")

  call UAM_XFER_start(ok=isOK)
  if (.not. isOK) then
     call UAM_write_error()
     call stop_gitm("error in xfer_start")
  endif

  !
  ! End the message passing
  !

  call UAM_XFER_finish(ok=isOK)
  if (.not. isOK) then
     call UAM_write_error()
     call stop_gitm("error in xfer_finish")
  endif

  do iBlock = 1, nBlocks

     call report("Starting to check a new block",3)

     if (Latitude(0,iBlock) < -pi/2.0) then

        Velocity(:,-1:0,:,iNorth_,iBlock) = &
             -Velocity(:,-1:0,:,iNorth_,iBlock)
        Velocity(:,-1:0,:,iEast_,iBlock) = &
             -Velocity(:,-1:0,:,iEast_,iBlock)

     endif

     if (Latitude(nLats+1,iBlock) > pi/2.0) then

        Velocity(:,nLats+1:nLats+2,:,iNorth_,iBlock) = &
             -Velocity(:,nLats+1:nLats+2,:,iNorth_,iBlock)

        Velocity(:,nLats+1:nLats+2,:,iEast_,iBlock) = &
             -Velocity(:,nLats+1:nLats+2,:,iEast_,iBlock)

     endif

     if (minval(temperature(1:nLons,1:nLats,1:nAlts,iBlock)) < 100.0) then
        write(*,*) "Low Temperature : ",iBlock, &
             minval(temperature(1:nLons,1:nLats,1:nAlts,iBlock))
        call stop_gitm("Temperature < 100.0")
     endif

     if (minval(rho(1:nLons,1:nLats,1:nAlts,iBlock)) < 0.0) then
        write(*,*) "Low mass density : ",iBlock, &
             minval(rho(1:nLons,1:nLats,1:nAlts,iBlock))
        call stop_gitm("mass density < 0.0")
     endif

  enddo

  call end_timing("Message Pass")

end subroutine exchange_messages_sphere
