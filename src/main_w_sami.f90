!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-----------------------------------------------------------------------------
! $Id: main.f90,v 1.10 2013/10/24 18:36:35 agburr Exp $
!
! Author: Aaron Ridley
!
! Comments: GITM main
!
! AGB 3/31/13: Added call for RCMR data assimilation
! AGB 10/23/13: Adapted RCMR call to new format
!-----------------------------------------------------------------------------

program GITM

  use ModInputs
  use ModTime
  use ModGITM
  use ModMpi
  use ModRCMR
  use ModSatellites, only: SatCurrentDat, SatAltDat, nRCMRSat
  use ModEUV, only: sza
  use advance_mod ,only: hrut,denit
  use ModCoupSAMI3
  use ModSamiInterp
  use coup_mod, only: IsCoupGITM
  use namelist_mod,only:hrinit

  use message_passing_mod, only: taskid, master_sami, iCommSami, iStatEdSimlat

  implicit none

  integer :: iBlock

  integer :: iError

  logical :: IsValidTime = .true.
  integer :: iprocglob,nprocsglob

  integer :: status(MPI_STATUS_SIZE)

  real :: xxx

  integer :: nIter= 1

  ! ------------------------------------------------------------------------
  ! initialize stuff
  ! ------------------------------------------------------------------------

  call init_mpi_coup

  call MPI_COMM_RANK(iCommGlobal, iprocglob, iError)
  call MPI_COMM_SIZE(iCommGlobal, nprocsglob, iError)

  ! initialize GITM

  if (iCommGITM /= MPI_COMM_NULL) then

     call init_mpi_gitm

     call start_timing("GITM")
     call delete_stop

     call init_planet
     call set_defaults

     call read_inputs(cInputFile)
     call set_inputs

     call initialize_gitm(CurrentTime)

     call write_output

     call report("Starting Main Time Loop",0)

  else if (iCommSAMI0 /= MPI_COMM_NULL) then

     call init_sami(iCommSAMI0,IsValidTime)

     if (taskid ==0) master_sami = iprocglob 

  endif

  ! broadcast CurrentTime to global procs
  call mpi_bcast(EndTime,1,MPI_REAL,0,iCommGlobal,iError)
  call mpi_bcast(CurrentTime,1,MPI_REAL,0,iCommGlobal,iError)
  call mpi_bcast(StartTime,1,MPI_REAL,0,iCommGlobal,iError)
  call mpi_bcast(DtCouple,1,MPI_REAL,0,iCommGlobal,iError)

  ! ------------------------------------------------------------------------
  ! Run for a few iterations
  ! ------------------------------------------------------------------------

  call MPI_BARRIER(iCommGlobal,iError)

  do while(CurrentTime < EndTime)

     if (iCommGITM /= MPI_COMM_NULL) then

        do while (CurrentTime < (StartTime+nIter*DtCouple))

           !write(*,*) '>>>GITM start iterations---',nIter,iproc ,&
           !    CurrentTime,EndTime
           call calc_pressure

!!! We may have to split cMax and Dt calculation!!!
           if(RCMRFlag) then
              Dt = 2
           else
              Dt = FixedDt
           end if

           call calc_timestep_vertical
           if (.not. Is1D) call calc_timestep_horizontal

           if(RCMRFlag) then
              call run_RCMR
           endif

           call advance

           if (.not.IsFramework) then
              call check_stop
           endif

           iStep = iStep + 1

           call write_output

        end do


     else if (iCommSAMI0 /= MPI_COMM_NULL) then

        if (taskid == 0) call advance_sami_pe0
        if (taskid>0) then

           do while ((IsValidTime .and. &
                ((hrut-hrinit)*3600.0 < (nIter*DtCouple-1.0e-6)))) 

              call Sami_run(IsValidTime,DtCouple)

           enddo

           if (taskid == 1) &
                print*, '---- SAMI3 completed one iter',(hrut-hrinit)*3600.0,nIter*DtCouple,IsValidTime,taskid
           ! take master out of advance
           xxx = 1.0
           call mpi_send(xxx, 1, MPI_REAL, 0,iStatEdSimlat ,&
                iCommSami, iError)
           if (taskid == 1) &
                write(*,*) '---- Done with sending xxx !!!',taskid
        endif
     endif

     ! turn on coupling flags and exchange data

     call MPI_BARRIER(iCommGlobal,iError)
     if (iCommGITM /= MPI_COMM_NULL) then

        StartCoup = .true.
        StartCoup_pot = .true.

        if (StartCoup) & 
             call gather_data_gitm

     else if (iCommSAMI0 /= MPI_COMM_NULL) then

        if (taskid > 0) IsCoupGITM = .true.

     endif

     call MPI_BARRIER(iCommGlobal,iError)

     call ExchangeData

     if (iprocglob == 0) print*,'-->> Done ExchangeData',iprocglob

     if (iCommSAMI0 /= MPI_COMM_NULL) then
        if (taskid >0 .and. IsCoupGITM) &
             call init_coup
     endif

     call MPI_BARRIER(iCommGlobal,iError)

     call mpi_bcast(CurrentTime,1,MPI_REAL,0,iCommGlobal,iError)

     nIter = nIter+1

  end do !(end while gobal)

  ! ------------------------------------------------------------------------
  ! Finish run
  ! ------------------------------------------------------------------------

  if (iCommSAMI0 /= MPI_COMM_NULL) call finalize_sami
  if (iCommGITM /= MPI_COMM_NULL) call finalize_gitm

  ! ---------------------------------------------------------------------
  ! MPI clean up
  ! ---------------------------------------------------------------------

  ! call MPI_BARRIER(iCommGlobal,iError)

  if (iCommGITM /= MPI_COMM_NULL) call MPI_Comm_free(iCommGITM,iError)
  if (iCommSAMI0 /= MPI_COMM_NULL) call MPI_Comm_free(iCommSAMI0,iError)

  call MPI_FINALIZE(iError)

end program GITM

!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================

subroutine CON_stop(StringError)

  implicit none
  character (len=*), intent(in) :: StringError
  call stop_gitm(StringError)

end subroutine CON_stop

subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe

  DoTest = .false.; DoTestMe = .false.

end subroutine CON_set_do_test

subroutine CON_io_unit_new(iUnit)

  implicit none
  integer, intent(in) :: iUnit

  return

end subroutine CON_io_unit_new

!---------------------------------------------------------------------------

