!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module UA_wrapper

  ! Wrapper for GITM Upper Atmosphere (UA) component

  implicit none

  private ! except

  public:: UA_set_param
  public:: UA_init_session
  public:: UA_run
  public:: UA_save_restart
  public:: UA_finalize

contains

  !============================================================================
  subroutine UA_set_param(CompInfo, TypeAction)

    use ModInputs, only: cInputText
    use ModReadParam, only: read_text, n_line_read

    use ModTime, ONLY: StartTime, tSimulation, CurrentTime
    use ModInputs, only: iStartTime, IsFramework, iOutputUnit_, set_defaults, &
         nInputLines
    use ModTimeConvert, ONLY: time_real_to_int
    use CON_physics,    ONLY: get_time
    use ModIoUnit
    use ModProcUA
    use ModGITM, only: iCommGITM, nProcs, iProcGITM => iProc
    use ModPlanet, only: init_planet
    use CON_comp_info
    use ModUtilities, ONLY: check_dir

    character (len=*), parameter :: NameSub='UA_set_param'

    ! Arguments
    type(CompInfoType), intent(inout):: CompInfo   ! Information for this comp.
    character (len=*), intent(in)    :: TypeAction ! What to do

    integer :: iError

    iError = 0

    !-------------------------------------------------------------------------
    select case(TypeAction)
    case('VERSION')

       call put(CompInfo,&
            Use=.true.,                                      &
            NameVersion='Global Iono-Thermo Model (Ridley)', &
            Version=2.0)

    case('MPI')

       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)

       iCommGITM = iComm
       iProcGITM = iProc
       nProcs    = nProc

       if(iProc==0)then
          call check_dir("UA/DataIn")
          call check_dir("UA/data")
          call check_dir("UA/RestartOUT")
       end if

       IsFramework = .true.

       call init_planet
       call set_defaults

    case('READ')

       call read_text(cInputText)
       cInputText(n_line_read()+1) = "#END"
       nInputLines=n_line_read()+1

       call set_inputs

    case('CHECK')

       iError = 0

    case('STDOUT')

       !     iUnitStdOut=STDOUT_
       !     if(nProc==1)then
       !        StringPrefix='UA:'
       !     else
       !        write(StringPrefix,'(a,i3.3,a)')'UA',iProc,':'
       !     end if

    case('FILEOUT')

       call get(CompInfo,iUnitOut=iOutputUnit_)
       !     StringPrefix=''

    case('GRID')

       call UA_set_grid

    case default

       call CON_stop(NameSub//' UA_ERROR: invalid TypeAction='//TypeAction)

    end select

    if (iError /= 0) &
         call CON_stop(NameSub//' UA_ERROR in TypeAction='//TypeAction)

  end subroutine UA_set_param

  !============================================================================

  subroutine UA_set_grid

    ! Set the grid descriptor for UA
    ! Since UA has a static grid the descriptor has to be set once.
    ! There can be many couplers that attempt to set the descriptor,
    ! so we must check IsInitialized.
    use ModProcUA
    use CON_Coupler
    use CON_comp_info
    use ModNumConst
    use ModSizeGitm
    use ModSphereInterface, only: iStartBLK
    use ModInputs, only: nBlocksLat, nBlocksLon
    use ModUtilities, ONLY: check_allocate

    character (len=*), parameter :: NameSub='UA_set_grid'
    logical :: IsInitialized=.false.
    integer, allocatable :: iProc_A(:), iProcPE_A(:)

    integer :: iError, iBlock, iBlockPE

    real, allocatable :: CoLat_I(:), Lon_I(:), Alt_I(:)
    real, allocatable :: LatPE_I(:), LonPE_I(:)

    logical :: DoTest, DoTestMe, Done

    !------------------------------------------------------
    !    call CON_set_do_test(NameSub,DoTest, DoTestMe)
    !    if(DoTest)write(*,*)NameSub,' IsInitialized=',IsInitialized
    if(IsInitialized) return

    IsInitialized=.true.

    if(iProc>=0)then
       allocate(CoLat_I(nBlocksLat*nLats), &
            Lon_I(nBlocksLon*nLons), &
            iProc_A(nBlocksLat*nBlocksLon), &
            LatPE_I(nBlocksLat*nBlocksLon*nLats), &
            LonPE_I(nBlocksLat*nBlocksLon*nLons), &
            iProcPE_A(nBlocksLat*nBlocksLon),&
            Alt_I(-1:nAlts+2),&
            stat = iError)

       if (iError /= 0) then
          write(*,*) NameSub, " Error in allocating variables"
          write(*,*) " Lat_I, Lon_I, iProc_A, LatPE_I, LonPE_I, iProcPE_A"
          call CON_stop(NameSub//' UA_ERROR')
       endif

       LatPE_I   = -1.0e32
       LonPE_I   = -1.0e32
       iProcPE_A = -1

       do iBlockPE = 1, nBlocks

          iBlock = iStartBLK + iBlockPE

          !        LatPE_I((iBlock-1)*nLats+1:iBlock*nLats) = &
          !             Latitude(1:nLats,iBlockPE)
          !        LonPE_I((iBlock-1)*nLons+1:iBlock*nLons) = &
          !             Longitude(1:nLons,iBlockPE)
          iProcPE_A(iBlock) = iProc

       enddo

       ! call MPI_allreduce( LatPE_I, CoLat_I, nBlocksLon*nBlocksLat*nLats, &
       !          MPI_REAL, MPI_MAX, iComm, iError)
       !     ! Save into colatitudes instead of latitude
       !     CoLat_I = cHalfPi - CoLat_I
       !
       !     call MPI_allreduce( LonPE_I, Lon_I, nBlocksLon*nBlocksLat*nLons, &
       !          MPI_REAL, MPI_MAX, iComm, iError)

       call MPI_allreduce( iProcPE_A, iProc_A, nBlocksLon*nBlocksLat, &
            MPI_INTEGER, MPI_MAX, iComm, iError)
       !     Alt_I=Altitude(:)
    else
       allocate( CoLat_I(1), Lon_I(1),iProc_A(1),Alt_I(1),stat=iError)
       call check_allocate(iError,NameSub)
    end if

    call set_grid_descriptor(                        &
         UA_,                                        &! component index
         nDim=3,                                     &! dimensionality
         nRootBlock_D=(/nBlocksLat,nBlocksLon,1/),     &! blocks
         nCell_D =(/nLats,nLons,nAlts/),             &! size of node based grid
         XyzMin_D=(/cHalf,cHalf,cHalf/),                   &! generalize coord
         XyzMax_D=(/nLats-cHalf,nLons-cHalf,nAlts-cHalf/), &! generalize coord
         TypeCoord='GEO',                            &! magnetic coordinates
         !Coord1_I= CoLat_I,                          &! colatitudes
         !Coord2_I= Lon_I,                            &! longitudes
         !Coord3_I= Alt_I,                            &! radial size in meters
         iProc_A = iProc_A,                          &! processor assigment
         IsPeriodic_D=(/.false.,.true.,.false./))     ! periodic in longitude

  end subroutine UA_set_grid

  !============================================================================

  subroutine UA_init_session(iSession, SWMFTime)

    use CON_physics,    ONLY: get_time
    use ModTime, only : StartTime, iTimeArray, CurrentTime

    real, intent(in)    :: SWMFTime
    integer, intent(in) :: iSession

    logical :: IsFirstTime = .true.

    if (IsFirstTime) then

       ! Set time related variables for UA
       call get_time(tStartOut = StartTime)

       CurrentTime = StartTime + SWMFTime
       call time_real_to_int(StartTime, iTimeArray)

       call fix_vernal_time

       call initialize_gitm(CurrentTime)
       call write_output

       IsFirstTime = .false.

    endif

  end subroutine UA_init_session

  !==========================================================================

  subroutine UA_run(SWMFTime, SWMFTimeLimit)

    use ModGITM
    use ModTime
    use ModInputs, only: iDebugLevel, Is1D
    use ModTimeConvert, ONLY: time_real_to_int, n_day_of_year

    save

    real, intent(in)    :: SWMFTimeLimit
    real, intent(inout) :: SWMFTime

    integer :: index,lat,long,s,i,j,a, k, jj
    logical :: done, status_ok
    integer :: ierror, CLAWiter, n
    real    :: maxi,tt
    integer :: time_array(7)

    logical :: exist, IsDone

    CurrentTime = StartTime + SWMFTime
    EndTime     = StartTime + SWMFTimeLimit

    if (iDebugLevel > 1) then
       call time_real_to_int(CurrentTime, time_array)
       write(*,"(a,i5,5i3,i3)") "> Running UA from time : ",time_array(1:7)
    endif

    call calc_pressure

    Dt = 1.e32

    call calc_timestep_vertical
    if (.not. Is1D) call calc_timestep_horizontal

    if (iDebugLevel > 1) write(*,"(a,f13.5)") "> UA_run Dt : ",Dt

    call advance

    iStep = iStep + 1

    call write_output

    SWMFTime = CurrentTime - StartTime

  end subroutine UA_run

  !==========================================================================

  subroutine UA_save_restart(TimeSimulation)

    use ModInputs

    real, intent(in) :: TimeSimulation

    call write_restart("UA/restartOUT/")

  end subroutine UA_save_restart

  !==========================================================================

  subroutine UA_finalize(TimeSimulation)

    real, intent(in) :: TimeSimulation

    call finalize_gitm

  end subroutine UA_finalize

end module UA_wrapper
