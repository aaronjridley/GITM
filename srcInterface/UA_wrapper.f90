!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module UA_wrapper

  ! Wrapper for GITM Upper Atmosphere (UA) component

  use ModUtilities, ONLY: CON_stop, CON_set_do_test

  implicit none
  save

  private ! except

  public:: UA_set_param
  public:: UA_init_session
  public:: UA_run
  public:: UA_save_restart
  public:: UA_finalize

  ! IE Coupler:
  public :: UA_get_info_for_ie
  public :: UA_get_for_ie
  public :: UA_put_from_ie

  ! GM Coupler (not implemented in this version of GITM):
  public :: UA_get_for_gm
  public :: UA_get_grid_info
  public :: UA_find_points

  ! Variables for UA-IE coupling:
  integer, save :: iSizeIeHemi, jSizeIeHemi ! Size of IE grid for 1 hemisphere.

  ! CON_coupler_points
  ! These are not called anywhere, should not be public.
  !public:: UA_find_points
  !public:: UA_get_grid_info

  ! GM coupler
  ! These are disabled for the time being:
  !public:: UA_get_for_gm
  !public:: UA_put_from_gm

  ! Are not called; should not be public:
  !public:: UA_put_from_gm_dt
  !public:: UA_put_from_gm_init

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
          call check_dir("UA/restartOUT")
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
    call CON_set_do_test(NameSub,DoTest, DoTestMe)
    if(DoTest) write(*,*)NameSub,' IsInitialized=',IsInitialized

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

                  !LatPE_I((iBlock-1)*nLats+1:iBlock*nLats) = &
                  !     Latitude(1:nLats,iBlockPE)
                  !LonPE_I((iBlock-1)*nLons+1:iBlock*nLons) = &
                  !     Longitude(1:nLons,iBlockPE)
                  ! write(*,*) 'iProc: ',iProc
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

    if(DoTest) write(*,*) NameSub//': nCell: ',nLats,nLons,nAlts

    call set_grid_descriptor(                        &
         UA_,                                        &! component index
         nDim=3,                                     &! dimensionality
         nRootBlock_D=(/nBlocksLat,nBlocksLon,1/),   &! blocks
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
    use ModTimeConvert, ONLY: TimeType
    use CON_physics,    ONLY: get_time
    use ModTime,        ONLY: StartTime, EndTime, iTimeArray, CurrentTime

    real, intent(in)    :: SWMFTime
    integer, intent(in) :: iSession

    type(TimeType) :: TimeSwmfEnd
    logical :: IsFirstTime = .true.

    ! Debug variables:
    integer, dimension(7) :: iTimeDebug
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'UA_init_session'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(DoTestMe) write(*,*) &
         NameSub//' called; IsFirstTime = ', IsFirstTime

    if (IsFirstTime) then

       ! Set time related variables for UA
       call get_time(tStartOut  = StartTime)
       call get_time(TimeEndOut = TimeSwmfEnd)

       ! Get times as real numbers:
       EndTime = TimeSwmfEnd % Time  ! End time as a real.
       CurrentTime = StartTime + SWMFTime
       call time_real_to_int(StartTime, iTimeArray) ! get time as integers

       if(DoTestMe) then
          write(*,*) NameSub//' Timing for UA (floating point):'
          write(*,*) NameSub//' Start time   = ', StartTime
          write(*,*) NameSub//' End time     = ', EndTime
          write(*,*) NameSub//' Current time = ', CurrentTime

          write(*,*) NameSub//' Timing for UA (integer datetime):'
          write(*,'(a, 7i5)') NameSub//' Start time   = ', iTimeArray
          call time_real_to_int(EndTime, iTimeDebug)
          write(*,'(a, 7i5)') NameSub//' End time     = ', iTimeDebug
          call time_real_to_int(CurrentTime, iTimeDebug)
          write(*,'(a, 7i5)') NameSub//' Current time = ', iTimeDebug
       end if

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

    character(len=*), parameter :: NameSub='UA_save_restart'
    !--------------------------------------------------------------------------
    call write_restart("UA/restartOUT/")

  end subroutine UA_save_restart

  !==========================================================================

  subroutine UA_finalize(TimeSimulation)

    real, intent(in) :: TimeSimulation

    character(len=*), parameter :: NameSub='UA_finalize'
    !--------------------------------------------------------------------------
    call finalize_gitm

  end subroutine UA_finalize

  !============================================================================
  subroutine UA_get_info_for_ie(nVar, NameVar_V, nMagLat, nMagLon)
    ! Get number and names of variables for IE to UA coupling.
    ! UA reports what variables it needs here.
    ! IE will use this info to fill buffers appropriately.

    use ModElectrodynamics, ONLY: MagLatRes, MagLonRes

    integer, intent(out) :: nVar
    integer, intent(out), optional :: nMagLat, nMagLon
    character(len=*), intent(out), optional :: NameVar_V(:)

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='UA_get_info_for_ie'
    !-------------------------------------------------------------------------

    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Right now, only 3 variables are needed for UA-IE coupling:
    ! potential, average energy, and energy flux.
    nVar = 3
    if(present(NameVar_V)) NameVar_V(1:3) = (/'pot','ave','tot'/)

    ! If more information is needed, PARAMs should be set to configure this
    ! behavior and changes made here, e.g.,
    ! if(DoCoupleThing) nVar=nVar+1
    ! if(present(NameVar_I)) NameVar_I(nVar) = 'new'

    ! Pass information on size of GITM's electrodynamic grid.
    ! Set following src/calc_electrodynamcs.f90::calc_electrodynamics
    ! Lines ~101-102 and 341-342 (approx.; see source for more details).
    if(present(nMagLat)) nMagLat =  88.0/MagLatRes ! 1 hemi ONLY, no equator.
    if(present(nMagLon)) nMagLon = 360.0/MagLonRes+1

    if(DoTestMe)then
       write(*,*) NameSub//': nVar=', nVar
       if(present(NameVar_V)) write(*,*) NameSub//': NameVar_V=',NameVar_V
    end if

  end subroutine UA_get_info_for_ie

  !============================================================================
  subroutine UA_put_from_ie(Buffer_IIV, iSizeIn, jSizeIn, nVarIn, &
       NameVarIn_V, iBlock)

    use ModNumConst, only:cPi
    use CON_coupler, ONLY: Grid_C, IE_, ncell_id
    use ModGITM, ONLY: iProcGITM=>iProc
    use ModEIE_Interface

  ! This gets called for each variable- external loop over all variable names.
  ! Namevar = Pot, Ave, and Tot.
  ! iBlock = 1:North, 2:South.

  ! Variables and what they do:
  ! EIEi_HavenLats = nLatsIE
  ! EIEi_HavenMlt  = nMltIE

  implicit none

  integer, intent(in)           :: iSizeIn, jSizeIn, nVarIn
  real, intent(in)              :: Buffer_IIV(iSizeIn,jSizeIn,nVarIn)
  character (len=*),intent(in)  :: NameVarIn_V(nVarIn)
  integer,intent(in)            :: iBlock

  logical :: IsInitialized = .false.

  integer :: i,j,ii, iVar, iError, nCells_D(2)

  logical :: DoTest, DoTestMe
  character (len=*), parameter :: NameSub='UA_put_from_ie'
  !--------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest,DoTestMe)

  if(.not. IsInitialized)then
     ! Gather information on size/shape of IE grid.
     ! Note that this is built for Ridley Serial, which stores its grid
     ! using "iSize - 1".  The +1 here compensates for that.
     nCells_D = ncell_id(IE_)
     iSizeIeHemi = nCells_D(1) + 1 ! Number of IE colats per hemisphere
     jSizeIeHemi = nCells_D(2) + 1 ! Number of IE mlts

     ! Build IE grid within UA infrastructure:
     call EIE_InitGrid(iSizeIeHemi, jSizeIeHemi, 2, iError)
     call EIE_FillLats(90.0-(Grid_C(IE_) % Coord1_I)*180.0/cPi,iError)
     call EIE_FillMltsOffset((Grid_C(IE_) % Coord2_I)*180.0/cPi,iError)
     IsInitialized = .true.

     ! Initial values for coupling:
     EIEr3_HavePotential = 0.0
     EIEr3_HaveAveE      = 0.0
     EIEr3_HaveEFlux     = 0.0
  end if

  ! Debug: print basic coupling info:
  if(DoTest .and. (iProcGITM==0)) then
     write(*,*)NameSub//' coupling info:'
     write(*,*) 'UA/GITM: IE grid set to nLats x nLons = ', &
          iSizeIeHemi, jSizeIeHemi
     write(*,'(a, 2i5)')' UA EIE nLats, nLons = ', EIEi_HavenLats, EIEi_HavenMlts
     write(*,*)'Buffer shape information:'
     write(*,'(a, 3i5)') ' iSizeIn, jSizeIn, nVarIn = ', iSizeIn, jSizeIn, nVarIn
     write(*,*)'Shape of Buffer_IIV = ', shape(Buffer_IIV)
     write(*,*)'ncell_id(IE_) = ', ncell_id(IE_)
     write(*,*)'iSizeIeHemi, jSizeIeHemi = ', iSizeIeHemi, jSizeIeHemi
  end if

  ! Debug statement: print max/mins of transferred variables.
  if(DoTest .and. (iProcGITM==0)) write(*,*) &
       'UA WRAPPER: Max/min of received variables (iBlock=',iBlock,'):'

  ! Put each variable where it belongs based on the name.x
  do iVar=1, nVarIn
     select case(NameVarIn_V(iVar))

     ! Electric Potential:
     case('pot')
        do i = 1, EIEi_HavenLats
           ii = i
           ! Flip southern hemisphere
           if (iBlock == 2) ii = EIEi_HavenLats - i + 1
           do j = 1, EIEi_HavenMlts
              EIEr3_HavePotential(j,i,iBlock) = Buffer_IIV(ii,j,iVar)
           enddo
        enddo
        if(DoTest .and. (iProcGITM==0)) write(*,'(a, 2(e12.5,1x))') 'UA pot: ', &
             maxval(EIEr3_HavePotential(:,:,iBlock)),   &
             minval(EIEr3_HavePotential(:,:,iBlock))

     ! Precipitation/average energy:
     case('ave')
        do i = 1, EIEi_HavenLats
           ii = i
           ! Flip southern hemisphere
           if (iBlock == 2) ii = EIEi_HavenLats - i + 1
           do j = 1, EIEi_HavenMlts
              EIEr3_HaveAveE(j,i,iBlock)  = Buffer_IIV(ii,j, iVar)
           enddo
        enddo
        if(DoTest .and. (iProcGITM==0)) write(*,'(a, 2(e12.5,1x))') 'UA ave: ', &
             maxval(EIEr3_HaveAveE(:,:,iBlock)),   &
             minval(EIEr3_HaveAveE(:,:,iBlock))

     ! Precipitation/Total energy flux:
     case('tot')
        do i = 1, EIEi_HavenLats
           ii = i
           ! Flip southern hemisphere
           if (iBlock == 2) ii = EIEi_HavenLats - i + 1
           do j = 1, EIEi_HavenMlts
              EIEr3_HaveEFlux(j,i,iBlock) = &
                   Buffer_IIV(ii,j,iVar) / (1.0e-7 * 100.0 * 100.0)
           enddo
        enddo
        if(DoTest .and. (iProcGITM==0)) write(*,'(a, 2(e12.5,1x))') 'UA tot: ', &
             maxval(EIEr3_HaveEFlux(:,:,iBlock)),   &
             minval(EIEr3_HaveEFlux(:,:,iBlock))

     case default
        call CON_stop(NameSub//' invalid NameVarIn='//NameVarIn_V(iVar))

     end select
  end do

  ! This bypasses some internal GITM issues.
  UAl_UseGridBasedEIE = .true.

end subroutine UA_put_from_ie

!============================================================================
subroutine UA_get_for_ie(BufferOut_IIBV, nMltIn, nLatIn, nVarIn, NameVarIn_V)

  use ModTime, ONLY: CurrentTime
  use ModPlotFile, ONLY: save_plot_file

  ! Input variables: size of grid/number and name of requested vars
  integer,          intent(in) :: nMltIn, nLatIn, nVarIn
  character(len=3), intent(in) :: NameVarIn_V(nVarIn)

  ! Buffer for the variables on the 2D IE grid
  real, intent(out) :: BufferOut_IIBV(nMltIn, nLatIn, 2, nVarIn) ! to fill

  ! Buffers for UA variables after they are calculated:
  real, dimension(:,:), allocatable :: UAr2_Mlts, UAr2_Lats
  real, dimension(:,:), allocatable :: UAr2_Hal, UAr2_Ped, UAr2_Fac

  integer :: UAi_nLats, UAi_nMlts, iError, iBlock, &
       iStartN, iEndN, iStartS, iEndS, iVar

  ! Debug related variables:
  character(len=35) :: NameFile
  character(len=3) :: NameVarSave_V(nVarIn)
  real :: BufferPlot_VII(nVarIn, nMltIn, nLatIn)
  integer :: time_array(7)
  logical :: DoTest, DoTestMe
  character (len=*), parameter :: NameSub='UA_get_for_ie'
  !-------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest,DoTestMe)
  ! Do we need to do this?  Uncomment or delete when we figure that out.
  !call initialize_ie_ua_buffers(iError)

  ! Calculate electrodynamics; get size of UA grid.
  call UA_calc_electrodynamics(UAi_nMlts, UAi_nLats)

  ! Ensure size of grid is same as buffer size for one hemisphere:
  if( (nMltIn/=UAi_nMlts).or.(2*nLatIn+1/=UAi_nLats) ) then
     write(*,*) NameSub//' UA electromagnetic grid does not match buffer:'
     write(*,*) NameSub//' UA E&M grid nLon, nLats = ', UAi_nMlts, UAi_nLats
     write(*,*) NameSub//' UA-IE coupling buffer nLon, nLats = ', nMltIn,nLatIn
     call CON_stop(NameSub//' Array size mismatch.')
  end if

  ! Collect relevant values using native GITM functions:
  allocate(UAr2_Fac(UAi_nMlts, UAi_nLats), &
       UAr2_Ped(UAi_nMlts, UAi_nLats), &
       UAr2_Hal(UAi_nMlts, UAi_nLats), &
       UAr2_Lats(UAi_nMlts, UAi_nLats), &
       UAr2_Mlts(UAi_nMlts, UAi_nLats), stat=iError)

  call UA_fill_electrodynamics(UAr2_Fac,UAr2_Ped,UAr2_Hal,UAr2_Lats, UAr2_Mlts)

  if(DoTest)write(*,*)NameSub//' nMlts,nLats= ',UAi_nMlts, UAi_nLats

  ! Set indices for northern/southern hemispheres.
  ! UA goes from the South pole to the north pole, while IE goes
  ! from the north pole to the south pole, so the blocks have to
  ! be reversed, basically.  Equatorial point is skipped.
  ! Northern hemisphere:
  iStartN = INT(CEILING(REAL(UAi_nLats)/2)) + 1 !rounding an odd number up
  iEndN   = UAi_nLats
  ! Southern Hemisphere:
  iStartS = 1
  iEndS   = UAi_nLats/2

  ! Only collect variables requested:
  do iVar=1, nVarIn
     select case (NameVarIn_V(iVar))
     case('lon')
        BufferOut_IIBV(:,:,1,iVar) = UAr2_Mlts(:,iStartN:iEndN)
        BufferOut_IIBV(:,:,2,iVar) = UAr2_Mlts(:,iStartS:iEndS)
     case('lat')
        BufferOut_IIBV(:,:,1,iVar) = UAr2_Lats(:,iStartN:iEndN)
        BufferOut_IIBV(:,:,2,iVar) = UAr2_Lats(:,iStartS:iEndS)
     case('hal')
        BufferOut_IIBV(:,:,1,iVar) = UAr2_Hal(:,iStartN:iEndN)
        BufferOut_IIBV(:,:,2,iVar) = UAr2_Hal(:,iStartS:iEndS)
     case('ped')
        BufferOut_IIBV(:,:,1,iVar) = UAr2_Ped(:,iStartN:iEndN)
        BufferOut_IIBV(:,:,2,iVar) = UAr2_Ped(:,iStartS:iEndS)
     case('fac')
        BufferOut_IIBV(:,:,1,iVar) = UAr2_Fac(:,iStartN:iEndN)
        BufferOut_IIBV(:,:,2,iVar) = UAr2_Fac(:,iStartS:iEndS)
     case default
        call CON_stop(NameSub//' Unrecognized coupling variable: '// &
             NameVarIn_V(iVar))
     end select

  end do

  ! Dump contents to file in debug mode:
  if(DoTestMe)then
     ! Get current time; create output file name:
     call time_real_to_int(CurrentTime, time_array)
     write(NameFile,'(a,i4.3,2i2.2,"_",3i2.2,a)') &
          'ua_ie_buffer_t', time_array(1:6), '.out'
     write(*,*) NameSub//': Saving buffer to ', NameFile

     ! Fill output array to match required ordering
     do iVar=1, nVarIn
        BufferPlot_VII(iVar,:,:) = BufferOut_IIBV(:,:,1,iVar)
     end do

     ! Create list of var names that includes dimensions:
     NameVarSave_V(1:2) = (/'mlt', 'lat'/)
     NameVarSave_V(3:nVarIn) = NameVarIn_V(3:nVarIn)

     ! Write file:
     call save_plot_file('UA/data/' //NameFile, &
          TypeFileIn='ascii', &
          StringHeaderIn='UA_get_for_ie BufferOut contents', &
          TimeIn=CurrentTime, &
          NameVarIn_I = NameVarSave_V, &
          NameUnitsIn = 'hours degrees Seimens Seimens', &
          nDimIn = 2, &
          Coord1In_I = BufferPlot_VII(1, :, 1), &
          Coord2In_I = BufferPlot_VII(2, 1, :), &
          VarIn_VII = BufferPlot_VII(3:nVarIn, :, :))
  end if

  ! Deallocate intermediate variables:
  deallocate(UAr2_Fac, UAr2_Ped, UAr2_Hal, UAr2_Lats, UAr2_Mlts)

end subroutine UA_get_for_ie

!============================================================================
subroutine UA_get_for_gm(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
  Data_VI)

  ! This version of GITM does not support coupling directly to GM.

  logical,          intent(in) :: IsNew   ! true for new point array
  character(len=*), intent(in) :: NameVar ! List of variables
  integer,          intent(in) :: nVarIn  ! Number of variables in Data_VI
  integer,          intent(in) :: nDimIn  ! Dimensionality of positions
  integer,          intent(in) :: nPoint  ! Number of points in Xyz_DI

  real, intent(in) :: Xyz_DI(nDimIn, nPoint)  ! Position vectors
  real, intent(out):: Data_VI(nVarIn, nPoint) ! Data array

  character(len=*), parameter :: NameSub='UA_get_for_gm'

  call CON_stop(NameSub// &
    ': UA_ERROR: This version of GITM does not support coupling to GM.')

end subroutine UA_get_for_gm
!============================================================================

subroutine UA_get_grid_info(nDimOut, iGridOut, iDecompOut)
   ! Not implemented in this version of GITM.

   integer, intent(out):: nDimOut    ! grid dimensionality
   integer, intent(out):: iGridOut   ! grid index
   integer, intent(out):: iDecompOut ! decomposition index

   character(len=*), parameter :: NameSub = 'UA_get_grid_info'
   call CON_stop(NameSub// &
     ': UA_ERROR: This version of GITM does not support coupling to GM.')

 end subroutine UA_get_grid_info
 !============================================================================!============================================================================

subroutine UA_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)
   ! Not implemented in this version of GITM.

   integer, intent(in) :: nDimIn                ! dimension of positions
   integer, intent(in) :: nPoint                ! number of positions
   real,    intent(in) :: Xyz_DI(nDimIn,nPoint) ! positions
   integer, intent(out):: iProc_I(nPoint)       ! processor owning position

   character(len=*), parameter:: NameSub = 'UA_find_points'
   call CON_stop(NameSub// &
     ': UA_ERROR: This version of GITM does not support coupling to GM.')

 end subroutine UA_find_points
 !============================================================================
end module UA_wrapper
