
module ModCoupSAMI3

  ! variables gathered in GITM
  real, allocatable:: GatLongitude(:,:)
  real, allocatable:: GatLatitude(:,:)
  real, allocatable:: GatAltitude(:,:) 
  real,allocatable :: GatGITMVars(:,:,:,:,:)

  ! GITM coupling flags
  logical :: StartCoup = .false.
  logical :: StartCoup_pot = .false.

  !logical :: CorotationAddedC = .false.

  real :: DtCouple = 300.0

  integer :: iProcGlobal

contains 

  subroutine gather_data_gitm

    use ModGITM
    use ModMpi
    use ModPlanet
    use coup_mod, only:nParams
    implicit none

    real,allocatable :: GITMLons(:)
    real,allocatable :: GITMLats(:)
    real,allocatable :: GITMAlts(:)
    real,allocatable :: GITMVars(:,:,:,:)

    real,allocatable :: Longitude0(:,:)

    integer :: iLon,iBlock,iLat,iAlt    
    integer :: iError
    integer :: loc(7)
    integer :: iParam,ip
    integer :: j

    ! allocate GITM data to exchange
    if (iproc==0) then
       allocate(GatLongitude(-1:(nLons+4)*nBlocks-2,nprocs))
       allocate(GatLatitude(-1:(nLats+4)*nBlocks-2,nprocs))
       allocate(GatAltitude(-1:(nAlts+4)*nBlocks-2,nprocs))

       allocate(GatGITMVars(nParams,&
            -1:(nLons+4)*nBlocks-2,&
            -1:(nLats+4)*nBlocks-2,&
            -1:(nAlts+4)*nBlocks-2,&
            nprocs))

    endif

    !! how to take the coratation into consideration when coupling
    ! This then implies that the grid is in local time, so we need to 
    ! feed SAMI the local time (still 0-360, though).

    !do iBlock = 1,nBlocks
    !   do iLon = -1,nLons+2
    !            Longitude0(iLon,iBlock) = Longitude(iLon,iBlock)

    !            if (CorotationAddedC) then 
    !                Longitude0(iLon,iBlock) = Longitude0(iLon,iBlock)+LocalTime(iLon)*15.0
    !                Longitude0(iLon,iBlock) = mod(Longitude0(iLon,iBlock),360.0)
    !            endif
    !   enddo
    !enddo

    allocate(GITMLons(-1:(nLons+4)*nBlocks-2))
    allocate(GITMLats(-1:(nLats+4)*nBlocks-2))
    allocate(GITMAlts(-1:(nAlts+4)*nBlocks-2))
    allocate(GITMVars(nParams,-1:(nLons+4)*nBlocks-2,&
         -1:(nLats+4)*nBlocks-2,&
         -1:(nAlts+4)*nBlocks-2))

    loc = (/iH_,  &
         iHe_,  &
         iN_4S_,&
         iO_3P_,&
         iN2_  ,&
         iO2_  ,&
         iNO_  /)

    do iBlock=1,nBlocks
       do iParam = 1,7
          GITMVars(iParam,&
               (iBlock-1)*(nLons+4)-1:iBlock*(nLons+4)-2,&
               (iBlock-1)*(nLats+4)-1:iBlock*(nLats+4)-2,&
               (iBlock-1)*(nAlts+4)-1:iBlock*(nAlts+4)-2) = &
               NDensityS(:,:,:,loc(iParam),iBlock)
       enddo

       GITMVars(8,&
            (iBlock-1)*(nLons+4)-1:iBlock*(nLons+4)-2,&
            (iBlock-1)*(nLats+4)-1:iBlock*(nLats+4)-2,&
            (iBlock-1)*(nAlts+4)-1:iBlock*(nAlts+4)-2) = &
            Temperature(:,:,:,iBlock)*TempUnit(:,:,:)

       GITMVars(9,&
            (iBlock-1)*(nLons+4)-1:iBlock*(nLons+4)-2,&
            (iBlock-1)*(nLats+4)-1:iBlock*(nLats+4)-2,&
            (iBlock-1)*(nAlts+4)-1:iBlock*(nAlts+4)-2) = &
            Velocity(:,:,:,iEast_,iBlock)
       GITMVars(10,&
            (iBlock-1)*(nLons+4)-1:iBlock*(nLons+4)-2,&
            (iBlock-1)*(nLats+4)-1:iBlock*(nLats+4)-2,&
            (iBlock-1)*(nAlts+4)-1:iBlock*(nAlts+4)-2) = &
            Velocity(:,:,:,iNorth_,iBlock)

       GITMLons((iBlock-1)*(nLons+4)-1:iBlock*(nLons+4)-2) = &
            Longitude(:,iBlock)

       GITMLats((iBlock-1)*(nLats+4)-1:iBlock*(nLats+4)-2) = &
            Latitude(:,iBlock)

       GITMAlts((iBlock-1)*(nAlts+4)-1:iBlock*(nAlts+4)-2) = &
            Altitude_GB(1,1,:,iBlock)

       do iLon = -1, nLons+2
          GITMVars(11,(iBlock-1)*(nLons+4)+iLon,:,:) = &
               Longitude(iLon,iBlock)
       enddo

       do iLat = -1, nLats+2
          GITMVars(12,:,(iBlock-1)*(nLats+4)+iLat,:) = &
               Latitude(iLat,iBlock)
       enddo
       do iAlt = -1, nAlts+2
          GITMVars(13,:,:,(iBlock-1)*(nAlts+4)+iAlt) = &
               Altitude_GB(1,1,iAlt,iBlock)
       enddo
    enddo

    call MPI_Gather(GITMLons,&
         (nLons+4)*nBlocks,&
         MPI_REAL,&
         GatLongitude,&
         (nLons+4)*nBlocks,&
         MPI_REAL,&
         0,iCommGITM,iError)

    call MPI_Gather(GITMLats,&
         (nLats+4)*nBlocks,&
         MPI_REAL,&
         GatLatitude,&
         (nLats+4)*nBlocks,&
         MPI_REAL,&
         0,iCommGITM,iError)

    call MPI_Gather(GITMAlts,&
         (nAlts+4)*nBlocks,&
         MPI_REAL,&
         GatAltitude,&
         (nAlts+4)*nBlocks,&
         MPI_REAL,&
         0,iCommGITM,iError)

    call MPI_Gather(GITMVars,&
         nParams*(nLons+4)*(nLats+4)*(nAlts+4)*nBlocks,&
         MPI_REAL,&
         GatGITMVars,&
         nParams*(nLons+4)*(nLats+4)*(nAlts+4)*nBlocks,&
         MPI_REAL,&
         0,iCommGITM,iError)

    deallocate(GITMLons,GITMLats,GITMAlts,GITMVars)

    if (iproc ==0) &
         print*,'---- Done with GITM data gathering!!!',iproc

  end subroutine gather_data_gitm

  subroutine overwrite_sami(iBlock)

    use advance_mod
    use ModConstants,only: pi
    use ModSamiInterp
    use ModGITM, only: MLongitude,MLatitude,Altitude_GB,&
         nLons,nLats,nAlts,&
         IDensityS,iHP_,ie_,iNOP_,iO_4SP_,iO2P_,iHeP_,&
         eTemperature,iTemperature,iproc,iCommGITM,nIons,&
         LocalTime

    implicit none

    save

    integer, intent(in) :: iBlock

    logical :: IsFirstTime = .true.
    real, allocatable :: InterpFactors(:,:,:,:)
    real, allocatable :: InterpIndex(:,:,:,:)

    integer :: iLon, iLat, iAlt,iIon
    real :: xd00,xd01,xd10,xd11,yd0,yd1,zd
    real :: data_inp(10),data_inp2(10)
    real :: alttmp
    integer :: iError
    integer :: locgitm(5)
    real :: lont,latt,altt

    call report("ionosphere_overwrite_sami",2)

    if (iproc == 0) &
         print*,'---- GITM Start Couping init!!!',IsFirstTime,iproc

    if (IsFirstTime) then

       allocate(InterpFactors(7,-1:nLons+2,-1:nLats+2,1:nAlts))
       allocate(InterpIndex(8,-1:nLons+2,-1:nLats+2,1:nAlts))

       do iLon = -1,nLons+2
          do iLat = -1,nLats+2
             do iAlt = 1, nAlts

                gitm_mlon = MLongitude(iLon,iLat,iAlt,iBlock)
                gitm_mlat = MLatitude(iLon,iLat,iAlt,iBlock)
                gitm_alt = Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0


                if (gitm_mlon<0) gitm_mlon = gitm_mlon+360.0

                alttmp = gitm_alt+re

                call get_index_8points_3

                InterpIndex(1 ,iLon,iLat,iAlt) = i0 
                InterpIndex(2 ,iLon,iLat,iAlt) = j0
                InterpIndex(3 ,iLon,iLat,iAlt) = k0
                InterpIndex(4 ,iLon,iLat,iAlt) = l0
                InterpIndex(5 ,iLon,iLat,iAlt) = i1
                InterpIndex(6 ,iLon,iLat,iAlt) = j1
                InterpIndex(7 ,iLon,iLat,iAlt) = k1
                InterpIndex(8 ,iLon,iLat,iAlt) = l1

                if ((j0>0).and.(k0>0).and.(l0>0)) then

                   call get_factors_from_8points_sph_3(gitm_mlon,gitm_mlat,alttmp,&
                        xd00,xd01,xd10,xd11,yd0,yd1,zd,lont,latt,altt)

                   if ((xd00*xd01*xd10*xd11*yd0*yd1*zd)<0) then
                      print*,'===>> Dist',iLon,iLat,iAlt,xd00,xd01,xd10,xd11,yd0,yd1,zd
                      stop
                   endif
                   InterpFactors(1,iLon,iLat,iAlt) = xd00
                   InterpFactors(2,iLon,iLat,iAlt) = xd01
                   InterpFactors(3,iLon,iLat,iAlt) = xd10
                   InterpFactors(4,iLon,iLat,iAlt) = xd11
                   InterpFactors(5,iLon,iLat,iAlt) = yd0
                   InterpFactors(6,iLon,iLat,iAlt) = yd1
                   InterpFactors(7,iLon,iLat,iAlt) = zd

                endif

             enddo
          enddo
       enddo

       IsFirstTime = .false.

    endif

    if (iproc == 0 ) &
         print*,'---- GITM Couping init Done!!!',iproc
!!! how to take the coratation into consideration when coupling
    !if (CorotationAdded) then 
    !   ! This then implies that the grid is in local time, so we need to 
    !   ! feed SAMI the local time (still 0-360, though).

    !   do iLon = -1,nLons+2
    !      do iLat = -1, nLats+2
    !         do iAlt = 1, nAlts
    !            GitmLons(iPoint) = LocalTime(iLon)*15.0
    !            if (iAlt == 1 .and. iLat == 1) write(*,*) iLon, LocalTime(iLon), GitmLons(iPoint) 
    !            iPoint = iPoint + 1
    !         enddo
    !      enddo
    !   enddo

    !endif

    locgitm = (/iHP_,iO_4SP_,iNOP_,iO2P_,iHeP_/)

    do iLon = -1,nLons+2
       do iLat = -1, nLats+2
          do iAlt = 1, nAlts

             !gitm_mlon = Longitude(iLon,iBlock)/pi*180.0
             !gitm_mlat = Latitude(iLat,iBlock)/pi*180.0
             !gitm_alt = Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0

             gitm_mlon = MLongitude(iLon,iLat,iAlt,iBlock)
             gitm_mlat = MLatitude(iLon,iLat,iAlt,iBlock)
             gitm_alt = Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0

             !print*,'===>> Magnetic Lon/Lat b',gitm_mlon,gitm_mlat,gitm_alt
             if (gitm_mlon<0) gitm_mlon = gitm_mlon+360.0

             !if (CorotationAddedC) then 
             !
             !    gitm_mlon =  gitm_mlon + LocalTime(iLon)*15.0
             !    gitm_mlon = mod(gitm_mlon,360.0)
             !endif

             alttmp = gitm_alt+re

             i0 = InterpIndex(1,iLon,iLat,iAlt)
             j0 = InterpIndex(2,iLon,iLat,iAlt)
             k0 = InterpIndex(3,iLon,iLat,iAlt)
             l0 = InterpIndex(4 ,iLon,iLat,iAlt) 

             if ((j0>0).and.(k0>0).and.(l0>0)) then

                i1 = InterpIndex(5 ,iLon,iLat,iAlt) 
                j1 = InterpIndex(6 ,iLon,iLat,iAlt) 
                k1 = InterpIndex(7 ,iLon,iLat,iAlt) 
                l1 = InterpIndex(8 ,iLon,iLat,iAlt) 

                call interp_with_8points(InterpFactors(1,iLon,iLat,iAlt), &
                     InterpFactors(2,iLon,iLat,iAlt),InterpFactors(3,iLon,iLat,iAlt),&
                     InterpFactors(4,iLon,iLat,iAlt),InterpFactors(5,iLon,iLat,iAlt),&
                     InterpFactors(6,iLon,iLat,iAlt),InterpFactors(7,iLon,iLat,iAlt),&
                     gitm_mlon,gitm_mlat,alttmp,data_inp)

                do iIon = 1, 5
                   if (data_inp(iIon+3)<=0) then

                      print*,'===>>2 Magnetic Lon/Lat ', iproc, gitm_mlon,gitm_mlat,gitm_alt
                      print*,'===>>2 The interpolated value', iproc, &
                           data_inp(1:2),data_inp(3)-re,data_inp(iIon+3),iIon+3

                      print*,'-->---loggitm,iIon ',locgitm(iIon),iIon
                      print*,'--> ',iproc,SAMIVars_g(1:2,k0,j0,i0),SAMIVars_g(3,k0,j0,i0)-re,SAMIVars_g(iIon+3,k0,j0,i0)
                      print*,'--> '
                      print*,'-->' &
                           ,SAMIVars_g(1:2,k1,j0,i0),SAMIVars_g(3,k1,j0,i0)-re,SAMIVars_g(iIon+3,k1,j0,i0)
                      print*,'--> '
                      print*,'--> ' &
                           ,SAMIVars_g(1:2,l0,j1,i0),SAMIVars_g(3,l0,j1,i0)-re,SAMIVars_g(iIon+3,l0,j1,i0)
                      print*,'--> '
                      print*,'--> ', &
                           SAMIVars_g(1:2,l1,j1,i0),SAMIVars_g(3,l1,j1,i0)-re,SAMIVars_g(iIon+3,l1,j1,i0)
                      print*,'--> '
                      print*,'--> ', &
                           SAMIVars_g(1:2,k0,j0,i1),SAMIVars_g(3,k0,j0,i1)-re,SAMIVars_g(iIon+3,k0,j0,i1)
                      print*,'--> '
                      print*,'--> ', &
                           SAMIVars_g(1:2,k1,j0,i1),SAMIVars_g(3,k1,j0,i1)-re,SAMIVars_g(iIon+3,k1,j0,i1)
                      print*,'--> '
                      print*,'--> ', &
                           SAMIVars_g(1:2,l0,j1,i1),SAMIVars_g(3,l0,j1,i1)-re,SAMIVars_g(iIon+3,l0,j1,i1)
                      print*,'--> '
                      print*,'--> ', &
                           SAMIVars_g(1:2,l1,j1,i1),SAMIVars_g(3,l1,j1,i1)-re,SAMIVars_g(iIon+3,l1,j1,i1)

                   else

                      IDensityS(iLon,iLat,iAlt,locgitm(iIon),iBlock) = data_inp(3+iIon)*1e6
                   endif
                enddo


                IDensityS(iLon,iLat,iAlt,ie_,iBlock) = &
                     sum(IDensityS(iLon,iLat,iAlt,1:nIons-1,iBlock))

                if (data_inp(9)<=0) then
                   print*,'===>>2 Magnetic Lon/Lat ', iproc, gitm_mlon,gitm_mlat,gitm_alt
                   print*,'===>>2 The interpolated value', iproc, &
                        data_inp(1:2),data_inp(3)-re,data_inp(9),9
                else
                   iTemperature(iLon,iLat,iAlt,iBlock) = data_inp(9)
                endif
                if (data_inp(10)<=0) then
                   print*,'===>>4 Magnetic Lon/Lat ', iproc, gitm_mlon,gitm_mlat,gitm_alt
                   print*,'===>>4 The interpolated value', iproc, &
                        data_inp(1:2),data_inp(3)-re,data_inp(10),10
                else
                   eTemperature(iLon,iLat,iAlt,iBlock) = data_inp(10)
                endif
             endif
          enddo
       enddo
    enddo

    if (iproc == 0) & 
         print*,'GITM overwrite Done!',iproc
  end subroutine overwrite_sami


  subroutine ExchangeData

    use ModGITM
    use ModMpi
    use coup_mod 
    ! use message_passing_mod
    use parameter_mod, only: nz,nf,nlt

    use ModSamiInterp,only: SAMIVars_g,SAMIPhi_g, &
         blonst_g,blatst_g,baltst_g
    use ModConstants,only : pi
    use ModInputs, only: nBlocksLon,nBlocksLat
    use ModGITM, only: SamiMaster

    implicit none

    integer :: intercomm1,iError,irank
    integer :: intercomm2
    integer :: tag_inter,tag_SandR,tag_SandR1,tag_SandR2
    integer :: status(MPI_STATUS_SIZE)
    real    :: testt =0.0

    logical :: flagsr
    integer :: i,j
    integer :: iLon,iLat,iAlt,ip,iParam
    real, allocatable :: GatLongitude_g0(:),GatLatitude_g0(:),&
         GatAltitude_g0(:),GatGITMVars_g0(:,:,:,:)
    real, allocatable :: GatGITMVar(:,:,:)
    integer :: iBlockLon,iBlockLat,iBlocksAlt
    integer :: nPointsTotal
    integer :: iiLon,iiLat,iiAlt,nLonsTotal,nLatsTotal,nAltsTotal
    logical :: IsFirstTime = .true.
    logical :: IsValidPoint
    real :: tmp(nz,nf,nlt), tmp2d(nf, nlt)
    integer :: iZ

    tag_inter = 11
    tag_SandR = 111
    tag_SandR1 = 222
    tag_SandR2 = 333

    if (iCommGITM /= MPI_COMM_NULL) then

       call MPI_INTERCOMM_CREATE(iCommGITM,0, &
            iCommGlobal,SamiMaster,tag_inter,intercomm1,iError)
       call MPI_COMM_RANK(intercomm1,irank,iError)

       if (irank == 0) then

          !recv from sami3
          call MPI_RECV(SAMIVars_g,&
               nf*nz*nlt*10,&
               MPI_REAL,0,tag_SandR,&
               intercomm1,status,iError)

          call MPI_RECV(SAMIPhi_g,&
               nf*nlt,&
               MPI_REAL,0,tag_SandR,&
               intercomm1,status,iError)

          blonst_g = SAMIVars_g(1,:,:,:)
          blatst_g = SAMIVars_g(2,:,:,:)
          baltst_g = SAMIVars_g(3,:,:,:)

          !print*,'SAMIPhi_g a',SAMIPhi_g(:,1)

       endif

       call report("GITM Broadcasting SAMI Data",0)

       ! write(*,*) "GITM Proc, before bcast : ", irank, iProcGlobal

       do i=1,10

          do iZ = 1, nz

             call MPI_BARRIER(iCommGITM,iError)
             ! if (iProcGlobal == 30) write(*,*) "Looping exchange start", iProcGlobal, i, iZ, nf*nlt
             if (iProc == 0) tmp2d = SAMIVars_g(i,iZ,:,:)
             call mpi_bcast(tmp2d, nf*nlt, MPI_REAL,0,iCommGITM,iError)
             ! if (iProcGlobal == 30) write(*,*) "Looping exchange after bcast", iProcGlobal, i, iZ

             if (iProc > 0) SAMIVars_g(i,iZ,:,:) = tmp2d

             ! write(*,*) "iError : ", iError, iProcGlobal, i 
             call MPI_BARRIER(iCommGITM,iError)

             !if (iProcGlobal == 30) write(*,*) "Looping exchange after barrier", iProcGlobal, i, iZ

          enddo

!          call MPI_BARRIER(iCommGITM,iError)
!          if (iProcGlobal == 30) write(*,*) "Looping exchange start", iProcGlobal, i, nf*nz*nlt
!          if (iProc == 0) tmp = SAMIVars_g(i,:,:,:)
!          call mpi_bcast(tmp, nf*nz*nlt, MPI_REAL,0,iCommGITM,iError)
!          if (iProcGlobal == 30) write(*,*) "Looping exchange after bcast", iProcGlobal, i
!
!          if (iProc > 0) SAMIVars_g(i,:,:,:) = tmp
!
!          write(*,*) "iError : ", iError, iProcGlobal, i 
!          call MPI_BARRIER(iCommGITM,iError)
!
!          if (iProcGlobal == 30) write(*,*) "Looping exchange after barrier", iProcGlobal, i

       enddo

       !write(*,*) "GITM Proc, after bcast : ", irank, iProcGlobal

       blonst_g = SAMIVars_g(1,:,:,:)
       blatst_g = SAMIVars_g(2,:,:,:)
       baltst_g = SAMIVars_g(3,:,:,:)

       call MPI_BARRIER(iCommGITM,iError)
       call mpi_bcast(SAMIPhi_g,nf*nlt,MPI_REAL,0,iCommGITM,iError)
       call MPI_BARRIER(iCommGITM,iError)

       call report("Done Broadcasting SAMI Data",0)

    else if (iCommSAMI0 /= MPI_COMM_NULL) then

       call MPI_INTERCOMM_CREATE(iCommSAMI0,0, &
            iCommGlobal,0,tag_inter,intercomm2,iError)
       call MPI_COMM_RANK(intercomm2,irank,iError)

       if(irank ==0) then

          !!send to GITM
          call MPI_send(SAMIVars,&
               nf*nz*nlt*10,&
               MPI_REAL,0,tag_SandR,&
               intercomm2,iError)

          call MPI_send(SAMIPhi,&
               nf*nlt,&
               MPI_REAL,0,tag_SandR,&
               intercomm2,iError)

          !print*,'SAMIPhi sami',SAMIPhi(:,1)

       endif
    endif


    if (iCommGITM /= MPI_COMM_NULL) then

       call MPI_COMM_RANK(intercomm1,irank,iError)

       if(irank ==0) then
          call MPI_send(nBlocks,1,MPI_INTEGER,&
               0,tag_SandR1,&
               intercomm1,iError)
          call MPI_send(nprocs,1,MPI_INTEGER,&
               0,tag_SandR1,&
               intercomm1,iError)
          call MPI_send(nBlocksLon,1,MPI_INTEGER,&
               0,tag_SandR1,&
               intercomm1,iError)
          call MPI_send(nBlocksLat,1,MPI_INTEGER,&
               0,tag_SandR1,&
               intercomm1,iError)

          nLonsTotal = nBlocksLon*nLons*nBlocks+4       
          nLatsTotal = nBlocksLat*nLats*nBlocks+4       
          nAltsTotal = nAlts+4       

          allocate(GatLongitude_g0(-1:nLonsTotal-2))
          allocate(GatLatitude_g0(-1:nLatsTotal-2))
          allocate(GatAltitude_g0(-1:nAltsTotal-2))
          allocate(GatGITMVars_g0(nParams,&
               -1:nLonsTotal-2,&
               -1:nLatsTotal-2,&
               -1:nAltsTotal-2))

          iBlocksAlt = 1
          do ip = 1, nprocs

             iBlockLat = (ip-1)/nBlocksLon
             iBlockLon = (ip-1)-(iBlockLat*nBlocksLon)

             iBlockLat = iBlockLat + 1
             iBlockLon = iBlockLon + 1

             do iLon = -1, nLons+2
                iiLon = (iBlockLon-1)*(nLons) + iLon
                GatLongitude_g0(iiLon) = GatLongitude(iLon,ip)/pi*180.0

                !if (GatLongitude_g0(iiLon)< 0) then
                !    GatLongitude_g0(iiLon) = GatLongitude_g0(iiLon)+360.0
                !    GatLongitude_g0(iiLon) = mod(GatLongitude_g0(iiLon),360.0)
                !endif
                !if (CorotationAddedC) then 
                !    GatLongitude_g0(iiLon) = GatLongitude_g0(iiLon) + LocalTime(iLon)*15.0
                !    GatLongitude_g0(iiLon) = mod(GatLongitude_g0(iiLon),360.0)
                !endif
                !print*,'-- ip,iLon,iBlockLon,iiLon',ip,iLon,iBlockLon,iiLon

             enddo

             do iLat = -1, nLats+2
                iiLat = (iBlockLat-1)*(nLats) + iLat
                GatLatitude_g0(iiLat) = GatLatitude(iLat,ip)/pi*180.0
             enddo

             do iLon = -1, nLons+2

                iiLon = (iBlockLon-1)*(nLons) + iLon

                IsValidPoint = .true.

                ! If point is ghostcell and not on the edge, it is not a good point
                if (iLon < 1 .and. iBlockLon > 1) IsValidPoint = .false.
                if (iLon > nLons .and. iBlockLon < nBlocksLon) IsValidPoint = .false.

                if (IsValidPoint) then

                   do iLat = -1, nLats+2
                      iiLat = (iBlockLat-1)*(nLats) + iLat

                      IsValidPoint = .true.

                      ! If point is ghostcell and not on the edge, it is not a good point
                      if (iLat < 1 .and. iBlockLat > 1) IsValidPoint = .false.
                      if (iLat > nLats .and. iBlockLat < nBlocksLat) IsValidPoint = .false.

                      if (IsValidPoint) then
                         do iAlt = -1, nAlts+2
                            GatGITMVars_g0(:,iiLon,iiLat,iAlt) = &
                                 GatGITMVars(:,iLon,iLat,iAlt,ip)
                         enddo
                      endif
                   enddo
                endif
             enddo

             !             do iLon = -1, nLons+2
             !                iiLon = (iBlockLon-1)*(nLons) + iLon
             !                do iLat = -1, nLats+2
             !                   iiLat = (iBlockLat-1)*(nLats) + iLat
             !                   do iAlt = -1, nAlts+2
             !                      iiAlt = iAlt
             !                      
             !                      GatGITMVars_g0(:,iiLon,iiLat,iiAlt) = &
             !                                         GatGITMVars(:,iLon,iLat,iAlt,ip)
             !                   enddo
             !                enddo
             !             enddo

          enddo

          GatAltitude_g0 = GatAltitude(-1:nAlts+2,1)/1000.0

          call MPI_send(GatLongitude_g0,&
               (nLonsTotal),&
               MPI_REAL,0,tag_SandR2,&
               intercomm1,iError)

          call MPI_send(GatLatitude_g0,&
               (nLatsTotal),&
               MPI_REAL,0,tag_SandR2,&
               intercomm1,iError)

          call MPI_send(GatAltitude_g0,&
               (nAltsTotal),&
               MPI_REAL,0,tag_SandR2,&
               intercomm1,iError)

          call MPI_send(GatGITMVars_g0,&
               nParams*(nLonsTotal)*(nLatsTotal)*(nAltsTotal),&
               MPI_REAL,0,tag_SandR2,&
               intercomm1,iError)

          deallocate (GatLongitude,GatLatitude,GatAltitude,GatGITMVars)

       endif

    else if (iCommSAMI0 /= MPI_COMM_NULL) then

       ! ----------------------------------------
       ! SAMI Receive Data
       ! ----------------------------------------

       call MPI_COMM_RANK(intercomm2,irank,iError)

       if (irank == 0) then

          ! -------------------------------
          ! Get all of the array sizes on PE 0
          ! -------------------------------

          ! First Get the sizes of the array:
          call MPI_RECV(nBlocks, 1, MPI_INTEGER, 0, tag_SandR1,&
               intercomm2, status, iError)

          call MPI_RECV(nprocs,  1, MPI_INTEGER, 0, tag_SandR1,&
               intercomm2,status,iError)

          call MPI_RECV(nBlocksLon, 1, MPI_INTEGER, 0, tag_SandR1,&
               intercomm2,status,iError)

          call MPI_RECV(nBlocksLat, 1, MPI_INTEGER, 0, tag_SandR1,&
               intercomm2,status,iError)

       endif

       ! if (irank==0) write(*,*) "> SAMI bcasting blocks"

       ! -------------------------------
       ! Pass Sizes to other PEs
       ! -------------------------------

       ! Use barriers to keep everyone in line:
       call MPI_BARRIER(iCommSAMI0,iError)
       call mpi_bcast(nBlocks,1,MPI_INTEGER,0,iCommSAMI0,iError)

       call MPI_BARRIER(iCommSAMI0,iError)
       call mpi_bcast(nprocs,1,MPI_INTEGER,0,iCommSAMI0,iError)

       call MPI_BARRIER(iCommSAMI0,iError)
       call mpi_bcast(nBlocksLon,1,MPI_INTEGER,0,iCommSAMI0,iError)

       call MPI_BARRIER(iCommSAMI0,iError)
       call mpi_bcast(nBlocksLat,1,MPI_INTEGER,0,iCommSAMI0,iError)

       call MPI_BARRIER(iCommSAMI0,iError)

       ! if (irank==0) write(*,*) "> SAMI bcasting blocks - done!"

       ! This is a bit troubling, since SAMI ~shouldn't~ know
       ! about nLons, nLats, nAlt
       nLonsTotal = nBlocksLon*nLons*nBlocks+4       
       nLatsTotal = nBlocksLat*nLats*nBlocks+4       
       nAltsTotal = nAlts+4       

       nPointsTotal = nLonsTotal * nLatsTotal * nAltsTotal

       if (IsFirstTime) then 
          allocate(GatLongitude_g(1:nLonsTotal))
          allocate(GatLatitude_g(1:nLatsTotal))
          allocate(GatAltitude_g(1:nAltsTotal))
          allocate(GatGITMVars_g(nParams,&
               1:nLonsTotal,&
               1:nLatsTotal,&
               1:nAltsTotal))
          IsFirstTime = .false.
       endif

       call MPI_COMM_RANK(intercomm2,irank,iError)

       ! -------------------------------
       ! Get all of the arrays on PE 0 from GITM
       ! -------------------------------

       if (irank == 0) then

          call MPI_RECV(GatLongitude_g, nLonsTotal, &
               MPI_REAL, 0, tag_SandR2, intercomm2, status, iError)

          call MPI_RECV(GatLatitude_g, nLatsTotal, &
               MPI_REAL, 0, tag_SandR2, intercomm2, status, iError)

          call MPI_RECV(GatAltitude_g, nAltsTotal, &
               MPI_REAL, 0, tag_SandR2, intercomm2, status, iError)

          call MPI_RECV(GatGITMVars_g, nParams * nPointsTotal,&
               MPI_REAL, 0, tag_SandR2, intercomm2, status, iError)

       endif

       ! -------------------------------
       ! Pass the data out to other SAMI PEs
       ! -------------------------------

       ! First spatial grid:
       call MPI_BARRIER(iCommSAMI0,iError)
       call mpi_bcast(GatLongitude_g, nLonsTotal, MPI_REAL, 0, iCommSAMI0, iError)

       call MPI_BARRIER(iCommSAMI0,iError)
       call mpi_bcast(GatLatitude_g, nLatsTotal, MPI_REAL, 0, iCommSAMI0, iError)

       call MPI_BARRIER(iCommSAMI0,iError)
       call mpi_bcast(GatAltitude_g, nAltsTotal, MPI_REAL, 0, iCommSAMI0, iError)

       call MPI_BARRIER(iCommSAMI0,iError)

       ! On Pleiades, passing a 4D array through mpi_bcast sometimes hangs.  So,
       ! loop through all of the parameters and send them one at a time:

       ! temp array:
       allocate(GatGITMVar(nLonsTotal,nLatsTotal,nAltsTotal))

       ! write(*,*) "SAMI Proc, before bcast : ", irank

       do i=1, nParams
          if (irank == 0) then
             ! write(*,*) "> var = ",i
             ! Move data into temp variable to bcast:
             GatGITMVar = GatGITMVars_g(i,:,:,:)
             ! write(*,*) "> done = ",i
          endif
          ! Pass temp variable to other PEs:
          call mpi_bcast(GatGITMVar,nPointsTotal,MPI_REAL,0,iCommSAMI0,iError)
          call MPI_BARRIER(iCommSAMI0,iError)
          if (irank > 0) then
             ! write(*,*) "> Moving bcast data: ",irank, i
             ! Move data from temp array to permenant array:
             GatGITMVars_g(i,:,:,:) = GatGITMVar
          endif
       enddo

       ! write(*,*) "SAMI Proc, after bcast : ", irank


       ! Deallocate temp array:
       deallocate(GatGITMVar)

       ! if (irank==0) write(*,*) "> SAMI bcasting - done!"

    endif

    if (iCommGITM /= MPI_COMM_NULL) then

       if (iproc == 0) &
            deallocate(GatLongitude_g0,GatLatitude_g0,GatAltitude_g0,GatGITMVars_g0)
    endif

    call MPI_BARRIER(iCommGlobal,iError)

    if (iCommGITM /= MPI_COMM_NULL) call MPI_COMM_FREE(intercomm1,iError)
    if (iCommSAMI0 /= MPI_COMM_NULL) call MPI_COMM_FREE(intercomm2,iError)

    call MPI_BARRIER(iCommGlobal,iError)

  endsubroutine ExchangeData

  subroutine overwrite_potential_sami(lons,lats,pot)

    use ModSamiInterp 
    use ModGITM, only: iproc,iCommGITM
    !use ModCoupSAMI3 
    use ModElectrodynamics

    implicit none

    real, dimension(-1:nLons+2,-1:nLats+2), intent(in) :: lons
    real, dimension(-1:nLons+2,-1:nLats+2), intent(in) :: lats
    real, dimension(-1:nLons+2,-1:nLats+2), intent(inout) :: pot

    save

    integer :: iLon, iLat
    real :: phi_xd00,phi_xd10,phi_yd0
    integer :: iError
    logical :: IsFirstTime = .true.
    real, allocatable :: phi_InterpFactors(:,:,:)
    real, allocatable :: phi_InterpIndex(:,:,:)
    real:: data_inp, data_inp2
    real :: phi_gitm_alt

    call report("potential_overwrite_sami",0)

    if (iproc == 0) &
         print*,'---- Gitm Start Overwrite Potential here!!!',iproc

    if (IsFirstTime) then

       allocate(phi_InterpIndex(5,-1:nLons+2,-1:nLats+2))
       allocate(phi_InterpFactors(3,-1:nLons+2,-1:nLats+2))

       do iLon = -1, nLons+2 
          do iLat = -1, nLats+2 

             phi_gitm_mlat = lats(iLon, iLat)

             if (phi_gitm_mlat < -90.0) phi_gitm_mlat = 180.0 + phi_gitm_mlat
             if (phi_gitm_mlat > 90.0)  phi_gitm_mlat = phi_gitm_mlat - 180.0

             phi_gitm_mlon = mod(lons(iLon, iLat)+360.0,360.0)

             !phi_gitm_mlon = MagLonMC(iLon,iLat)
             !phi_gitm_mlat = MagLatMC(iLon,iLat)
             phi_gitm_alt = 100.0

             !if (phi_gitm_mlon<0) phi_gitm_mlon = phi_gitm_mlon+360.0
             !print*,'===>> Magnetic Lon/Lat p',iproc,iLon,iLat,phi_gitm_mlon,phi_gitm_mlat,phi_gitm_alt

             call get_index_4points_2(phi_gitm_alt)

             phi_InterpIndex(1 ,iLon,iLat) = phi_i0
             phi_InterpIndex(2 ,iLon,iLat) = phi_j0
             phi_InterpIndex(3 ,iLon,iLat) = phi_i1
             phi_InterpIndex(4 ,iLon,iLat) = phi_j1
             phi_InterpIndex(5 ,iLon,iLat) = phi_flag

             !print*,phi_i0,phi_j0,phi_i1,phi_j1,phi_flag,phi_gitm_mlon,phi_gitm_mlat

             if (phi_flag == 1) then

                call get_factors_from_4points_sph(phi_gitm_mlon,phi_gitm_mlat,phi_xd00,phi_xd10,phi_yd0)

                if ((phi_xd00*phi_xd10*phi_yd0)<0) then
                   print*,'===>> Dist',phi_xd00,phi_xd10,phi_yd0
                   stop
                endif

                phi_InterpFactors(1 ,iLon,iLat) = phi_xd00
                phi_InterpFactors(2 ,iLon,iLat) = phi_xd10
                phi_InterpFactors(3 ,iLon,iLat) = phi_yd0

                !call  interp_with_4points(phi_InterpFactors(1,iLon,iLat),&
                !             phi_InterpFactors(2,iLon,iLat),phi_InterpFactors(3,iLon,iLat),&
                !             data_inp2)
             endif

          enddo
       enddo

       IsFirstTime = .false.

    endif

    if (iproc == 0) &
         print*,'-- GITM Couping Init Potential Done!!!',iproc

    do iLon = -1,nLons+2
       do iLat = -1,nLats+2

          phi_gitm_mlat = lats(iLon, iLat)
          if (phi_gitm_mlat < -90.0) phi_gitm_mlat = 180.0 + phi_gitm_mlat
          if (phi_gitm_mlat > 90.0)  phi_gitm_mlat = phi_gitm_mlat - 180.0

          phi_gitm_mlon = mod(lons(iLon, iLat)+360.0,360.0)

          ! phi_gitm_mlon = MagLonMC(iLon,iLat)
          ! phi_gitm_mlat = MagLatMC(iLon,iLat)

          phi_gitm_alt = 100.0

          !if (phi_gitm_mlon<0) phi_gitm_mlon = phi_gitm_mlon+360.0

          !print*,'===>> Magnetic Lon/Lat',iproc,iLon,iLat,phi_gitm_mlon,phi_gitm_mlat,phi_gitm_alt

          phi_flag = phi_InterpIndex(5 ,iLon,iLat) 

          if (phi_flag == 1) then

             call  interp_with_4points(phi_InterpFactors(1,iLon,iLat),&
                  phi_InterpFactors(2,iLon,iLat), &
                  phi_InterpFactors(3,iLon,iLat),&
                  data_inp)


             !print*,'---- DynamoPotentialMC,data_inp' , &
             !        pot(iLon,iLat),data_inp, &
             !        iLon,iLat,phi_gitm_mlon,phi_gitm_mlat

             pot(iLon,iLat) = data_inp

          endif
       enddo
    enddo

    if (iproc == 0) &
         print*,'GITM Overwrite_potential Done!',iproc
  end subroutine overwrite_potential_sami

end module ModCoupSAMI3
