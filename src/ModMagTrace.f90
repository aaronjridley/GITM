!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModMagTrace

  integer, parameter :: MaxMMTPoints = 1500
  integer, allocatable :: MMTblk(:,:,:,:,:)

  logical, parameter :: MMTSaveInterp = .false., MMTDebug=.false.

  real, parameter :: MMTlen=1000.
  real, allocatable :: MMTalt(:,:,:,:,:), MMTlat(:,:,:,:,:), MMTlon(:,:,:,:,:)
  real, allocatable :: MMTaltLoc(:,:,:,:,:), MMTlatLoc(:,:,:,:,:), MMTlonLoc(:,:,:,:,:)

contains

  !-\
  !- This routine traces the magnetic field for every block starting from the lowest altitude.
  !- It stores the alt/lat/lon for every step along the trace in MMTlen increments.
  !-
  subroutine MMT_Init
    use ModGITM
    use ModMpi
    use ModElectrodynamics, ONLY: LengthFieldLine
    implicit none

    integer :: i,j,k,m,jLon,jLat, iBlock, iLon, iLat, iAlt, iLoop, iSize, iError, iCount=0

    logical :: IsFound, VerticalOnly=.false.

    real :: GeoLat, GeoLon, GeoAlt, xAlt, signz
    real :: GeoLatInitial, GeoLonInitial, GeoAltInitial
    real :: xmag, ymag, zmag, bmag
    real, dimension(MaxMMTPoints, -1:nLons+2, -1:nLats+2, nBlocksMax) :: LOCalt,LOClat,LOClon

    logical :: IsDone
    !--------------------

    ! Arrays to hold the alt/lat/lon of fieldline traces, synced to all processors
    if(.not.allocated(MMTalt)) &
         allocate(MMTalt(MaxMMTPoints, -1:nLons+2, -1:nLats+2, nBlocks, nProcs))
    if(.not.allocated(MMTlat)) &
         allocate(MMTlat(MaxMMTPoints, -1:nLons+2, -1:nLats+2, nBlocks, nProcs))
    if(.not.allocated(MMTlon)) &
         allocate(MMTlon(MaxMMTPoints, -1:nLons+2, -1:nLats+2, nBlocks, nProcs))
    MMTalt = 0.;
    MMTlat = 0.;
    MMTlon = 0.;

    if(MMTSaveInterp)then
       ! Arrays to hold the block and location values of point along fieldline for integrals
       if(.not.allocated(MMTblk)) &
            allocate(MMTblk(MaxMMTPoints, -1:nLons+2, -1:nLats+2, nBlocks, nProcs))
       if(.not.allocated(MMTaltLoc)) &
            allocate(MMTaltLoc(MaxMMTPoints, -1:nLons+2, -1:nLats+2, nBlocks, nProcs))
       if(.not.allocated(MMTlatLoc)) &
            allocate(MMTlatLoc(MaxMMTPoints, -1:nLons+2, -1:nLats+2, nBlocks, nProcs))
       if(.not.allocated(MMTlonLoc)) &
            allocate(MMTlonLoc(MaxMMTPoints, -1:nLons+2, -1:nLats+2, nBlocks, nProcs))
       MMTblk = 0;
       MMTaltLoc = -9.;
       MMTlatLoc = -9.;
       MMTlonLoc = -9.;
    end if

    if(MMTDebug) write(*,*)iProc,' MMT_Init:  ',nProcs,'  ', &
         nAlts,nLons,NLats,'  ',nBlocksMax,'  ',MaxMMTPoints

    do iBlock=1,nBlocks;  do iLat=-1,nLats+2;  do iLon=-1,nLons+2

       GeoLat = Latitude(iLat, iBlock)
       GeoLon = Longitude(iLon,iBlock)

       if (GeoLat > pi/2.) then
          GeoLat = pi - GeoLat
          GeoLon = mod(GeoLon + pi,twopi)
       endif
       if (GeoLat < -pi/2.) then
          GeoLat = -pi - GeoLat
          GeoLon = mod(GeoLon + pi,twopi)
       endif
       GeoLon = mod(GeoLon, twopi)
       if(GeoLon<0.) GeoLon=GeoLon+twopi

       GeoAlt = Altitude_GB(iLon,iLat,1,iBlock)
       IsDone = .false.
       xAlt = 1.0
       iAlt = 1

       GeoLatInitial = GeoLat
       GeoLonInitial = GeoLon
       GeoAltInitial = GeoAlt

       CALL get_magfield(GeoLat*180.0/pi,GeoLon*180.0/pi,GeoALT/1000.0,XMAG,YMAG,ZMAG)
       signz = sign(1.0,zmag)

       iLoop = 0
       do while (.not. IsDone)
          iLoop = iLoop +1
          if(iLoop > MaxMMTPoints)then
             write(*,*)'ERROR: increase size of MaxMMTPoints, ',MaxMMTPoints
             stop
          end if

          ! save values
          LOCalt(iLoop,iLon,iLat,iBlock) = GeoAlt
          if(VerticalOnly)then
             LOClat(iLoop,iLon,iLat,iBlock) = GeoLatInitial
             LOClon(iLoop,iLon,iLat,iBlock) = GeoLonInitial
          else
             LOClat(iLoop,iLon,iLat,iBlock) = GeoLat
             LOClon(iLoop,iLon,iLat,iBlock) = GeoLon
          end if

          ! this "integral" is computed now as the length is the integral of 1.0 along fieldline
          LengthFieldLine(iLon,iLat) = LengthFieldLine(iLon,iLat) + MMTlen

          CALL get_magfield(GeoLat*180.0/pi,GeoLon*180.0/pi,GeoALT/1000.0,XMAG,YMAG,ZMAG)

          if (sign(1.0,zmag)*signz < 0) then
             IsDone = .true.
          else
             bmag = sqrt(xmag*xmag + ymag*ymag + zmag*zmag)
             GeoAlt = GeoAlt + abs(zmag)/bmag * MMTlen
             if (GeoAlt > Altitude_GB(iLon,iLat,nAlts,iBlock)) then
                IsDone = .true.
             else
                if (GeoAlt > Altitude_GB(iLon,iLat,iAlt+1,iBlock)) iAlt = iAlt+1
                xAlt = (GeoAlt - Altitude_GB(iLon,iLat,iAlt,iBlock)) / &
                     ( Altitude_GB(iLon,iLat,iAlt+1,iBlock) &
                     - Altitude_GB(iLon,iLat,iAlt  ,iBlock))
                GeoLat = GeoLat + signz*xmag/bmag * MMTlen/(RBody + GeoAlt)*pi
                GeoLon = GeoLon + signz*ymag/bmag * MMTlen/(RBody + GeoAlt)*pi &
                     /cos(GeoLon)
                !    /(sign(1.,GeoLon)*max(0.01,abs(cos(GeoLon))))

                if (GeoLat > pi/2.) then
                   GeoLat = pi - GeoLat
                   GeoLon = mod(GeoLon + pi,twopi)
                endif
                if (GeoLat < -pi/2.) then
                   GeoLat = -pi - GeoLat
                   GeoLon = mod(GeoLon + pi,twopi)
                endif
                GeoLon = mod(GeoLon, twopi)
                if(GeoLon<0.) GeoLon=GeoLon+twopi

             endif
          endif

       enddo  !! while

    end do;  end do;  end do

    iSize = MaxMMTPoints * (nLons+4) * (nLats+4) * nBlocks
    call MPI_ALLGATHER(LOCalt, iSize, MPI_REAL, MMTalt, iSize, MPI_REAL, iCommGITM, iError)
    call MPI_ALLGATHER(LOClat, iSize, MPI_REAL, MMTlat, iSize, MPI_REAL, iCommGITM, iError)
    call MPI_ALLGATHER(LOClon, iSize, MPI_REAL, MMTlon, iSize, MPI_REAL, iCommGITM, iError)

    if(MMTSaveInterp)then
       do k=1,nProcs; do iBlock=1,nBlocks
          do iLat=-1,nLats+2;  do iLon=-1,nLons+2; do i=1,MaxMMTPoints
             GeoLat = MMTlat(i,iLon,iLat,iBlock,k)
             GeoLon = MMTlon(i,iLon,iLat,iBlock,k)
             GeoAlt = MMTalt(i,iLon,iLat,iBlock,k)

             IsFound = .false.
             if(GeoAlt > 1.)then
                do j=1,nBlocks
                   if(  GeoLat >=  ( Latitude(    0,j)+ Latitude(      1,j))/2. .and. &
                        GeoLat <   ( Latitude(nLats,j)+ Latitude(nLats+1,j))/2. .and. &
                        GeoLon >=  (Longitude(0    ,j)+Longitude(      1,j))/2. .and. &
                        GeoLon <   (Longitude(nLons,j)+Longitude(nLons+1,j))/2. )then
                      if(IsFound)then
                         write(*,*)'WARNING, BLOCK FOUND TWICE=] ',i,iLon,iLat,iBlock,k, &
                              '  ',GeoLat*180.0/pi,GeoLon*180.0/pi,GeoALT/1000.0
                      end if
                      IsFound = .true.
                      iCount = iCount +1

                      MMTblk(i,iLon,iLat,iBlock,k) = j
                      do m=0,nLats
                         if(  GeoLat>=Latitude(m  ,j) .and. &
                              GeoLat< Latitude(m+1,j))then
                            MMTlatLoc(i,iLon,iLat,iBlock,k) = m + &
                                 (GeoLat-Latitude(m,j))/&
                                 (Latitude(m+1,j)-Latitude(m,j))
                            jLat = m
                            exit
                         end if
                      end do
                      do m=0,nLons
                         if(  GeoLon>=Longitude(m  ,j) .and. &
                              GeoLon< Longitude(m+1,j))then
                            MMTlonLoc(i,iLon,iLat,iBlock,k) = m + &
                                 (GeoLon-Longitude(m,j))/&
                                 (Longitude(m+1,j)-Longitude(m,j))
                            jLon=m
                            exit
                         end if
                      end do
                      do m=1,nAlts
                         if(  GeoAlt>=Altitude_GB(jLon,jLat,m  ,j) .and. &
                              GeoAlt< Altitude_GB(jLon,jLat,m+1,j))then
                            MMTaltLoc(i,iLon,iLat,iBlock,k) = m + &
                                 (GeoAlt-Altitude_GB(jLon,jLat,m,j))/ &
                                 (Altitude_GB(jLon,jLat,m+1,j)-Altitude_GB(jLon,jLat,m,j))
                            exit
                         end if
                      end do
                      if(  MMTlatLoc(i,iLon,iLat,iBlock,k) == -9 .or. &
                           MMTlonLoc(i,iLon,iLat,iBlock,k) == -9 .or. &
                           MMTaltLoc(i,iLon,iLat,iBlock,k) == -9 )then
                         write(*,*)'WARNING, BLOCK VALUE NOT VALID=] ',i,iLon,iLat,iBlock,k, &
                              '  ',GeoLat*180.0/pi,GeoLon*180.0/pi,GeoALT/1000.0
                      end if

                   end if
                end do
             end if
          end do;  end do
       end do;  end do;  end do
       if(MMTDebug) write(*,*)iProc,' MMT_Init total processor fieldline length=',iCount*MMTlen
    end if

    if(MMTDebug) call MMT_Test
    if(MMTDebug) stop

  end subroutine MMT_Init
  !-/


  !-\
  !- This routine uses the previously stored trace positions to interpolate values for
  !- integrated quantities on that fieldline.
  !-
  subroutine MMT_Integrate(InValues,OutIntegral)
    use ModGITM
    use ModMpi
    implicit none

    integer :: i,j,k,m,jLat,jLon, jBlk, iBlock, iLon, iLat, iAlt, iSize, iError, iCount=0
    integer :: i1, i2, j1, j2, k1, k2

    logical :: IsFound

    real :: Dx1, Dx2, Dy1, Dy2, Dz1, Dz2
    real :: GeoLat, GeoLon, GeoAlt, InterpValue, lonLoc, latLoc, altLoc
    real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocksMax), intent(in) :: InValues
    real, dimension(-1:nLons+2, -1:nLats+2, nBlocksMax), intent(out) :: OutIntegral
    real, allocatable :: PartialIntegral(:,:,:,:), FullIntegral(:,:,:,:)

    !--------------------

    if(.not.allocated(PartialIntegral)) &
         allocate(PartialIntegral(-1:nLons+2, -1:nLats+2, nBlocks, nProcs))
    if(.not.allocated(FullIntegral)) &
         allocate(FullIntegral(-1:nLons+2, -1:nLats+2, nBlocks, nProcs))
    PartialIntegral = 0.
    FullIntegral = 0.
    do k=1,nProcs; do iBlock=1,nBlocks
       do iLat=-1,nLats+2;  do iLon=-1,nLons+2; do i=1,MaxMMTPoints
          GeoLat = MMTlat(i,iLon,iLat,iBlock,k)
          GeoLon = MMTlon(i,iLon,iLat,iBlock,k)
          GeoAlt = MMTalt(i,iLon,iLat,iBlock,k)
          if(MMTSaveInterp)then
             jBlk = MMTblk(i,iLon,iLat,iBlock,k)
             lonLoc = MMTlonLoc(i,iLon,iLat,iBlock,k)
             latLoc = MMTlatLoc(i,iLon,iLat,iBlock,k)  
             altLoc = MMTaltLoc(i,iLon,iLat,iBlock,k)
          else
             jBlk=-9;  lonLoc=-9.;  latLoc=-9.;  altLoc=-9.
             IsFound = .false.
             if(GeoAlt > 1.)then
                do j=1,nBlocks
                   if(  GeoLat >=  ( Latitude(    0,j)+ Latitude(      1,j))/2. .and. &
                        GeoLat <   ( Latitude(nLats,j)+ Latitude(nLats+1,j))/2. .and. &
                        GeoLon >=  (Longitude(0    ,j)+Longitude(      1,j))/2. .and. &
                        GeoLon <   (Longitude(nLons,j)+Longitude(nLons+1,j))/2. )then
                      if(IsFound)then
                         write(*,*)'WARNING, BLOCK FOUND TWICE=] ',i,iLon,iLat,iBlock,k, &
                              '  ',GeoLat*180.0/pi,GeoLon*180.0/pi,GeoALT/1000.0
                      end if
                      IsFound = .true.
                      iCount = iCount +1

                      jBlk = j
                      do m=0,nLats
                         if(  GeoLat>=Latitude(m  ,j) .and. &
                              GeoLat< Latitude(m+1,j))then
                            latLoc = m + &
                                 (GeoLat-Latitude(m,j))/&
                                 (Latitude(m+1,j)-Latitude(m,j))
                            jLat = m
                            exit
                         end if
                      end do
                      do m=0,nLons
                         if(  GeoLon>=Longitude(m  ,j) .and. &
                              GeoLon< Longitude(m+1,j))then
                            lonLoc = m + &
                                 (GeoLon-Longitude(m,j))/&
                                 (Longitude(m+1,j)-Longitude(m,j))
                            jLon=m
                            exit
                         end if
                      end do
                      do m=1,nAlts
                         if(  GeoAlt>=Altitude_GB(jLon,jLat,m  ,j) .and. &
                              GeoAlt< Altitude_GB(jLon,jLat,m+1,j))then
                            altLoc = m + &
                                 (GeoAlt-Altitude_GB(jLon,jLat,m,j))/ &
                                 (Altitude_GB(jLon,jLat,m+1,j)-Altitude_GB(jLon,jLat,m,j))
                            exit
                         end if
                      end do
                      if(  latLoc == -9 .or. &
                           lonLoc == -9 .or. &
                           altLoc == -9 )then
                         write(*,*)'WARNING, BLOCK VALUE NOT VALID=] ',i,iLon,iLat,iBlock,k, &
                              '  ',GeoLat*180.0/pi,GeoLon*180.0/pi,GeoALT/1000.0
                      end if

                   end if
                end do
             end if
          end if
          if(jBlk > 0)then
             !Set location
             i1 =   floor(lonLoc)
             i2 = ceiling(lonLoc)

             j1 =   floor(latLoc)  
             j2 = ceiling(latLoc)

             k1 =   floor(altLoc)
             k2 = ceiling(altLoc)

             !Set interpolation weights
             Dx1= lonLoc - i1;  Dx2 = 1.0 - Dx1
             Dy1= latLoc - j1;  Dy2 = 1.0 - Dy1
             Dz1= altLoc - k1;  Dz2 = 1.0 - Dz1

             !Perform interpolation
             InterpValue = Dz2*( Dy2*( Dx2*InValues(i1,j1,k1,jBlk)   &
                  +                    Dx1*InValues(i2,j1,k1,jBlk))  &
                  +              Dy1*( Dx2*InValues(i1,j2,k1,jBlk)   &
                  +                    Dx1*InValues(i2,j2,k1,jBlk))) &
                  +        Dz1*( Dy2*( Dx2*InValues(i1,j1,k2,jBlk)   &
                  +                    Dx1*InValues(i2,j1,k2,jBlk))  &
                  +              Dy1*( Dx2*InValues(i1,j2,k2,jBlk)   &
                  +                    Dx1*InValues(i2,j2,k2,jBlk)))
             InterpValue = InterpValue * MMTlen

             PartialIntegral(iLon,iLat,iBlock,k) = &
                  PartialIntegral(iLon,iLat,iBlock,k) + InterpValue
          end if
       end do;  end do
    end do;  end do;  end do

    if(.not.MMTSaveInterp)then
       if(MMTDebug) write(*,*)iProc,' MMT_Integrate total processor fieldline length=',iCount*MMTlen
    end if

    iSize = (nLons+4) * (nLats+4) * nBlocks * nProcs
    call MPI_AllREDUCE(PartialIntegral, FullIntegral, iSize, MPI_REAL, MPI_SUM, iCommGITM, iError)
    OutIntegral = FullIntegral(:,:,:,iProc+1)

  end subroutine MMT_Integrate
  !-/


  !-\
  !- This routine tests the MMT_Integrate routine
  !- Turn MMTDebug to T in MMT_Init and this will get called.
  !- The total fieldline points should equal the sum of the integral values
  !-
  subroutine MMT_Test
    use ModGITM
    implicit none

    integer :: i,j,k

    real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocksMax) :: Values
    real, dimension(-1:nLons+2, -1:nLats+2, nBlocksMax) :: Results

    !--------------------
    write(*,*)iProc,'Starting MMT_Test ...'

    Values=1.
    Results=0.
    call MMT_Integrate(Values,Results)

    write(*,*)iProc,'Integral min/max/sum values: ',minval(Results),maxval(Results),sum(Results)
    write(*,*)iProc,'                  ... finished MMT_Test'
  end subroutine MMT_Test
  !-/

end module ModMagTrace
