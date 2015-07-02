!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

program PostProcess

  implicit none

  character (len=80) :: FileName, cLine
  character (len=6)  :: cBlock
  integer :: iStart, iError

  integer, parameter :: iInputUnitD_  = 13
  integer, parameter :: iInputUnitH_  = 14
  integer, parameter :: iOutputUnit_  = 15
  integer, parameter :: nMaxVars = 200
  integer :: Position = 0, PositionSave = 0

  logical :: IsThere, DoWrite, UseGhostCells

  integer :: iYear,iMonth,iDay,iHour,iMinute,iSecond,iMilli
  integer :: nBlocksLon, nBlocksLat, nBlocksAlt
  integer :: nVars, nLons, nLats, nAlts, nLonsTotal, nLatsTotal, nAltsTotal
  integer :: nGhostLats, nGhostLons, nGhostAlts
  integer :: iLon, iLat, iAlt, iBlockLon, iBlockLat, iBlockAlt
  integer :: iiLon, iiLat, iiAlt
  real    :: Version
  real    :: Data(nMaxVars)

  real, allocatable :: AllData(:,:,:,:)

  character (len=40) :: Variables(nMaxVars)

  integer :: i, iVar, iBlock

  logical :: IsEnd = .false., IsDone = .false., IsFirstTime = .true.

  nGhostLats = 2
  nGhostLons = 2
  nGhostAlts = 2

  write(*,*) "Enter file group to process (can include .header) : "
  read(5,*) FileName

  iStart = index(FileName,".header")-1
  if (iStart > 0) then 
  else
     iStart = index(FileName," ")-1
  endif

  inquire(file=FileName(1:iStart)//".header",EXIST=IsThere)
  if (.not.IsThere) then
     write(*,*) "Could not find header file : ",FileName(1:iStart)//".header"
     stop
  endif

  open(iInputUnitH_, file=FileName(1:iStart)//".header",status="old")

  do while (.not. IsDone)

     iError = 0

     nBlocksAlt = 1
     nBlocksLat = 1
     nBlocksLon = 1

     nVars = 0
     nAlts = 0

     UseGhostCells = .true.
     IsEnd = .false.

     do while (iError == 0)

        read(iInputUnitH_,'(a)',iostat=iError) cLine

        if (iError /= 0) IsDone = .true.

        if (index(cLine,"BLOCKS") > 0) then
           read(iInputUnitH_,*) nBlocksAlt
           read(iInputUnitH_,*) nBlocksLat
           read(iInputUnitH_,*) nBlocksLon
        endif

        if (index(cLine,"NUMERICAL") > 0) then
           read(iInputUnitH_,*) nVars
           read(iInputUnitH_,*) nAlts
           read(iInputUnitH_,*) nLats
           read(iInputUnitH_,*) nLons
        endif

        if (index(cLine,"TIME") > 0) then
           read(iInputUnitH_,*) iYear
           read(iInputUnitH_,*) iMonth
           read(iInputUnitH_,*) iDay
           read(iInputUnitH_,*) iHour
           read(iInputUnitH_,*) iMinute
           read(iInputUnitH_,*) iSecond
           read(iInputUnitH_,*) iMilli
        endif

        if (index(cLine,"VERSION") > 0) then
           read(iInputUnitH_,*) Version
        endif

        if (index(cLine,"VARIABLE") > 0) then
           do iVar = 1, nVars
              read(iInputUnitH_,"(i7,a40)") i, Variables(iVar)
           enddo
        endif

        if (index(cLine,"NGHOSTCELLS") > 0) then
           read(iInputUnitH_,*) nGhostAlts
           read(iInputUnitH_,*) nGhostLats
           read(iInputUnitH_,*) nGhostLons
        endif

        if (index(cLine,"NO GHOSTCELLS") > 0) then
           UseGhostCells = .false.
           nGhostLats = 0
           nGhostLons = 0
           nGhostAlts = 0
        endif

        if (index(cLine,"END") > 0) then
           iError = 1
           IsEnd = .true.
        endif

     enddo

     if (nAlts > 0) then 

        ! If it is 1D, then we know that we have to have 1,1,1,1.....
        if (nLons == 1 .and. nLats == 1) then
           nBlocksLon = 1
           nBlocksLat = 1
        endif

        write(*,*) "Inputs :"
        write(*,*) "  nBlocksLon, nBlocksLat, nBlocksAlt : ", &
             nBlocksLon, nBlocksLat, nBlocksAlt
        write(*,*) "  nLons,      nLats,      nAlts      : ", nLons, nLats, nAlts
        if (UseGhostCells) then
           nLonsTotal = nBlocksLon*(nLons-nGhostLons*2)+nGhostLons*2
           !     write(*,*) nBlocksLon, nLons, nGhostLons
           nLatsTotal = nBlocksLat*(nLats-nGhostLats*2)+nGhostLats*2
           !     write(*,*) nBlocksLat, nLats, nGhostLats, nLatsTotal
           nAltsTotal = nBlocksAlt*(nAlts-nGhostAlts*2)+nGhostAlts*2
        else
           nLonsTotal = nBlocksLon*nLons
           nLatsTotal = nBlocksLat*nLats
           nAltsTotal = nBlocksAlt*nAlts
        endif
        write(*,*) "  nLonsTotal, nLatsTotal, nAltsTotal : ", &
             nLonsTotal, nLatsTotal, nAltsTotal, " (Predicted Values)"

        allocate(AllData(nLonsTotal,nLatsTotal,nAltsTotal,nVars))

        if (.not. IsEnd .or. IsFirstTime) then
           write(*,*) "Not Appending...."
           open(iOutputUnit_,file=FileName(1:iStart)//".bin",&
                status="unknown",form="unformatted")
        else
           write(*,*) "Appending Satellite Files...."
           open(iOutputUnit_,file=FileName(1:iStart)//".bin",&
                status="old",form="unformatted",position="append")
        endif

        write(iOutputUnit_) Version
        write(iOutputUnit_) nLonsTotal, nLatsTotal, nAltsTotal
        write(iOutputUnit_) nVars
        do iVar=1,nVars
           write(iOutputUnit_) Variables(iVar)
        enddo

        write(iOutputUnit_) iYear, iMonth, iDay, iHour, iMinute, iSecond, iMilli

        nAltsTotal = 0
        nLatsTotal = 0
        nLonsTotal = 0

        iBlock = 1
        do iBlockAlt = 1, nBlocksAlt
           do iBlockLat = 1, nBlocksLat
              do iBlockLon = 1, nBlocksLon

                 write(cBlock,'(a2,i4.4)') ".b",iBlock

                 inquire(file=FileName(1:iStart)//cBlock,EXIST=IsThere)
                 if (.not.IsThere) then
                    write(*,*) "Must be a satellite file...."
                    cBlock = ".sat"
                 endif

                 !write(*,*) "Reading File : ",FileName(1:iStart)//cBlock

                 if (IsFirstTime .or. .not.IsEnd .or. iBlock > 1) then
                    !write(*,*) "Opening File!!!!"
                    open(iInputUnitD_, file=FileName(1:iStart)//cBlock, &
                         status="old", form="unformatted")
                    IsFirstTime = .false.
                 endif

                 do iAlt = 1, nAlts
                    do iLat = 1, nLats
                       do iLon = 1, nLons

                          read(iInputUnitD_) Data(1:nVars)

                          DoWrite = .true.

                          !! Only write ghostcells on the edges

                          if (UseGhostCells) then

                             if (iBlockLon /= 1 .and. iLon <= nGhostLons) &
                                  DoWrite = .false.
                             if (iBlockLon /= nBlocksLon .and. &
                                  iLon > nLons-nGhostLons) DoWrite = .false.

                             if (iBlockLat /= 1 .and. iLat <= nGhostLats) &
                                  DoWrite = .false.
                             if (iBlockLat /= nBlocksLat .and. &
                                  iLat > nLats-nGhostLats) DoWrite = .false.

                             if (iBlockAlt /= 1 .and. iAlt <= nGhostAlts) &
                                  DoWrite = .false.
                             if (iBlockAlt /= nBlocksAlt .and. &
                                  iAlt > nAlts-nGhostAlts) DoWrite = .false.

                          endif

                          if (DoWrite) then

                             if (UseGhostCells) then

                                iiLon = (iBlockLon-1)*(nLons-nGhostLons*2) + iLon
                                iiLat = (iBlockLat-1)*(nLats-nGhostLats*2) + iLat
                                iiAlt = (iBlockAlt-1)*(nAlts-nGhostAlts*2) + iAlt

                             else

                                iiLon = (iBlockLon-1)*(nLons) + iLon
                                iiLat = (iBlockLat-1)*(nLats) + iLat
                                iiAlt = (iBlockAlt-1)*(nAlts) + iAlt

                             endif

                             AllData(iiLon,iiLat,iiAlt,1:nVars) = Data(1:nVars)

                             if (iiLon > nLonsTotal) nLonsTotal = iiLon
                             if (iiLat > nLatsTotal) nLatsTotal = iiLat
                             if (iiAlt > nAltsTotal) nAltsTotal = iiAlt


                          endif

                       enddo
                    enddo
                 enddo

                 if (.not. IsEnd .or. nBlocksLat > 1) close(iInputUnitD_)

                 iBlock = iBlock + 1

              enddo
           enddo
        enddo

        do iVar = 1, nVars
           write(iOutputUnit_) AllData(:,:,:,iVar)
        enddo

        write(*,*) "Output :"
        write(*,*) "  nLons, nLats, nAlts : ", nLonsTotal, nLatsTotal, nAltsTotal
     
        close(iOutputUnit_)
        deallocate(AllData)

     endif

  enddo

  close(iInputUnitH_)
  close(iInputUnitD_)

end program PostProcess
