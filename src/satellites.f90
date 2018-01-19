!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!----------------------------------------------------------------------
! $Id: satellites.f90,v 1.17 2017/07/17 19:52:02 lregoli Exp $
!
! Author: Aaron Ridley, UMichigan
!
! Modified: 
!           AGB, Oct 2013 - Fixed RCMR flag useage and treatment of sat data
!           AGB, Sept 2013 - when you are not doing RCAC the data is not being
!                            read in, so you should not try to save it
!           Aaron, July 2013 - when you are not doing RCAC, then you shouldn't
!                              try to read data.  This was causing only every 
!                              other line to be read.
!           Alexey Morozov (alex), UM, May 2013 - fixed interpolation for the 
!                    case when satellites cross the prime meridian (lon=0)
!           Angeline Burrell (AGB), UMichigan, Feb 2013 - Modified variable
!                                                         names for consistency
!                                                         and added comments
!           Asad, UMichigan, Feb 2013 - Added data reading capability to allow
!                                       data assimilation
!           Aaron, May 2015 - fixed some error with cLine. Needed iCharLen_ to
!                             be used from ModInputs
!
! Comments: Routines to handle satellite data in GITM.
!
! Includes: read_satellites: reads space/time (and now data) from the satellite
!                            files specified in the UAM.in file.
!           move_satellites: Finds the current position and data value for
!                            each satellite at the current universal time
!----------------------------------------------------------------------

subroutine read_satellites(iError)

  use ModRCMR, only: RCMRFlag
  use ModSatellites
  use ModGITM, only: iUp_, iEast_, iNorth_
  use ModInputs, only: iDebugLevel, iCharLen_
  use ModConstants

  implicit none

  integer, intent(out)      :: iError
  integer                   :: iSat, iMax, i
  logical                   :: IsStartFound, IsFine
  character (len=iCharLen_) :: cLine
  real (kind = dblprec)     :: OldTime, NewTime
  integer                   :: itime(7)
  real                      :: pos(3), dat
  !---------------------------------------------------------------------
  if (RCMRFlag) then
     call init_mod_satellites_rcmr
  else
     call init_mod_satellites
  end if

  iError = 0
  nSatLines = 0
  nSatPos = 0
  iSatCurrentIndex = 0

  do iSat = 1, nSats

     if (iDebugLevel > 2) &
          write(*,*) "Reading Satellite File : ",cSatFileName(iSat), iSat

     open(unit=iSatUnit,file=cSatFileName(iSat),status="old",iostat = iError)

     if (iError /= 0) then
        write(*,*) "Error opening satellite file : ",cSatFileName(iSat)
        return
     endif

     IsStartFound = .false.

     do while (.not. IsStartFound)
        cLine = ""
        read(iSatUnit, *, iostat = iError) cLine
        if (iError /= 0) IsStartFound = .true.
        if (index(cline,"#START") > 0) IsStartFound = .true.
     enddo

     if (iError /= 0) then
        write(*,*) "Error finding #START in satellite file : ",&
             cSatFileName(iSat)
        close(iSatUnit)
        return
     endif

     OldTime = 0.0
     IsFine = .false.

     do while (iError == 0)

        ! Asad Note: Read in the time, 3D position, and a single data column
        if (RCMRFlag) then
           read(iSatUnit,*,iostat=iError) iTime, Pos, dat
        else
           read(iSatUnit,*,iostat=iError) iTime, Pos
        endif

        if (iError == 0) then

           IsFine = .true.

           ! Convert position data from degrees and km to radians and m
           if (Pos(iEast_) < 0.0) Pos(iEast_) = Pos(iEast_) + 360.0
           Pos(iEast_)  = Pos(iEast_)*pi/180.0
           Pos(iNorth_) = Pos(iNorth_)*pi/180.0
           Pos(iUp_)    = Pos(iUp_)*1000.0

           call time_int_to_real(iTime, NewTime)

           if (NewTime /= OldTime) then
              nSatLines(iSat) = nSatLines(iSat) + 1
              if (nSatLines(iSat) > nMaxSatInputLines) then
                 write(*,*) "Too Many Lines in satfile : ",cSatFileName(iSat)
                 iError = 1
              else
                 nSatPos(iSat, nSatLines(iSat)) = &
                      nSatPos(iSat, nSatLines(iSat)) + 1
                 SatPos(iSat, 1, nSatPos(iSat,nSatLines(iSat)), &
                      nSatLines(iSat)) = Pos(1)
                 SatPos(iSat, 2, nSatPos(iSat, nSatLines(iSat)), &
                      nSatLines(iSat)) = Pos(2)
                 SatPos(iSat, 3, nSatPos(iSat, nSatLines(iSat)), &
                      nSatLines(iSat)) = Pos(3)
                 SatTime(iSat, nSatLines(iSat)) = NewTime

                 if (RCMRFlag) SatDat(iSat, nSatLines(iSat)) = dat
              endif
           else
              nSatPos(iSat,nSatLines(iSat)) = nSatPos(iSat,nSatLines(iSat)) + 1
              if (nSatPos(iSat,nSatLines(iSat)) > nMaxSatPos-2) then
                 write(*,*) "Too Many Positions in satfile : ", &
                      cSatFileName(iSat)
                 write(*,*) "Line : ",nSatLines(iSat)
                 iError = 1
              else
                 SatPos(iSat,:,nSatPos(iSat,nSatLines(iSat)),nSatLines(iSat))=&
                      Pos
                 if (RCMRFlag) SatDat(iSat,nSatLines(iSat)) = dat
              endif
           endif

           OldTime = NewTime

        endif

     enddo

     close(iSatUnit)
     if (IsFine) iError = 0

  enddo

  if (iDebugLevel > 1) then

     write(*,*) "Number of Satellite Files Read : ",nSats
     do iSat = 1, nSats
        write(*,*) iSat,". Number of Different times : ",nSatLines(iSat)
        iMax = 1
        do i=1,nSatLines(iSat)
           if (nSatPos(iSat,i) > iMax) iMax = nSatPos(iSat,i)
        enddo
        write(*,*) iSat,". Max number of Positions : ",iMax
     enddo

  endif

end subroutine read_satellites

!----------------------------------------------------------------------
! Find locations for current time
!
! AGB/Asad 3/31/13: Adapted to use satellites in RCMR data assimilation
! Alexey 5/17/13: fixed the interpolation when satellite goes over 0lon
!----------------------------------------------------------------------

subroutine move_satellites

  use ModRCMR, only: RCMRFlag
  use ModSatellites
  use ModGITM, only: nBlocks, dt, iProc
  use ModTime, only: CurrentTime, tSimulation
  use ModConstants, only: pi !alex
  implicit none

  integer :: iSat, iPos, iLine = 1, i, iBlock, iOut
  real    :: r
  character (len=8)         :: cName
  character (len=3)         :: cPos
  character (len=iCharLen_) :: cTmp
  integer :: l1, l2

  iOut = -1
  if (RCMRFlag) iOut = -2

  do iSat = 1, nSats
if (SatDtPlot(iSat) < 0) then
!Then we only write the sat file if it is the time listed in the input file
iLine = iSatCurrentIndex(iSat)

if (iSatCurrentIndex(iSat) == 0) then

if ( CurrentTime >= SatTime(iSat,1) .and. &
CurrentTime <= SatTime(iSat,nSatLines(iSat))) then
iLine = 1

if (CurrentTime == SatTime(iSat,nSatLines(iSat))) then
iLine = nSatLines(iSat)
else
do while (SatTime(iSat,iLine) < CurrentTime)
iLine = iLine + 1
enddo
iLine = iLine - 1

endif
iSatCurrentIndex(iSat) = iLine
iLine = 0
endif
else
if (CurrentTime > SatTime(iSat,iSatCurrentIndex(iSat)+1) .and. &
(CurrentTime-dt) < SatTime(iSat,iSatCurrentIndex(iSat)+1)) then

iSatCurrentIndex(iSat) = iSatCurrentIndex(iSat) + 1

endif

endif

if (iSatCurrentIndex(iSat) /= iLine) then
!only plot if the index has changed!

iLine = iSatCurrentIndex(iSat)

do iPos = 1, nSatPos(iSat,iLine)
do i=1,3
if (i == 1 .and. &
SatPos(iSat, i, iPos, iLine) > 300 .and. &
SatPos(iSat, i, iPos, iLine+1) < 60) then
SatCurrentPos(iSat, i, iPos) = &
SatPos(iSat, i, iPos, iLine)
else
SatCurrentPos(iSat, i, iPos) = &
SatPos(iSat, i, iPos, iLine)

endif
enddo

do iBlock = 1, nBlocks

write(cPos,'(i3.3)') iPos
cTmp  = cSatFileName(iSat)
cName = cTmp(1:4)//"_"//cPos
CurrentSatellitePosition = SatCurrentPos(iSat,:,iPos)
CurrentSatelliteName     = cName
call output("UA/data/",iBlock, -1)
enddo
enddo
endif

else

     !Asad added this so that each satellite would be updated in output()
     CurrSat = iSat

     if (floor((tSimulation-dt)/SatDtPlot(iSat)) /= &
          floor((tsimulation)/SatDtPlot(iSat))) then

        if (iSatCurrentIndex(iSat) == 0) then
           if ( CurrentTime >= SatTime(iSat,1) .and. &
                CurrentTime <= SatTime(iSat,nSatLines(iSat))) then
              iLine = 1
              if (CurrentTime == SatTime(iSat,nSatLines(iSat))) then
                 iLine = nSatLines(iSat)
              else
                 do while (SatTime(iSat,iLine) < CurrentTime)
                    iLine = iLine + 1
                 enddo
                 iLine = iLine - 1
              endif
              iSatCurrentIndex(iSat) = iLine
           endif
        else
           iLine = iSatCurrentIndex(iSat)
           if (SatTime(iSat,nSatLines(iSat)) >= CurrentTime) then
              if (CurrentTime == SatTime(iSat,nSatLines(iSat))) then
                 iLine = nSatLines(iSat)
              else
                 do while (SatTime(iSat,iLine) < CurrentTime)
                    iLine = iLine + 1
                 enddo
                 iLine = iLine - 1
              endif
           else
              iLine = 0
           endif
           iSatCurrentIndex(iSat) = iLine
        endif

        if (iSatCurrentIndex(iSat) > 0) then

           if (iSatCurrentIndex(iSat) < nSatLines(iSat)) then
              iLine = iSatCurrentIndex(iSat)
              r =  1-(CurrentTime         - SatTime(iSat,iLine)) / &
                   (SatTime(iSat,iLine+1) - SatTime(iSat,iLine))
           else
              iLine = iSatCurrentIndex(iSat) - 1
              r = 1-(CurrentTime          - SatTime(iSat,iLine)) / &
                   (SatTime(iSat,iLine+1) - SatTime(iSat,iLine))
           endif

           if (RCMRFlag) then
              ! Asad: Estimate the current time's data value for each satellite
              SatCurrentDat(iSat) = r * SatDat(iSat, iLine) + &
                   (1 - r) * SatDat(iSat, iLine + 1)
           end if

           do iPos = 1, nSatPos(iSat,iLine)
              do i=1,3
                 if (i == 1 .and. &
                      SatPos(iSat, i, iPos, iLine) > 260*pi/180 .and. & !alex
                      SatPos(iSat, i, iPos, iLine+1) < 100*pi/180) then !alex
                    SatCurrentPos(iSat, i, iPos) = &
                        mod( (  r)*SatPos(iSat, i, iPos, iLine) + & !alex
                         (1-r)*(SatPos(iSat, i, iPos, iLine+1)+2*pi) , 2*pi) !alex
                 else
                    SatCurrentPos(iSat, i, iPos) = &
                         (  r)*SatPos(iSat, i, iPos, iLine) + &
                         (1-r)*SatPos(iSat, i, iPos, iLine+1)
                 endif

!                write(*,*) "Satellite : ",iSat, iSatCurrentIndex(iSat), &
!                     CurrentTime - SatTime(iSat,iSatCurrentIndex(iSat)), &
!                     iPos, i, SatCurrentPos(iSat, i, iPos)
              enddo
              write(cPos,'(i3.3)') iPos

              l1 = index(cSatFileName(iSat), '/', back=.true.) + 1
              l2 = index(cSatFileName(iSat), '.') - 1
              if (l1-1 <= 0) l1 = 1
              if (l2+1 <= 0) l2 = len_trim(cSatFileName(iSat))

              CurrentSatellitePosition = SatCurrentPos(iSat,:,iPos)
              CurrentSatelliteName     = cSatFileName(iSat)(l1:l2)//"_"//cPos
              do iBlock = 1, nBlocks
!                 call output_1d("UA/data/",cName,iBlock, &
!                      SatCurrentPos(iSat,:,iPos))
                 call output("UA/data/",iBlock, iOut)
              enddo
           enddo

        endif
    endif
     endif

  enddo

end subroutine move_satellites

