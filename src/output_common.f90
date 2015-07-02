!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!----------------------------------------------------------------------------
! $Id: output_common.f90,v 1.55 2014/10/29 22:23:28 xingm Exp $
!
! Author: Aaron Ridley, UMichigan
!
! Comments: Routines to output binary files
!
! AGB 3/31/13: Added 1D routine to output data at a specific altitude
! AJR 8/28/13: The code was outputting data on all processors for satellite
!              files.  I corrected this to make it so that if the linear
!              interpolation routine returns -1, the processor returns. 
! AGB 10/18/13: Added gravity, collision frequency, and pressure gradient
!               to 3DION output
! AGB 12/20/13: Removed gravity from 3DION output
!----------------------------------------------------------------------------

integer function bad_outputtype()

  use ModInputs, only : OutputType, nOutputTypes

  implicit none

  integer :: iOutputType
  logical :: IsFound

  do iOutputType = 1, nOutputTypes

     IsFound = .false.

     if (OutputType(iOutputType) == '3DALL')     IsFound = .true.
     if (OutputType(iOutputType) == '3DNEU')     IsFound = .true.
     if (OutputType(iOutputType) == '3DION')     IsFound = .true.
     if (OutputType(iOutputType) == '3DTHM')     IsFound = .true.
     if (OutputType(iOutputType) == '3DCHM')     IsFound = .true.
     if (OutputType(iOutputType) == '3DUSR')     IsFound = .true.
     if (OutputType(iOutputType) == '3DGLO')     IsFound = .true.
     if (OutputType(iOutputType) == '3DMAG')     IsFound = .true.

     if (OutputType(iOutputType) == '2DGEL')     IsFound = .true.
     if (OutputType(iOutputType) == '2DMEL')     IsFound = .true.
     if (OutputType(iOutputType) == '2DUSR')     IsFound = .true.
     if (OutputType(iOutputType) == '2DTEC')     IsFound = .true.

     if (OutputType(iOutputType) == '1DALL')     IsFound = .true.
     if (OutputType(iOutputType) == '0DALL')     IsFound = .true.
     if (OutputType(iOutputType) == '1DGLO')     IsFound = .true.
     if (OutputType(iOutputType) == '1DTHM')     IsFound = .true.
     if (OutputType(iOutputType) == '1DNEW')     IsFound = .true.
     if (OutputType(iOutputType) == '1DCHM')     IsFound = .true.
     if (OutputType(iOutputType) == '1DCMS')     IsFound = .true.
     if (OutputType(iOutputType) == '1DUSR')     IsFound = .true.

     if (.not. IsFound) then
        bad_outputtype = iOutputType
        return
     endif

  enddo

  bad_outputtype = 0
  return

end function bad_outputtype


!----------------------------------------------------------------
! Comments: Asad added data to allow output from RCAC
!----------------------------------------------------------------

subroutine output(dir, iBlock, iOutputType)

  use ModSphereInterface, only:iStartBlk
  use ModSatellites, only : CurrentSatellitePosition, CurrentSatelliteName, &
       CurrSat, SatAltDat
  use ModGITM
  use ModEUV
  use ModTime
  use ModInputs
  use ModSources
  use ModUserGITM, only: nVarsUser2d, nVarsUser3d, nVarsUser1d
  use ModRCMR, only: RCMRFlag

  implicit none

  character (len=*), intent(in) :: dir
  integer, intent(in) :: iBlock
  integer, intent(in) :: iOutputType

  character (len=5) :: proc_str,cBlock, cType
  character (len=24) :: cTime='', cTimeSave=''
  integer :: iiLat, iiLon, iiAlt, nGCs, cL=0
  integer :: iLon,iLat,iAlt, nVars_to_Write, nlines, iBLK,iSpecies
  logical :: done, IsFirstTime = .true., IsThere

  real :: LatFind, LonFind, AltFind
  real :: rLon, rLat, rAlt

  character (len=2) :: cYear, cMonth, cDay, cHour, cMinute, cSecond
  character (len=4) :: cYearL

  !! construct naming strings

  if (iOutputType == -1) then
     cType = "1DALL"
  else if(iOutputType == -2) then
     cType = "0DALL"
  else
     cType = OutputType(iOutputType)
     if (cType(1:2) == "3D" .and. Is1D) then 
        cType(1:2) = "1D"
     endif
  endif

  if (Is1D) then
     iiLat = 1
     iiLon = 1
     rLon = 1.0
     rLat = 1.0
  endif

  ! If there are satellites, initialize the current satellite so that
  ! the maximum from all processors will contain the real value.  This is
  ! done by setting the current value to something rediculously small for
  ! all currently known satellite input data types.
  
  if(CurrSat > 0 .and. RCMRFlag) then
     SatAltDat(CurrSat) = -1.0e32
  end if

  if (iOutputType <= -1) then
     LatFind = CurrentSatellitePosition(iNorth_)
     LonFind = CurrentSatellitePosition(iEast_)
     call BlockLocationIndex(LonFind,LatFind,iBlock,iiLon,iiLat,rLon,rLat)

     if(iOutputType == -2) then
        AltFind = CurrentSatellitePosition(iUp_)
        call BlockAltIndex(AltFind,iBlock,iiLon,iiLat,iiAlt,rAlt)

        if (iiAlt < 0) return
     end if
     
     if(iDebugLevel > 2)then
        write(*,*) 'For BlockLocationIndex:'
        write(*,*) 'LonFind, LatFind = ', LonFind, LatFind 
        write(*,*) 'Found iBlock, iiLon, iiLat, rLon, rLat =', &
             iBlock, iiLon, iiLat, rLon, rLat
     endif

     if (iiLon < 0 .or. iiLat < 0) return
  endif

  if((iProc == 0.and.iBlock == 1).and.(iOutputType /= -1)) &
       write(*,'(a,i7,i5,5i3)') &
       "Writing Output files ("//cType//") at iStep : ",&
       iStep, iTimeArray(1:6)

  if (iOutputType <= -1) &
       write(*,'(a,i7,i5,5i3)') &
       "Writing satellite file ("//trim(CurrentSatelliteName)//") at iStep : ",&
       iStep, iTimeArray(1:6)

  call calc_physics(iBlock)
  call calc_rates(iBlock)
  call calc_collisions(iBlock)
  call chapman_integrals(iBlock)
  call set_horizontal_bcs(iBlock)
  if (.not. Is1D) call calc_efield(iBlock)

  iBLK = iStartBLK + iBlock

  write(cBlock,'(a1,i4.4)') "b",iBLK

  call i2s(mod(iTimeArray(1),100), cYear, 2)
  call i2s(iTimeArray(1), cYearL, 4)
  call i2s(iTimeArray(2), cMonth, 2)
  call i2s(iTimeArray(3), cDay, 2)
  call i2s(iTimeArray(4), cHour, 2)
  call i2s(iTimeArray(5), cMinute, 2)
  call i2s(iTimeArray(6), cSecond, 2)

  if (.not. UseSecondsInFilename) cSecond='00'  !xianjing

  !-----------
  ! New feature - we want to be able to write to the same file over and
  ! over and over again, so we don't get 100,000,000 satellite files.
  ! So, here we are going to name the file the first time, then open
  ! that same file over and over again with an append.
  !-----------

  if ( IsFirstTime         .or. &
       .not. DoAppendFiles .or. &
       (iOutputType /= -1 .and. .not. Is1D .and. iOutputType /= -2)) then
     if (UseCCMCFileName) then
        cTime = "GITM_"//cYearL//"-"//cMonth//"-"//cDay//"T" &
             //cHour//"-"//cMinute//"-"//cSecond
        cL = 24
     else
        cTime = "t"//cYear//cMonth//cDay//"_"//cHour//cMinute//cSecond
        cL = 14
     endif
        
     if (IsFirstTime) cTimeSave = cTime
  else
     cTime = cTimeSave
  endif

  !! ---------------------------------------------
  !! Write the binary data files
  !! ---------------------------------------------

  if (iOutputType <= -1) then
     inquire(file=dir//"/"//trim(CurrentSatelliteName)//"_"//cTime(1:cL)//".sat", &
          EXIST=IsThere)
     if (.not. DoAppendFiles .or. tSimulation < 0.1 .or. .not. IsThere) then
        open(unit=iOutputUnit_, form="unformatted", &
             file=dir//"/"//trim(CurrentSatelliteName)//"_"//cTime(1:cL)//".sat",&
             status="unknown")
     else
        open(unit=iOutputUnit_, form="unformatted", &
             file=dir//"/"//trim(CurrentSatelliteName)//"_"//cTime(1:cL)//".sat",&
             status="unknown",position='append')
     endif
  else
     if (cType /= '2DMEL' .or. iBLK == 1) then
        inquire(file=dir//"/"//cType//"_"//cTime(1:cL)//"."//cBlock,&
             EXIST=IsThere)
        if (.not. DoAppendFiles .or. tSimulation < 0.1 .or. .not. IsThere) then
           open(unit=iOutputUnit_, form="unformatted", &
                file=dir//"/"//cType//"_"//cTime(1:cL)//"."//cBlock,&
                status="unknown")
        else
           open(unit=iOutputUnit_, form="unformatted", &
                file=dir//"/"//cType//"_"//cTime(1:cL)//"."//cBlock,&
                status="unknown",position='append')
        endif
     endif
  endif

  nGCs = 2

  select case (cType)

  case ('3DALL')

     nvars_to_write = 13+nSpeciesTotal+nSpecies+nIons
     call output_3dall(iBlock)

  case ('3DNEU')

     nvars_to_write = 8+nSpeciesTotal+nSpecies
     call output_3dneu(iBlock)

  case ('3DION')

     nvars_to_write = 8+nIons+6+4+4+1+4
     ! AGB: added nu_in (1) + pressure gradient (3)
     call output_3dion(iBlock)

  case ('3DTHM')

     nvars_to_write = 14
     call output_3dthm(iBlock)

  case ('1DCHM')

     nGCs = 0
     nvars_to_write = 30
     call output_1dchm(iBlock)

  case ('3DCHM')

     nvars_to_write = 30
     call output_3dchm(iBlock)

  case ('3DGLO')

     nvars_to_write = 3 + 3
     call output_3dglo(iBlock)

  case ('3DUSR')

     if (iBlock == 1) call set_nVarsUser3d
     nvars_to_write = nVarsUser3d
     call output_3duser(iBlock, iOutputUnit_)

  case ('3DMAG')

     nvars_to_write = 5+4
     call output_3dmag(iBlock)

  case ('2DGEL')

     nvars_to_write = 13
     call output_2dgel(iBlock)

  case ('2DMEL')

     nvars_to_write = 29
     if (iBLK == 1) call output_2dmel(iBlock)

  case ('2DUSR')

     if (iBlock == 1) call set_nVarsUser2d
     nvars_to_write = nVarsUser2d
     call output_2duser(iBlock, iOutputUnit_)

  case ('2DTEC')

     if (iBlock == 1) call set_nVarsUser2d
     nvars_to_write = 5
     call output_2dtec(iBlock)

  case('1DALL')

     nGCs = 0
     nvars_to_write = 13+nSpeciesTotal+nSpecies+nIons+nSpecies+5
     call output_1dall(iiLon, iiLat, iBlock, rLon, rLat, iOutputUnit_)

  case ('0DALL')
     ! AGB: added output type used by Asad to allow satellite output at the
     !      exact orbit location

     nGCs = 0
     nvars_to_write = 13+nSpeciesTotal+nSpecies+nIons+nSpecies+5
     call output_0dall(iiLon, iiLat, iiAlt, iBlock, rLon, rLat, rAlt, &
          iOutputUnit_)

  case ('1DGLO')

     nGCs = 0
     nvars_to_write = 6
     call output_1dglo

 case ('1DTHM')
     
     nGCs = 0
     nvars_to_write = 14 + (nspeciestotal*2)
     call output_1dthm

  case ('1DNEW')
     nGCs = 0
     nvars_to_write = 15 + nSpeciesTotal + nSpecies + nIons + nSpecies
     call output_1dnew(iiLon, iiLat, iBlock, rLon, rLat, iOutputUnit_)

  end select

  close(unit=iOutputUnit_)

  !! Now write the header file

  if ((iProc == 0 .and. iBlock == nBlocks) .or. iOutputType <= -1) then 

     if (iOutputType <= -1) then
        inquire(file=dir//"/"//trim(CurrentSatelliteName)//"_"//cTime(1:cL)//".header", EXIST=IsThere)
        if (.not. DoAppendFiles .or. tSimulation < 0.1 .or. .not. IsThere) then
           open(unit=iOutputUnit_, &
                file=dir//"/"//trim(CurrentSatelliteName)//"_"//cTime(1:cL)//".header",&
                status="unknown") 
        else
           open(unit=iOutputUnit_, &
                file=dir//"/"//trim(CurrentSatelliteName)//"_"//cTime(1:cL)//".header",&
                status="unknown",position='append')
        endif
     else
        inquire(file=dir//"/"//cType//"_"//cTime(1:cL)//".header",&
             EXIST=IsThere)
        if (.not. DoAppendFiles .or. tSimulation < 0.1 .or. .not. IsThere) then
           open(unit=iOutputUnit_, &
                file=dir//"/"//cType//"_"//cTime(1:cL)//".header",&
                status="unknown")
        else
           open(unit=iOutputUnit_, &
                file=dir//"/"//cType//"_"//cTime(1:cL)//".header",&
                status="unknown",position='append')
        endif
     endif

     call write_head_blocks
     call write_head_time
     call write_head_version

     if (cType(3:5) == 'USR') then
        call output_header_user(cType, iOutputUnit_)
     elseif (cType(3:5) == 'NEW') then
        call output_header_new
     else
        call output_header
     endif

     write(iOutputUnit_,*) ""
     write(iOutputUnit_,*) "END"
     write(iOutputUnit_,*) ""

     close(unit=iOutputUnit_)

  endif

  IsFirstTime = .false.

contains

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine output_header

    use ModElectrodynamics, only : nMagLats, nMagLons

    integer :: iOff, iSpecies, iIon

    write(iOutputUnit_,*) "NUMERICAL VALUES"

    write(iOutputUnit_,"(I7,6A)") nvars_to_write, " nvars"
    if (cType(1:2) /= "2D" .and. cType(1:2) /= '0D') then 
       write(iOutputUnit_,"(I7,7A)") nAlts+4, " nAltitudes"
    else
       write(iOutputUnit_,"(I7,7A)") 1, " nAltitudes"
    endif
    if (cType(1:2) == "1D" .or. cType(1:2) == '0D') then 
       write(iOutputUnit_,"(I7,7A)") 1, " nLatitudes"
       write(iOutputUnit_,"(I7,7A)") 1, " nLongitudes"
    else
       if (cType(3:5) =="MEL") then
          write(iOutputUnit_,"(I7,A)") nMagLats, " nLatitude"
          write(iOutputUnit_,"(I7,A)") nMagLons+1, " nLongitudes"
          write(iOutputUnit_,*) " "
          write(iOutputUnit_,*) "NO GHOSTCELLS"
       elseif (cType(3:5) =="GEL".or.cType(3:5)=="TEC") then
          write(iOutputUnit_,"(I7,A)") nLats, " nLatitude"
          write(iOutputUnit_,"(I7,A)") nLons, " nLongitudes"
          write(iOutputUnit_,*) " "
          write(iOutputUnit_,*) "NO GHOSTCELLS"
       else
          write(iOutputUnit_,"(I7,7A)") nLats+nGCs*2, " nLatitudes"
          write(iOutputUnit_,"(I7,7A)") nLons+nGCs*2, " nLongitudes"
       endif
    endif
    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "VARIABLE LIST"
    write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
    write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
    write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"

    if (cType(3:5) == "MAG") then
       write(iOutputUnit_,"(I7,A1,a)") iOff+5, " ", "Magnetic Latitude"
       write(iOutputUnit_,"(I7,A1,a)") iOff+6, " ", "Magnetic Longitude"
       write(iOutputUnit_,"(I7,A1,a)") iOff+8, " ", "B.F. East"
       write(iOutputUnit_,"(I7,A1,a)") iOff+9, " ", "B.F. North"
       write(iOutputUnit_,"(I7,A1,a)") iOff+10, " ", "B.F. Vertical"
       write(iOutputUnit_,"(I7,A1,a)") iOff+11, " ", "B.F. Magnitude"
    end if

    if (cType(3:5) == "GEL") then

       write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Potential"
       write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Pedersen Conductance"
       write(iOutputUnit_,"(I7,A1,a)")  6, " ", "Hall Conductance"
       write(iOutputUnit_,"(I7,A1,a)")  7, " ", "Electron_Average_Energy"
       write(iOutputUnit_,"(I7,A1,a)")  8, " ", "Electron_Energy_Flux"
       write(iOutputUnit_,"(I7,A1,a)")  9, " ", "DivJuAlt"
       write(iOutputUnit_,"(I7,A1,a)") 10, " ", "Pedersen FL Conductance"
       write(iOutputUnit_,"(I7,A1,a)") 11, " ", "Hall FL Conductance"
       write(iOutputUnit_,"(I7,A1,a)") 12, " ", "DivJu FL"
       write(iOutputUnit_,"(I7,A1,a)") 13, " ", "FL Length"

    endif

    if(cType(3:5) == "TEC") then

       write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Solar Zenith Angle"
       write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Vertical TEC"

    endif
    
    if (cType(3:5) == "THM") then

       write(iOutputUnit_,"(I7,A1,a)")  4, " ", "EUV Heating"
       write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Conduction"
       write(iOutputUnit_,"(I7,A1,a)")  6, " ", "Molecular Conduction"
       write(iOutputUnit_,"(I7,A1,a)")  7, " ", "Eddy Conduction"
       write(iOutputUnit_,"(I7,A1,a)")  8, " ", "Eddy Adiabatic Conduction"
       write(iOutputUnit_,"(I7,A1,a)")  9, " ", "Chemical Heating"
       write(iOutputUnit_,"(I7,A1,a)")  10, " ", "Auroral Heating"
       write(iOutputUnit_,"(I7,A1,a)")  11, " ", "Joule Heating"
       write(iOutputUnit_,"(I7,A1,a)")  12, " ", "NO Cooling"
       write(iOutputUnit_,"(I7,A1,a)")  13, " ", "O Cooling"
       write(iOutputUnit_,"(I7,A1,a)")  14, " ", "Total Abs EUV"
       if (cType(1:2) == "1D") then
          do iSpecies = 1, nSpeciesTotal
             write(iOutputUnit_,"(I7,A1,a,a)") 11 + iSpecies, " ", &
                  "Production Rate ",cSpecies(iSpecies)
          enddo
          do iSpecies = 1, nSpeciesTotal
             write(iOutputUnit_,"(I7,A1,a,a)") 11 + nSpeciesTotal + iSpecies, " ", &
                  "Loss Rate ",cSpecies(iSpecies)
             
          enddo
       endif
       
    endif

    if (cType(3:5) == "CHM") then

       write(iOutputUnit_,"(I7,A1,a)") 4, " ", "N!D2!U+!N + e"
       write(iOutputUnit_,"(I7,A1,a)") 5, " ", "O!D2!U+!N + e"
       write(iOutputUnit_,"(I7,A1,a)") 6, " ", "N!D2!U+!N + O"
       write(iOutputUnit_,"(I7,A1,a)") 7, " ", "NO!U+!N + e"
       write(iOutputUnit_,"(I7,A1,a)") 8, " ", "N!U+!N + O!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 9, " ", "NO + N"
       write(iOutputUnit_,"(I7,A1,a)") 10, " ","O!U+!N + O!D2!N" 
       write(iOutputUnit_,"(I7,A1,a)") 11, " ", "N + O!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 12, " ", "O!D2!U+!N + N"
       write(iOutputUnit_,"(I7,A1,a)") 13, " ", "O!D2!U+!N + NO"
       write(iOutputUnit_,"(I7,A1,a)") 14, " ", "O!D2!U+!N + N2"
       write(iOutputUnit_,"(I7,A1,a)") 15, " ", "N!D2!U+!N + O!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 16, " ", "N!U+!N + O"
       write(iOutputUnit_,"(I7,A1,a)") 17, " ", "O!+!N + N!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 18, " ", "O(1D) + N!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 19, " ", "O(1D) + O!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 20, " ", "O(1D) + O"
       write(iOutputUnit_,"(I7,A1,a)") 21, " ", "O(1D) + e"
       write(iOutputUnit_,"(I7,A1,a)") 22, " ", "N(2D) + O!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 23, " ", "O!U+!N(2D)+e"
       write(iOutputUnit_,"(I7,A1,a)") 24, " ", "N(2D) + O"
       write(iOutputUnit_,"(I7,A1,a)") 25, " ", "N(2D) + e"
       write(iOutputUnit_,"(I7,A1,a)") 26, " ", "O!U+!N(2D + N!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 27, " ", "O!U+!N(2P) + e"
       write(iOutputUnit_,"(I7,A1,a)") 28, " ", "O!U+!N(2P) + O"
       write(iOutputUnit_,"(I7,A1,a)") 29, " ", "O!U+!N(2P) + N!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 30, " ", "Chemical Heating Rate"
       
       
    endif

    if (cType(3:5) == "GLO") then

       write(iOutputUnit_,"(I7,A1,a)")  4, " ", "6300 A Emission"
       write(iOutputUnit_,"(I7,A1,a)")  5, " ", "PhotoElectronUp"
       write(iOutputUnit_,"(I7,A1,a)")  6, " ", "PhotoElectronDown"

    endif
       
    if (cType(3:5) == "MEL") then

       write(iOutputUnit_,"(I7,A1,a)")  4, " ", "MLT"
       write(iOutputUnit_,"(I7,A1,a)")  5, " ", "GeoLat"
       write(iOutputUnit_,"(I7,A1,a)")  6, " ", "GeoLon"
       write(iOutputUnit_,"(I7,A1,a)")  7, " ", "Pedersen Conductance"
       write(iOutputUnit_,"(I7,A1,a)")  8, " ", "Hall Conductance"
       write(iOutputUnit_,"(I7,A1,a)")  9, " ", "DivJuAlt"
       write(iOutputUnit_,"(I7,A1,a)") 10, " ", "Field Line Length"
       write(iOutputUnit_,"(I7,A1,a)") 11, " ", "Sigma PP"
       write(iOutputUnit_,"(I7,A1,a)") 12, " ", "Sigma LL"
       write(iOutputUnit_,"(I7,A1,a)") 13, " ", "Sigma H"
       write(iOutputUnit_,"(I7,A1,a)") 14, " ", "Sigma C"
       write(iOutputUnit_,"(I7,A1,a)") 15, " ", "Sigma PL"
       write(iOutputUnit_,"(I7,A1,a)") 16, " ", "Sigma LP"
       write(iOutputUnit_,"(I7,A1,a)") 17, " ", "K^D_{m\phi}"
       write(iOutputUnit_,"(I7,A1,a)") 18, " ", "K^D_{m\lamda}"
       write(iOutputUnit_,"(I7,A1,a)") 19, " ", "Solver A"
       write(iOutputUnit_,"(I7,A1,a)") 20, " ", "Solver B"
       write(iOutputUnit_,"(I7,A1,a)") 21, " ", "Solver C"
       write(iOutputUnit_,"(I7,A1,a)") 22, " ", "Solver D"
       write(iOutputUnit_,"(I7,A1,a)") 23, " ", "Solver E"
       write(iOutputUnit_,"(I7,A1,a)") 24, " ", "Solver S"
       write(iOutputUnit_,"(I7,A1,a)") 25, " ", "DynamoPotential"
       write(iOutputUnit_,"(I7,A1,a)") 26, " ", "Ed1new"
       write(iOutputUnit_,"(I7,A1,a)") 27, " ", "Ed2new"
       write(iOutputUnit_,"(I7,A1,a)") 28, " ", "Kphi"
       write(iOutputUnit_,"(I7,A1,a)") 29, " ", "Klamda"

    endif

    if (cType(3:5) == "ALL" .or. cType(3:5) == "NEU") then

       write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Rho"
   
       iOff = 4
       do iSpecies = 1, nSpeciesTotal
          write(iOutputUnit_,"(I7,A1,a)")  iOff+iSpecies, " ", &
               "["//cSpecies(iSpecies)//"]"
       enddo
    
       iOff = 4+nSpeciesTotal
       write(iOutputUnit_,"(I7,A1,a)")  iOff+1, " ", "Temperature"
       write(iOutputUnit_,"(I7,A1,a)")  iOff+2, " ", "V!Dn!N (east)"
       write(iOutputUnit_,"(I7,A1,a)")  iOff+3, " ", "V!Dn!N (north)"
       write(iOutputUnit_,"(I7,A1,a)")  iOff+4, " ", "V!Dn!N (up)"

       iOff = 8+nSpeciesTotal
       do iSpecies = 1, nSpecies
          write(iOutputUnit_,"(I7,A1,a)")  iOff+iSpecies, " ",&
               "V!Dn!N (up,"//cSpecies(iSpecies)//")"
       enddo

    endif

    if (cType(3:5) == "ALL" .or. cType(3:5) == "ION") then

       iOff = 3
       if (cType(3:5) == "ALL") iOff = 8+nSpeciesTotal+nSpecies
       do iIon = 1, nIons
          write(iOutputUnit_,"(I7,A1,a)") iOff+iIon, " ", "["//cIons(iIon)//"]"
       enddo

       iOff = iOff+nIons

       write(iOutputUnit_,"(I7,A1,a)") iOff+1, " ", "eTemperature"
       write(iOutputUnit_,"(I7,A1,a)") iOff+2, " ", "iTemperature"
       write(iOutputUnit_,"(I7,A1,a)") iOff+3, " ", "V!Di!N (east)"
       write(iOutputUnit_,"(I7,A1,a)") iOff+4, " ", "V!Di!N (north)"
       write(iOutputUnit_,"(I7,A1,a)") iOff+5, " ", "V!Di!N (up)"

       iOff = iOff + 5

       if (cType(3:5) == "ALL") then

          write(iOutputUnit_,"(I7,A1,a)") iOff+1, " ", "N2 Mixing Ratio"
          write(iOutputUnit_,"(I7,A1,a)") iOff+2, " ", "CH4 Mixing Ratio"
          write(iOutputUnit_,"(I7,A1,a)") iOff+3, " ", "Ar Mixing Ratio"
          write(iOutputUnit_,"(I7,A1,a)") iOff+4, " ", "HCN Mixing Ratio"
          write(iOutputUnit_,"(I7,A1,a)") iOff+5, " ", "H2 Mixing Ratio"

!       write(iOutputUnit_,"(I7,A1,a)") iOff+6, " ", "15N2 Mixing Ratio"
!       write(iOutputUnit_,"(I7,A1,a)") iOff+7, " ", "13CH4 Mixing Ratio"

          iOff = iOff + nSpecies
          write(iOutputUnit_,"(I7,A1,a)") iOff+1, " ", "RadCooling"
          write(iOutputUnit_,"(I7,A1,a)") iOff+2, " ", "EuvHeating"
          write(iOutputUnit_,"(I7,A1,a)") iOff+3, " ", "Conduction"
          write(iOutputUnit_,"(I7,A1,a)") iOff+4, " ", "Heat Balance Total"
          write(iOutputUnit_,"(I7,A1,a)") iOff+5, " ", "Heaing Efficiency"

       else

          write(iOutputUnit_,"(I7,A1,a)") iOff+1, " ", "Ed1"
          write(iOutputUnit_,"(I7,A1,a)") iOff+2, " ", "Ed2"
          write(iOutputUnit_,"(I7,A1,a)") iOff+3, " ", "Je1"
          write(iOutputUnit_,"(I7,A1,a)") iOff+4, " ", "Je2"
          write(iOutputUnit_,"(I7,A1,a)") iOff+5, " ", "Magnetic Latitude"
          write(iOutputUnit_,"(I7,A1,a)") iOff+6, " ", "Magnetic Longitude"
          write(iOutputUnit_,"(I7,A1,a)") iOff+7, " ", "B.F. East"
          write(iOutputUnit_,"(I7,A1,a)") iOff+8, " ", "B.F. North"
          write(iOutputUnit_,"(I7,A1,a)") iOff+9, " ", "B.F. Vertical"
          write(iOutputUnit_,"(I7,A1,a)") iOff+10, " ", "B.F. Magnitude"
          write(iOutputUnit_,"(I7,A1,a)") iOff+11, " ", "Potential"
          write(iOutputUnit_,"(I7,A1,a)") iOff+12, " ", "E.F. East"
          write(iOutputUnit_,"(I7,A1,a)") iOff+13, " ", "E.F. North"
          write(iOutputUnit_,"(I7,A1,a)") iOff+14, " ", "E.F. Vertical"
          write(iOutputUnit_,"(I7,A1,a)") iOff+15, " ", "E.F. Magnitude"

          ! AGB: 10/18/17: Add Collision Frequency and Pressure Gradient
          ! to output

          write(iOutputUnit_,"(I7,A1,a)") iOff+16, " ", "IN Collision Freq"
          write(iOutputUnit_,"(I7,A1,a)") iOff+17, " ", "PressGrad (east)"
          write(iOutputUnit_,"(I7,A1,a)") iOff+18, " ", "PressGrad (north)"
          write(iOutputUnit_,"(I7,A1,a)") iOff+19, " ", "PressGrad (up)"
       endif

    endif

    write(iOutputUnit_,*) ""

  end subroutine output_header


 subroutine output_header_new

    use ModElectrodynamics, only : nMagLats, nMagLons

    integer :: iOff, iSpecies, iIon

    write(iOutputUnit_,*) "NUMERICAL VALUES"

    write(iOutputUnit_,"(I7,6A)") nvars_to_write, " nvars"
    write(iOutputUnit_,"(I7,7A)") nAlts, " nAltitudes"
    write(iOutputUnit_,"(I7,7A)") 1, " nLatitudes"
    write(iOutputUnit_,"(I7,7A)") 1, " nLongitudes"
    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "VARIABLE LIST"
    write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
    write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Local Time"
    write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Latitude"
    write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Altitude"
    write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Solar Zenith Angle"
    write(iOutputUnit_,"(I7,A1,a)")  6, " ", "Rho"
    iOff = 6
    do iSpecies = 1, nSpeciesTotal
       write(iOutputUnit_,"(I7,A1,a)")  iOff+iSpecies, " ", &
            "["//cSpecies(iSpecies)//"]"
    enddo
    
    iOff = iOff + nSpeciesTotal
    write(iOutputUnit_,"(I7,A1,a)")  iOff+1, " ", "Temperature"
    write(iOutputUnit_,"(I7,A1,a)")  iOff+2, " ", "V!Dn!N (east)"
    write(iOutputUnit_,"(I7,A1,a)")  iOff+3, " ", "V!Dn!N (north)"
    write(iOutputUnit_,"(I7,A1,a)")  iOff+4, " ", "V!Dn!N (up)"

    iOff = iOff + 4
    do iSpecies = 1, nSpecies
       write(iOutputUnit_,"(I7,A1,a)")  iOff+iSpecies, " ",&
            "V!Dn!N (up,"//cSpecies(iSpecies)//")"
    enddo

    iOff = iOff + nSpecies
    do iIon = 1, nIons
       write(iOutputUnit_,"(I7,A1,a)") iOff+iIon, " ", "["//cIons(iIon)//"]"
    enddo

    iOff = iOff+nIons
    write(iOutputUnit_,"(I7,A1,a)") iOff+1, " ", "eTemperature"
    write(iOutputUnit_,"(I7,A1,a)") iOff+2, " ", "iTemperature"
    write(iOutputUnit_,"(I7,A1,a)") iOff+3, " ", "V!Di!N (east)"
    write(iOutputUnit_,"(I7,A1,a)") iOff+4, " ", "V!Di!N (north)"
    write(iOutputUnit_,"(I7,A1,a)") iOff+5, " ", "V!Di!N (up)"

    iOff = iOff + 5
    do iSpecies = 1, nSpecies
       write(iOutputUnit_,"(I7,A1,a)")  iOff+iSpecies, " ", &
            " "//cSpecies(iSpecies)//"Mixing Ratio"
    enddo


    write(iOutputUnit_,*) ""

  end subroutine output_header_new

  !----------------------------------------------------------------

  subroutine write_head_blocks

    if (cType(1:2) == "1D") return

    write(iOutputUnit_,*) "BLOCKS"
    write(iOutputUnit_,"(I7,A)") 1, " nBlocksAlt"
    if (cType /= "2DMEL") then
       write(iOutputUnit_,"(I7,A)") nBlocksLat, " nBlocksLat"
       write(iOutputUnit_,"(I7,A)") nBlocksLon, " nBlocksLon"
    else
       write(iOutputUnit_,"(I7,A)") 1, " nBlocksLat"
       write(iOutputUnit_,"(I7,A)") 1, " nBlocksLon"
    endif
    write(iOutputUnit_,*) ""

  end subroutine write_head_blocks


  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine write_head_time

    write(iOutputUnit_,*) "TIME"
    write(iOutputUnit_,"(I7,A)") iTimeArray(1), " Year"
    write(iOutputUnit_,"(I7,A)") iTimeArray(2), " Month"
    write(iOutputUnit_,"(I7,A)") iTimeArray(3), " Day"
    write(iOutputUnit_,"(I7,A)") iTimeArray(4), " Hour"
    write(iOutputUnit_,"(I7,A)") iTimeArray(5), " Minute"
    write(iOutputUnit_,"(I7,A)") iTimeArray(6), " Second"
    write(iOutputUnit_,"(I7,A)") iTimeArray(7), " Millisecond"
    write(iOutputUnit_,*) ""

  end subroutine write_head_time

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine write_head_version

    write(iOutputUnit_,*) "VERSION"
    write(iOutputUnit_,*) 3.1+PlanetNum
    write(iOutputUnit_,*) ""

  end subroutine write_head_version

end subroutine output

!!  !----------------------------------------------------------------
!!  !
!!  !----------------------------------------------------------------
!!  
!!  subroutine output_1d(dir, cName, iBlock, Position)
!!  
!!    use ModGITM
!!    use ModEUV
!!    use ModTime
!!    use ModInputs
!!    use ModSources
!!    use ModConstants
!!  
!!    implicit none
!!  
!!    character (len=*), intent(in) :: dir
!!    character (len=*), intent(in) :: cName
!!    integer, intent(in)           :: iBlock
!!    real, intent(in)              :: Position(3)
!!  
!!    character (len=14) :: cTime
!!    integer :: iLon,iLat,iAlt, nvars_to_write, nlines, iBLK, i
!!    integer :: iiLon, iiLat, iiAlt
!!    logical :: done
!!  
!!    real :: LatFind, LonFind
!!    real :: rLon, rLat
!!  
!!    character (len=2) :: cYear, cMonth, cDay, cHour, cMinute, cSecond
!!  
!!    LatFind = Position(iNorth_)
!!    LonFind = Position(iEast_)
!!  
!!    if ((Longitude(0,iBlock)+Longitude(1,iBlock))/2 <=LonFind .and. &
!!         (Longitude(nLons,iBlock)+Longitude(nLons+1,iBlock))/2 >LonFind) then
!!       if ((Latitude(0,iBlock)+Latitude(1,iBlock))/2 <=LatFind .and. &
!!            (Latitude(nLats,iBlock)+Latitude(nLats+1,iBlock))/2 >LatFind) then
!!  
!!          iiLat = -1
!!          iiLon = -1
!!          do iLon = 0,nLons
!!             if (Longitude(iLon,iBlock) <= LonFind .and. &
!!                  Longitude(iLon+1,iBlock) > LonFind) then
!!                iiLon = iLon
!!                rLon = 1.0 - (LonFind - Longitude(iLon,iBlock)) / &
!!                     (Longitude(iLon+1,iBlock) - Longitude(iLon,iBlock))
!!             endif
!!          enddo
!!  
!!          do iLat = 0,nLats
!!             if (Latitude(iLat,iBlock) <= LatFind .and. &
!!                  Latitude(iLat+1,iBlock) > LatFind) then
!!                iiLat = iLat
!!                rLat = 1.0 - (LatFind - Latitude(iLat,iBlock)) / &
!!                     (Latitude(iLat+1,iBlock) - Latitude(iLat,iBlock))
!!             endif
!!          enddo
!!  
!!       else
!!          return
!!       endif
!!    else 
!!       return
!!    endif
!!  
!!    if (iProc == 0 .and. iBlock == 1) &
!!         write(*,'(a,i7,i5,5i3)') &
!!         "Writing Output files at iStep : ",iStep, iTimeArray(1:6)
!!  
!!    call calc_physics(iBlock)
!!    call chapman_integrals(iBlock)
!!    call calc_rates(iBlock)
!!    if (.not. Is1D) call calc_efield(iBlock)
!!  
!!    !! construct naming strings
!!  
!!    call i2s(mod(iTimeArray(1),100), cYear, 2)
!!    call i2s(iTimeArray(2), cMonth, 2)
!!    call i2s(iTimeArray(3), cDay, 2)
!!    call i2s(iTimeArray(4), cHour, 2)
!!    call i2s(iTimeArray(5), cMinute, 2)
!!    call i2s(iTimeArray(6), cSecond, 2)
!!  
!!    cTime = "t"//cYear//cMonth//cDay//"_"//cHour//cMinute//cSecond
!!  
!!    !! open file
!!    open(unit=iOutputUnit_, &
!!         file=dir//"/"//cName//"_"//cTime(1:cL)//"."//"dat",&
!!         status="unknown")
!!  
!!    write(iOutputUnit_,*) ""
!!    write(iOutputUnit_,*) "TIME"
!!    do i=1,7
!!       write(iOutputUnit_,*) iTimeArray(i)
!!    enddo
!!    write(iOutputUnit_,*) ""
!!  
!!    call output_header(.true., nVars_to_Write)
!!  
!!    call output_1dall(iiLon, iiLat, iBlock, rLon, rLat, iOutputUnit_)
!!  
!!    close(unit=iOutputUnit_)
!!  
!!  end subroutine output_1d

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_3dall(iBlock)

  use ModGITM
  use ModInputs

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon, iiAlt, iiLat, iiLon, i

  do iAlt=-1,nAlts+2
     !!! Why ???
     iiAlt = max(min(iAlt,nAlts),1)
     do iLat=-1,nLats+2
        !!! Why ???
        iiLat = min(max(iLat,1),nLats)
        do iLon=-1,nLons+2
           !!! Why ???
           iiLon = min(max(iLon,1),nLons)
           write(iOutputUnit_)       &
                Longitude(iLon,iBlock), &
                Latitude(iLat,iBlock), &
                Altitude_GB(iLon,iLat,iAlt,iBlock),&
                Rho(iLon,iLat,iAlt,iBlock),&
                (NDensityS(iLon,iLat,iAlt,i,iBlock),i=1,nSpeciesTotal), &
                Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt),&
                (Velocity(iLon,iLat,iAlt,i,iBlock),i=1,3), &
                (VerticalVelocity(iLon,iLat,iAlt,i,iBlock),i=1,nSpecies), &
                (IDensityS(iLon,iLat,iAlt,i,iBlock),i=1,nIons), &
                eTemperature(iLon,iLat,iAlt,iBlock)  ,&
                ITemperature(iLon,iLat,iAlt,iBlock)  ,&
                (Ivelocity(iLon,iLat,iAlt,i,iBlock),i=1,3)
        enddo
     enddo
  enddo

end subroutine output_3dall

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_3dneu(iBlock)

  use ModGITM
  use ModInputs

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon, iiAlt, iiLat, iiLon

  do iAlt=-1,nAlts+2
     iiAlt = max(min(iAlt,nAlts),1)
     do iLat=-1,nLats+2
        iiLat = min(max(iLat,1),nLats)
        do iLon=-1,nLons+2
           iiLon = min(max(iLon,1),nLons)
           write(iOutputUnit_)       &
                Longitude(iLon,iBlock), &
                Latitude(iLat,iBlock), &
                Altitude_GB(iLon,iLat,iAlt,iBlock),&
                Rho(iLon,iLat,iAlt,iBlock),&
                NDensityS(iLon,iLat,iAlt,:,iBlock), &
                Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt),&
                velocity(iLon,iLat,iAlt,:,iBlock) , &
                VerticalVelocity(iLon,iLat,iAlt,:,iBlock)
        enddo
     enddo
  enddo

end subroutine output_3dneu

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_3dion(iBlock)

  use ModGITM
  use ModInputs
  use ModElectrodynamics

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon, iiAlt, iiLat, iiLon

  do iAlt=-1,nAlts+2
     do iLat=-1,nLats+2
        do iLon=-1,nLons+2
           write(iOutputUnit_)          &
                Longitude(iLon,iBlock),               &
                Latitude(iLat,iBlock),                &
                Altitude_GB(iLon,iLat,iAlt,iBlock),   &
                IDensityS(iLon,iLat,iAlt,:,iBlock),   &
                eTemperature(iLon,iLat,iAlt,iBlock),  &
                ITemperature(iLon,iLat,iAlt,iBlock),  &
                Ivelocity(iLon,iLat,iAlt,:,iBlock),   &
                ed1(iLon,iLat,iAlt), &
 		ed2(iLon,iLat,iAlt), &
		je1(iLon,iLat,iAlt), &
		je2(iLon,iLat,iAlt), &
                mLatitude(iLon,iLat,iAlt,iBlock), &
                mLongitude(iLon,iLat,iAlt,iBlock), &
                B0(iLon,iLat,iAlt,:,iBlock), &  !Geomagnetic B0(nLons,nLats,nAlts,4[iEast_,iNorth_,iUp_,iMag_],nBlocks)
                potential(iLon,iLat,iAlt,iBlock), &
                EField(iLon,iLat,iAlt,:), &  ! EField(Lon,lat,alt,3)
                sqrt(sum(EField(iLon,iLat,iAlt,:)**2)), & ! magnitude of E.F.
                Collisions(iLon,iLat,iAlt,iVIN_), & ! AGB: nu_in
                PressureGradient(iLon,iLat,iAlt,:,iBlock) ! AGB: 3D Grad P
        enddo
     enddo
  enddo

end subroutine output_3dion

!----------------------------------------------------------------
!
!----------------------------------------------------------------
subroutine output_3dthm(iBlock)

  use ModGITM
  use ModInputs
  use ModSources
  use ModEuv, only : EuvTotal
  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon, iiAlt, iiLat, iiLon


  do iAlt=-1,nAlts+2
     iiAlt = max(min(iAlt,nAlts),1)
     do iLat=-1,nLats+2
        iiLat = min(max(iLat,1),nLats)
        do iLon=-1,nLons+2
           iiLon = min(max(iLon,1),nLons)

           write(iOutputUnit_)          &
                Longitude(iLon,iBlock),               &
                Latitude(iLat,iBlock),                &
                Altitude_GB(iLon,iLat,iAlt,iBlock),   &
                EuvHeating(iiLon,iiLat,iiAlt,iBlock)*dt*TempUnit(iiLon,iiLat,iiAlt),    &
                Conduction(iiLon,iiLat,iiAlt)*TempUnit(iiLon,iiLat,iiAlt),              &
                MoleConduction(iiLon,iiLat,iiAlt),                             &
                EddyCond(iiLon,iiLat,iiAlt),                                   &
                EddyCondAdia(iiLon,iiLat,iiAlt),                               &
                ChemicalHeatingRate(iiLon,iiLat,iiAlt)*TempUnit(iiLon,iiLat,iiAlt),     &
                AuroralHeating(iiLon,iiLat,iiAlt)*dt*TempUnit(iiLon,iiLat,iiAlt),       &
                JouleHeating(iiLon,iiLat,iiAlt)*dt*TempUnit(iiLon,iiLat,iiAlt),         &
                -NOCooling(iiLon,iiLat,iiAlt)*dt*TempUnit(iiLon,iiLat,iiAlt),           &
                -OCooling(iiLon,iiLat,iiAlt)*dt*TempUnit(iiLon,iiLat,iiAlt),            &
                EuvTotal(iiLon,iiLat,iiAlt,iBlock) * dt
           
        enddo
     enddo
  enddo
     
end subroutine output_3dthm

!----------------------------------------------------------------
!
!----------------------------------------------------------------


subroutine output_1dthm

  use ModGITM
  use ModInputs
  use ModSources
  use ModEuv, only : EuvTotal
  implicit none

  integer :: iAlt, iLat, iLon, iiAlt,iSpecies
  real    :: varsS(nSpeciesTotal),varsL(nSpeciesTotal)

  
  do iAlt=-1,nAlts+2
     iiAlt = max(min(iAlt,nAlts),1)
 
     do iSpecies = 1, nSpeciesTotal 
        varsS(iSpecies) = NeutralSourcesTotal(iialt,iSpecies)
        varsL(iSpecies) = NeutralLossesTotal(iialt,iSpecies)
     enddo

     write(iOutputUnit_) &
          Longitude(1,1),               &
          Latitude(1,1),                &
          Altitude_GB(1,1,iAlt,1),   &
          EuvHeating(1,1,iiAlt,1)*dt*TempUnit(1,1,iiAlt),    &
          Conduction(1,1,iiAlt)*TempUnit(1,1,iiAlt),              &
          MoleConduction(1,1,iiAlt),                             &
          EddyCond(1,1,iiAlt),                                   &
          EddyCondAdia(1,1,iiAlt),                               &
          ChemicalHeatingRate(1,1,iiAlt)*TempUnit(1,1,iiAlt),     &
          AuroralHeating(1,1,iiAlt)*dt*TempUnit(1,1,iiAlt),       &
          JouleHeating(1,1,iiAlt)*dt*TempUnit(1,1,iiAlt),         &
          -RadCooling(1,1,iiAlt,1)*dt*TempUnit(1,1,iiAlt),           &
          -OCooling(1,1,iiAlt)*dt*TempUnit(1,1,iiAlt),            &
          EuvTotal(1,1,iiAlt,1) * dt,                             &
          varsS, varsL
                   
  enddo

end subroutine output_1dthm

!----------------------------------------------------------------
!
!----------------------------------------------------------------
subroutine output_1dchm(iBlock)

  use ModGITM
  use ModInputs
  use ModSources
  use ModConstants
  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iiAlt, iReact
  real :: vars(nReactions)

  do iAlt=-1,nAlts+2
     iiAlt = max(min(iAlt,nAlts),1)
     do iReact = 1, nReactions
        vars(iReact) = ChemicalHeatingSpecies(1,1,iiAlt,iReact) / &
             Element_Charge
     enddo
              
     write(iOutputUnit_) &
          Longitude(1,iBlock),               &
          Latitude(1,iBlock),                &
          Altitude_GB(1,1,iAlt,iBlock),   &
          Vars, &
          ChemicalHeatingRate(1,1,iiAlt) * &
          cp(1,1,iiAlt,iBlock) *   &
          Rho(1,1,iiAlt,iBlock)*TempUnit(1,1,iiAlt) / &
          Element_Charge
  enddo

end subroutine output_1dchm


!----------------------------------------------------------------
!
!----------------------------------------------------------------
subroutine output_3dchm(iBlock)

  use ModGITM
  use ModInputs
  use ModSources
  use ModConstants
  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon, iiAlt, iiLat, iiLon, iReact
  real :: vars(nReactions)

  do iAlt=-1,nAlts+2
     iiAlt = max(min(iAlt,nAlts),1)
     do iLat=-1,nLats+2
        iiLat = min(max(iLat,1),nLats)
        do iLon=-1,nLons+2
           iiLon = min(max(iLon,1),nLons)
           do iReact = 1, nReactions
              
              vars(iReact) = ChemicalHeatingSpecies(iiLon,iiLat,iiAlt,iReact) / &
                   Element_Charge
              
              enddo
              
              write(iOutputUnit_) &
                   Longitude(iLon,iBlock),               &
                   Latitude(iLat,iBlock),                &
                   Altitude_GB(iLon,iLat,iAlt,iBlock),   &
                   Vars, &
                   ChemicalHeatingRate(iiLon,iiLat,iiAlt) * &
                   cp(iilon,iiLat,iiAlt,iBlock) *   &
                   Rho(iilon,iiLat,iiAlt,iBlock)*TempUnit(iilon,iiLat,iiAlt) / &
                   Element_Charge
        enddo
     enddo
  enddo

end subroutine output_3dchm

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_3dglo(iBlock)

  use ModGITM
  use ModInputs

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon, iiAlt, iiLat, iiLon

  return

!  do iAlt=-1,nAlts+2
!     iiAlt = max(min(iAlt,nAlts),1)
!     do iLat=-1,nLats+2
!        iiLat = min(max(iLat,1),nLats)
!        do iLon=-1,nLons+2
!           iiLon = min(max(iLon,1),nLons)
!              
!              write(iOutputUnit_) &
!                   Longitude(iLon,iBlock),               &
!                   Latitude(iLat,iBlock),                &
!                   Altitude_GB(iLon,iLat,iAlt,iBlock),   &
!                   vEmissionRate(iiLon,iiLat,iiAlt,i6300_,iBlock), &
!                   PhotoEFluxTotal(iiLon,iiLat,iiAlt,iBlock,1),    &
!                   PhotoEFluxTotal(iiLon,iiLat,iiAlt,iBlock,2)
!                   
!        enddo
!     enddo
!  enddo
!
end subroutine output_3dglo

!----------------------------------------------------------------
!
!----------------------------------------------------------------
subroutine output_1dglo

  use ModGITM
  use ModInputs

  implicit none

  integer :: iAlt, iLat, iLon, iiAlt

!  do iAlt=-1,nAlts+2
!     iiAlt = max(min(iAlt,nAlts),1)
! 
!     write(iOutputUnit_) &
!          Longitude(1,1),               &
!          Latitude(1,1),                &
!          Altitude_GB(1,1,iAlt,1),   &
!          vEmissionRate(1,1,iiAlt,i6300_,1), &
!          PhotoEFluxTotal(1,1,iiAlt,1,1),    &
!          PhotoEFluxTotal(1,1,iiAlt,1,2)
!                   
!  enddo

  return

end subroutine output_1dglo

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_2dgel(iBlock)

  use ModElectrodynamics
  use ModGITM
  use ModInputs

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon, iiAlt, iiLat, iiLon

  iAlt = 1
  do iLat=1,nLats
     do iLon=1,nLons
        write(iOutputUnit_)       &
             Longitude(iLon,iBlock), &
             Latitude(iLat,iBlock),&
             Altitude_GB(iLon,iLat,iAlt,iBlock), &
             Potential(iLon,iLat,iAlt,iBlock), &
             PedersenConductance(iLon,iLat,iBlock), &
             HallConductance(iLon,iLat,iBlock), &
             ElectronAverageEnergy(iLon,iLat), &
             ElectronEnergyFlux(iLon,iLat), &
             DivJuAlt(iLon,iLat), &
             PedersenFieldLine(iLon, iLat), &
             HallFieldLine(iLon, iLat), &
             DivJuFieldLine(iLon, iLat), &
             LengthFieldLine(iLon, iLat)
     enddo
  enddo

end subroutine output_2dgel

!-------------------------------------------------------------------------------
! AGB: Routine to output a 2D TEC file that includes the lat, lon, SZA, and VTEC
!-------------------------------------------------------------------------------

subroutine output_2dtec(iBlock)

  use ModGITM
  use ModInputs
  use ModEUV, only : Sza

  implicit none

  integer, intent(in) :: iBlock
  integer :: iLat, iLon, iAlt, iiLat, iiLon

  call calc_vtec(iBlock)

  iAlt = 1
  do iLat=1,nLats
     do iLon=1,nLons
        write(iOutputUnit_)       &
             Longitude(iLon,iBlock), &
             Latitude(iLat,iBlock),&
             Altitude_GB(iLon,iLat,iAlt,iBlock), &
             Sza(iLon,iLat,iBlock), &
             VTEC(iLon,iLat,iBlock)
     enddo
  enddo

end subroutine output_2dtec

!-------------------------------------------------------------------------------
! AGB: Routine to output a 3D Mag file that includes the geomagnetic field info
!-------------------------------------------------------------------------------

subroutine output_3dmag(iBlock)

  use ModGITM
  use ModInputs

  implicit none

  integer, intent(in) :: iBlock
  integer :: iLat, iLon, iAlt

  do iAlt=-1,nAlts+2
     do iLat=-1,nLats+2
        do iLon=-1,nLons+2
           write(iOutputUnit_)                      &
                Longitude(iLon,iBlock),             &
                Latitude(iLat,iBlock),              &
                Altitude_GB(iLon,iLat,iAlt,iBlock), &
                mLatitude(iLon,iLat,iAlt,iBlock),   &
                mLongitude(iLon,iLat,iAlt,iBlock),  &
                B0(iLon,iLat,iAlt,:,iBlock)
           ! B0 has 4 elements of the geomag field, the E, N, Up, and magnitude
        enddo
     enddo
  end do

end subroutine output_3dmag

!----------------------------------------------------------------

subroutine output_2dmel(iBlock)
  
  use ModElectrodynamics
  use ModConstants, only:Pi
  use ModGITM
  use ModInputs

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon
  !--------------------------------------------------------------------------

  iAlt = 1
  do iLat=1,nMagLats
     do iLon=1,nMagLons+1
        write(iOutputUnit_)       &
             MagLonMC(iLon,iLat)*Pi/180.0, &
             MagLatMC(iLon,iLat)*Pi/180.0, &
             Altitude_GB(iLon,iLat,iAlt,iBlock), &
             MagLocTimeMC(iLon,iLat)*Pi/180.0, &
             GeoLatMC(iLon,iLat), &
             GeoLonMC(iLon,iLat), &
             SigmaPedersenMC(iLon,iLat), &
             SigmaHallMC(iLon,iLat), &
             DivJuAltMC(iLon,iLat), &
             LengthMC(iLon,iLat), &
             SigmaPPMC(iLon,iLat), &
             SigmaLLMC(iLon,iLat), &
             SigmaHHMC(iLon,iLat), &
             SigmaCCMC(iLon,iLat), &
             SigmaPLMC(iLon,iLat), &
             SigmaLPMC(iLon,iLat), &
             KDpmMC(iLon,iLat), &
             KdlmMC(iLon,iLat), &
             solver_a_mc(iLon,iLat), &
             solver_b_mc(iLon,iLat), &
             solver_c_mc(iLon,iLat), &
             solver_d_mc(iLon,iLat), &
             solver_e_mc(iLon,iLat), &
             solver_s_mc(iLon,iLat), &
             DynamoPotentialMC(iLon,iLat), &
             Ed1new(iLon,iLat), &
             Ed2new(iLon,iLat), &
             kpmMC(iLon,iLat), &
             klmMC(iLon,iLat)
     enddo
  enddo

end subroutine output_2dmel


!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_1dall(iiLon, iiLat, iBlock, rLon, rLat, iUnit)

  use ModGITM
  use ModEUV, only: HeatingEfficiency_CB
  use ModSources, only: JouleHeating, RadCooling, EuvHeating, Conduction
                        
  use ModInputs, only: iOutputUnit_
  implicit none

  integer, intent(in) :: iiLat, iiLon, iBlock, iUnit
  real, intent(in)    :: rLon, rLat

  integer, parameter :: nVars = 13+nSpeciesTotal+nSpecies+nIons+nSpecies+5
  real :: Vars(nVars)
  real :: Tmp(0:nLons+1,0:nLats+1)
  integer :: iAlt, iiAlt, iOff, iIon, iSpecies, iDir

  do iAlt=-1,nAlts+2

     iiAlt = max(min(iAlt,nAlts),1)


     Vars(1) = &
          rLon*Longitude(iiLon,iBlock)+(1-rLon)*Longitude(iiLon+1,iBlock)
     Vars(2) = &
          rLat*Latitude(iiLat,iBlock)+(1-rLat)*Latitude(iiLat+1,iBlock)

     Vars(3) = Altitude_GB(iiLon, iiLat, iAlt, iBlock)

     Tmp = Rho(0:nLons+1,0:nLats+1,iAlt,iBlock)
     Vars(4) = inter(Tmp,iiLon,iiLat,rlon,rlat)

     iOff = 4
     do iSpecies = 1, nSpeciesTotal
        Tmp = NDensityS(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)
!        Tmp = NDensityS(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)/NDensity(0:nLons+1,0:nLats+1,iAlt,iBlock)
        Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     Tmp = Temperature(0:nLons+1,0:nLats+1,iAlt,iBlock) * &
          TempUnit(0:nLons+1,0:nLats+1,iAlt)
     iOff = 5+nSpeciesTotal
     Vars(iOff) = inter(Tmp,iiLon,iiLat,rlon,rlat)

     do iDir = 1, 3
        Tmp = Velocity(0:nLons+1,0:nLats+1,iAlt,iDir,iBlock)
        Vars(iOff+iDir) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     iOff = 8+nSpeciesTotal
     do iSpecies = 1, nSpecies
        Tmp = VerticalVelocity(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)
        Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     iOff = 8+nSpeciesTotal+nSpecies
     do iIon = 1, nIons
        Tmp = IDensityS(0:nLons+1,0:nLats+1,iAlt,iIon,iBlock)
        Vars(iOff+iIon) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     iOff = 8+nSpeciesTotal+nSpecies+nIons+1
     Tmp = eTemperature(0:nLons+1,0:nLats+1,iAlt,iBlock)
     Vars(iOff)   = inter(Tmp,iiLon,iiLat,rlon,rlat)

     iOff = iOff+1
     Tmp = iTemperature(0:nLons+1,0:nLats+1,iAlt,iBlock)
     Vars(iOff) = inter(Tmp,iiLon,iiLat,rlon,rlat)

!     do iDir = 1, 3
!        Tmp = IVelocity(0:nLons+1,0:nLats+1,iAlt,iDir,iBlock)
!        Vars(iOff+iDir) = inter(Tmp,iiLon,iiLat,rlon,rlat)
!     enddo

        iOff = iOff + 1
        Tmp = IVelocity(0:nLons+1,0:nLats+1,iAlt,iEast_,iBlock)
        Vars(iOff) = inter(Tmp,iiLon,iiLat,rlon,rlat)

        iOff = iOff + 1
        Tmp = IVelocity(0:nLons+1,0:nLats+1,iAlt,iNorth_,iBlock)
        Vars(iOff) = inter(Tmp,iiLon,iiLat,rlon,rlat)

        iOff = iOff + 1
        Tmp = IVelocity(0:nLons+1,0:nLats+1,iAlt,iUp_,iBlock)
        Vars(iOff) = inter(Tmp,iiLon,iiLat,rlon,rlat)

        do iSpecies = 1, nSpecies
           Tmp = NDensityS(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)/NDensity(0:nLons+1,0:nLats+1,iAlt,iBlock)
           Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,rlon,rlat)
        enddo

!        iOff = iOff + 1
!        Tmp = NDensityS(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)/NDensity(0:nLons+1,0:nLats+1,iAlt,iBlock)
!        Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,rlon,rlat)
! AGB: Fixed this to use interpolation like all the other variables.
!      Also the porper lat and lon
        iOff = iOff + nSpecies
        Tmp = Dt*RadCooling(0:nLons+1,0:nLats+1,iAlt,iBlock) * &
             TempUnit(0:nLons+1,0:nLats+1,iAlt)
        Vars(iOff+1) = inter(Tmp,iiLon,iiLat,rlon,rlat)

        Tmp = Dt*EuvHeating(0:nLons+1,0:nLats+1,iAlt,iBlock) * &
             TempUnit(0:nLons+1,0:nLats+1,iAlt)
        Vars(iOff+2) = inter(Tmp,iiLon,iiLat,rlon,rlat)

        Tmp(1:nLons,1:nLats) = Conduction(1:nLons,1:nLats,iAlt) * &
             TempUnit(1:nLons,1:nLats,iAlt)
        Vars(iOff+3) = inter(Tmp,iiLon,iiLat,rlon,rlat)
 
        Vars(iOff+4) = Vars(iOff+2) - Vars(iOff+1) + Vars(iOff+3)

        Tmp = HeatingEfficiency_CB(0:nLons+1,0:nLats+1,iAlt,iBlock)
        Vars(iOff+5) = inter(Tmp,iiLon,iiLat,rlon,rlat)
! AGB: End of corrections

     write(iOutputUnit_) Vars

  enddo

contains

  real function inter(variable, iiLon, iiLat, rLon, rLat) result(PointValue)

    implicit none

    real :: variable(0:nLons+1, 0:nLats+1), rLon, rLat
    integer :: iiLon, iiLat

    PointValue = &  
          (  rLon)*(  rLat)*Variable(iiLon  ,iiLat  ) + &
          (1-rLon)*(  rLat)*Variable(iiLon+1,iiLat  ) + &
          (  rLon)*(1-rLat)*Variable(iiLon  ,iiLat+1) + &
          (1-rLon)*(1-rLat)*Variable(iiLon+1,iiLat+1)

  end function inter

end subroutine output_1dall

!----------------------------------------------------------------
! output_0dall: outputs GITM data at a specified 3D satellite location
!----------------------------------------------------------------

subroutine output_0dall(iiLon, iiLat, iiAlt, iBlock, rLon, rLat, rAlt, iUnit)

  use ModGITM
  use ModEUV, only: HeatingEfficiency_CB
  use ModSources, only: JouleHeating, RadCooling, EuvHeating, Conduction
  use ModMpi
                        
  use ModInputs, only: iOutputUnit_
  implicit none

  integer, intent(in) :: iiLat, iiLon, iiAlt, iBlock, iUnit
  real, intent(in)    :: rLon, rLat, rAlt

  integer :: ierr
  integer, parameter :: nVars = 13+nSpeciesTotal+nSpecies+nIons+nSpecies+5
  real :: Vars(nVars)
  real :: Tmp(0:nLons+1,0:nLats+1,0:nAlts+1)
  integer :: jAlt, iOff, iIon, iSpecies, iDir

  jAlt = max(min(iiAlt,nAlts),1)

  !Determine the satellite position using linear interpolation
  Vars(1) = rLon*Longitude(iiLon,iBlock)+(1-rLon)*Longitude(iiLon+1,iBlock)
  Vars(2) = rLat*Latitude(iiLat,iBlock)+(1-rLat)*Latitude(iiLat+1,iBlock)
  Vars(3) = rAlt*Altitude_GB(iiLon, iiLat, iiAlt, iBlock) + &
       (1-rAlt)*Altitude_GB(iiLon+1, iiLat+1, iiAlt+1, iBlock)

  ! Get the species characteristics at this location
  Tmp     = Rho(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock)
  Vars(4) = inter(Tmp,iiLon,iiLat,iiAlt,rLon,rLat,rAlt)

  iOff = 4
  do iSpecies = 1, nSpeciesTotal
     Tmp = NDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iSpecies,iBlock)
     Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,iiAlt,rLon,rLat,rAlt)
  enddo

  Tmp = Temperature(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock) * &
       TempUnit(0:nLons+1,0:nLats+1,0:nAlts+1)
  iOff = 5+nSpeciesTotal
  Vars(iOff) = inter(Tmp,iiLon,iiLat,iiAlt,rLon,rLat,rAlt)

  do iDir = 1, 3
     Tmp = Velocity(0:nLons+1,0:nLats+1,0:nAlts+1,iDir,iBlock)
     Vars(iOff+iDir) = inter(Tmp,iiLon,iiLat,iiAlt,rLon,rLat,rAlt)
  enddo

  iOff = 8+nSpeciesTotal
  do iSpecies = 1, nSpecies
     Tmp = VerticalVelocity(0:nLons+1,0:nLats+1,0:nAlts+1,iSpecies,iBlock)
     Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,iiAlt,rLon,rLat,rAlt)
  enddo

  iOff = 8+nSpeciesTotal+nSpecies
  do iIon = 1, nIons
     Tmp = IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iIon,iBlock)
     Vars(iOff+iIon) = inter(Tmp,iiLon,iiLat,iiAlt,rLon,rLat,rAlt)
  enddo

  iOff = 8+nSpeciesTotal+nSpecies+nIons+1
  Tmp = eTemperature(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock)
  Vars(iOff) = inter(Tmp,iiLon,iiLat,iiAlt,rLon,rLat,rAlt)

  iOff = iOff+1
  Tmp = iTemperature(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock)
  Vars(iOff) = inter(Tmp,iiLon,iiLat,iiAlt,rLon,rLat,rAlt)

  iOff = iOff + 1
  Tmp = IVelocity(0:nLons+1,0:nLats+1,0:nAlts+1,iEast_,iBlock)
  Vars(iOff) = inter(Tmp,iiLon,iiLat,iiAlt,rLon,rLat,rAlt)

  iOff = iOff + 1
  Tmp = IVelocity(0:nLons+1,0:nLats+1,0:nAlts+1,iNorth_,iBlock)
  Vars(iOff) = inter(Tmp,iiLon,iiLat,iiAlt,rLon,rLat,rAlt)

  iOff = iOff + 1
  Tmp = IVelocity(0:nLons+1,0:nLats+1,0:nAlts+1,iUp_,iBlock)
  Vars(iOff) = inter(Tmp,iiLon,iiLat,iiAlt,rLon,rLat,rAlt)

  do iSpecies = 1, nSpecies
     Tmp = NDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,iSpecies,iBlock) &
          / NDensity(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock)
     Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,iiAlt,rLon,rLat,rAlt)
  enddo

  iOff = iOff + nSpecies
  Vars(iOff+1) = Dt*RadCooling(1,1,jAlt,iBlock)*TempUnit(1,1,jAlt)

  Vars(iOff+2) = Dt*EuvHeating(1,1,jAlt,iBlock)*TempUnit(1,1,jAlt)

  Vars(iOff+3) = Conduction(1,1,jAlt)*TempUnit(1,1,jAlt)

  Vars(iOff+4) = Dt*EuvHeating(1,1,jAlt,iBlock)*TempUnit(1,1,jAlt) - &
       Dt*RadCooling(1,1,jAlt,iBlock)*TempUnit(1,1,jAlt) + &
       Conduction(1,1,jAlt)*TempUnit(1,1,jAlt)

  Vars(iOff+5) = HeatingEfficiency_CB(1,1,jAlt,iBlock)

  ! Write the output data
  write(iOutputUnit_) Vars

contains

  real function inter(variable, iiLon, iiLat, iiAlt, rLon, rLat, rAlt) &
       result(PointValue)

    implicit none

    real :: variable(0:nLons+1, 0:nLats+1, 0:nAlts+1), rLon, rLat, rAlt
    integer :: iiLon, iiLat, iiAlt

    PointValue = &  
          (  rLon)*(  rLat)*(  rAlt)*Variable(iiLon  ,iiLat  ,iiAlt  ) + &
          (1-rLon)*(  rLat)*(  rAlt)*Variable(iiLon+1,iiLat  ,iiAlt  ) + &
          (  rLon)*(1-rLat)*(  rAlt)*Variable(iiLon  ,iiLat+1,iiAlt  ) + &
          (1-rLon)*(1-rLat)*(  rAlt)*Variable(iiLon+1,iiLat+1,iiAlt  ) + &
          (  rLon)*(  rLat)*(1-rAlt)*Variable(iiLon  ,iiLat  ,iiAlt+1) + &
          (1-rLon)*(  rLat)*(1-rAlt)*Variable(iiLon+1,iiLat  ,iiAlt+1) + &
          (  rLon)*(1-rLat)*(1-rAlt)*Variable(iiLon  ,iiLat+1,iiAlt+1) + &
          (1-rLon)*(1-rLat)*(1-rAlt)*Variable(iiLon+1,iiLat+1,iiAlt+1)

  end function inter
end subroutine output_0dall

subroutine output_1dnew(iiLon, iiLat, iBlock, rLon, rLat, iUnit)

  use ModGITM
  use ModEUV, only : Sza, HeatingEfficiency_CB
  use ModSources, only: JouleHeating
  use ModInputs, only: iOutputUnit_
  implicit none

  integer, intent(in) :: iiLat, iiLon, iBlock, iUnit
  real, intent(in)    :: rLon, rLat

  integer, parameter :: nVars = 13+nSpeciesTotal+nSpecies+nIons+nSpecies+2
  real :: Vars(nVars)
  real :: Tmp(0:nLons+1,0:nLats+1)
  integer :: iAlt, iiAlt, iOff, iIon, iSpecies, iDir

  do iAlt=-1,nAlts+2

     iiAlt = max(min(iAlt,nAlts),1)

     Vars(1) = &
          rLon*Longitude(iiLon,iBlock)+(1-rLon)*Longitude(iiLon+1,iBlock)
     Vars(2) = &
          rLon*LocalTime(iiLon)+(1-rLon)*LocalTime(iiLon+1)
     Vars(3) = &
          rLat*Latitude(iiLat,iBlock)+(1-rLat)*Latitude(iiLat+1,iBlock)

     Vars(4) = Sza(iiLon,iiLat,iBlock)

     Vars(5) = Altitude_GB(iiLon, iiLat, iAlt, iBlock)

     Tmp = Rho(0:nLons+1,0:nLats+1,iAlt,iBlock)
     Vars(6) = inter(Tmp,iiLon,iiLat,rlon,rlat)

     iOff = 6
     do iSpecies = 1, nSpeciesTotal
        Tmp = NDensityS(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)
        Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     Tmp = Temperature(0:nLons+1,0:nLats+1,iAlt,iBlock) * &
          TempUnit(0:nLons+1,0:nLats+1,iAlt)
     iOff = 7+nSpeciesTotal
     Vars(iOff) = inter(Tmp,iiLon,iiLat,rlon,rlat)

     do iDir = 1, 3
        Tmp = Velocity(0:nLons+1,0:nLats+1,iAlt,iDir,iBlock)
        Vars(iOff+iDir) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     iOff = 10+nSpeciesTotal
     do iSpecies = 1, nSpecies
        Tmp = VerticalVelocity(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)
        Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     iOff = 10+nSpeciesTotal+nSpecies
     do iIon = 1, nIons
        Tmp = IDensityS(0:nLons+1,0:nLats+1,iAlt,iIon,iBlock)
        Vars(iOff+iIon) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     iOff = iOff+1
     Tmp = eTemperature(0:nLons+1,0:nLats+1,iAlt,iBlock)
     Vars(iOff)   = inter(Tmp,iiLon,iiLat,rlon,rlat)

     iOff = iOff+1
     Tmp = iTemperature(0:nLons+1,0:nLats+1,iAlt,iBlock)
     Vars(iOff) = inter(Tmp,iiLon,iiLat,rlon,rlat)

     iOff = iOff
     do iDir = 1, 3
        Tmp = IVelocity(0:nLons+1,0:nLats+1,iAlt,iDir,iBlock)
        Vars(iOff+iDir) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     iOff = iOff
     do iSpecies = 1, nSpecies
        Tmp = NDensityS(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)/NDensity(0:nLons+1,0:nLats+1,iAlt,iBlock)
        Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     write(iOutputUnit_) Vars

  enddo

contains

  real function inter(variable, iiLon, iiLat, rLon, rLat) result(PointValue)

    implicit none

    real :: variable(0:nLons+1, 0:nLats+1), rLon, rLat
    integer :: iiLon, iiLat

    PointValue = &  
          (  rLon)*(  rLat)*Variable(iiLon  ,iiLat  ) + &
          (1-rLon)*(  rLat)*Variable(iiLon+1,iiLat  ) + &
          (  rLon)*(1-rLat)*Variable(iiLon  ,iiLat+1) + &
          (1-rLon)*(1-rLat)*Variable(iiLon+1,iiLat+1)

  end function inter

end subroutine output_1dnew



