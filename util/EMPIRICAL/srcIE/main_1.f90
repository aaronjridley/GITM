!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program Interface

  use ModEIE_Interface
  use ModTimeConvert, ONLY: time_int_to_real

  implicit none

  character (len=100) :: inFileName
  integer             :: iError
  real*8              :: rTime
  integer, dimension(7) :: iTime_i
  real, dimension(4,2) :: templat, tempmlt, temppot

  iDebugLevel = 100

  write(6,*) 'Enter file name :'
  read(5,'(A100)') inFileName

  iError = 0

  call AMIE_SetFileName(inFileName)

  call readAMIEOutput(iError)

  call AMIE_GetnLats(EIEi_HavenLats)
  call AMIE_GetnMLTs(EIEi_HavenMLTs)
  EIEi_HavenBLKs = 2

  if (iDebugLevel > 1) then
     write(*,*) "EIEi_HavenBLKs : ", EIEi_HavenBLKs
     write(*,*) "EIEi_HavenLats : ", EIEi_HavenLats
     write(*,*) "EIEi_HavenMLTs : ", EIEi_HavenMLTs
  endif

  allocate(EIEr3_HaveLats(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs), &
       stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array EIEr3_HaveLats in Interface"
     stop
  endif

  allocate(EIEr3_HaveMlts(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs), &
       stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array EIEr3_HaveMlts in Interface"
     stop
  endif

  allocate(EIEr3_HavePotential(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs), &
       stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array EIEr3_HavePotential in Interface"
     stop
  endif

  allocate(EIEr3_HaveEFlux(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs), &
       stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array EIEr3_HaveEFlux in Interface"
     stop
  endif

  allocate(EIEr3_HaveAveE(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs), &
       stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array EIEr3_HaveAveE in Interface"
     stop
  endif

  call AMIE_GetLats(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs,&
       EIEr3_HaveLats,iError)

  call AMIE_GetMLTs(EIEi_HavenMlts,EIEi_HavenLats,EIEi_HavenBLKs,&
       EIEr3_HaveMLTs,iError)

  itime_i(1) = 1998
  itime_i(2) = 05
  itime_i(3) = 01
  itime_i(4) = 12
  itime_i(5) = 00
  itime_i(6) = 0
  itime_i(7) = 0
  call time_int_to_real(itime_i,rtime)

  call AMIE_GetPotential(rtime, EIE_Interpolate_, &
       EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs, EIEr3_HavePotential, iError)

  call AMIE_GetAveE(rtime, EIE_Closest_, &
       EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs, EIEr3_HaveAveE, iError)

  call AMIE_GetEFlux(rtime, EIE_Closest_, &
       EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs, EIEr3_HaveEFlux, iError)

  call IO_SetnMLTs(4)
  call IO_SetnLats(2)

  templat(:,1) = -60.0
  templat(:,2) = 70.0
  tempmlt(1,:) =  0.0
  tempmlt(2,:) =  6.0
  tempmlt(3,:) = 12.0
  tempmlt(4,:) = 18.0

  call IO_SetGrid(tempmlt, templat, iError)

  call IO_GetPotential(temppot,iError)

  write(*,*) templat(:,2)
  write(*,*) tempmlt(:,2)
  write(*,*) temppot(:,2)

  call EIE_End

end program Interface
