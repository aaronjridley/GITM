!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine readMHDoutput(iBLK, iError)

  use ModMHD_Interface
  use ModEIEFiles
  use ModTimeConvert, ONLY: time_int_to_real
  use ModIoUnit, ONLY : UNITTMP_
  
  implicit none

  integer, intent(out) :: iError
  integer, intent(in)  :: iBLK

  integer :: iTime
  integer :: nfields
  integer :: ntemp, iyr, imo, ida, ihr, imi
  integer :: i,j, iField, iPot_, iAveE_, iEFlux_
  real*4  :: swv,bx,by,bz,aei,ae,au,al,dsti,dst,hpi,sjh,pot
  real*8  :: rtime
  integer, dimension(7) :: itime_i

  real*4, allocatable, dimension(:,:,:) :: AllData

  integer, parameter :: nFieldsMax = 100
  character (len=30), dimension(nFieldsMax) :: Fields

  logical :: IsBinary

  real :: dPotential
  integer :: nCellsPad

  integer :: strpos
  character (len=100) :: FileName, cLine
  character (len=6)   :: cDate
  character (len=6)   :: cTime
  logical :: IsThere, IsDone
  logical :: IsFile(1440)
  character (len=100) :: FileNameList(1440)
  integer :: iFirstFile, iMLT, iLat

  real*4, allocatable,dimension(:,:,:,:) :: SigmaP, SigmaH

  iError = 0

  strpos = index(MHD_FileName,"*")

  iyr = mod(MHD_Year,100)
  imo = MHD_Month
  ida = MHD_Day

  write(cDate,'(3i2.2)') iyr,imo,ida

  iFirstFile = -1

  do ihr=0,23
     do imi=0,59

        write(cTime,'(3i2.2)') ihr,imi,0
        FileName = MHD_FileName(1:strpos-1)//cDate//"_"//cTime//".idl"
        inquire(file=FileName,EXIST=IsThere)

        if (IsThere) write(*,*) FileName,'--> ',IsThere  
        IsFile(ihr*60+imi+1) = IsThere
        if (IsThere) then
           if (iFirstFile < 0) iFirstFile = ihr*60+imi+1
           MHD_ntimes = MHD_ntimes+1
           FileNameList(ihr*60+imi+1) = FileName
        endif

     enddo
  enddo

  if (iFirstFile < 0) then
     write(*,*) "Could not find any MHD files which you specified!!!"
     stop
  endif

  FileName = FileNameList(iFirstFile)

  open(UnitTmp_, file=FileName, status='old',iostat=iError)
  if (iError.ne.0) then
     write(*,*) "Error opening file:", FileName
     stop
  endif
  MHD_nLats = 0

  IsDone = .false.

  i = 1

  do while (.not.IsDone)

     read(UnitTmp_,*) cLine
     write(*,*) cLine

     if (i > 100) IsDone = .true.
     i =i + 1

     if (index(cLine,"NUMERICAL") > 0) then

        read(UnitTmp_,*) nFields
        read(UnitTmp_,*) MHD_nLats
        read(UnitTmp_,*) MHD_nMlts

        write(*,*) nFields, MHD_nLats, MHD_nMlts

        if (allocated(MHD_Lats)) deallocate(MHD_Lats)
        allocate(MHD_Lats(MHD_nLats), stat=iError)
        if (iError /= 0) then
           write(*,*) "Error in allocating array MHD_Lats in "
           stop
        endif

        if (allocated(MHD_Mlts)) deallocate(MHD_Mlts)
        allocate(MHD_Mlts(MHD_nMlts), stat=iError)
        if (iError /= 0) then
           write(*,*) "Error in allocating array Mlts in "
           stop
        endif

        if (.not.allocated(MHD_Potential)) then

           allocate(MHD_Potential(MHD_nMlts,MHD_nLats, &
                MHD_nTimes,2), stat=iError)
           if (iError /= 0) then
              write(*,*) "Error in allocating array MHD_Potential in "
              stop
           endif

           allocate(MHD_EFlux(MHD_nMlts,MHD_nLats, &
                MHD_nTimes,2), stat=iError)
           if (iError /= 0) then
              write(*,*) "Error in allocating array MHD_EFlux in "
              stop
           endif

           allocate(MHD_AveE(MHD_nMlts,MHD_nLats, &
                MHD_nTimes,2), stat=iError)
           if (iError /= 0) then
              write(*,*) "Error in allocating array MHD_AveE in "
              stop
           endif

           allocate(MHD_Value(MHD_nMlts,MHD_nLats, &
                MHD_nTimes,2), stat=iError)
           if (iError /= 0) then
              write(*,*) "Error in allocating array MHD_Value in "
              stop
           endif

           allocate(SigmaP(MHD_nMlts,MHD_nLats, &
                MHD_nTimes,2), stat=iError)
           if (iError /= 0) then
              write(*,*) "Error in allocating array SigmaP in "
              stop
           endif

           allocate(SigmaH(MHD_nMlts,MHD_nLats, &
                MHD_nTimes,2), stat=iError)
           if (iError /= 0) then
              write(*,*) "Error in allocating array SigmaH in "
              stop
           endif

           allocate(MHD_Time(MHD_nTimes,2), stat=iError)
           if (iError /= 0) then
              write(*,*) "Error in allocating array MHDTimes in "
              stop
           endif

        endif

     endif

     if (index(cLine,"BEGIN") > 0) then

        MHD_Mlts = -1.0e32

        do iMLT = 1, MHD_nMlts
           i = mod(iMLT+(MHD_nMlts-1)/2-1, (MHD_nMlts-1))+1
           do iLat = 1, MHD_nLats
              read(UnitTmp_,*) MHD_Lats(iLat), MHD_Mlts(i)
           enddo
        enddo

        IsDone = .true.

     endif

  end do

  MHD_Mlts(MHD_nMLTs) = MHD_Mlts(1)

  MHD_Mlts = MHD_Mlts * 24.0 / 360.0 + 12.0
  where (MHD_Mlts > 24) MHD_Mlts = MHD_Mlts - 24.0

  do iMLT = 1, MHD_nMlts-1
     if (MHD_Mlts(iMLT)-12 > MHD_Mlts(iMLT+1))  &
          MHD_Mlts(iMLT) = MHD_Mlts(iMLT)-24
  enddo

  MHD_Lats = 90.0 - MHD_Lats

  if (nFields > nFieldsMax) then
     write(*,*) "Maximum number of fields in MHD is ",nFieldsMax
     stop
  endif

  allocate(AllData(MHD_nMlts,MHD_nLats,nFields), stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array AllData in "
     stop
  endif

  MHD_iDebugLevel = 0

  do iField=1,nfields
     if (IsBinary) then
        read(UnitTmp_) Fields(iField)
     else
        read(UnitTmp_,'(a)') Fields(iField)
     endif

     if (MHD_iDebugLevel > 1) write(*,*) Fields(iField)

     if ((index(Fields(iField),"Potential") > 0).and. &
         (index(Fields(iField),"odel") < 1)) then
        iPot_ = iField
        if (MHD_iDebugLevel > 1) write(*,*) "<--- Potential Found", iPot_
     endif

     if ((index(Fields(iField),"Mean Energy") > 0) .and. &
         (index(Fields(iField),"odel") < 1)) then
        iAveE_ = iField
        if (MHD_iDebugLevel > 1) write(*,*) "<--- Mean Energy Found", iAveE_
     endif

     if ((index(Fields(iField),"Energy Flux") > 0) .and. &
         (index(Fields(iField),"odel") < 1)) then
        iEFlux_ = iField
        if (MHD_iDebugLevel > 1) write(*,*) "<--- Energy Flux Found", iEFlux_
     endif

  enddo

  do iTime=1,MHD_ntimes

     if (IsBinary) then

        read(UnitTmp_) ntemp,iyr,imo,ida,ihr,imi
        read(UnitTmp_) swv,bx,by,bz,aei,ae,au,al,dsti,dst,hpi,sjh,pot

        do iField=1,nfields
           read(UnitTmp_) ((AllData(j,i,iField),j=1,MHD_nMlts),i=1,MHD_nLats)
        enddo

     else

        read(UnitTmp_,*) ntemp,iyr,imo,ida,ihr,imi
        read(UnitTmp_,*) swv,bx,by,bz,aei,ae,au,al,dsti,dst,hpi,sjh,pot

        do iField=1,nfields
           read(UnitTmp_,*) ((AllData(j,i,iField),j=1,MHD_nMlts),i=1,MHD_nLats)
        enddo

     endif

     itime_i(1) = iyr
     itime_i(2) = imo
     itime_i(3) = ida
     itime_i(4) = ihr
     itime_i(5) = imi
     itime_i(6) = 0
     itime_i(7) = 0
     call time_int_to_real(itime_i,rtime)
     MHD_Time(iTime,iBLK) = rtime

     ! We need Potential to be in Volts
     !         AveE to be in keV
     !         EFlux to be in W/m2

     MHD_Potential(:,1:MHD_nLats,iTime,iBLK) = AllData(:,1:MHD_nLats,iPot_)
     MHD_AveE(:,1:MHD_nLats,iTime,iBLK)      = AllData(:,1:MHD_nLats,iAveE_)
     ! Need to convert from erg/cm2/s to W/m2
     MHD_EFlux(:,1:MHD_nLats,iTime,iBLK)     = &
          AllData(:,1:MHD_nLats,iEFlux_) !* 1.0e-7 * 100.0 * 100.0

     do i=1,MHD_nMlts

        dPotential = MHD_Potential(i,MHD_nLats,iTime,iBLK)/nCellsPad

        do j=MHD_nLats+1, MHD_nLats+nCellsPad

           MHD_AveE(i,j,iTime,iBLK)  = MHD_AveE(i,MHD_nLats,iTime,iBLK)
           MHD_EFlux(i,j,iTime,iBLK) = MHD_EFlux(i,MHD_nLats,iTime,iBLK)

           MHD_Potential(i,j,iTime,iBLK) = &
                MHD_Potential(i,j-1,iTime,iBLK) - dPotential

        enddo
     enddo

  enddo

  MHD_nLats = MHD_nLats + nCellsPad

  close(UnitTmp_)

  deallocate(AllData, stat=iError)

end subroutine readMHDoutput
