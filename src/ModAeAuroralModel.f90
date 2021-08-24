! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

module ModAeAuroralModel

  implicit none

  integer, parameter :: nMlts = 96
  integer, parameter :: nMLats = 19
  integer, parameter :: nLevels = 7

  integer, parameter :: MaxAe = 350
  integer, parameter :: AeBinSize = 50
  
  real, dimension(nMLats) :: mlts
  real, dimension(nLevels, nMlts, nMLats) :: mLats, eFlux, AveE

  character (len=*), parameter :: cDirectory = "UA/DataIn/Aurora/"

  character (len=*), parameter :: cFileStart = "AE_"
  character (len=*), parameter :: cFileEnd = ".old.dat"
  
contains

  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------

  subroutine run_ae_model(GitmCurrentTime, iBlock)

    use ModSizeGITM
    use ModGITM
    use ModInputs

    use ModKind, ONLY: Real8_
    use ModIndicesInterfaces

    real(Real8_), intent(in) :: GitmCurrentTime
    integer, intent(in) :: iBlock

    logical :: IsFirstTime=.true.

    integer :: iError, iAe, iMLat, i, iLat, iLon, iMlt
    real :: ae, gMlat, mindist
    real :: dist(nMLats)
    
    call start_timing("run_ae_model")

    if (IsFirstTime) then
       iError = 0
       IsFirstTime = .false.
    endif
    
    call get_ae(GitmCurrentTime+TimeDelayHighLat, ae, iError)

    if (ae > MaxAe) ae = MaxAe
    iAe = ae/AeBinSize

    ElectronEnergyFlux = 0.01
    ElectronAverageEnergy = 0.1

    do iLon = -1, nLons+2
       do iLat = -1, nLats+2

          iMlt = mod(floor(mod(MLT(iLon, iLat, nAlts+1)+24.0,24.0)*4),nMlts)
          if (iMlt == 0) iMlt = nMlts

          gMlat = abs(MLatitude(iLon, iLat, nAlts+1, iBlock))
          
          if (gMlat >= mLats(iAe,iMlt,1) .and.   &
             gMlat <= mLats(iAe,iMlt,nMlats)) then
             
             dist = abs(gMlat - mLats(iAe,iMlt,:))
             mindist = minval(dist)
             ! use minloc?
             iMLat = 1
             do i = 1,nMlats
                if (dist(i) == mindist) iMLat = i
             enddo             

             ElectronEnergyFlux(iLon, iLat) = eFlux(iAe, iMlt, iMlat)
             ElectronAverageEnergy(iLon, iLat) = AveE(iAe, iMlt, iMlat)

          endif

       enddo

    enddo

    call end_timing("run_ae_model")

  end subroutine run_ae_model

  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------

  subroutine read_ae_model_files(iError)

    use ModInputs
    
    integer, intent(out) :: iError

    integer :: ae, iAe, iLat, iMlt

    character (len=3) :: cAe1, cAe2
    character (len=iCharLen_) :: cFile
    character (len=iCharLen_) :: cDummyLine
    logical :: IsThere
    
    do ae = 0, MaxAe, AeBinSize

       iAe = ae/AeBinSize
       call i2s(ae, cAe1, 3)
       call i2s(ae+AeBinSize, cAe2, 3)
       cFile = cDirectory//cFileStart//cAe1//"_"//cAe2//cFileEnd

       inquire(file=cFile,EXIST=IsThere)
       if (.not.IsThere) then
          write(*,*) cFile//" cannot be found by read_ovationsm_files"
          call stop_gitm("must stop!!!")
       endif

       iError = 0

       open(iInputUnit_,file=cFile,status="old", iostat=iError)

       ! Read in the header:
       read(iInputUnit_,*,iostat=iError) cDummyLine

       ! Read in mlts, only once:
       read(iInputUnit_,*,iostat=iError) mlts

       ! Read in mlats for all bins:
       do iLat = 1, nMLats 
          read(iInputUnit_,*,iostat=iError) mlats(iAe, :, iLat)
       enddo

       ! Read in eflux for all bins:
       do iLat = 1, nMLats 
          read(iInputUnit_,*,iostat=iError) eFlux(iAe, :, iLat)
       enddo

       ! Read in avee for all bins:
       do iLat = 1, nMLats 
          read(iInputUnit_,*,iostat=iError) AveE(iAe, :, iLat)
       enddo

    enddo
       
    close(iInputUnit_)

  end subroutine read_ae_model_files

end module ModAeAuroralModel
