!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModNewell

  implicit none

  integer, parameter        :: iCharLen_     = 100

  character (len=iCharLen_) :: dir="UA/DataIn/Aurora/"

  character (len=iCharLen_) :: cFileDiffef="diff.txt"
  character (len=iCharLen_) :: cFileDiffnf="diff_n.txt"
  character (len=iCharLen_) :: cFileDiffp ="prob_b_diff.txt"
  character (len=iCharLen_) :: cFileMonoef="mono.txt"
  character (len=iCharLen_) :: cFileMononf="mono_n.txt"
  character (len=iCharLen_) :: cFileMonop ="prob_b_mono.txt"
  character (len=iCharLen_) :: cFileWaveef="wave.txt"
  character (len=iCharLen_) :: cFileWavenf="wave_n.txt"
  character (len=iCharLen_) :: cFileWavep ="prob_b_wave.txt"
  character (len=iCharLen_) :: cFileIonsef="ions.txt"
  character (len=iCharLen_) :: cFileIonsnf="ions_n.txt"

  integer, parameter :: nMlts  = 96
  integer, parameter :: nMLats = 160
  integer, parameter :: nProbs = 3
  integer, parameter :: ndF    = 12

  real, dimension(nMlts, nMlats) :: &
       B1aDiff,B2aDiff,B1aDiffn,B2aDiffn,rFaDiff,rFaDiffn,B1pDiff,B2pDiff, &
       B1aMono,B2aMono,B1aMonon,B2aMonon,rFaMono,rFaMonon,B1pMono,B2pMono, &
       B1aWave,B2aWave,B1aWaven,B2aWaven,rFaWave,rFaWaven,B1pWave,B2pWave, &
       B1aIons,B2aIons,B1aIonsn,B2aIonsn,rFaIons,rFaIonsn, &
       ProbDiffTotal, ProbMonoTotal, ProbWaveTotal, ProbIonsTotal, &
       NumberFluxDiff, EnergyFluxDiff, &
       NumberFluxMono, EnergyFluxMono, &
       NumberFluxWave, EnergyFluxWave, &
       NumberFluxIons, EnergyFluxIons, &
       Area

  real, dimension(ndF, nMlts, nMlats) :: &
       ProbDiff, ProbMono, ProbWave

  real    :: dFdt = 1000.0
  integer :: dFBin

contains

  ! -------------------------------------------------------------------

  subroutine init_newell

    use ModConstants, only: pi
    integer :: iLat

    call report("Newell Aurora Initializing",2)

    call read_all_regression_files
    call read_all_probability_files

    ! This assumes about 500 km altitude and 1/2 deg and 1/4 hour resolution
    area = 120.0 * 120.0 * 0.5 * 3.75

    do iLat = 1,nMLats/2-1
       area(:,iLat)         = area(:,iLat) * cos((50.0 + 0.5*(iLat-1))*pi/180.0)
       area(:,iLat+nMLats/2)= area(:,iLat) * cos((50.0 + 0.5*(iLat-1))*pi/180.0)
    enddo

    call report("Newell Aurora Initialized",2)

  end subroutine init_newell

  ! -------------------------------------------------------------------

  subroutine run_newell(iBlock)

    use ModSizeGITM
    use ModGITM
    use ModInputs

    integer, intent(in) :: iBlock
    integer :: iLon, iLat, iMlat, iMlt, iMlat2
    real :: numflux, hps, hpn

    call start_timing("run_newell")

    ElectronEnergyFlux = 0.01
    ElectronAverageEnergy = 0.1

    if (iBlock == 1) then

       call calc_probability(ProbDiff, B1pDiff, B2pDiff, ProbDiffTotal)
       call calc_probability(ProbMono, B1pMono, B2pMono, ProbMonoTotal)
       call calc_probability(ProbWave, B1pWave, B2pWave, ProbWaveTotal)
       ProbIonsTotal = 1.0

       call calc_flux(ProbDiffTotal, B1aDiff,  B2aDiff,  EnergyFluxDiff)
       call calc_flux(ProbDiffTotal, B1aDiffn, B2aDiffn, NumberFluxDiff)
       
       call calc_flux(ProbMonoTotal, B1aMono,  B2aMono,  EnergyFluxMono)
       call calc_flux(ProbMonoTotal, B1aMonon, B2aMonon, NumberFluxMono)

!       if (iDebugLevel > -1) then
!          do iMlat = 1, nMLats
!             write(*,*) "from pat: ",iMlat,dfdt,&
!                  ProbMonoTotal(8,iMlat), B1aMono(8,iMlat), B2aMono(8,iMlat), &
!                  EnergyFluxMono(8,iMlat), NumberFluxMono(8,iMlat)
!          enddo
!       endif

       call calc_flux(ProbWaveTotal, B1aWave,  B2aWave,  EnergyFluxWave)
       call calc_flux(ProbWaveTotal, B1aWaven, B2aWaven, NumberFluxWave)
       
       call calc_flux(ProbIonsTotal, B1aIons,  B2aIons,  EnergyFluxIons)
       call calc_flux(ProbIonsTotal, B1aIonsn, B2aIonsn, NumberFluxIons)

       where (EnergyFluxDiff > 10.0) EnergyFluxDiff = 0.5
       where (EnergyFluxDiff > 5.0)  EnergyFluxDiff = 5.0
       where (EnergyFluxMono > 10.0) EnergyFluxMono = 0.5
       where (EnergyFluxMono > 5.0)  EnergyFluxMono = 5.0
       where (EnergyFluxWave > 10.0) EnergyFluxWave = 0.5
       where (EnergyFluxWave > 5.0)  EnergyFluxWave = 5.0

       where (EnergyFluxIons > 4.0)  EnergyFluxIons = 0.25
       where (EnergyFluxIons > 2.0)  EnergyFluxIons = 2.0

       where (NumberFluxDiff > 2.0e10) NumberFluxDiff = 0.0
       where (NumberFluxDiff > 2.0e9)  NumberFluxDiff = 1.0e9
       where (NumberFluxMono > 2.0e10) NumberFluxMono = 0.0
       where (NumberFluxMono > 2.0e9)  NumberFluxMono = 1.0e9
       where (NumberFluxWave > 2.0e10) NumberFluxWave = 0.0
       where (NumberFluxWave > 2.0e9)  NumberFluxWave = 1.0e9

       where (NumberFluxIons > 5.0e8)  NumberFluxIons = 0.0
       where (NumberFluxIons > 1.0e8)  NumberFluxIons = 1.0e8

       if (DoNewellRemoveSpikes .or. DoNewellAverage) then

          call calc_hp(EnergyFluxDiff, hps, hpn)
!          write(*,*) "Before Smooth, Diffuse : ", hps, hpn
          call smooth(EnergyFluxDiff)
          call smooth(NumberFluxDiff)
          call calc_hp(EnergyFluxDiff, hps, hpn)
!          write(*,*) "After Smooth, Diffuse : ", hps, hpn

          call calc_hp(EnergyFluxMono, hps, hpn)
!          write(*,*) "Before Smooth, Mono : ", hps, hpn
          call smooth(EnergyFluxMono)
          call smooth(NumberFluxMono)
          call calc_hp(EnergyFluxMono, hps, hpn)
!          write(*,*) "After Smooth, Mono : ", hps, hpn

          call calc_hp(EnergyFluxWave, hps, hpn)
!          write(*,*) "Before Smooth, Wave : ", hps, hpn
          call smooth(EnergyFluxWave)
          call smooth(NumberFluxWave)
!          write(*,*) "After Smooth, Wave : ", hps, hpn

          !call smooth(EnergyFluxIons)
          !call smooth(NumberFluxIons)

       endif

    endif

    do iLon = -1, nLons+2
       do iLat = -1, nLats+2
          if (abs(MLatitude(iLon, iLat, nAlts+1, iBlock)) > 50.0) then

             iMlat = floor(abs(MLatitude(iLon, iLat, nAlts+1, iBlock)) - 50.0)*2
             iMlat = min(max(iMlat,1),nMlats/2)
             if (MLatitude(iLon, iLat, nAlts+1, iBlock) > 0.0) &
                  iMlat = iMlat + nMlats/2
             iMlt = mod(floor(mod(MLT(iLon, iLat, nAlts+1)+24.0,24.0)*4),nMlts)
             if (iMlt == 0) iMlt = nMlts

             if (iMlat > nMlats/2) then
                iMlat2 = iMlat-nMlats/2
             else
                iMlat2 = iMlat+nMlats/2
             endif

             iMlat2 = min(max(iMlat2,1),nMlats)

             ! Diffuse Energy Flux

             if (UseNewellAveraged .or. EnergyFluxDiff(iMlt, iMlat)==0) then

                ! Add North and South together
                ElectronEnergyFlux(iLon, iLat) = &
                     EnergyFluxDiff(iMlt, iMlat) + EnergyFluxDiff(iMlt, iMlat2)

                ! If there are values in both hemisphere, then divide by 2
                if ( EnergyFluxDiff(iMlt, iMlat) * &
                     EnergyFluxDiff(iMlt, iMlat2) /= 0) &
                     ElectronEnergyFlux(iLon, iLat) = &
                     ElectronEnergyFlux(iLon, iLat)/2.0
             else
                ElectronEnergyFlux(iLon, iLat) = &
                     EnergyFluxDiff(iMlt, iMlat)
             endif

             ! Diffuse Number Flux

             if (UseNewellAveraged .or. NumberFluxDiff(iMlt, iMlat)==0) then

                ! Add North and South together
                numflux = &
                     NumberFluxDiff(iMlt, iMlat) + NumberFluxDiff(iMlt, iMlat2)

                ! If there are values in both hemisphere, then divide by 2
                if ( NumberFluxDiff(iMlt, iMlat) * &
                     NumberFluxDiff(iMlt, iMlat2) /= 0) &
                     numflux = numflux/2

             else
                numflux = NumberFluxDiff(iMlt, iMlat)
             endif

             if (numflux /= 0) then
                ElectronAverageEnergy(iLon,iLat) = &
                     ElectronEnergyFlux(iLon, iLat)/numflux * &
                     6.242e11 / 1000.0 ! ergs -> keV
             endif

             ! Mono Energy Flux

             if (UseNewellAveraged .or. EnergyFluxMono(iMlt, iMlat)==0) then

                ! Add North and South together
                ElectronEnergyFluxMono(iLon, iLat) = &
                     EnergyFluxMono(iMlt, iMlat) + &
                     EnergyFluxMono(iMlt, iMlat2)

                ! If there are values in both hemisphere, then divide by 2
                if ( EnergyFluxMono(iMlt, iMlat) * &
                     EnergyFluxMono(iMlt, iMlat2) /= 0) &
                     ElectronEnergyFluxMono(iLon, iLat) = &
                     ElectronEnergyFluxMono(iLon, iLat)/2.0
             else
                ElectronEnergyFluxMono(iLon, iLat) = &
                     EnergyFluxMono(iMlt, iMlat)
             endif

             ! Mono Number Flux

             if (UseNewellAveraged .or. NumberFluxMono(iMlt, iMlat)==0) then

                ! Add North and South together
                ElectronNumberFluxMono(iLon, iLat) = &
                     NumberFluxMono(iMlt, iMlat) + NumberFluxMono(iMlt, iMlat2)

                ! If there are values in both hemisphere, then divide by 2
                if ( NumberFluxMono(iMlt, iMlat) * &
                     NumberFluxMono(iMlt, iMlat2) /= 0) &
                     ElectronNumberFluxMono(iLon, iLat) = &
                     ElectronNumberFluxMono(iLon, iLat)/2.0

             else
                ElectronNumberFluxMono(iLon, iLat) = &
                     NumberFluxMono(iMlt, iMlat)
             endif

             ! Wave Energy Flux

             if (UseNewellAveraged .or. EnergyFluxWave(iMlt, iMlat)==0) then

                ! Add North and South together
                ElectronEnergyFluxWave(iLon, iLat) = &
                     EnergyFluxWave(iMlt, iMlat) + &
                     EnergyFluxWave(iMlt, iMlat2)

                ! If there are values in both hemisphere, then divide by 2
                if ( EnergyFluxWave(iMlt, iMlat) * &
                     EnergyFluxWave(iMlt, iMlat2) /= 0) &
                     ElectronEnergyFluxWave(iLon, iLat) = &
                     ElectronEnergyFluxWave(iLon, iLat)/2.0
             else
                ElectronEnergyFluxWave(iLon, iLat) = &
                     EnergyFluxWave(iMlt, iMlat)
             endif

             ! Wave Number Flux

             if (UseNewellAveraged .or. NumberFluxWave(iMlt, iMlat)==0) then

                ! Add North and South together
                ElectronNumberFluxWave(iLon, iLat) = &
                     NumberFluxWave(iMlt, iMlat) + NumberFluxWave(iMlt, iMlat2)

                ! If there are values in both hemisphere, then divide by 2
                if ( NumberFluxWave(iMlt, iMlat) * &
                     NumberFluxWave(iMlt, iMlat2) /= 0) &
                     ElectronNumberFluxWave(iLon, iLat) = &
                     ElectronNumberFluxWave(iLon, iLat)/2.0

             else
                ElectronNumberFluxWave(iLon, iLat) = &
                     NumberFluxWave(iMlt, iMlat)
             endif

          endif
       enddo
    enddo

    call end_timing("run_newell")

  end subroutine run_newell

  ! -------------------------------------------------------------------

  subroutine calc_flux(ProbTotal, b1a, b2a, Flux)

    real, dimension(nMlts, nMlats), intent(in)  :: b1a, b2a
    real, dimension(nMlts, nMlats), intent(in)  :: ProbTotal
    real, dimension(nMlts, nMlats), intent(out) :: Flux

    where (ProbTotal > 0) &
         Flux = (dFdt * b2a + b1a) * ProbTotal
    where (Flux < 0) Flux = 0.0

  end subroutine calc_flux


  ! -------------------------------------------------------------------

  subroutine calc_probability(Prob, b1p, b2p, ProbTotal)

    real, dimension(nMlts, nMlats), intent(in)      :: b1p, b2p
    real, dimension(ndF, nMlts, nMlats), intent(in) :: Prob
    real, dimension(nMlts, nMlats), intent(out)     :: ProbTotal

    integer :: idfp, idfm
    integer :: iMlt, iMLat

    ProbTotal = b1p + b2p*dFdt
    where (ProbTotal < 0.0) ProbTotal = 0.0
    where (ProbTotal > 1.0) ProbTotal = 1.0

    do iMlt = 1, nMlts
       do iMlat = 1, nMlats
          if (b1p(iMlt, iMlat) == 0 .and. b2p(iMlt, iMlat) == 0) then
             ProbTotal(iMlt, iMlat) = Prob(dfBin, iMlt, iMlat)
             if (ProbTotal(iMlt, iMlat) == 0.0) then
                idfp = dfBin + 1
                idfm = dfBin - 1
                if (idfm < 1)   idfm = dfBin + 2
                if (idfp > ndF) idfp = dfBin - 2
                ProbTotal(iMlt, iMlat) = &
                     (Prob(idfm, iMlt, iMlat)+Prob(idfp, iMlt, iMlat))/2
             endif
          endif
       enddo
    enddo

  end subroutine calc_probability

  ! -------------------------------------------------------------------

  subroutine calc_dfdt(by, bz, vx)

    use ModConstants, only : pi

    real, intent(in) :: by, bz, vx
    real :: dFAve, dFStep
    real :: v, sintc, bt, tc, bzt

    bzt = bz
    bt = sqrt(by**2 + bz**2)
    v = abs(vx)

    if (bzt == 0.0) bzt = 0.001
    tc = atan2(by,bzt)
    if (bt*cos(tc)*bz < 0.0) tc = tc + pi

    sintc = abs(sin(tc/2.))

    dFdt = (v**1.33333)*(sintc**2.66667)*(BT**0.66667)

    dfave = 4421.0
    dfStep = dfAve/(ndF-1)
    dfbin = min(max(floor(dFdt/dfStep),1),ndF)

  end subroutine calc_dfdt

  ! -------------------------------------------------------------------

  subroutine calc_hp(value,outs,outn)

    real, dimension(nMlts, nMlats), intent(in) :: value
    real, intent(out) :: outs, outn

    outs = sum(value(:,1:nMLats/2-1) * Area(:,1:nMLats/2-1))
    outn = sum(value(:,nMLats/2:nMLats) * Area(:,nMLats/2:nMLats))

  end subroutine calc_hp

  ! -------------------------------------------------------------------

  subroutine smooth(value)

    use ModSizeGITM
    use ModInputs

    implicit none

    real, dimension(nMlts, nMlats) :: value, valueout

    integer :: nPL = 2, nPM = 2, nMin = 2
    integer :: iMlt, iLat, iM, iL, n, iMa, iLa, iHem
    integer :: iLatStart, iLatEnd
    real    :: ave, std

    ! for removing bad points, we want the zone of consideration to be at 
    ! least 2 cells on each side (i.e., 5x5).  For averaging, allow only 1 cell.
    if (DoNewellAverage) nMin = 1

    ! How many points to average over in Lat and MLT.
    ! Newell is at 1/2 deg resolution in lat, so averaging would be over 
    nPL = max(1,floor(180.0 / float(nBlocksLat*nLats) + 0.499))
    ! no *2, because this is for each side

    ! Newell is at 1/4 hour MLT, which is 3.75 deg so averaging would be over
    nPM = max(1,floor(360.0 / float(nBlocksLon*nLons)/3.75/2 + 0.499))

!    write(*,*) "nPL, nPM in smooth : ", nPL, nPM

    valueout = value*0.0
    do iMlt = 1, nMlts
       do iHem = 1,2
          if (iHem == 1) then
             iLatStart = nPl+1
             iLatEnd = nMlats/2-nPl-1
          endif
          if (iHem == 2) then
             iLatStart = nMlats/2+nPl+1
             iLatEnd = nMlats-nPl-1
          endif

          do iLat = iLatStart, iLatEnd
             if (value(iMlt, iLat) > 0.0) then
                n = 0
                ave = 0.0
                do iM = iMlt-nPM, iMlt+nPM
                   iMa = iM
                   iMa = mod(iMa + nMlts, nMlts)
                   if (iMa == 0) iMa = nMlts
                   do iLa = iLat-nPL, iLat+nPL
                      if (value(iMa,iLa) > 0.0) then
                         ave = ave + value(iMa,iLa)
                         n   = n + 1
                      endif
                   enddo
                enddo
                if (n > (2*nPL+1)*(2*nPM+1)/2) then 
                   ave = ave/n
                   std = 0.0
                   do iM = iMlt-nPM, iMlt+nPM
                      iMa = iM
                      iMa = mod(iMa + nMlts, nMlts)
                      if (iMa == 0) iMa = nMlts
                      do iLa = iLat-nPL, iLat+nPL
                         if (value(iMa,iLa) > 0.0) then
                            std = std + (ave - value(iMa,iLa))**2
                         endif
                      enddo
                   enddo
                   std = sqrt(std/n)
                ! We only want to kill points that are 2 stdev ABOVE the average
                ! value.
                   if (abs(value(iMlt,iLat)-ave) > 2*std .or. &
                        DoNewellAverage) then
                   !write(*,*) "ave : ", valueout(iMlt,iLat), ave, &
                   !   std, abs(ave - value(iMlt,iLat)), 2*std
                      valueout(iMlt,iLat) = ave
                   else
                      valueout(iMlt,iLat) = value(iMlt,iLat)
                   endif
                endif
             endif
          enddo
       enddo
    enddo

    value = valueout

  end subroutine smooth

  ! -------------------------------------------------------------------

  subroutine read_single_regression_file(cFile, rFa, b1a, b2a)

    use ModInputs, only : iInputUnit_

    character (len=*), intent(in) :: cFile
    real, dimension(nMlts, nMlats), intent(out) :: rFa, b1a, b2a
    integer :: year0, day0, year1, day1, nFiles, sf0
    integer :: iMlt, iMlat, i, j, iError

    iError = 0

    open(iInputUnit_,file=cFile,status="old",iostat=iError)

    if (iError /= 0) then
       write(*,*) "Error in read_single_regression_file"
       call stop_gitm(cFile//" cannot be opened")
    endif

    read(iInputUnit_, *, iostat=iError) year0, day0, year1, day1, nFiles, sf0

    if (iError /= 0) then
       write(*,*) "Error in read_single_regression_file"
       call stop_gitm(cFile//" cannot read first line")
    endif

    do iMlt = 1, nMlts
       do iMlat = 1, nMlats
          read(iInputUnit_,*,iostat=iError) &
               i,j,b1a(iMlt,iMlat),b2a(iMlt,iMlat),rfa(iMlt,iMlat)
          if (iError /= 0) then
             write(*,*) "Error in read_single_regression_file:", iMlt, iMlat
             call stop_gitm(cFile//" error reading file")
          endif
       enddo
    enddo

    close(iInputUnit_)

  end subroutine read_single_regression_file

  ! -------------------------------------------------------------------

  subroutine read_single_probability_file(cFile, b1p, b2p, Prob)

    use ModInputs, only : iInputUnit_

    character (len=*), intent(in) :: cFile
    real, dimension(nMlts, nMlats), intent(out) :: b1p, b2p
    real, dimension(ndF, nMlts, nMlats), intent(out) :: Prob
    integer :: year0, day0, year1, day1, nFiles, sf0
    integer :: iMlt, iMlat, idF, iError

    iError = 0

    open(iInputUnit_,file=cFile,status="old",iostat=iError)

    if (iError /= 0) then
       write(*,*) "Error in read_single_probability_file"
       call stop_gitm(cFile//" cannot be opened")
    endif

    read(iInputUnit_, *, iostat=iError) year0, day0, year1, day1, nFiles, sf0

    if (iError /= 0) then
       write(*,*) "Error in read_single_probability_file"
       call stop_gitm(cFile//" cannot read first line")
    endif

    do iMlt = 1, nMlts
       do iMlat = 1, nMlats
          read(iInputUnit_,*,iostat=iError) &
               b1p(iMlt,iMlat),b2p(iMlt,iMlat)
          if (iError /= 0) then
             write(*,*) "Error in read_single_probability_file:", iMlt, iMlat
             call stop_gitm(cFile//" error reading file")
          endif
       enddo
    enddo

    do iMlt = 1, nMlts
       do iMlat = 1, nMlats
          do idF = 1, ndF
             read(iInputUnit_,*,iostat=iError) Prob(idF,iMlt,iMlat)
             if (iError /= 0) then
                write(*,*) "Error in read_single_probability_file:", idF,iMlt, iMlat
                call stop_gitm(cFile//" error reading file")
             endif
          enddo
       enddo
    enddo

    close(iInputUnit_)

  end subroutine read_single_probability_file

  ! -------------------------------------------------------------------

  subroutine read_all_regression_files

    character (len=iCharLen_) :: cFile

    cFile = cFileDiffef
    call merge_str(dir,cFile)
    call read_single_regression_file(cFile, rFaDiff, b1aDiff, b2aDiff)
    cFile = cFileDiffnf
    call merge_str(dir,cFile)
    call read_single_regression_file(cFile, rFaDiffn, b1aDiffn, b2aDiffn)

    cFile = cFileMonoef
    call merge_str(dir,cFile)
    call read_single_regression_file(cFile, rFaMono, b1aMono, b2aMono)
    cFile = cFileMononf
    call merge_str(dir,cFile)
    call read_single_regression_file(cFile, rFaMonon, b1aMonon, b2aMonon)

    cFile = cFileWaveef
    call merge_str(dir,cFile)
    call read_single_regression_file(cFile, rFaWave, b1aWave, b2aWave)
    cFile = cFileWavenf
    call merge_str(dir,cFile)
    call read_single_regression_file(cFile, rFaWaven, b1aWaven, b2aWaven)

    cFile = cFileIonsef
    call merge_str(dir,cFile)
    call read_single_regression_file(cFile, rFaIons, b1aIons, b2aIons)
    cFile = cFileIonsnf
    call merge_str(dir,cFile)
    call read_single_regression_file(cFile, rFaIonsn, b1aIonsn, b2aIonsn)

  end subroutine read_all_regression_files

  ! -------------------------------------------------------------------

  subroutine read_all_probability_files

    character (len=iCharLen_) :: cFile

    cFile = cFileDiffp
    call merge_str(dir,cFile)
    call read_single_probability_file(cFile, b1pDiff, b2pDiff, ProbDiff)

    cFile = cFileMonop
    call merge_str(dir,cFile)
    call read_single_probability_file(cFile, b1pMono, b2pMono, ProbMono)

    cFile = cFileWavep
    call merge_str(dir,cFile)
    call read_single_probability_file(cFile, b1pWave, b2pWave, ProbWave)

  end subroutine read_all_probability_files

end module ModNewell
