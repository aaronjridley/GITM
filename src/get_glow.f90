!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine get_glow(iLon,iLat,iBlock)

  use ModPlanet, only : nEmissionWavelengths,nPhotoBins
  use ModTime
  use ModGITM
  use ModEUV, only    : WAVES,WAVEL,Num_WaveLengths_High,Flux_of_EUV

  implicit none

  integer, intent(in)    :: iLon, iLat, iBlock
  real, dimension(nAlts) :: ZRHO,ZO,ZN2,ZO2,ZNO,ZNS,ZND
  real, dimension(nAlts) :: ZTN,ZTI,ZTE,ZE,Zo4sp,Zo2p,Znop,Zn2p,Znp
  real, dimension(nEmissionWavelengths,nAlts) :: emissionrates
  real, dimension(nPhotoBins,nAlts) :: peup, pedown, pespec
  integer :: iAlt,iError,iBin
 
  integer, parameter :: iO3P  = 1
  integer, parameter :: iO2 = 2
  integer, parameter :: iN2 = 3
  integer, parameter :: iN4S =  4
  integer, parameter :: iNO   =  5
  integer, parameter :: iN2D =  6
  integer, parameter  :: iO4SP = 1
  integer, parameter  :: iO2P   = 2
  integer, parameter  :: iN2P   = 3
  integer, parameter  :: iNP    = 4
  integer, parameter  :: iNOP   = 5
  integer, parameter  :: iO2DP = 6
  integer, parameter  :: ie    = 10

  call GL_settime(iTimeArray,utime)
  call GL_setloc(Altitude_GB(iLon,iLat,1:nAlts,iBlock),nAlts,Latitude(iLat,iBlock),&
       Longitude(iLon,iBlock),nPhotoBins,iError)
  if (iError.gt.0) then
     write(*,*) "nAlts too large for glow!.. in get_glow"
     call stop_GITM
  end if

  if (isFirstGlow) then
     call GL_init
     isFirstGlow = .False.
  endif

  if (isInitialGlow) then
     call GL_interp_flux(WaveS,WaveL,Flux_of_EUV,Num_WaveLengths_High)
     isInitialGlow = .False.
  endif

  ZRHO(1:nAlts)  = Rho(iLon,iLat,1:nAlts,iBlock)*0.000001*1000.0
  ZO(1:nAlts)    = NDensityS(iLon,iLat,1:nAlts,iO3P,iBlock)*0.000001
  ZN2(1:nAlts)   = NDensityS(iLon,iLat,1:nAlts,iN2,iBlock)*0.000001
  ZO2(1:nAlts)   = NDensityS(iLon,iLat,1:nAlts,iO2,iBlock)*0.000001
  ZNO(1:nAlts)   = NDensityS(iLon,iLat,1:nAlts,iNO,iBlock)*0.000001
  ZNS(1:nAlts)   = NDensityS(iLon,iLat,1:nAlts,iN4S,iBlock)*0.000001
  ZND(1:nAlts)   = NDensityS(iLon,iLat,1:nAlts,iN2D,iBlock)*0.000001
 
  call GL_setND(ZRHO,ZO,ZN2,ZO2,ZNO,ZNS,ZND)

  ZTN(1:nAlts)   = Temperature(iLon,iLat,1:nAlts,iBlock) * &
             TempUnit(iLon,iLat,1:nAlts)
  ZTI(1:nAlts)   = iTemperature(iLon,iLat,1:nAlts,iBlock)
  ZTE(1:nAlts)   = eTemperature(iLon,iLat,1:nAlts,iBlock)
  
  call GL_settemp(ZTN,ZTI,ZTE)

  Zo4sp = IDensityS(iLon,iLat,1:nAlts,iO4SP,iBlock)*0.000001
  Zo2p =  IDensityS(iLon,iLat,1:nAlts,iO2P,iBlock)*0.000001
  Znop = IDensityS(iLon,iLat,1:nAlts,iNOP,iBlock)*0.000001
  ZE = IDensityS(iLon,iLat,1:nAlts,ie,iBlock)*0.000001
  Zn2p = IDensityS(iLon,iLat,1:nAlts,iN2P,iBlock)*0.000001
  Znp =  IDensityS(iLon,iLat,1:nAlts,iNP,iBlock)*0.000001

  call GL_setID(Zo4sp, Zo2p, Znop, Zn2p, Znp, ZE)

  call GL_getvals(iProc,emissionrates,peup,pedown,pespec,PhotoEBins)

  do iAlt = 1, nAlts
     vEmissionRate(iLon,iLat,iAlt,:,iBlock)     = emissionrates(:,iAlt) * 100000.0
     PhotoEFluxU(iLon,iLat,iAlt,:,iBlock)       = peup(:,iAlt)   * 1000000.0
     PhotoEFluxD(iLon,iLat,iAlt,:,iBlock)       = pedown(:,iAlt) * 1000000.0
     PhotoElectronRate(iLon,iLat,iAlt,:,iBlock) = pespec(:,iAlt) * 1000000.0
  enddo


     PhotoEFluxTotal = 0

     !Sum over all Energy Bins...
     do iBin = 1, nPhotoBins
        PhotoEFluxTotal(:,:,:,iBlock,1) = PhotoEFluxTotal(:,:,:,iBlock,1) + &
             PhotoEFluxU(:,:,:,ibin,iBlock)
        PhotoEFluxTotal(:,:,:,iBlock,2) = PhotoEFluxTotal(:,:,:,iBlock,2) + &
             PhotoEFluxD(:,:,:,iBin,iBlock)
     enddo


end subroutine get_glow
