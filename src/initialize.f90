!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine initialize_gitm(TimeIn)

  use ModGITM
  use ModInputs
  use ModConstants
  use ModPlanet
  use ModRates
  use ModSphereInterface
  use ModTime
  use ModEUV
  use ModIndicesInterfaces
  implicit none

  type (UAM_ITER) :: r_iter

  real(Real8_), intent(in) :: TimeIn
  integer :: iLat, iAlt, iBlock, iSpecies, iLon, iError,nAltsMean,iilon,iilat

  real :: TempAve
  real :: TempDiff
  real :: InvScaleHeightS(-1:nLons+2,-1:nLats+2)

  real :: LogRho(-1:nLons+2,-1:nLats+2), NewSumRho(-1:nLons+2,-1:nLats+2)
  real :: GradAlt_CD(nLons, nLats, nAlts, 3)

  real :: TempUnit_const, t, h,meanaltzero,altdiff, dAlt

  logical :: IsThere, IsOk, IsDone, IsFirstTime = .true.

  real :: DistM, DistP, Ratio2, InvDenom
  !----------------------------------------------------------------------------

!! Sorry but this is a double-negative
!! This checks to see if the planet is not, not Titan.
!! That is to say, it checks to see if this is actually Titan.

  if (   .not. (index(cPlanet,"Titan") == 0)  ) then 
     call init_radcooling
  endif

  if (.not.IsFirstTime) return

  IsFirstTime = .false.

  call start_timing("initialize")
  call report("initialize",1)

  ! ------------------------------------------------------------------------
  ! initialize stuff
  ! ------------------------------------------------------------------------

  if (TimeIn /= CurrentTime) then
     CurrentTime = TimeIn
     call time_real_to_int(CurrentTime, iTimeArray)
     call fix_vernal_time
  endif

  call get_f107(CurrentTime, f107, iError)
  call get_f107a(CurrentTime, f107a, iError)

  call init_grid

  if (DoRestart) then
     call read_inputs("UA/restartIN/header.rst")
     call set_inputs
     call read_restart("UA/restartIN")
     call init_msis
!     if (UsePerturbation) call user_create_perturbation
  endif

  call set_RrTempInd

  inquire(file='GITM.STOP',EXIST=IsThere)
  if (IsThere .and. iProc == 0) then
     open(iOutputUnit_, file = 'GITM.STOP', status = 'OLD')
     close(iOutputUnit_, status = 'DELETE')
  endif

  !\
  ! Initialize the EUV wave spectrum
  !/

  call init_euv

  Gravity_GB = 0.0

  if (.not. DoRestart) then

     if (UseStretchedAltitude) then
        call init_altitude
     else
        if (UseTopography) then
           if (AltMin > 1.0) then 
              write(*,*) 'When using topography, the minimum altitude'
              write(*,*) 'must be zero.  Stopping...'
              call stop_gitm('Incorrect minimum altitude')
           endif
           altzero = 0.0
           call init_topography

           meanAltZero = sum(altzero)/(nlons*nlats*nblocks)
           altdiff = AltMinUniform - meanAltZero
           naltsmean = floor(nalts - nalts*((AltMax-AltMinUniform)/(AltMax)))

           do iBlock = 1, nBlocks
              do iLon = -1,nLons +2
                 do iLat = -1, nLats +2
                     iiLon = min(max(iLon,1),nLons)
                     iilat = min(max(iLat,1),nLats)

                    do iAlt = -1,nAltsmean
                       
                       dAlt = (AltMinUniform-AltZero(iiLon,iiLat,iBlock))/nAltsmean
                       Altitude_GB(iLon,iLat,iAlt,iBlock) = &
                            AltZero(iiLon,iiLat,iBlock) + ialt*dAlt
                    enddo
!                  if (ilon ==1) then
!                     write(93,*) ilat,altzero(iilon,iilat,iblock)
!                     write(93,*)ilat, iblock,latitude(ilat,iblock)*180/pi,altitude_GB(1,ilat,0,iblock)
!                  endif
                                                         
                 enddo
              enddo

              dAlt = (AltMax-AltMinUniform)/(nAlts-naltsmean)

              do iAlt = 1, (nAlts+2)-naltsmean
                 Altitude_GB(:,:,iAlt+naltsmean,1:nBlocks) = &
                      AltMinUniform + iAlt * dAlt
              enddo
           enddo

        else
           ! Uniform grid
           do iAlt=-1,nAlts+2
              Altitude_GB(:,:,iAlt,1:nBlocks) = &
                   AltMin + (iAlt-0.5)*(AltMax-AltMin)/nAlts
           enddo
        endif
        
     end if
  endif

  ! Calculate vertical cell sizes
  do iAlt = 0,nAlts+1
     ! Cell interface is taken to be half way between cell centers
     ! so the cell size is half of the cell center distance
     ! between cells i-1 and i+1: 
     dAlt_GB(:,:,iAlt,1:nBlocks) = 0.5* &
          ( Altitude_GB(:,:,iAlt+1,1:nBlocks) &
          - Altitude_GB(:,:,iAlt-1,1:nBlocks))
  enddo
  dAlt_GB(:,:,-1,1:nBlocks)      = dAlt_GB(:,:,0,1:nBlocks)
  dAlt_GB(:,:,nAlts+2,1:nBlocks) = dAlt_GB(:,:,nAlts+1,1:nBlocks)

  RadialDistance_GB(:,:,:,1:nBlocks)    = rBody + Altitude_GB(:,:,:,1:nBlocks)
  InvRadialDistance_GB(:,:,:,1:nBlocks) = 1/RadialDistance_GB(:,:,:,1:nBlocks)

  Gravity_GB(:,:,:,1:nBlocks) = -Gravitational_Constant &
       *(rBody*InvRadialDistance_GB(:,:,:,1:nBlocks)) ** 2

  if (UseStretchedAltitude) then
     do iAlt=1,nAlts
        Gravity_GB(:,:,iAlt,:) = -Gravitational_Constant &
             *(rBody/RadialDistance_GB(:,:,iAlt,:))**2 
     enddo
  endif

  if (iDebugLevel > 2) then
     do iAlt=-1,nAlts+2
        write(*,*) "===>Altitude : ", &
             iAlt, Altitude_GB(1,1,iAlt,1), RadialDistance_GB(1,1,iAlt,1), &
             Gravity_GB(1,1,iAlt,1)
     end do
  endif

  if (Is1D) then
     Latitude(0,1)  = Latitude(1,1) - 1.0 * pi/180.0
     Latitude(-1,1) = Latitude(0,1) - 1.0 * pi/180.0
     Latitude(2,1)  = Latitude(1,1) + 1.0 * pi/180.0
     Latitude(3,1)  = Latitude(2,1) + 1.0 * pi/180.0

     Longitude(0,1)  = Longitude(1,1) - 1.0 * pi/180.0
     Longitude(-1,1) = Longitude(0,1) - 1.0 * pi/180.0
     Longitude(2,1)  = Longitude(1,1) + 1.0 * pi/180.0
     Longitude(3,1)  = Longitude(2,1) + 1.0 * pi/180.0
  endif

  ! Precalculate the limited tangent and unlimited cosine of the latitude
  TanLatitude(:,1:nBlocks) = min(abs(tan(Latitude(:,1:nBlocks))),100.0) * &
       sign(1.0,Latitude(:,1:nBlocks))

  CosLatitude(:,1:nBlocks) = cos(Latitude(:,1:nBlocks))

  ! This is done so we don't get a /0 below.

  dLatDist_GB = 1.0
  dLatDist_FB = 1.0
  dLonDist_GB = 1.0
  dLonDist_FB = 1.0

  ! Precalculate the physical size of cells in the Lat and Lon directions
  do iBlock = 1, nBlocks
     do iAlt = -1, nAlts+2
        do iLat = 0, nLats+1
           do iLon = 0, nLons+1
              ! This is the cell size assuming that cell interface is half way
              dLatDist_GB(iLon, iLat, iAlt, iBlock) = 0.5 * &
                (Latitude(iLat+1, iBlock) - Latitude(iLat-1, iBlock)) * &
                RadialDistance_GB(iLon, iLat, iAlt, iBlock)

              ! This is the distance between neighboring cells
              ! Note that face(i) is between cells i and i-1 (like in BATSRUS)
              dLatDist_FB(iLon, iLat, iAlt, iBlock) = &
                   (Latitude(iLat, iBlock) - Latitude(iLat-1, iBlock)) * &
                   0.5*(RadialDistance_GB(iLon, iLat  , iAlt, iBlock) &
                   +    RadialDistance_GB(iLon, iLat-1, iAlt, iBlock))
              
              ! This is the cell size assuming that cell interface is half way
              dLonDist_GB(iLon, iLat, iAlt, iBlock) = 0.5 * &
                   (Longitude(iLon+1,iBlock) - Longitude(iLon-1,iBlock)) * &
                   RadialDistance_GB(iLon, iLat, iAlt, iBlock)* &
                   max(abs(CosLatitude(iLat,iBlock)),0.01)

              ! This is the distance between neighboring cells
              dLonDist_FB(iLon, iLat, iAlt, iBlock) = &
                   (Longitude(iLon,iBlock) - Longitude(iLon-1,iBlock)) * &
                   0.5*(RadialDistance_GB(iLon,   iLat, iAlt, iBlock) &
                   +    RadialDistance_GB(iLon-1, iLat, iAlt, iBlock)) &
                   *max(abs(CosLatitude(iLat,iBlock)),0.01)

              CellVolume(iLon,iLat,iAlt,iBlock) = &
                   dLonDist_GB(iLon, iLat, iAlt, iBlock) * &
                   dLatDist_GB(iLon, iLat, iAlt, iBlock) * &
                   dAlt_GB(iLon,iLat,iAlt,iBlock) 

           enddo

           ! Fill in longitude ghost cells
           dLonDist_FB(-1, iLat, iAlt, iBlock) = &
                dLonDist_FB(0, iLat, iAlt, iBlock)
           dLonDist_FB(nLons+2, iLat, iAlt, iBlock) = &
                dLonDist_FB(nLons+1, iLat, iAlt, iBlock)
           dLatDist_FB(-1, iLat, iAlt, iBlock) = &
                dLatDist_FB(0, iLat, iAlt, iBlock)
           dLatDist_FB(nLons+2, iLat, iAlt, iBlock) = &
                dLatDist_FB(nLons+1, iLat, iAlt, iBlock)
        enddo
        ! Fill in latitude ghost cells
        dLonDist_FB(:, -1, iAlt, iBlock) = &
             dLonDist_FB(:, 0, iAlt, iBlock)
        dLonDist_FB(:, nLats+2, iAlt, iBlock) = &
             dLonDist_FB(:, nLats+1, iAlt, iBlock)
        dLatDist_FB(:, -1, iAlt, iBlock) = &
             dLatDist_FB(:, 0, iAlt, iBlock)
        dLatDist_FB(:, nLats+2, iAlt, iBlock) = &
             dLatDist_FB(:, nLats+1, iAlt, iBlock)
     enddo
  enddo

  InvDLatDist_GB = 1.0/dLatDist_GB
  InvDLatDist_FB = 1.0/dLatDist_FB
  InvDLonDist_GB = 1.0/dLonDist_GB
  InvDLonDist_FB = 1.0/dLonDist_FB

  ! Precalculate the coefficients for the gradient calculation
  do iBlock = 1, nBlocks
     do iLon = 1, nLons
        DistM    = Longitude(iLon,iBlock) - Longitude(iLon-1,iBlock)
        DistP    = Longitude(iLon+1,iBlock) - Longitude(iLon,iBlock)
        Ratio2   = (DistM / DistP)**2
        InvDenom = 1.0/(Ratio2*DistP + DistM)

        GradLonP_CB(iLon, iBlock) =  InvDenom*Ratio2
        GradLon0_CB(iLon, iBlock) =  InvDenom*(1-Ratio2)
        GradLonM_CB(iLon, iBlock) = -InvDenom
     enddo

     do iLat = 1, nLats
        DistM    = Latitude(iLat,iBlock) - Latitude(iLat-1,iBlock)
        DistP    = Latitude(iLat+1,iBlock) - Latitude(iLat,iBlock)
        Ratio2   = (DistM / DistP)**2
        InvDenom = 1.0/(Ratio2*DistP + DistM)

        GradLatP_CB(iLat, iBlock) =  InvDenom*Ratio2
        GradLat0_CB(iLat, iBlock) =  InvDenom*(1-Ratio2)
        GradLatM_CB(iLat, iBlock) = -InvDenom
     enddo
  end do

  if(UseTopography)then
     !!! What about the maxi tricks ??? Why is that there???
     do iBlock = 1, nBlocks
        call UAM_gradient(Altitude_GB, GradAlt_CD, iBlock)
        dAltDLon_CB(:,:,:,iBlock) = GradAlt_CD(:,:,:,iEast_)
        dAltDLat_CB(:,:,:,iBlock) = GradAlt_CD(:,:,:,iNorth_)
     end do
  else
     dAltDLon_CB = 0.
     dAltDLat_CB = 0.
  end if

  call init_heating_efficiency

!! Some Titan-Specific Startup Routines here (Regardless of Restart or Not)

  if (   .not. (index(cPlanet,"Titan") == 0)  ) then 
     call init_magheat
     call init_isochem
     call init_aerosol
  endif

  if (   .not. (index(cPlanet,"Mars") == 0)  ) then 
     call init_isochem
  endif

  if (.not. DoRestart) then

     Potential = 0.0
     Velocity = 0.0
     IVelocity = 0.0
     VerticalVelocity = 0.0

     if (UseMsis) then
        call init_msis
     else

        TempUnit_const = 1. * Mass(1) / Boltzmanns_Constant
        TempAve  = (TempMax+TempMin)/2/TempUnit_const
        TempDiff = (TempMax-TempMin)/2/TempUnit_const

        do iBlock = 1, nBlocks

           do iAlt=-1,nAlts+2
              call get_temperature(0.0, 0.0, Altitude_GB(:,:,iAlt,iBlock), t, h)
              Temperature(:,:,iAlt,iBlock)  = t/TempUnit_const
              eTemperature(:,:,iAlt,iBlock) = t
              iTemperature(:,:,iAlt,iBlock) = t
           enddo
           
           do iAlt=-1,nAlts+2

              InvScaleHeight(:,:,iAlt,iBlock)  =  &
                   -Gravity_GB(:,:,iAlt,iBlock) / &
                   Temperature(:,:,iAlt,iBlock)

              Rho(:,:,iAlt,iBlock) = 0.0

              NewSumRho = 0.0

              do iSpecies = 1, nSpecies

                 InvScaleHeightS = -Gravity_GB(:,:,iAlt,iBlock) * &
                      Mass(iSpecies) / &
                      (Temperature(:,:,iAlt,iBlock)*TempUnit_const* &
                      Boltzmanns_Constant)

                 if(iAlt < 2)then
                    LogNS(:,:,iAlt,iSpecies,iBlock) = &
                         - (Altitude_GB(:,:,iAlt,iBlock)-AltMin) &
                         * InvScaleHeightS + LogNS0(iSpecies)
                 elseif(iAlt > 1 .and. iAlt < nAlts+2)then
                    LogNS(:,:,iAlt,iSpecies,iBlock) = &
                         LogNS(:,:,iAlt-1,iSpecies,iBlock) &
                         -(Altitude_GB(:,:,iAlt  ,iBlock) &
                         - Altitude_GB(:,:,iAlt-1,iBlock))*InvScaleHeightS &
                         -(Temperature(:,:,iAlt+1,iBlock)   &
                         - Temperature(:,:,iAlt-1,iBlock))  &
                         /(2.0*Temperature(:,:,iAlt,iBlock))
                 else
                    LogNS(:,:,iAlt,iSpecies,iBlock) = &
                         LogNS(:,:,iAlt-1,iSpecies,iBlock) &
                         -(Altitude_GB(:,:,iAlt  ,iBlock) &
                         - Altitude_GB(:,:,iAlt-1,iBlock))*InvScaleHeightS &
                         -(Temperature(:,:,iAlt  ,iBlock)  &
                         - Temperature(:,:,iAlt-1,iBlock)) &
                         /(Temperature(:,:,iAlt,iBlock))
                 endif
              
                 NewSumRho      = NewSumRho + &
                      Mass(iSpecies)*exp(LogNS(:,:,iAlt,iSpecies,iBlock))
              
              enddo
              
              do iSpecies=1,nSpecies

                 NDensityS(:,:,iAlt,iSpecies,iBlock) = &
                      exp(LogNS(:,:,iAlt,iSpecies,iBlock))

                 NDensity(:,:,iAlt,iBlock) = NDensity(:,:,iAlt,iBlock) + &
                      NDensityS(:,:,iAlt,iSpecies,iBlock)

                 Rho(:,:,iAlt,iBlock) = Rho(:,:,iAlt,iBlock) + &
                      Mass(iSpecies) * NDensityS(:,:,iAlt,iSpecies,iBlock)

              enddo

           enddo
           
        enddo

     endif

     if (UseIRI .and. IsEarth) then
        call init_iri
     else
        if (IsEarth) then
           do iBlock = 1, nBlocks
              IDensityS(:,:,:,:,iBlock)    = 1.00e8
              IDensityS(:,:,:,ie_,iBlock)  = 1.00e8*(nIons-1)
           enddo
        endif
     endif

  else

     do iBlock = 1, nBlocks
        InvScaleHeight(:,:,:,iBlock)  =  &
             -Gravity_GB(:,:,:,iBlock) / &
             Temperature(:,:,:,iBlock)
     enddo

  endif

  if (UseWACCMTides) then
     call read_waccm_tides
     call update_waccm_tides
  endif

  if (UseGSWMTides) then
     call read_tides
     call update_tides
  endif

  call init_b0
  if (IsEarth) call init_energy_deposition

  if (UseApex .and. IsEarth) then
     call report("subsolr",2)
     call SUBSOLR(iTimeArray(1),iJulianDay,iTimeArray(4),&
          iTimeArray(5),iTimeArray(6),SubsolarLatitude, &
          SubsolarLongitude)
  endif

  if (.not.Is1D) call exchange_messages_sphere

  call calc_pressure

  ! The iLon and iLat are dummy variables...
  call UA_calc_electrodynamics(iLon,iLat)

  do iBlock = 1, nBlocks
    call calc_eddy_diffusion_coefficient(iBlock)
    call calc_rates(iBlock)
    call calc_viscosity(iBlock)
  enddo

  call end_timing("initialize")

end subroutine initialize_gitm

