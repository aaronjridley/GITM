!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine user_create_perturbation

  use ModGITM
  use ModInputs

  NDensityS(nLons/2,nLats/2,2,1:3,1) = NDensityS(nLons/2,nLats/2,2,1:3,1)*50.0
  temperature(nLons/2,nLats/2,2,1) = temperature(nLons/2,nLats/2,2,1)*100.0

end subroutine user_create_perturbation


subroutine user_perturbation

  use ModGITM
  use ModInputs
  use ModNumConst
  use ModKind, ONLY: Real8_
  use ModTime
  use ModSources, only: UserHeatingRate

  implicit none

  real         :: latcenter, loncenter, amp
  real         :: latwidth, lonwidth
  real         :: f(-1:nLons+2,-1:nLats+2)
  integer      :: iBlock, iSpecies, iLat, iLon, iAlt
  real(Real8_) :: PerturbTimeStart, PerturbTimeEnd, MidTime
  real         :: tsave = 0.0
  real         :: dla, dlo, lac, loc
  
  UserHeatingRate = 0.0

  ! Start 1 hour1 after start of the simulation
  PerturbTimeStart = StartTime + 1.0*3600.0+60.0
  ! End 1 minute after that
  PerturbTimeEnd   = PerturbTimeStart + 60.0

  dla = latitude(2,1) - latitude(1,1)
  lac = (LatEnd+LatStart)/2.0
  dlo = longitude(2,1) - longitude(1,1)
  loc = (LonEnd+LonStart)/2.0

  latcenter = lac + dla/2.0
  loncenter = loc + dlo/2.0
  latwidth  = dla/2.0
  lonwidth  = dlo/2.0

  if  (CurrentTime < PerturbTimeStart) then
     tsave = sum(temperature(1:nLons,1:nLats,1,1:nBlocks)) / &
          (nLons*nLats*nBlocks)
  endif

  DuringPerturb = .false.
  
  if  (CurrentTime >= PerturbTimeStart .and. &
       CurrentTime <  PerturbTimeEnd) then

     MidTime = (PerturbTimeStart + PerturbTimeEnd)/2.0

     DuringPerturb = .true.

     amp = exp(-((CurrentTime-MidTime)/(PerturbTimeEnd-PerturbTimeStart)*5)**2)
!     write(*,*) "Perturbing!", tsave, amp,latcenter/3.1415*180.0,&
!          loncenter/3.1415*180.0, latwidth, lonwidth

     do iBlock = 1, nBlocks
        do iLon = 1, nLons
           do iLat = 1, nLats
              if ((abs(latitude(iLat,iBlock)-latcenter) < 4*latwidth) .and. &
                  ( abs(longitude(iLon,iBlock)-loncenter) < 4*lonwidth)) then
                 f(iLon,iLat) = & !amp*&
                      exp(-((latitude(iLat,iBlock)-latcenter)/latwidth)**2) * &
                      exp(-((longitude(iLon,iBlock)-loncenter)/lonwidth)**2)
              else
                 f(iLon,iLat) = 0.0
              endif

              iAlt = 1
              UserHeatingRate(iLon,iLat,iAlt,iBlock) =  &
                   1000.0 * 4.184e6 * f(iLon,iLat) / &
                   cellvolume(iLon,iLat,iAlt,iBlock) / &
                   TempUnit(iLon,iLat,iAlt) / &
                   cp(iLon,iLat,iAlt,iBlock) / &
                   rho(iLon,iLat,iAlt,iBlock)
           enddo
        enddo

!        temperature(:,:,-1,iBlock) = tsave + 5.0 * f * tsave
!        temperature(:,:, 0,iBlock) = tsave + 5.0 * f * tsave
!        temperature(:,:, 1,iBlock) = tsave + 5.0 * f * tsave
!        do iSpecies = 1, nSpecies
!           VerticalVelocity(:, :, -1, iSpecies, iBlock) = 2*f
!           VerticalVelocity(:, :,  0, iSpecies, iBlock) = 2*f
!           VerticalVelocity(:, :,  1, iSpecies, iBlock) = 2*f
!        enddo
!        Velocity(:,:,-1, iUp_, iBlock) = 2*f
!        Velocity(:,:, 0, iUp_, iBlock) = 2*f
!        Velocity(:,:, 1, iUp_, iBlock) = 2*f
     enddo

  endif

!  NDensityS(nLons/2,nLats/2,2,1:3,1) = NDensityS(nLons/2,nLats/2,2,1:3,1)*50.0
!  temperature(nLons/2,nLats/2,2,1) = temperature(nLons/2,nLats/2,2,1)*100.0

end subroutine user_perturbation


!=================================================================
subroutine user_bc_perturbation(LogRhoBc, LogNSBc, VelBc_GD, TempBc)

  ! PURPOSE: 
  !   This subroutine contains the source code for Wave Perturbation (WP) 
  !   model for simulating tsunami or earthquake generated atmospheric
  !   acoustic-gravity waves. The waves perturb the neutral atmosphere, 
  !   propagate upward, and enter GITM through GITM's lower boundary. 
  !
  ! USAGE:
  !   This subroutine is called by specifying UseBcPerturbation = .true.
  !   and providing surface perturbation characteristics under
  !   #USEBCPERTURBATION in UAM.in. For details please refer to 
  !   srcDoc/manual_WPGITM.pdf. 
  !   
  ! NOTE:
  !   The time of perturbation is set to the time delay since StartTime
  !   ("CurrentTime-StartTime-PerturbTimeDelay" below). One needs to make 
  !   sure 'PerturbTimeDelay' is set properly so the perturbation is added 
  !   at the desired time, especially for restart runs.
  !   This can be resolved by implementing the start time (year month day
  !   hour minute second) of the perturbation in the future.
  !
  ! AUTHOR: 
  !   Xing Meng (Xing.Meng@jpl.nasa.gov or xingm@umich.edu)
  !
  ! REVISION HISTORY:
  !   Initial version: March 19, 2021
  !
  ! LICENSE:
  !   Copyright 2021 California Institute of Technology
  !   Licensed under the Apache License, Version 2.0 (the "License");
  !   you may not use this code except in compliance with the License.
  !   You may obtain a copy of the License at
  !   http://www.apache.org/licenses/LICENSE-2.0

  use ModGITM, ONLY: iEast_, iNorth_, iUp_, iProc, MeanTempBc, MeanVelBc_D
  use ModVertical, ONLY: Lat, Lon, Altitude_G, Gravity_G
  use ModPlanet, ONLY: nSpecies, RBody
  use ModConstants, ONLY: pi, Univ_Gas_Constant, Gamma_const
  use ModTime, only: CurrentTime, StartTime
  use ModInputs, ONLY: iTypeBcPerturb, RefLon, RefLat, EpicenterLon, &
       EpicenterLat, PerturbTimeDelay, PerturbDuration, SeisWaveTimeDelay, &
       PerturbWaveDirection, PerturbWaveSpeed, PerturbWavePeriod, &
       PerturbWaveHeight, EpiDistance, FFTReal, FFTImag, &
       PerturbWaveFreq, nPerturbFreq, nMaxPerturbFreq, idebugLevel

  implicit none

  real, intent(inout):: LogRhoBc(-1:0), LogNSBc(-1:0,nSpecies), &
       VelBc_GD(-1:0,iEast_:iUp_), TempBc(-1:0)
  ! Temp0 and Rho0 -- neutral temperature[K] and density[kg/m^3] at 0km 
  ! altitude, can be modified according to the specific region/time to model.  
  real, parameter :: Temp0 = 288.15, Rho0 = 1.225, dLat = 1.11e5, &
       AvgGravity = 9.65674
  real :: ScaleHeight = 0.0, BouyancyFreq2 = 0.0, SoundSpeed2 = 0.0, &
       OldTime = 0.0, Slope = 0.0, WaveSpeedLat = 0.0
  real :: dLon, RadialD, Kr, Kh, Kh2, Kx, Ky, Kz, WaveFreqIntri, WaveFreqIntri2
  real :: WaveForm, CosineWave, SineWave, WaveCombo1, WaveCombo2, &
       Coeff1, Coeff2, Coeff3, CoeffS, CoeffC, dRho, dVel
  integer :: iFreq, iAlt
  logical :: IsFirstCall=.true., DoPerturb
  logical, dimension(1:nMaxPerturbFreq) :: IsForbidden = .false.

! ----------------------------------------------------------------

  if (IsFirstCall .and. iTypeBcPerturb < 2) then
     ! For tsunami plane waves, these are constants and only need to 
     ! be calculated once
     Slope = tan(pi/2.+PerturbWaveDirection*pi/180.)
     WaveSpeedLat = PerturbWaveSpeed/sin(PerturbWaveDirection*pi/180.)
     IsFirstCall = .false.
  endif

  if(CurrentTime /= OldTime)then
     ! These are the same for all grid cells within each time step
     ScaleHeight = Univ_Gas_Constant/2.9e-2*(Temp0 + MeanTempBc)/2./AvgGravity
     BouyancyFreq2 = (Gamma_const-1)/Gamma_const*AvgGravity/ScaleHeight
     SoundSpeed2 = Gamma_const*AvgGravity*ScaleHeight
     OldTime = CurrentTime
  endif

  ! longitude degree-meter conversion factor of the current Lat
  dLon = pi*Rbody*cos(Lat*pi/180.)/180.

  do iFreq = 1, nPerturbFreq

     ! skip zero frequency and non-propagating frequency components
     if (PerturbWaveFreq(iFreq) == 0 .or. IsForbidden(iFreq)) cycle

     if (iTypeBcPerturb < 2) then 
        Kh = PerturbWaveFreq(iFreq)/PerturbWaveSpeed
        Kh2 = Kh**2
        Kx = Kh*cos(PerturbWaveDirection*pi/180.)
        Ky = Kh*sin(PerturbWaveDirection*pi/180.)
        WaveFreqIntri = PerturbWaveFreq(iFreq) - &
             Kx*MeanVelBc_D(iEast_)/2. - Ky*MeanVelBc_D(iNorth_)/2.
        WaveFreqIntri2 = WaveFreqIntri**2

        if(WaveFreqIntri < 0 .or. &
             (WaveFreqIntri2 >= SoundSpeed2*(Kh2+1/ScaleHeight**2/4.)/2. &
             - sqrt(SoundSpeed2**2*(Kh2+1/ScaleHeight**2/4.)**2 &
             - 4.*Kh2*SoundSpeed2*BouyancyFreq2)/2. &
             .and. WaveFreqIntri2 <= SoundSpeed2*(Kh2+1/ScaleHeight**2/4.)/2. &
             + sqrt(SoundSpeed2**2*(Kh2+1/ScaleHeight**2/4.)**2 &
             - 4.*Kh2*SoundSpeed2*BouyancyFreq2)/2.))then
           if (.not. IsForbidden(iFreq) .and. iProc == 0)then
              write(*,*) 'No longer a upward-propagating waves at Time = ', &
                   CurrentTime-StartTime, 'iFreq =', iFreq
              IsForbidden(iFreq) = .true.
           endif
           cycle
        endif

        Kz = sqrt(WaveFreqIntri2/SoundSpeed2 + &
             (BouyancyFreq2-WaveFreqIntri2)*Kh2/WaveFreqIntri2 - &
             1/(4.*ScaleHeight**2))

        ! find region to perturb
        if (WaveSpeedLat > 0) then
           DoPerturb = (Lat-RefLat)*dLat - Slope*(Lon-RefLon)*dLon <= &
                WaveSpeedLat*(CurrentTime-StartTime-PerturbTimeDelay) &
                .and. (Lat-RefLat)*dLat - Slope*(Lon-RefLon)*dLon >= &
                WaveSpeedLat*(CurrentTime-StartTime-PerturbTimeDelay &
                -PerturbDuration)
        else
           DoPerturb = (Lat-RefLat)*dLat - Slope*(Lon-RefLon)*dLon >= &
                WaveSpeedLat*(CurrentTime-StartTime-PerturbTimeDelay) &
                .and. (Lat-RefLat)*dLat - Slope*(Lon-RefLon)*dLon <= &
                WaveSpeedLat*(CurrentTime-StartTime-PerturbTimeDelay &
                -PerturbDuration)
        endif
     endif

     do iAlt = -1, 0

        if (iTypeBcPerturb == 2) then           
           RadialD = sqrt(((Lon-EpicenterLon)*dLon)**2 + &
                ((Lat-EpicenterLat)*dLat)**2 + Altitude_G(iAlt)**2)
           ! unlike plane waves, neutral wind effect is negelected here
           WaveFreqIntri = PerturbWaveFreq(iFreq)
           WaveFreqIntri2 = WaveFreqIntri**2
          
           if((WaveFreqIntri2-SoundSpeed2/(4.*ScaleHeight**2))* &
                (WaveFreqIntri2-BouyancyFreq2*(((Lon-EpicenterLon)*dLon)**2 + &
                ((Lat-EpicenterLat)*dLat)**2)/RadialD**2) < 0)then
              if (.not. IsForbidden(iFreq))then
                 write(*,*) 'No longer a upward-propagating AGW at Time = ', &
                      CurrentTime-StartTime, 'iFreq =', iFreq
                 IsForbidden(iFreq) = .true.
              endif
              cycle
           endif

           Kr = sqrt(WaveFreqIntri2*(WaveFreqIntri2 - SoundSpeed2/ &
                (4.*ScaleHeight**2))/SoundSpeed2/(WaveFreqIntri2 - &
                BouyancyFreq2*(((Lon-EpicenterLon)*dLon)**2 + &
                ((Lat-EpicenterLat)*dLat)**2)/RadialD**2))
           Kx = Kr*(Lon-EpicenterLon)*dLon/RadialD
           Ky = Kr*(Lat-EpicenterLat)*dLat/RadialD
           Kz = Kr*Altitude_G(iAlt)/RadialD
           Kh2 = Kx**2 + Ky**2

           ! find region to perturb
           DoPerturb = (RadialD <= PerturbWaveFreq(iFreq)/Kr* &
                (CurrentTime-StartTime-PerturbTimeDelay) &
                .and. RadialD >= PerturbWaveFreq(iFreq)/Kr* &
                (CurrentTime-StartTime-PerturbTimeDelay-PerturbDuration))
        endif

        if (DoPerturb) then 

           Coeff1 = WaveFreqIntri2/ScaleHeight/2.-Kh2*AvgGravity
           Coeff2 = Kh2+Kz**2-1/4./ScaleHeight**2
           Coeff3 = WaveFreqIntri2**2*Kz**2 + Coeff1**2

           if (iTypeBcPerturb < 2) then
              ! plane waves
              WaveForm = Kx*(Lon-RefLon)*dLon + Ky*(Lat-RefLat)*dLat + &
                   Kz*Altitude_G(iAlt) - PerturbWaveFreq(iFreq)* &
                   (CurrentTime-StartTime-PerturbTimeDelay)
              CosineWave = cos(WaveForm)
              SineWave = sin(WaveForm)
              ! FFTImag and FFTReal are for ocean surface vertical displacement
              CoeffS = -2.*PerturbWaveFreq(iFreq)*FFTImag(iFreq)* &
                   exp(Altitude_G(iAlt)/(2.*ScaleHeight))
              CoeffC = 2.*PerturbWaveFreq(iFreq)*FFTReal(iFreq)* &
                   exp(Altitude_G(iAlt)/(2.*ScaleHeight))         
              WaveCombo1 = CoeffS*CosineWave + CoeffC*SineWave
              WaveCombo2 = CoeffS*SineWave + CoeffC*CosineWave

           else if (iTypeBcPerturb == 2) then
              ! spherical waves
               WaveForm = Kr*RadialD - PerturbWaveFreq(iFreq)* &
                    (CurrentTime-StartTime-PerturbTimeDelay)
              CosineWave = cos(WaveForm)
              SineWave = sin(WaveForm)
              ! FFTImag and FFTReal are for the surface vertical velocity
              CoeffS = 2*EpiDistance*exp(pi*0.005*SeisWaveTimeDelay)* &
                   FFTImag(iFreq)*exp(Altitude_G(iAlt)/(2.*ScaleHeight))
              CoeffC = 2*EpiDistance*exp(pi*0.005*SeisWaveTimeDelay)* &
                   FFTReal(iFreq)*exp(Altitude_G(iAlt)/(2.*ScaleHeight))
              WaveCombo1 = (CoeffS*SineWave + CoeffC*CosineWave)/RadialD
              WaveCombo2 = (CoeffS*CosineWave + CoeffC*SineWave)/RadialD
           endif
           
           dRho = exp(-Altitude_G(iAlt)/ScaleHeight)*Rho0*WaveFreqIntri/Coeff3* &
                ((Coeff2*WaveFreqIntri2 + Coeff1/ScaleHeight)*Kz*WaveCombo1 - &
                (WaveFreqIntri2*Kz**2/ScaleHeight - Coeff2*Coeff1)*WaveCombo2)
           ! assume all species are perturbed by the same precentage
           LogNSBc(iAlt,:) = alog(exp(LogNSBc(iAlt,:))*(1 + dRho/exp(LogRhoBc(iAlt))))
           LogRhoBc(iAlt) = alog(exp(LogRhoBc(iAlt)) + dRho)

           dVel = 1/Coeff3* &
                (((WaveFreqIntri2-AvgGravity/ScaleHeight/2.)*WaveFreqIntri2 + &
                AvgGravity*Coeff1)*Kz*WaveCombo1 - &
                (WaveFreqIntri2*Kz**2*AvgGravity - &
                (WaveFreqIntri2-AvgGravity/ScaleHeight/2.)*Coeff1)*WaveCombo2)
           VelBc_GD(iAlt,iEast_) = VelBc_GD(iAlt,iEast_) + Kx*dVel
           VelBc_GD(iAlt,iNorth_) = VelBc_GD(iAlt,iNorth_) + Ky*dVel
           VelBc_GD(iAlt, iUp_) = VelBc_GD(iAlt,iUp_) + WaveCombo1           

           TempBc(iAlt) = TempBc(iAlt) + &
                (Gamma_const-1)*(Temp0 + MeanTempBc)/2.*WaveFreqIntri/ &
                ((Coeff1-(Gamma_const-2)*Kh2*AvgGravity)**2 + &
                WaveFreqIntri2**2*Kz**2)*((Coeff2*WaveFreqIntri2 + &
                (Coeff1-(Gamma_const-2)*Kh2*AvgGravity)/ScaleHeight)*Kz*WaveCombo1 - &
                (Coeff2*(Coeff1-(Gamma_const-2)*Kh2*AvgGravity) - &
                WaveFreqIntri2*Kz**2/ScaleHeight)*WaveCombo2)

           if(iDebugLevel > 0) &
                write(*,*) 'user_bc_perturbation: CurrentTime = ', CurrentTime, &
                'Lat, Lon = ', Lat, Lon, 'iFreq = ', iFreq, 'iAlt = ', iAlt, &
                'Kx, Ky, Kz = ', Kx, Ky, Kz, 'WaveFreqIntri=', WaveFreqIntri, &
                'WaveForm = ', WaveForm, 'CoeffS =', CoeffS,  'CoeffC =', CoeffC, &
                'dRho = ', dRho, 'dVel = ', dVel
                
        endif
     enddo
  enddo

end subroutine user_bc_perturbation


! ----------------------------------------------------------------
! If you want to output some specific variables, then do that here.
! In ModUserGITM, there are two variables defined, UserData2D and UserData3D.
! To output variables:
! 1. Figure out which variable you want to output.
! 2. Go into the code where the variable is set and copy it into
!    UserData3D or UserData2D.
! 3. Do this for each variable you want to output.
! 4. Edit output_header_user below, making a list of all the variables
!    that you have added.  Make sure you leave longitude, latitude, and
!    altitude in the list of variables.
! 5. Count the number of variables that you output (including 
!    Lon, Lat, and Alt). Change nVarsUser3d or nVarsUser2d in the 
!    subroutines towards the top of this file.
! 6. If you add more than 40 variables, you probably should check 
!    nUserOutputs in ModUserGITM.f90 and make sure that this number is
!    larger than the number of variables that you added.
! 7. Recompile and run. Debug. Repeat 7.
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine set_nVarsUser3d

  use ModUserGITM

  ! Make sure to include Lat, Lon, and Alt

  nVarsUser3d = 5

  if (nVarsUser3d-3 > nUserOutputs) &
       call stop_gitm("Too many user outputs!! Increase nUserOutputs!!")

end subroutine set_nVarsUser3d

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine set_nVarsUser2d

  use ModUserGITM
  use ModSources, only:ED_N_Energies

  ! Make sure to include Lat, Lon, and Alt

  nVarsUser2d = 10+ED_N_Energies

  if (nVarsUser2d-3 > nUserOutputs) &
       call stop_gitm("Too many user outputs!! Increase nUserOutputs!!")

end subroutine set_nVarsUser2d

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine set_nVarsUser1d

  use ModUserGITM

  ! Make sure to include Lat, Lon, and Alt

  nVarsUser1d = 4

  if (nVarsUser1d-3 > nUserOutputs) &
       call stop_gitm("Too many user outputs!! Increase nUserOutputs!!")

end subroutine set_nVarsUser1d

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine set_nVarsUser0d

  use ModUserGITM

  nVarsUser0d = 6

  if (nVarsUser0d-3 > nUserOutputs) &
       call stop_gitm("Too many user outputs!! Increase nUserOutputs!!")

end subroutine set_nVarsUser0d

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine output_header_user(cType, iOutputUnit_)

  use ModUserGITM
  use ModSources, only:ED_Energies, ED_N_Energies

  implicit none

  character (len=5), intent(in) :: cType
  integer, intent(in)           :: iOutputUnit_
  integer :: n

  ! ------------------------------------------
  ! 3D Output Header
  ! ------------------------------------------

  if (cType(1:2) == '3D') then 

     write(iOutputUnit_,*) "NUMERICAL VALUES"
     write(iOutputUnit_,"(I7,6A)") nVarsUser3d, " nvars"
     write(iOutputUnit_,"(I7,7A)") nAlts+4, " nAltitudes"
     write(iOutputUnit_,"(I7,7A)") nLats+4, " nLatitudes"
     write(iOutputUnit_,"(I7,7A)") nLons+4, " nLongitudes"

     write(iOutputUnit_,*) ""

     write(iOutputUnit_,*) "VARIABLE LIST"
     write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
     write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
     write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
     write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Joule Heating"
     write(iOutputUnit_,"(I7,A1,a)")  5, " ", "JPara"

  endif

  ! ------------------------------------------
  ! 2D Output Header
  ! ------------------------------------------

  if (cType(1:2) == '2D') then 

     write(iOutputUnit_,*) "NUMERICAL VALUES"
     write(iOutputUnit_,"(I7,6A)") nVarsUser2d, " nvars"
     write(iOutputUnit_,"(I7,7A)")     1, " nAltitudes"
     write(iOutputUnit_,"(I7,7A)") nLats, " nLatitudes"
     write(iOutputUnit_,"(I7,7A)") nLons, " nLongitudes"

     write(iOutputUnit_,*) ""
     write(iOutputUnit_,*) "NO GHOSTCELLS"
     write(iOutputUnit_,*) ""

     write(iOutputUnit_,*) "VARIABLE LIST"
     write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
     write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
     write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
     write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Potential (kV)"
     write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Average Energy (keV)"
     write(iOutputUnit_,"(I7,A1,a)")  6, " ", "Total Energy (ergs)"
     write(iOutputUnit_,"(I7,A1,a)")  7, " ", "Discrete Average Energy (keV)"
     write(iOutputUnit_,"(I7,A1,a)")  8, " ", "Discrete Total Energy (ergs)"
     write(iOutputUnit_,"(I7,A1,a)")  9, " ", "Wave Average Energy (keV)"
     write(iOutputUnit_,"(I7,A1,a)") 10, " ", "Wave Total Energy (ergs)"
     do n=1,ED_N_Energies
        write(iOutputUnit_,"(I7,A6,1P,E9.3,A11)") 10+n, " Flux@",ED_energies(n), "eV (/cm2/s)"
     enddo
  endif

  ! ------------------------------------------
  ! 1D Output Header
  ! ------------------------------------------

  if (cType(1:2) == '1D') then

     write(iOutputUnit_,*) "NUMERICAL VALUES"
     write(iOutputUnit_,"(I7,6A)") nVarsUser1d, " nvars"
     write(iOutputUnit_,"(I7,7A)") nAlts+4, " nAltitudes"
     write(iOutputUnit_,"(I7,7A)") 1, " nLatitudes"
     write(iOutputUnit_,"(I7,7A)") 1, " nLongitudes"

     write(iOutputUnit_,*) "VARIABLE LIST"
     write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
     write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
     write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
     write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Electron Density"
  endif

  ! ------------------------------------------
  ! 0D Output Header
  ! ------------------------------------------

  if (cType(1:2) == '0D') then

     write(iOutputUnit_,*) "NUMERICAL VALUES"
     write(iOutputUnit_,"(I7,6A)") nVarsUser0d, " nvars"
     write(iOutputUnit_,"(I7,7A)") 1, " nAltitudes"
     write(iOutputUnit_,"(I7,7A)") 1, " nLatitudes"
     write(iOutputUnit_,"(I7,7A)") 1, " nLongitudes"

     write(iOutputUnit_,*) "VARIABLE LIST"
     write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
     write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
     write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
     write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Electron Density"
     write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Electron Temperature"
     write(iOutputUnit_,"(I7,A1,a)")  6, " ", "Ion Temperature"
  endif

  write(iOutputUnit_,*) ""

end subroutine output_header_user

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_3dUser(iBlock, iOutputUnit_)

  use ModGITM
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock, iOutputUnit_
  integer :: iAlt, iLat, iLon

  do iAlt=-1,nAlts+2
     do iLat=-1,nLats+2
        do iLon=-1,nLons+2
           write(iOutputUnit_)       &
                Longitude(iLon,iBlock), &
                Latitude(iLat,iBlock), &
                Altitude_GB(iLon, iLat, iAlt, iBlock),&
                UserData3D(iLon,iLat,iAlt,1:nVarsUser3d-3,iBlock)
        enddo
     enddo
  enddo

end subroutine output_3dUser

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_2dUser(iBlock, iOutputUnit_)

  use ModGITM
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock, iOutputUnit_
  integer :: iAlt, iLat, iLon

  iAlt = 1
  do iLat=1,nLats
     do iLon=1,nLons
        write(iOutputUnit_)       &
             Longitude(iLon,iBlock), &
             Latitude(iLat,iBlock), &
             Altitude_GB(iLon, iLat, iAlt, iBlock),&
             UserData2D(iLon,iLat,iAlt,1:nVarsUser2d-3,iBlock)
     enddo
  enddo

end subroutine output_2dUser

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_1dUser(iiLon, iiLat, iBlock, rLon, rLat, iOutputUnit_)

  use ModGITM
  use ModUserGITM
  use ModPlanet, only: ie_

  implicit none

  integer, intent(in) :: iiLat, iiLon, iBlock, iOutputUnit_
  real, intent(in)    :: rLon, rLat
  integer :: iAlt

  do iAlt=-1,nAlts+2
     write(iOutputUnit_)       &
          rLon*Longitude(iiLon,iBlock)+(1-rLon)*Longitude(iiLon+1,iBlock), &
          rLat*Latitude(iiLat,iBlock)+(1-rLat)*Latitude(iiLat+1,iBlock), &
          Altitude_GB(iiLon, iiLat, iAlt, iBlock), &
          inter(IDensityS(0:nLons+1,0:nLats+1,iAlt,ie_,iBlock), &
          iiLon, iiLat, rLon, rLat)
  enddo

  contains

  real function inter(variable, iiLon, iiLat, rLon, rLat) &
       result(PointValue)

    implicit none

    real :: variable(:,:), rLon, rLat
    integer :: iiLon, iiLat

    PointValue = &
         (  rLon)*(  rLat)*Variable(iiLon  ,iiLat  ) + &
         (1-rLon)*(  rLat)*Variable(iiLon+1,iiLat  ) + &
         (  rLon)*(1-rLat)*Variable(iiLon  ,iiLat+1) + &
         (1-rLon)*(1-rLat)*Variable(iiLon+1,iiLat+1)

  end function inter

end subroutine output_1dUser

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_0dUser(iiLon, iiLat, iiAlt, iBlock, rLon, rLat, rAlt, iOutputUnit_)

  use ModGITM
  use ModUserGITM
  use ModPlanet, only: ie_

  implicit none

  integer, intent(in) :: iiLat, iiLon, iiAlt, iBlock, iOutputUnit_
  real, intent(in)    :: rLon, rLat, rAlt

  write(iOutputUnit_)       &
       rLon*Longitude(iiLon,iBlock)+(1-rLon)*Longitude(iiLon+1,iBlock), &
       rLat*Latitude(iiLat,iBlock)+(1-rLat)*Latitude(iiLat+1,iBlock), &
       rAlt*Altitude_GB(iiLon, iiLat, iiAlt, iBlock)+ &
       (1-rAlt)*Altitude_GB(iiLon, iiLat, iiAlt+1, iBlock), &
       inter(IDensityS(0:nLons+1,0:nLats+1,0:nAlts+1,ie_,iBlock), &
       iiLon, iiLat, iiAlt, rLon, rLat, rAlt), &
       inter(eTemperature(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock), &
       iiLon, iiLat, iiAlt, rLon, rLat, rAlt), &
       inter(ITemperature(0:nLons+1,0:nLats+1,0:nAlts+1,iBlock), &
       iiLon, iiLat, iiAlt, rLon, rLat, rAlt)

  contains

    real function inter(variable, iiLon, iiLat, iiAlt, rLon, rLat, rAlt) &
       result(PointValue)

    implicit none

    real :: variable(:,:,:), rLon, rLat, rAlt
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

end subroutine output_0dUser
