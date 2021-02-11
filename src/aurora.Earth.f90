!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine aurora(iBlock)

  use ModGITM
  use ModSources
  use ModTime, only : tSimulation, CurrentTime
  use ModInputs
  use ModConstants
  use ModUserGITM
  use ModMpi
  use ModIndicesInterfaces

  implicit none

  integer, intent(in) :: iBlock

  real :: alat, hpi, ped, hal, av_kev, eflx_ergs, a,b, maxi
  real :: Factor,temp_ED, avee, eflux, p, E0, Q0
  integer :: i, j, k, n, iError, iED, iErr, iEnergy
  logical :: IsDone, IsTop, HasSomeAurora, UseMono, UseWave

  real, dimension(nLons,nLats,nAlts) :: temp, AuroralBulkIonRate, &
       IonPrecipitationBulkIonRate, IonPrecipitationHeatingRate

  real :: fac(nAlts)
  
  logical :: IsFirstTime(nBlocksMax) = .true.

  real :: f1, f2, f3, f4, f5, power
  real :: de1, de2, de3, de4, de5, detotal, h

  real :: LocalVar, HPn, HPs, avepower, ratio

  real :: Fang_Pij(8,4), Ci(8), Fang_de = 0.035
  data Fang_Pij(1,:) /1.25E+00,1.45903,-2.42E-01,5.95E-02/
  data Fang_Pij(2,:) /2.24E+00,-4.23E-07,1.36E-02,2.53E-03/
  data Fang_Pij(3,:) /1.42E+00,1.45E-01,1.70E-02,6.40E-04/
  data Fang_Pij(4,:) /0.248775,-1.51E-01,6.31E-09,1.24E-03/
  data Fang_Pij(5,:) /-0.465119,-1.05E-01,-8.96E-02,1.22E-02/
  data Fang_Pij(6,:) /3.86E-01,1.75E-03,-7.43E-04,4.61E-04/
  data Fang_Pij(7,:) /-6.45E-01,8.50E-04,-4.29E-02,-2.99E-03/
  data Fang_Pij(8,:) /9.49E-01,1.97E-01,-2.51E-03,-2.07E-03/

  real BulkScaleHeight1d(nAlts)

  HPn = 0.0
  HPs = 0.0

  if (IsFirstTime(iBlock)) then
     
     IsFirstTime(iBlock) = .false.

     if (iBlock == 1 .and. UseIonPrecipitation) then
        call ReadIonHeat(IonIonizationFilename,  .true.) 
        call ReadIonHeat(IonHeatingRateFilename, .false.) 
     endif

     if (UseFangEnergyDeposition) then
        allocate(Fang_Ci(ED_N_Energies,8), stat=iErr)
        allocate(Fang_y(ED_N_Energies,nAlts), stat=iErr)
        allocate(Fang_f(ED_N_Energies,nAlts), stat=iErr)
        if (iErr /= 0) then
           call stop_gitm("Error allocating Fang arrays in aurora")
        endif

        do iEnergy = 1, ED_N_Energies
           do i=1,8
              Fang_Ci(iEnergy,i) = 0.0
              do j=0,3
                 Fang_Ci(iEnergy,i) = Fang_Ci(iEnergy,i) + &
                      Fang_Pij(i,j+1) * log(ED_Energies(iEnergy)/1000.0)**j
              enddo
           enddo
        enddo
        Fang_Ci = exp(Fang_Ci)
     endif

  else
     if (floor((tSimulation - dT)/dTAurora) == &
          floor(tSimulation/dTAurora)) return
  endif

  AuroralBulkIonRate               = 0.0
  AuroralHeatingRate(:,:,:,iBlock) = 0.0

  call report("Aurora",1)
  call start_timing("Aurora")

  if (UseIonPrecipitation) call interpolate_ions( &
       nLons, nLats, nAlts, &
       Longitude(1:nLons,iBlock), Latitude(1:nLats,iBlock), &
       Altitude_GB(1:nLons, 1:nLats, 1:nAlts, iBlock),&
       IonPrecipitationBulkIonRate, IonPrecipitationHeatingRate)

  if (iBlock == 1) then
     HemisphericPowerNorth = 0.0
     HemisphericPowerSouth = 0.0
  endif

  ! Let's scale our hemispheric power so it is roughly the same as what
  ! is measured.

  if (NormalizeAuroraToHP) then
  
     do i=1,nLats
        do j=1,nLons

           eflx_ergs = ElectronEnergyFlux(j,i) !/ (1.0e-7 * 100.0 * 100.0)

           if (eflx_ergs > 0.1) then
              eflux = eflx_ergs * 6.242e11  ! ergs/cm2/s -> eV/cm2/s

              !(eV/cm2/s -> J/m2/s)
              power = eflux * Element_Charge*100.0*100.0 * & 
                   dLatDist_FB(j, i, nAlts, iBlock) * &
                   dLonDist_FB(j, i, nAlts, iBlock)

              if (latitude(i,iBlock) < 0.0) then
                 HemisphericPowerSouth = HemisphericPowerSouth + power
              else
                 HemisphericPowerNorth = HemisphericPowerNorth + power
              endif

           endif

        enddo
     enddo

     ! Collect all of the powers by summing them together

     LocalVar = HemisphericPowerNorth/1.0e9
     call MPI_REDUCE(LocalVar, HPn, 1, MPI_REAL, MPI_SUM, &
          0, iCommGITM, iError)

     LocalVar = HemisphericPowerSouth/1.0e9
     call MPI_REDUCE(LocalVar, HPs, 1, MPI_REAL, MPI_SUM, &
          0, iCommGITM, iError)

     ! Average north and south together

     avepower = (HPn+HPs)/2.0

     ! If we are only have one hemisphere or the other, assign to avepower
     if (HPs < 0.1*HPn) avepower = HPn
     if (HPn < 0.1*HPs) avepower = HPs

     call MPI_Bcast(avepower,1,MPI_Real,0,iCommGITM,ierror)

     call get_hpi(CurrentTime,Hpi,iError)
     ratio = Hpi/avepower

     write(*,*) 'Auroral normalizing ratio : ', Hpi, avepower, ratio
     
     do i=1,nLats
        do j=1,nLons
           if (ElectronEnergyFlux(j,i)>0.1) then
              ElectronEnergyFlux(j,i) = ElectronEnergyFlux(j,i)*ratio
           endif
        enddo
     enddo

  endif

  ! Reset the hemispheric power

  if (iBlock == 1) then
     HemisphericPowerNorth = 0.0
     HemisphericPowerSouth = 0.0
  endif

  AveEFactor = 1.0
  if (iProc == 0 .and. AveEFactor /= 1.0) then
     write(*,*) "Auroral Experiments!!!!"
     write(*,*) "AveEFactor : ", AveEFactor
  endif
  if (iProc == 0 .and. IsKappaAurora) then
     write(*,*) "Auroral Experiments!!!!"
     write(*,*) "kappa : ", AuroraKappa
  endif
  
  do i=1,nLats
     do j=1,nLons

        UserData2d(j,i,1,2:nUserOutputs,iBlock) = 0.0

        eflx_ergs = ElectronEnergyFlux(j,i) !/ (1.0e-7 * 100.0 * 100.0)
        av_kev    = ElectronAverageEnergy(j,i)

        ! For diffuse auroral models

        ED_Flux = 0.0
        HasSomeAurora = .false.

        if (eflx_ergs > 0.1) then

           UserData2d(j,i,1,2,iBlock) = av_kev
           UserData2d(j,i,1,3,iBlock) = eflx_ergs

           HasSomeAurora = .true.
           avee = av_kev * 1000.0        ! keV -> eV
           eflux = eflx_ergs * 6.242e11  ! ergs/cm2/s -> eV/cm2/s

           power = eflux * Element_Charge*100.0*100.0 * & !(eV/cm2/s -> J/m2/s)
                dLatDist_FB(j, i, nAlts, iBlock) * &
                dLonDist_FB(j, i, nAlts, iBlock)

           if (latitude(i,iBlock) < 0.0) then
              HemisphericPowerSouth = HemisphericPowerSouth + power
           else
              HemisphericPowerNorth = HemisphericPowerNorth + power
           endif

           ! Looking at other papers (Fang et al., [2010]), a
           ! Maxwellian is defined as:
           ! DifferentialNumberFlux = Q0/2/E0**3 * E * exp(-E/E0),
           ! where:
           ! Q0 = Total Energy Flux
           ! E0 = Characteristic Energy (0.5*avee)
           ! E = mid-point of energy bin
           !
           Q0 = eflux
           E0 = avee/2
           a = Q0/2/E0**3

           do n=1,ED_N_Energies


              if (IsKappaAurora) then
                 ! This is a Kappa Function from Fang et al. [2010]:
                 ED_Flux(n) = a * (AuroraKappa-1) * (AuroraKappa-2) / &
                      (AuroraKappa**2) * &
                      ed_energies(n) * &
                      (1 + ed_energies(n) / (AuroraKappa * E0)) ** (-k-1)
              else
                 ! This is a Maxwellian from Fang et al. [2010]:
                 ED_flux(n) = a * ed_energies(n) * exp(-ed_energies(n)/E0)
              endif

              ED_EnergyFlux(n) = &
                   ED_flux(n) * &
                   ED_Energies(n) * &
                   ED_delta_energy(n)
              
           enddo

        endif

        UseMono = .false.
        if (UseNewellAurora .and. UseNewellMono    ) UseMono = .true.
        if (UseOvationSME   .and. UseOvationSMEMono) UseMono = .true.

        if (UseMono .and. ElectronNumberFluxMono(j,i) > 0.0) then

           av_kev = ElectronEnergyFluxMono(j, i) / &
                    ElectronNumberFluxMono(j, i) * 6.242e11 ! eV

           power = ElectronNumberFluxMono(j, i) * &
                Element_Charge * 100.0 * 100.0 * &    ! (eV/cm2/s -> J/m2/s)
                dLatDist_FB(j, i, nAlts, iBlock) * &
                dLonDist_FB(j, i, nAlts, iBlock)

           if (latitude(i,iBlock) < 0.0) then
              HemisphericPowerSouth = HemisphericPowerSouth + power
           else
              HemisphericPowerNorth = HemisphericPowerNorth + power
           endif

           UserData2d(j,i,1,4,iBlock) = av_kev / 1000.0
           UserData2d(j,i,1,5,iBlock) = ElectronEnergyFluxMono(j, i)

           ! Mono-Energetic goes into one bin only!
           do n=2,ED_N_Energies-1
              if (av_kev < ED_energies(n-1) .and. av_kev >= ED_energies(n)) then
                 ED_flux(n) = ED_Flux(n) + &
                      ElectronNumberFluxMono(j, i) / &
                      (ED_Energies(n-1) - ED_Energies(n))
                 HasSomeAurora = .true.
              endif
           enddo

        endif

        UseWave = .false.
        if (UseNewellAurora .and. UseNewellWave    ) UseWave = .true.
        if (UseOvationSME   .and. UseOvationSMEWave) UseWave = .true.

        if (UseWave .and. ElectronNumberFluxWave(j,i) > 0.0) then

           av_kev = ElectronEnergyFluxWave(j, i) / &
                    ElectronNumberFluxWave(j, i) * 6.242e11 ! eV

           power = ElectronNumberFluxWave(j, i) * &
                Element_Charge * 100.0 * 100.0 * &    ! (eV/cm2/s -> J/m2/s)
                dLatDist_FB(j, i, nAlts, iBlock) * dLonDist_FB(j, i, nAlts, iBlock)

           if (latitude(i,iBlock) < 0.0) then
              HemisphericPowerSouth = HemisphericPowerSouth + power
           else
              HemisphericPowerNorth = HemisphericPowerNorth + power
           endif

           UserData2d(j,i,1,6,iBlock) = av_kev / 1000.0
           UserData2d(j,i,1,7,iBlock) = ElectronEnergyFluxWave(j, i)

           ! Waves goes into five bins only!
           k = 0
           do n=3,ED_N_Energies-3
              if (av_kev < ED_energies(n-1) .and. av_kev >= ED_energies(n)) then
                 k = n
              endif
           enddo
           if (k > 3) then 
              f1 = 1.0
              f2 = 1.2
              f3 = 1.3
              f4 = f2
              f5 = f1
              de1 = ED_energies(k-3)-ED_energies(k-2)
              de2 = ED_energies(k-2)-ED_energies(k-1)
              de3 = ED_energies(k-1)-ED_energies(k)  
              de4 = ED_energies(k)  -ED_energies(k+1)
              de5 = ED_energies(k+1)-ED_energies(k+2)
!              detotal = (de1+de2+de3+de4+de5) * (f1+f2+f3+f4+f5) / 5
              detotal = (f1+f2+f3+f4+f5)
              ED_flux(k-2) = ED_Flux(k-2)+f1*ElectronNumberFluxWave(j, i)/detotal/de1
              ED_flux(k-1) = ED_Flux(k-1)+f2*ElectronNumberFluxWave(j, i)/detotal/de2
              ED_flux(k  ) = ED_Flux(k  )+f3*ElectronNumberFluxWave(j, i)/detotal/de3
              ED_flux(k+1) = ED_Flux(k+1)+f4*ElectronNumberFluxWave(j, i)/detotal/de4
              ED_flux(k+2) = ED_Flux(k+2)+f5*ElectronNumberFluxWave(j, i)/detotal/de5
!              ED_flux(k-2) = ED_Flux(k-2) + f1*ElectronNumberFluxWave(j, i) / detotal
!              ED_flux(k-1) = ED_Flux(k-1) + f2*ElectronNumberFluxWave(j, i) / detotal
!              ED_flux(k  ) = ED_Flux(k  ) + f3*ElectronNumberFluxWave(j, i) / detotal
!              ED_flux(k+1) = ED_Flux(k+1) + f4*ElectronNumberFluxWave(j, i) / detotal
!              ED_flux(k+2) = ED_Flux(k+2) + f5*ElectronNumberFluxWave(j, i) / detotal
           endif

        endif

        if (HasSomeAurora) then

           if (UseFangEnergyDeposition) then

              BulkScaleHeight1d = &
                      Temperature(j,i,1:nAlts,iBlock) &
                      * TempUnit(j,i,1:nAlts) * Boltzmanns_Constant &
                      / (-Gravity_GB(j,i,1:nAlts,iBlock) * &
                      MeanMajorMass(j,i,1:nAlts))*100.0 ! Convert to cm
              
              do iEnergy = 1,ED_N_Energies

                 ! /10.0 in this statement is for kg/m2 to g/cm2
                 ! /1000.0 is conversion from eV to keV
                 ! Fang doesn't include the dip angle, be we do.
                 Fang_y(iEnergy,:) = 2.0 / (ED_Energies(iEnergy)/1000.0) * &
                      (ColumnIntegralRho(j,i,:) / 10.0 / 6e-6) ** 0.7
                      !sinDipAngle(j,i,1:nAlts,iBlock) / 6e-6) ** 0.7

                 Ci = Fang_Ci(iEnergy,:)
                 ! write(*,*) iEnergy, ED_Energies(iEnergy), Ci
                 Fang_f(iEnergy,:) = &
                      Ci(1) * Fang_y(iEnergy,:) ** Ci(2) * &
                      exp(-Ci(3) * Fang_y(iEnergy,:) ** Ci(4)) + &
                      Ci(5) * Fang_y(iEnergy,:) ** Ci(6) * & 
                      exp(-Ci(7) * Fang_y(iEnergy,:) ** Ci(8)) 

                 ! Energy flux is in eV/cm2/s and Fang needs keV/cm2/s:
                 fac = ED_energyflux(iEnergy)/1000.0 / &
                      Fang_de / &
                      BulkScaleHeight1d

                 AuroralBulkIonRate(j,i,1:nAlts) = &
                      AuroralBulkIonRate(j,i,1:nAlts) + 1e6*Fang_f(iEnergy,:) * fac

              enddo
              !do k = 1, nAlts
              !   write(*,*) k,  AuroralBulkIonRate(j,i,k), fac(k), fang_f(25,k), ColumnIntegralRho(j,i,k)
              !enddo
              AuroralHeatingRate(j,i,1:nAlts,iBlock) = 0.0

              !write(*,*) AuroralBulkIonRate(j,i,1:nAlts)
              
           else

              call R_ELEC_EDEP (ED_Flux, 15, ED_Energies, 3, ED_Ion, 7)
              call R_ELEC_EDEP (ED_Flux, 15, ED_Energies, 3, ED_Heating, 11)

              iED = 1

              factor = 1.0

              do k = 1, nAlts

                 p = alog(Pressure(j,i,k,iBlock)*factor)

                 IsDone = .false.
                 IsTop = .false.
                 do while (.not.IsDone)
                    if (ED_grid(iED) >= p .and. ED_grid(iED+1) <= p) then
                       IsDone = .true.
                       ED_Interpolation_Index(k) = iED
                       ED_Interpolation_Weight(k) = (ED_grid(iED) - p) /  &
                            (ED_grid(iED) - ED_grid(iED+1))
                    else
                       if (iED == ED_N_Alts-1) then
                          IsDone = .true.
                          IsTop = .true.
                       else
                          iED = iED + 1
                       endif
                    endif
                 enddo

                 if (.not.IsTop) then
                    n = ED_Interpolation_Index(k)
                    AuroralBulkIonRate(j,i,k) = ED_Ion(n) - &
                         (ED_Ion(n) - ED_Ion(n+1))*ED_Interpolation_Weight(k)
                    AuroralHeatingRate(j,i,k,iBlock) = ED_Heating(n) - &
                         (ED_Heating(n) - ED_Heating(n+1))*ED_Interpolation_Weight(k)
                 else

                    ! Decrease after top of model
                    AuroralBulkIonRate(j,i,k) = ED_Ion(ED_N_Alts) * &
                         factor*Pressure(j,i,k,iBlock) / exp(ED_grid(ED_N_Alts)) 
                    AuroralHeatingRate(j,i,k,iBlock) = ED_Heating(ED_N_Alts) * &
                         factor*Pressure(j,i,k,iBlock) / exp(ED_grid(ED_N_Alts))

                 endif

              enddo

           endif

        endif

     enddo
  enddo

  ! From Rees's book:

  temp = 0.92 * NDensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock) + &
         1.00 * NDensityS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock) + &
         0.56 * NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)

  if (UseIonPrecipitation) then

     IonPrecipIonRateS(:,:,:,iO_3P_,iBlock)  = &
          0.56*IonPrecipitationBulkIonRate*&
          NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)/temp
     IonPrecipIonRateS(:,:,:,iO2_,iBlock) = &
          1.00*IonPrecipitationBulkIonRate*&
          NDensityS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock)/temp
     IonPrecipIonRateS(:,:,:,iN2_,iBlock) = &
          0.92*IonPrecipitationBulkIonRate*&
          NDensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock)/temp

     IonPrecipHeatingRate(:,:,:,iBlock) = IonPrecipitationHeatingRate

  else

     IonPrecipIonRateS(:,:,:,:,iBlock) = 0.0
     IonPrecipHeatingRate(:,:,:,iBlock) = 0.0

  endif

  AuroralIonRateS(:,:,:,iO_3P_,iBlock)  = &
       0.56*AuroralBulkIonRate*&
       NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)/temp
  AuroralIonRateS(:,:,:,iO2_,iBlock) = &
       1.00*AuroralBulkIonRate*&
       NDensityS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock)/temp
  AuroralIonRateS(:,:,:,iN2_,iBlock) = &
       0.92*AuroralBulkIonRate*&
       NDensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock)/temp

  call end_timing("Aurora")

end subroutine aurora




