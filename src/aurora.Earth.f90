! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

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
  real :: ion_av_kev, ion_eflx_ergs, ion_eflux, ion_avee
  real :: Factor,temp_ED, avee, eflux, p, E0, Q0, E0i, Q0i, ai
  integer :: i, j, k, n, iAlt, iError, iED, iErr, iEnergy
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

  ! This is from:
  ! Fang, X., D. Lummerzheim, and C. H. Jackman (2013),
  !           Proton impact ionization and a fast calculation method,
  !           J. Geophys. Res. Space Physics, 118, 5369–5378,
  !           doi:10.1002/jgra.50484:

  real :: Fang_Ion_Pij(12,4), Ion_Ci(12)
  data Fang_Ion_Pij( 1, :) / 2.55050E+00,  2.69476e-01, -2.58425E-01,  4.43190E-02/
  data Fang_Ion_Pij( 2, :) / 6.39287E-01, -1.85817e-01, -3.15636E-02,  1.01370E-02/
  data Fang_Ion_Pij( 3, :) / 1.63996E+00,  2.43580e-01,  4.29873E-02,  3.77803E-02/
  data Fang_Ion_Pij( 4, :) /-2.13479E-01,  1.42464e-01,  1.55840E-02,  1.97407E-03/
  data Fang_Ion_Pij( 5, :) /-1.65764E-01,  3.39654e-01, -9.87971E-03,  4.02411E-03/
  data Fang_Ion_Pij( 6, :) /-3.59358E-02,  2.50330e-02, -3.29365E-02,  5.08057E-03/
  data Fang_Ion_Pij( 7, :) /-6.26528E-01,  1.46865e+00,  2.51853E-01, -4.57132E-02/
  data Fang_Ion_Pij( 8, :) / 1.01384E+00,  5.94301e-02, -3.27839E-02,  3.42688E-03/
  data Fang_Ion_Pij( 9, :) /-1.29454E-06, -1.43623e-01,  2.82583E-01,  8.29809E-02/
  data Fang_Ion_Pij(10, :) /-1.18622E-01,  1.79191e-01,  6.49171E-02, -3.99715E-03/
  data Fang_Ion_Pij(11, :) / 2.94890E+00, -5.75821e-01,  2.48563E-02,  8.31078E-02/
  data Fang_Ion_Pij(12, :) /-1.89515E-01,  3.53452e-02,  7.77964E-02, -4.06034E-03/
  
  real BulkScaleHeight1d(nAlts)

  HPn = 0.0
  HPs = 0.0

  if (IsFirstTime(iBlock)) then
     
     if (UseFangEnergyDeposition) then

        ! Electrons
        allocate(Fang_Ci(ED_N_Energies,8), stat=iErr)
        allocate(Fang_y(ED_N_Energies,nAlts), stat=iErr)
        allocate(Fang_f(ED_N_Energies,nAlts), stat=iErr)

        ! Ions
        if (UseIonPrecipitation) then
           allocate(Fang_Ion_Ci(ED_N_Energies,12), stat=iErr)
           allocate(Fang_Ion_y(ED_N_Energies,nAlts), stat=iErr)
           allocate(Fang_Ion_f(ED_N_Energies,nAlts), stat=iErr)
        endif
        
        if (iErr /= 0) then
           call stop_gitm("Error allocating Fang arrays in aurora")
        endif

        ! Electrons
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

        ! Ions
        if (UseIonPrecipitation) then
           do iEnergy = 1, ED_N_Energies
              do i=1,12
                 Fang_Ion_Ci(iEnergy,i) = 0.0
                 do j=0,3
                    Fang_Ion_Ci(iEnergy,i) = Fang_Ion_Ci(iEnergy,i) + &
                         Fang_Ion_Pij(i,j+1) * log(ED_Energies(iEnergy)/1000.0)**j
                 enddo
              enddo
           enddo
           Fang_Ion_Ci = exp(Fang_Ion_Ci)
        endif
        
     endif

  else
     if (floor((tSimulation - dT)/dTAurora) == &
          floor(tSimulation/dTAurora)) return
  endif

  AuroralBulkIonRate               = 0.0
  AuroralHeatingRate(:,:,:,iBlock) = 0.0
  AuroralIonRateS = 0.0

  call report("Aurora",1)
  call start_timing("Aurora")

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

     if (iDebugLevel >= 0) then
        if ((iDebugLevel == 0) .and. IsFirstTime(iBlock)) then
           write(*,*) '---------------------------------------------------'
           write(*,*) 'Using auroral normalizing ratios!!! '
           write(*,*) 'no longer reporting!'
           write(*,*) '---------------------------------------------------'
        else
           if (iDebugLevel >= 1) &
              write(*,*) 'auroral normalizing ratio: ', Hpi, avepower, ratio
        endif
     endif
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

        eflx_ergs = ElectronEnergyFlux(j,i)
        av_kev    = ElectronAverageEnergy(j,i)

        if (UseIonPrecipitation) then
           ion_eflx_ergs = IonEnergyFlux(j,i)
           ion_av_kev = IonAverageEnergy(j,i)
        else
           ion_eflx_ergs = 0.001
           ion_av_kev = 10.0
        endif
        
        ! For diffuse auroral models

        ED_Flux = 0.0
        HasSomeAurora = .false.

        if (eflx_ergs > 0.1) then

           UserData2d(j,i,1,2,iBlock) = av_kev
           UserData2d(j,i,1,3,iBlock) = eflx_ergs

           HasSomeAurora = .true.
           avee = av_kev * 1000.0        ! keV -> eV
           eflux = eflx_ergs * 6.242e11  ! ergs/cm2/s -> eV/cm2/s

           ion_avee = ion_av_kev * 1000.0        ! keV -> eV
           ion_eflux = ion_eflx_ergs * 6.242e11  ! ergs/cm2/s -> eV/cm2/s
           
           ! 100 * 100 is for (eV/cm2/s -> J/m2/s)
           power = (eflux+ion_eflux) * Element_Charge*100.0*100.0 * &
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

           Q0i = ion_eflux
           E0i = ion_avee/2
           ai = Q0i/2/E0i**3

           do n=1,ED_N_Energies

              if (IsKappaAurora) then
                 ! This is a Kappa Function from Fang et al. [2010]:
                 ED_Flux(n) = a * (AuroraKappa-1) * (AuroraKappa-2) / &
                      (AuroraKappa**2) * &
                      ed_energies(n) * &
                      (1 + ed_energies(n) / (AuroraKappa * E0)) ** &
                      (-AuroraKappa-1)
              else
                 ! This is a Maxwellian from Fang et al. [2010]:
                 ED_flux(n) = a * ed_energies(n) * exp(-ed_energies(n)/E0)
                 ED_ion_flux(n) = ai * ed_energies(n) * exp(-ed_energies(n)/E0i)
              endif

              ED_EnergyFlux(n) = &
                   ED_flux(n) * &
                   ED_Energies(n) * &
                   ED_delta_energy(n)
              ED_Ion_EnergyFlux(n) = &
                   ED_Ion_flux(n) * &
                   ED_Energies(n) * &
                   ED_delta_energy(n)
              
           enddo

        endif

        UseMono = .false.
        if (UseNewellAurora .and. UseNewellMono    ) UseMono = .true.
        if (UseOvationSME   .and. UseOvationSMEMono) UseMono = .true.

        if ( UseMono .and. &
             ElectronNumberFluxMono(j,i) > 1.0e4 .and. &
             ElectronEnergyFluxMono(j, i) > 0.1) then

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
                 ED_EnergyFlux(n) = &
                      ED_flux(n) * &
                      ED_Energies(n) * &
                      ED_delta_energy(n)
                 HasSomeAurora = .true.
              endif
           enddo

        endif

        UseWave = .false.
        if (UseNewellAurora .and. UseNewellWave    ) UseWave = .true.
        if (UseOvationSME   .and. UseOvationSMEWave) UseWave = .true.

        if ( UseWave .and. &
             ElectronNumberFluxWave(j,i) > 1.0e4 .and. &
             ElectronEnergyFluxWave(j,i) > 0.1) then

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
              do n = k-2, n+2
                 ED_EnergyFlux(n) = ED_flux(n) * ED_Energies(n) * ED_delta_energy(n)
              enddo
              HasSomeAurora = .true.
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
                      (ColumnIntegralRho(j,i,1:nAlts) / 10.0 / 6e-6) ** 0.7
                      !sinDipAngle(j,i,1:nAlts,iBlock) / 6e-6) ** 0.7

                 Ci = Fang_Ci(iEnergy,:)
                 Fang_f(iEnergy,:) = &
                      Ci(1) * Fang_y(iEnergy,:) ** Ci(2) * &
                      exp(-Ci(3) * Fang_y(iEnergy,:) ** Ci(4)) + &
                      Ci(5) * Fang_y(iEnergy,:) ** Ci(6) * & 
                      exp(-Ci(7) * Fang_y(iEnergy,:) ** Ci(8)) 

                 ! Energy flux is in eV/cm2/s and Fang needs keV/cm2/s:
                 fac = ED_energyflux(iEnergy)/1000.0 / &
                      Fang_de / &
                      BulkScaleHeight1d

                 ! I think that the 1e6 is cm3 to m3
                 AuroralBulkIonRate(j,i,1:nAlts) = &
                      AuroralBulkIonRate(j,i,1:nAlts) + 1e6*Fang_f(iEnergy,:) * fac

                 if (UseIonPrecipitation) then

                    ! /10.0 in this statement is for kg/m2 to g/cm2
                    ! /1000.0 is conversion from eV to keV
                    Fang_Ion_y(iEnergy,:) = 7.5 / (ED_Energies(iEnergy)/1000.0) * &
                         (ColumnIntegralRho(j,i,1:nAlts) / 10.0 / 1e-4) ** 0.9

                    Ion_Ci = Fang_Ion_Ci(iEnergy,:)

                    Fang_Ion_f(iEnergy,:) = &
                         Ion_Ci(1) * Fang_Ion_y(iEnergy,:) ** Ion_Ci(2) * &
                         exp(-Ion_Ci(3) * Fang_Ion_y(iEnergy,:) ** Ion_Ci(4)) + &
                         Ion_Ci(5) * Fang_Ion_y(iEnergy,:) ** Ion_Ci(6) * & 
                         exp(-Ion_Ci(7) * Fang_Ion_y(iEnergy,:) ** Ion_Ci(8)) + &
                         Ion_Ci(9) * Fang_Ion_y(iEnergy,:) ** Ion_Ci(10) * & 
                         exp(-Ion_Ci(11) * Fang_Ion_y(iEnergy,:) ** Ion_Ci(12)) 

                    fac = ED_ion_energyflux(iEnergy)/1000.0 / &
                         Fang_de / &
                         BulkScaleHeight1d

                    AuroralBulkIonRate(j,i,1:nAlts) = &
                         AuroralBulkIonRate(j,i,1:nAlts) + 1e6*Fang_Ion_f(iEnergy,:) * fac
                    
                 endif
                 
              enddo

              AuroralHeatingRate(j,i,1:nAlts,iBlock) = 0.0
              
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

  AuroralIonRateS(:,:,:,iO_3P_,iBlock)  = &
       0.56*AuroralBulkIonRate*&
       NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)/temp
  AuroralIonRateS(:,:,:,iO2_,iBlock) = &
       1.00*AuroralBulkIonRate*&
       NDensityS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock)/temp
  AuroralIonRateS(:,:,:,iN2_,iBlock) = &
       0.92*AuroralBulkIonRate*&
       NDensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock)/temp

  if (UseAuroralHeating) then
     AuroralHeating = AuroralHeatingRate(:,:,:,iBlock) / &
          TempUnit(1:nLons,1:nLats,1:nAlts) / cp(:,:,1:nAlts,iBlock) / &
          rho(1:nLons,1:nLats,1:nAlts, iBlock)
  else
     AuroralHeating = 0.0
  endif

  IsFirstTime(iBlock) = .false.

  call end_timing("Aurora")

end subroutine aurora




