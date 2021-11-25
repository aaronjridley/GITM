! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine calc_ir_heating(iBlock)

  use ModInputs
  use ModSources
  use ModGITM
  use ModTime, only: UTime
  use ModEUV, only: sza
  implicit none

  integer, intent(in) :: iBlock
  integer :: iErr, iAlt, iLon, iLat, jAlt, iSza
  real :: m !slope for IR heating

  !BP: Proportionality for IR interpolation in sza and in altitude
  real :: r_sza
  real :: x03, x12
  real :: LST
  real :: RP_hours = Rotation_Period/3600.0
  integer :: iLatTable, iAltTable

  !BP: sza data for IR table                                                                 
  real, dimension(11) :: sza_table = (/0.0,25.0,36.0,45.0,53.0,60.0, &
                                      66.0,72.0,78.0,84.0,89.0 /)*pi/180.0
  real, dimension(40) :: altitude_table = (/80.0,82.05,84.10,86.15,88.21, &
                                            90.26,92.31,94.36,96.41,98.46, &
                                            100.51,102.56,104.62,106.67,108.72, &
                                            110.77,112.82,114.87,116.92,118.97, &
                                            121.03,123.08,125.13,127.18,129.23, &
                                            131.28,133.33,135.38,137.44,139.49, &
                                            141.54,143.59,145.64,147.69,149.74, &
                                            151.79,153.85,155.90,157.95,160.00/)*1000.0

  !BP (5/29/2020)                                                             
  !Compute the IR heating at every location in GITM.                             
  !Inputs:   SZA & Altitude                                                             
  !Outputs:  Q_IR (iLon, iLat, iAlt, iProc)
  call start_timing("IR Heating")
  do iAlt = 1, nAlts
     do iLon = 1, nLons
        do iLat = 1, nLats
           !Everywhere we don't want heating                                        
           !1.) Below 80 km                                                          
           !2.) Above 160 km                                                        
           !3.) SZA < 0                                                             
           !4.) SZA > 90                                                      
           if (Altitude_GB(iLon,iLat,iAlt,iBlock) < 80.0*1000.0 .or. &
                Altitude_GB(iLon,iLat,iAlt,iBlock) > 160.0*1000.0 .or. &
                sza(iLon,iLat,iBlock) < 0.0 .or. &
                sza(iLon,iLat,iBlock) .ge. 89.0*pi/180.0) then
              QnirTOT(iLon,iLat,iAlt,iBlock) = 0.0
           else
              do iSza = 1, 11
                 if (sza(iLon, iLat, iBlock) .le. sza_table(iSza+1) .and. &
                      sza(iLon, iLat, iBlock) .ge. sza_table(iSza)) then
                    exit
                 endif
              enddo

              do jAlt = 1,40
                 if (Altitude_GB(iLon,iLat,iAlt,iBlock) .le. altitude_table(jAlt+1) .and. &
                      Altitude_GB(iLon,iLat,iAlt,iBlock) .ge. altitude_table(jAlt)) then
                    exit
                 endif
              enddo

              !Debugging stuff                                                      
              !if (iProc == 0) then                                                
              ! write(*,*) iLon, iLat, iAlt, nLons                                    
              !write(*,*) "GITM SZA", sza(iLon,iLat, iBlock)                               
              !write(*,*) "SZA Found in the table:", sza_table(i)                          
              !write(*,*) "Altitude Found in table:", altitude_table(j)              
              !write(*,*) "Longitude", Longitude(iLon, iBlock) * 180.0/pi             
              !write(*,*) "Latitude", Latitude(iLat, iBlock) * 180.0/pi                
              !write(*,*) "Altitude",  Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0      
              !endif                                                  
              r_sza = (sza(iLon, iLat, iBlock) - sza_table(iSza)) / &
                   (sza_table(iSza+1) - sza_table(iSza))

              m = (altitude_table(jAlt) - altitude_table(jAlt+1))/&
                   (qIR_NLTE_table(jAlt,iSza) - qIR_NLTE_table(jAlt+1,iSza))

              x03 = (Altitude_GB(iLon,iLat,iAlt,iBlock) - &
                   (altitude_table(jAlt+1) - m*qIR_NLTE_table(jAlt+1,iSza))) / m
              m = (altitude_table(jAlt)- altitude_table(jAlt+1))/&
                   (qIR_NLTE_table(jAlt,iSza+1) - qIR_NLTE_table(jAlt+1,iSza+1))

              x12 = (Altitude_GB(iLon,iLat,iAlt,iBlock) - &
                   (altitude_table(jAlt+1) - m*qIR_NLTE_table(jAlt+1,iSza+1))) / m

              if (r_sza > 1.0) then
                 write(*,*) "r is too big...", iLon, iLat, iAlt
                 write(*,*) "GITM SZA:", sza(iLon, iLat, iBlock)
                 write(*,*) "Table SZA:", sza_table(iSza), sza_table(iSza+1)
              endif
                 
              if (x03 < 0.0 .or. x12 < 0.0) then
                 write(*,*) "First interpolated value is negative...this is wrong."
                 write(*,*) "x03", x03
                 write(*,*) "x12", x12
              endif

              QnirTOT(iLon, iLat, iAlt, iBlock) = &
                   x03*(1 - r_sza) + x12*r_sza

              !Given in K/s. No need for rho or cp because that             
              !converts J/s -> K/s. Need TempUnit to normalize it to               
              !GITM units                                                        
              QnirTOT(iLon, iLat, iAlt, iBlock) = QnirTOT(iLon, iLat, iAlt, iBlock)  &
                   / TempUnit(iLon,iLat,iAlt)! &        
              ! / Rho(iLon,iLat,iAlt, iBlock) &       
              ! / cp(iLon,iLat,iAlt,iBlock)        
           endif

           !LTE IR Heating Section                                                    
           if (Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0 < 100.0) then
              LST = mod(UTime/3600.0 + Longitude(iLon, iBlock)*180.0/(360*pi/RP_hours), &
                   RP_hours)

              !Convert to a 0 - 24 scale instead of the 0 - VenusHoursPerDay             
              LST = 24*LST/RP_hours

              !0 - 90 in 5 deg increments                            
              do iLatTable = 1,19
                 if (ABS(Latitude(iLat,iBlock)*180.0/pi) .ge. 5*(iLatTable - 1) .and. &
                      ABS(Latitude(iLat,iBlock)*180.0/pi) .le. 5*(iLatTable)) then
                    exit
                 endif
              enddo

              do iAltTable = 1,16
                 if (Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0 .ge. 70 + 2*(iAltTable - 1) &
                      .and. &
                      Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0 .le. 70 + 2*(iAltTable)) then
                    exit
                 endif
              enddo

              QnirTOT(iLon,iLat,iAlt,iBlock) = diurnalHeating(iLatTable,iAltTable)* &
                   cos(2*pi/24*(LST - 12)) + &
                   semidiurnalHeating(iLatTable,iAltTable)* &
                   cos(4*pi/24*(LST - 12))
              
              if (QnirTOT(iLon,iLat,iAlt,iBlock) < 0.0) then
                 QnirTOT(iLon,iLat,iAlt,iBlock) = 0.0
              endif

              !Given in K/s. No need for rho or cp because that                              
              !converts J/s -> K/s. Need TempUnit to normalize it to              
              !GITM units                                                       
              QnirTOT(iLon, iLat, iAlt, iBlock) = QnirTOT(iLon, iLat, iAlt, iBlock)  &
                   / TempUnit(iLon,iLat,iAlt)! &      
              ! / Rho(iLon,iLat,iAlt, iBlock) &       
              ! / cp(iLon,iLat,iAlt,iBlock)           
           endif
        enddo
     enddo
  enddo

  call end_timing("IR Heating")
  
end subroutine calc_ir_heating
