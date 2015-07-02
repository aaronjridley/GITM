!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine calc_avesza(iLon,iLat,iBlock, SinDec, CosDec)

  use ModGITM
  use ModTime
  use ModInputs
  use ModEUV

  implicit none

  integer, intent(in) :: iLon, iLat, iBlock
  real, intent(in) :: SinDec, CosDec

  !################### Nelli, April 07 ##Mars only##################
  !variables needed for calculation of the cosine of the zenith angle
  real :: X5,X6a,X6,X3,X4,XT1,XT2,NDT,TOFDAY,COSZ1,COSZ2,TR1,TRR
  !################################################################

  !################ Nelli, April 07 ##Mars only#####################
  !Calculating the cosine of the zenith angle to input into lower
  !atmosphere radiation code

  if (tSimulation == 0.0 .or. Dt == 0.0) then
     AveCosSza(iLon,iLat,iBlock) = max(0.0,cos(sza(iLon,iLat,iBlock)))
  else

     if( floor((tSimulation - dT)/DtLTERadiation) /= &
          floor(tSimulation/DtLTERadiation)) then         

        X5 = 2.0*PI/Rotation_Period

        ! Longitude needs to be multiplied by 180/pi to 
        ! convert from radians to degree

        if(Longitude(iLon,iBlock)*180.0/pi.LE.180.0) then
           X6a=Longitude(iLon,iBlock)/pi + 1.0
        else
           X6a=Longitude(iLon,iBlock)/pi - 1.0
        endif
        X6 = pi*X6a !same as 2*pi*X6a/2
        X3   = SinDec*sin(Latitude(iLat,iBlock))
        X4    = CosDec*cos(Latitude(iLat,iBlock))
        !   TOFDAY = localtime(ilon)-((Longitude(iLon,iBlock)/360.0)&
             !             *HoursPerDay)
        TOFDAY = localtime(ilon)-((Longitude(iLon,iBlock)/(2.0*pi))&
             *HoursPerDay)
        if(TOFDAY.LT.0.0) TOFDAY=TOFDAY + HoursPerDay
        XT1   = (TOFDAY*Rotation_Period/HoursPerDay)-0.5*DT
        NDT   = DtLTERadiation 
        XT2   = XT1 + NDT
        COSZ1 = X3 + X4*COS(X5*XT1+X6)
        COSZ2 = X3 + X4*COS(X5*XT2+X6)

        IF(COSZ1.GE.0.0.AND.COSZ2.GE.0.0) THEN

           !C                 SUN ABOVE THE HORIZON FOR THE ENTIRE PERIOD.

           AveCosSza(iLon,iLat,iBlock) =&
                X3+X4*(SIN(X5*XT2+X6)-SIN(X5*XT1+X6))/&
                (NDT*X5)

        ELSE IF(COSZ1.LE.0.0.AND.COSZ2.LE.0.0) THEN

           !C                 SUN BELOW THE HORIZON FOR THE ENTIRE PERIOD.

           AveCosSza(iLon,iLat,iBlock) = 0.0

        ELSE

           !C                 SUN RISES OR SETS DURING THE TIME INTERVAL.

           TR1 = ACOS(-X3/X4)
           TRR = (TR1-X6)/X5

           !C                 GET 'TRR' THAT IS BETWEEN XT1 AND XT2.

           IF(TRR.GE.XT1.AND.TRR.LE.XT2) GOTO 999
           IF(TRR.GT.XT2) GOTO 998
           TRR = TRR+Rotation_Period
           IF(TRR.GE.XT1.AND.TRR.LE.XT2) GOTO 999
           TRR = TRR+Rotation_Period
           IF(TRR.GE.XT1.AND.TRR.LE.XT2) GOTO 999

998        CONTINUE

           TR1 = 2.0*PI-TR1
           TRR = (TR1-X6)/X5

           IF(TRR.GE.XT1.AND.TRR.LE.XT2) GOTO 999
           IF(TRR.LT.XT1) TRR = TRR+Rotation_Period
           IF(TRR.GT.XT2) TRR = TRR-Rotation_Period

999        CONTINUE

           IF(COSZ1.LT.0.0.AND.COSZ2.GT.0.0) THEN

              !C                    SUN RISES AFTER TIME XT1.

              IF(TRR.LT.XT1.OR.TRR.GT.XT2) THEN
                 TRR = XT2
              END IF

              AveCosSza(iLon,iLat,iBlock) =&
                   (X3*(XT2-TRR)+X4*(SIN(X5*XT2+X6)-&
                   SIN(X5*TRR+X6))/X5)/NDT

           ELSE

              !C                    SUN SETS BEFORE TIME XT2.

              IF(TRR.LT.XT1.OR.TRR.GT.XT2) THEN
                 TRR = XT1
              END IF

              AveCosSza(iLon,iLat,iBlock) =& 
                   (X3*(TRR-XT1)+X4*(SIN(X5*TRR+X6)-&
                   SIN(X5*XT1+X6))/X5)/NDT

           END IF

        END IF

     endif

  endif
  !##########################################################

end subroutine calc_avesza
