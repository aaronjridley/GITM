C Subroutine RCOLUM
C
C This software is part of the GLOW model.  Use is governed by the Open Source
C Academic Research License Agreement contained in the file glowlicense.txt.
C For more information see the file glow.txt.
C
C Stan Solomon, 1988, 1991
C
C Calculates the column density ZCOL for each species ZMAJ above height
C ZZ at zenith angle CHI.  Uses a Chapman function fit [Smith and Smith,
C JGR 77, 3592, 1972].  If CHI is less than 90 degrees, column
C densities are calculated directly; if CHI is greater than 90 degrees
C the column density at grazing height for 90 degrees is calculated and
C doubled and the column density above ZZ(J) is subtracted; if CHI is
C such that the grazing height is less than the radius of the earth the
C column densities are set to 'infinity', i.e., 1.0E30.  Densities
C supplied in array ZMAJ are used in the calculation except where
C grazing height is below the lowest level specified, in this case
C values are interpolated logarithmically from a US Standard Atmosphere
C at sea level, the tropopause, the stratopause, and the mesopause.
C
      SUBROUTINE RCOLUM (CHI, ZZ, ZMAJ, TN, ZCOL, ZVCD, JMAX, NMAJ)
C
      PARAMETER (NM=3)
      PARAMETER (NU=4)
C
      DIMENSION ZZ(JMAX), ZMAJ(NMAJ,JMAX), TN(JMAX), ZCOL(NMAJ,JMAX),
     >          ZVCD(NMAJ,JMAX), ZCG(NM), ZUS(NU), TNUS(NU), ZCUS(NM,NU)
C
C    DAVE PAWLOWSKI CHANGED ZUS(4) TO 1.e7 from 9.e6!!!!!!!!!!!!!! 
      DATA PI/3.1415926535/, RE/6.37E8/
      DATA ZUS/0., 1.5E6, 5.E6, 1.E7/, TNUS/288., 217., 271., 187./
      DATA ZCUS/8.00E17, 4.54E24, 1.69E25,
     >          8.00E17, 5.46E23, 2.03E24,
     >          8.00E17, 3.63E21, 1.35E22,
     >          7.80E17, 8.48E18, 3.16E19/
C
C
      CALL VCD (ZZ, ZMAJ, ZVCD, JMAX, NMAJ)
C     
      IF (CHI .GE. 2.) THEN 
         DO 40 I=1,NMAJ
            DO 40 J=1,JMAX
               ZCOL(I,J) = 1.0E30
 40         CONTINUE
            RETURN
      ENDIF
C     
      IF (CHI .LE. PI/2.) THEN
         DO 60 I=1,NMAJ
            DO 60 J=1,JMAX
               ZCOL(I,J) = ZVCD(I,J) * CHAP(CHI,ZZ(J),TN(J),I)
 60         CONTINUE
      ELSE
         DO 220 J=1,JMAX
            GHRG=(RE+ZZ(J))*SIN(CHI) 
            GHZ=GHRG-RE 
C            write(*,*) "GHRG: ", GHRG,ZZ(J),CHI  
            IF (GHZ .LE. 0.) THEN
               DO 80 I=1,NMAJ
                  ZCOL(I,J) = 1.0E30
 80            CONTINUE
               GOTO 220
            ENDIF

            IF (GHZ .GE. ZZ(1)) THEN
               DO 100 JG=1,J-1
                  IF (ZZ(JG) .LE. GHZ .AND. ZZ(JG+1) .GT. GHZ) GOTO 120
 100           CONTINUE
 120           TNG = TN(JG)+(TN(JG+1)-TN(JG))*(GHZ-ZZ(JG))/(ZZ(JG+1)-ZZ(JG))
               DO 140 I=1,NMAJ
                  ZCG(I) = ZVCD(I,JG) * (ZVCD(I,JG+1) / ZVCD(I,JG)) **
     >                 ((GHZ-ZZ(JG)) / (ZZ(JG+1)-ZZ(JG)))
 140           CONTINUE
            ELSE
               DO 160 JG=1,3
                  IF (ZUS(JG) .LT. GHZ .AND. ZUS(JG+1) .GT. GHZ) GOTO 180
C                  IF (J.eq.1) write(*,*) JG, ZUS(JG), GHZ, ZUS(JG+1)
 160           CONTINUE
 180           TNG = TNUS(JG)
     >              + (TNUS(JG+1)-TNUS(JG))*(GHZ-ZUS(JG))/(ZUS(JG+1)-ZUS(JG))
               DO 200 I=1,NMAJ
                  ZCG(I) = ZCUS(I,JG) * (ZCUS(I,JG+1) / ZCUS(I,JG)) **
     >                 ((GHZ-ZUS(JG)) / (ZUS(JG+1)-ZUS(JG)))
C                  if (J.eq.JMAX-1) write(*,*) "ZCG: ", ZCG(I), ZCUS(I,JG),GHZ
C                  if (J.eq.JMAX-1) write(*,*) "-------> ",ZUS(JG)
 200           CONTINUE 
c               if (J.eq.1)  then
c                     write(*,*) 'ZCG:   ', ZCG, ZCUS(:,JG), JG
c                     write(*,*) 'ZCUSL  ', ZCUS(:,JG+1)
c                     write(*,*) 'ZUS:   ',ZUS(JG), ZUS(JG+1)
c                  endif
            ENDIF

            DO 210 I=1,NMAJ
               ZCOL(I,J) = 2. * ZCG(I) * CHAP(PI/2.,GHZ,TNG,I)
     >              - ZVCD(I,J) * CHAP(CHI,ZZ(J),TN(J),I)
               
C     if (J.eq.1 .and. I.eq.3) write(*,*) ZCOL(I,J), ZCG(I),CHAP(PI/2.,GHZ,TNG,I)
C     if (J.eq.1 .and. I.eq.3) write(*,*) ZVCD(I,J),CHAP(CHI,ZZ(J),TN(J),I)

 210        CONTINUE
 220     CONTINUE
               
       
 
      ENDIF
C
      RETURN 
      END 

C
C
C
C
      FUNCTION CHAP (CHI, Z, T, I)
      PARAMETER (NMAJ=3)
      DIMENSION AM(NMAJ)
      DATA AM/16., 32., 28./, PI/3.1415926535/, RE/6.37E8/, G/978.1/
      GR=G*(RE/(RE+Z))**2 
      HN=1.38E-16*T/(AM(I)*1.662E-24*GR)
      HG=(RE+Z)/HN 
      HF=0.5*HG*(COS(CHI)**2) 
      SQHF=SQRT(HF) 
      CHAP=SQRT(0.5*PI*HG)*SPERFC(SQHF) 
      RETURN
      END
C
C
C
C
      FUNCTION SPERFC(DUMMY) 
      IF (DUMMY .LE. 8.) THEN
        SPERFC = (1.0606963+0.55643831*DUMMY) /
     >           (1.0619896+1.7245609*DUMMY+DUMMY*DUMMY)
      ELSE
        SPERFC=0.56498823/(0.06651874+DUMMY) 
      ENDIF 
      RETURN 
      END 
C
C
C
C
      SUBROUTINE VCD(ZZ,ZMAJ,ZVCD,JMAX,NMAJ)
      DIMENSION ZZ(JMAX), ZMAJ(NMAJ,JMAX), ZVCD(NMAJ,JMAX)
C
      DO 200 I=1,NMAJ
      ZVCD(I,JMAX) =   ZMAJ(I,JMAX)
     >               * (ZZ(JMAX)-ZZ(JMAX-1))
     >               / ALOG(ZMAJ(I,JMAX-1)/ZMAJ(I,JMAX))
      DO 200 J=JMAX-1,1,-1
         RAT = ZMAJ(I,J+1) / ZMAJ(I,J)
         ZVCD(I,J) =   ZVCD(I,J+1)
     >        + ZMAJ(I,J) * (ZZ(J)-ZZ(J+1)) / ALOG(RAT) * (1.-RAT)
  200 CONTINUE
C        write(*,*) ZVCD(3,:)
C     stop
      RETURN
      END
