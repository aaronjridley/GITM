C Subroutine QBACK
C
C This software is part of the GLOW model.  Use is governed by the Open Source
C Academic Research License Agreement contained in the file glowlicense.txt.
C For more information see the file glow.txt.
C
C Stan Solomon, 11/88, 11/92
C Comment updated 3/05
C
C Estimates background ("nighttime") ionization rates.
C Intended to conform to the TIGCM background (Roble et al., 1989).
C Three fluxes, from the stellar background and multiple scattering
C of solar atomic hydrogen emissions in the geocorona, are considered:
C The sum of ionizing flux in the 100-1000 A range (FIONT)
C H Lyman beta (FLYBT)
C H Ly alpha (FLYAT)
C
      SUBROUTINE QBACK (ZMAJ, ZNO, ZVCD, PHOTOI, PHONO, JM, NMAJ, NST)
C
      DIMENSION ZMAJ(NMAJ,JM), ZNO(JM), ZVCD(NMAJ,JM),
     +          PHOTOI(NST,NMAJ,JM), PHONO(NST,JM)
C
      DATA FIONT/5.0E7/, FLYBT/1.0E7/, FLYAT/1.0E9/
      DATA SIGIO/1.0E-17/, SIGIO2/2.0E-17/, SIGIN2/2.0E-17/,
     +     SLBAO2/1.6E-18/, SLBIO2/1.0E-18/, SLAAO2/1.0E-20/,
     +     SLAINO/2.0E-18/
C
C Calculate ionization rates at each altitude:
C
      DO 200 J=1,JM
        TAUI = 2.*(SIGIO*ZVCD(1,J)+SIGIO2*ZVCD(2,J)+SIGIN2*ZVCD(3,J))
        IF (TAUI .GT. 60.) TAUI = 60.
        TAULYB = 2.*SLBAO2*ZVCD(2,J)
        IF (TAULYB .GT. 60.) TAULYB = 60.
        TAULYA = 2.*SLAAO2*ZVCD(2,J)
        IF (TAULYA .GT. 60.) TAULYB = 60.
        FION = FIONT * EXP(-TAUI)
        PHOTOI(1,1,J) = PHOTOI(1,1,J) + FION * ZMAJ(1,J) * SIGIO
        PHOTOI(1,2,J) = PHOTOI(1,2,J) + FION * ZMAJ(2,J) * SIGIO2
     +                  + FLYBT * EXP(-TAULYB) * ZMAJ(2,J) * SLBIO2
        PHOTOI(1,3,J) = PHOTOI(1,3,J) + FION * ZMAJ(3,J) * SIGIN2
        PHONO(1,J) = PHONO(1,J) + FLYAT * EXP(-TAULYA) * ZNO(J) * SLAINO
  200 CONTINUE
C
      RETURN
      END
