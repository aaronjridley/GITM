C Subroutine MAXT
C
C This software is part of the GLOW model.  Use is governed by the Open Source
C Academic Research License Agreement contained in the file glowlicense.txt.
C For more information see the file glow.txt.
C
C Stan Solomon, 11/89, 9/91, 1/94, 3/05
C
C Generates Maxwellian electron spectra with, optionally, a low energy tail
C of the form used by Meier et al., JGR 94, 13541, 1989.
C
C Supplied by calling routine:
C     EFLUX  total energy flux in erg cm-2 s-1
C     EZER   characteristic energy in eV
C     ENER   energy grid in eV
C     DEL    energy bin width in eV
C     NBINS  number of energy bins (dimension of ENER, DEL, and PHI)
C     ITAIL  1 = Maxwellian with low-energy tail, 0 = regular Maxwellian
C     FMONO  additional monoenergetic energy flux in erg cm-2 s-1
C     EMONO  characteristic enerngy of FMONO in eV
C
C Returned by subroutine:
C     PHI    Hemispherical flux in cm-2 s-1 eV-1
C
      SUBROUTINE MAXT (EFLUX, EZER, ENER, DEL, NBINS, ITAIL, 
     >                 FMONO, EMONO, PHI)
      DIMENSION ENER(NBINS), DEL(NBINS), PHI(NBINS)
C
      TE = 0.
C
      IF (EZER .LT. 500.) THEN
        B = 0.8*EZER
      ELSE
        B = 0.1*EZER + 350.
      ENDIF
C
      PHIMAX = EXP(-1.)
C
      DO 300 K=1,NBINS
        ERAT = ENER(K) / EZER
        IF (ERAT .GT. 60.) ERAT = 60.
        PHI(K) = ERAT * EXP(-ERAT)
        IF (ITAIL .GT. 0)
     >    PHI(K) = PHI(K) + 0.4*PHIMAX*(EZER/ENER(K))*EXP(-ENER(K)/B)
        TE = TE + PHI(K) * DEL(K) * ENER(K) * 1.6022E-12
  300 CONTINUE
C
      DO 400 K=1,NBINS
      PHI(K) = PHI(K) * EFLUX / TE
  400 CONTINUE
c
c
c
      if (fmono .gt. 0.) then
      do 500 k=1,nbins
      if (emono.gt.ener(k)-del(k)/2..and.emono.lt.ener(k)+del(k)/2.)
     >  phi(k)=phi(k)+fmono/(1.6022E-12*del(k)*ener(k))
  500 continue
      endif
C
      RETURN
      END
