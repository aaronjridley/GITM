C Subroutine GCHEM
C
C This software is part of the GLOW model.  Use is governed by the Open Source
C Academic Research License Agreement contained in the file glowlicense.txt.
C For more information see the file glow.txt.
C
C Stan Solomon, 1988, 1989, 1992, 1999, 2005
C
C Includes quartic solution to electron density equation.
C Electron density must be supplied above 200 km in array ZE, a priori
C values are expected but not strictly required at and below 200 km,
C depending on the value of KCHEM.  An initial guess of the N(2D)
C density in array ZND is also expected but not strictly required.
C Other neutral species (O, N2, O2, NO, N(4S)) must be supplied,
C see subroutine GLOW.
C
C Chemical calculations are controlled by switch KCHEM:
C 0 = no calculations at all are performed.
C 1 = electron density, O+(4S), N+, N2+, O2+, NO+ supplied at all
C     altitudes; O+(2P), O+(2D), excited neutrals, and emission rates
C     are calculated.
C 2 = electron density, O+(4S), O2+, NO+ supplied at all altitudes;
C     O+(2P), O+(2D), N+, N2+, excited neutrals, emissions calculated.
C 3 = electron density supplied at all altitudes; everthing else
C     calculated.
C 4 = electron density supplied above 200 km; electron density below
C     200 km is calculated, everything else calculated at all altitudes.
C     An initial guess of the electron density below 200 km should be
C     supplied in array ZE.  Electron density for the next two levels
C     above J200 is log interpolated between E(J200) and E(J200+3).
C
C For definitions of common block /CGLOW/ see subroutine GLOW
C
C Other definitions:
C A        Einstein coefficients; s-1
C B        Branching ratios
C BZ       Altitude-dependent branching ratios
C G        Resonant scattering g-factors at each altitude; s-1
C KZ       Temperature dependent rate coeffs at each altitude; cm3s-1
C OEI      O electron impact ionization rates; cm-3s-1
C O2EI     O2   "       "        "       "      "
C RN2EI    N2   "       "        "       "      "
C O2PI     O2 photoionization rates; cm-3s-1
C RN2PI    N2       "          "      "
C RN2ED    N2 electron impact dissociation rate; cm-3s-1
C SRCED    O2     "      "         "         " (SR continuum); cm-3s-1
C P        Volume production rate for each species, altitude; cm-3s-1
C L        Loss rate for each species, altitude; s-1
C T1       Effective temperature divided by 300 for O+ + N2; K
C T2          "          "        "                 O+ + O2; K
C T3          "          "        "                 N2+ + O; K
C T4          "          "        "                 N2+ + O2; K
C T5          "          "        "                 O+ + NO; K
C QQ, RR, SS, TT, UU, VV, WW, XX:  Combined terms for calculation of
C                              O+(4S) given e
C AA, BB, CC, DD, EE, FF, GG, HH: Combined terms for solution of
C                        electron density equation, Roble & Ridley, 1988
C COEF     Coefficients of quartic electron density equation, see "
C ROOT     Roots of          "       "         "       "
C
C
C References for rate coefficients, transition coefficients, branching
C ratios, and g-factors:
C
C      k1  O+(4S) + N2             St. Maurice & Torr, 1978
C                                  (from Albritton et al., 1977)
C      k2  O+(4S) + O2             Combination of (Chen et al, 1978
C                                  and St. Maurice & Torr, 1978)
C      k3  N2+ + O -> NO+ + O      New fit to McFarland et al, 1974
C      k4  N2+ + O2                McFarland et al, 1973
C      k5  N(2D) + O2              Lin & Kaufman, 1971,
C                                  cf. Piper et al, 1987
C      k6  N(2D) + O               Fell et al., 1990
C      k7  N(2D) + e               Frederick & Rusch, 1977;
C                                  Queffelec et al, 1985
C      k8  O(1D) + N2              Streit et al, 1976
C      k9  O(1D) + O2              Streit et al, 1976
C      k10 O(1D) + e               Link, 1982
C      k11 O(1S) + O               Slanger & Black, 1981
C      k12 O+(2D) + N2             Johnsen & Biondi, 1980
C      k13 O+(2D) + O2             Johnsen & Biondi, 1980
C      k14 O+(2D) + e              Henry et al., 1969
C      k15 O+(2D) + O              Torr & Torr, 1980
C      k16 O+(2P) + N2             Rusch et al, 1977
C      k17 O+(2P) + O2             Link, 1982
C      k18 O+(2P) + e              Henry et al., 1969
C      k19 O+(2P) + O              Rusch et al, 1977
C      k20 O2(c) + O               Solheim & Llewellyn, 1979
C      k21 O2(c) + N2              Solheim & Llewellyn, 1979
C      k22 NO+ + e                 Walls & Dunn, 1974; Torr et al, 1977;
C                                  Alge et al, 1983; Dulaney et al,1987;
C                                  Davidson & Hobson, 1987
C      k23 N2+ + e                 Mehr & Biondi, 1969
C      k24 O2+ + e                 Mehr & Biondi; Walls & Dunn, 1974;
C                                  Torr et al, 1976; Alge et al, 1983
C      k25 N+ + O2                 Langford et al, 1985
C      k26 N2(A) + O               Piper et al, 1981b (av. v=1,2)
C      k27 O(1D) + O               Abreu et al, 1986; Yee, pc, 1991
C      k28 O + et                  Link, 1982
C      k29 N2(A) + O2              Piper et al, 1981a
C      k30 O2+ + NO                Lindeger & Ferguson, 1974; G&R, 1979
C      k31 N(2D) + NO              Black et al, 1969; fr. Roble, 1986
C      k32 N+ + O                  Torr, 1985 (Levine, Photochemistry)
C      k33 N(2P) + O               Zipf et al, 1980; Young & Dunn, 1975
C      k34 N(2P) + O2              Zipf et al, 1980; Rawlins, 1988
C                                  (cf Ianuzzi & Kaufman, 1980, 3.5E-12)
C      k35 N(2P) + NO              Rees & Jones, 1973, Zipf et al, 1980
C      k36 O(1S) + O2              Slanger et al, 1972; fr. Bates, 1978
C      k37 O2+ + N                 Fehsenfeld (1977)
C      k38 O+ + N(2D)              Bates, 1989 (PSS 37, 363)
C      k39 N2+ + O -> N2 + O+      Torr, 1985; Torr et al, 1988;
C                                  Knutsen et al, 1988
C      k40 O+ + NO                 St. Maurice & Torr, 1978
C      A1  5200       N(4S-2D)     Wiese et al, 1966
C      A2  6300       O(3P-1D)     Baluja and Zeippen, 1988
C      A3  6364       O(3P-1D)     Baluja and Zeippen, 1988
C      A4  2972       O(3P-1S)     Kernahan & Pang, 1975
C      A5  5577       O(1D-1S)     Kernahan & Pang, 1975
C      A6  3726       O+(4S-2D)    Kernahan & Pang, 1975
C      A7  2470       O+(4S-2P)    Weise et al, 1966
C      A8  7319-30    O+(2D-2P)    Weise et al, 1966
C      A9  (Hertz II) O2(X-c)      Solheim & Llewellyn, 1978
C      A10 (Veg-Kap)  N2(X-A)      Shemansky, 1969
C      A11 3466       N(4S-2P)     Chamberlain, 1961
C      A12 10400      N(2D-2P)     Chamberlain, 1961
C      B1  O(1S) from O2+ + e      Yee et al, 1988
C      B2  O(1D) from O2+ + e      Abreu et al, 1986
C      B3  N(2D) from NO+ + e      Kley et al, 1976
C      B4  N(2D) from N2+ + e      Queffelec et al, 1985
C      B5  N(2D) from N2+ + O      Frederick & Rusch, 1977
C      B6  O(1D) from N(2D) + O2   Link, 1983; Langford et al, 1985
C      B7  O(1D) from O+(2D) + O   ?
C      B8  O+(2D) from O+(2P) + e  Link, 1982
C      B9  O+(2P) from O + e*      Gerard & Rusch, 1979; Jones, 1975
C      B10 O+(2D) from O + e*      Gerard & Rusch, 1979; Jones, 1975
C      B11 O+(4S) from O + e*      Gerard & Rusch, 1979; Jones, 1975
C      B12 O+(2P) from O2 + e*     Link, 1982; guess of .3/3
C      B13 O+(2D) from O2 + e*     Link, 1982; guess of .3/3
C      B14 O+(4S) from O2 + e*     Link, 1982; guess of .3/3
C      B15 N+ from N2 + e*         Richards & Torr, 1985
C      B16 N(2D) from above        Zipf et al, 1980
C      B17 O(1D) from N+ + O2      Langford et al, 1985
C      B18 O(1S) from N2(A) + O    Sharp & Torr, 1979
C      B19 O(1S) from O2(*) + O    ? (= 0 at present)
C      B20 O(1D) from N(2D) + O    ?
C      B21 NO+ from N+ + O2        Langford et al, 1985
C      B22 O2+ from N+ + O2        Langford et al, 1985
C      B23 N(2P) from N2+ + e      Queffelec et al, 1985
C      B24 N2 + e* pri. dis/ion    xsect ratio > 250 eV
C      B25 N(2D) from N2 + e* dis  Zipf et al, 1980
C      B26 N(2P) from N2 + e* dis  Zipf et al, 1980
C      B27 N(2D) from N2 + hv      Richards et al, 1981 (add to B28)
C      B28 N(2P) from N2 + hv      ? (cf Zipf & McGlaughlin, 1978)
C      B29 N(2D) from N(2P) + O    ?
C      B30 O+(2P) from O2 + hv     ?
C      B31 O+(2D)  "               ?
C      B32 O+(4S)  "               ?
C      B33 O(1S) from O2+ + N      Frederick et al, 1976; Kopp ea, 1977
C      B34 O(1S) from N(2D) + NO   Frederick et al, 1976; Kopp ea, 1977
C      B35 O2 + e* pri SRC/ion     xsect ratio > 250 eV
C      B36 N2+(B) from N2 + e*     Borst & Zipf, 1970;
C                                  Shemansky & Broadfoot, 1971
C      B37 (0,0) (3914) fr. N2+(B) Shemansky & Broadfoot, 1971
C      B38 (0,1) (4278) fr. N2+(B) Shemansky & Broadfoot, 1971
C      B39 (0,0) (3371) fr. N2(C)  Conway, 1983; Benesch et al, 1966
C      B40 (0,9) (3352) fr. N2(A)  Cartwright, 1978; Shemansky, 1969
C      B41 O+(2Po) fr. O+(2Pe)     Kirby et al, 1979
C      B42 O+(2Do) fr. O+(2Pe)     Kirby et al, 1979
C      B43 N2(C) bound fraction    ?
C      B44 7990 fr. O(3s'3D)       appx. fr. Hecht, p.c.
C      G1  N2+B(0,0) (3914)        Broadfoot, 1967
C      G2  N2+B(0,1) (4278)        Broadfoot, 1967
C
C
C Array dimensions:
C JMAX    number of altitude levels
C NBINS   number of energetic electron energy bins
C LMAX    number of wavelength intervals for solar flux
C NMAJ    number of major species
C NEX     number of ionized/excited species
C NW      number of airglow emission wavelengths
C NC      number of component production terms for each emission
C NST     number of states produced by photoionization/dissociation
C NEI     number of states produced by electron impact
C NR      number of rate coefficients, branching ratios, A and G factors
C NF      number of available types of auroral fluxes
C
C
      SUBROUTINE GCHEM
C
      use Mod_GLOW

      PARAMETER (NR=50)
      REAL KZ, L
C
      DIMENSION A(NR), B(NR), BZ(NR,JMAX), G(NR,JMAX), KZ(NR,JMAX),
     >          OEI(JMAX), O2EI(JMAX), RN2EI(JMAX), O2PI(JMAX),
     >          RN2PI(JMAX), RN2ED(JMAX), SRCED(JMAX),
     >          P(NEX,JMAX), L(NEX,JMAX),
     >          T1(JMAX), T2(JMAX), T3(JMAX), T4(JMAX), T5(JMAX),
     >          QQ(JMAX), RR(JMAX), SS(JMAX), TT(JMAX), UU(JMAX),
     >          VV(JMAX), WW(JMAX), XX(JMAX),
     >          AA(JMAX), BB(JMAX), CC(JMAX), DD(JMAX), EE(JMAX),
     >          FF(JMAX), GG(JMAX), HH(JMAX)
C     
     
      DOUBLE PRECISION COEF(JMAX,5), ROOT(JMAX)
C
      DATA RE/6.37E8/
      DATA A / 1.07E-5, 0.00585, 0.00185, 0.04500, 1.06000,
     >         9.70E-5, 0.04790, 0.17120, 0.00100, 0.77000,
     >         0.00540, 0.07900, 38*0.0 /
      DATA B/0.07, 1.20, 0.76, 1.85, 1.00, 0.10, 0.50, 0.81, 0.20, 0.32,
     >       0.48, 0.10, 0.10, 0.10, 0.16, 0.50, 0.30, 0.19, 0.00, 0.10,
     >       0.43, 0.51, 0.10, 0.60, 0.54, 0.44, 0.80, 0.20, 1.00, 0.33,
     >       0.33, 0.34, 0.21, 0.20, 0.10, 0.11, 0.65, 0.20, 0.24, 0.02,
     >       0.18, 0.72, 0.75, 0.10, 6*0./
C
C
      
      IF (KCHEM .EQ. 0) RETURN
C
C
C Zero airglow and density arrays:
C
      DO 10 I=1,JMAX
      DO 10 IW=1,NW
      ZETA(IW,I) = 0.0
   10 CONTINUE
C
      DO 20 I=1,JMAX
      DO 20 IW=1,NW
      DO 20 IC=1,NC
      ZCETA(IC,IW,I) = 0.0
   20 CONTINUE
C
      IF (KCHEM .GE. 3) THEN
      DO 30 I=1,JMAX
      DO 30 IX=1,NEX
      ZXDEN(IX,I)=0.
   30 CONTINUE
      ENDIF

C
C
C Assign g-factors at altitudes which are sunlit:
C
      DO 40 I=1,JMAX
      DO 40 N=1,NR
      G(N,I) = 0.0
   40 CONTINUE
C
      DO 50 I=1,JMAX
      GH = (RE+ZZ(I)) * SIN(SZA)
      IF (SZA .LT. 1.6 .OR. GH .GT. RE) THEN
        G(1,I) = 0.041
        G(2,I) = 0.013
      ENDIF
   50 CONTINUE


C
C
C Calculate rate coefficients as a function of altitude:
C

      DO 70 I=1,JMAX
      T1(I) = (16.*ZTN(I)+28.*ZTI(I)) / (16.+28.) / 300.
      T2(I) = (16.*ZTN(I)+32.*ZTI(I)) / (16.+32.) / 300.
      T3(I) = (28.*ZTN(I)+16.*ZTI(I)) / (28.+16.) / 300.
      T4(I) = (28.*ZTN(I)+32.*ZTI(I)) / (28.+32.) / 300.
      T5(I) = (16.*ZTN(I)+30.*ZTI(I)) / (16.+30.) / 300.
      IF (T1(I) .LT. 5.6667) THEN
        KZ(1,I) = 1.533E-12 - 5.92E-13*T1(I) + 8.6E-14*T1(I)**2
      ELSE
        KZ(1,I) = 2.73E-12 - 1.155E-12*T1(I) + 1.483E-13*T1(I)**2
      ENDIF
      IF (T2(I) .LT. 6.6667) THEN
        KZ(2,I) = 3.53E-11 - 1.84E-11*T2(I) + 4.62E-12*T2(I)**2
     >                     - 4.95E-13*T2(I)**3 + 2.00E-14*T2(I)**4
      ELSE
        KZ(2,I) = 2.82E-11 - 7.74E-12*T2(I) + 1.073E-12*T2(I)**2
     >                     - 5.17E-14*T2(I)**3 + 9.65E-16*T2(I)**4
      ENDIF
      KZ(3,I) = 1.793e-10 - 6.242e-11*T3(I) + 1.225e-11*T3(I)**2
     >                    - 1.016e-12*T3(I)**3 + 3.087e-14*T3(I)**4
      KZ(4,I) = 5.0E-11 * (1./T4(I)) ** 0.8
      KZ(5,I) = 6.0E-12
      KZ(6,I) = 6.9E-13
      KZ(7,I) = 5.5E-10 * (ZTE(I)/300.) ** 0.5
      KZ(8,I) = 2.0E-11 * EXP(107.8/ZTN(I))
      KZ(9,I) = 2.9E-11 * EXP(67.5 /ZTN(I))
      KZ(10,I) = 8.1E-10 * (ZTE(I)/300.) ** 0.5
      KZ(11,I) = 2.0E-14
      KZ(12,I) = 8.0E-10
      KZ(13,I) = 7.0E-10
      KZ(14,I) = 6.6E-08 * (300./ZTE(I)) ** 0.5
      KZ(15,I) = 1.0E-11
      KZ(16,I) = 4.8E-10
      KZ(17,I) = 4.8E-10
      KZ(18,I) = 1.7E-07 * (300./ZTE(I)) ** 0.5
      KZ(19,I) = 5.2E-11
      KZ(20,I) = 2.1E-11 * EXP(-1136./ZTN(I))
      KZ(21,I) = 1.0E-13
      KZ(22,I) = 4.2E-07 * (300./ZTE(I)) ** 0.85
      KZ(23,I) = 1.8E-07 * (300./ZTE(I)) ** 0.39
      KZ(24,I) = 1.95E-07 * (300./ZTE(I)) ** 0.70
      IF (ZTE(I) .GE. 1200.) KZ(24,I) = 1.6E-07 * (300./ZTE(I)) ** 0.55
      KZ(25,I) = 6.0E-10
      KZ(26,I) = 3.1E-11
      KZ(27,I) = 3.0E-12
      IF (ZTE(I) .LT. 500.) THEN
        KZ(28,I) = 1.0E-29
      ELSE
        KZ(28,I) = 2.6E-11 * ZTE(I)**0.5 * EXP(-22740./ZTE(I))
      ENDIF
      KZ(29,I) = 4.1E-12
      KZ(30,I) = 4.4E-10
      KZ(31,I) = 7.0E-11
      KZ(32,I) = 1.0E-12
      KZ(33,I) = 1.2E-11
      KZ(34,I) = 2.0E-12
      KZ(35,I) = 1.8E-10
      KZ(36,I) = 4.0E-12 * EXP(-865./ZTN(I))
      KZ(37,I) = 1.2E-10
      KZ(38,I) = 1.3E-10
      KZ(39,I) = 2.0E-11
      IF (T5(I) .LT. 5) THEN
        KZ(40,I) = 8.36E-13 - 2.02E-13*T5(I) + 6.95E-14*T5(I)**2
      ELSE
        KZ(40,I) = 5.33E-13 - 1.64E-14*T5(I) + 4.72E-14*T5(I)**2
     >                      - 7.05E-16*T5(I)**3
      ENDIF
   70 CONTINUE
C
C
C Calculate Electron impact ionization, photoionization, and electron
C impact dissociation rates at each altitude; put a priori electron
C density in calculated electron density array: put rough estimate of
C O+ and a priori N(2D) in DEN array:
      
C
      DO 100 I=1,JMAX
      OEI(I)   = SION(1,I) + PIA(1,I)
      O2EI(I)  = SION(2,I) + PIA(2,I)
      RN2EI(I) = SION(3,I) + PIA(3,I)
      O2PI(I)  = PHOTOI(1,2,I) + PHOTOI(2,2,I) + PHOTOI(3,2,I)
      RN2PI(I) = PHOTOI(1,3,I) + PHOTOI(2,3,I) + PHOTOI(3,3,I) +
     >           PHOTOI(4,3,I) + PHOTOI(5,3,I)
      RN2ED(I) = AGLW(5,3,I) + AGLW(6,3,I) + AGLW(7,3,I)
     >          + B(24)*PIA(3,I)
      SRCED(I) = AGLW(4,2,I) + B(35)*PIA(2,I)
      ECALC(I)     = ZE(I)
      ZXDEN(10,I)= ZND(I)
  100 CONTINUE
C
C
C Find level below which electron density will be calculated:
C


      IF (KCHEM .GE. 4) THEN
        DO 120 I=JMAX,1,-1
           IF (ZZ(I) .GT. 2.0001E7) J200=I-1
  120   CONTINUE
      ELSE
        J200=0
      ENDIF
      
C
C
C Iterative loop assures that feedback reactions (O+(2P,2D)+e,
C O+(4S)+N(2D), N2++O) are correctly computed:
C

      DO 350 ITER=1,3
C
C
C Calculate atomic ion densities at each altitude:
C
      DO 150 I=1,JMAX
C
C
C O+(2P):
C
      P(1,I)= PHOTOI(3,1,I)
     >      + B(41) * PHOTOI(5,1,I)
     >      + B(30) * PHOTOI(4,2,I)
     >      + B(9)  * OEI(I)
     >      + B(12) * O2EI(I)
      L(1,I)= KZ(16,I) * ZN2(I)
     >      + KZ(17,I) * ZO2(I)
     >      + KZ(19,I) * ZO(I)
     >      + KZ(18,I) * ECALC(I) 
     >      + A(8)
     >      + A(7)
      ZXDEN(1,I) = P(1,I) / L(1,I)
C
C
C O+(2D):
C
      P(2,I)= PHOTOI(2,1,I)
     >      + B(42) * PHOTOI(5,1,I)
     >      + B(31) * PHOTOI(4,2,I)
     >      + B(10) * OEI(I)
     >      + B(13) * O2EI(I)
     >      + B(8)  * KZ(18,I) * ZXDEN(1,I) * ECALC(I)
     >      + A(8)  * ZXDEN(1,I)
      L(2,I)= KZ(12,I) * ZN2(I) 
     >      + KZ(13,I) * ZO2(I) 
     >      + KZ(15,I) * ZO(I)
     >      + KZ(14,I) * ECALC(I)
     >      + A(6)
      ZXDEN(2,I) = P(2,I) / L(2,I)
C
C
C N+:
C
      
      IF (KCHEM .GE. 2) THEN
        P(4,I) = PHOTOI(6,3,I)
     >         + B(15) * RN2EI(I)
     >         + KZ(38,I) * ZXDEN(3,I) * ZXDEN(10,I)
        L(4,I) = KZ(25,I) * ZO2(I)
     >         + KZ(32,I) * ZO(I)
        ZXDEN(4,I) = P(4,I) / L(4,I)
      ENDIF
C
C
C O+(4S):
C
      IF (KCHEM .GE. 3) THEN
        P(3,I)= PHOTOI(1,1,I) + PHOTOI(4,1,I)
     >        + B(32) * PHOTOI(4,2,I)
     >        + B(11) * OEI(I)
     >        + B(14) * O2EI(I) 
     >        + KZ(14,I) * ZXDEN(2,I) * ECALC(I) 
     >        + KZ(15,I) * ZXDEN(2,I) * ZO(I) 
     >        + A(6) * ZXDEN(2,I)
     >        + (1.-B(8)) * KZ(18,I) * ZXDEN(1,I) * ECALC(I)
     >        + KZ(19,I) * ZXDEN(1,I) * ZO(I) 
     >        + A(7) * ZXDEN(1,I)
     >        + KZ(32,I) * ZXDEN(4,I) * ZO(I)
     >        + KZ(39,I) * ZXDEN(5,I) * ZO(I)
        L(3,I)= KZ(1,I) * ZN2(I)
     >        + KZ(2,I) * ZO2(I)
     >        + KZ(38,I) * ZXDEN(10,I)
        ZXDEN(3,I) = P(3,I) / L(3,I)
      ENDIF

C
  150 CONTINUE
      
C
C
C Above 200 km, (or at all altitudes if KCHEM=3) use a priori
C electron density to calculate O+(4S):
C
      IF (KCHEM .GE. 3) THEN
C
      DO 170 I=J200+1,JMAX
      P(5,I)= RN2PI(I)
     >      + (1.-B(15)) * RN2EI(I)
     >      + KZ(12,I) * ZXDEN(2,I) * ZN2(I)
     >      + KZ(16,I) * ZXDEN(1,I) * ZN2(I)
      L(5,I)= KZ(3,I)  * ZO(I)
     >      + KZ(4,I)  * ZO2(I)
     >      + KZ(23,I) * ECALC(I)
     >      + KZ(39,I) * ZO(I)
      ZXDEN(5,I) = P(5,I) / L(5,I)
      QQ(I) = PHONO(1,I)
     >      + KZ(3,I)  * ZXDEN(5,I) * ZO(I)
     >      + B(21) * KZ(25,I) * ZXDEN(4,I) * ZO2(I)
      RR(I) = KZ(30,I) * ZNO(I)
     >      + KZ(37,I) * ZNS(I)
      SS(I) = KZ(1,I) * ZN2(I)
     >      + KZ(40,I) * ZNO(I)
      TT(I) = KZ(22,I) * ECALC(I)
      UU(I) = O2PI(I)
     >      + (1.-B(12)-B(13)-B(14)) * O2EI(I)
     >      + KZ(13,I) * ZXDEN(2,I) * ZO2(I)
     >      + KZ(17,I) * ZXDEN(1,I) * ZO2(I)
     >      + KZ(4,I)  * ZXDEN(5,I) * ZO2(I)
     >      + B(22) * KZ(25,I) * ZXDEN(4,I) * ZO2(I)
      VV(I) = KZ(2,I) * ZO2(I)
      WW(I) = KZ(24,I) * ECALC(I)
     >      + KZ(30,I) * ZNO(I)
     >      + KZ(37,I) * ZNS(I)
      XX(I) = ZXDEN(1,I) + ZXDEN(2,I) + ZXDEN(4,I) + ZXDEN(5,I)
      ZXDEN(3,I) = (TT(I)*WW(I)*ECALC(I) - TT(I)*WW(I)*XX(I) - TT(I)*UU(I)
     >            - QQ(I)*WW(I) - RR(I)*UU(I) ) /
     >           (TT(I)*WW(I) + TT(I)*VV(I) + RR(I)*VV(I) + SS(I)*WW(I))
  170 CONTINUE
C
      ENDIF
      
C
C
C If KCHEM=4, calculate electron density using quartic equation method
C below 200 km:
C
      
      IF (KCHEM .GE. 4) THEN
C
      DO 200 I=1,J200
      AA(I) =  PHONO(1,I)
     >       + KZ(1,I)  * ZXDEN(3,I) * ZN2(I)
     >       + KZ(40,I) * ZXDEN(3,I) * ZNO(I)
     >       + B(21) * KZ(25,I) * ZXDEN(4,I) * ZO2(I)
      BB(I) =  O2PI(I)
     >       + (1.-B(12)-B(13)-B(14)) * O2EI(I)
     >       + KZ(2,I)  * ZXDEN(3,I) * ZO2(I)
     >       + KZ(13,I) * ZXDEN(2,I) * ZO2(I)
     >       + KZ(17,I) * ZXDEN(1,I) * ZO2(I)
     >       + B(22) * KZ(25,I) * ZXDEN(4,I) * ZO2(I)
      CC(I) =  KZ(30,I) * ZNO(I)
     >       + KZ(37,I) * ZNS(I)
      DD(I) =  RN2PI(I)
     >       + (1.-B(15)) * RN2EI(I)
     >       + KZ(12,I) * ZXDEN(2,I) * ZN2(I)
     >       + KZ(16,I) * ZXDEN(1,I) * ZN2(I)
      EE(I) = KZ(3,I) * ZO(I)
      FF(I) = ZXDEN(1,I) + ZXDEN(2,I) + ZXDEN(3,I) + ZXDEN(4,I)
      HH(I) = KZ(4,I)  * ZO2(I)
      GG(I) = KZ(39,I) * ZO(I)
      COEF(I,5) = KZ(22,I) * KZ(23,I) * KZ(24,I)
      COEF(I,4) =(KZ(22,I) * (KZ(24,I)*(EE(I)+HH(I)+GG(I))
     >                        + KZ(23,I)*CC(I))
     >            - KZ(22,I) * KZ(23,I) * KZ(24,I) * FF(I)) / 4.0
      COEF(I,3) =(KZ(22,I) * CC(I) * (EE(I)+HH(I)+GG(I))
     >            - KZ(22,I) * FF(I)
     >                       * (KZ(24,I)*(EE(I)+HH(I)+GG(I))
     >                          + KZ(23,I)*CC(I))
     >            - KZ(24,I) * KZ(23,I) * AA(I)
     >            - KZ(22,I) * KZ(23,I) * BB(I)
     >            - KZ(22,I) * KZ(24,I) * DD(I)) / 6.0
      COEF(I,2) =(- KZ(22,I) * (CC(I)*FF(I)*(EE(I)+HH(I)+GG(I))
     >                         + DD(I)*HH(I) + BB(I)*(EE(I)+HH(I)+GG(I))
     >                          + CC(I)*DD(I))
     >            - KZ(24,I) * (AA(I)*(EE(I)+HH(I)+GG(I)) + DD(I)*EE(I))
     >            - KZ(23,I) * CC(I) * (AA(I)+BB(I))) / 4.0
      COEF(I,1) = - CC(I) * (EE(I)+HH(I)+GG(I)) * (AA(I)+BB(I)+DD(I))
      COEF(I,5) = COEF(I,5) * 1.D10
      COEF(I,4) = COEF(I,4) * 1.D10
      COEF(I,3) = COEF(I,3) * 1.D10
      COEF(I,2) = COEF(I,2) * 1.D10
      COEF(I,1) = COEF(I,1) * 1.D10
  200 CONTINUE
C
      
      CALL VQUART (COEF, ROOT, J200)
C
      
      DO 250 I=1,J200
      ECALC(I) = ROOT(I)
  250 CONTINUE
C
      ECALC(J200+1) = ECALC(J200) * ( ECALC(J200+3) / ECALC(J200) )
     >            ** ( (ZZ(J200+1)-ZZ(J200)) / (ZZ(J200+3)-ZZ(J200)) )
      ECALC(J200+2) = ECALC(J200) * (ECALC(J200+3)/ECALC(J200))
     >            ** ( (ZZ(J200+2)-ZZ(J200)) / (ZZ(J200+3)-ZZ(J200)) )
C
      
      ENDIF
C
C
C Calculate molecular ion densities and excited species densites:

     
      DO 300 I=1,JMAX
C
C
C N2+:
C
      IF (KCHEM .GE. 2) THEN
        P(5,I)= RN2PI(I)
     >        + (1.-B(15)) * RN2EI(I)
     >        + KZ(12,I) * ZXDEN(2,I) * ZN2(I)
     >        + KZ(16,I) * ZXDEN(1,I) * ZN2(I)
        L(5,I)= KZ(3,I)  * ZO(I)
     >        + KZ(4,I)  * ZO2(I)
     >        + KZ(23,I) * ECALC(I)
     >        + KZ(39,I) * ZO(I)
          ZXDEN(5,I) = P(5,I) / L(5,I)
      ENDIF
C
C
C O2+:
C
      IF (KCHEM .GE. 3) THEN
        P(6,I)= O2PI(I)
     >        + (1.-B(12)-B(13)-B(14)) * O2EI(I)
     >        + KZ(2,I)  * ZXDEN(3,I) * ZO2(I)
     >        + KZ(13,I) * ZXDEN(2,I) * ZO2(I)
     >        + KZ(17,I) * ZXDEN(1,I) * ZO2(I)
     >        + KZ(4,I)  * ZXDEN(5,I) * ZO2(I)
     >        + B(22) * KZ(25,I) * ZXDEN(4,I) * ZO2(I)
        L(6,I)= KZ(24,I) * ECALC(I)
     >        + KZ(30,I) * ZNO(I)
     >        + KZ(37,I) * ZNS(I)
        ZXDEN(6,I) = P(6,I)/ L(6,I)
      ENDIF
C
C
C NO+:
C
      IF (KCHEM .GE. 3) THEN
        P(7,I)= PHONO(1,I)
     >        + KZ(1,I)  * ZXDEN(3,I) * ZN2(I)
     >        + KZ(40,I) * ZXDEN(3,I) * ZNO(I)
     >        + KZ(3,I)  * ZXDEN(5,I) * ZO(I)
     >        + B(21) * KZ(25,I) * ZXDEN(4,I) * ZO2(I)
     >        + KZ(30,I) * ZXDEN(6,I) * ZNO(I)
     >        + KZ(37,I) * ZXDEN(6,I) * ZNS(I)
        L(7,I)= KZ(22,I) * ECALC(I)
        ZXDEN(7,I) = P(7,I) / L(7,I)
      ENDIF
      

C
C
C N2(A):
C
      P(8,I)= AGLW(1,3,I) + AGLW(2,3,I) + B(43)*AGLW(3,3,I)
      L(8,I)= KZ(26,I) * ZO(I)
     >      + KZ(29,I) * ZO2(I)
     >      + A(10)
      ZXDEN(8,I) = P(8,I) / L(8,I)
C
C
C N(2P):
C
      P(9,I)= B(28) * PHOTOD(1,3,I)
     >      + B(28) * PHOTOI(6,3,I)
     >      + B(26) * RN2ED(I)
     >      + B(23) * KZ(23,I) * ZXDEN(5,I) * ECALC(I)
      L(9,I)= KZ(33,I) * ZO(I)
     >      + KZ(34,I) * ZO2(I)
     >      + KZ(35,I) * ZNO(I)
     >      + A(11)
     >      + A(12)
      ZXDEN(9,I) = P(9,I) / L(9,I)
C
C
C N(2D):
C
      P(10,I)= B(27) * PHOTOD(1,3,I)
     >       + B(27) * PHOTOI(6,3,I)
     >       + B(25) * RN2ED(I)
     >       + B(16) * B(15) * RN2EI(I)
     >       + B(3)  * KZ(22,I) * ZXDEN(7,I) * ECALC(I)
     >       + B(4)  * KZ(23,I) * ZXDEN(5,I) * ECALC(I)
     >       + B(5)  * KZ(3,I)  * ZXDEN(5,I) * ZO(I)
     >       + B(29) * KZ(33,I) * ZXDEN(9,I) * ZO(I)
     >       + A(12) * ZXDEN(9,I)
      L(10,I)= KZ(5,I)  * ZO2(I)
     >       + KZ(6,I)  * ZO(I)
     >       + KZ(7,I)  * ECALC(I)
     >       + KZ(31,I) * ZNO(I)
     >       + KZ(38,I) * ZXDEN(3,I)
     >       + A(1)
      ZXDEN(10,I) = P(10,I) / L(10,I)

C
C
C O(1S):
C
      BZ(1,I) = 0.12 + 0.02 * ALOG10 (ECALC(I)/ZO(I)*(300./ZTE(I))**0.7)
      IF (BZ(1,I) .LT. 0.03) BZ(1,I)=0.03
      P(11,I)= AGLW(2,1,I)
     >       + BZ(1,I) * KZ(24,I) * ZXDEN(6,I)  * ECALC(I)
     >       + B(18) * KZ(26,I) * ZXDEN(8,I) * ZO(I)
     >       + B(33) * KZ(37,I) * ZXDEN(6,I) * ZNS(I)
     >       + B(34) * KZ(31,I) * ZXDEN(10,I) * ZNO(I)
     >       + PHOTOD(2,2,I)
      L(11,I)= KZ(11,I) * ZO(I)
     >       + KZ(36,I) * ZO2(I)
     >       + A(5)
     >       + A(4)
      ZXDEN(11,I) = P(11,I) / L(11,I)
C
C
C O(1D):
C
      P(12,I)= AGLW(1,1,I)
     >       + KZ(28,I) * ECALC(I)  * ZO(I)
     >       + B(2)  * KZ(24,I) * ZXDEN(6,I)  * ECALC(I)
     >       + B(6)  * KZ(5,I)  * ZXDEN(10,I) * ZO2(I)
     >       + B(20) * KZ(6,I)  * ZXDEN(10,I) * ZO(I)
     >       + B(17) * KZ(25,I) * ZXDEN(4,I)  * ZO2(I)
     >       + B(7)  * KZ(15,I) * ZXDEN(2,I)  * ZO(I)
     >       + SRCED(I)
     >       + PHOTOD(1,2,I)
     >       + A(5)  * ZXDEN(11,I)
      L(12,I)= KZ(8,I)  * ZN2(I) 
     >       + KZ(9,I)  * ZO2(I)
     >       + KZ(10,I) * ECALC(I)
     >       + KZ(27,I) * ZO(I)
     >       + A(2)
     >       + A(3)
      ZXDEN(12,I) = P(12,I) / L(12,I)
C
  300 CONTINUE
C
  350 CONTINUE
C
C
C Calculate airglow emission rates; fill ZCETA array with partial rates
C from each source; fill ZETA array with total rate for each emission:
C
    
      DO 400 I=1,JMAX
C
      ZCETA(1,1,I) = B(39) * AGLW(3,3,I)
      ZCETA(2,1,I) = B(40) * A(10) * P(8,I) / L(8,I)
C
      ZCETA(1,2,I) = B(38) * B(36) * RN2EI(I)
      ZCETA(2,2,I) = B(38) * PHOTOI(3,3,I)
      ZCETA(3,2,I) = G(2,I) * ZXDEN(5,I)
C
      ZCETA(1,3,I) = A(1) * B(27) * PHOTOD(1,3,I) / L(10,I)
      ZCETA(2,3,I) = A(1) * B(27) * PHOTOI(6,3,I) / L(10,I)
      ZCETA(3,3,I) = A(1) * B(25) * RN2ED(I) / L(10,I)
      ZCETA(4,3,I) = A(1) * B(16) * B(15) * RN2EI(I) / L(10,I)
      ZCETA(5,3,I) = A(1) * B(3)  * KZ(22,I) * ZXDEN(7,I) * ECALC(I) /L(10,I)
      ZCETA(6,3,I) = A(1) * B(4)  * KZ(23,I) * ZXDEN(5,I) * ECALC(I) /L(10,I)
      ZCETA(7,3,I) = A(1) * B(5)  * KZ(3,I)  * ZXDEN(5,I) * ZO(I) /L(10,I)
      ZCETA(8,3,I) = A(1) * B(29) * KZ(33,I) * ZXDEN(9,I) * ZO(I) /L(10,I)
      ZCETA(9,3,I) = A(1) * A(12) * ZXDEN(9,I) / L(10,I)
C
      ZCETA(1,4,I) = A(5) * AGLW(2,1,I) / L(11,I)
      ZCETA(2,4,I) = A(5) * BZ(1,I)*KZ(24,I) * ZXDEN(6,I) * ECALC(I)  /L(11,I)
      ZCETA(3,4,I) = A(5) * B(18) * KZ(26,I) * ZXDEN(8,I) * ZO(I) /L(11,I)
      ZCETA(4,4,I) = A(5) * B(33) * KZ(37,I) * ZXDEN(6,I) * ZNS(I)/L(11,I)
      ZCETA(5,4,I) = A(5) * B(34) * KZ(31,I) * ZXDEN(10,I)* ZNO(I)/L(11,I)
      ZCETA(6,4,I) = PHOTOD(2,2,I) / L(11,I)
C
      ZCETA(1,5,I) = A(2) * AGLW(1,1,I) / L(12,I)
      ZCETA(2,5,I) = A(2) * KZ(28,I) * ECALC(I)  * ZO(I) / L(12,I)
      ZCETA(3,5,I) = A(2) * B(2)  * KZ(24,I) * ZXDEN(6,I)  * ECALC(I)/L(12,I)
      ZCETA(4,5,I) = A(2) * B(6)  * KZ(5,I)  * ZXDEN(10,I) *ZO2(I)/L(12,I)
      ZCETA(5,5,I) = A(2) * B(20) * KZ(6,I)  * ZXDEN(10,I) * ZO(I)/L(12,I)
      ZCETA(6,5,I) = A(2) * B(17) * KZ(25,I) * ZXDEN(4,I) * ZO2(I)/L(12,I)
      ZCETA(7,5,I) = A(2) * B(7)  * KZ(15,I) * ZXDEN(2,I)  * ZO(I)/L(12,I)
      ZCETA(8,5,I) = A(2) * SRCED(I) / L(12,I)
      ZCETA(9,5,I) = A(2) * PHOTOD(1,2,I) / L(12,I)
      ZCETA(10,5,I)= A(2) * A(5)  * ZXDEN(11,I) / L(12,I)
C
      ZCETA(1,6,I) = A(8) * (PHOTOI(3,1,I)+B(41)*PHOTOI(5,1,I)) / L(1,I)
      ZCETA(2,6,I) = A(8) * B(30) * PHOTOI(4,2,I) / L(1,I)
      ZCETA(3,6,I) = A(8) * B(9) * OEI(I) / L(1,I)
      ZCETA(4,6,I) = A(8) * B(12) * O2EI(I) / L(1,I)
C
      ZCETA(1,7,I) = A(12) * B(28) * PHOTOD(1,3,I) / L(9,I)
      ZCETA(2,7,I) = A(12) * B(28) * PHOTOI(6,3,I) / L(9,I)
      ZCETA(3,7,I) = A(12) * B(26) * RN2ED(I) / L(9,I)
      ZCETA(4,7,I) = A(12) * B(23) * KZ(23,I) * ZXDEN(5,I) * ECALC(I) / L(9,I)
C
      ZCETA(1,8,I) = A(11) * B(28) * PHOTOD(1,3,I) / L(9,I)
      ZCETA(2,8,I) = A(11) * B(28) * PHOTOI(6,3,I) / L(9,I)
      ZCETA(3,8,I) = A(11) * B(26) * RN2ED(I) / L(9,I)
      ZCETA(4,8,I) = A(11) * B(23) * KZ(23,I) * ZXDEN(5,I) * ECALC(I) / L(9,I)
C
      ZCETA(1,9,I) = AGLW(5,1,I)
C
      ZCETA(1,10,I) = AGLW(6,1,I)
      ZCETA(2,10,I) = AGLW(7,1,I)
      ZCETA(3,10,I) = B(44) * AGLW(8,1,I)
C
      ZCETA(1,11,I) = A(6) * (PHOTOI(2,1,I)+B(42)*PHOTOI(5,1,I))/ L(2,I)
      ZCETA(2,11,I) = A(6) * B(31) * PHOTOI(4,2,I) / L(2,I)
      ZCETA(3,11,I) = A(6) * B(10) * OEI(I) / L(2,I)
      ZCETA(4,11,I) = A(6) * B(13) * O2EI(I) / L(2,I)
      ZCETA(5,11,I) = A(6) * B(8)  * KZ(18,I) * ZXDEN(1,I) * ECALC(I) / L(2,I)
      ZCETA(6,11,I) = A(6) * A(8)  * ZXDEN(1,I) / L(2,I)
C
      ZETA(1,I)  = ZCETA(1,1,I)+ZCETA(2,1,I)
      ZETA(2,I)  = ZCETA(1,2,I)+ZCETA(2,2,I)+ZCETA(3,2,I)
      ZETA(3,I)  = ZCETA(1,3,I)+ZCETA(2,3,I)+ZCETA(3,3,I)
     >            +ZCETA(4,3,I)+ZCETA(5,3,I)+ZCETA(6,3,I)
     >            +ZCETA(7,3,I)+ZCETA(8,3,I)+ZCETA(9,3,I)
      ZETA(4,I)  = ZCETA(1,4,I)+ZCETA(2,4,I)+ZCETA(3,4,I)
     >            +ZCETA(4,4,I)+ZCETA(5,4,I)+ZCETA(6,4,I)
      ZETA(5,I)  = ZCETA(1,5,I)+ZCETA(2,5,I)+ZCETA(3,5,I)
     >            +ZCETA(4,5,I)+ZCETA(5,5,I)+ZCETA(6,5,I)
     >            +ZCETA(7,5,I)+ZCETA(8,5,I)+ZCETA(9,5,I)
     >            +ZCETA(10,5,I)
      ZETA(6,I)  = ZCETA(1,6,I)+ZCETA(2,6,I)+ZCETA(3,6,I)+ZCETA(4,6,I)
      ZETA(7,I)  = ZCETA(1,7,I)+ZCETA(2,7,I)+ZCETA(3,7,I)+ZCETA(4,7,I)
      ZETA(8,I)  = ZCETA(1,8,I)+ZCETA(2,8,I)+ZCETA(3,8,I)+ZCETA(4,8,I)
      ZETA(9,I)  = ZCETA(1,9,I)
      ZETA(10,I) = ZCETA(1,10,I)+ZCETA(2,10,I)+ZCETA(3,10,I)
      ZETA(11,I) = ZCETA(1,11,I)+ZCETA(2,11,I)+ZCETA(3,11,I)
     >            +ZCETA(4,11,I)+ZCETA(5,11,I)+ZCETA(6,11,I)
      
     
C        if(ipc.eq.8) then
C         if (GLAT .le. -2 .and.GLAT.gt.-3) then
C            if(GLONG .gt. 2 .and.  GLONG .lt. 3) then 
C               if(ZETA(5,I).gt.1000) then
C                  write(*,*) 'in GLLIB:', GLAT, GLONG
CC     !do 610 J = 1,JMAX
C                  write(*,*)  ZXDEN(6,I),ZXDEN(4,I),ZXDEN(10,I), ZXDEN(2,I)
C                  write(*,*) ZCETA(1,5,I),ZCETA(12,5,I),ZCETA(3,5,I),ZCETA(4,5,I),
C     >             ZCETA(5,5,I),ZCETA(6,5,I),ZCETA(7,5,I),ZCETA(8,5,I),
C     >             ZCETA(9,5,I),ZCETA(10,5,I)
C                  write(*,*) 'srced: ',SRCED(I), AGLW(4,2,I),B(35),PIA(2,I)
C               endif
CC 610           continue
C            endif
C         endif
C      endif

  400 CONTINUE


C      write(*,*) ZN2(22) 
C     >        , ZO2(22)
C     >        , ECALC(22)
C     >        , ZO(22)


!A(2) , B(2)  , KZ(24,22) , ZXDEN(6,22)  , ECALC(22),L(12,22)

! ZCETA(1,5,22),ZCETA(2,5,22),ZCETA(3,5,22)
!     >            ,ZCETA(4,5,22),ZCETA(5,5,22),ZCETA(6,5,22)
!     >            ,ZCETA(7,5,22),ZCETA(8,5,22),ZCETA(9,5,22)
!     >            ,ZCETA(10,5,22),ZETA(5,22)
!      stop

C
C Calculate vertical column brightnesses:
C
      DO 450 IW=1,NW
      VCB(IW) = 0.
  450 CONTINUE
C
      DO 550 I=1,JMAX
      IF (I .EQ. JMAX) THEN
        DZ = (ZZ(I) - ZZ(I-1))
      ELSE
        IF (I .EQ. 1) THEN
          DZ = (ZZ(I+1) - ZZ(I))
        ELSE
          DZ = (ZZ(I+1) - ZZ(I-1)) / 2.0
        ENDIF
      ENDIF
      DO 500 IW=1,NW
      VCB(IW) = VCB(IW) + ZETA(IW,I) * DZ
  500 CONTINUE
  550 CONTINUE
       
C
C
C Convert brightnesses to Rayleighs:
C
      DO 600 IW=1,NW
        VCB(IW) = VCB(IW) / 1.E6
  600 CONTINUE
     

C
C
      RETURN
      END
