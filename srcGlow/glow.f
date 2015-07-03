C Subroutine GLOW
C
C This software is part of the GLOW model.  Use is governed by the Open Source
C Academic Research License Agreement contained in the file glowlicense.txt.
C For more information see the file glow.txt.
C
C Version 0.97
C
C Stan Solomon, 1988, 1989, 1990, 1991, 1992, 1994, 2000, 2002, 2005
C
C Subroutine GLOW is the master routine of the /glow package.  It
C receives input parameters from the calling program in common block
C /CGLOW/, calls the other subroutines, and returns results in /CGLOW/.
C Header file glow.h supplies array sizes for the altitude, electron,
C and solar spectrum energy grids.
C
C Subroutines called by GLOW are:
C   EGRID   sets up electron energy grid
C   FIELDM  calculates magnetic dip angle
C   SOLZEN  calculates solar zenith angle
C   SSFLUX  scales solar flux for activity level
C   RCOLUM  calculates slant column density of major species
C   EPHOTO  calculates photoionization and photoelectron production
C   QBACK   estimates background ionization
C   ETRANS  computes electron transport, ionization, excitation
C             calls EXSECT for cross-sections, first call only
C   GCHEM   finds electron/ion/metastable densities, airglow emissions
C             uses VQUART to solve electron density equation
C
C Supplied to subroutine in labeled common /CGLOW/:
C IDATE   Date, in form yyddd
C UT      Universal Time; seconds
C GLAT    Geographic latitude; degrees
C GLONG   Geographic longitude; degrees
C ISCALE  Solar flux scaling switch, see subroutine SSFLUX
C JLOCAL  =0 for electron transport calculation, =1 for local calc only
C KCHEM   Ion/electron chemistry switch, see subroutine GCHEM
C F107    Solar 10.7 cm flux for day being modeled, 1.E-22 W m-2 Hz-1
C F107A   Solar 10.7 cm flux 81-day centered average
C HLYBR   H Ly-b (1026A) enhancement ratio
C FEXVIR  Fe XVI (335A) enhancement ratio
C HLYA    H Ly-a flux; photons cm-2 s-1
C HEIEW   He I 10830 equivalent width; milliAngstroms (obsolete)
C XUVFAC  Factor by which to multiply to solar flux 16-250 A or 16-50 A.
C ZZ      altitude array; cm
C ZO      O number density at each altitude; cm-3
C ZN2     N2  "      "      "   "     "       "
C ZO2     O2         "
C ZNO     NO         "
C ZNS     N(4S)      "
C ZND     N(2D)      "
C ZRHO    mass density at each altitude; gm cm-3
C ZE      electron density at each alt; cm-3
C ZTN     neutral temperature at each alt; K
C ZTI     ion temperature at each alt; K
C ZTE     electron temp at each alt; K
C PHITOP  energetic electron flux into top of atmosphere; cm-2 s-1 eV-1
C EFLUX   obsolete
C EZERO   obsolete
C
C Calculated by subroutine:
C SZA     solar zenith angle; radians
C DIP     magnetic field dip angle; radians
C EFRAC   energy conservation check from ETRANS, (out-in)/in
C IERR    error code returned from ETRANS:
C           0=normal, 1=local problem, 2=transport problem
C ZMAJ    major species density array, O, O2, N2; cm-3
C ZCOL    major species slant column density array, O, O2, N2; cm-2
C WAVE1   longwave edge of solar flux wavelength range; A
C WAVE2   shortwave edge of solar flux wavelength range; A
C SFLUX   scaled solar flux in each wavelength range; photons cm-2 s-1
C ENER    electron energy grid; eV
C DEL     width of each bin in electron energy grid; eV
C PESPEC  photoelectron production rate at energy, altitude; cm-3 s-1
C SESPEC  background electron production rate (obsolete); cm-3 s-1
C PHOTOI  photoionization rates for state, species, altitude; cm-3 s-1
C           O+ states: 4S, 2Do, 2Po, 4Pe, 2Pe
C           O2+ states: X, a+A, b, dissoc.
C           N2+ states: X, A, B, C, F, dissoc.
C PHOTOD  photodissoc. & exc. rates for state, species, alt.; cm-3 s-1
C           (1,2,J) = O2 -> O(3P) + O(1D))
C           (2,2,J) = O2 -> O(3P) + O(1S)
C           (1,3,J) = N2 -> N + N
C PHONO   photoionization/dissociation/excitation rates for NO, cm-3 s-1
C         (1,J) = NO+ from H Ly-a ionization
C QTI     obsolete
C AURI    obsolete
C PIA     obsolete
C SION    electron impact ioniz. rates calculated by ETRANS; cm-3 s-1
C UFLX    upward hemispherical electron flux; cm-2 s-1 eV-1
C DFLX    downward hemispherical electron flux; cm-2 s-1 eV-1
C AGLW    Electron impact exc. rates; state, species, alt.; cm-3 s-1
C           O states: 1D, 1S, 5S, 3S, 3p5P, 3p3P, 3d3D, 3s'3D
C           O2 states: a, b, (A+A'+c), B(SRC), 9.9eV, Ryds., vib.
C           N2 states: (A+B+W), B', C, (a+a'+w), 1Pu, b', Ryds., vib.
C EHEAT   ambient electron heating rate, eV cm-3 s-1
C TEZ     total energetic electron energy deposition, eV cm-3 s-1
C ECALC   electron density, calculated below 200 km, cm-3
C ZXDEN   array of excited and and/or ionized state densities:
C           O+(2P), O+(2D), O+(4S), N+, N2+, O2+, NO+, N2(A), N(2P),
C           N(2D), O(1S), O(1D), 8 spares, at each altitude; cm-3
C ZETA    array of volume emission rates:
C           3371A, 4278A, 5200A, 5577A, 6300A, 7320A, 10400A, 3466A,
C           7774A, 8446A, 3726A, 9 spares; cm-3 s-1
C ZCETA   array of contributions to each v.e.r at each alt; cm-3 s-1
C VCB     array of vertical column brightnesses (as above); Rayleighs
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
C NF      obsolete
C
C
      SUBROUTINE GLOW
C

      USE Mod_GLOW
     
C
      DIMENSION ZVCD(NMAJ,JMAX)
C
C      DATA IFIRST/1/
C
C
C First call only: set up energy grid:
C
      
      IF (GFIRST .EQ. 1) THEN
         CALL EGRID (ENER, DEL, NBINS)
      ENDIF
       
C
C
C Find magnetic dip angle and solar zenith angle (radians):
C
      
      CALL FIELDM (GLAT, GLONG, 300., XF, YF, ZF, FF, DIP, DEC, SDIP)
      DIP = ABS(DIP) * PI/180.
C
      
      CALL SOLZEN (IDATE, UT, GLAT, GLONG, SZA)
      
      
      SZA = SZA * PI/180.
      
      
C
C !!!!!! The GITM SOLAR FLUX IS NOW USED !!!!!!!!!
C Scale solar flux:  
C
C
     
C      CALL SSFLUX (F107, F107A, HLYBR, FEXVIR, HLYA,
C     >               HEIEW, XUVFAC, WAVE1, WAVE2, SFLUX)

C
C
C Pack major species density array:
C
     
      DO 100 J=1,JMAX
        ZMAJ(1,J) = ZO(J)
        ZMAJ(2,J) = ZO2(J)
        ZMAJ(3,J) = ZN2(J)
  100 CONTINUE


C
C
C Calculate slant path column densities of major species in the
C direction of the sun:
C


      CALL RCOLUM (SZA, ZZ, ZMAJ, ZTN, ZCOL, ZVCD, JMAX, NMAJ)

C
C
C Call subroutine EPHOTO to calculate the photoelectron production
C spectrum and photoionization rates as a function of altitude,
C unless all altitudes are dark, in which case zero arrays:
C

      IF (SZA .LT. 2.) THEN
         CALL EPHOTO
         EPFIRST = 0
      ELSE
        DO 240 J=1,JMAX
          DO 200 I=1,NMAJ
          DO 200 IST=1,NST
            PHOTOI(IST,I,J) = 0.0
            PHOTOD(IST,I,J) = 0.0
  200     CONTINUE
          DO 210 IST=1,NST
            PHONO(IST,J) = 0.0
  210     CONTINUE
          DO 220 N=1,NBINS
            PESPEC(N,J) = 0.0
  220     CONTINUE
 240   CONTINUE
      ENDIF

C
C
C Zero obsolete secondary spectrum:
C


      DO 245 J=1,JMAX
      DO 245 N=1,NBINS
        SESPEC(N,J) = 0.0
  245 CONTINUE

C
C 
C Add background ionization to photoionization:
C
      CALL QBACK (ZMAJ, ZNO, ZVCD, PHOTOI, PHONO, JMAX, NMAJ, NST)
C
C
C Call subroutine ETRANS to calculate photoelectron and auroral
C electron transport and electron impact excitation rates, unless
C there are no energetic electrons, in which case zero arrays:
C

      TEFLUX = 0.
      DO 260 N=1,NBINS
        TEFLUX = TEFLUX + PHITOP(N)
  260 CONTINUE

C

      IF (TEFLUX .GT. 0.001 .OR. SZA .LT. 2.) THEN
         CALL ETRANS
         ETFIRST = 0
      ELSE
        DO 350 J=1,JMAX
          DO 320 N=1,NBINS
             UFLX(N,J) = 0.0
             DFLX(N,J) = 0.0
  320     CONTINUE
          DO 340 I=1,NMAJ
            SION(I,J) = 0.0
            DO 340 IEI=1,NEI
              AGLW(IEI,I,J) = 0.0
  340     CONTINUE
          EHEAT(J) = 0.0
          TEZ(J) = 0.0
 350   CONTINUE
        EFRAC = 0.0
        IERR = 0
      ENDIF

C
C
C Call subroutine GCHEM to calculate the densities of excited and
C ionized consituents, airglow emission rates, and vertical column
C brightnesses:
C



      CALL GCHEM
             
C
C
      GFIRST = 0
      RETURN
      END
