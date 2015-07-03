C Subroutine SSFLUX
C
C This software is part of the GLOW model.  Use is governed by the Open Source
C Academic Research License Agreement contained in the file glowlicense.txt.
C For more information see the file glow.txt.
C
C Subroutine SSFLUX calculates the solar EUV and FUV flux in the range
C 0.5 to 1750 Angstroms for a specified level of solar activity.
C
C The calling routine supplies a scaling switch ISCALE, the daily 10.7
C cm flux F107, its 81-day centered average F107A, and, optionally,
C the H Lyman-Beta 1026A ratio to its solar minimum value HLYBR, the
C Fe XVI 335A ratio to its solar minimum value FEXVIR, the H Lyman-alpha
C 1216A flux HLYA, the He I 10830A equivalent width HEIEW, and an XUV
C enhancement factor XUVFAC.  Any optional calling parameters not used
C should be set to zero.
C
C XUVFAC is applied from 18-250 A for the Hinteregger model (ISCALE=0),
C from 18-50 A for the EUVAC model (ISCALE=1), and not at all for
C user-supplied data (ISCALE=2)
C
C The subroutine returns the longwave boundary WAVE1 and shortwave
C boundary WAVE2 of the wavelenth bins, and the solar flux in each bin
C SFLUX.  Bins are arranged in energy multiples of 2 from 0.5 to 8 A,
C aligned with k-shell boundaries at 23, 32, and 44 A from 18 to 44 nm,
C 10 A in width from 60 to 1050 A, and 50 A in width from 1050 to
C 1750 A with the exception of Lyman-alpha which has its own bin from
C 1210 to 1220 A. 
C
C Methods used:
C   If ISCALE=0 the flux is scaled using parameterization methods based
C on F107 and F107A.  For ionizing EUV, Hinteregger's contrast ratio
C method (Hinteregger et al., GRL, 8, 1147, 1981) is used, based on the
C reference spectrum SC#21REFW at 1 nm bin resolution.  If the
C H Lyman-Beta (1026A) or Fe XVI (335A) enhancement ratios are provided
C (>0) as calling arguments, they are used to scale the EUV spectrum.
C Otherwise, enhancement ratios for H Ly-B and Fe XVI are calculated
C from F107 and F107A using Hinteregger's formula, employing
C coefficients which reduce to the reference values at F107=67.6,
C F107A=71.5.  The 'best fit' coefficients are not used as they produce
C some negative values at low solar activity, but remain in a
C 'commented out' data statement for reference.  The EUV spectrum is
C then scaled from these modeled ratios.  Scaling factors were
C calculated from contrast ratios in the SC#21REFW data file.
C   If ISCALE=1, the EUV flux (50-1050A) is scaled using the EUVAC model
C (Richards et al., JGR 99, 8981, 1994) re-binned onto ~1 nm intervals.
C The Hinteregger spectrum, scaled using the EUVAC algorithm, is used
C from 18 to 50A.
C   Neither of these models extends shortward of 18A, so from 1-18 A
C an amalgam of sources are used to derive an estimated flux, e.g.,
C DeJager, in Astronomical Observations from Space Vehicles, Steinberg,
C ed., 1964; Smith & Gottlieb, SSR 16, 771, 1974; Manson, in The Solar
C Output and its Variation, White, ed., 1977; Kreplin et al, ibid;
C Horan & Kreplin, Solar Physics 74, 265, 1981; Wagner, Adv. Space Res.
C 8, (7)67, 1988.
C    For FUV from 1050A-1750A, 50A interval bins from the Woods and
C Rottman [2002] reference spectrum and scale factors based on
C UARS SOLSTICE data are used.  The scaling method follows the
C Hinteregger or EUVAC algorithm, whichever is selected, so as to
C linearly scale the spectrum between the reference value and maximum
C value calculated with F10.7=F10.7A=200.  If a value for Lyman-alpha
C (HLYA>0) is provided by the calling program, it is subsituted into
C the spectrum.
C   If ISCALE=2, the solar flux (0-1750A) is read from a file named
C ssflux_user.dat in the current working directory.  The file must
C contain three columns:  WAVES, WAVEL, SFLUX (Angstroms and cm-2 s-1)
C in order of increasing wavelength.  The number of lines in the file
C must match the value of LMAX in glow.h.
C
C Modification history:
C   Stan Solomon, 12/88  Basic Hinteregger EUV, approx. SME FUV
C   Chris Gaskill, 7/89  Added early Tobiska model
C   Stan Solomon,  8/89  Corrections to above
C   Stan Solomon,  1/90  Tobiska SERF2; added W & R spectra
C   Stan Solomon,  6/91  Tobiska EUV 91; Hntggr Ly-B, Fe XVI scaling
C   Stan Solomon,  2/92  Updated Tobiska EUV91; corrected SME FUV
C   Scott Bailey, 12/93  Initial one-nm bins version
C   Stan Solomon,  6/04  Added EUVAC option, cleaned up artifacts
C   Stan Solomon,  9/04  Added ability to specify input data file
C   Stan Solomon,  3/05  Changed all to photon units
C
C Calling parameters:
C ISCALE   =0 for Hinteregger contrast ratio method
C          =1 for EUVAC
C          =2 for user-supplied data
C F107     daily 10.7 cm flux (1.E-22 W m-2 Hz-1)
C F107A    81-day centered average 10.7 cm flux
C HLYBR    ratio of H Ly-b 1026A flux to solar minimum value (optional)
C FEXVIR   ratio of Fe XVI 335A flux to solar minimum value (optional)
C HLYA     H Lyman-alpha flux (photons cm-2 s-1) (optional)
C HEIEW    He I 10830A equivalent width (mAngstroms) (obsolete)
C XUVFAC   factor for scaling flux 18-250A or 18-50A (optional)

C Returned parameters:
C WAVE1    longwave bound of spectral intervals (Angstroms)
C WAVE2    shortwave bound of intervals
C SFLUX    scaled solar flux returned by subroutine (photons cm-2 s-1)
C
C Other definitions:
C LMAX     dimension of flux and scaling arrays, currently = 123
C WAVEL    = WAVE1
C WAVES    = WAVE2
C RFLUX    low solar activity flux
C XFLUX    high solar activity flux
C SCALE1   scaling factors for H LyB-keyed chromospheric emissions
C SCALE2   scaling factors for FeXVI-keyed coronal emissions
C B1       fit coefficients for H LyB
C B2       fit coefficients for FeXVI
C R1       enhancement ratio for H LyB
C R2       enhancement ratio for FeXVI
C P107     average of F107 and F107A
C A        scaling factor for EUVAC model
C
C
      SUBROUTINE SSFLUX (F107, F107A, HLYBR, FEXVIR, HLYA,
     >                   HEIEW, XUVFAC, WAVE1, WAVE2, SFLUX)
C
      use Mod_GLOW, only: JMAX,NBINS,LMAX,ISCALE

C
      DIMENSION WAVE1(LMAX), WAVE2(LMAX), SFLUX(LMAX),
     >          WAVEL(LMAX), WAVES(LMAX), RFLUX(LMAX),
     >          SCALE1(LMAX), SCALE2(LMAX), A(LMAX), B1(3), B2(3) 
      DATA EPSIL/1.0E-6/

C     
C regression coefficients which reduce to solar min. spectrum:
      DATA B1/1.0, 0.0138, 0.005/, B2/1.0, 0.59425, 0.3811/
C
C 'best fit' regression coefficients, commented out, for reference:
C     DATA B1/1.31, 0.01106, 0.00492/, B2/-6.618, 0.66159, 0.38319/
C
C
C Hinteregger contrast ratio method:
C
C      write(*,*) ISCALE,F107,F107A,FEXVIR,HLYA,HEIEW,xuvfac
C      stop

      IF (ISCALE .EQ. 0) THEN
        
         open(unit=1,file='../srcGlow/ssflux_hint.dat',status='old')
         read(1,*)
         do 40,l=lmax,1,-1
            read(1,*) waves(l),wavel(l),rflux(l),scale1(l),scale2(l)
 40      continue
         close(unit=1)
C     
         IF (HLYBR .GT. EPSIL) THEN
            R1 = HLYBR
         ELSE
            R1 =  B1(1) + B1(2)*(F107A-71.5) + B1(3)*(F107-F107A+3.9)
         ENDIF
         IF (FEXVIR .GT. EPSIL) THEN
            R2 = FEXVIR
         ELSE
            R2 =  B2(1) + B2(2)*(F107A-71.5) + B2(3)*(F107-F107A+3.9)
         ENDIF
C     
         DO 100 L=1,LMAX
            SFLUX(L) = RFLUX(L) + (R1-1.)*SCALE1(L) + (R2-1.)*SCALE2(L)
            IF (SFLUX(L) .LT. 0.0) SFLUX(L) = 0.0
            IF (XUVFAC .GT. EPSIL .AND.
     >           WAVEL(L).LT.251.0 .AND. WAVES(L).GT.17.0)
     >           SFLUX(L)=SFLUX(L)*XUVFAC
 100     CONTINUE
      ENDIF
C     
C EUVAC Method:
C
      IF (ISCALE .EQ. 1) THEN
        open(unit=1,file='../srcGlow/ssflux_euvac.dat',status='old')
        read(1,*)
        do 200,l=lmax,1,-1
          read(1,*) waves(l),wavel(l),rflux(l),a(l)
  200   continue
        close(unit=1)
C
      P107 = (F107+F107A)/2.
        DO 300 L=1,LMAX
          SFLUX(L) = RFLUX(L) * (1. + A(L)*(P107-80.))
          IF (SFLUX(L) .LT. 0.8*RFLUX(L)) SFLUX(L) = 0.8*RFLUX(L)
          IF (XUVFAC .GT. EPSIL .AND.
     >        WAVEL(L).LT.51.0 .AND. WAVES(L).GT.17.0)
     >        SFLUX(L)=SFLUX(L)*XUVFAC
  300   CONTINUE

      ENDIF
C
C User-supplied data:
C
      IF (ISCALE .EQ. 2) THEN
         open(unit=1,file='~/GLOWGITM/srcGlow/ssflux_user.dat',status='old')
         read(1,*)
         do 400,l=lmax,1,-1
            read(1,*) waves(l),wavel(l),sflux(l)
 400     continue
         close(unit=1)
      ENDIF

     
C     
C Fill wavelength arrays, substitute in H Lyman-alpha if provided:
C
      DO 600 L=1,LMAX
        WAVE1(L) = WAVEL(L)
        WAVE2(L) = WAVES(L)
        IF (HLYA .GT. EPSIL .AND.
     >      WAVEL(L).LT.1221. .AND. WAVES(L).GT.1209.)
     >      SFLUX(L) = HLYA
  600 CONTINUE
C

      RETURN
      END
