      SUBROUTINE MHEMODL (RMLAT,RMLT,HP,BY,BZ,MODL,ET,EP,EPOT)
C          Millstone Hill electric potential models.  This subroutine
C          determines the electric potential for a specified location
C          poleward of 55 or -55 degrees magnetic latitude.  One of two
C          models (MHI or MHS) may be used for the estimation.  Also
C          determined are the electic field components, since they are
C          -grad(electrostatic potential).
C
C          A prerequisite:  call MHINIT to initialize constants and load
C          model coefficients into memory prior calling MHEMODL.
C
C          MHI model:
C          Millstone Hill incoherent scatter scatter (IS) radar data
C          combined with the Hemispheric Power (HP) defined by Fuller-Rowell
C          and Evans (JGR, 92, pp 7606-7618, 1987).  An integer index (HPI)
C          from 1 to 10 is approximately related to the hemispheric power
C          from auroral electrons, HP in GW by
C
C               HP = 4.2**(HPI/3)
C
C          where HP is the median value in the index level HPI.  The
C          inverse of this formula is:
C
C               HPI = 2.09 * ln(HP)
C
C          However, these are approximations to the actual boundaries
C
C               2.50, 3.94, 6.22, 9.82, 15.49, 24.44, 38.56, 60.85, 96.00
C
C          in GW (David Evans, private communication, 1988).  The minumum
C          value in level 2 is 2.50 GW; the minimum in level 10 is 96.00 GW.
C          5% of the hemispheric power dataset is in level 1 and 1%-1.5%
C          is in level 10.
C
C          If the above table is not used, a better formula for HPI is:
C
C               IHPI = MAX0(1,MIN0(10,INT( 1.0475 * 2.09 * ln(HP))))
C
C          Because of the sparse statistics in level 10, electric potential
C          patterns were developed from Millstone Hill data for HPI levels
C          1-9 (Foster, Holt, Musgrove and Evans, GRL, 13, pp 656-659, 1986).
C          There is no IMF By dependence so the patterns are identical in
C          both hemispheres.
C
C          MHS model:
C          Millstone Hill and Sondestrom IS radar data were combined with
C          the sign of the By and Bz IMF components to produce four potential
C          patterns (described in Foster "Proceedings of the International
C          Symposium on Quantitative Modeling of Magnetosphere-Ionosphere
C          Coupling Processes, editors Kamide and Wolf, Kyoto, Japan,
C          pp 71-76, 1987).  The patterns were derived for the northern
C          hemisphere, so a southern hemisphere pattern can be inferred by
C          assuming an opposite sign for By.
C
C          INPUTS:
C            RMLAT   = magnetic latitude in degrees.
C            RMLTLON = magnetic local time in degrees (0=midnight).
C            RMLT    = magnetic local time in degrees (0=midnight).
C            HP      = The Hemispheric Power in GW
C            BY      = IMF GSM By in nT
C            BZ      = IMF GSM Bz in nT
C            MODL    = Model selection switch:
C                    = 1 Use HP model from Millstone
C                    = 2 Use Radar model based on Millstone/Sonde for 4 By/Bz
C                        If Bz=0 use Bz+, and if By=0 use By- (NH) or By+ (SH)
C          RETURNS:
C            ET    = Etheta (magnetic antipoleward*) E field component in V/m
C            EP    = Ephi   (magnetic eastward)      E field component in V/m
C            EPOT  = Electrostatic potential in kV.
C
C          * ET direction is along the magnetic meridian away from the
C            current hemisphere; i.e., the direction is
C                 southward when RMLAT > 0
C                 northward when RMLAT < 0
C            The MHI potential patterns and E components are otherwise
C            the same in both hemispheres.  The MHS model recognizes the
C            inverse relationship between the sign of BY and the affected
C            hemisphere by taking both into account when selecting the
C            potential pattern; so, when RMLAT < 0, MHS chooses the pattern
C            for the opposite sign of By.
      SAVE

C          MHRFCM contains model coefficients and is assigned in MHINIT:
      PARAMETER (MXNBETA=143,MXNNDX=9,MXNF=196)
      REAL PI,RE,TMN,TMX,RMLAMN,RMLAMX,TX,TY,BETA,F
      COMMON /MHRFCM/ NX,KX,NY,KY,NBETA,KFIT,PI,RE,TMN,TMX,RMLAMN(2),
     + RMLAMX(2),TX(100),TY(100,2),BETA(MXNBETA,MXNNDX,2),F(3,3,MXNF)

C          ADDPOT(NDX,MODL) = Electrostatic potential (kV) additive offset
C          required make a zero average at the lowest model latitude.
C          When MODL=1, NDX 1-9 corresponds to the 9 Hemispheric Power Input
C          patterns of the MHI model.  When MODL=2, NDX 1-4 corresponds to
C          the four MHS model By/Bz patterns.
      DIMENSION ADDPOT(9,2)
      DATA  ADDPOT /1.04, 2.08,7.54,10.11,4.11,10.37,6.29,5.39,5.00,
     +             14.62,22.66,0.00, 0.09,     999.,999.,999.,999.,999./

C          Declarations for GETFIT
      REAL  TIME,RMLA,XLA,DFIT,ESFIT,EEFIT

      IF (MODL .NE. 1 .AND. MODL .NE. 2) THEN
	WRITE (6,'(''MHEMODL: MODL .ne. 1 or 2; MODL='',I3)') MODL
	STOP
      ENDIF

      TIME = RMLTLON/15.
      TIME = RMLT
   20 IF (TIME .GT. TMX) THEN
	TIME = TIME - 24.
	GO TO 20
      ENDIF
   30 IF (TIME .LT. TMN) THEN
	TIME = TIME + 24.
	GO TO 30
      ENDIF

      H = SIGN (1.,RMLAT)
      RMLA = MIN(MAX(ABS(RMLAT),RMLAMN(MODL)), RMLAMX(MODL))

C         Millstone HP model for any By
      IF (MODL .EQ. 1) THEN
C          Avoid error with HP=0. when determining integer hemi pwr index
	P = AMAX1 (0.01,HP)
	NDX = MIN0 (9, MAX0 (1 , INT(1.0475*2.09*ALOG(P))))
      ENDIF

C         Millstone/Sonde 4 By/Bz patterns
      IF (MODL .EQ. 2) THEN
C          Now revert to Bz+ for Bz=0, and to By- (NH) or By+ (SH) for By=0
	IF (H*BY .GT. 0.0 .AND. BZ .GE. 0.0) NDX = 1
	IF (H*BY .GT. 0.0 .AND. BZ .LT. 0.0) NDX = 2
	IF (H*BY .LE. 0.0 .AND. BZ .GE. 0.0) NDX = 3
	IF (H*BY .LE. 0.0 .AND. BZ .LT. 0.0) NDX = 4
      ENDIF

C          Adjust sign for hemisphere, convert from mV/m -> V/m (with
C          factor 1.E-3) and maintain constant field equatorward of
C          minimum model magnetic latitude (~55.06 deg).
      IF (RMLA .GE. RMLAMN(MODL)) THEN
	XLA = RMLA
	FET = 1.E-3
	FEP = 1.E-3
      ELSE IF (RMLA .LE. -RMLAMN(MODL)) THEN
	XLA = -RMLA
	FET = 1.E-3
	FEP = 1.E-3
      ELSE
	XLA = RMLAMN(MODL)
	FET = 0.
	FEP = 1.E-3*COS(RMLAMN(MODL)*PI/180.)/SIN((90.-RMLAT)*PI/180.)
      ENDIF
      CALL GETFIT (NDX,MODL,TIME,XLA,DFIT,ESFIT,EEFIT,IST)
      IF (IST .NE. 0) WRITE (6,'(''MHEMODL: GETFIT returned IST='',I3)')
     +                IST

      EPOT = DFIT - ADDPOT(NDX,MODL)
      ET = ESFIT*FET
      EP = EEFIT*FEP

      RETURN
      END

      SUBROUTINE MHINIT (MODL,IUN,IPR,ISTAT)
C          Initialize Millstone Hill hemispheric electric field model.
C          Read 'FITVEC' output (Millstone Hill terminology) parameters
C          into common MHRFCM.  This is specifically for two models:
C          Millstone Hill electric potential patterns binned by nine
C          Hemispheric Power Input levels (MHI) and Millstone Hill/
C          Sondrestrom electric potential patterns binned by 4 IMF By/Bz
C          conditions (MHS).
C
C          INPUTS:
C            MODL = Model selection flag, where
C                 = 1 for MHI model
C                 = 2 for MHS model
C            IUN  = Logical unit number of FITVEC parameter file which has
C                   already been opened; e.g.,
C                     OPEN (IUN,FILE='mhi.cofcnts')   when MODL=1
C                     OPEN (IUN,FILE='mhs.cofcnts')   when MODL=2
C            IPR  = 1,0  if do,not print diagnostics during initialization.
C          RETURNS:
C            ISTAT  = 0 - Normal return.
C                    -2 - END-OF-FILE when more data are expected.
C                    -1 - Error reading parameter file.
C                     1 - Illegal KFIT
C                     2 - Inconsistency in parameters defining time nodes.
C                     3 - Inconsistency in parameters defining latitude node
C                     4 - Inconsistency in total number of parameters.
C                     5 - Too many parameters for BETA, F arrays. The array
C                         dimensions should be increased.
C
C          John M. Holt, Millstone Hill, jmh@chaos.haystack.edu - 1/86
C
C          Mod Oct 96 (Roy Barnes NCAR): omit reference to third MH model.
C          Removed reads here, first sets of coefficients from 'mhs.cofcnts',
C          and reduced last array dimension from 4 to 2 for TY and BETA.
C          Also, formal argument MODL was added and two remaining read loops
C          (for MHI and MHS model coefficients) were combined and fit limits
C          parameters (TMN-RMLAMX) were moved to MHRFCM.  Reset variables to
C          REAL*8 to conform with current Millstone Hill version and eliminate
C          underflow in BASPRC.  Also eliminated return of ion-drift components
C          which were based on a simple approximation (B =.5) and the
C          electostatic field components.

C          MHRFCM contains model coefficients and is assigned here:
      PARAMETER (MXNBETA=143,MXNNDX=9,MXNF=196)
      REAL PI,RE,TMN,TMX,RMLAMN,RMLAMX,TX,TY,BETA,F
      COMMON /MHRFCM/ NX,KX,NY,KY,NBETA,KFIT,PI,RE,TMN,TMX,RMLAMN(2),
     + RMLAMX(2),TX(100),TY(100,2),BETA(MXNBETA,MXNNDX,2),F(3,3,MXNF)
C             NX+KX  = No. assigned elements in TX
C             NY+KY  = No. assigned elements in TY   (1st dim)
C             NBETA  = No. assigned elements in BETA (1st dim)
C             KFIT   = Kind of fit (currently always 2)
C             RE     = Earth radius
C             TMN    = Minimum time of fit (hours).
C             TMX    = Maximum time of fit.
C             RMLAMN = Minimum magnetic latitude of fit (degrees).
C             RMLAMX = Maximum magnetic latitude of fit (degrees).
C             TX     = spline fit independent variable MLT (radians)
C             TY     = spline fit independent variable mag colat (km)
C             BETA   = fit parameter electric potential
C                      BETA 3rd dimension is MODL (1=MHI, 2=MHS)
C                      BETA 2nd dimension is hemispheric power input index (1-9)
C                      for MHI model or By/Bz sign (1-4) for MHS, where
C                               NDX  By    Bz
C                                1    +     +
C                                2    +     -
C                                3    -     +
C                                4    -     -
C                      (B_ is considered + or - when |B_| > .5 nT)

C          Local declarations
      CHARACTER CMODL*3, CBYZ*16, LABELS*80

C          Initialize common constants
      PI = 3.141592653589793123846
      RE = 6368.0

      IF (MODL .EQ. 1) THEN
	NSET  = 1
	CMODL ='MHI'
      ELSE IF (MODL .EQ. 2) THEN
	NSET  = 4
	CMODL ='MHS'
      ELSE
	WRITE (6,'(''MHINIT:  MODL may be 1 or 2, not'',I8)') MODL
	STOP
      ENDIF

      DO 100 ISET = 1,NSET
      CBYZ  = ' '
      IF (MODL .EQ. 2) THEN
	READ (IUN,'(9X,F3.1,5X,F3.1)') BY,BZ
	WRITE (CBYZ,'(2F8.1)') BY,BZ
      ENDIF

      READ (IUN,*) NLABS
      DO 10 I=1,NLABS
   10 READ (IUN,'(A)',END=998,ERR=999) LABELS

      READ (IUN,*,END=998,ERR=999) KFIT
      IF (KFIT .NE. 2) THEN
	ISTAT = 1
	RETURN
      ENDIF

      READ (IUN,*,END=998,ERR=999) NKTSX,NX,NBETAX,KX
      IF (NBETAX .NE. NX-(KX-1)) THEN
	ISTAT = 2
	RETURN
      ENDIF
      READ (IUN,*,END=998,ERR=999) (TX(I),I=1,NKTSX+2*(KX-1))

      READ (IUN,*,END=998,ERR=999) NKTSY,NY,NBETAY,KY
      IF (NBETAY .NE. NY-1) THEN
	ISTAT = 3
	RETURN
      ENDIF
      READ (IUN,*,END=998,ERR=999) (TY(I,MODL),I=1,NKTSY+2*(KY-1))

      IBEG = 1
      IEND = MXNNDX
      IF (MODL .EQ. 2) THEN
	IBEG = ISET
	IEND = ISET
      ENDIF
      DO NDX=IBEG,IEND
      READ (IUN,*,END=998,ERR=999) NBETA
      IF (NBETA .NE. NBETAX*NBETAY) THEN
	ISTAT = 4
	RETURN
      ELSE IF (NBETA .GT. MXNBETA) THEN
	ISTAT = 5
	RETURN
      ENDIF
      READ (IUN,*,END=998,ERR=999) (BETA(I,NDX,MODL),I=1,NBETA)

      IF (MODL .EQ. 2) THEN
C          Reverse the sign of BETA for MHS to conform to MHI convention
	DO 20 I=1,NBETA
   20   BETA(I,NDX,MODL) = -BETA(I,NDX,MODL)
      ENDIF
      ENDDO
C  30 IF (IPR .EQ. 1) WRITE (6,'(''MHINIT: Read ''A,'' pars: I MODL BY B
C    +Z BETA1 ='',2I3,A,E12.3)') CMODL,NDX,MODL,CBYZ,BETA(1,NDX,MODL)


C          Establish independent fit parameter (MLT and mag lat) model extremes
C          Model inputs should be filtered by these limits, even though
C          inputs are constrained in GETFIT.
      TMN    =  0.
      TMX    = 24.
      RMLAMN(MODL) = 90.0-TY(NY+1,MODL)*180./(PI*RE)
      RMLAMX(MODL) = 90.0-TY(1   ,MODL)*180./(PI*RE)

  100 CONTINUE

      ISTAT = 0
      RETURN

C          Bad file format traps
  998 ISTAT = -1
      RETURN
  999 ISTAT = -2
      RETURN
      END

      SUBROUTINE GETFIT (NDX,MODL,TIME,RMLA,DFIT,ESFIT,EEFIT,ISTAT)
C          GETFIT recovers the electrostatic potential and electric field
C          components from the FITVEC parameters stored in COMMON MHRFCM
C          for the point TIME,RMLA.    The potential is a B-splines fit and
C          the electic field is -grad(potential).
C
C          INPUTS:
C            NDX   = Index of potential pattern.  If MODL=1, NDX may be 1-9
C                    corresponding to Hemispheric Power Input Index.  If MODL=2,
C                    index may be 1-4 corresponding to the 4 By/Bz paterns.
C            MODL  = model to be used (1) MHI or (2) MHS
C            TIME  = Magnetic local time (hours); must be in domain TMN to TMX.
C            RMLA  = Magnetic latitude (degrees); must be in domain RMLAMN to
C                    RMLAMX.
C          RETURNS:
C            DFIT  = Electrostatic potential (kV)
C            ESFIT = Southward electric field (mV/m).
C            EEFIT = Eastward  electric field (mV/m).
C            ISTAT = Status:  0 = Normal return
C                             1 = Bad NDX,MODL input
C          Author: John M. Holt, Millstone Hill, jmh@chaos.haystack.edu - 1/86
C          Mod Oct 96 (Roy Barnes):  Reset variables to REAL*8 to conform with
C          current Millstone Hill version and eliminate (32-bit) underflow
C          when computing F in BASPRC.

C          Formal argument declarations
      REAL TIME,RMLA,DFIT,ESFIT,EEFIT

C          MHRFCM contains model coefficients and is assigned in MHINIT
      PARAMETER (MXNBETA=143,MXNNDX=9,MXNF=196)
      REAL PI,RE,TMN,TMX,RMLAMN,RMLAMX,TX,TY,BETA,F
      COMMON /MHRFCM/ NX,KX,NY,KY,NBETA,KFIT,PI,RE,TMN,TMX,RMLAMN(2),
     + RMLAMX(2),TX(100),TY(100,2),BETA(MXNBETA,MXNNDX,2),F(3,3,MXNF)
      DIMENSION XY(NY+KY)

C          Local declarations
      REAL X,Y,DDXFIT,DDYFIT
      SAVE

C          Check valid model indices
C          Note:  MXNDX ought to be added to MHRFCM as a 2 element array
      ISTAT = 1
      MXNDX = 0
      IF (MODL .EQ. 1) MXNDX = 9
      IF (MODL .EQ. 2) MXNDX = 4
      IF (MXNDX .EQ. 0) GO TO 100
      IF (NDX .LT. 1 .OR. NDX .GT. MXNDX) GO TO 100
      ISTAT  = 0

C          Convert hours MLT to radians and degrees magnetic latitude to ground
C          distance from the pole (km), which are the independent variables
C          in the spline fits.  Constrain the results to be slightly within
C          the maximum fit extreme (thus, avoiding a divide by zero (Y) below
C          and trouble when INTERV finds the interval for the biggest X or Y.
      X = 2.*PI*(1.0-TIME/24.0)
      Y = RE*PI*(90.0-RMLA)/180.0
      X = MIN(MAX(X,.00000001),TX(NX+1)-.00000001)
      Y = MIN(MAX(Y,.0000001),TY(NY+1,MODL)-.0000001)

C          Compute the fit
C     CALL BASPRC (TX,TY(1,MODL),NX,NY,KX,KY,X,Y,F)
      XY(1:NY+KY) = TY(1:NY+KY,MODL)
      CALL BASPRC (TX,XY,NX,NY,KX,KY,X,Y,F)

      DFIT = 0.0
      DDXFIT = 0.
      DDYFIT = 0.
      DO 30 I=1,NBETA
      DFIT   = DFIT   + BETA(I,NDX,MODL)*F(1,1,I)
      DDXFIT = DDXFIT + BETA(I,NDX,MODL)*F(2,1,I)/Y
   30 DDYFIT = DDYFIT + BETA(I,NDX,MODL)*F(1,2,I)
      DDYFIT = -DDYFIT
      DFIT   = -DFIT/1000.
      ESFIT  = -DDYFIT
      EEFIT  = -DDXFIT

  100 RETURN
      END

      SUBROUTINE BASPRC (TX,TY,NX,NY,KX,KY,X,Y,F)
C          BASPRC is a two-dimensional tensor product B-spline basis function
C          routine.  The splines generated by these basis functions are
C          periodic in the X direction. this version of BASPRC is also zero
C          at the first Y knot. this function is a useful representation of
C          the electrostatic potential determined from an IS radar azimuth scan.
C
C          Author: John M. Holt, Millstone Hill, jmh@chaos.haystack.edu - 8/83

C          Formal argument declarations
      REAL TX(*),TY(*),X,Y,F(9,*)

C          Local declarations
      REAL XP,YP,VALX(16),VALY(16),A(4,4)
      SAVE

      NDERIV = MIN0 (3,KX-1,KY-1)
      N = NX*NY
      DO 10 I=1,N
      DO 10 J=1,9
   10 F(J,I) = 0.

      XP = MOD (X,TX(NX+1))
      CALL INTERV (TX, NX+KX, XP, ILEFTX, MFLAG)
      YP = Y
      CALL INTERV (TY, NY+KY, YP, ILEFTY, MFLAG)
      IF (ILEFTX .LT. KX .OR. ILEFTX .GT. NX) THEN
C	WRITE (6,'(''BASPRC: out of range of spline:  ILEFTX KX NX X='',
C    +  3I3,E12.4)') ILEFTX,KX,NX,X
        ILEFTX = NX+KX
C	STOP
      ENDIF
      IF (ILEFTY .LT. KY .OR. ILEFTY .GT. NY) THEN
C	WRITE (6,'(''BASPRC: out of range of spline:  ILEFTY KY NY Y='',
C    +  3I3,E12.4)') ILEFTY,KY,NY,Y
        ILEFTY = NY+KY
C	STOP
      ENDIF

      LFTMKX = ILEFTX - KX
      CALL BSPLVD (TX, KX, XP, ILEFTX, A, VALX, NDERIV)
      LFTMKY = ILEFTY - KY
      CALL BSPLVD (TY, KY, YP, ILEFTY, A, VALY, NDERIV)
      NPX = NX - (KX-1)

      DO 20 MX=1,KX
      IX = LFTMKX + MX
      IF (IX .GE. NX-(KX-2)) IX = IX-NX+(KX-1)

      DO 30 MY=1,KY
      IY = LFTMKY + MY - 1
      IF (IY .GT. 0) THEN
	DO 40 JX=1,3
	LX = MX + KX*(JX-1)
	DO 50 JY=1,3
	J = JX + (JY-1)*3
	LY = MY + KY*(JY-1)
	I = IX + (IY-1)*NPX
   50   F(J,I) = VALX(LX)*VALY(LY)
   40   CONTINUE
      ENDIF
   30 CONTINUE
   20 CONTINUE

      RETURN
      END

      SUBROUTINE INTERV ( XT, LXT, X, LEFT, MFLAG )
C          Find the adjacent array elements in XT that bound X.
C           INPUTS:
C             XT  = Array of nondecreasing values
C             LXT = Number of assigned elements in XT
C             X   = Value whose location with respect to XT is to be determined.
C
C           RETURNS:
C             LEFT  = Index of XT, such that XT(LEFT) <= X < XT(LEFT+1) except
C                     when MFLAG is not 0.
C             MFLAG = Return status
C                   =  0  (normal)
C                   = -1  (when X <  XT(1)  ; in this case LEFT = 1)
C                   =  1  (when X >= XT(LXT); in this case LEFT = LXT)
C
C          Non-zero MFLAG indicates that X lies outside the halfopen interval
C          XT(1) <= Y < XT(LXT).  The asymmetric treatment of the interval
C          is due to the decision to make all PP functions continuous from
C          the right.
C
C          ALGORITHM:
C          From:  "A PRACTICAL GUIDE TO SPLINES"  by C. De Boor.
C          This is designed to be efficient in the common situation that it
C          is called repeatedly, with X taken from an increasing or decreasing
C          sequence.  This will happen, e.g., when a PP function is to be
C          graphed.  The first guess for LEFT is therefore taken to be the val-
C          ue returned at the previous call and stored in the  local  variable
C          ILO.  A first check ascertains that ILO .LT. LXT (this is necessary
C          since the present call may have nothing to do with the previous
C          call) then, if XT(ILO) .LE. X .LT. XT(ILO+1), we set  LEFT = ILO
C          and are done after just three comparisons.  Otherwise, we
C          repeatedly double the difference (ISTEP = IHI - ILO) while also
C          moving ILO and IHI in the direction of X until
C               XT(ILO) .LE. X .LT. XT(IHI)
C          after which we use bisection to get, in addition, ILO+1 = IHI.
C          LEFT = ILO is then returned.

C          Formal argument declarations
      REAL X,XT(LXT)

C          Local declarations
      DATA ILO /1/
      SAVE

      IHI = ILO + 1
      IF (IHI .LT.    LXT ) GO TO 20
      IF (X   .GE. XT(LXT)) GO TO 110
      IF (LXT .LE.      1 ) GO TO 90
      ILO = LXT - 1
      IHI = LXT

   20 IF (X .GE. XT(IHI)) GO TO 40
      IF (X .GE. XT(ILO)) GO TO 100

C          Now X < XT(ILO), decrease ILO to capture X
      ISTEP = 1
   31 IHI = ILO
      ILO = IHI - ISTEP
      IF (ILO .LE.      1 ) GO TO 35
      IF (X   .GE. XT(ILO)) GO TO 50
      ISTEP = ISTEP*2
      GO TO 31

   35 ILO = 1
      IF (X .LT. XT(1)) GO TO 90
      GO TO 50

C          Now X >= XT(IHI), increase  IHI to capture X
   40 ISTEP = 1
   41 ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI .GE.    LXT ) GO TO 45
      IF (X   .LT. XT(IHI)) GO TO 50
      ISTEP = ISTEP*2
      GO TO 41

   45 IF (X .GE. XT(LXT)) GO TO 110
      IHI = LXT

C          Now XT(ILO) <= X < XT(IHI), narrow the interval
   50 MIDDLE = (ILO + IHI)/2
      IF (MIDDLE .EQ. ILO) GO TO 100
C          Note: It is assumed that MIDDLE = ILO in case IHI = ILO+1
      IF (X .LT. XT(MIDDLE)) GO TO 53
      ILO = MIDDLE
      GO TO 50
   53 IHI = MIDDLE
      GO TO 50

C          Set output and return.
   90 MFLAG = -1
      LEFT = 1
      RETURN

  100 MFLAG = 0
      LEFT = ILO
      RETURN

  110 MFLAG = 1
      LEFT = LXT
      RETURN
      END

      SUBROUTINE BSPLVD ( T, K, X, LEFT, A, DBIATX, NDERIV )
C  FROM  * A PRACTICAL GUIDE TO SPLINES *  BY C. DE BOOR
CALLS BSPLVB
CALCULATES VALUE AND DERIV.S OF ALL B-SPLINES WHICH DO NOT VANISH AT X
C
C******  I N P U T  ******
C  T     THE KNOT ARRAY, OF LENGTH LEFT+K (AT LEAST)
C  K     THE ORDER OF THE B-SPLINES TO BE EVALUATED
C  X     THE POINT AT WHICH THESE VALUES ARE SOUGHT
C  LEFT  AN INTEGER INDICATING THE LEFT ENDPOINT OF THE INTERVAL OF
C        INTEREST. THE  K  B-SPLINES WHOSE SUPPORT CONTAINS THE INTERVAL
C               (T(LEFT), T(LEFT+1))
C        ARE TO BE CONSIDERED.
C  A S S U M P T I O N  - - -  IT IS ASSUMED THAT
C               T(LEFT) .LT. T(LEFT+1)
C        DIVISION BY ZERO WILL RESULT OTHERWISE (IN  B S P L V B ).
C        ALSO, THE OUTPUT IS AS ADVERTISED ONLY IF
C               T(LEFT) .LE. X .LE. T(LEFT+1) .
C  NDERIV   AN INTEGER INDICATING THAT VALUES OF B-SPLINES AND THEIR
C        DERIVATIVES UP TO BUT NOT INCLUDING THE  NDERIV-TH  ARE ASKED
C        FOR. ( NDERIV  IS REPLACED INTERNALLY BY THE INTEGER  M H I G H
C        IN  (1,K)  CLOSEST TO IT.)
C
C******  W O R K   A R E A  ******
C  A     AN ARRAY OF ORDER (K,K), TO CONTAIN B-COEFF.S OF THE DERIVAT-
C        IVES OF A CERTAIN ORDER OF THE  K  B-SPLINES OF INTEREST.
C
C******  O U T P U T  ******
C  DBIATX   AN ARRAY OF ORDER (K,NDERIV). ITS ENTRY  (I,M)  CONTAINS
C        VALUE OF  (M-1)ST  DERIVATIVE OF  (LEFT-K+I)-TH  B-SPLINE OF
C        ORDER  K  FOR KNOT SEQUENCE  T , I=M,...,K, M=1,...,NDERIV.
C
C******  M E T H O D  ******
C  FROM  * A PRACTICAL GUIDE TO SPLINES *  BY C. DE BOOR
C  VALUES AT  X  OF ALL THE RELEVANT B-SPLINES OF ORDER K,K-1,...,
C  K+1-NDERIV  ARE GENERATED VIA  BSPLVB  AND STORED TEMPORARILY IN
C  DBIATX .  THEN, THE B-COEFFS OF THE REQUIRED DERIVATIVES OF THE B-
C  SPLINES OF INTEREST ARE GENERATED BY DIFFERENCING, EACH FROM THE PRE-
C  CEDING ONE OF LOWER ORDER, AND COMBINED WITH THE VALUES OF B-SPLINES
C  OF CORRESPONDING ORDER IN  DBIATX  TO PRODUCE THE DESIRED VALUES .
C
      SAVE
      INTEGER K,LEFT,NDERIV,   I,IDERIV,IL,J,JLOW,JP1MID,KP1,KP1MM
     *                        ,LDUMMY,M,MHIGH

      REAL A(K,K),DBIATX(K,NDERIV),T(LEFT+K),X,FACTOR,FKP1MM,SUM

      MHIGH = MAX0(MIN0(NDERIV,K),1)
C     MHIGH IS USUALLY EQUAL TO NDERIV.
      KP1 = K+1
      CALL BSPLVB(T,KP1-MHIGH,1,X,LEFT,DBIATX)
      IF (MHIGH .EQ. 1)                 GO TO 99
C     THE FIRST COLUMN OF  DBIATX  ALWAYS CONTAINS THE B-SPLINE VALUES
C     FOR THE CURRENT ORDER. THESE ARE STORED IN COLUMN K+1-CURRENT
C     ORDER  BEFORE  BSPLVB  IS CALLED TO PUT VALUES FOR THE NEXT
C     HIGHER ORDER ON TOP OF IT.
      IDERIV = MHIGH
      DO 15 M=2,MHIGH
         JP1MID = 1
         DO 11 J=IDERIV,K
            DBIATX(J,IDERIV) = DBIATX(JP1MID,1)
   11       JP1MID = JP1MID + 1
         IDERIV = IDERIV - 1
         CALL BSPLVB(T,KP1-IDERIV,2,X,LEFT,DBIATX)
   15    CONTINUE
C
C     AT THIS POINT,  B(LEFT-K+I, K+1-J)(X) IS IN  DBIATX(I,J) FOR
C     I=J,...,K AND J=1,...,MHIGH ('=' NDERIV). IN PARTICULAR, THE
C     FIRST COLUMN OF  DBIATX  IS ALREADY IN FINAL FORM. TO OBTAIN COR-
C     RESPONDING DERIVATIVES OF B-SPLINES IN SUBSEQUENT COLUMNS, GENE-
C     RATE THEIR B-REPR. BY DIFFERENCING, THEN EVALUATE AT  X.
C
      JLOW = 1
      DO 20 I=1,K
         DO 19 J=JLOW,K
   19       A(J,I) = 0.
         JLOW = I
   20    A(I,I) = 1.
C     AT THIS POINT, A(.,J) CONTAINS THE B-COEFFS FOR THE J-TH OF THE
C     K  B-SPLINES OF INTEREST HERE.
C
      DO 40 M=2,MHIGH
         KP1MM = KP1 - M
         FKP1MM = FLOAT(KP1MM)
         IL = LEFT
         I = K
C
C        FOR J=1,...,K, CONSTRUCT B-COEFFS OF  (M-1)ST  DERIVATIVE OF
C        B-SPLINES FROM THOSE FOR PRECEDING DERIVATIVE BY DIFFERENCING
C        AND STORE AGAIN IN  A(.,J) . THE FACT THAT  A(I,J) = 0  FOR
C        I .LT. J  IS USED.
         DO 25 LDUMMY=1,KP1MM
            FACTOR = FKP1MM/(T(IL+KP1MM) - T(IL))
C           THE ASSUMPTION THAT T(LEFT).LT.T(LEFT+1) MAKES DENOMINATOR
C           IN  FACTOR  NONZERO.
            DO 24 J=1,I
   24          A(I,J) = (A(I,J) - A(I-1,J))*FACTOR
            IL = IL - 1
   25       I = I - 1
C
C        FOR I=1,...,K, COMBINE B-COEFFS A(.,I) WITH B-SPLINE VALUES
C        STORED IN DBIATX(.,M) TO GET VALUE OF  (M-1)ST  DERIVATIVE OF
C        I-TH B-SPLINE (OF INTEREST HERE) AT  X , AND STORE IN
C        DBIATX(I,M). STORAGE OF THIS VALUE OVER THE VALUE OF A B-SPLINE
C        OF ORDER M THERE IS SAFE SINCE THE REMAINING B-SPLINE DERIVAT-
C        IVES OF THE SAME ORDER DO NOT USE THIS VALUE DUE TO THE FACT
C        THAT  A(J,I) = 0  FOR J .LT. I .
   30    DO 40 I=1,K
            SUM = 0.
            JLOW = MAX0(I,M)
            DO 35 J=JLOW,K
   35          SUM = A(J,I)*DBIATX(J,M) + SUM
   40       DBIATX(I,M) = SUM
   99                                   RETURN
      END

      SUBROUTINE BSPLVB ( T, JHIGH, INDEX, X, LEFT, BIATX )
C  FROM  * A PRACTICAL GUIDE TO SPLINES *  BY C. DE BOOR
CALCULATES THE VALUE OF ALL POSSIBLY NONZERO B-SPLINES AT  X  OF ORDER
C
C               JOUT  =  MAX( JHIGH , (J+1)*(INDEX-1) )
C
C  WITH KNOT SEQUENCE  T .
C
C******  I N P U T  ******
C  T.....KNOT SEQUENCE, OF LENGTH  LEFT + JOUT  , ASSUMED TO BE NONDE-
C        CREASING.  A S S U M P T I O N . . . .
C                       T(LEFT)  .LT.  T(LEFT + 1)   .
C   D I V I S I O N  B Y  Z E R O  WILL RESULT IF  T(LEFT) = T(LEFT+1)
C  JHIGH,
C  INDEX.....INTEGERS WHICH DETERMINE THE ORDER  JOUT = MAX(JHIGH,
C        (J+1)*(INDEX-1))  OF THE B-SPLINES WHOSE VALUES AT  X  ARE TO
C        BE RETURNED.  INDEX  IS USED TO AVOID RECALCULATIONS WHEN SEVE-
C        RAL COLUMNS OF THE TRIANGULAR ARRAY OF B-SPLINE VALUES ARE NEE-
C        DED (E.G., IN  BVALUE  OR IN  BSPLVD ). PRECISELY,
C                     IF  INDEX = 1 ,
C        THE CALCULATION STARTS FROM SCRATCH AND THE ENTIRE TRIANGULAR
C        ARRAY OF B-SPLINE VALUES OF ORDERS 1,2,...,JHIGH  IS GENERATED
C        ORDER BY ORDER , I.E., COLUMN BY COLUMN .
C                     IF  INDEX = 2 ,
C        ONLY THE B-SPLINE VALUES OF ORDER  J+1, J+2, ..., JOUT  ARE GE-
C        NERATED, THE ASSUMPTION BEING THAT  BIATX , J , DELTAL , DELTAR
C        ARE, ON ENTRY, AS THEY WERE ON EXIT AT THE PREVIOUS CALL.
C           IN PARTICULAR, IF  JHIGH = 0, THEN  JOUT = J+1, I.E., JUST
C        THE NEXT COLUMN OF B-SPLINE VALUES IS GENERATED.
C
C  W A R N I N G . . .  THE RESTRICTION   JOUT .LE. JMAX (= 20)  IS IM-
C        POSED ARBITRARILY BY THE DIMENSION STATEMENT FOR  DELTAL  AND
C        DELTAR  BELOW, BUT IS  N O W H E R E  C H E C K E D  FOR .
C
C  X.....THE POINT AT WHICH THE B-SPLINES ARE TO BE EVALUATED.
C  LEFT.....AN INTEGER CHOSEN (USUALLY) SO THAT
C                  T(LEFT) .LE. X .LE. T(LEFT+1)  .
C
C******  O U T P U T  ******
C  BIATX.....ARRAY OF LENGTH  JOUT , WITH  BIATX(I)  CONTAINING THE VAL-
C        UE AT  X  OF THE POLYNOMIAL OF ORDER  JOUT  WHICH AGREES WITH
C        THE B-SPLINE  B(LEFT-JOUT+I,JOUT,T)  ON THE INTERVAL (T(LEFT),
C        T(LEFT+1)) .
C
C******  M E T H O D  ******
C  THE RECURRENCE RELATION
C
C                       X - T(I)              T(I+J+1) - X
C     B(I,J+1)(X)  =  -----------B(I,J)(X) + ---------------B(I+1,J)(X)
C                     T(I+J)-T(I)            T(I+J+1)-T(I+1)
C
C  IS USED (REPEATEDLY) TO GENERATE THE (J+1)-VECTOR  B(LEFT-J,J+1)(X),
C  ...,B(LEFT,J+1)(X)  FROM THE J-VECTOR  B(LEFT-J+1,J)(X),...,
C  B(LEFT,J)(X), STORING THE NEW VALUES IN  BIATX  OVER THE OLD. THE
C  FACTS THAT
C            B(I,1) = 1  IF  T(I) .LE. X .LT. T(I+1)
C  AND THAT
C            B(I,J)(X) = 0  UNLESS  T(I) .LE. X .LT. T(I+J)
C  ARE USED. THE PARTICULAR ORGANIZATION OF THE CALCULATIONS FOLLOWS AL-
C  GORITHM  (8)  IN CHAPTER X OF THE TEXT.
C
      PARAMETER (JMAX=20)
      REAL BIATX(JHIGH),T(LEFT+JHIGH),X,DELTAL(JMAX),DELTAR(JMAX),
     +       SAVED,TERM

C     DIMENSION BIATX(JOUT), T(LEFT+JOUT)
CURRENT FORTRAN STANDARD MAKES IT IMPOSSIBLE TO SPECIFY THE LENGTH OF
C  T  AND OF  BIATX  PRECISELY WITHOUT THE INTRODUCTION OF OTHERWISE
C  SUPERFLUOUS ADDITIONAL ARGUMENTS.
      DATA J/1/
C     SAVE J,DELTAL,DELTAR
      SAVE
C
                                        GO TO (10,20), INDEX
   10 J = 1
      BIATX(1) = 1.
      IF (J .GE. JHIGH)                 GO TO 99
C
   20    JP1 = J + 1
         DELTAR(J) = T(LEFT+J) - X
         DELTAL(J) = X - T(LEFT+1-J)
         SAVED = 0.
         DO 26 I=1,J
            TERM = BIATX(I)/(DELTAR(I) + DELTAL(JP1-I))
            BIATX(I) = SAVED + DELTAR(I)*TERM
   26       SAVED = DELTAL(JP1-I)*TERM
         BIATX(JP1) = SAVED
         J = JP1
         IF (J .LT. JHIGH)              GO TO 20
C
   99                                   RETURN
      END

