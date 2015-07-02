CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HMREPOT (RMLAT,RMLT,BY,BZ,RKP,IFRST,ETHETA,EPHI,
     | POTKV)
      SAVE

C HMREPOT computes the electric potential in kV poleward of +50/-50 
C  degrees magnetic latitude from the Heppner/Maynard/Rich model
C  described in Rich and Maynard (JGR, 94, 3687-3701, 1989) from patterns
C  described by Heppner and Maynard (JGR, 92, 4467-4489, 1987).
C The original patterns are for the Northern Hemisphere and are found at
C  a particular magnetic latitude and longitude (0=0000MLT) in HMRMOD.
C  The patterns are only valid for magnetic latitudes between 50 and 90
C  degrees (CMIN and CMAX in GETHMR).  However, the Southern Hemisphere
C  can be found by assuming an opposite sign for By.  The patterns are 
C  based on the iabc model index:
C              IABC=1 IS FOR MODEL A   (IMF BZ < 0, BY < 0, NORTHERN HEM
C              IABC=2 IS FOR MODEL BC  (IMF BZ < 0, BY > 0, NORTHERN HEM
C              IABC=3 IS FOR MODEL DE  (IMF BZ < 0, BY < 0, NORTHERN HEM
C              IABC=4 IS FOR MODEL BCP (IMF BZ > 0, BY > 0, NORTHERN HEM
C              IABC=5 IS FOR MODEL BCPP(IMF BZ >>0, BY > 0, NORTHERN HEM
C              IABC=6 IS FOR MODEL DEP (IMF BZ > 0, BY < 0, NORTHERN HEM
C              IABC=7 IS FOR MODEL DEPP(IMF BZ >>0, BY < 0, NORTHERN HEM
C    For iabc=3, model DE is for more negative By than model A (iabc=1),
C       and is interpreted as By<=-5., while A is 0<=By<5.
C    For iabc=7, model DEPP, Bz >> 0 is interpreted as Bz>5.
C    If Bz=By=0, the model chosen is A (iabc=1)
C    For Bz .le. 0 models (iabc=1-3), the coefficients are adjusted
C       in the subroutine SETHM3 by the value of Kp.  The default
C       coefficients are for Kp=3.5 (between 3+ and 4-), and the
C       maximum Kp used is 7.0.
C  A particular model is chosen only if IFIRST=1 (which is then changed).
C  In addition to finding the electric potential, HMREPOT calculates the
C  electric field in V/m by finding the difference at 300 km in the
C  potential at nearby points 10 km away (STEP in /EHMR) along the same
C  magnetic latitude or longitude.
C HMREPOT is used in the AMIE (Assimilative Mapping of Ionospheric 
C  Electrodynamics) procedure to define the statistical electric 
C  potential if this model is chosen as the initial condition.
C HMREPOT solves for the E fields using the space centered technique.
C  Is based on the routine FLOWV2.  (B. Emery, 10/89)
C
C RMLAT,RMLTLON = mag lat and lon (0 midnight) in degrees
C BY,BZ         = GSM IMF Y and Z components in nT (determine model used)
C RKP           = Kp index, used if Bz is negative or zero.
C ETHETA,EPHI   = E field components in V/m.  (THETA=+S, PHI=+E)
C  NOTE:  SH uses NH RMLAT, but changes sign of By, so +S is +equatorward)
C POTKV         = electrostatic potential in kV.
C  We use the formula:  E_field = -del(electrostatic potential) to find E
C
C    /EHMR/
C      ALTE = 6671 KM = MEAN EARTH RADIUS (6371 KM) + 300 KM = ALTITUDE 
C             AT WHICH E FIELDS ARE COMPUTED FROM CENTER OF EARTH IN KM 
C             I.E. COMPUTE E FIELDS AT 300 KM ABOVE THE EARTH.
C      STEP,... - VARIOUS STEPS (OFTEN 10 KM) AWAY FROM THE LOCATION 
C             USED IN ORDER TO FIND THE ELECTRIC FIELDS IN VOLTS/M 
C             USING THE SPACE CENTERED TECHNIQUE
      COMMON /EHMR/ ALTE,STEP,STF2,STF3,STF4,STF5,STF6
      COMMON /RUNCON/ PII,RAD,RE
      COMMON/CHGCOM/IABC,ICHGHM,CHGV,AHM(3),BHM(3),DXHM(3),DYHM(3),
     |   ACHG(3),BCHG(3),DXCHG(3),DYCHG(3)
      COMMON/SCOM/CMIN,CMAX

      IF (IFRST .EQ. 1) THEN
C  Choose model
      IHSOLN = 1.
      IF (RMLAT .LT. 0.) IHSOLN = -1.
      IF (BZ .LE. 0. .AND. BY*IHSOLN .LE. 0.) IABC = 1
      IF (BZ .LE. 0. .AND. BY*IHSOLN .GT. 0.) IABC = 2
      IF (BZ .LE. 0. .AND. BY*IHSOLN .LE. -5.) IABC = 3
      IF (BZ .GT. 0. .AND. BY*IHSOLN .GE. 0.) IABC = 4
      IF (BZ .GT. 5. .AND. BY*IHSOLN .GE. 0.) IABC = 5
      IF (BZ .GT. 0. .AND. BY*IHSOLN .LT. 0.) IABC = 6
      IF (BZ .GT. 5. .AND. BY*IHSOLN .LT. 0.) IABC = 7
      ICHGHM = 1
      CHGV = 1.
      IF (BZ .LE. 0.) CALL SETHM3 (RKP)
      IFRST = IFRST + 1
      ENDIF
C
C  Keep potential constant equatorward of CMIN (50 degrees)
      DLAT = AMAX1(CMIN,ABS(RMLAT))/RAD
C  DLAT = magnetic latitude in radians
c     DLON = RMLTLON/RAD
      DLON = RMLT*15./RAD
C  Compute potential POTEN in V
      CALL HMRMOD (DLAT,DLON,1.0,POTEN)
      POTKV = POTEN * 1.E-3
C
      if (dlat+stf3 .lt. pii/2.) go to 100
C  Special computation near the pole
C  Compute ETHETA (+S, or POTEN(N-S)) at pole
      RATIO = STF4
      THETA = DLAT - STF3
      CALL HMRMOD (THETA,DLON,RATIO,PHB)
      PHI = DLON + 180./RAD
      CALL HMRMOD (THETA,PHI,RATIO,PHA)
      ETHETA = (PHA - PHB) * STF2
C  Compute EPHI (+E, or POTEN(W-E)) at pole
      THETA = DLAT - STF3
      CC = COS(THETA)
      R = SQRT(STF5 + CC**2*STF6)
      RS = ALTE*SIN(THETA)
      ALAM = ATAN(RS/R)
      R = SQRT(RS**2 + R**2)
      RATIO = ALTE / R
      PHI = DLON + ATAN(STEP/(CC*ALTE))
      CALL HMRMOD (ALAM,PHI,RATIO,PHA)
      PHI = DLON - ATAN(STEP/(CC*ALTE))
      CALL HMRMOD (ALAM,PHI,RATIO,PHB)
      EPHI1 = (PHB - PHA) * STF2
C  Go 180 degrees across the pole
      PHI = DLON + 180./RAD + ATAN(STEP/(CC*ALTE))
      CALL HMRMOD (ALAM,PHI,RATIO,PHA)
      PHI = DLON + 180./RAD - ATAN(STEP/(CC*ALTE))
      CALL HMRMOD (ALAM,PHI,RATIO,PHB)
      EPHI2 = (PHB - PHA) * STF2
C  Average by subtraction since sign should change over the pole
      EPHI = (EPHI1 - EPHI2)/2.

      return

100   continue

C  Compute EPHI (+E, or POTEN(W-E)) at DLAT
      CC = COS(DLAT)
      R = SQRT(STF5 + CC**2*STF6)
      RS = ALTE*SIN(DLAT)
      ALAM = ATAN(RS/R)
      R = SQRT(RS**2 + R**2)
      RATIO = ALTE / R
      PHI = DLON + ATAN(STEP/(CC*ALTE))
      CALL HMRMOD (ALAM,PHI,RATIO,PHA)
      PHI = DLON - ATAN(STEP/(CC*ALTE))
      CALL HMRMOD (ALAM,PHI,RATIO,PHB)
      EPHI = (PHB - PHA) * STF2
C
C  Compute ETHETA (+S, or POTEN(N-S)) at DLON
      THETA = DLAT + STF3
      RATIO = STF4
      CALL HMRMOD (THETA,DLON,RATIO,PHA)
      THETA = DLAT - STF3
      CALL HMRMOD (THETA,DLON,RATIO,PHB)
      ETHETA = (PHA - PHB) * STF2
C
      if (ABS(RMLAT) .lt. CMIN) then
C  Special computation for E fields equatorward of CMIN
      ETHETA = 0.
      EPHI = EPHI * COS(CMIN/RAD)/SIN((90.-RMLAT)/RAD)
      endif

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GETHMR (NTAPE)
      SAVE
C
C  GETHMR gets the coefficients for the HMR (Heppner-Maynard-Rich) model
C   and sets initial values for some parameters.
C  The coefficient file is hmr89.cofcnts
C
C  INPUTS - CALLING SEQUENCE
C     NTAPE  - UNIT NUMBER FOR DATA FILE OF COEFFICIENTS
C
C   COMMON BLOCKS USED
C     /PCOM/
C       COEF - ARRAY OF COEFFICIENTS TO LEGENDRE POLYNOMIAL TO BE READ I
C       P    - ARRAY OF LEGENDRE POLYNOMIAL TO BE CALCULATED
C       NMAX - MAXIMUM ORDER OF COEFFICIENTS
C       PKV51 - DIURNAL AVE OF POTENTIAL IN KV AT 51 DEG FOR KP=3.5
C               (IS SUBTRACTED AND IS ACCURATE TO WITHIN 1.0 KV FOR
C                KP=7 AND WITHIN 0.3 KV FOR KP=0 IF PKV51*CHGV)
C     /SCOM/
C       CMIN - MINIMUM LATITUDE USED FOR FITTING
C       CMAX - MAXIMUM LATITUDE USED FOR FITTING
C     /RUNCON/
C       PII - VALUE OF PI
C       RAD - VALUE OF RATIO OF DEG/RADIAN
C       RE  - RADIUS OF EARTH IN M - RE MEAN = 6371
C     /EHMR/
C       ALTE = 6671 KM = MEAN EARTH RADIUS (6371 KM) + 300 KM = ALTITUDE 
C              AT WHICH E FIELDS ARE COMPUTED FROM CENTER OF EARTH IN KM 
C              I.E. COMPUTE E FIELDS AT 300 KM ABOVE THE EARTH.
C       STEP,... - VARIOUS STEPS (OFTEN 10 KM) AWAY FROM THE LOCATION 
C              USED IN ORDER TO FIND THE ELECTRIC FIELDS IN VOLTS/M 
C              USING THE SPACE CENTERED TECHNIQUE
C
C  INPUTS - DATA FILE
C     NNN
C     MMM
C     I      - INDEX OF COEFFICIENT
C     J      - INDEX OF COEFFICIENT
C     CCF    - COEFFICIENT OF SPHERICAL HARMONIC EXPANSION
C
      COMMON/PCOM/COEF(18,18,7),P(18,18),NNMAX(7),NMAX,PKV51(7)
     |,XCO(18,18),DP(18,18),CONST(18,18),SP(18),CP(18),FN(18),FM(18)
      COMMON/SCOM/CMIN,CMAX
      COMMON /RUNCON/ PII,RAD,RE
      COMMON /EHMR/ ALTE,STEP,STF2,STF3,STF4,STF5,STF6
C
      CHARACTER*10 MLBL(5)
C
      CMIN=50. 
      CMAX=90.
      NNMAX(1:7) = (/12,12,12,13,13,13,13/)
      PKV51(1:7) = (/ 2.012, 23.536, -13.525, 
     |       0.064, 0.074, -0.176, -0.268/)
C   CONSTANTS FOR LEGENDRE POLYNOMIAL
c     DATA CONST/324*0./
      XCO(1,1:18)=(/
     | .282095E+00, .488603E+00, .109255E+01, .228523E+01,
     | .468333E+01, .951188E+01, .192265E+02, .387523E+02,
     | .779645E+02, .156658E+03, .314501E+03, .630964E+03,
     | .126523E+04, .253611E+04, .508196E+04, .101809E+05,
     |   .203918E+05,   .408366E+05/)
      XCO(2,1:18)=(/
     | .488603E+00, .488603E+00, .546274E+00, .144531E+01,
     | .331161E+01, .719031E+01, .151999E+02, .316411E+02,
     | .652298E+02, .133599E+03, .272366E+03, .553392E+03,
     | .112151E+04, .226837E+04, .458082E+04, .923904E+04,
     |   .186151E+05,   .374743E+05/)
      XCO( 3,1:18)=(/
     | .946175E+00, .109255E+01, .546274E+00, .590044E+00,
     | .177013E+01, .440314E+01, .101333E+02, .223736E+02,
     | .481754E+02, .102038E+03, .213661E+03, .443701E+03,
     | .915709E+03, .188083E+04, .384866E+04, .785168E+04,
     |   .159791E+05,   .324537E+05/)
      XCO( 4,1:18)=(/
     | .186588E+01, .228523E+01, .144531E+01, .590044E+00,
     | .625836E+00, .207566E+01, .555021E+01, .134918E+02,
     | .310971E+02, .693209E+02, .151081E+03, .324033E+03,
     | .686782E+03, .144253E+04, .300864E+04, .623988E+04,
     |   .128827E+05,   .264983E+05/)
      XCO( 5,1:18)=(/
     | .370249E+01, .468333E+01, .331161E+01, .177013E+01,
     | .625836E+00, .656382E+00, .236662E+01, .674590E+01,
     | .172496E+02, .414272E+02, .955522E+02, .214328E+03,
     | .471128E+03, .102002E+04, .218269E+04, .462762E+04,
     |   .973844E+04,   .203694E+05/)
      XCO( 6,1:18)=(/
     | .736787E+01, .951188E+01, .719031E+01, .440314E+01,
     | .207566E+01, .656382E+00, .683184E+00, .264596E+01,
     | .798499E+01, .213929E+02, .534153E+02, .127330E+03,
     | .293800E+03, .661878E+03, .146420E+04, .319336E+04,
     |   .688612E+04,   .147131E+05/)
      XCO( 7,1:18)=(/
     | .146845E+02, .192265E+02, .151999E+02, .101333E+02,
     | .555021E+01, .236662E+01, .683184E+00, .707163E+00,
     | .291571E+01, .926339E+01, .259102E+02, .671087E+02,
     | .165101E+03, .391572E+03, .903721E+03, .204248E+04,
     |   .454057E+04,   .996084E+04/)
      XCO( 8,1:18)=(/
     | .292940E+02, .387523E+02, .316411E+02, .223736E+02,
     | .134918E+02, .674590E+01, .264596E+01, .707163E+00,
     | .728927E+00, .317732E+01, .105778E+02, .307916E+02,
     | .825507E+02, .209304E+03, .509767E+03, .120459E+04,
     |   .278052E+04,   .629979E+04/)
      XCO( 9,1:18)=(/
     | .584734E+02, .779645E+02, .652298E+02, .481754E+02,
     | .310971E+02, .172496E+02, .798499E+01, .291571E+01,
     | .728927E+00, .748901E+00, .343190E+01, .119255E+02,
     | .360281E+02, .997819E+02, .260366E+03, .650553E+03,
     |   .157290E+04,   .370647E+04/)
      XCO(10,1:18)=(/
     | .116766E+03, .156658E+03, .133599E+03, .102038E+03,
     | .693209E+02, .414272E+02, .213929E+02, .926339E+01,
     | .317732E+01, .748901E+00, .767395E+00, .368030E+01,
     | .133043E+02, .416119E+02, .118840E+03, .318704E+03,
     |   .816138E+03,   .201755E+04/)
      XCO(11,1:18)=(/
     | .233240E+03, .314501E+03, .272366E+03, .213661E+03,
     | .151081E+03, .955522E+02, .534153E+02, .259102E+02,
     | .105778E+02, .343190E+01, .767395E+00, .784642E+00,
     | .392321E+01, .147120E+02, .475361E+02, .139761E+03,
     |   .384731E+03,   .100877E+04/)
      XCO(12,1:18)=(/
     | .465998E+03, .630964E+03, .553392E+03, .443701E+03,
     | .324033E+03, .214328E+03, .127330E+03, .671087E+02,
     | .307916E+02, .119255E+02, .368030E+01, .784642E+00,
     | .800822E+00, .416119E+01, .161472E+02, .537941E+02,
     |   .162579E+03,   .458849E+03/)
      XCO(13,1:18)=(/
     | .931187E+03, .126523E+04, .112151E+04, .915709E+03,
     | .686782E+03, .471128E+03, .293800E+03, .165101E+03,
     | .825507E+02, .360281E+02, .133043E+02, .392321E+01,
     | .800822E+00, .816077E+00, .439471E+01, .176082E+02,
     |   .603802E+02,   .187325E+03/)
      XCO(14,1:18)=(/
     | .186100E+04, .253611E+04, .226837E+04, .188083E+04,
     | .144253E+04, .102002E+04, .661878E+03, .391572E+03,
     | .209304E+03, .997819E+02, .416119E+02, .147120E+02,
     | .416119E+01, .816077E+00, .830522E+00, .462415E+01,
     |   .190939E+02,   .672889E+02/)
      XCO(15,1:18)=(/
     | .371962E+04, .508196E+04, .458082E+04, .384866E+04,
     | .300864E+04, .218269E+04, .146420E+04, .903721E+03,
     | .509767E+03, .260366E+03, .118840E+03, .475361E+02,
     | .161472E+02, .439471E+01, .830522E+00, .844251E+00,
     |   .484985E+01,   .206029E+02/)
      XCO(16,1:18)=(/
     | .743510E+04, .101809E+05, .923904E+04, .785168E+04,
     | .623988E+04, .462762E+04, .319336E+04, .204248E+04,
     | .120459E+04, .650553E+03, .318704E+03, .139761E+03,
     | .537941E+02, .176082E+02, .462415E+01, .844251E+00,
     |   .857341E+00,   .507210E+01/)
      XCO(17,1:18)=(/
     | .148629E+05, .203918E+05, .186151E+05, .159791E+05,
     | .128827E+05, .973844E+04, .688612E+04, .454057E+04,
     | .278052E+04, .157290E+04, .816138E+03, .384731E+03,
     | .162579E+03, .603802E+02, .190939E+02, .484985E+01,
     |   .857341E+00,   .869857E+00/)
      XCO(18,1:18)=(/
     | .297130E+05, .408366E+05, .374743E+05, .324537E+05,
     | .264983E+05, .203694E+05, .147131E+05, .996084E+04,
     | .629979E+04, .370647E+04, .201755E+04, .100877E+04,
     | .458849E+03, .187325E+03, .672889E+02, .206029E+02,
     | .507210E+01,   .869857E+00/)
C
C   INITIALIZE ARRAY OF COEFFICIENTS FOR FIT FROM UNIT NTAPE
       P(:,:) = 0.
       CONST(:,:) = 0.

       REWIND NTAPE
       DO 80 IA=1,7
        READ(NTAPE,400,END=   82)MLBL
  400   FORMAT(5A10)
   82   READ(NTAPE,*,END=10000) NNN,MMM,I,J,CCF
10000   IF(NNN.EQ.-1) GO TO 81
        COEF(I,J,IA)=CCF
        GO TO 82
   81  CONTINUE

C
      WRITE(6,"(1X,5A10)") MLBL
      WRITE(6,"(1X,'USING MODEL COEFFICIENTS UP TO NMAX=',I2)")NNMAX(IA)
   80 CONTINUE
      WRITE(6,"(1X,'MODEL VALID FOR LATITUDES ',F5.1,'  TO  ',F5.1)")
     | CMIN,CMAX
C
C  Set constants in RUNCON and EHMR
      RE = 6.49E6
      PII = 4 * ATAN(1.)
      RAD = 180. / PII
      ALTE = 6671.
      STEP = 10.
      STF2 = 1.0 / (2.0 * STEP * 1000.)
      STF3 = ATAN(STEP/ALTE)
      STF4 = ALTE / SQRT(STEP**2+ALTE**2)
      STF5 = STEP**2
      STF6 = ALTE**2
C  Set low order coefficients
      DP(1,1)=0.0
      SP(1)=0.0
      CP(1)=1.0
      DO 2 N=2,18
      FN(N)=N
      DO 2 M=1,N
      FM(M)=M-1
    2 CONST(N,M)=FLOAT((N-2)**2-(M-1)**2)/
     |   FLOAT((2*N-3)*(2*N-5))
C
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SETHM3 (RKP)
      SAVE

C RKP = The real value for Kp, the log planetary magnetic index between 
C       0. and 9. Can be any real between these values, but the maximum 
C       value used is 7.
C.......................................................................
C                                                                      .
C            SETHM3  - INTERFACE TO ALLOW USER TO CHANGE THE HEPPNER-  .
C                      MAYNARD MODELS. (ONLY IMF SOUTHWARD MODELS AT   .
C                      THIS TIME.)                                     .
C                                                                      .
C  INPUT - CALLING SEQUENCE                                            .
C     MODLE  - INDEX FOR HEPPNER-MAYNARD MODELS                        .
C              (1 = A, 2 = BC, 3 = DE)                                 .
C                                                                      .
C  OUTPUT - CALLING SEQUENCE                                           .
C     ACHG,BCHG   - SEMI-MAJOR,MINOR AXIS OF NEW COORD. SYSTEM FOR H-M M
C     DXCHG,DYCHG - DISPLACEMENT OF CENTER OF ELIPSE FROM POLE FOR NEW S
C     AHM,BHM     - SEMI-MAJOR,MINOR AXIS OF OLD COORD. SYSTEM FOR H-M M
C     DYHM,DYHM   - DISPLACEMENT OF CENTER OF ELIPSE FROM POLE FOR OLD S
C     CHGV   -  PERCENT CHANGE IN POTENTIAL LEVELS                     .
C                                                                      .
C.......................................................................
      COMMON/CHGCOM/IABC,ICHGHM,CHGV,AHM(3),BHM(3),DXHM(3),DYHM(3),
     |   ACHG(3),BCHG(3),DXCHG(3),DYCHG(3)
      REAL KP
      real, parameter :: A(1:3) = (/17.88, 14.35, 15.96/)
      real, parameter :: B(1:3) = (/15.28, 14.15, 15.09/)
      real, parameter :: DX(1:3) = (/-3.27, -3.83, -2.50/)
      real, parameter :: DY(1:3) = (/ 0.57, -0.19,  0.76/)
C
C  VALID MODEL - SET STANDARD VALUES.
C  Default value for KP is 3.5 and maximum is 7.
      KP = AMIN1(RKP,7.)
      ICHGHM = 1
C  Change coefficients for all 3 Bz negative models
      CHGV = 0.223 + 0.222 * KP
      DO 1000 N=1,3
      AHM(N) = A(N)
      BHM(N) = B(N)
      DXHM(N) = DX(N)
      DYHM(N) = DY(N)
       ACHG(N) = AHM(N) * (0.721 + 0.087 * KP)
       BCHG(N) = BHM(N) * (0.735 + 0.082 * KP)
       DXCHG(N) = DXHM(N)
       DYCHG(N) = DYHM(N)
 1000 CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HMRMOD (DLAT,DLON,RATIO,POTEN)
C
C  DLAT,DLON = Magnetic latitude and longitude in radians, where
C              midnight is 0 degrees.
C  POTEN = electric potential in V
C
C  RATIO = A multiplicative factor for the latitude.  Is usually=1 
C          unless derivatives are sought to find the electric field as 
C          the difference in potential at 300 km.  It accounts for the 
C          fact that the magnetic field lines curve away from the Earth.
C          (Comment added by Barbara Emery)
C
C     IABC   - MODEL INDEX
C              IABC=1 IS FOR MODEL A   (IMF BZ < 0, BY < 0, NORTHERN HEM
C              IABC=2 IS FOR MODEL BC  (IMF BZ < 0, BY > 0, NORTHERN HEM
C              IABC=3 IS FOR MODEL DE  (IMF BZ < 0, BY < 0, NORTHERN HEM
C              IABC=4 IS FOR MODEL BCP (IMF BZ > 0, BY > 0, NORTHERN HEM
C              IABC=5 IS FOR MODEL BCPP(IMF BZ >>0, BY > 0, NORTHERN HEM
C              IABC=6 IS FOR MODEL DEP (IMF BZ > 0, BY < 0, NORTHERN HEM
C              IABC=7 IS FOR MODEL DEPP(IMF BZ >>0, BY < 0, NORTHERN HEM
C     NNMAX  - MAXIMUM ORDER OF POLYNOMIAL TO EVALUATE
C     ACHG,BCHG   - SEMI-MAJOR,MINOR AXIS OF NEW COORD. SYSTEM FOR H-M M
C     DXCHG,DYCHG - DISPLACEMENT OF CENTER OF ELIPSE FROM POLE FOR NEW S
C     AHM,BHM     - SEMI-MAJOR,MINOR AXIS OF OLD COORD. SYSTEM FOR H-M M
C     DYHM,DYHM   - DISPLACEMENT OF CENTER OF ELIPSE FROM POLE FOR OLD S
C     ICHGHM  - FLAG TO INDICATE WHETHER TO USE DEFAULT H-M PATTERN (=0)
C               TO USE MODIFIED PATTERN (=1) WITH ACHG,BCHG,DXCHG,DYCHG
C
C  INPUTS - COMMON/RUNCON/
C     RAD   - DEGREES PER RADIAN (=57.28...)
C
C  OUTPUTS - CALLING SEQUENCE
C    VALUE  - POTENTIAL AT (LAT, LON)  (KV)
C
C  SUBROUTINES CALLED
C    NONE (EXCEPT INTRINSIC FORTRAN FUNCTIONS)
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      COMMON /RUNCON/ PII,RAD,RE
C  COEFFICIENTS AND LIMITS FOR LEGENDRE POLYNOMIAL FIT TO H-M MODEL
C  (PLACED IN COMMON TO PREVENT LOADER FROM ERASING ARRAY.)
      COMMON/PCOM/COEF(18,18,7),P(18,18),NNMAX(7),NMAX,PKV51(7)
     |,XCO(18,18),DP(18,18),CONST(18,18),SP(18),CP(18),FN(18),FM(18)
      COMMON/SCOM/CMIN,CMAX
      COMMON/CHGCOM/IABC,ICHGHM,CHGV,AHM(3),BHM(3),DXHM(3),DYHM(3),
     |   ACHG(3),BCHG(3),DXCHG(3),DYCHG(3)
C
C
C   IABC=1 IS FOR MODEL A
C   IABC=2 IS FOR MODEL BC
C   IABC=3 IS FOR MODEL DE
C   IABC=4 IS FOR MODEL BCP
C   IABC=5 IS FOR MODEL BCPP
C   IABC=6 IS FOR MODEL DEP
C   IABC=7 IS FOR MODEL DEPP
C
      real, parameter :: ICHEM(1:7) = (/1,3,2,6,7,4,5/)
C
      NMAX = NNMAX(IABC)
C  Set IHEM=1,2 for S,N hemisphere
      IHEM = IFIX(DLAT*2./PII + 2.)
C
C  Change model with hemisphere except for model 1 (A) for now
      JABC = IABC
      IF (IHEM .EQ. 1) JABC = ICHEM(IABC)
C
      CMINRAD = CMIN/RAD
      CMN5RAD = (CMIN+5.)/RAD
      CMINP5 = CMIN+5.
      XLAT = ABS(DLAT) * SQRT(RATIO)
      XLON = DLON
C  SHOULD COORDINATE SYSTEM BE CHANGED?
      IF( ICHGHM.EQ.1 .AND. IABC .LE. 3) THEN
C                              YES.
       COLAT = PII/2. - DLAT
       TLONG = DLON - PII
       XX = DXHM(JABC)+AHM(JABC)*(COLAT*RAD*COS(TLONG)-DXCHG(JABC))/
     |  ACHG(JABC)
       YY = DYHM(JABC)+BHM(JABC)*(COLAT*RAD*SIN(TLONG)-DYCHG(JABC))/
     |  BCHG(JABC)
       XCOL = SQRT(XX**2 + YY**2)
       XLAT = (90.-XCOL) / RAD
       XLON = ATAN2(YY,XX) + PII
      ENDIF
C  CONVERT LATITUDE,LONGITUDE TO DIMENSIONLESS PARAMETER FOR EVALUATION
C  OF POLYNOMIAL
      ALPHA=2./(CMAX-CMIN)
      BETA=1.-ALPHA*CMAX
      CT=AMIN1( 1.0, AMAX1(XLAT*RAD,CMIN)*ALPHA+BETA)
      ST=SQRT(1.-CT*CT)
      SPH=SIN(XLON)
      CPH=COS(XLON)
C
      P(1,1) = 1.
      SP(2)=SPH
      CP(2)=CPH
      DO 4 M=3,NMAX
      SP(M)=SP(2)*CP(M-1)+CP(2)*SP(M-1)
    4 CP(M)=CP(2)*CP(M-1)-SP(2)*SP(M-1)
C
      P(2,1) = CT*P(1,1)
      P(2,2) = ST*P(1,1)
      DO 200 N=3,NMAX
      DO 100 M=1,N-1
      P(N,M)=CT*P(N-1,M)-CONST(N,M)*P(N-2,M)
C     DP(N,M)=CT*DP(N-1,M)-ST(I)*P(N-1,M)-
C    X   CONST(N,M)*DP(N-2,M)
 100  CONTINUE
      P(N,N)=ST*P(N-1,N-1)
C     DP(N,N)=ST*DP(N-1,N-1)+CT*P(N-1,N-1)
 200  CONTINUE
C
      P(1,1)=P(1,1)*XCO(1,1)
C
      POTEN=COEF(1,1,JABC)*P(1,1)
      DO 300 N=2,NMAX
      P(N,1)=P(N,1)*XCO(N,1)
      POTEN=POTEN+P(N,1)*COEF(N,1,JABC)
      DO 300 M=2,N
      POL=P(N,M)*XCO(N,M)
      P(M-1,N)=CP(M)*POL
      P(N,M)=SP(M)*POL
      POTEN=POTEN+P(M-1,N)*COEF(M-1,N,JABC) +
     |   P(N,M)*COEF(N,M,JABC)
  300 CONTINUE
      POTEN = (POTEN - PKV51(JABC)) * EXP(AMIN1(0.,XLAT*RAD-CMINP5)) *
     |         1.E+3 * CHGV
C
      RETURN
      END
