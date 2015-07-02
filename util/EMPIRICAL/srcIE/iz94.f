CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE IZEPOT (IMON,RMLAT,RMLT,BY,BZ,ET,EP,EPOT)
C          Obtain the electrostatic potential for the specified location
C          and conditions (hemisphere, season, IMF) and compute electric
C          field components from E = - grad (electrostatic potential).
C
C          INPUTS:
C            IMON  = Month of year (e.g., IMON=1 for Jan).  Season is
C                    determined by the combination of month and sign of
C                    RMLAT; it may be equinox (IMON = 3,4,9,10) summer or
C                    winter.
C            RMLAT = Magnetic latitude (deg) from 58 to 90 or -58 to -90
C            RMLT  = Magnetic Local Time (hr) from 0 to 24
C            BY,BZ = Interplanetary magnetic field components (nT)
C          RETURNS:
C            ET   = Etheta (magnetic equatorward*) E field component (V/m)
C            EP   = Ephi   (magnetic eastward)     E field component (V/m)
C            EPOT = Electostatic potential (kV)
C
C          * ET direction is along the magnetic meridian away from the
C            current hemisphere; i.e., when ET > 0, the direction is
C              southward when RMLAT > 0
C              northward when RMLAT < 0
C
C          Based on the PC routine, AM_SYNT.FOR by V. O. Papitashvili.
C
C          Barbara Emery (emery@ncar.ucar.edu) and Wm Golesorkhi, NCAR (6/95)
C          Mod Dec 96:  Correct model cofficient logic and electric field
C          component calculation and (R.Barnes, bozo@ncar.ucar.edu).

      COMMON /IZC1/ HIZA(24,1632), IZMODE, ALTE, ALAMN, STEPA

C          Local declarations (degrees to radians and hrs to radians):
      PARAMETER ( D2R = 0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147   ,
     +            H2R = 0.261799387799148538154743922256)

C          Establish hemsiphere and filter latitude to 58-90 deg
      IHEM  = 0
      IF (RMLAT .LT. 0) IHEM = 1
      ABLA = ABS(RMLAT)
      FMLA = AMAX1 (ALAMN, ABLA)

C          Obtain the potential at the requested point
      CALL IZNTRP (IMON, RMLT,FMLA,IHEM, BY,BZ, EPOT)

C          Calculate -(latitude gradient) by stepping 1 deg (model
C          resolution) along the meridian in each direction (flipping
C          coordinates when going over pole to keep lat <= 90).
      KPOL = 0
      XMLT  = RMLT
   10 XMLT1 = XMLT
      FMLA1 = FMLA + 1.
      IF (FMLA1 .GT. 90.) THEN
	FMLA1 = 180.  - FMLA1
	XMLT1 = XMLT1 + 12.
      ENDIF
      CALL IZNTRP (IMON, XMLT1, FMLA1  , IHEM, BY,BZ, P1)
      CALL IZNTRP (IMON, XMLT , FMLA-1., IHEM, BY,BZ, P2)
      IF (KPOL .EQ. 1) GO TO 20
      ET = (P1 - P2) / STEPA

C          Calculate -(lon gradient).  For most latitudes, step 15 degrees
C          in longitude (1 hr MLT = model resolution) along a great circle.
C          However, limit minimum latitude to the model minimum, distorting
C          the path onto a latitude line.  The step shrinks as latitude
C          increases and would become zero at the pole, but a different
C          scheme is used near the pole:  Assume lat=90 degrees and use
C          Art's trick where Ephi(90,lon) = Etheta(90,lon+90.)
      IF (FMLA .LT. 89.9) THEN
	SL = SIN (FMLA*D2R)
	CL = SQRT (1.-SL*SL)
	SP = SIN (15.*D2R)
	SA = SQRT (1.-SL*SL*SP*SP)
	FMLA1 = AMAX1 (ALAMN , ACOS(CL/SA)*R2D)
	CALL IZNTRP (IMON, XMLT+1.,FMLA1, IHEM, BY,BZ, P1)
	CALL IZNTRP (IMON, XMLT-1.,FMLA1, IHEM, BY,BZ, P2)
	STEP2 = 2.*ALTE*ASIN(CL*SP/SA)
      ELSE
	STEP2 = STEPA
	XMLT = XMLT + 6.
	FMLA = 90.
	KPOL = 1
	GO TO 10
      ENDIF
   20 EP = (P2 - P1) / STEP2
      IF (KPOL .EQ. 1) EP = -EP

C          Below model minimum lat, the potential is value at min lat
C  2/98:  RLAMN changed to ALAMN in next 3 lines
      IF (ABLA .LT. ALAMN) THEN
	ET = 0.
	EP = EP * COS(ALAMN*D2R)/SIN((90.-ABLA)*D2R)
      ENDIF

      RETURN
      END   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE IZINIT (IUN)
C          Initialize the IZMEM 94 model:  Load coefficients file and
C          model constants into common block IZC1.
C          INPUTS:
C             IUN = Fortran unit no. of (previously opened) coefficients file

      COMMON /IZC1/ HIZA(24,1632), IZMODE, ALTE, ALAMN, STEPA
C          HIZA   = coefficients array
C          IZMODE = computation mode.  Initialized here, but may be revised
C                   (before calling IZEPOT), to change contributing terms
C                   of the electrostatic potential formula:
C
C                          EPOT = AF*(BZLIM*HZ + BY*HY + H0)
C
C                 = 0,1  EPOT is based on corrective multiplier (AF), the
C                        Bz effect (BZLIM*HZ), the By effect (BY*HY), and
C                        the constant (H0).  This is the default, because
C                        the full potential with all corrective multipliers
C                        is the comon choice.  However, it can be changed
C                        to see how the potential is partitioned.
C                 = 2    the corrective multiplier (AF) is omitted (so the
C                        potential is generally too large).
C                 = 3    EPOT includes only the constant (H0).
C                 = 4    EPOT includes only the BZ effect (BZLIM*HZ).
C                 = 5    EPOT includes only the BY effect (BY*HY).
C          ALTE   = Mean Earth radius (6371 km) + 300 km, the altitude at
C                   which the electric field is computed.
C          ALAMN  = Absolute value (degrees) of model minimum magnetic
C                   latitude (with non-zero coefficients).
C          STEPA  = Twice the lat step used in Etheta calc (2 deg in radians)

      PARAMETER (D2R = 0.0174532925199432957692369076847 ,
     +           R2D = 57.2957795130823208767981548147)

C  Write out the valid limits
       write (6,"(1x,'IZMEM model has valid limits between 58-90 mlat',
     |  ' and between 0-24 MLT.'/1x,
     |  'SH values are derived by inverting the sign of By from NH.')")

      IZMODE = 1
      ALTE   = 6671.
      ALAMN  = 58.
      STEPA  = 2.*ALTE*D2R

C          Skip header text records when reading the combined coefficients file
      IROW = 0
      JROW = 0
   10 IF (IROW .EQ. 0) THEN
	READ (IUN,'(1X)',END=100)
      ELSE
	JROW = JROW + 1
	READ (IUN,'(14X,24F10.4)',END=100) (HIZA(I, JROW),I=1,24)
      ENDIF
      IROW = IROW+1
      IF (IROW .EQ. 35) IROW = 0
      GO TO 10

  100 RETURN
      END   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE IZNTRP (IMON, RMLT,AMLAT,IHEM, BY,BZ, EPOT)
C          Determine the electrostatic potential for the given location by
C          interpolating from the 4 adjacent (integer) model grid locations.
C
C          INPUTS:
C            IMON   = Month of year (e.g., IMON=1 for Jan).  Season is
C                     determined by the combination of IMON and IHEM.
C                     Season may be equinox (IMON = 3,4,9,10), summer
C                     or winter.
C            RMLT   = Magnetic local time (hr)
C            AMLAT  = Absolute value of magnetic latitude (degrees)
C            IHEM   = 0 (northern hemisphere) or 1 (southern hemisphere)
C            BY, BZ = IMF components (nT)
C          RETURNS:
C            EPOT   = electric potential (kV)
C
C          William Golesorkhi, HAO/NCAR (6/95)

      DIMENSION ISNS(24)
      SAVE      ISNS
      DATA      ISNS /2,2,0,0,1,1,1,1,0,0,2,2,1,1,0,0,2,2,2,2,0,0,1,1/

C          Define season, where 0=Equinox, 1=Summer, 2=Winter
      ISEA = ISNS(IMON+12*IHEM)

C          Filter MLT to range 1-24 hr
      N24 = INT (RMLT/24.)
      IF (RMLT .LT. 0.) N24 = N24 - 1
      FMLT = RMLT - REAL(N24)*24.
      IF (FMLT .LT. 1.) FMLT = FMLT + 24.

C          Establish coordinate indices required by IZMOD
      MLT  = INT (FMLT)
      LAT  = INT (AMLAT)
      MLT1 = MLT+1
      LAT1 = LAT+1
      IF (LAT1 .GT. 90) THEN
	LAT1 = 180 - LAT1
	MLT1 = MLT1 + 12
      ENDIF
      IF (MLT1 .GT. 24) MLT1 = MLT1 - 24

      FRACT = FMLT  - REAL (MLT)
      FRACL = AMLAT - REAL (LAT)

      IF (FRACT .EQ. 0. .AND. FRACL .EQ. 0.) THEN
	CALL IZMOD (MLT , LAT , IHEM, ISEA, BY,BZ, EPOT)
      ELSE
	CALL IZMOD (MLT , LAT1, IHEM, ISEA, BY,BZ, PNW)
	CALL IZMOD (MLT1, LAT1, IHEM, ISEA, BY,BZ, PNE)
	CALL IZMOD (MLT , LAT , IHEM, ISEA, BY,BZ, PSW)
	CALL IZMOD (MLT1, LAT , IHEM, ISEA, BY,BZ, PSE)
	PN = PNW + (PNE-PNW)*FRACT
	PS = PSW + (PSE-PSW)*FRACT
	EPOT = PS + (PN-PS) *FRACL
      ENDIF

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE IZMOD (MLT,LAT,IHEM, ISEA, BY,BZ, EPOT)
C          Determine the electrostatic potential for an integer MLT and
C          magnetic latitude.  This is known as the "IZMEM 94" model by
C          Papitashvili (JGR, 99, 13251 - 13262, 1994).
C
C          INPUTS:
C            MLT    = Magnetic local time (integer hr: 1,2,3 ... 24)
C            LAT    = Magnetic latitude (integer deg: 58,59,60 ... 90)
C            IHEM   = 0 (northern hemisphere) or 1 (southern hemisphere)
C            ISEA   = Season:  (0) equinox, (1) summer, or (2) winter
C            BY, BZ = IMF components (nT)
C          RETURNS:
C            EPOT  = electrostatic potential (kV)
C
C          Based on the PC routine, am_synt.for by V. O. Papitashvili.
C
C          William Golesorkhi, HAO/NCAR (6/95);
C          Mod Dec 96:  correct indexing into coefficient array (Two
C          conditions were reversed (By neg+Bz pos was swapped with
C          By pos+Bz neg) R Barnes (bozo@ncar.ucar.edu).

      COMMON /IZC1/ HIZA(24,1632), IZMODE, ALTE, ALAMN, STEPA

C          Local declarations
C          Originally, the coeffecients were placed in 48 seperate files
C          that the routine accessed by putting together different pieces of
C          the file name (eg. 'n'+'e'+'bh'+'1'+'.pot'=nebh1.pot) depending
C          on the hemisphere, season, component, and if BY and/or BZ was
C          positive or negative.  The coefficient files were combined into
C          one file which was read into HIZA.  An additional header line
C          was added to each block (which is omitted here) and a 90 degree
C          line containing the average of the values at lat=89.  Each
C          coefficient block contains 34 rows (90-57 deg lat) and 24 columns
C          (1 to 24 hrs MLT).  Latitude 57 contains 24 zeros.  The files were
C          split into two hemispheres (north and south), each hemisphere was
C          split into three seasons (equinox, summer, winter), each season
C          was split into three categories (Bh, By, Bz), and finally, each
C          category was split into two Bh blocks (H0 coef. for Bz-,Bz+),
C          4 By blocks (By coef for By-Bz-,By+Bz-,By-Bz+,By+Bz+), and 2 Bz
C          blocks (Bz coef for Bz-,Bz+).  The order of the original files
C          within HIZA is as follows:
C            nebh1.pot  nsbh1.pot  nwbh1.pot  sebh1.pot  ssbh1.pot  swbh1.pot
C            nebh2.pot  nsbh2.pot  nwbh2.pot  sebh2.pot  ssbh2.pot  swbh2.pot
C            neby3.pot  nsby3.pot  nwby3.pot  seby3.pot  ssby3.pot  swby3.pot
C            neby4.pot  nsby4.pot  nwby4.pot  seby4.pot  ssby4.pot  swby4.pot
C            neby5.pot  nsby5.pot  nwby5.pot  seby5.pot  ssby5.pot  swby5.pot
C            neby6.pot  nsby6.pot  nwby6.pot  seby6.pot  ssby6.pot  swby6.pot
C            nebz1.pot  nsbz1.pot  nwbz1.pot  sebz1.pot  ssbz1.pot  swbz1.pot
C            nebz2.pot  nsbz2.pot  nwbz2.pot  sebz2.pot  ssbz2.pot  swbz2.pot
C          The following constants are used to locate a certain block:
C            NROW = 34 = No. rows in a block (also no. magnetic latitudes)
C            NSEA = 3  = No. seasons (winter/summer/equinox)
C            NCAT = 8  = No. coefficients categories: 2(H0) + 4(HY) + 2(HZ)
C            IOFH =      No. of rows in 1 hemisphere
C            IOFS =      No. of rows in 1 season
C          Thus, to read the coefficients from file number 37, 'ssby5.pot',
C          which corresponds to BY<0, BZ>=0, and a southern summer:
C            IHEM=1, ISEA=1, ISBZ=2
C            LOC = NROW * NSEA * NCAT * IHEM  (LOC is in the 2nd half of files)
C            LOC = LOC  + NCAT * NROW * ISEA  (LOC is in the summer season
C            LOC = LOC  + NROW * 2            (LOC is in the By category
C            LOC = LOC  + NROW * ISBZ/Y       (LOC is at coeff block ssby5.pot
      PARAMETER (NROW = 34, NSEA = 3, NCAT = 8,
     +           IOFH = NROW*NCAT*NSEA, IOFS = NROW*NCAT)

C          Determine no blks to skip to reach desired category
      IF (BZ .LE. 0.) THEN
	ISBZ = 0
	ISBY = 0
	IF (BY .GT. 0.) ISBY = 1
      ELSE
	ISBZ = 1
	ISBY = 2
	IF (BY .GT. 0.) ISBY = 3
      ENDIF
C          previous code reversed order of two category offsets:
C     IF (BY .LE. 0 .AND. ISBZ .EQ. 0) ISBY = 0
C     IF (BY .LE. 0 .AND. ISBZ .EQ. 1) ISBY = 1    s.b. 2
C     IF (BY .GT. 0 .AND. ISBZ .EQ. 0) ISBY = 2    s.b. 1
C     IF (BY .GT. 0 .AND. ISBZ .EQ. 1) ISBY = 3

C          Determine the amplification factor, AF
      IF (IHEM .EQ. 0 .AND. ISEA .EQ. 1) AF = 3.2
      IF (IHEM .EQ. 1 .AND. ISEA .EQ. 1) AF = 3.0
      IF (IHEM .EQ. 0 .AND. ISEA .EQ. 0) AF = 2.3
      IF (IHEM .EQ. 1 .AND. ISEA .EQ. 0) AF = 2.7
      IF (IHEM .EQ. 0 .AND. ISEA .EQ. 2) AF = 1.5
      IF (IHEM .EQ. 1 .AND. ISEA .EQ. 2) AF = 2.5

C          Location within each block is ordered by colatitude
      KOLA = 90-LAT+1

C          LH0 = location for H0 coefficient
C          LHY = location for HY coefficient
C          LHZ = location for HZ coefficient
      I = IOFH*IHEM + IOFS*ISEA + KOLA
      LH0 = I + NROW* ISBZ
      LHY = I + NROW*(ISBY + 2)
      LHZ = I + NROW*(ISBZ + 6)
      
      H0 = HIZA (MLT,LH0)
      HY = HIZA (MLT,LHY)
      HZ = HIZA (MLT,LHZ)

      BZLIM = AMAX1 (-12., BZ)

      EPOT = 99999.
      IF (IZMODE .EQ. 0) EPOT = AF*(BZLIM*HZ + BY*HY + H0)
      IF (IZMODE .EQ. 1) EPOT = AF*(BZLIM*HZ + BY*HY + H0)
      IF (IZMODE .EQ. 2) EPOT = BZLIM*HZ + BY*HY + H0
      IF (IZMODE .EQ. 3) EPOT = H0
      IF (IZMODE .EQ. 4) EPOT = BZLIM*HZ
      IF (IZMODE .EQ. 5) EPOT = BY*HY
      IF (EPOT .EQ. 99999.) WRITE (6,'(''IZMOD:  IZMODE must be 0 to 5,
     + not='',I5)') IZMODE

      RETURN
      END   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
