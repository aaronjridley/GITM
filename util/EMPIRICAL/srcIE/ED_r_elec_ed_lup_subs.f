
      subroutine init_all

c      BLOCK DATA

      IMPLICIT NONE
      INCLUDE 'ED_R_elec_ed_lup_subs.inc'
      integer I,J
      integer LOADED
      real    PMAT(ENERGY_LEVELS,ALTITUDES)
      real    ENGMAT(ENERGY_LEVELS)
      real    EDEP_ALTS(ALTITUDES)
      real    EDEP_PRES(ALTITUDES)
      COMMON /R_EDEP/ LOADED,PMAT,ENGMAT,EDEP_ALTS,EDEP_PRES
      LOADED = 0
      DO I=1,ENERGY_LEVELS
        ENGMAT(I) = 0.0
        DO J=1,ALTITUDES
          PMAT(I,J) = 0.0
        ENDDO
      ENDDO
c      STOP
      END

C------------------------------------------------------------------------------

      FUNCTION R_LOAD_EDEP_LOOKUP ()
      use ModIoUnit, ONLY: io_unit_new,UnitTmp_
C******************************************************************************
C
C  THIS FUNCTION READS THE ASCII ENERGY DEPOSITION VALUES AND CREATES AN ENERGY
C  DEPOSITION LOOKUP TABLE WHICH WILL BE USED THROUGHOUT THE CODE.  
C  COMMUNICATION OF THE LOOKUP TABLE IS THROUGH COMMON MEMORY.
C
C  INPUTS - DIRECTLY PASSED
C           NONE
C
C  INPUTS - PASSED BY COMMON
C           NONE
C
C  OUTPUTS - DIRECTLY PASSED
C           R_LOAD_EDEP_LOOKUP              integer   STATUS OF LOAD
C                                                       -1 LOAD FAILED
C                                                       +1 LOAD SUCCEED
C
C  OUTPUTS - PASSED BY NAMED COMMON /R_EDEP/
C           LOADED                          LOGICAL*4   INDICATING STATe
C                                                       OF DATA FILE
C                                                       1  FILE LOADED
C                                                       0 FILE NOT LOADED
C           PMAT(ENERGY_LEVELS,ALTITUDES)   real     ENERGY DEPOSITION VALUEs
C                                                       FOR INCIDENT ISOTROPIC 
C                                                       FLUX [eV/cm]
C           ENGMAT(ENERGY_LEVELS)           real      ENERGY VALUES [MeV]
C           EDEP_ALTS(ALTITUDES)            real      ALTITUDE VALUES USED BY
C                                                       THE ENERGY DEPOSITION 
C                                                       PROFILES [km]
C           EDEP_PERS(ALTITUDES)            real      PRESSURE VALUES USED BY
C                                                       THE ENERGY DEPOSITION 
C                                                       PROFILES [mbar]
C
C******************************************************************************

      IMPLICIT NONE

C  Include Electron Energy Deposition Parameters

      INCLUDE 'ED_R_elec_ed_lup_subs.inc'

C  DEFINE LOCAL VARIABLES

      INTEGER   STORE
      integer I,J
      integer ITERME
      integer IOERR
      integer LEVELS
      integer LIST
      integer PTR
      integer R_LOAD_EDEP_LOOKUP 
      real    ENGMEV
      real    ALPHA
      real    FLUX

C Define common blocks

      INTEGER   LOADED, iUnit
      real    PMAT(ENERGY_LEVELS,ALTITUDES)
      real    ENGMAT(ENERGY_LEVELS)
      real    EDEP_ALTS(ALTITUDES)
      real    EDEP_PRES(ALTITUDES)
      COMMON /R_EDEP/ LOADED,PMAT,ENGMAT,EDEP_ALTS,EDEP_PRES

C******************************************************************************
C
C  ENERGY_LEVELS     PARAMETER  NUMBER OF ENERGY LEVELS
C  ALTITUDES         PARAMETER  NUMBER OF ALTITUDE LEVELS
C  STORE             LOGICAL*4  FLAG TO STORE THE ENERGY DEPOSITION PROFILE
C                               1  EVERYTHING CHECKS OUT, SAVE PROFILE
C                               0 ERRORS ENCOUNTERED, DO NOT SAVE PROFILE
C  I                 integer  INDEX COUNTER
C  J                 integer  INDEX COUNTER
C  ITERME            integer  NUMBER OF ENERGY LEVELS RECORDED
C  IOERR             integer  STATUS OF I/O OPERATION
C  R_LOAD_EDEP_LOOKUP
C                    integer  STATUS OF LOAD
C                               -1 LOAD FAILED
C                               +1 LOAD SUCCEED
C  LEVELS            integer  PREVIOUS ENTRY LOCATION IN ENERGY MATRIX -- 
C                               SHOULD ALWAYS BE ZERO TO INDICATE EACH 
C                               DEPOSITION TABLE IS UNIQUE
C                               DEPOSITION TABLE
C  LIST              integer  STATUS OF I/O OPERATION
C  PTR               integer  ALTITUDE LEVEL INDEX OF PROFILE
C  ENGMEV            real     TABLE ENERGY [MeV]
C  ALPHA             real     PARAMETER DESCRIBING DISTRIBUTION PITCH ANGLE 
C                               SHAPE -- NOTE THAT THIS SHOULD ALWAYS BE ZERO
C                               FOR THE ISOTROPIC CASE
C  FLUX              real     ENERGY DEPOSITION TABLE FLUX [erg-cm-eV]
C  LOADED            LOGICAL*4  /R_EDEP/  LOGICAL INDICATING STATE OF DATA FILE
C                                         1  FILE HAS BEEN LOADED
C                                         0 FILE HAS NOT BEEN LOADED
C  PMAT(ENERGY_LEVELS,ALTITUDES)   
C                    real     /R_EDEP/  ENERGY DEPOSITION VALUES FOR INCIDENT
C                                         ISOTROPIC FLUX [eV/cm]
C  ENGMAT(ENERGY_LEVELS)           
C                    real     /R_EDEP/  ENERGY VALUES [MeV]
C  EDEP_ALTS(ALTITUDES)            
C                    real     /R_EDEP/  ALTITUDE VALUES USED BY THE ENERGY 
C                                         DEPOSITION PROFILES [km]
C  EDEP_PRES(ALTITUDES)            
C                    real     /R_EDEP/  PRESSURE VALUES USED BY THE ENERGY 
C                                         DEPOSITION PROFILES [mbar]
C
C******************************************************************************

C******************************************************************************
C  FORMAT STATEMENTS
C******************************************************************************

  800 FORMAT(1X,'WARNING: Energy ',1PE12.5,' MeV is not for an ',
     1       'isotropic plasma distribution (alpha =',E10.3,')',/,12X,
     2       'Isotropic distribution assumed!')

  810 FORMAT(1X,'WARNING: Energy ',1PE12.5,' MeV was already ',
     1       'processed, profile skipped')

  820 FORMAT(1X,'WARNING: Expected',I3,'profiles for these energy ',
     1       'levels, but',I3,' were found')

C******************************************************************************

C  Indicate start

      R_LOAD_EDEP_LOOKUP = -1

C  Initialize energy level counter

      ITERME = 0

C Open the file of energy deposition values

C     Must set this before calling CON_io_unit_new, so if it is running
C     stand-alone, this routine doesn't have to do anything.

      iUnit = io_unit_new()

      OPEN(UNIT=iUnit,FILE=
     |     'EIE/isoelec_edep.dat',FORM='FORMATTED',
     +       STATUS='OLD', IOSTAT=IOERR)

C     Loop through the number of energy values

      DO I=1,ENERGY_LEVELS

C  Read the Electron lookup table and associated values

        READ(iUnit,*,IOSTAT=IOERR) ENGMEV,ALPHA,LEVELS

C  Handle any errors

        IF (IOERR .NE. 0) THEN
           WRITE(*,'(1X,A25,/,8X,A50,I5)',IOSTAT = IOERR) 
     1           'FATAL: Program aborted --',
     2           'Unable to read energy deposition look-up table.',
     3           IOERR
           GOTO 999
        ENDIF

C  Print warning message for non-isotropic distribltion profiles

        IF (ALPHA .NE. 0.0) WRITE(*,800) ENGMEV,ALPHA,LEVELS

C  Search for energy to determine if it is already in the list

        IF (ITERME .EQ. 0) THEN

C  No energies recorded yet, so add the first one to the list

          ITERME = 1
          ENGMAT(ITERME) = ENGMEV
          STORE = 1

        ELSE

C  At least one energy is in the list already
C  Set flag to indicate where it appears in the list and search the list

          LIST = 0

          DO J=1,ITERME
            IF (ENGMAT(J) .EQ. ENGMEV) LIST = J
          ENDDO

C  If the energy does not show up in the list, append it and indicate that
C  we want to acquire the profile

          IF (LIST .EQ. 0) THEN
            ITERME = ITERME + 1
            ENGMAT(ITERME) = ENGMEV
            STORE = 1

C  If the energy does show up in the list, write out a warning message and
C  do not acquire the profile

          ELSE
            WRITE(*,810) ENGMEV
            STORE = 0
          ENDIF

        ENDIF

C  Read in the energy deposition values and handle any errors

        DO J = 1, LEVELS

          READ(iUnit,*,IOSTAT=IOERR) PTR,FLUX
          IF (IOERR .NE. 0) THEN
             WRITE(*,'(5X,A55,I5)',IOSTAT = IOERR)
     +  'FATAL: Program aborted -- Unable to read EAR.DAT file.',IOERR
             GOTO 999
          ENDIF

C  Store the data in the energy deposition matrix only if requested

          IF (STORE.eq.1) PMAT(ITERME,PTR) = FLUX

        ENDDO
      ENDDO

C  Do we have all of the expected energies?

      IF (ENERGY_LEVELS .NE. ITERME) THEN
        WRITE(*,820) ENERGY_LEVELS,ITERME
      ENDIF

C  Indicate that the data file has been loaded

      LOADED = 1

 999  CLOSE(iUnit)

C  Set the altitudes in km

      DO I=1,12
        EDEP_ALTS(I) = 5.0 * I
      ENDDO

      DO I=13,32
        EDEP_ALTS(I) = 60.0 + (I - 12) * 3.0
      ENDDO

      DO I=33,88
        EDEP_ALTS(I) = 120.0 + (I - 32) * 5.0
      ENDDO

C  Set the Pressure -- the 1976 US Standard Atmosphere is used
C  At ground, EDEP_PRES(0) = 1013.25 mbar

      EDEP_PRES(1) = 540.479
      EDEP_PRES(2) = 264.989
      EDEP_PRES(3) = 121.110
      EDEP_PRES(4) = 55.2929
      EDEP_PRES(5) = 25.4920
      EDEP_PRES(6) = 11.9700
      EDEP_PRES(7) = 5.74590
      EDEP_PRES(8) = 2.87140
      EDEP_PRES(9) = 1.49100
      EDEP_PRES(10) = 797789.0E-6
      EDEP_PRES(11) = 425249.0E-6
      EDEP_PRES(12) = 219579.0E-6
      EDEP_PRES(13) = 145150.0E-6
      EDEP_PRES(14) = 946089.0E-7
      EDEP_PRES(15) = 607360.0E-7
      EDEP_PRES(16) = 338619.0E-7
      EDEP_PRES(17) = 238809.0E-7
      EDEP_PRES(18) = 146730.0E-7
      EDEP_PRES(19) = 889229.0E-8
      EDEP_PRES(20) = 531050.0E-8
      EDEP_PRES(21) = 312589.0E-8
      EDEP_PRES(22) = 183590.0E-8
      EDEP_PRES(23) = 108009.0E-8
      EDEP_PRES(24) = 637650.0E-9
      EDEP_PRES(25) = 379480.0E-9
      EDEP_PRES(26) = 231440.0E-9
      EDEP_PRES(27) = 144770.0E-9
      EDEP_PRES(28) = 931879.0E-10
      EDEP_PRES(29) = 626140.0E-10
      EDEP_PRES(30) = 444729.0E-10
      EDEP_PRES(31) = 330219.0E-10
      EDEP_PRES(32) = 253819.0E-10
      EDEP_PRES(33) = 173539.0E-10
      EDEP_PRES(34) = 125050.0E-10
      EDEP_PRES(35) = 935679.0E-11
      EDEP_PRES(36) = 720280.0E-11
      EDEP_PRES(37) = 566910.0E-11
      EDEP_PRES(38) = 454219.0E-11
      EDEP_PRES(39) = 369300.0E-11
      EDEP_PRES(40) = 303949.0E-11
      EDEP_PRES(41) = 252780.0E-11
      EDEP_PRES(42) = 212100.0E-11
      EDEP_PRES(43) = 179359.0E-11
      EDEP_PRES(44) = 152710.0E-11
      EDEP_PRES(45) = 130809.0E-11
      EDEP_PRES(46) = 112659.0E-11
      EDEP_PRES(47) = 974909.0E-12
      EDEP_PRES(48) = 847360.0E-12
      EDEP_PRES(49) = 739419.0E-12
      EDEP_PRES(50) = 647560.0E-12
      EDEP_PRES(51) = 569000.0E-12
      EDEP_PRES(52) = 501489.0E-12
      EDEP_PRES(53) = 443239.0E-12
      EDEP_PRES(54) = 392760.0E-12
      EDEP_PRES(55) = 348880.0E-12
      EDEP_PRES(56) = 310589.0E-12
      EDEP_PRES(57) = 277080.0E-12
      EDEP_PRES(58) = 247669.0E-12
      EDEP_PRES(59) = 221779.0E-12
      EDEP_PRES(60) = 198940.0E-12
      EDEP_PRES(61) = 178739.0E-12
      EDEP_PRES(62) = 160829.0E-12
      EDEP_PRES(63) = 144919.0E-12
      EDEP_PRES(64) = 130760.0E-12
      EDEP_PRES(65) = 118129.0E-12
      EDEP_PRES(66) = 106850.0E-12
      EDEP_PRES(67) = 967510.0E-13
      EDEP_PRES(68) = 877040.0E-13
      EDEP_PRES(69) = 795851.0E-13
      EDEP_PRES(70) = 722849.0E-13
      EDEP_PRES(71) = 657166.0E-13
      EDEP_PRES(72) = 597960.0E-13
      EDEP_PRES(73) = 544555.0E-13
      EDEP_PRES(74) = 496299.0E-13
      EDEP_PRES(75) = 452688.0E-13
      EDEP_PRES(76) = 413199.0E-13
      EDEP_PRES(77) = 377428.0E-13
      EDEP_PRES(78) = 344980.0E-13
      EDEP_PRES(79) = 315540.0E-13
      EDEP_PRES(80) = 288780.0E-13
      EDEP_PRES(81) = 264469.0E-13
      EDEP_PRES(82) = 242339.0E-13
      EDEP_PRES(83) = 222196.0E-13
      EDEP_PRES(84) = 203840.0E-13
      EDEP_PRES(85) = 187107.0E-13
      EDEP_PRES(86) = 171840.0E-13
      EDEP_PRES(87) = 157907.0E-13
      EDEP_PRES(88) = 145180.0E-13

C  Indicate success

      R_LOAD_EDEP_LOOKUP = 0

      RETURN
      END

C-------------------------------------------------------------------------------

      FUNCTION R_EDEP_ALT_INDEX (GALT,UNIT)

C******************************************************************************
C
C  THIS FUNCTION DETERMINES THE CLOSEST BIN NUMBER TO THE REQUESTED ALTITUDE
C
C  INPUTS - DIRECTLY PASSED
C           GALT                  real      THE REQUESTED ALTITUDE
C           UNIT                  integer   UNIT OF ALTITUDE
C                                             1 = cm
C                                             2 = m
C                                             3 = km
C                                             4 = mbar
C
C  INPUTS - PASSED BY NAMED COMMON /R_EDEP/
C           EDEP_ALTS(ALTITUDES)  real      ALTITUDE VALUES USED BY THE 
C                                             ENERGY DEPOSITION PROFILES [km]
C
C  OUTPUTS - DIRECTLY PASSED
C           R_EDEP_ALT_INDEX      integer   THE CLOSEST ALTITUDE INDEX
C
C  OUTPUTS - PASSED BY COMMON
C           NONE
C
C******************************************************************************

      IMPLICIT NONE

C  Include Electron Energy Deposition Parameters

      INCLUDE 'ED_R_elec_ed_lup_subs.inc'

C  DEFINE LOCAL VARIABLES

      integer R_EDEP_ALT_INDEX
      integer I
      integer UNIT
      integer IOERR
      real    GALT
      real    ALT
      real    UFACT

C Define common blocks

      INTEGER   LOADED
      real    PMAT(ENERGY_LEVELS,ALTITUDES)
      real    ENGMAT(ENERGY_LEVELS)
      real    EDEP_ALTS(ALTITUDES)
      real    EDEP_PRES(ALTITUDES)
      COMMON /R_EDEP/ LOADED,PMAT,ENGMAT,EDEP_ALTS,EDEP_PRES

C******************************************************************************
C
C  ENERGY_LEVELS     PARAMETER  NUMBER OF ENERGY LEVELS
C  ALTITUDES         PARAMETER  NUMBER OF ALTITUDE LEVELS
C  R_EDEP_ALT_INDEX  integer  CLOSEST ALTITUDE INDEX
C  I                 integer  INDEX COUNTER
C  IOERR             integer  STATUS OF I/O OPERATION
C  UNIT              integer  UNIT FOR ALTITUDE
C                               1 = cm
C                               2 = m
C                               3 = km
C  GALT              real     ALTITUDE FOR WHICH INDEX HAS BEEN REQUESTED
C  ALT               real     CONVERTED ALTITUDE FOR WHICH INDEX HAS BEEN 
C                               REQUESTED [km]
C  UFACT             real     UNIT FACTOR FOR CONVERSION INTO [km]
C                               UFACT = 1.0E-5 FOR cm
C                               UFACT = 1.0E-3 FOR m
C                               UFACT = 1.0 FOR km
C                               UFACT = 1.0 FOR mbar
C  LOADED            LOGICAL*4  /R_EDEP/  LOGICAL INDICATING STATE OF DATA FILE
C                                         1  FILE HAS BEEN LOADED
C                                         0 FILE HAS NOT BEEN LOADED
C  PMAT(ENERGY_LEVELS,ALTITUDES)   
C                    real     /R_EDEP/  ENERGY DEPOSITION VALUES FOR INCIDENT
C                                         ISOTROPIC FLUX [eV/cm]
C  ENGMAT(ENERGY_LEVELS)           
C                    real     /R_EDEP/  ENERGY VALUES [MeV]
C  EDEP_ALTS(ALTITUDES)            
C                    real     /R_EDEP/  ALTITUDE VALUES USED BY THE ENERGY 
C                                         DEPOSITION PROFILES [km]
C  EDEP_PRES(ALTITUDES)            
C                    real     /R_EDEP/  PRESSURE VALUES USED BY THE ENERGY 
C                                         DEPOSITION PROFILES [mbar]
C
C******************************************************************************

C  Set indicator indicating off the altitude grid

      R_EDEP_ALT_INDEX = -1

C  Have the energy deposition tables been loaded?

      IF (LOADED.eq.0) THEN
        WRITE (*,'(1X,A55,I5)',IOSTAT=IOERR)
     +         'SKIP: Energy Deposition Files Not Loaded.',IOERR
        RETURN
      ENDIF

C  Set the unit factor

      UFACT = 1.0
      IF (MOD(UNIT,4) .EQ. 1) UFACT = 1.0E-5
      IF (MOD(UNIT,4) .EQ. 2) UFACT = 1.0E-3
      ALT = GALT*UFACT

C  Handle the mbar unit separately

      IF (UNIT .EQ. 4) THEN

C  Handle pressure distance -- log
C  Signify bad value if altitues is above or below altitude grid

        IF (ALT .LT. EDEP_PRES(ALTITUDES)) RETURN
        IF (ALT .GT. EDEP_PRES(1)) RETURN
      
C  MINIMUM INDEX

        IF (ALT .GT. 10.0**(0.5*(LOG10(EDEP_PRES(1))+
     +                                     LOG10(EDEP_PRES(2))))) THEN
          R_EDEP_ALT_INDEX = 1

C  MAXIMUM INDEX

        ELSE IF (ALT .LE. 10.0**(0.5*(LOG10(EDEP_PRES(ALTITUDES-1))+
     +                             LOG10(EDEP_PRES(ALTITUDES))))) THEN
          R_EDEP_ALT_INDEX = ALTITUDES

C  ALL OTHER INDICIES

        ELSE
          DO I=2,ALTITUDES-1
            IF ((ALT .LE. 10.0**(0.5*(LOG10(EDEP_PRES(I-1))+
     1                           LOG10(EDEP_PRES(I))))) .AND. (ALT .GT. 
     2       10.0**(0.5*(LOG10(EDEP_PRES(I+1))+LOG10(EDEP_PRES(I))))))
     3          R_EDEP_ALT_INDEX = I
          ENDDO
        ENDIF

      ELSE

C  Handle altitude distance -- linear
C  Signify bad value if altitues is above or below altitude grid

        IF (ALT .GT. EDEP_ALTS(ALTITUDES)) RETURN
        IF (ALT .LT. EDEP_ALTS(1)) RETURN
      
C  MINIMUM INDEX

        IF (ALT .LT. 0.5*(EDEP_ALTS(1)+EDEP_ALTS(2))) THEN
          R_EDEP_ALT_INDEX = 1

C  MAXIMUM INDEX

        ELSE IF (ALT .GE. 0.5*(EDEP_ALTS(ALTITUDES-1)+
     +                                      EDEP_ALTS(ALTITUDES))) THEN
          R_EDEP_ALT_INDEX = ALTITUDES

C  ALL OTHER INDICIES

        ELSE
          DO I=2,ALTITUDES-1
            IF ((ALT .GE. 0.5*(EDEP_ALTS(I-1)+EDEP_ALTS(I))) .AND.
     1          (ALT .LT. 0.5*(EDEP_ALTS(I+1)+EDEP_ALTS(I))))
     2          R_EDEP_ALT_INDEX = I
          ENDDO
        ENDIF

      ENDIF

      RETURN
      END

C------------------------------------------------------------------------------

      FUNCTION R_EDEP_ENG_INDEX (GENG,UNIT)

C******************************************************************************
C
C  THIS FUNCTION DETERMINES THE CLOSEST BIN NUMBER TO THE REQUESTED ENERGY
C
C  INPUTS - DIRECTLY PASSED
C           GENG                  real      THE REQUESTED ENERGY
C           UNIT                  integer   UNIT OF ENERGY
C                                             1 = erg
C                                             2 = Joule
C                                             3 = eV
C                                             4 = keV
C                                             5 = MeV
C
C  INPUTS - PASSED BY NAMED COMMON /R_EDEP/
C           ENGMAT(ENERGY_LEVELS) real      ENERGY VALUES USED BY THE ENERGY
C                                             DEPOSITION PROFILES [MeV]
C
C  OUTPUTS - DIRECTLY PASSED
C           R_EDEP_ENG_INDEX      integer   THE CLOSEST ALTITUDE INDEX
C
C  OUTPUTS - PASSED BY COMMON
C           NONE
C
C******************************************************************************

      IMPLICIT NONE

C  Include Electron Energy Deposition Parameters

      INCLUDE 'ED_R_elec_ed_lup_subs.inc'

C  DEFINE LOCAL VARIABLES

      integer R_EDEP_ENG_INDEX
      integer I
      integer UNIT
      integer IOERR
      real    GENG
      real    ENG
      real    UFACT

C Define common blocks

      INTEGER LOADED
      real    PMAT(ENERGY_LEVELS,ALTITUDES)
      real    ENGMAT(ENERGY_LEVELS)
      real    EDEP_ALTS(ALTITUDES)
      real    EDEP_PRES(ALTITUDES)
      COMMON /R_EDEP/ LOADED,PMAT,ENGMAT,EDEP_ALTS,EDEP_PRES

C******************************************************************************
C
C  ENERGY_LEVELS     PARAMETER  NUMBER OF ENERGY LEVELS
C  ALTITUDES         PARAMETER  NUMBER OF ALTITUDE LEVELS
C  R_EDEP_ENG_INDEX  integer  CLOSEST ENERGY INDEX
C  I                 integer  INDEX COUNTER
C  IOERR             integer  STATUS OF I/O OPERATION
C  UNIT              integer  UNIT FOR ALTITUDE
C                               1 = erg
C                               2 = Joule
C                               3 = eV
C                               4 = keV
C                               5 = MeV
C  GENG              real     ENERGY FOR WHICH INDEX HAS BEEN REQUESTED
C  ENG               real     ENERGY FOR WHICH INDEX HAS BEEN REQUESTED [km]
C  UFACT             real     UNIT FACTOR FOR CONVERSION INTO [MeV]
C                               UFACT = 6.242E+5 FOR erg
C                               UFACT = 6.242E+12 FOR Joule
C                               UFACT = 1.0E-6 FOR eV
C                               UFACT = 1.0E-3 FOR keV
C                               UFACT = 1.0 FOR MeV
C  LOADED            LOGICAL*4  /R_EDEP/  LOGICAL INDICATING STATE OF DATA FILE
C                                         1  FILE HAS BEEN LOADED
C                                         0 FILE HAS NOT BEEN LOADED
C  PMAT(ENERGY_LEVELS,ALTITUDES)   
C                    real     /R_EDEP/  ENERGY DEPOSITION VALUES FOR INCIDENT
C                                         ISOTROPIC FLUX [eV/cm]
C  ENGMAT(ENERGY_LEVELS)           
C                    real     /R_EDEP/  ENERGY VALUES [MeV]
C  EDEP_ALTS(ALTITUDES)            
C                    real     /R_EDEP/  ALTITUDE VALUES USED BY THE ENERGY 
C                                         DEPOSITION PROFILES [km]
C  EDEP_PRES(ALTITUDES)            
C                    real     /R_EDEP/  PRESSURE VALUES USED BY THE ENERGY 
C                                         DEPOSITION PROFILES [mbar]
C
C******************************************************************************

C  Set indicator indicating off the energy grid

      R_EDEP_ENG_INDEX = -1

C  Have the energy deposition tables been loaded?

      IF (LOADED.eq.0) THEN
        WRITE (*,'(1X,A55,I5)',IOSTAT=IOERR)
     +         'SKIP: Energy Deposition Files Not Loaded.',IOERR
        RETURN
      ENDIF

C  Set the unit factor

      UFACT = 1.0
      IF (MOD(UNIT,5) .EQ. 1) UFACT = 6.242E+5
      IF (MOD(UNIT,5) .EQ. 2) UFACT = 6.242E+12
      IF (MOD(UNIT,5) .EQ. 3) UFACT = 1.0E-6
      IF (MOD(UNIT,5) .EQ. 4) UFACT = 1.0E-3
      ENG = GENG*UFACT

C  Signify bad value if energy is above or below energy grid

      IF (ENG .LT. ENGMAT(ENERGY_LEVELS)) RETURN
      IF (ENG .GT. ENGMAT(1)) RETURN
      
C  MINIMUM INDEX

      IF (ENG .GT. 10.0**(0.5*(LOG10(ENGMAT(1))+
     +                                       LOG10(ENGMAT(2))))) THEN
        R_EDEP_ENG_INDEX = 1

C  MAXIMUM INDEX

      ELSE IF (ENG .LE. 10.0**(0.5*(LOG10(ENGMAT(ENERGY_LEVELS-1))+
     +                           LOG10(ENGMAT(ENERGY_LEVELS))))) THEN
        R_EDEP_ENG_INDEX = ENERGY_LEVELS

C  ALL OTHER INDICIES

      ELSE
        DO I=2,ENERGY_LEVELS-1
          IF ((ENG .LE. 10.0**(0.5*(LOG10(ENGMAT(I-1))+
     1         LOG10(ENGMAT(I))))) .AND.  (ENG .GT. 10.0**(0.5*(
     2         LOG10(ENGMAT(I+1))+LOG10(ENGMAT(I))))))
     3                                           R_EDEP_ENG_INDEX = I
        ENDDO
      ENDIF

      RETURN
      END

C------------------------------------------------------------------------------

      FUNCTION R_EDEP_ALT_VALUE (ALT_INDEX,UNIT)

C******************************************************************************
C
C  THIS FUNCTION DETERMINES THE ALTITUDE ASSOCIATED WITH THE BIN INDEX IN THE
C  UNITS REQUESTED
C
C  INPUTS - DIRECTLY PASSED
C           ALT_INDEX             integer   THE REQUESTED BIN INDEX
C           UNIT                  integer   UNIT OF ALTITUDE REQUESTED
C                                             1 = cm
C                                             2 = m
C                                             3 = km
C                                             4 = mbar
C
C  INPUTS - PASSED BY NAMED COMMON /R_EDEP/
C           EDEP_ALTS(ALTITUDES)  real      ALTITUDE VALUES USED BY THE ENERGY
C                                             DEPOSITION PROFILES [km]
C
C  OUTPUTS - DIRECTLY PASSED
C           R_EDEP_ALT_VALUE      real      THE ALTITUDE OF THE REQUESTED 
C                                             INDEX EXPRESSED IN THE REQUESTED
C                                             UNITS
C
C  OUTPUTS - PASSED BY COMMON
C           NONE
C
C******************************************************************************

      IMPLICIT NONE

C  Include Electron Energy Deposition Parameters

      INCLUDE 'ED_R_elec_ed_lup_subs.inc'

C  DEFINE LOCAL VARIABLES

      integer UNIT
      integer ALT_INDEX
      integer IOERR
      real    R_EDEP_ALT_VALUE
      real    UFACT

C Define common blocks

      INTEGER   LOADED
      real    PMAT(ENERGY_LEVELS,ALTITUDES)
      real    ENGMAT(ENERGY_LEVELS)
      real    EDEP_ALTS(ALTITUDES)
      real    EDEP_PRES(ALTITUDES)
      COMMON /R_EDEP/ LOADED,PMAT,ENGMAT,EDEP_ALTS,EDEP_PRES

C******************************************************************************
C
C  ENERGY_LEVELS     PARAMETER  NUMBER OF ENERGY LEVELS
C  ALTITUDES         PARAMETER  NUMBER OF ALTITUDE LEVELS
C  IOERR             integer  STATUS OF I/O OPERATION
C  UNIT              integer  UNIT FOR ALTITUDE
C                               1 = cm
C                               2 = m
C                               3 = km
C  ALT_INDEX         integer  REQUESTED ALTITUDE BIN
C  R_EDEP_ALT_VALUE  real     ALTITUDE OF BIN INDEX
C  UFACT             real     UNIT FACTOR FOR CONVERSION INTO [km]
C                               UFACT = 1.0E+5 FOR cm
C                               UFACT = 1.0E+3 FOR m
C                               UFACT = 1.0 FOR km
C  LOADED            LOGICAL*4  /R_EDEP/  LOGICAL INDICATING STATE OF DATA FILE
C                                         1  FILE HAS BEEN LOADED
C                                         0 FILE HAS NOT BEEN LOADED
C  PMAT(ENERGY_LEVELS,ALTITUDES)   
C                    real     /R_EDEP/  ENERGY DEPOSITION VALUES FOR INCIDENT
C                                         ISOTROPIC FLUX [eV/cm]
C  ENGMAT(ENERGY_LEVELS)           
C                    real     /R_EDEP/  ENERGY VALUES [MeV]
C  EDEP_ALTS(ALTITUDES)            
C                    real     /R_EDEP/  ALTITUDE VALUES USED BY THE ENERGY 
C                                         DEPOSITION PROFILES [km]
C  EDEP_PRES(ALTITUDES)            
C                    real     /R_EDEP/  PRESSURE VALUES USED BY THE ENERGY 
C                                         DEPOSITION PROFILES [mbar]
C
C******************************************************************************

C  Set indicator indicating off the altitude grid

      R_EDEP_ALT_VALUE = -1.0

C  Have the energy deposition tables been loaded?

      IF (LOADED.eq.0) THEN
        WRITE (*,'(1X,A55,I5)',IOSTAT=IOERR)
     +         'SKIP: Energy Deposition Files Not Loaded.',IOERR
        RETURN
      ENDIF

C  Set the unit factor

      UFACT = 1.0
      IF (MOD(UNIT,4) .EQ. 1) UFACT = 1.0E+5
      IF (MOD(UNIT,4) .EQ. 2) UFACT = 1.0E+3

C  Signify bad value if altitues is above or below altitude grid

      IF (ALT_INDEX .LT. 1) RETURN
      IF (ALT_INDEX .GT. ALTITUDES) RETURN
      
C  Determine value of altitude and shift to requested units

      IF (UNIT .EQ. 4) THEN
        R_EDEP_ALT_VALUE = UFACT*EDEP_PRES(ALT_INDEX)
      ELSE
        R_EDEP_ALT_VALUE = UFACT*EDEP_ALTS(ALT_INDEX)
      ENDIF

      RETURN
      END

C------------------------------------------------------------------------------

      FUNCTION R_EDEP_ENG_VALUE (ENG_INDEX,UNIT)

C******************************************************************************
C
C  THIS FUNCTION DETERMINES THE ENERGY ASSOCIATED WITH THE BIN INDEX IN THE
C  UNITS REQUESTED
C
C  INPUTS - DIRECTLY PASSED
C           ENG_INDEX             integer   THE REQUESTED BIN INDEX
C           UNIT                  integer   UNIT OF ENERGY REQUESTED
C                                             1 = erg
C                                             2 = Joule
C                                             3 = eV
C                                             4 = keV
C                                             5 = MeV
C
C  INPUTS - PASSED BY NAMED COMMON /R_EDEP/
C           ENGMAT(ENERGY_LEVELS) real      ENERGY VALUES USED BY THE ENERGY
C                                             DEPOSITION PROFILES [MeV]
C
C  OUTPUTS - DIRECTLY PASSED
C           R_EDEP_ENG_VALUE      real      THE CLOSEST ALTITUDE INDEX
C
C  OUTPUTS - PASSED BY COMMON
C           NONE
C
C******************************************************************************

      IMPLICIT NONE

C  Include Electron Energy Deposition Parameters

      INCLUDE 'ED_R_elec_ed_lup_subs.inc'

C  DEFINE LOCAL VARIABLES

      integer ENG_INDEX
      integer IOERR
      integer UNIT
      real    R_EDEP_ENG_VALUE
      real    UFACT

C Define common blocks

      INTEGER   LOADED
      real    PMAT(ENERGY_LEVELS,ALTITUDES)
      real    ENGMAT(ENERGY_LEVELS)
      real    EDEP_ALTS(ALTITUDES)
      real    EDEP_PRES(ALTITUDES)
      COMMON /R_EDEP/ LOADED,PMAT,ENGMAT,EDEP_ALTS,EDEP_PRES

C******************************************************************************
C
C  ENERGY_LEVELS     PARAMETER  NUMBER OF ENERGY LEVELS
C  ALTITUDES         PARAMETER  NUMBER OF ALTITUDE LEVELS
C  ENG_INDEX         integer  REQUESTED ENERGY BIN
C  IOERR             integer  STATUS OF I/O OPERATION
C  UNIT              integer  UNIT FOR ENERGY
C                               1 = erg
C                               2 = Joule
C                               3 = eV
C                               4 = keV
C                               5 = MeV
C  R_EDEP_ENG_VALUE  real     ENERGY OF THE INDEX REQUESTED
C  UFACT             real     UNIT FACTOR FOR CONVERSION INTO [MeV]
C                               UFACT = 1.602E-6 FOR erg
C                               UFACT = 1.602E-13 FOR Joule
C                               UFACT = 1.0E+6 FOR eV
C                               UFACT = 1.0E+3 FOR keV
C                               UFACT = 1.0 FOR MeV
C  LOADED            LOGICAL*4  /R_EDEP/  LOGICAL INDICATING STATE OF DATA FILE
C                                         1  FILE HAS BEEN LOADED
C                                         0 FILE HAS NOT BEEN LOADED
C  PMAT(ENERGY_LEVELS,ALTITUDES)   
C                    real     /R_EDEP/  ENERGY DEPOSITION VALUES FOR INCIDENT
C                                         ISOTROPIC FLUX [eV/cm]
C  ENGMAT(ENERGY_LEVELS)           
C                    real     /R_EDEP/  ENERGY VALUES [MeV]
C  EDEP_ALTS(ALTITUDES)            
C                    real     /R_EDEP/  ALTITUDE VALUES USED BY THE ENERGY 
C                                         DEPOSITION PROFILES [km]
C  EDEP_PRES(ALTITUDES)            
C                    real     /R_EDEP/  PRESSURE VALUES USED BY THE ENERGY 
C                                         DEPOSITION PROFILES [mbar]
C
C******************************************************************************

C  Set indicator indicating off the energy grid

      R_EDEP_ENG_VALUE = -1.0

C  Have the energy deposition tables been loaded?

      IF (LOADED.eq.0) THEN
        WRITE (*,'(1X,A55,I5)',IOSTAT=IOERR)
     +         'SKIP: Energy Deposition Files Not Loaded.',IOERR
        RETURN
      ENDIF

C  Set the unit factor

      UFACT = 1.0
      IF (MOD(UNIT,5) .EQ. 1) UFACT = 1.602E-6
      IF (MOD(UNIT,5) .EQ. 2) UFACT = 1.602E-13
      IF (MOD(UNIT,5) .EQ. 3) UFACT = 1.0E+6
      IF (MOD(UNIT,5) .EQ. 4) UFACT = 1.0E+3

C  Signify bad value if energy is above or below altitude grid

      IF (ENG_INDEX .LT. 1) RETURN
      IF (ENG_INDEX .GT. ENERGY_LEVELS) RETURN
      
C  Determine value of altitude and shift to requested units

      R_EDEP_ENG_VALUE = UFACT*ENGMAT(ENG_INDEX)

      RETURN
      END

C------------------------------------------------------------------------------

      SUBROUTINE R_ELEC_EDEP (FLUXMAT,FUNIT,ENRGMAT,EUNIT,PROFMAT,PUNIT)

C******************************************************************************
C
C  THIS FUNCTION COMPUTES THE ENERGY DEPOSITION PROFILE FROM AN INPUT 
C  ENERGY-NUMBER FLUX SPECTRA IN UNITS REQUESTED
C
C  WHEN THE sr UNIT IS INCLUDED IN THE INPUT, THE RESULT IS MULTIPLIED BY
C  PI TO ACCOUNT FOR AN ISOTROPIC PITCH ANGLE DISTRIBUTION.  IF FLUX UNITS DO
C  NOT HAVE sr, THEN IT IS ASSUMED THAT THE INPUT SPECTRA IS ALREADY AN 
C  ISOTROPIC FLUX.  ENERGIES STORED ARE IN [MeV].  THE PROFILES ARE COMPUTED 
C  FOR A FLUX UNIT OF #/(cm**2-sr-s-eV) AT EACH ALTITUDE IN [km].  A SUM OF
C  THE INDIVIDUAL PROFILES IS THE TOTAL ENERGY DEPOSITION RATE.
C
C  INPUTS - DIRECTLY PASSED
C           FLUXMAT(SPECTRA_LEVELS) real     FLUX VALUES
C           FUNIT                  integer   UNIT OF FLUX
C                                              1 = #/(km**2-s-eV)
C                                              2 = #/(km**2-s-keV)
C                                              3 = #/(km**2-s-MeV)
C                                              4 = #/(km**2-s-Joule)
C                                              5 = #/(km**2-s-erg)
C                                              6 = #/(m**2-s-eV)
C                                              7 = #/(m**2-s-keV)
C                                              8 = #/(m**2-s-MeV)
C                                              9 = #/(m**2-s-Joule)
C                                              10 = #/(m**2-s-erg)
C                                              11 = #/(cm**2-s-eV)
C                                              12 = #/(cm**2-s-keV)
C                                              13 = #/(cm**2-s-MeV)
C                                              14 = #/(cm**2-s-Joule)
C                                              15 = #/(cm**2-s-erg)
C                                              16 = #/(km**2-s-sr-erg)
C                                              17 = #/(km**2-s-sr-Joule)
C                                              18 = #/(km**2-s-sr-MeV)
C                                              19 = #/(km**2-s-sr-keV)
C                                              20 = #/(km**2-s-sr-eV)
C                                              21 = #/(m**2-s-sr-erg)
C                                              22 = #/(m**2-s-sr-Joule)
C                                              23 = #/(m**2-s-sr-MeV)
C                                              24 = #/(m**2-s-sr-keV)
C                                              25 = #/(m**2-s-sr-eV)
C                                              26 = #/(cm**2-s-sr-erg)
C                                              27 = #/(cm**2-s-sr-Joule)
C                                              28 = #/(cm**2-s-sr-MeV)
C                                              29 = #/(cm**2-s-sr-keV)
C                                              30 = #/(cm**2-s-sr-eV)
C           ENRGMAT(SPECTRA_LEVELS) real     ENERGY VALUES
C           EUNIT                  integer   UNIT OF ENERGY
C                                              1 = erg
C                                              2 = Joule
C                                              3 = eV
C                                              4 = keV
C                                              5 = MeV
C           PUNIT                  integer   REQUESTED UNIT OF PROFILE
C                                              1 = Ionizations/(km**3-s)
C                                              2 = eV/(km**3-s)
C                                              3 = keV/(km**3-s)
C                                              4 = MeV/(km**3-s)
C                                              5 = Joule/(km**3-s)
C                                              6 = erg/(km**3-s)
C                                              7 = Ionizations/(m**3-s)
C                                              8 = eV/(m**3-s)
C                                              9 = keV/(m**3-s)
C                                              10 = MeV/(m**3-s)
C                                              11 = Joule/(m**3-s)
C                                              12 = erg/(m**3-s)
C                                              13 = Ionizations/(cm**3-s)
C                                              14 = eV/(cm**3-s)
C                                              15 = keV/(cm**3-s)
C                                              16 = MeV/(cm**3-s)
C                                              17 = Joule/(cm**3-s)
C                                              18 = erg/(cm**3-s)
C
C  INPUTS - PASSED BY NAMED COMMON /R_EDEP/
C           ENGMAT(ENERGY_LEVELS) real      ENERGY VALUES USED BY THE ENERGY
C                                             DEPOSITION PROFILES [MeV]
C
C  OUTPUTS - DIRECTLY PASSED
C           PROFMAT(ALTITUDES)    real      ENERGY DEPOSITION PROFILE
C
C  OUTPUTS - PASSED BY COMMON
C           NONE
C
C******************************************************************************

      IMPLICIT NONE

C  Include Electron Energy Deposition Parameters

      INCLUDE 'ED_R_elec_ed_lup_subs.inc'

C  DEFINE LOCAL VARIABLES

      integer IOERR
      integer FUNIT
      integer EUNIT
      integer PUNIT
      integer I,J,K
      real    FUFACT
      real    EUFACT
      real    PUFACT
      real    PI
      real    AFACT
      real    ADJUSTMENT
      real    R_LOG_INTERP
      real    FLUXMAT(SPECTRA_LEVELS)
      real    ENRGMAT(SPECTRA_LEVELS)
      real    FMAT(SPECTRA_LEVELS)
      real    EMAT(SPECTRA_LEVELS)
      real    PROFMAT(ALTITUDES)
      real    DEMAT(SPECTRA_LEVELS)

C Define common blocks

      INTEGER   LOADED
      real    PMAT(ENERGY_LEVELS,ALTITUDES)
      real    ENGMAT(ENERGY_LEVELS)
      real    EDEP_ALTS(ALTITUDES)
      real    EDEP_PRES(ALTITUDES)
      COMMON /R_EDEP/ LOADED,PMAT,ENGMAT,EDEP_ALTS,EDEP_PRES

C Define data values

      DATA PI / 3.14159826 /

C******************************************************************************
C
C  SPECTRA_LEVELS    PARAMETER  NUMBER OF ENERGY LEVELS IN USER SPECTRUM
C  ENERGY_LEVELS     PARAMETER  NUMBER OF ENERGY LEVELS
C  ALTITUDES         PARAMETER  NUMBER OF ALTITUDE LEVELS
C  IOERR             integer  STATUS OF I/O OPERATION
C  FUNIT             integer  UNIT OF FLUX
C                               1 = #/(km**2-s-eV)
C                               2 = #/(km**2-s-keV)
C                               3 = #/(km**2-s-MeV)
C                               4 = #/(km**2-s-Joule)
C                               5 = #/(km**2-s-erg)
C                               6 = #/(m**2-s-eV)
C                               7 = #/(m**2-s-keV)
C                               8 = #/(m**2-s-MeV)
C                               9 = #/(m**2-s-Joule)
C                               10 = #/(m**2-s-erg)
C                               11 = #/(cm**2-s-eV)
C                               12 = #/(cm**2-s-keV)
C                               13 = #/(cm**2-s-MeV)
C                               14 = #/(cm**2-s-Joule)
C                               15 = #/(cm**2-s-erg)
C                               16 = #/(km**2-s-sr-erg)
C                               17 = #/(km**2-s-sr-Joule)
C                               18 = #/(km**2-s-sr-MeV)
C                               19 = #/(km**2-s-sr-keV)
C                               20 = #/(km**2-s-sr-eV)
C                               21 = #/(m**2-s-sr-erg)
C                               22 = #/(m**2-s-sr-Joule)
C                               23 = #/(m**2-s-sr-MeV)
C                               24 = #/(m**2-s-sr-keV)
C                               25 = #/(m**2-s-sr-eV)
C                               26 = #/(cm**2-s-sr-erg)
C                               27 = #/(cm**2-s-sr-Joule)
C                               28 = #/(cm**2-s-sr-MeV)
C                               29 = #/(cm**2-s-sr-keV)
C                               30 = #/(cm**2-s-sr-eV)
C  EUNIT             integer  UNIT OF ENERGY
C                               1 = erg
C                               2 = Joule
C                               3 = eV
C                               4 = keV
C                               5 = MeV
C  PUNIT             integer  REQUESTED UNIT FOR PROFILE
C                               1 = Ionizations/(km**3-s)
C                               2 = eV/(km**3-s)
C                               3 = keV/(km**3-s)
C                               4 = MeV/(km**3-s)
C                               5 = Joule/(km**3-s)
C                               6 = erg/(km**3-s)
C                               7 = Ionizations/(m**3-s)
C                               8 = eV/(m**3-s)
C                               9 = keV/(m**3-s)
C                               10 = MeV/(m**3-s)
C                               11 = Joule/(m**3-s)
C                               12 = erg/(m**3-s)
C                               13 = Ionizations/(cm**3-s)
C                               14 = eV/(cm**3-s)
C                               15 = keV/(cm**3-s)
C                               16 = MeV/(cm**3-s)
C                               17 = Joule/(cm**3-s)
C                               18 = erg/(cm**3-s)
C  I                 integer  INDEX COUNTER
C  J                 integer  INDEX COUNTER
C  K                 integer  INDEX COUNTER
C  FUFACT            real     UNIT FACTOR FOR CONVERSION INTO PER eV UNIT
C                               FUFACT = 1.602E-22 FOR #/(km**2-s-erg)
C                               FUFACT = 1.602E-29 FOR #/(km**2-s-Joule)
C                               FUFACT = 1.0E-16 FOR #/(km**2-s-MeV)
C                               FUFACT = 1.0E-13 FOR #/(km**2-s-keV)
C                               FUFACT = 1.0E-10 FOR #/(km**2-s-eV)
C                               FUFACT = 1.602E-16 FOR #/(m**2-s-erg)
C                               FUFACT = 1.602E-23 FOR #/(m**2-s-Joule)
C                               FUFACT = 1.0E-10 FOR #/(m**2-s-MeV)
C                               FUFACT = 1.0E-7 FOR #/(m**2-s-keV)
C                               FUFACT = 1.0E-4 FOR #/(m**2-s-eV)
C                               FUFACT = 1.602E-12 FOR #/(cm**2-s-erg)
C                               FUFACT = 1.602E-19 FOR #/(cm**2-s-Joule)
C                               FUFACT = 1.0E-6 FOR #/(cm**2-s-MeV)
C                               FUFACT = 1.0E-3 FOR #/(cm**2-s-keV)
C                               FUFACT = 1.0 FOR #/(cm**2-s-eV)
C                               FUFACT = 1.602E-22 FOR #/(km**2-s-sr-erg)
C                               FUFACT = 1.602E-29 FOR #/(km**2-s-sr-Joule)
C                               FUFACT = 1.0E-16 FOR #/(km**2-s-sr-MeV)
C                               FUFACT = 1.0E-13 FOR #/(km**2-s-sr-keV)
C                               FUFACT = 1.0E-10 FOR #/(km**2-s-sr-eV)
C                               FUFACT = 1.602E-16 FOR #/(m**2-s-sr-erg)
C                               FUFACT = 1.602E-23 FOR #/(m**2-s-sr-Joule)
C                               FUFACT = 1.0E-10 FOR #/(m**2-s-sr-MeV)
C                               FUFACT = 1.0E-7 FOR #/(m**2-s-sr-keV)
C                               FUFACT = 1.0E-4 FOR #/(m**2-s-sr-eV)
C                               FUFACT = 1.602E-12 FOR #/(cm**2-s-sr-erg)
C                               FUFACT = 1.602E-19 FOR #/(cm**2-s-sr-Joule)
C                               FUFACT = 1.0E-6 FOR #/(cm**2-s-sr-MeV)
C                               FUFACT = 1.0E-3 FOR #/(cm**2-s-sr-keV)
C                               FUFACT = 1.0 FOR #/(cm**2-s-sr-eV)
C  EUFACT            real     UNIT FACTOR FOR CONVERSION INTO [MeV]
C                               EUFACT = 6.242E+5 FOR erg
C                               EUFACT = 6.242E+12 FOR Joule
C                               EUFACT = 1.0E-6 FOR eV
C                               EUFACT = 1.0E-3 FOR keV
C                               EUFACT = 1.0 FOR MeV
C  PUFACT            real     UNIT FACTOR FOR CONVERSION INTO REQUESTED UNITS
C                               PUFACT = 1.783E+25 FOR Ionizations/(km**3-s)
C                               PUFACT = 6.242E+26 FOR eV/(km**3-s)
C                               PUFACT = 6.242E+23 FOR keV/(km**3-s)
C                               PUFACT = 6.242E+20 FOR MeV/(km**3-s)
C                               PUFACT = 1.0E+8 FOR Joule/(km**3-s)
C                               PUFACT = 1.0E+15 FOR erg/(km**3-s)
C                               PUFACT = 1.783E+16 FOR Ionizations/(m**3-s)
C                               PUFACT = 6.242E+17 FOR eV/(m**3-s)
C                               PUFACT = 6.242E+14 FOR keV/(m**3-s)
C                               PUFACT = 6.242E+11 FOR MeV/(m**3-s)
C                               PUFACT = 1.0E-1 FOR Joule/(m**3-s)
C                               PUFACT = 1.0E+6 FOR erg/(m**3-s)
C                               PUFACT = 1.783E+10 FOR Ionizations/(cm**3-s)
C                               PUFACT = 6.242E+11 FOR eV/(cm**3-s)
C                               PUFACT = 6.242E+8 FOR keV/(cm**3-s)
C                               PUFACT = 6.242E+5 FOR MeV/(cm**3-s)
C                               PUFACT = 1.0E-7 FOR Joule/(cm**3-s)
C                               PUFACT = 1.0 FOR erg/(cm**3-s)
C  PI                real     NUMBER pi = 3.14159826
C  R_LOG_INTERP      real     FUNCTION FOR LOG INTERPOLATION
C  AFACT             real     ANGULAR FACTOR
C                               AFACT = 1.0 FOR sr INCLUDED IN FLUX
C                               AFACT = PI FOR sr NOT INCLUDED IN FLUX
C  ADJUSTMENT        real     INTERPOLATED PROFILE VALUE
C  FLUXMAT(SPECTRA_LEVELS)
C                    real     ARRAY HOLDING FLUX VALUES
C  ENRGMAT(SPECTRA_LEVELS)
C                    real     ARRAY HOLDING INPUT ENERGY LEVELS
C  FMAT(SPECTRA_LEVELS)
C                    real     INTERNAL ARRAY HOLDING FLUX VALUES 
C                               [#/(cm**2-sr-s-eV) OR #/(cm**2-s-eV)]
C  EMAT(SPECTRA_LEVELS)
C                    real     INTERNAL ARRAY HOLDING INPUT ENERGY LEVELS[MeV]
C  DEMAT(SPECTRA_LEVELS)
C                    real     WIDTH OF EACH DESCRETE SPECTRAL ENERGY LEVEL
C  PROFMAT(ALTITUDES)
C  LOADED            LOGICAL*4  /R_EDEP/  LOGICAL INDICATING STATE OF DATA FILE
C                                         1  FILE HAS BEEN LOADED
C                                         0 FILE HAS NOT BEEN LOADED
C  PMAT(ENERGY_LEVELS,ALTITUDES)   
C                    real     /R_EDEP/  ENERGY DEPOSITION VALUES FOR INCIDENT
C                                         ISOTROPIC FLUX [eV/cm]
C  ENGMAT(ENERGY_LEVELS)           
C                    real     /R_EDEP/  ENERGY VALUES [MeV]
C  EDEP_ALTS(ALTITUDES)            
C                    real     /R_EDEP/  ALTITUDE VALUES USED BY THE ENERGY 
C                                         DEPOSITION PROFILES [km]
C  EDEP_PRES(ALTITUDES)            
C                    real     /R_EDEP/  PRESSURE VALUES USED BY THE ENERGY 
C                                         DEPOSITION PROFILES [mbar]
C
C******************************************************************************

C Set the output profile all zeros

      DO I=1,ALTITUDES
        PROFMAT(I) = 0.0
      ENDDO

C  Have the energy deposition tables been loaded?

      IF (LOADED.eq.0) THEN
        WRITE (*,'(1X,A55,I5)',IOSTAT=IOERR)
     +         'SKIP: Energy Deposition Files Not Loaded.',IOERR
        RETURN
      ENDIF

C  Set the flux unit factor and angular integral factor for isotropy

C  Set the energy unit factor

      EUFACT = 1.0
      I = MOD(EUNIT,5)
      IF (I .EQ. 1) EUFACT = 6.242E+5
      IF (I .EQ. 2) EUFACT = 6.242E+12
      IF (I .EQ. 3) EUFACT = 1.0E-6
      IF (I .EQ. 4) EUFACT = 1.0E-3

C  Adjust energy for internal use to MeV

      DO I=1,SPECTRA_LEVELS
        EMAT(I) = EUFACT*ENRGMAT(I)
      ENDDO

C  Set the angular multiplication factor for isotropic angular distribution

      AFACT = PI
      I = MOD(FUNIT,30)
      IF ((I .GE. 1) .AND. (I .LE. 15)) AFACT = 1.0

C  Set Flux conversion to #/(cm**2-s-sr-eV) or #/(cm**2-s-eV)

      FUFACT = 1.0
      I = MOD(FUNIT,15)
      IF (I .EQ. 1) EUFACT = 1.602E-22
      IF (I .EQ. 2) FUFACT = 1.602E-29
      IF (I .EQ. 3) FUFACT = 1.0E-16
      IF (I .EQ. 4) FUFACT = 1.0E-13
      IF (I .EQ. 5) FUFACT = 1.0E-10
      IF (I .EQ. 6) FUFACT = 1.602E-16
      IF (I .EQ. 7) FUFACT = 1.602E-23
      IF (I .EQ. 8) FUFACT = 1.0E-10
      IF (I .EQ. 9) FUFACT = 1.0E-7
      IF (I .EQ. 10) FUFACT = 1.0E-4
      IF (I .EQ. 11) FUFACT = 1.602E-12
      IF (I .EQ. 12) FUFACT = 1.602E-19
      IF (I .EQ. 13) FUFACT = 1.0E-6
      IF (I .EQ. 14) FUFACT = 1.0E-3

C  Adjust flux for internal use

      DO I=1,SPECTRA_LEVELS
        FMAT(I) = FUFACT*FLUXMAT(I)
      ENDDO

C  Set output profile conversion factor

      PUFACT = 1.0
      I = MOD(PUNIT,18)
      IF (I .EQ. 1) PUFACT = 1.783485E+25
      IF (I .EQ. 2) PUFACT = 6.242197E+26
      IF (I .EQ. 3) PUFACT = 6.242197E+23
      IF (I .EQ. 4) PUFACT = 6.242197E+20
      IF (I .EQ. 5) PUFACT = 1.0E+8
      IF (I .EQ. 6) PUFACT = 1.0E+15
      IF (I .EQ. 7) PUFACT = 1.783485E+16
      IF (I .EQ. 8) PUFACT = 6.242197E+17
      IF (I .EQ. 9) PUFACT = 6.242197E+14
      IF (I .EQ. 10) PUFACT = 6.242197E+11
      IF (I .EQ. 11) PUFACT = 1.0E-1
      IF (I .EQ. 12) PUFACT = 1.0E+6
      IF (I .EQ. 13) PUFACT = 1.783485E+10
      IF (I .EQ. 14) PUFACT = 6.242197E+11
      IF (I .EQ. 15) PUFACT = 6.242197E+8
      IF (I .EQ. 16) PUFACT = 6.242197E+5
      IF (I .EQ. 17) PUFACT = 1.0E-7

C  Compute the spectral channel widths

      CALL R_FIND_DE (EMAT,DEMAT)

C  Adjust energy widths from MeV to eV

      DO I=1,SPECTRA_LEVELS
        DEMAT(I) = 1.0E6*DEMAT(I)
      ENDDO

C  Look through the spectra and collect individual profiles

      DO I=1,SPECTRA_LEVELS

C  Process for valid energies only

        IF (EMAT(I) .GT. 0.0) THEN

C  Handle the first profile separately

          IF (EMAT(I) .EQ. ENGMAT(1)) THEN
            DO K=1,ALTITUDES
              PROFMAT(K) = PROFMAT(K) + FMAT(1)*PMAT(1,K)*DEMAT(1)
            ENDDO

C  Check to see if the energy is greater than the maximum energy profile

          ELSE IF (EMAT(I) .GT. ENGMAT(1)) THEN

C  Extrapolate a profile

            DO K=1,ALTITUDES
              ADJUSTMENT = R_LOG_INTERP (EMAT(I),ENGMAT(1),ENGMAT(2),
     1                                            PMAT(1,K),PMAT(2,K))
              PROFMAT(K) = PROFMAT(K) + FMAT(I)*ADJUSTMENT*DEMAT(I)
            ENDDO

C  Check to see if the energy is less than the minimum energy profile

          ELSE IF (EMAT(I) .LT. ENGMAT(ENERGY_LEVELS)) THEN

C  Extrapolate a profile

            DO K=1,ALTITUDES
              ADJUSTMENT = R_LOG_INTERP (EMAT(I),ENGMAT(ENERGY_LEVELS),
     1                  ENGMAT(ENERGY_LEVELS-1),PMAT(ENERGY_LEVELS,K),
     2                                        PMAT(ENERGY_LEVELS-1,K))
              PROFMAT(K) = PROFMAT(K) + FMAT(I)*ADJUSTMENT*DEMAT(I)
            ENDDO
          ELSE

C  The requested energy lies in the profile grid
C  Find the spectral energy in the energy lookup table
C  Loop through the remainder of the energy levels to find the appropriate 
C  profile(s)

            DO J=2,ENERGY_LEVELS

C  Handle cases where the energy is exact -- no interpolation needed

              IF (EMAT(I) .EQ. ENGMAT(J)) THEN
                DO K=1,ALTITUDES
                  PROFMAT(K) = PROFMAT(K) + FMAT(I)*PMAT(J,K)*DEMAT(I)
                ENDDO

              ELSE IF ((EMAT(I) .GT. ENGMAT(J)) .AND. (EMAT(I) .LT. 
     1                                              ENGMAT(J-1))) THEN

C  Energy requested lies between precomputed profiles - interpolation is needed

                DO K=1,ALTITUDES
                  ADJUSTMENT = R_LOG_INTERP (EMAT(I),ENGMAT(J-1),
     1                                 ENGMAT(J),PMAT(J-1,K),PMAT(J,K))

                  PROFMAT(K) = PROFMAT(K) + FMAT(I)*ADJUSTMENT*DEMAT(I)

                ENDDO

              ENDIF
            ENDDO  ! INDIVIDUAL PROFILES
          ENDIF  ! TYPE OF PROFILE
        ENDIF  ! VALLID ENERGIES
      ENDDO  ! SPECTRAL VALUES

C  Adjust total profile by angular distribution factor and 
C  convert to output units

      DO I=1,ALTITUDES
        PROFMAT(I) = PUFACT*AFACT*PROFMAT(I)
      ENDDO

      RETURN
      END

C------------------------------------------------------------------------------

      FUNCTION R_LOG_INTERP (X,X2,X1,Y2,Y1)

C******************************************************************************
C
C  THIS FUNCTION DETERMINES AN INDEPENDENT VALUE AT A GIVEN DEPENDENT VALUE 
C  LOCATION GIVEN TWO COORDINATE POINTS.  THE VALUE OF THE INDEPENDENT VARIABLE
C  IS DETERMINED IN LOG-LOG SPACE.  ALL VALUES ARE GIVEN LINEARLY AND THIS
C  FUNCTION DETERMINES THE LOG VALUES AS NEEDED.  IF ANY ARGUMENT IS LESS
C  THAN OR EQUAL TO ZERO, LINEAR INTERPOLATION IS INVOLKED.  IF THE GIVEN
C  X1 AND X2 CORRDINATES ARE THE SAME, THEN THE OUTPUT IS THE LOG AVERAGE OF
C  THE Y1 AND Y2 POINT REGARDLESS ON THE DEPENDENT VARIABLE.
C
C  INPUTS - DIRECTLY PASSED
C           X                     real      REQUEST AT THIS INDEAENDENT VALUE
C           X1                    real      KNOWN X COORDINATE OF POINT 1
C           X2                    real      KNOWN X COORDINATE OF POINT 2
C           Y1                    real      KNOWN Y COORDINATE OF POINT 1
C           Y2                    real      KNOWN Y COORDINATE OF POINT 2
C
C  INPUTS - PASSED BY COMMON
C           NONE
C
C  OUTPUTS - DIRECTLY PASSED
C           R_LOG_INTERP          real      THE DEPENDENT VALUE
C
C  OUTPUTS - PASSED BY COMMON
C           NONE
C
C******************************************************************************

      real    X,X1,X2,Y1,Y2
      real    U,U1,U2
      real    V,V1,V2
      real    R_LOG_INTERP
      real    R_LINEAR_INTERP
      real    SLOPE
      real    INTERCEPT

C******************************************************************************
C
C  X                 real     INTERPOLATION REQUESTED AT THIS DEPENDENT VALUE
C  X1                real     KNOWN LINEAR X COORDINATE OF POINT 1
C  X2                real     KNOWN LINEAR X COORDINATE OF POINT 2
C  Y1                real     KNOWN LINEAR Y COORDINATE OF POINT 1
C  Y2                real     KNOWN LINEAR Y COORDINATE OF POINT 2
C  U                 real     LOG REQUESTED DEPENDENT VALUE
C  U1                real     KNOWN LOG X COORDINATE OF POINT 1
C  U2                real     KNOWN LOG X COORDINATE OF POINT 2
C  V                 real     LOG DEPENDENT VALUE
C  V1                real     KNOWN LOG Y COORDINATE OF POINT 1
C  V2                real     KNOWN LOG Y COORDINATE OF POINT 2
C  R_LOG_INTERP      real     FUNCTION FOR LOG INTERPOLATION
C  R_LINEAR_INTERP   real     FUNCTION FOR LINEAR INTERPOLATION
C  SLOPE             real     Y VALUE CHANGE WITH RESPECT TO X VALUE CHANGE
C  INTERCEPT         real     Y VALUE LOCATION AT AT X OF ZERO
C
C******************************************************************************

C  All values must be non=zero and positive, if not, you can not use this 
C  function -- Linearly Interpolate

      IF ((X .LE. 0.0) .OR. (X1 .LE. 0.0) .OR. (X2 .LE. 0.0)) THEN
        R_LOG_INTERP = R_LINEAR_INTERP (X,X2,X1,Y2,Y1)
      ELSE IF ((Y1 .LE. 0.0) .OR. (Y2 .LE. 0.0)) THEN
        R_LOG_INTERP = R_LINEAR_INTERP (X,X2,X1,Y2,Y1)
      ELSE

C  Determine the log values and degeneracy

        V1 = LOG10(Y1)
        V2 = LOG10(Y2)
        U = LOG10(X)
        U1 = LOG10(X1)
        U2 = LOG10(X2)
        SLOPE = U2 - U1

C  Handle the case of degenerate dependent values separately

        IF (SLOPE .EQ. 0.0) THEN
          V = 0.5*(V2 + V1)
          R_LOG_INTERP = 10.0**V
        ELSE

C  Interpolate in log-log space

          SLOPE = (V2 - V1)/SLOPE
          INTERCEPT = V2 - SLOPE*U2
          V = SLOPE*U + INTERCEPT
          R_LOG_INTERP = 10.0**V
        ENDIF
      ENDIF

      RETURN
      END

C------------------------------------------------------------------------------

      FUNCTION R_SEMILOG_INTERP (X,X2,X1,Y2,Y1)

C******************************************************************************
C
C  THIS FUNCTION DETERMINES AN INDEPENDENT VALUE AT A GIVEN DEPENDENT VALUE 
C  LOCATION GIVEN TWO COORDINATE POINTS.  THE VALUE OF THE INDEPENDENT VARIABLE
C  IS DETERMINED IN LOG-LINEAR SPACE.  ALL VALUES ARE GIVEN LINEARLY AND THIS
C  FUNCTION DETERMINES THE LOG VALUES AS NEEDED.  IF ANY X ARGUMENT IS LESS
C  THAN OR EQUAL TO ZERO, LINEAR INTERPOLATION IS INVOLKED.  IF THE GIVEN
C  X1 AND X2 CORRDINATES ARE THE SAME, THEN THE OUTPUT IS THE LINEAR AVERAGE OF
C  THE Y1 AND Y2 POINT REGARDLESS ON THE DEPENDENT VARIABLE.
C
C  INPUTS - DIRECTLY PASSED
C           X                     real      REQUEST AT THIS INDEAENDENT VALUE
C           X1                    real      KNOWN X COORDINATE OF POINT 1
C           X2                    real      KNOWN X COORDINATE OF POINT 2
C           Y1                    real      KNOWN Y COORDINATE OF POINT 1
C           Y2                    real      KNOWN Y COORDINATE OF POINT 2
C
C  INPUTS - PASSED BY COMMON
C           NONE
C
C  OUTPUTS - DIRECTLY PASSED
C           R_SEMILOG_INTERP      real      THE DEPENDENT VALUE
C
C  OUTPUTS - PASSED BY COMMON
C           NONE
C
C******************************************************************************

      real    X,X1,X2,Y1,Y2
      real    U,U1,U2
      real    R_SEMILOG_INTERP
      real    R_LINEAR_INTERP
      real    SLOPE
      real    INTERCEPT

C******************************************************************************
C
C  X                 real     INTERPOLATION REQUESTED AT THIS DEPENDENT VALUE
C  X1                real     KNOWN LINEAR X COORDINATE OF POINT 1
C  X2                real     KNOWN LINEAR X COORDINATE OF POINT 2
C  Y1                real     KNOWN LINEAR Y COORDINATE OF POINT 1
C  Y2                real     KNOWN LINEAR Y COORDINATE OF POINT 2
C  U                 real     LOG REQUESTED DEPENDENT VALUE
C  U1                real     KNOWN LOG X COORDINATE OF POINT 1
C  U2                real     KNOWN LOG X COORDINATE OF POINT 2
C  R_SEMILOG_INTERP  real     FUNCTION FOR SEMILOG INTERPOLATION
C  R_LINEAR_INTERP   real     FUNCTION FOR LINEAR INTERPOLATION
C  SLOPE             real     Y VALUE CHANGE WITH RESPECT TO X VALUE CHANGE
C  INTERCEPT         real     Y VALUE LOCATION AT AT X OF ZERO
C
C******************************************************************************

C  All values must be non=zero and positive, if not, you can not use this 
C  function -- Linearly Interpolate

      IF ((X .LE. 0.0) .OR. (X1 .LE. 0.0) .OR. (X2 .LE. 0.0)) THEN
        R_SEMILOG_INTERP = R_LINEAR_INTERP (X,X2,X1,Y2,Y1)
      ELSE

C  Determine the log values and degeneracy

        U = LOG10(X)
        U1 = LOG10(X1)
        U2 = LOG10(X2)
        SLOPE = U2 - U1

C  Handle the case of degenerate dependent values separately

        IF (SLOPE .EQ. 0.0) THEN
          R_SEMILOG_INTERP = 0.5*(Y2 + Y1)
        ELSE

C  Interpolate in log-linear space

          SLOPE = (Y2 - Y1)/SLOPE
          INTERCEPT = Y2 - SLOPE*U2
          R_SEMILOG_INTERP = SLOPE*U + INTERCEPT
        ENDIF
      ENDIF

      RETURN
      END

C------------------------------------------------------------------------------

      FUNCTION R_SEMILINEAR_INTERP (X,X2,X1,Y2,Y1)

C******************************************************************************
C
C  THIS FUNCTION DETERMINES AN INDEPENDENT VALUE AT A GIVEN DEPENDENT VALUE 
C  LOCATION GIVEN TWO COORDINATE POINTS.  THE VALUE OF THE INDEPENDENT VARIABLE
C  IS DETERMINED IN LINEAR-LOG SPACE.  ALL VALUES ARE GIVEN LINEARLY AND THIS
C  FUNCTION DETERMINES THE LOG VALUES AS NEEDED.  IF ANY Y ARGUMENT IS LESS
C  THAN OR EQUAL TO ZERO, LINEAR INTERPOLATION IS INVOLKED.  IF THE GIVEN
C  X1 AND X2 CORRDINATES ARE THE SAME, THEN THE OUTPUT IS THE LOG AVERAGE OF
C  THE Y1 AND Y2 POINT REGARDLESS ON THE DEPENDENT VARIABLE.
C
C  INPUTS - DIRECTLY PASSED
C           X                     real      REQUEST AT THIS INDEAENDENT VALUE
C           X1                    real      KNOWN X COORDINATE OF POINT 1
C           X2                    real      KNOWN X COORDINATE OF POINT 2
C           Y1                    real      KNOWN Y COORDINATE OF POINT 1
C           Y2                    real      KNOWN Y COORDINATE OF POINT 2
C
C  INPUTS - PASSED BY COMMON
C           NONE
C
C  OUTPUTS - DIRECTLY PASSED
C           R_SEMILINEAR_INTERP   real      THE DEPENDENT VALUE
C
C  OUTPUTS - PASSED BY COMMON
C           NONE
C
C******************************************************************************

      real    X,X1,X2,Y1,Y2
      real    V,V1,V2
      real    R_SEMILINEAR_INTERP
      real    R_LINEAR_INTERP
      real    SLOPE
      real    INTERCEPT

C******************************************************************************
C
C  X                 real     INTERPOLATION REQUESTED AT THIS DEPENDENT VALUE
C  X1                real     KNOWN LINEAR X COORDINATE OF POINT 1
C  X2                real     KNOWN LINEAR X COORDINATE OF POINT 2
C  Y1                real     KNOWN LINEAR Y COORDINATE OF POINT 1
C  Y2                real     KNOWN LINEAR Y COORDINATE OF POINT 2
C  V                 real     LOG DEPENDENT VALUE
C  V1                real     KNOWN LOG Y COORDINATE OF POINT 1
C  V2                real     KNOWN LOG Y COORDINATE OF POINT 2
C  R_SEMILINEAR_INTERP real   FUNCTION FOR SEMILINEAR INTERPOLATION
C  R_LINEAR_INTERP   real     FUNCTION FOR LINEAR INTERPOLATION
C  SLOPE             real     Y VALUE CHANGE WITH RESPECT TO X VALUE CHANGE
C  INTERCEPT         real     Y VALUE LOCATION AT AT X OF ZERO
C
C******************************************************************************

C  All values must be non=zero and positive, if not, you can not use this 
C  function -- Linearly Interpolate

      IF ((Y1 .LE. 0.0) .OR. (Y2 .LE. 0.0)) THEN
        R_SEMILINEAR_INTERP = R_LINEAR_INTERP (X,X2,X1,Y2,Y1)
      ELSE

C  Determine the log values and degeneracy

        V1 = LOG10(Y1)
        V2 = LOG10(Y2)
        SLOPE = X2 - X1

C  Handle the case of degenerate dependent values separately

        IF (SLOPE .EQ. 0.0) THEN
          V = 0.5*(V2 + V1)
          R_SEMILINEAR_INTERP = 10.0**V
        ELSE

C  Interpolate in linear-log space

          SLOPE = (V2 - V1)/SLOPE
          INTERCEPT = V2 - SLOPE*X2
          V = SLOPE*X + INTERCEPT
          R_SEMILINEAR_INTERP = 10.0**V
        ENDIF
      ENDIF

      RETURN
      END

C------------------------------------------------------------------------------

      FUNCTION R_LINEAR_INTERP (X,X2,X1,Y2,Y1)

C******************************************************************************
C
C  THIS FUNCTION DETERMINES AN INDEPENDENT VALUE AT A GIVEN DEPENDENT VALUE 
C  LOCATION GIVEN TWO COORDINATE POINTS.  THE VALUE OF THE INDEPENDENT VARIABLE
C  IS DETERMINED IN LINEAR SPACE.  IF THE GIVEN X1 AND X2 CORRDINATES ARE THE 
C  SAME, THEN THE OUTPUT IS THE AVERAGE OF THE Y1 AND Y2 POINT REGARDLESS OF 
C  THE DEPENDENT VARIABLE.
C
C  INPUTS - DIRECTLY PASSED
C           X                     real      REQUEST AT THIS INDEAENDENT VALUE
C           X1                    real      KNOWN X COORDINATE OF POINT 1
C           X2                    real      KNOWN X COORDINATE OF POINT 2
C           Y1                    real      KNOWN Y COORDINATE OF POINT 1
C           Y2                    real      KNOWN Y COORDINATE OF POINT 2
C
C  INPUTS - PASSED BY COMMON
C           NONE
C
C  OUTPUTS - DIRECTLY PASSED
C           R_LINEAR_INTERP       real      THE DEPENDENT VALUE
C
C  OUTPUTS - PASSED BY COMMON
C           NONE
C
C******************************************************************************

      real    X,X1,X2,Y1,Y2
      real    R_LINEAR_INTERP
      real    SLOPE
      real    INTERCEPT

C******************************************************************************
C
C  X                 real     INTERPOLATION REQUESTED AT THIS DEPENDENT VALUE
C  X1                real     KNOWN LINEAR X COORDINATE OF POINT 1
C  X2                real     KNOWN LINEAR X COORDINATE OF POINT 2
C  Y1                real     KNOWN LINEAR Y COORDINATE OF POINT 1
C  Y2                real     KNOWN LINEAR Y COORDINATE OF POINT 2
C  R_LINEAR_INTERP   real     FUNCTION FOR LINEAR INTERPOLATION
C  SLOPE             real     Y VALUE CHANGE WITH RESPECT TO X VALUE CHANGE
C  INTERCEPT         real     Y VALUE LOCATION AT AT X OF ZERO
C
C******************************************************************************

C  Determine the degeneracy

      SLOPE = X2 - X1

C  Handle the case of degenerate dependent values separately

      IF (SLOPE .EQ. 0.0) THEN
        R_LINEAR_INTERP = 0.5*(Y2 + Y1)
      ELSE

C  Interpolate in linear space

        SLOPE = (Y2 - Y1)/SLOPE
        INTERCEPT = Y2 - SLOPE*X2
        R_LINEAR_INTERP = SLOPE*X + INTERCEPT
      ENDIF

      RETURN
      END

C------------------------------------------------------------------------------

      SUBROUTINE R_FIND_DE (EMAT,DEMAT)

C******************************************************************************
C
C  THIS SUBROUTINE DETERMINES THE CHANNEL WIDTH OF THE INPUT SPECTRA.  ENERGY
C  WIDTHS ARE RETAINED FOR FUTURE CALLS TO THIS ROUTINE.
C
C  INPUTS - DIRECTLY PASSED
C           EMAT(SPECTRA_LEVELS)   real      ENERGIES OF THE SPECTRA
C
C  INPUTS - PASSED BY NAMED COMMON /R_DE/
C  EMAT_SAVE(SPECTRA_LEVELS)           
C                    real     /R_DE/    PREVIOUS ENERGY
C  DEMAT_SAVE(SPECTRA_LEVELS)           
C                    real     /R_DE/    PREVIOUS ENERGY WIDTH FOR EACH ENERGY
C
C  OUTPUTS - DIRECTLY PASSED
C           DEMAT(SPECTRA_LEVELS)  real      ENERGY WIDTHS OF THE SPECTRA
C
C  OUTPUTS - PASSED BY NAMED COMMON /R_DE/
C  EMAT_SAVE(SPECTRA_LEVELS)           
C                    real     /R_DE/    CURRENT ENERGY
C  DEMAT_SAVE(SPECTRA_LEVELS)           
C                    real     /R_DE/    CURRENT ENERGY WIDTH FOR EACH ENERGY
C
C******************************************************************************

      IMPLICIT NONE

C  Include Electron Energy Deposition Parameters

      INCLUDE 'ED_R_elec_ed_lup_subs.inc'

C  DEFINE LOCAL VARIABLES

      INTEGER   DUPLICATE
      integer I,J,K
      integer FIRST_GOOD_POINT
      integer LAST_GOOD_POINT
      real    EMAT(SPECTRA_LEVELS)
      real    DEMAT(SPECTRA_LEVELS)
      real    EMAT_SAVE(SPECTRA_LEVELS)
      real    DEMAT_SAVE(SPECTRA_LEVELS)
      COMMON /R_DE/ EMAT_SAVE,DEMAT_SAVE

C******************************************************************************
C
C  DUPLICATE         LOGICAL*2  LOGICAL FOR TEST OF OLD ENERGY MATRIX
C                               1 MEANS THE ENERGY MATRIX HAS NOT CHANGED
C                               0 MEANS THAT THERE IS A NEW ENERGY MATRIX
C  I                 integer  INDEX COUNTER
C  J                 integer  INDEX COUNTER
C  K                 integer  INDEX COUNTER
C  FIRST_GOOD_POINT  integer  INDEX OF FIRST NON-ZERO ENERGY
C  LAST_GOOD_POINT   integer  INDEX OF LAST NON-ZERO ENERGY
C  FMAT(SPECTRA_LEVELS)
C                    real     INTERNAL ARRAY HOLDING FLUX VALUES 
C                               [#/(cm**2-sr-s-eV) OR #/(cm**2-s-eV)]
C  EMAT(SPECTRA_LEVELS)
C                    real     REQUESTED ENERGY CENTER VALUES
C  DEMAT(SPECTRA_LEVELS)
C                    real     WIDTH OF EACH DESCRETE SPECTRAL ENERGY LEVEL
C  EMAT_SAVE(SPECTRA_LEVELS)
C                    real     /R_DE/    PREVIOUS REQUESTED ENERGY VALUES
C  DEMAT_SAVE(SPECTRA_LEVELS)
C                    real     /R_DE/    PREVIOUS ENERGY WIDTHS
C
C******************************************************************************

C  Check to see if channel widths have already been computed

      DUPLICATE = 1
      DO I=1,SPECTRA_LEVELS
        IF (EMAT(I) .NE. EMAT_SAVE(I)) DUPLICATE = 0
      ENDDO

C  If you have computed them previously, just reload the previous values

      IF (DUPLICATE.eq.1) THEN

        DO I=1,SPECTRA_LEVELS
          DEMAT(I) = DEMAT_SAVE(I)
        ENDDO

C  If the previous energy sample has changed, you need to recompute 
C  channelwidths

      ELSE

C  Initialize energy width matrix

        DO I=1,spectra_LEVELS
          DEMAT(I) = 0.0
        ENDDO

C  Initialize search indicies

        J = 0
        K = 0
        FIRST_GOOD_POINT = 0
        LAST_GOOD_POINT = 0

C  Search for the first good point in the energy matrix

        DO I=spectra_LEVELS,1,-1
          IF (EMAT(I) .GT. 0.0) THEN
            K = J
            J = FIRST_GOOD_POINT
            FIRST_GOOD_POINT = I
          ENDIF
        ENDDO

C  No valid data values in energy spectrum

        IF (FIRST_GOOD_POINT .EQ. 0) RETURN

C  Only one valid data values in energy spectrum

        IF (J .EQ. 0) RETURN

C  Only two data points

        IF (K .EQ. 0) THEN
          I = FIRST_GOOD_POINT
          DEMAT(I) = ABS(EMAT(I) - EMAT(J))
          DEMAT(J) = DEMAT(I)
          RETURN
        ENDIF

C  Search for the last good point in the energy matrix

        DO I=1,SPECTRA_LEVELS
          IF (EMAT(I) .GT. 0.0) LAST_GOOD_POINT = I
        ENDDO

C  Compute the initial energy width

        i = 1
        j = 2
        k = 3

        DEMAT(I) = ABS(EMAT(I) - EMAT(J))

        i = 0
        j = 1
        k = 2

C  Do three point stepping ignoring zero's in energy matrix

  100   I = J
        J = K
  200   K = K + 1
        IF (K .GT. LAST_GOOD_POINT) go to 300
        IF (EMAT(K) .LE. 0.0) GO TO 200

C  Compute the energy width

        DEMAT(J) = ABS(0.5*(EMAT(I)+EMAT(J)) - 0.5*(EMAT(K)+EMAT(J)))

        GO TO 100

C  Compute the final energy width 

  300   DEMAT(J) = ABS(EMAT(I) - EMAT(J))

C  Save the newly generated channel widths

        DO I=1,SPECTRA_LEVELS
          DEMAT_SAVE(I) = DEMAT(I)
          EMAT_SAVE(I) = EMAT(I)
        ENDDO

      ENDIF

      RETURN
      END
