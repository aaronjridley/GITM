!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModPlanetConst

  use ModNumConst, ONLY: cDegToRad
  use ModConst, ONLY: cAU, cHour => cSecondPerHour, cDay => cSecondPerDay
  use ModKind

  implicit none

  save

  ! All astronomical bodies other than the Sun are defined below.  
  ! Solar constants are defined in ModConst.
  !
  ! The maximum number of astronomical bodies.  This is set at 200, it can be
  ! increased if necessary

  integer,parameter :: MaxPlanet = 200

  integer, parameter :: lNamePlanet = 40
  integer, parameter :: lTypeBField = 40


  ! Declarations for the variables that we are storing to define each body.
  ! 
  ! NOTE THAT ALL VARIABLES IN THIS FILE SHOULD BE IN SI UNITS. 
  !  (m,s,kg,m/s, ... )
  !
  ! NOTE THE THE PRECISE DEFINITIONS OF WHAT THE VARIABLES MEAN CAN BE
  ! FOUND AT THE END OF THE FILE AND IN THE CODE DOCUMENTATION (we hope).

  real,dimension(0:MaxPlanet+1):: rPlanet_I, mPlanet_I, rOrbitPlanet_I
  real,dimension(0:MaxPlanet+1):: OrbitalPeriodPlanet_I, RotationPeriodPlanet_I

  integer,dimension(0:MaxPlanet+1) :: &
       iYearEquinoxPlanet_I, iMonthEquinoxPlanet_I, iDayEquinoxPlanet_I, &
       iHourEquinoxPlanet_I,iMinuteEquinoxPlanet_I,iSecondEquinoxPlanet_I
  real,dimension(0:MaxPlanet+1)    :: FracSecondEquinoxPlanet_I
  real,dimension(0:MaxPlanet+1)    :: TiltPlanet_I

  character (len=lTypeBField)   :: TypeBFieldPlanet_I(0:MaxPlanet+1)
  real,dimension(0:MaxPlanet+1) :: DipoleStrengthPlanet_I
  real,dimension(0:MaxPlanet+1) :: bAxisThetaPlanet_I, bAxisPhiPlanet_I

  real,dimension(0:MaxPlanet+1) :: IonoHeightPlanet_I

  character (len=lNamePlanet)   :: NamePlanet_I(0:MaxPlanet+1)

  integer :: Planet_ 

  ! Below are defining constants for all astronomical bodies.  They are
  ! grouped using a system similar to JPL's naif/spice toolkit although
  ! the definitions are not quite the same.

  ! First define the storage location for all bodies.  This is so that you
  ! can easily find the index and can also see the naming system

  ! No Planet (in other words, no body)
  integer,parameter :: NoPlanet_  =  0

  ! New Planet (a body that is not in the database below)
  integer,parameter :: NewPlanet_  =  MaxPlanet+1

  ! Sun, planets + Pluto
  integer,parameter :: Sun_       =  1
  integer,parameter :: Mercury_   = 10
  integer,parameter :: Venus_     = 20
  integer,parameter :: Earth_     = 30
  integer,parameter :: Mars_      = 40
  integer,parameter :: Jupiter_   = 50
  integer,parameter :: Saturn_    = 60
  integer,parameter :: Uranus_    = 70
  integer,parameter :: Neptune_   = 80
  integer,parameter :: Pluto_     = 90

  ! Moons of planets (the order of the moons is not in radial distance)
  integer,parameter :: Moon_      = 31
  integer,parameter :: Io_        = 51
  integer,parameter :: Europa_    = 52
  integer,parameter :: Titan_     = 61
  integer,parameter :: Enceladus_ = 62

  ! For other solar system bodies (comets, asteroids, extra solar planets) 
  ! the index is 100 or above
  integer,parameter :: Halley_               = 100
  integer,parameter :: Comet1P_              = 100
  integer,parameter :: Borrelly_             = 101
  integer,parameter :: Comet19P_             = 101
  integer,parameter :: ChuryumovGerasimenko_ = 102
  integer,parameter :: Comet67P_             = 102
  integer,parameter :: HaleBopp_             = 103 

contains

   subroutine init_planet_const

     use ModUtilities,     ONLY: upper_case   

     implicit none

     save

     integer :: i
     !-----------------------------------------------------------------------
     ! Initialize all values - below set only the non-default values
     NamePlanet_I                     = ''

     rPlanet_I                        = 0.0                     ! [ m]
     mPlanet_I                        = 0.0                     ! [kg]
     rOrbitPlanet_I                   = 0.0                     ! [ m]
     OrbitalPeriodPlanet_I            = 0.0                     ! [ s]
     RotationPeriodPlanet_I           = 0.0                     ! [ s]
                                            
     iYearEquinoxPlanet_I             =2000                     ! [yr]
     iMonthEquinoxPlanet_I            =   1                     ! [mo]
     iDayEquinoxPlanet_I              =   1                     ! [dy]
     iHourEquinoxPlanet_I             =   0                     ! [hr]
     iMinuteEquinoxPlanet_I           =   0                     ! [mn]
     iSecondEquinoxPlanet_I           =   0                     ! [ s]
     FracSecondEquinoxPlanet_I        = 0.0                     ! [ s]
     TiltPlanet_I                     = 0.0                     ! [rad]
                         
     TypeBFieldPlanet_I               = "NONE"                
     DipoleStrengthPlanet_I           = 0.0                     ! [ T]
     bAxisThetaPlanet_I               = 0.0                     ! [rad]
     bAxisPhiPlanet_I                 = 0.0                     ! [rad]
                                          
     IonoHeightPlanet_I               = 0.0                     ! [ m]

     !\
     ! Mercury (10)                        
     !/ 
     NamePlanet_I(Mercury_)                = 'MERCURY'

     rPlanet_I(Mercury_)                   = 2439.0e+3               ! [ m]
     mPlanet_I(Mercury_)                   = 3.3022e+23               ! [kg]
     OrbitalPeriodPlanet_I(Mercury_)       = 87.969* cDay          ! [ s]
     RotationPeriodPlanet_I(Mercury_)      = 175.942* cDay          ! [ s]

     TypeBFieldPlanet_I(Mercury_)         = "DIPOLE"
     DipoleStrengthPlanet_I(Mercury_)     = -200.0e-9              ! [ T]

     !\                                
     ! Venus (20)                          
     !/                                
     NamePlanet_I(Venus_)                = 'VENUS'

     rPlanet_I(Venus_)                   = 6052.0e+3               ! [ m]
     mPlanet_I(Venus_)                   = 4.865e+24               ! [kg]
     OrbitalPeriodPlanet_I(Venus_)       = 224.7   * cDay          ! [ s]
     RotationPeriodPlanet_I(Venus_)      = 116.750 * cDay          ! [ s]
                                       
     !\                                
     ! Earth (30)                         
     !/                                
     NamePlanet_I(Earth_)                = 'EARTH'

     rPlanet_I(Earth_)                   = 6378.0e+3               ! [ m]
     mPlanet_I(Earth_)                   = 5.976e+24               ! [kg]
     rOrbitPlanet_I(Earth_)              = cAU                     ! [ m]
     OrbitalPeriodPlanet_I(Earth_)       = 365.24218967 * cDay     ! [ s]
     RotationPeriodPlanet_I(Earth_)      = cDay                    ! [ s]
                                       
     iYearEquinoxPlanet_I(Earth_)        = 2000                    ! [yr]
     iMonthEquinoxPlanet_I(Earth_)       =    3                    ! [mo]
     iDayEquinoxPlanet_I(Earth_)         =   20                    ! [dy]
     iHourEquinoxPlanet_I(Earth_)        =    7                    ! [hr]
     iMinuteEquinoxPlanet_I(Earth_)      =   35                    ! [mn]
     iSecondEquinoxPlanet_I(Earth_)      =    0                    ! [ s]
     FracSecondEquinoxPlanet_I(Earth_)   =  0.0                    ! [ s]
     TiltPlanet_I(Earth_)                = 23.5 * cDegToRad        ! [rad]
   
     TypeBFieldPlanet_I(Earth_)          = "DIPOLE"                
     DipoleStrengthPlanet_I(Earth_)      = -31100.0e-9             ! [ T]
     bAxisThetaPlanet_I(Earth_)          =  11.0 * cDegToRad       ! [rad]
     bAxisPhiPlanet_I(Earth_)            = 289.1 * cDegToRad       ! [rad]
                                       
     IonoHeightPlanet_I(Earth_)          = 110000.0                ! [ m]
   
     !\                               
     ! Moon (31)                         
     !/                               
     NamePlanet_I(Moon_)                 = 'MOON'

     rPlanet_I(Moon_)                    = 1737.0e+3               ! [ m]
     mPlanet_I(Moon_)                    = 7.3477e+22              ! [kg]
     OrbitalPeriodPlanet_I(Moon_)        = 27.321582 * cDay        ! [ s]
     RotationPeriodPlanet_I(Moon_)       = 27.321582 * cDay        ! [ s]

     !\                               
     ! Mars (40)                         
     !/                               
     NamePlanet_I(Mars_)                 = 'MARS'

     rPlanet_I(Mars_)                    = 3396.0e+3               ! [ m]
     mPlanet_I(Mars_)                    = 0.6436e+24              ! [kg]
     OrbitalPeriodPlanet_I(Mars_)        = 686.98* cDay            ! [ s]
     RotationPeriodPlanet_I(Mars_)       = 1.0275 * cDay            ! [ s]
                                        
     !\
     ! Jupiter (50)
     !/
     NamePlanet_I(Jupiter_)              = 'JUPITER'

     rPlanet_I(Jupiter_)                 = 71492.0e+3              ! [ m]
     mPlanet_I(Jupiter_)                 = 1.8980e+27              ! [kg]
     OrbitalPeriodPlanet_I(Jupiter_)     = 4330.6 * cDay           ! [ s]
     RotationPeriodPlanet_I(Jupiter_)    = 9.9259 * cHour           ! [ s]

     TypeBFieldPlanet_I(Jupiter_)        = "DIPOLE"                
     DipoleStrengthPlanet_I(Jupiter_)    = 428000.0e-9             ! [ T]
                                       
     IonoHeightPlanet_I(Jupiter_)        = 1000.0e+3               ! [ m]
   
     !\                               
     ! Saturn (60)                        
     !/                               
     NamePlanet_I(Saturn_)               = 'SATURN'

     rPlanet_I(Saturn_)                  = 60268.0e+3              ! [ m]
     mPlanet_I(Saturn_)                  = 0.5685e+27              ! [kg]
     OrbitalPeriodPlanet_I(Saturn_)      = 10746.94 * cDay         ! [ s]
     RotationPeriodPlanet_I(Saturn_)     = 10.656 * cHour            ! [ s]
                                       
     TypeBFieldPlanet_I(Saturn_)         = "DIPOLE"                
     DipoleStrengthPlanet_I(Saturn_)     = 20800.0e-9              ! [ T]
                                       
     IonoHeightPlanet_I(Saturn_)         = 1000.0e+3               ! [ m]
   
     !\
     ! Uranus (70)
     !/
   
     !\
     ! Neptune (80)
     !/
   
     !\
     ! Pluto (90)
     !/
   
     !\
     ! Io (51)
     !/
     NamePlanet_I(Io_)                   = 'IO'

     rPlanet_I(Io_)                      = 1821.0e+3               ! [ m]
     mPlanet_I(Io_)                      = 0.0                     ! [kg]
     OrbitalPeriodPlanet_I(Io_)          = 0.0                     ! [ s]
     RotationPeriodPlanet_I(Io_)         = 0.0                     ! [ s]

     !\
     ! Europa (52)
     !/
     NamePlanet_I(Europa_)               = 'EUROPA'

     rPlanet_I(Europa_)                  = 1569.0e+3               ! [ m]
     mPlanet_I(Europa_)                  = 4.80e22                 ! [kg]
     OrbitalPeriodPlanet_I(Europa_)      = 3.551 * cDay            ! [ s]
     RotationPeriodPlanet_I(Europa_)     = 3.551 * cDay            ! [ s]

     iYearEquinoxPlanet_I(Europa_)       = 2000                    ! [yr]
     iMonthEquinoxPlanet_I(Europa_)      =    1                    ! [mo]
     iDayEquinoxPlanet_I(Europa_)        =    1                    ! [dy]
     iHourEquinoxPlanet_I(Europa_)       =    0                    ! [hr]
     iMinuteEquinoxPlanet_I(Europa_)     =    0                    ! [mn]
     iSecondEquinoxPlanet_I(Europa_)     =    0                    ! [ s]
     FracSecondEquinoxPlanet_I(Europa_)  =  0.0                    ! [ s]
     TiltPlanet_I(Europa_)               =  0.0 * cDegToRad        ! [rad]

     TypeBFieldPlanet_I(Europa_)         = "DIPOLE"                
     DipoleStrengthPlanet_I(Europa_)     =    100.0e-9             ! [ T]
     bAxisThetaPlanet_I(Europa_)         =  90.0 * cDegToRad       ! [rad]
     bAxisPhiPlanet_I(Europa_)           =   0.0 * cDegToRad       ! [rad]
     
     !\
     ! Titan (61)
     !/
     NamePlanet_I(Titan_)                = 'TITAN'

     rPlanet_I(Titan_)                   = 2575.0e+3               ! [ m]
     mPlanet_I(Titan_)                   = 0.1346e+24              ! [kg]
     rOrbitPlanet_I(Titan_)              = 1.222e+9                ! [ m]
     OrbitalPeriodPlanet_I(Titan_)       = 15.945 * cDay           ! [ s]
     RotationPeriodPlanet_I(Titan_)      = 15.945 * cDay           ! [ s]

     !\
     ! Enceladus (62)
     !/
     NamePlanet_I(Enceladus_)                = 'ENCELADUS'

     rPlanet_I(Enceladus_)                   = 252.0e+3            ! [ m]
     mPlanet_I(Enceladus_)                   = 8.4e+19             ! [kg]
     rOrbitPlanet_I(Enceladus_)              = 2.3802e+8           ! [ m]
     OrbitalPeriodPlanet_I(Enceladus_)       = 1.370218 * cDay     ! [ s]
     RotationPeriodPlanet_I(Enceladus_)      = 1.370218 * cDay     ! [ s]

     ! This is the field of Saturn but using a different distance unit !!!
     TypeBFieldPlanet_I(Enceladus_)         = "DIPOLE"
     DipoleStrengthPlanet_I(Enceladus_) = &
          DipoleStrengthPlanet_I(Saturn_)* &
          (rPlanet_I(Saturn_)/rPlanet_I(Enceladus_))**3

     !\
     ! Comets
     !/
     NamePlanet_I(Halley_) = 'HALLEY'

     rPlanet_I(Halley_)                   = 1.0E9    ! [ m] for distance units
     mPlanet_I(Halley_)                   = 1.0      ! [kg]
     rOrbitPlanet_I(Halley_)              = cAU      ! [ m]
     OrbitalPeriodPlanet_I(Halley_)       = 1.0      ! [ s]
     RotationPeriodPlanet_I(Halley_)      = 1.0      ! [ s]

     !\
     ! No Planet (0)
     !     - No Planet and no body - defaults for everything, just set name.
     !/
     NamePlanet_I(NoPlanet_)             = 'NONE'

     !\
     ! New Planet (MaxPlanet+1)
     !     - A planet whose parameters are not defined in the database above.
     !       A few values are set to clearly meaningless values to prevent
     !       a users from using the values incorrectly. 
     !/
     NamePlanet_I(NewPlanet_)            = 'NEW/UNKNOWN'
     rPlanet_I(NewPlanet_)               = -1.0
     mPlanet_I(NewPlanet_)               = -1.0

     ! make all the planet names upper case
     do i=NoPlanet_,MaxPlanet+1
        call upper_case(NamePlanet_I(i))  ! make all the names upper case
     end do

   end subroutine init_planet_const

end module ModPlanetConst
   

!====================================================================
! Documentation
!====================================================================


! The orbital period belongs to the TROPICAL YEAR, which is
! relative to the vernal equinox which is slowly moving 
! due to the precession of the Earth's rotation axis.

! The rotational angular velocity is relative to an inertial frame

! Reference equinox time taken from 
! http://aa.usno.navy.mil/data/docs/Earth_Seasons.html

! The angle between the zero meridian and the eqinox direction at 
! equinox time. For Earth this can be calculated from the time of day.
! For other planets there is no analogous method to calculate this angle.
