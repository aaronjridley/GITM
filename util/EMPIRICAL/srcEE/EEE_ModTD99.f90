!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
module EEE_ModTD99
  use EEE_ModCommonVariables
  implicit none
  save

  private

  public :: set_parameters_TD99,get_transformed_TD99fluxrope

  !----------------------------------------------------------------------
  ! Logical variables related to the magnetic field computation::
  logical, public :: DoBqField=.false.

  logical :: UseVariedCurrent=.false.

  logical :: DoEquilItube=.false.
  logical :: DoRevCurrent=.false.
  logical :: DoMaintainEpot=.false.
  real :: CurrentRiseTime,CurrentStartTime

  !----------------------------------------------------------------------
  ! Variables related to the position of the flux rope::
  real, parameter :: Li_TD99=0.5

  !----------------------------------------------------------------------
  ! Variables related to the flux rope properties::
  real :: Itube_TD99=0.0, Rtube_TD99=0.0, atube_TD99=0.0
  real :: d_TD99=0.0, aratio_TD99=0.0
  real :: Mass_TD99=0.0
  real :: InvH0_TD99=0.0, Rho0_TD99=0.0
  real :: ItubeSaved=0.0

  !----------------------------------------------------------------------
  ! Variables related to the properties of the strapping field, Bq::
  integer :: nStepSaved=-1
  real, parameter :: AlphaRamp=9.52381E-04         !in [-]
  real, parameter :: VTransX=1.500E+03             !in [m/s]
  real, parameter :: VTransY=-2.900E+04            !in [m/s]
  real, parameter :: UVorCMax0=2.5                 !in units of 100 km/s
  real, parameter :: BqZMax0=3.768210E+01          !in [Gauss]
  real :: BqZMax=0.0, BqZMaxSaved=0.0, UVorCMax=0.0
  real :: q_TD99=0.0, L_TD99=0.0

  !----------------------------------------------------------------------
  ! Declare the rotational matrix of coordinate transformation::
  real, dimension(3,3) :: RotateTD99_DD

contains

  !============================================================================

  subroutine set_parameters_TD99(NameCommand)
    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'set_parameters_TD99'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#TD99FLUXROPE")
       call read_var('UseCme',              UseCme)
       call read_var('UseVariedCurrent'   , UseVariedCurrent)
       call read_var('CurrentStartTime'   , CurrentStartTime)
       call read_var('CurrentRiseTime '   , CurrentRiseTime)
       call read_var('DoAddFluxRope'      , DoAddFluxRope)
       call read_var('DoEquilItube'       , DoEquilItube)
       call read_var('DoRevCurrent'       , DoRevCurrent)
       call read_var('aratio_TD99'        , aratio_TD99)
       call read_var('Itube_TD99'         , Itube_TD99)
       call read_var('Rtube_TD99'         , Rtube_TD99)
       call read_var('atube_TD99'         , atube_TD99)
       call read_var('d_TD99'             , d_TD99)
       call read_var('Mass_TD99'          , Mass_TD99)
       call read_var('LongitudeCme'       , LongitudeCme)
       call read_var('LatitudeCme'        , LatitudeCme)
       call read_var('OrientationCme'     , OrientationCme)
       call read_var('DoBqField'          , DoBqField)
       call read_var('q_TD99'             , q_TD99)
       call read_var('L_TD99'             , L_TD99)
    case("#CME")
       call read_var('Current',     Itube_TD99)
       call read_var('RadiusMajor', Rtube_TD99)
       call read_var('RadiusMinor', atube_TD99)
       call read_var('Depth',       d_TD99)
       call read_var('Mass',        Mass_TD99)
    case default
       call CON_stop(NameSub//' unknown NameCommand='//NameCommand)
    end select

  end subroutine set_parameters_TD99

  !============================================================================

  real function varied_current(Time)
    real, intent(in):: Time
    !--------------------------------------------------------------------------
    varied_current = min(&
         max(Time-CurrentStartTime,0.0)/CurrentRiseTime,1.0)*&
         ItubeSaved 
  end function varied_current

  !============================================================================

  real function time_increment(Time)
    real, intent(in):: Time
    !--------------------------------------------------------------------------
    time_increment = min(&
         max(Time-CurrentStartTime,0.0)/CurrentRiseTime,1.0)
  end function time_increment

  !============================================================================

  subroutine get_transformed_TD99fluxrope(RFace_D,BFRope_D,UVorT_D,&
       n_step,Iteration_Number,RhoFRope,Time)

    real, dimension(3), intent(in) :: RFace_D
    real, dimension(3), intent(out) :: BFRope_D,UVorT_D
    integer, intent(in) :: n_step,Iteration_Number
    real, intent(out), optional :: RhoFRope
    real, intent(in), optional :: Time

    logical:: DoFirstCallTD99=.true.
    real:: atemp,Itemp
    real:: UVorR
    real, dimension(3):: B1FRopeTemp_D
    real, dimension(3):: R1Face_D,B1FRope_D,B1qField_D
    real, dimension(3):: UVorC_D,U1VorC_D
    !--------------------------------------------------------------------------
    ! Initialize the TD99 model parameters once::

    if (DoFirstCallTD99) then
       call init_TD99_parameters
       DoFirstCallTD99=.false.
    endif
    if (present(Time).and.UseVariedCurrent) &
         Itube_TD99 = varied_current(Time)

    ! Check if the potential electric field needs to be applied.
    ! Maintain Epot for (CurrentStartTime < Time < CurrentRiseTime)::

    if (present(Time).and.DoBqField) then
       DoMaintainEpot = (&
            (Time.gt.CurrentStartTime).and.&
            (Time.lt.CurrentRiseTime))
    else
       DoMaintainEpot = .false.
    endif

    ! Compute the flux rope, and transform coordinates and vectors to
    ! position the flux rope in the desired way::

    R1Face_D = matmul(RotateTD99_DD,RFace_D)
    if (DoAddFluxRope) then
       call compute_TD99_FluxRope(R1Face_D,B1FRope_D,RhoFRope)
       if (DoRevCurrent) then
          Itemp = Itube_TD99; Itube_TD99 = -Itemp
          atemp = atube_TD99; atube_TD99 = aratio_TD99*atemp
          call compute_TD99_FluxRope(R1Face_D,B1FRopeTemp_D)
          B1FRope_D = B1FRope_D+B1FRopeTemp_D
          Itube_TD99 = Itemp
          atube_TD99 = atemp
       endif
    else
       B1FRope_D = 0.0
       RhoFRope  = 0.0
    endif
    U1VorC_D = 0.0
    if (DoBqField) then
       if (present(Time)) then
          call compute_TD99_BqField(R1Face_D,B1qField_D,&
               U1VorC_D,DoMaintainEpot,n_step,Iteration_Number,Time)
       else
          call compute_TD99_BqField(R1Face_D,B1qField_D,&
               U1VorC_D,DoMaintainEpot,n_step,Iteration_Number)
       endif
       B1FRope_D = B1FRope_D+B1qField_D
    endif
    BFRope_D = matmul(B1FRope_D,RotateTD99_DD)
    UVorC_D  = matmul(U1VorC_D,RotateTD99_DD)

    ! Compute the tangential component of the velocity field, UVorT_D,
    ! associated with the potential electric field, Epot::

    UVorR   = dot_product(RFace_D,UVorC_D)
    UVorT_D = UVorC_D-RFace_D*UVorR

    BFRope_D = BFRope_D*No2Si_V(UnitB_)
    UVorT_D = UVorT_D*No2Si_V(UnitU_)
    RhoFRope = RhoFRope*No2Si_V(UnitRho_)

  end subroutine get_transformed_TD99fluxrope

  !============================================================================

  subroutine init_TD99_parameters

    use ModCoordTransform, ONLY: rot_matrix_x,rot_matrix_y,rot_matrix_z

    real:: AlphaRope,LInduct,WFRope,FootSepar,ItubeDim
    !--------------------------------------------------------------------------

    ! Compute the magnetic energy, WFRope, associated with the portion
    ! of the flux rope current that is above the solar surface::

    InvH0_TD99 = cGravitation*Msun/Rsun*Si2No_V(UnitU_)**2   ! in [-]
    AlphaRope  = 2.0*acos(d_TD99/Rtube_TD99)                 ! in [rad]
    FootSepar  = Rtube_TD99*sin(0.5*AlphaRope)/1.0e6         ! in [Mm]
    LInduct    = cMu*(0.5*AlphaRope/cPi)*Rtube_TD99*log(2.0**3 &
         *(Rtube_TD99-d_TD99)/atube_TD99-2.0+0.25)           ! in [H]
    WFRope     = 0.5*LInduct*Itube_TD99**2*1.0e7             ! in [ergs]

    ! Compute the average density inside the flux rope assuming that the
    ! total amount of prominence mass is Mass_TD99 (=10^16g=10^13kg)::

    Rho0_TD99  = Mass_TD99/(AlphaRope*Rtube_TD99*cPi*atube_TD99**2)
    ! in [kg/m^3]

    ! Define the normalized model parameters here::

    ! Flux rope::
    Rtube_TD99 = Rtube_TD99*Si2No_V(UnitX_)
    atube_TD99 = atube_TD99*Si2No_V(UnitX_)
    ItubeDim   = Itube_TD99
    Itube_TD99 = ItubeDim*Si2No_V(UnitJ_)*Si2No_V(UnitX_)**2 ! in [A]
    Rho0_TD99  = Rho0_TD99*Si2No_V(UnitRho_)

    ! Save the maximum value of the current for possible use in
    ! varied_current case::

    ItubeSaved = Itube_TD99 

    ! Strapping field::

    d_TD99     = d_TD99*Si2No_V(UnitX_)
    L_TD99     = L_TD99*Si2No_V(UnitX_)
    q_TD99     = q_TD99*Si2No_V(UnitB_)*Si2No_V(UnitX_)**2

    ! Construct the rotational matrix, RotateTD99_DD, to position the
    ! flux rope in the desired way on the solar surface::

    RotateTD99_DD = matmul(rot_matrix_y(-0.5*cPi),&
         rot_matrix_x(-OrientationCme*cDegToRad))
    RotateTD99_DD = matmul(RotateTD99_DD,          &
         rot_matrix_y(LatitudeCme*cDegToRad))
    RotateTD99_DD = matmul(RotateTD99_DD,          &
         rot_matrix_z(-LongitudeCme*cDegToRad))

    if (iProc==0) then
       write(*,*) prefix
       write(*,*) prefix,'>>>>>>>>>>>>>>>>>>>                   <<<<<<<<<<<<<<<<<<<<<'
       write(*,*) prefix
       write(*,*) prefix,'    Twisted Flux Rope Model by Titov & Demoulin, 1999.     '
       write(*,*) prefix
       write(*,*) prefix,'>>>>>>>>>>>>>>>>>>>                   <<<<<<<<<<<<<<<<<<<<<'
       write(*,*) prefix
       write(*,*) prefix,'d_TD99      = ',d_TD99*No2Si_V(UnitX_)/1.0E6,'[Mm]'
       write(*,*) prefix,'Rtube_TD99  = ', &
            Rtube_TD99*No2Si_V(UnitX_)/1.0E6,'[Mm]'
       write(*,*) prefix,'atube_TD99  = ', &
            atube_TD99*No2Si_V(UnitX_)/1.0E6,'[Mm]'
       write(*,*) prefix,'atube/Rtube = ',atube_TD99/Rtube_TD99,'[-]'
       write(*,*) prefix,'Itube_TD99  = ',ItubeDim,'[A]'
       write(*,*) prefix,'aratio_TD99 = ',aratio_TD99,'[-]'
       write(*,*) prefix,'Mass_TD99   = ',Mass_TD99*1.0e3,'[g] '
       write(*,*) prefix,'Rho0_TD99   = ',Rho0_TD99*No2Io_V(UnitRho_),'[g/cm^3]'
       write(*,*) prefix
       write(*,*) prefix,'q_TD99      = ', &
            q_TD99*No2Si_V(UnitB_)*No2Si_V(UnitX_)**2,'[T m^2]'
       write(*,*) prefix,'L_TD99      = ',L_TD99*No2Si_V(UnitX_)/1.0e6,'[Mm]'
       write(*,*) prefix
       write(*,*) prefix,'Free energy of flux rope is ',WFRope,'Ergs.'
       write(*,*) prefix,'Separation of flux rope ends is ',FootSepar,'Mm,'
       write(*,*) prefix,'   or ',cPi*FootSepar*1.0e6/(2.0*Rsun)*cRadToDeg,'deg.'
       write(*,*) prefix
       if (UseVariedCurrent) then
          write(*,*) prefix,'>>>>>       UseVariedCurrent is set to .true.!!!      <<<<<'
          write(*,*) prefix,'CurrentStartTime = ',CurrentStartTime,'[s]'
          write(*,*) prefix,'CurrentRiseTime  = ',CurrentRiseTime,'[s]'
          write(*,*) prefix
       endif
    endif
    if (DoEquilItube) then

       ! Compute the equilibrium toroidal current, Itube_TD99, based
       ! on the force balance in direction normal to the surface of
       ! the flux tube.

       Itube_TD99 = 8.0*cPi*q_TD99*L_TD99*Rtube_TD99 &
            *(L_TD99**2+Rtube_TD99**2)**(-1.5) &
            /(alog(8.0*Rtube_TD99/atube_TD99) &
            -1.5+Li_TD99/2.0)                           ! in [-]
       WFRope    = 0.5*LInduct*(ItubeDim)**2*1.0e7      ! in [ergs]
    endif
    if (DoEquilItube.and.iProc==0) then
       write(*,*) prefix,'The strapping field, Bq, is added and the EQUILIBRIUM value'
       write(*,*) prefix,'of Itube_TD99 is computed!!!'
       write(*,*) prefix
       write(*,*) prefix,'The value of Itube_TD99 is reset to :: ',Itube_TD99
       write(*,*) prefix,'The free energy of the flux rope is :: ',WFRope,'Ergs.'
       write(*,*) prefix
    endif

  end subroutine init_TD99_parameters

  !=====================================================================!

  subroutine compute_TD99_BqField(RFace_D,BqField_D,UVorC_D,&
       DoMaintainEpot,n_step,Iteration_Number,TimeNow)

    real, intent(in), dimension(3) :: RFace_D
    real, intent(out), dimension(3) :: BqField_D,UVorC_D
    logical, intent(in) :: DoMaintainEpot
    integer, intent(in) :: n_step,Iteration_Number
    real, intent(in), optional :: TimeNow

    ! Variables related to coordinates::
    real:: R2Plus,R2Mins
    real, dimension(3):: RPlus_D,RMins_D

    ! Variables related to computations of potential electric field::
    real:: BqZOverBqZ0,BqZFunction
    real, dimension(3):: EpotC_D,UTranC_D
    real, dimension(3):: GradBqZ_D,GradPsiC_D
    !--------------------------------------------------------------------

    ! Compute the locations, RMins_D and RPlus_D, of the two magnetic
    ! charges, -/+q::

    if (present(TimeNow)) then
       RPlus_D(x_) = RFace_D(x_)-L_TD99 &
            - VTransX*(TimeNow-CurrentStartTime)*Si2No_V(UnitX_)
       RMins_D(x_) = RFace_D(x_)+L_TD99 &
            + VTransX*(TimeNow-CurrentStartTime)*Si2No_V(UnitX_)
       RPlus_D(y_) = RFace_D(y_) &
            - VTransY*(TimeNow-CurrentStartTime)*Si2No_V(UnitX_)
       RMins_D(y_) = RFace_D(y_) &
            + VTransY*(TimeNow-CurrentStartTime)*Si2No_V(UnitX_)
    else
       RPlus_D(x_) = RFace_D(x_)-L_TD99
       RMins_D(x_) = RFace_D(x_)+L_TD99
       RPlus_D(y_) = RFace_D(y_)
       RMins_D(y_) = RPlus_D(y_)
    endif
    RPlus_D(z_) = RFace_D(z_)+d_TD99-1.0
    RMins_D(z_) = RPlus_D(z_)
    R2Plus = sqrt(dot_product(RPlus_D,RPlus_D))
    R2Mins = sqrt(dot_product(RMins_D,RMins_D))

    ! Compute the field of the strapping magnetic field, BqField_D::

    BqField_D = q_TD99*(RPlus_D/R2Plus**3 - RMins_D/R2Mins**3)

    ! Update the values of BqZMax and BqZMaxSaved once at each iteration::

    if (n_step/=nStepSaved) then
       nStepSaved  = n_step
       BqZMaxSaved = BqZMax
       BqZMax  = 0.0
    end if
    if (iteration_number==1) &
         BqZMaxSaved = BqZMax0/No2Io_V(UnitB_)
    BqZMax = max(abs(BqField_D(z_)),BqZMax)

    ! Apply Epot only if DoMaintainEpot=.true.!!!

    if (DoMaintainEpot) then

       ! Compute the gradient of the z-component of the Bq field on the
       ! solar surface in Cartesian geometry -- GradBqZ_D(x_:z_)::

       GradBqZ_D = 3.0*q_TD99*RMins_D(z_)*&
            (RMins_D/R2Mins**5 - RPlus_D/R2Plus**5)
       GradBqZ_D(z_)    = 0.0
       !    GradBqZ_D(x_:z_) = GradBqZ_D(x_:z_)  - &
       !         RFace_D(x_:z_)*dot_product(RFace_D,GradBqZ_D)

       ! Compute the gradient of the scalar potential in Cartesian
       ! geometry -- GradPsiC_D::

       BqZOverBqZ0 = min(1.0,abs(BqField_D(z_))/BqZMaxSaved)
       BqZFunction = max(0.0,1.0-BqZOverBqZ0**2)
       if (BqZOverBqZ0.gt.0.0) then
          GradPsiC_D = GradBqZ_D*2.0 &
               *AlphaRamp*BqZFunction*BqZOverBqZ0*exp(-BqZFunction)
       else
          GradPsiC_D = 0.0
       endif

       ! Compute the potential electric field on the solar surface to
       ! be applied at one of the spots -- EpotC_D(x_:z_).
       ! This is given by:
       ! EpotC_D(x_:z_) = sign(BqField_D(z_))*BqZOverBqZ0*GradPsiC_D::

       EpotC_D = sign(1.0,BqField_D(z_))*GradPsiC_D*BqZOverBqZ0

       ! Compute the plasma velocity, UVorC_D, associated with the static
       ! eletric field, EpotC_D.
       ! This is given by:
       ! UVorC_D(x_:z_) = sign(BqField_D(z_))*GradPsiC_D(x_:z_) X Ez,
       ! where Ez = (0,0,1)

       UVorC_D(x_)  =  GradPsiC_D(y_)*sign(1.0,BqField_D(z_))
       UVorC_D(y_)  = -GradPsiC_D(x_)*sign(1.0,BqField_D(z_))
       UVorC_D(z_)  =  0.0
       UVorC_D      = UVorC_D*UVorCMax0
       if (iteration_number.gt.2) &
            UVorC_D = UVorC_D*UVorCMax0/UVorCMax

       ! Compute the translational velocity, UTranC_D, at which the two
       ! q-sources are moved with respect to each other.

       if (present(TimeNow)) then
          UTranC_D(x_) = (VTransX*Si2No_V(UnitU_))*BqZOverBqZ0 &
               *sign(1.0,BqField_D(z_))*exp(-BqZFunction)
          UTranC_D(y_) = (VTransY*Si2No_V(UnitU_))*BqZOverBqZ0 &
               *sign(1.0,BqField_D(z_))*exp(-BqZFunction)
          UTranC_D(z_) = 0.0
       else
          UTranC_D = 0.0
       endif
       ! Add the translational velocity, UTranC_D, to UVorC_D::
       UVorC_D = UVorC_D+UTranC_D
    else
       EpotC_D = 0.0
       UVorC_D = 0.0
    endif

  end subroutine compute_TD99_BqField

  !=====================================================================!

  subroutine compute_TD99_FluxRope(RFace_D,BFRope_D,RhoFRope)

    !\__                                                             __/!
    !    Twisted Magnetic Field Configuration by Titov & Demoulin '99   !
    !                                                                   !
    ! An instability that causes a CME eruption is expected to occur at !
    ! R > L*sqrt(2). For a detailed description of the initial state    !
    ! refer to A&A, 1999, v.351, pp.707-720                             !
    !                                                                   !
    ! ___  This module was written by Ilia Roussev on June 10, 2002 ___ !
    !/                                                                 \!

    real, intent(in), dimension(3):: RFace_D
    real, intent(out), dimension(3):: BFRope_D

    real, intent(out), optional:: RhoFRope

    real:: xxx,yyy,zzz,R2Face
    real:: RhoTB_TD99,Rperp_TD99
    real:: CHIin_TD99,CHIex_TD99
    real:: xUVx_TD99,xUVy_TD99,xUVz_TD99
    real:: ThetaUVy_TD99,ThetaUVz_TD99
    real:: RperpUVx_TD99,RperpUVy_TD99,RperpUVz_TD99
    real:: Kappa_TD99,dKappadx_TD99,dKappadr_TD99
    real:: KappaA_TD99,dKappaAdx_TD99,dKappaAdr_TD99

    ! Complete elliptic integrals related variables::
    real:: K_elliptic, E_elliptic
    real:: K_ellipticA, E_ellipticA

    ! Vector potential related variables::
    real:: Ak_TD99,dAkdk_TD99
    real:: AkA_TD99,dAkdkA_TD99,d2Akdk2A_TD99
    real:: AI_TD99,dAIdx_TD99,dAIdr_TD99
    real:: AIin_TD99,dAIindx_TD99,dAIindr_TD99
    real:: AIex_TD99,dAIexdx_TD99,dAIexdr_TD99
    ! Flux-rope related variables::
    real:: BIphix_TD99,BIphiy_TD99,BIphiz_TD99
    !--------------------------------------------------------------------
    ! Assign X,Y,Z coordinates at which to compute the magnetic field::

    xxx = RFace_D(x_)
    yyy = RFace_D(y_)
    zzz = RFace_D(z_)
    R2Face = sqrt(dot_product(RFace_D,RFace_D))

    ! Compute Rperp_TD99 and RhoTB_TD99::

    Rperp_TD99 = sqrt(yyy**2+(zzz+d_TD99-1.0)**2)
    RhoTB_TD99 = sqrt(xxx**2+(Rperp_TD99-Rtube_TD99)**2)

    ! Define the Heaviside step function in the internal region 
    ! (RhoTB_TD99<atube_TD99), CHIin_TD99, and the external one
    ! (RhoTB_TD99>atube_TD99), CHIex_TD99::

    if (RhoTB_TD99.lt.atube_TD99) then
       CHIin_TD99 = 1.0
       CHIex_TD99 = 0.0
    else
       CHIin_TD99 = 0.0
       CHIex_TD99 = 1.0
    endif

    ! Add the prominence material inside the flux rope, assuming that the
    ! total amount mass is 10^13kg, and that the desnity scale-height is
    ! the same as the pressure scale-height, 1/InvH0 (i.e., iso-thermal
    ! atmoshpere)::

    if (present(RhoFRope)) &
         RhoFRope = Rho0_TD99*exp(-10.0*(RhoTB_TD99/atube_TD99)**6) &
         *exp(-InvH0_TD99*abs(R2Face-1.0))    

    ! Compute the field produced by the ring current, Itube_TD99, both
    ! inside and outside the torus, BI_TD99 = BFRope_D(x_:z_)::

    ThetaUVy_TD99 = -(zzz+d_TD99-1.0)/Rperp_TD99
    ThetaUVz_TD99 = yyy/Rperp_TD99

    ! Compute the toroidal field (BIphix_TD99, BIphiy_TD99, BIphiz_TD99)
    ! produced by the azimuthal current Iphi. This is needed to ensure
    ! that the flux rope configuration is force free. 

    BIphix_TD99 = 0.0
    BIphiy_TD99 = abs(Itube_TD99)/(2.0*cPi*atube_TD99**2) &
         *sqrt(CHIin_TD99*2.0*(atube_TD99**2-RhoTB_TD99**2)) &
         *ThetaUVy_TD99
    BIphiz_TD99 = abs(Itube_TD99)/(2.0*cPi*atube_TD99**2) &
         *sqrt(CHIin_TD99*2.0*(atube_TD99**2-RhoTB_TD99**2)) &
         *ThetaUVz_TD99

    ! Compute the components of the unit vector in the plane of symmetry
    ! x=0::

    RperpUVx_TD99 = 0.0
    RperpUVy_TD99 = yyy/Rperp_TD99
    RperpUVz_TD99 = (zzz+d_TD99-1.0)/Rperp_TD99

    ! Compute the components of the unit vector pointing in the positive 
    ! x-direction::

    xUVx_TD99 = 1.0
    xUVy_TD99 = 0.0
    xUVz_TD99 = 0.0

    ! Define two model parameters, Kappa_TD99 and KappaA_TD99::

    Kappa_TD99 = 2.0*sqrt(Rperp_TD99*Rtube_TD99 &
         /((Rperp_TD99+Rtube_TD99)**2+xxx**2))
    KappaA_TD99 = 2.0*sqrt(Rperp_TD99*Rtube_TD99 &
         /(4.0*Rperp_TD99*Rtube_TD99+atube_TD99**2))

    ! Truncate the value of Kappa_TD99::

    if (abs(1.0-Kappa_TD99).lt.cTiny/10.0) &
         Kappa_TD99 = 1.0-cTiny/10.0

    ! Compute the vector potential in the internal, AIin_TD99, and
    ! external (outside the current torus), AIex_TD99, regions::   

    call calc_elliptic_int_1kind(Kappa_TD99,K_elliptic)
    call calc_elliptic_int_2kind(Kappa_TD99,E_elliptic)
    Ak_TD99       = ((2.0-Kappa_TD99**2)*K_elliptic &
         - 2.0*E_elliptic)/Kappa_TD99
    dAkdk_TD99    = (2.0-Kappa_TD99**2)*E_elliptic &
         /(Kappa_TD99**2*(1.0-Kappa_TD99**2)) &
         - 2.0*K_elliptic/Kappa_TD99**2
    call calc_elliptic_int_1kind(KappaA_TD99,K_ellipticA)
    call calc_elliptic_int_2kind(KappaA_TD99,E_ellipticA)
    AkA_TD99      = ((2.0-KappaA_TD99**2)*K_ellipticA &
         - 2.0*E_ellipticA)/KappaA_TD99
    dAkdkA_TD99   = (2.0-KappaA_TD99**2)*E_ellipticA &
         /(KappaA_TD99**2*(1.0-KappaA_TD99**2)) &
         - 2.0*K_ellipticA/KappaA_TD99**2
    d2Akdk2A_TD99 = ((7.0*KappaA_TD99**2-4.0 &
         - KappaA_TD99**4)*E_ellipticA/(1.0-KappaA_TD99**2) &
         + (4.0-5.0*KappaA_TD99**2)*K_ellipticA) &
         /(KappaA_TD99**3*(1.0-KappaA_TD99**2))

    ! Define AIin_TD99 and AIex_TD99::

    AIex_TD99     = Itube_TD99/(2.0*cPi)*sqrt(Rtube_TD99 &
         /Rperp_TD99)*Ak_TD99
    AIin_TD99     = Itube_TD99/(2.0*cPi)*sqrt(Rtube_TD99 &
         /Rperp_TD99)*(AkA_TD99+dAkdkA_TD99*(Kappa_TD99 &
         - KappaA_TD99))

    ! Compute the vector potential, AI_TD99, of the magnetic field 
    ! produced by the ring current Itube_TD99 in the whole space::

    AI_TD99 = CHIin_TD99*AIin_TD99+CHIex_TD99*AIex_TD99

    ! Derive the BI_TD99 field from the corresponding vector potential,
    ! AI_TD99 (this involves the comp. of some nasty derivatives)::

    dKappadx_TD99  = -xxx*Kappa_TD99/(xxx**2+(Rperp_TD99+Rtube_TD99)**2)
    dKappadr_TD99  = Kappa_TD99*(Rtube_TD99**2-Rperp_TD99**2+xxx**2) &
         /(2.0*Rperp_TD99*((Rtube_TD99+Rperp_TD99)**2+xxx**2))
    dKappaAdx_TD99 = 0.0 
    dKappaAdr_TD99 = KappaA_TD99*atube_TD99**2/(2.0*Rperp_TD99 &
         *(4.0*Rperp_TD99*Rtube_TD99+atube_TD99**2))

    ! Derivative of AIin_TD99 with respect to `x` and `rperp`:: 

    dAIindx_TD99   = Itube_TD99/(2.0*cPi)*sqrt(Rtube_TD99/Rperp_TD99) &
         *(dAkdkA_TD99*dKappadx_TD99)
    dAIindr_TD99   = Itube_TD99/(2.0*cPi)*sqrt(Rtube_TD99/Rperp_TD99) &
         *(dAkdkA_TD99*dKappadr_TD99+d2Akdk2A_TD99*dKappaAdr_TD99 &
         *(Kappa_TD99-KappaA_TD99))-AIin_TD99/(2.0*Rperp_TD99)

    ! Derivative of AIex_TD99 with respect to `x` and `rperp`::

    dAIexdx_TD99   = Itube_TD99/(2.0*cPi)*sqrt(Rtube_TD99/Rperp_TD99) &
         *(dAkdk_TD99*dKappadx_TD99)
    dAIexdr_TD99   = Itube_TD99/(2.0*cPi)*sqrt(Rtube_TD99/Rperp_TD99) &
         *(dAkdk_TD99*dKappadr_TD99)-AIex_TD99/(2.0*Rperp_TD99)

    ! Derivatives of AI with respect to `x` and `rperp`::

    dAIdx_TD99 = CHIin_TD99*dAIindx_TD99+CHIex_TD99*dAIexdx_TD99
    dAIdr_TD99 = CHIin_TD99*dAIindr_TD99+CHIex_TD99*dAIexdr_TD99

    ! Obtain the BI_TD99 field in the whole space from the corresponding
    ! vector potential, AI_TD99 -->
    ! BI_TD99 = curl(AI_TD99*ThetaUV_TD99) = BFRope_D(x_:z_)::

    BFRope_D(x_) = -dAIdx_TD99*RperpUVx_TD99 &
         + (dAIdr_TD99+AI_TD99/Rperp_TD99)*xUVx_TD99
    BFRope_D(y_) = -dAIdx_TD99*RperpUVy_TD99 &
         + (dAIdr_TD99+AI_TD99/Rperp_TD99)*xUVy_TD99
    BFRope_D(z_) = -dAIdx_TD99*RperpUVz_TD99 &
         + (dAIdr_TD99+AI_TD99/Rperp_TD99)*xUVz_TD99

    ! Add the field of the azimuthal current, Iphi::

    BFRope_D(x_) = BFRope_D(x_)+BIphix_TD99
    BFRope_D(y_) = BFRope_D(y_)+BIphiy_TD99
    BFRope_D(z_) = BFRope_D(z_)+BIphiz_TD99

  contains
    !====================================================================
    subroutine calc_elliptic_int_1kind(Kappa,K_elliptic)

      real, intent(in):: Kappa
      real, intent(out):: K_elliptic
      !------------------------------------------------------------------
      integer:: iN
      real,parameter:: pK_LIMIT1 = 0.5*cSqrtTwo
      real,parameter:: pK_LIMIT2 = 0.993
      real:: DESIRED_CEI_ACCURACY
      real:: pK,pK1,K_ell_sum,K_ell_sum_old
      real:: TwoN_1FactOverNFact

      ! Compute the complete elliptic integral of 1st kind from the series
      ! representations given by ...
      ! see formulae 8.113.1 (for 0<k<0.701) and 8.113.2 (for 0.701=<k<1)
      ! therein::
      !
      ! The stability is ensured up to pK = 0.9935 (sin**-1=83.5deg)!!!
      !
      ! Set the desired accuracy for the integral computation::KappaA_TD99

      if (iRealPrec==1) then
         DESIRED_CEI_ACCURACY = 1.0E-15
      else
         DESIRED_CEI_ACCURACY = 1.0E-7
      endif

      ! Initialize some variables::

      iN                  = 1
      TwoN_1FactOverNFact = 1.0 
      pK                  = Kappa
      pK1                 = sqrt(1.0-pK**2)

      ! Compute the CEI of 1st kind::

      if (abs(pK).lt.pK_LIMIT1) then
         K_ell_sum_old = 0.0
         K_ell_sum     = (1.0+pK**2/4.0)*cPi/2.0
         do while (abs(K_ell_sum-K_ell_sum_old).gt.DESIRED_CEI_ACCURACY)
            iN                  = iN+1
            TwoN_1FactOverNFact = TwoN_1FactOverNFact*(2.0*iN-1.0)/iN
            K_ell_sum_old       = K_ell_sum
            K_ell_sum           = K_ell_sum+(TwoN_1FactOverNFact/2.0**iN)**2 &
                 *pK**(2*iN)*cPi/2.0
         enddo
      else
         if (abs(pK).lt.pK_LIMIT2) then
            K_ell_sum_old = 0.0
            K_ell_sum     = (1.0+((1.0-pK1)/(1.0+pK1))**2/4.0) &
                 *cPi/(1.0+pK1)
            do while (abs(K_ell_sum-K_ell_sum_old).gt.DESIRED_CEI_ACCURACY)
               iN                  = iN+1
               TwoN_1FactOverNFact = TwoN_1FactOverNFact*(2.0*iN-1.0)/iN
               K_ell_sum_old       = K_ell_sum
               K_ell_sum           = K_ell_sum+(TwoN_1FactOverNFact/2.0**iN)**2 &
                    *((1.0-pK1)/(1.0+pK1))**(2*iN)*cPi/(1.0+pK1)
            enddo
         else
            K_ell_sum = alog(4.0/pK1)+(alog(4.0/pK1)-1.0)*pK1**2/4.0 &
                 + (alog(4.0/pK1)-1.0-1.0/3.0/2.0)*pK1**4*(3.0/2.0/4.0)**2 &
                 + (alog(4.0/pK1)-1.0-1.0/3.0/2.0-1.0/5.0/3.0) &
                 *pK1**6*(3.0*5.0/2.0/4.0/6.0)**2 &
                 + (alog(4.0/pK1)-1.0-1.0/3.0/2.0-1.0/5.0/3.0 &
                 - 1.0/7.0/4.0)*pK1**8*(3.0*5.0*7.0/2.0/4.0/6.0/8.0)**2
         endif
      endif
      K_elliptic = K_ell_sum

    end subroutine calc_elliptic_int_1kind

    !====================================================================

    subroutine calc_elliptic_int_2kind(Kappa,E_elliptic)

      real, intent(in):: Kappa
      real, intent(out):: E_elliptic

      integer:: iN
      real,parameter:: pK_LIMIT1 = 0.5*cSqrtTwo
      real,parameter:: pK_LIMIT2 = 0.999
      real:: DESIRED_CEI_ACCURACY
      real:: pK,pK1,E_ell_sum,E_ell_sum_old
      real:: TwoN_1FactOverNFact,TwoN_3FactOverNFact
      !------------------------------------------------------------------
      ! Compute the complete elliptic integral of 2nd kind from the series
      ! representations given by ...
      ! see formulae 8.114.1 (for 0<k<0.701) and 8.114.2 (for 0.701=<k<1)
      ! therein::
      !
      ! The stability is ensured up to pK = 0.9993 (sin**-1=88.0deg)!!!
      !
      ! Set the desired accuracy for the integral computation::

      if (iRealPrec==1) then
         DESIRED_CEI_ACCURACY = 1.0E-15
      else
         DESIRED_CEI_ACCURACY = 1.0E-7
      endif

      ! Initialize some variables::

      iN                  = 1
      TwoN_1FactOverNFact = 1.0
      TwoN_3FactOverNFact = 1.0
      pK                  = Kappa
      pK1                 = sqrt(1.0-pK**2)

      ! Compute the CEI of 2nd kind::

      if (abs(pK).lt.pK_LIMIT1) then
         E_ell_sum_old = 0.0
         E_ell_sum     = (1.0-pK**2/4.0)*cPi/2.0
         do while (abs(E_ell_sum-E_ell_sum_old).gt.DESIRED_CEI_ACCURACY)
            iN                  = iN+1
            TwoN_1FactOverNFact = TwoN_1FactOverNFact*(2.0*iN-1.0)/iN
            E_ell_sum_old       = E_ell_sum
            E_ell_sum           = E_ell_sum-(TwoN_1FactOverNFact/2.0**iN)**2 &
                 /(2.0*iN-1.0)*pK**(2*iN)*cPi/2.0
         enddo
      else
         if (abs(pK).lt.pK_LIMIT2) then
            E_ell_sum_old = 0.0
            E_ell_sum     = (1.0+((1.0-pK1)/(1.0+pK1))**2/4.0) &
                 *cPi*(1.0+pK1)/4.0
            do while (abs(E_ell_sum-E_ell_sum_old).gt.DESIRED_CEI_ACCURACY)
               iN                  = iN+1
               TwoN_3FactOverNFact = TwoN_3FactOverNFact*(2.0*iN-3.0)/iN
               E_ell_sum_old       = E_ell_sum
               E_ell_sum           = E_ell_sum+(TwoN_3FactOverNFact/2.0**iN)**2 &
                    *((1.0-pK1)/(1.0+pK1))**(2*iN)*cPi*(1.0+pK1)/4.0
            enddo
         else
            E_ell_sum = 1.0+(alog(4.0/pK1)-1.0/2.0)*pK1**2/2.0 &
                 + (alog(4.0/pK1)-1.0-1.0/3.0/4.0)*pK1**4*(3.0/4.0/4.0) &
                 + (alog(4.0/pK1)-1.0-1.0/3.0/2.0-1.0/5.0/6.0) &
                 *pK1**6*(3.0/2.0/4.0)**2*5.0/6.0 &
                 + (alog(4.0/pK1)-1.0-1.0/3.0/2.0-1.0/5.0/3.0 &
                 - 1.0/7.0/8.0)*pK1**8*(3.0*5.0/2.0/4.0/6.0)**2*7.0/8.0
         endif
      endif
      E_elliptic = E_ell_sum

    end subroutine calc_elliptic_int_2kind

  end subroutine compute_TD99_FluxRope

end module EEE_ModTD99
