!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModMagnetogram
  use ModNumConst
  use ModMpi
  use ModIoUnit,   ONLY: io_unit_new
  use CON_axes,    ONLY: dLongitudeHgrDeg
  implicit none
  save

  private

  !\
  ! PFSSM related control parameters
  !/

  ! Maximum order of spherical harmonics
  integer,parameter::nHarmonicsMax=180 !90
  ! Weights of the spherical harmonics
  real, dimension(nHarmonicsMax+1,nHarmonicsMax+1):: g_nm,h_nm

  integer:: N_PFSSM=nHarmonicsMax

  ! Number of header lines in the file
  integer:: iHead_PFSSM=12

  ! Name of the input file
  character (LEN=32):: File_PFSSM='mf.dat'

  ! Name of output directory
  character(len=32):: NameOutDir

  ! Rs - radius of outer source surface where field is taken to be 0.
  ! Ro - radius of inner boundary for the potential
  ! H  - height of ??
  !  

  real, public :: Rs_PFSSM=2.50,Ro_PFSSM=1.0,H_PFSSM=0.0

  ! Units of the magnetic field in the file including corrections
  ! relative to the magnetic units in the SC (Gauss)
  real, public :: UnitB=1.0

  ! Rotation angle around z-axis, in degrees,
  ! from the coordinate system of the component
  ! towards the coordinate system of the magnetic 
  ! map in the positive direction - hence, is positive. 
  ! Is set automatically to be equal to the 
  ! H(eliographic) L(ongitude) of C(entral) M(eridian)
  ! of the M(ap) minus 180 Deg, if in the input file 
  ! Phi_Shift is negative
  !
  real, public :: Phi_Shift=-1.0


  !Global arrays at the magnetogram grid
  integer, parameter, public :: nDim=3,R_=1,Phi_=2,Theta_=3

  real,allocatable,dimension(:,:,:,:)::B_DN

  ! logical for reading the potential field B_DN
  logical :: DoSavePotentialField = .false.
  character(len=100) :: NamePotentialFieldFile

  real, public :: dR=1.0,dPhi=1.0,dSinTheta=1.0,dInv_D(nDim)=1.0
  integer, public :: nThetaPerProc,nRExt=2
  integer, public :: nR=29, nPhi=72, nTheta=29


  ! private variables for MPI
  integer :: iProc, nProc, iComm

  ! write prefix for magnetogram processing
  character(len=*), parameter :: prefix=''


  ! public available procedures
  public :: read_magnetogram_file
  !Then, calls set_magnetogram

  public :: set_parameters_magnetogram

  public :: set_magnetogram
  !Reads the file of magnetic field harmonics  
  !and recovers the spatial distribution of the 
  !potential mganetic field at the spherical grid 
  !nR*nPhi*nTheta

  public :: get_hlcmm
  ! Read H(eliographic) L(ongitude) of the C(entral) M(eridian) of 
  ! the M(ap) from the file header. Assign Phi_Shift=HLCMM-180

  public :: get_magnetogram_field
  !Gives the interpolated values of the Cartesian components of
  !the macnetic vector in HGR system, input parameters
  !being the cartesian coordinates in the HGR system

  public :: sin_latitude, r_latitude, colatitude
  public :: correct_angles
  public :: interpolate_field

  ! read the potential field source surface solution
  public :: read_potential_field

contains
  !=================================================================
  real function sin_latitude(iTheta)
    !\
    !Uniform in sin(latitude) grid enumerated by index iTheta
    !iTheta varies from 0 to nTheta
    !dSinTheta = 2.0/(nTheta+1)
    !/ 
    integer,intent(in)::iTheta
    sin_latitude=(real(iTheta)+0.5)*dSinTheta-1.0
  end function sin_latitude
  !=================================================================
  real function r_latitude(iTheta)
    integer,intent(in)::iTheta
    r_latitude=asin(sin_latitude(iTheta))
  end function r_latitude
  !=================================================================
  real function colatitude(iTheta)
    integer,intent(in)::iTheta
    colatitude=0.5*cPi-r_latitude(iTheta)
  end function colatitude
  !=================================================================
  ! SUBROUTINE get_hlcmm
  ! Read H(eliographic) L(ongitude) of the C(entral) M(eridian) of 
  ! the M(ap) from the header line. Assign Phi_Shift=HLCMM-180
  subroutine get_hlcmm(Head_PFSSM,Shift)
    character (LEN=80),intent(inout):: Head_PFSSM
    real,intent(inout)::Shift
    real::HLCMM     !Heliographic Longitudee of the central meridian of map
    integer::iHLCMM !The same, but HLCMM is integer at WSO magnetograms
    integer::iErrorRead,iPosition
    iPosition=index(Head_PFSSM,'Centered')	
    if (iPosition>0.and.Shift<0.0)then	
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       iPosition=index(Head_PFSSM,':')
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       iPosition=index(Head_PFSSM,':')
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       read(Head_PFSSM,'(i3)',iostat=iErrorRead)iHLCMM
       if(iErrorRead>0)call Con_stop(&
            'Can nod find HLCMM, '//File_PFSSM//&
            ' is not a true WSO magnetogram')
       Shift=modulo(iHLCMM-180-dLongitudeHgrDeg, 360.0) 
       if(iProc==0)then
          write(*,*) prefix, 'Phi_Shift=',Shift
       end if
       return
    end if
    iPosition=index(Head_PFSSM,'Central')
    if(iPosition>0.and.Shift<0.0)then
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       iPosition=index(Head_PFSSM,':')
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       read(Head_PFSSM,*,iostat=iErrorRead)HLCMM
       if(iErrorRead>0)call CON_stop(&
            'Can nod find HLCMM, '//File_PFSSM//&
            ' is not a true MDI magnetogram')
       Shift=modulo(HLCMM-180-dLongitudeHgrDeg, 360.0) 
       if(iProc==0)then
          write(*,*) prefix, 'Phi_Shift=',Shift
       end if
    end if
  end subroutine get_hlcmm

  !============================================================================

  subroutine set_parameters_magnetogram(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand

    character(len=*), parameter :: &
         NameSub = 'ModMagnetogram::set_parameters_magnetogram'
    !--------------------------------------------------------------------------

    select case(NameCommand)
    case("#MAGNETOGRAM")
       call read_var('rMagnetogram',        Ro_PFSSM)
       call read_var('rSourceSurface',      Rs_PFSSM)
       call read_var('HeightInnerBc',       H_PFSSM)
       call read_var('NameMagnetogramFile', File_PFSSM)
       call read_var('nHeaderLine',         iHead_PFSSM)
       call read_var('PhiShift',            Phi_Shift)
       call read_var('UnitB',               UnitB)

    case("#B0GRID")
       call read_var('nR', nR)

    case("#SAVEPOTENTIALFIELD")
       call read_var('DoSavePotentialField', DoSavePotentialField)

    case("#READPOTENTIALFIELD")
       call read_var('NamePotentialFieldFile', NamePotentialFieldFile)
       call read_var('HeightInnerBc',          H_PFSSM)
       call read_var('UnitB',                  UnitB)

    case default
       call CON_stop(NameSub//' invalid NameCommand='//NameCommand)
    end select

  end subroutine set_parameters_magnetogram

  !============================================================================

  subroutine read_magnetogram_file(NamePlotDir)
    implicit none

    character(len=*),intent(in) :: NamePlotDir
    integer :: iError
    !--------------------------------------------------------------------------
    iComm = MPI_COMM_WORLD
    call MPI_COMM_SIZE(iComm,nProc,iError)
    call MPI_COMM_RANK(iComm,iProc,iError)

    if (iProc==0) then
       write(*,*) prefix, 'Norder = ',N_PFSSM
       write(*,*) prefix, 'Entered coefficient file name :: ',File_PFSSM
       write(*,*) prefix, 'Entered number of header lines:: ',iHead_PFSSM
    endif
    !\
    ! Initialize once g(n+1,m+1) & h(n+1,m+1) by reading a file
    ! created from Web data::
    !/ 
    NameOutDir = NamePlotDir

    call read_harmonics

    call set_magnetogram

  end subroutine read_magnetogram_file

  !============================================================================

  subroutine read_harmonics
    integer :: iUnit,nOrderIn,iPosition,iPosition1,n,m,i,iError
    character (LEN=80):: Head_PFSSM=''
    real::gtemp,htemp,stuff1,stuff
    !\
    ! Formats adjusted for wso CR rad coeffs::
    !/
    iUnit = io_unit_new()
    open(iUnit,file=File_PFSSM,status='old',iostat=iError)
    if(iError>0)call CON_stop('Cannot open '//File_PFSSM)
    if (iHead_PFSSM /= 0) then
       do i=1,iHead_PFSSM
          read(iUnit,'(a)') Head_PFSSM
          iPosition=index(Head_PFSSM,'rder')
          iPosition1=0
          if(iPosition>0)&
               iPosition1=max(&
               index(Head_PFSSM(iPosition+4:iPosition+9),'='),&
               index(Head_PFSSM(iPosition+4:iPosition+9),':'))
          if(iPosition1>0)then
             read(Head_PFSSM(iPosition+4+iPosition1:len(Head_PFSSM)),&
                  '(i3)',iostat=iError)nOrderIn
             if(iError>0)call CON_stop('Cannot figure out nOrder')
             if(nOrderIn<N_PFSSM)then
                if(iProc==0)then
                   write(*,*) prefix, 'Reset N_PFSSM=',nOrderIn
                end if
                N_PFSSM=nOrderIn
             end if
          end if
          if(Phi_Shift<0.0)call get_hlcmm(Head_PFSSM,Phi_Shift)
       enddo
    endif
    nR=max(nR,N_PFSSM)
    nPhi=max(nPhi,N_PFSSM)
    nTheta=max(nTheta,N_PFSSM)
    nRExt=min(10,1+floor(real(nR)*(0.1*Ro_PFSSM)/(Rs_PFSSM-Ro_PFSSM)))
    if(iProc==0)then
       write(*,*) prefix, &
            'Magnetogram is extended by ',nRExt,' nodes towards the Sun'
    end if
    !\
    ! Initialize all coefficient arrays::
    !/
    g_nm = 0.0; h_nm = 0.0
    !\
    ! Read file with coefficients, g_nm and h_nm::
    !/
    do
       read(iUnit,*,iostat=iError) n,m,gtemp,htemp
       if (iError /= 0) EXIT
       if (n > N_PFSSM .or. m > N_PFSSM) CYCLE
       g_nm(n+1,m+1) = gtemp
       h_nm(n+1,m+1) = htemp
    enddo
    close(iUnit)
    !\
    ! Add correction factor for radial, not LOS, coefficients::
    ! Note old "coefficients" file are LOS, all new coeffs and 
    ! files are radial)
    !/
    stuff=Rs_PFSSM
    do n=0,N_PFSSM
       stuff=stuff/Rs_PFSSM**2
       stuff1 = 1.0/(real(n+1)+real(n)*stuff)
       !       stuff1 = 1.0/real(n+1+(n/(Rs_PFSSM**(2*n+1))))
       g_nm(n+1,1:n+1) = g_nm(n+1,1:n+1)*stuff1
       h_nm(n+1,1:n+1) = h_nm(n+1,1:n+1)*stuff1
    enddo
    !\
    ! Leave out monopole (n=0) term::
    !/
    g_nm(1,1) = 0.0
  end subroutine read_harmonics
  !===========================================================================
  subroutine set_magnetogram
    !
    !---------------------------------------------------------------------------
    ! This subroutine computes PFSS (Potential Field Source Surface)
    ! field model components in spherical coordinates at user-specified
    ! r (solar radii units, r>1) and theta, phi (both radians).
    ! The subroutine requires a file containing the spherical
    ! harmonic coefficients g(n,m) & h(n.m) obtained from a separate analysis 
    ! of a photospheric synoptic map, in a standard ("radial")
    ! format and normalization used by Stanford.
    !
    ! The PFSS field model assumes no currents in the corona and
    ! a pure radial field at the source surface, here R=2.5 Rsun
    !
    ! Get solar coefficients from Todd Hoeksema's files:
    !    1. Go to http://solar.stanford.edu/~wso/forms/prgs.html
    !    2. Fill in name and email as required
    !    3. Chose Carrington rotation (leave default 180 center longitude)
    ! For most requests of integer CRs with order < 20, result will come back
    ! immediately on the web.
    !    4. Count header lines before 1st (0,0) coefficient -this will be asked!
    !---------------------------------------------------------------------------
    ! Notes:
    !
    ! In the calling routine you must initialize one variable: istart=0 (it is a 
    ! flag used to tell the subroutine to read the coefficient file the first 
    ! time only). The first time around (DoFirst=0), the subroutine will ask for
    ! the coefficient file name, the order of the expansion to use (N_PFSSM=40 or 
    ! less*, but the coeff file can contain more orders than you use), and the 
    ! number of lines in the coefficient file header. (*note computation time 
    ! increases greatly with order used).
    !
    ! The source surface surface radius has been set at Rs=2.5*Ro in the 
    ! subroutine. PFSS fields at R>Rs are radial.(br,bthet,bphi) are the resulting
    ! components. Note the units of the B fields will differ with observatory used
    ! for the coefficients. Here we assume use of the wso coefficients so units are
    ! microT. The computation of the B fields is taken mainly from Altschuler, 
    ! Levine, Stix, and Harvey, "High Resolutin Mapping of the Magnetic Field of
    ! the Solar Corona," Solar Physics 51 (1977) pp. 345-375. It uses Schmidt
    ! normalized Legendre polynomials and the normalization is explained in the 
    ! paper. The field expansion in terms of the Schmidt normalized Pnm and dPnm's
    ! is best taken from Todd Hoeksema's notes which can be downloaded from the Web
    ! http://quake.stanford.edu/~wso/Description.ps
    ! The expansions  used to get include radial factors to make the field become
    ! purely radial at the source surface. The expans. in Altschuler et al assumed
    ! that the the coefficient g(n,m) and h(n,m) were the LOS coefficients -- but 
    ! the g(n,m) and h(n,m) now available are radial (according to Janet Luhman). 
    ! Therefore, she performs an initial correction to the g(n,m) and h(n,m) to 
    ! make them consistent with the the expansion. There is no reference for this
    ! correction.
    !--------------------------------------------------------------------------  
    integer:: i,n,m,iTheta,iPhi,iR
    real:: c_n
    real:: SinPhi,CosPhi
    real:: SinPhi_I(0:N_PFSSM),CosPhi_I(0:N_PFSSM)
    real:: CosTheta,SinTheta
    real:: stuff1,stuff2,stuff3
    real:: SumR,SumT,SumP,SumPsi
    real:: Theta,Phi,R_PFSSM
    real:: Psi_PFSSM
    integer::iBcast, iStart, nSize, iError
    integer::iNorth,isouth

    ! Temporary variable
    real, dimension(N_PFSSM+1):: FactRatio1
    real, dimension(N_PFSSM+1,N_PFSSM+1):: p_nm,dp_nm

    integer, parameter:: MaxInt=100000
    real:: Sqrt_I(MaxInt)
    !\
    !
    !/

    real, dimension(-1:N_PFSSM+2,-nRExt:nR):: &
         RoRsPower_I, RoRPower_I, rRsPower_I


    !Introduce a spherical grid with the resolution, depending on the
    !magnetogram resolution (N_PFSSM)
    dR=(Rs_PFSSM-Ro_PFSSM)/real(nR)
    dPhi=cTwoPi/real(nPhi)
    dSinTheta=2.0/real(nTheta+1) 
    !sin of latitude is taken in 29+1 points in the MWO file
    dInv_D=1.0/(/dR,dPhi,dSinTheta/)

    !Calculate the radial part of spherical functions
    call calc_radial_functions

    call set_auxiliary_arrays

    !Allocate the magnetic field array, at the spherical grid.
    if(allocated(B_DN))deallocate(B_DN)
    allocate(B_DN(R_:Theta_,-nRExt:nR,0:nPhi,0:nTheta))

    B_DN=0.0


    !Parallel computation of the magnetic field at the grid
    nThetaPerProc=nTheta/nProc+1 


    !Loop by theta, each processor treats a separate part of the grid
    do iTheta=iProc*nThetaPerProc,(iProc+1)*nThetaPerProc-1

       if(iTheta>nTheta)EXIT !Some processors may have less amount of work
       Theta=colatitude(iTheta)
       CosTheta=cos(Theta)
       SinTheta=max(sin(Theta), 1.0E-10)

       !Calculate the set of Legandre polynoms, for given CosTheta,SinTheta
       call calc_Legandre_polynoms

       !Start loop by Phi
       do iPhi=0,nPhi
          Phi=real(iPhi)*dPhi

          !Calculate azymuthal harmonics, for a given Phi
          do m=0,N_PFSSM
             CosPhi_I(m)=cos(m*Phi)
             SinPhi_I(m)=sin(m*phi)
          end do

          !Loop by radius
          do iR=-nRExt,nR
             R_PFSSM=Ro_PFSSM+dR*iR
             !\
             ! Initialize the values of SumR,SumT,SumP, and SumPsi::
             !/
             SumR = 0.0; SumT   = 0.0
             SumP = 0.0; SumPsi = 0.0

             !\
             ! Calculate B for (R_PFSSM,Phi_PFSSM)::
             ! Also calculate magnetic potential Psi_PFSSM
             !/
             do m=0,N_PFSSM
                !if((iTheta==0.or.iTheta==nTheta).and.m>1)CYCLE
                do n=m,N_PFSSM
                   !\
                   ! c_n corresponds to Todd's c_l::
                   !/
                   c_n    = -RoRsPower_I(n+2,iR)
                   !\
                   ! Br_PFSSM = -d(Psi_PFSSM)/dR_PFSSM::
                   !/
                   stuff1 = (n+1)*RoRPower_I(n+2,iR)-c_n*n*rRsPower_I(n-1,iR)
                   stuff2 = g_nm(n+1,m+1)*CosPhi_I(m)+h_nm(n+1,m+1)*SinPhi_I(m)
                   SumR   = SumR + p_nm(n+1,m+1)*stuff1*stuff2
                   !\
                   ! Bt_PFSSM = -(1/R_PFSSM)*d(Psi_PFSSM)/dTheta_PFSSM::
                   !/
                   stuff1 = RoRPower_I(n+2,iR)+c_n*rRsPower_I(n-1,iR)
                   SumT   = SumT-dp_nm(n+1,m+1)*stuff1*stuff2
                   !\
                   ! Psi_PFSSM::
                   !/
                   SumPsi = SumPsi+R_PFSSM*p_nm(n+1,m+1)*stuff1*stuff2
                   !\
                   ! Bp_PFSSM = -(1/R_PFSSM)*d(Psi_PFSSM)/dPhi_PFSSM::
                   !/
                   stuff2 = g_nm(n+1,m+1)*SinPhi_I(m)-h_nm(n+1,m+1)*CosPhi_I(m)
                   SumP   = SumP + p_nm(n+1,m+1)*m/SinTheta*stuff1*stuff2
                enddo
             enddo
             !\
             ! Compute (Br_PFSSM,Btheta_PFSSM,Bphi_PFSSM) and Psi_PFSSM::
             !/
             Psi_PFSSM    = SumPsi
             B_DN(R_,iR,iPhi,iTheta)     = SumR
             B_DN(Phi_,iR,iPhi,iTheta)   = SumP 
             B_DN(Theta_,iR,iPhi,iTheta) = SumT
          end do
       end do
    end do

    if(nProc>1)then
       do iBcast=0,nProc-1
          iStart=iBcast*nThetaPerProc
          if(iStart>nTheta)EXIT
          nSize=min(nThetaPerProc,nTheta+1-iStart)*(nPhi+1)*&
               (nR+1+nRExt)*3
          call MPI_bcast(B_DN(1,-nRExt,0,iStart),nSize,MPI_REAL,iBcast,iComm,iError)
       end do
    end if

    if(iProc==0)call write_Br_plot
    if(iProc==0 .and. DoSavePotentialField)call save_potential_field
  contains
    subroutine set_auxiliary_arrays

      !\
      ! Optimization by G. Toth::
      !/
      !\
      ! Calculate sqrt(integer) from 1 to 10000::
      !/
      do m=1,MaxInt
         Sqrt_I(m) = sqrt(real(m))
      end do
      !\
      ! Calculate the ratio sqrt(2m!)/(2^m*m!)::
      !/
      factRatio1(:) = 0.0; factRatio1(1) = 1.0
      do m=1,N_PFSSM
         factRatio1(m+1) = factRatio1(m)*Sqrt_I(2*m-1)/Sqrt_I(2*m)
      enddo
    end subroutine set_auxiliary_arrays

    subroutine calc_Legandre_polynoms
      real:: SinThetaM, SinThetaM1
      integer:: delta_m0
      !\
      ! Calculate polynomials with appropriate normalization
      ! for Theta_PFSSMa::
      !/
      SinThetaM  = 1.0
      SinThetaM1 = 1.0
      p_nm(:,:) = 0.0; dp_nm(:,:) = 0.0

      do m=0,N_PFSSM
         if (m == 0) then
            delta_m0 = 1
         else
            delta_m0 = 0
         endif
         !\
         ! Eq.(27) from Altschuler et al. 1976::
         !/
         p_nm(m+1,m+1) = factRatio1(m+1)*Sqrt_I((2-delta_m0)*(2*m+1))* &
              SinThetaM
         !\
         ! Eq.(28) from Altschuler et al. 1976::
         !/
         if (m < N_PFSSM) p_nm(m+2,m+1) = p_nm(m+1,m+1)*Sqrt_I(2*m+3)* &
              CosTheta
         !\
         ! Eq.(30) from Altschuler et al. 1976::
         !/
         dp_nm(m+1,m+1) = factRatio1(m+1)*Sqrt_I((2-delta_m0)*(2*m+1))*&
              m*CosTheta*SinThetaM1
         !\
         ! Eq.(31) from Altschuler et al. 1976::
         !/
         if (m < N_PFSSM) &
              dp_nm(m+2,m+1) = Sqrt_I(2*m+3)*(CosTheta*&
              dp_nm(m+1,m+1)-SinTheta*p_nm(m+1,m+1))

         SinThetaM1 = SinThetaM
         SinThetaM  = SinThetaM*SinTheta

      enddo
      do m=0,N_PFSSM-2; do n=m+2,N_PFSSM
         !\
         ! Eq.(29) from Altschuler et al. 1976::
         !/
         stuff1         = Sqrt_I(2*n+1)/Sqrt_I(n**2-m**2)
         stuff2         = Sqrt_I(2*n-1)
         stuff3         = Sqrt_I((n-1)**2-m**2)/Sqrt_I(2*n-3)
         p_nm(n+1,m+1)  = stuff1*(stuff2*CosTheta*p_nm(n,m+1)-  &
              stuff3*p_nm(n-1,m+1))
         !\         ! Eq.(32) from Altschuler et al. 1976::
         !/
         dp_nm(n+1,m+1) = stuff1*(stuff2*(CosTheta*dp_nm(n,m+1)-&
              SinTheta*p_nm(n,m+1))-stuff3*dp_nm(n-1,m+1))
      enddo; enddo
      !\
      ! Apply Schmidt normalization::
      !/
      do m=0,N_PFSSM; do n=m,N_PFSSM
         !\
         ! Eq.(33) from Altschuler et al. 1976::
         !/
         stuff1 = 1.0/Sqrt_I(2*n+1)
         !\
         ! Eq.(34) from Altschuler et al. 1976::
         !/
         p_nm(n+1,m+1)  = p_nm(n+1,m+1)*stuff1
         dp_nm(n+1,m+1) = dp_nm(n+1,m+1)*stuff1
      enddo; enddo
    end subroutine calc_Legandre_polynoms
    subroutine calc_radial_functions
      do iR=-nRExt,nR
         !\ 
         ! Calculate powers of the ratios of radii
         !/
         R_PFSSM=Ro_PFSSM+dR*iR
         rRsPower_I(-1,iR) = Rs_PFSSM/R_PFSSM 
         ! This one can have negative power.
         rRsPower_I(0,iR)  = 1.0
         RoRsPower_I(0,iR) = 1.0
         RoRPower_I(0,iR)  = 1.0
         do m=1,N_PFSSM+2
            RoRsPower_I(m,iR) = RoRsPower_I(m-1,iR) * (Ro_PFSSM/Rs_PFSSM)
            RoRPower_I(m,iR)  = RoRPower_I(m-1,iR)  * (Ro_PFSSM/R_PFSSM)
            rRsPower_I(m,iR)  = rRsPower_I(m-1,iR)  * (R_PFSSM /Rs_PFSSM)
         end do
      end do
    end subroutine calc_radial_functions

  end subroutine set_magnetogram

  !=====================================================================

  subroutine write_Br_plot
    use ModPlotFile, ONLY: save_plot_file

    integer :: iError,iPhi,iTheta,iUnit
    real,dimension(2,0:nPhi,0:nTheta):: Coord_DII,State_VII
    character(len=32)                :: FileNameDat
    character(len=32)                :: FileNameOut
    !-----------------------------------------------------
    FileNameDat = trim(NameOutDir)//'PFSSM_Br.dat'
    FileNameOut = trim(NameOutDir)//'PFSSM_Br.outs'

    write(*,*) prefix, 'Writing PFSSM_Br output file, named'
    write(*,*) prefix, FileNameDat
    iUnit = io_unit_new()
    open ( unit = iUnit, &
         file = FileNameDat, &
         form = 'formatted', &
         access = 'sequential', &
         status = 'replace', iostat = iError )

    if ( iError /= 0 ) then
       write (*,*) prefix, ' '
       write (*,*) prefix, 'TECPLOT_WRITE_OPEN - Fatal error!'
       write (*,*) prefix, '  Could not open the PFSSM_Br output file.'
       call CON_stop('')
    end if

    write ( iUnit, '(a)' ) 'Title = "'     // trim ('PFSSM_Br') // '"'
    write ( iUnit, '(a)' ) &
         'Variables = ' // trim (&
         '"Longitude [Deg]", "Latitude [Deg]", "Br_0 [G]","Br_SS [G]"')
    write ( iUnit, '(a)' ) ' '
    write ( iUnit, '(a,i6,a,i6,a)' ) 'Zone I = ', nPhi+1, ', J=', nTheta+1,&
         ', F=point' 

    do iTheta=0,nTheta
       do iPhi=0,nPhi
          Coord_DII(1,iPhi,iTheta) = real(iPhi)*dPhi/cDegToRad
          Coord_DII(2,iPhi,iTheta) = r_latitude(iTheta)/cDegToRad

          State_VII(1,iPhi,iTheta) = UnitB*B_DN(R_,0,iPhi,iTheta)
          State_VII(2,iPhi,iTheta) = UnitB*B_DN(R_,nR,iPhi,iTheta)

          write ( iUnit, '(4f10.3)' )Coord_DII(:,iPhi,iTheta),&
               State_VII(:,iPhi,iTheta)
       end do
    end do
    close(iUnit)
    call save_plot_file(FileNameOut,TypeFileIn='ascii',&
         StringHeaderIn=&
         'Longitude [Deg], Latitude [Deg], Br_0 [G], Br_SS [G]',&
         nDimIn=2, CoordIn_DII=Coord_DII,VarIn_VII=State_VII)

  end subroutine Write_Br_plot

  !==========================================================================

  subroutine save_potential_field

    use ModPlotFile, ONLY: save_plot_file

    real, allocatable :: Radius_I(:), Theta_I(:), Phi_I(:)

    integer :: iR, iTheta, iPhi
    character(len=32) :: FileNameOut
    !------------------------------------------------------------------------

    allocate(Radius_I(-nRExt:nR), Phi_I(0:nPhi), Theta_I(0:nTheta))

    do iR = -nRExt, nR
       Radius_I(iR) = Ro_PFSSM + dR*iR
    end do
    do iPhi = 0, nPhi
       Phi_I(iPhi) = real(iPhi)*dPhi
    end do
    do iTheta = 0, nTheta
       Theta_I(iTheta) = r_latitude(iTheta)
    end do

    FileNameOut = trim(NameOutDir)//'potentialfield.outs'
    call save_plot_file(FileNameOut, TypeFileIn='real8', StringHeaderIn = &
         'Radius [Rs] Longitude [Deg] Latitude [Deg] B [G]', &
         nameVarIn = 'Radius Longitude Latitude Br Bphi Btheta' &
         //' Ro_PFSSM Rs_PFSSM', &
         ParamIn_I = (/ Ro_PFSSM, Rs_PFSSM /), &
         nDimIn=3, VarIn_VIII=B_DN, &
         Coord1In_I=Radius_I, &
         Coord2In_I=Phi_I, &
         Coord3In_I=Theta_I)

    deallocate(Radius_I, Theta_I, Phi_I)

  end subroutine save_potential_field

  !============================================================================

  subroutine read_potential_field(NamePlotDir)

    use ModPlotFile, ONLY: read_plot_file

    character(len=*),intent(in) :: NamePlotDir

    integer :: iError
    integer :: n_D(3), nParam
    real :: Param_I(2)

    character(len=*), parameter :: &
         NameSub = 'ModMagnetogram::read_potential_field'
    !--------------------------------------------------------------------------
    iComm = MPI_COMM_WORLD
    call MPI_COMM_SIZE(iComm,nProc,iError)
    call MPI_COMM_RANK(iComm,iProc,iError)

    n_D = 1
    call read_plot_file(NamePotentialFieldFile, TypeFileIn='real8', &
         nOut_D=n_D, nParamOut=nParam, ParamOut_I=Param_I)

    !!!write(*,*)'nParam = ',nParam

    Ro_PFSSM = Param_I(1)
    Rs_PFSSM = Param_I(2)
    Phi_Shift = 0.0
    nRExt = 0

    nR = n_D(1) - 1
    nPhi = n_D(2) - 1
    nTheta = n_D(3) - 1

    allocate(B_DN(3,-nRExt:nR,0:nPhi,0:nTheta))

    call read_plot_file(NamePotentialFieldFile, TypeFileIn='real8', &
         VarOut_VIII=B_DN)

    ! recover public vaiables
    dR = (Rs_PFSSM - Ro_PFSSM)/nR
    dPhi = cTwoPi/nPhi
    dSinTheta = 2.0/(nTheta+1)
    dInv_D = 1.0/(/dR,dPhi,dSinTheta/)
    nThetaPerProc = nTheta/nProc + 1

    NameOutDir = NamePlotDir

    if(iProc==0)call write_Br_plot
    if(iProc==0 .and. DoSavePotentialField)call save_potential_field

  end subroutine read_potential_field

  !==========================================================================
  ! This subroutine corrects the angles Phi and Theta after every 
  ! step of the field-alligned integration 
  subroutine correct_angles(TR_D)
    real,dimension(nDim),intent(inout) :: TR_D

    if(TR_D(Theta_) < 0.0)then
       TR_D(Theta_) = -TR_D(Theta_)
       TR_D(Phi_)=TR_D(Phi_)+cPi
    end if
    if(TR_D(Theta_) > cPi)then
       TR_D(Theta_)=cTwoPi-TR_D(Theta_)
       TR_D(Phi_)=TR_D(Phi_)+cPi
    end if
    TR_D(Phi_)=modulo(TR_D(Phi_),cTwoPi)

  end subroutine correct_angles
  !==========================================================================
  !\
  ! Calculate B_R, B_Phi, B_Theta in Gauss for given
  ! R, Phi, Theta, Theta being the colatitude.
  !/
  subroutine interpolate_field(R_D,BMap_D)
    real,intent(in)::R_D(nDim)
    real,intent(out)::BMap_D(nDim)
    integer::Node_D(nDim)
    real::Res_D(nDim)
    integer::iDim,i,j,k  
    real::Weight_III(0:1,0:1,0:1)
    logical::DoCorrection
    real::ReductionCoeff,BRAvr
    Res_D=R_D
    !Limit a value of R:
    Res_D(R_)=max(min(Res_D(R_),Rs_PFSSM-cTiny),Ro_PFSSM-nRExt*dR+cTiny)

    Res_D(R_)=Res_D(R_)-Ro_PFSSM

    call correct_angles(Res_D)
    DoCorrection=Res_D(Theta_)/=R_D(Theta_)
    Res_D(Theta_)=cos(Res_D(Theta_)) & !This is sin(latitude)
         -sin_latitude(0)             !This is sin(latitude) for iTheta=0 
    Res_D=Res_D*dInv_D
    Node_D=floor(Res_D)
    if(Node_D(R_)==nR)Node_D(R_)=Node_D(R_)-1
    Res_D=Res_D-real(Node_D)
    ReductionCoeff=1.0
    BRAvr=0.0
    !Near poles reduce the Phi and Theta components of the field and
    !take a weithed average for the R component from using the  value in
    !the last grid node and that averaged over longitude 
    !(iTheta=0 or iTheta=nTheta)
    if(Node_D(Theta_)==nTheta)then
       Node_D(Theta_)=Node_D(Theta_)-1
       ReductionCoeff=sqrt(max(1.0-2.0*Res_D(Theta_),0.0))
       Res_D(Theta_)=1.0
       BRAvr=sum(&
            (1.0-Res_D(R_))*B_DN(R_,Node_D(R_)  ,:,nTheta)&
            +Res_D(R_) *B_DN(R_,Node_D(R_)+1,:,nTheta))/(nPhi+1)
    elseif(Node_D(Theta_)==-1)then
       Node_D(Theta_)= Node_D(Theta_)+1
       ReductionCoeff=sqrt(max(0.0,2.0*Res_D(Theta_)-1.0))
       Res_D(Theta_) = 0.0
       BRAvr=sum(&
            (1.0-Res_D(R_))*B_DN(R_,Node_D(R_)  ,:,     0)&
            +Res_D(R_) *B_DN(R_,Node_D(R_)+1,:,     0))/(nPhi+1)
    end if

    if(Node_D(Phi_)==nPhi)Node_D(Phi_)=0

    Weight_III(0,:,:)=1.0-Res_D(R_)
    Weight_III(1,:,:)=Res_D(R_)
    Weight_III(:,0,:)=Weight_III(:,0,:)*(1.0-Res_D(Phi_))
    Weight_III(:,1,:)=Weight_III(:,1,:)*Res_D(Phi_)
    Weight_III(:,:,0)=Weight_III(:,:,0)*(1.0-Res_D(Theta_))
    Weight_III(:,:,1)=Weight_III(:,:,1)*Res_D(Theta_)
    do iDim=1,nDim
       BMap_D(iDim)=&
            sum(Weight_III*B_DN(iDim,&
            Node_D(R_):Node_D(R_)+1,&
            Node_D(Phi_):Node_D(Phi_)+1,&
            Node_D(Theta_):Node_D(Theta_)+1))
    end do
    if(DoCorrection)then
       BMap_D(Phi_:Theta_)=-ReductionCoeff*BMap_D(Phi_:Theta_)
    else
       BMap_D(Phi_:Theta_)=ReductionCoeff*BMap_D(Phi_:Theta_)
    end if
    BMap_D(R_)=ReductionCoeff*BMap_D(R_)+(1.0-ReductionCoeff)*BRAvr
  end subroutine interpolate_field
  !==========================================================================
  subroutine get_magnetogram_field(xInput,yInput,zInput,B0_D)
    !\
    ! For Cartesian xInput,yInput,zInput, all normalized per the
    ! solar radius, calculate the magnetic field vector 
    ! in Cartesian components, in SI units.
    !/
    real, intent(in):: xInput,yInput,zInput
    real, intent(out), dimension(3):: B0_D

    real:: Rin_PFSSM,Theta_PFSSM,Phi_PFSSM
    real:: BMap_D(nDim)
    real:: SinPhi,CosPhi,XY,R_PFSSM
    real:: CosTheta,SinTheta

    ! Variables to obtain FDIPS field below Ro_PFSSM
    real :: BMapo_D(nDim), BMapi_D(nDim), Ri_PFSSM
    !--------------------------------------------------------------------------
    !\
    ! Calculate cell-centered spherical coordinates::
    !/
    Rin_PFSSM   = sqrt(xInput**2+yInput**2+zInput**2)
    !\
    ! Avoid calculating B0 inside a critical radius = 0.9*Rsun
    !/
    if(nRExt/=0 .and. Rin_PFSSM < max(Ro_PFSSM-dR*nRExt,0.9*Ro_PFSSM) &
         .or. nRExt==0 .and. Rin_PFSSM < 0.9*Ro_PFSSM)then
       B0_D= 0.0
       RETURN
    end if
    Theta_PFSSM = acos(zInput/Rin_PFSSM)
    Xy          = sqrt(xInput**2+yInput**2+cTiny**2)
    Phi_PFSSM   = atan2(yInput,xInput)
    SinTheta    = Xy/Rin_PFSSM
    CosTheta    = zInput/Rin_PFSSM
    SinPhi      = yInput/Xy
    CosPhi      = xInput/Xy
    !\
    ! Set the source surface radius::
    ! The inner boundary in the simulations starts at a height
    ! H_PFSSM above that of the magnetic field measurements!
    !/
    R_PFSSM =min(Rin_PFSSM+H_PFSSM, Rs_PFSSM)

    !\
    ! Transform Phi_PFSSM from the component's frame to the magnetogram's frame
    !/
    Phi_PFSSM = Phi_PFSSM - Phi_Shift*cDegToRad

    if(nRExt == 0 .and. R_PFSSM < Ro_PFSSM)then
       ! linearly extrapolate FDIPS field for locations below Ro_PFSSM
       Ri_PFSSM = 2.0*Ro_PFSSM - R_PFSSM

       call interpolate_field((/Ri_PFSSM,Phi_PFSSM,Theta_PFSSM/),BMapi_D)
       call interpolate_field((/Ro_PFSSM,Phi_PFSSM,Theta_PFSSM/),BMapo_D)

       BMap_D = 2.0*BMapo_D - BMapi_D
    else
       call interpolate_field((/R_PFSSM,Phi_PFSSM,Theta_PFSSM/),BMap_D)
    end if
    !\
    ! Magnetic field components in global Cartesian coordinates::
    ! Set B0x::
    !/
    B0_D(1) = BMap_D(R_)*SinTheta*CosPhi+    &
         BMap_D(Theta_)*CosTheta*CosPhi-&
         BMap_D(Phi_)*SinPhi
    !\
    ! Set B0y::
    !/
    B0_D(2) =  BMap_D(R_)*SinTheta*SinPhi+    &
         BMap_D(Theta_)*CosTheta*SinPhi+&
         BMap_D(Phi_)*CosPhi
    !\
    ! Set B0z::
    !/
    B0_D(3) = BMap_D(R_)*CosTheta-           &
         BMap_D(Theta_)*SinTheta
    !\
    ! Apply field strength normalization::
    ! UnitB contains the units of the CR file relative to 1 Gauss
    ! and a possible correction factor (e.g. 1.8 or 1.7).
    !/
    B0_D=B0_D*UnitB

    if (Rin_PFSSM > Rs_PFSSM) &
         B0_D  =  B0_D*(Rs_PFSSM/Rin_PFSSM)**2

    ! Transform from Gauss to Tesla
    B0_D = B0_D*1.0E-4
  end subroutine get_magnetogram_field
end module ModMagnetogram
