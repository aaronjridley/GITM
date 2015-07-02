!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf


module CRASH_ModMolecularData
  use ModConst,               ONLY: cRyToEv, cElectronMass, cAtomicMass
  use CRASH_ModAtomicDataMix, ONLY: nMixMax, nMix   !Number of components in the mixture
  use CRASH_ModAtomicDataMix, ONLY: Concentration_I !Relative atomic concentrations
  use CRASH_ModAtomicDataMix, ONLY: nZ_I            !Z for each component
  implicit none

!Dissociation energies of the diatomic molecules, [eV], (for H, C, N, O only):
!The contribution to the internal energy from the dissociation energy is
! - 1eV * N_{ij} * cDissEnergy

real, parameter, dimension(8,8) :: cDissEnergy_II = reshape(  (/   &
  !   H   !   He    !   Li   !    Be  !   B   !   C   !   N   !   O   ! 
  4.5202  , 0.      , 0.     , 0.     , 0.    , 0.    , 0.    , 0.    ,&  ! 1 - H
  0.      , 0.      , 0.     , 0.     , 0.    , 0.    , 0.    , 0.    ,&  ! 2 - He
  0.      , 0.      , 0.     , 0.     , 0.    , 0.    , 0.    , 0.    ,&  ! 3 - Li
  0.      , 0.      , 0.     , 0.     , 0.    , 0.    , 0.    , 0.    ,&  ! 4 - Be
  0.      , 0.      , 0.     , 0.     , 0.    , 0.    , 0.    , 0.    ,&  ! 5 - B
  3.5141  , 0.      , 0.     , 0.     , 0.    , 6.2907, 0.    , 0.    ,&  ! 6 - C
  3.6876  , 0.      , 0.     , 0.     , 0.    , 7.5922, 9.8030, 0.    ,&  ! 7 - N
  3.8829  , 0.      , 0.     , 0.     , 0.    , 11.163, 0.    , 5.1673 /),&  ! 8 - O
                                                                        (/8,8/))
!Moment of inertia of diatomic molecules, M*(a/cBohrRadius)**2

real, parameter, dimension(8,8) :: cMomentOfInertiaTab_II = reshape(  (/   &
  !   H   !    He   !    Li  !    Be  !   B   !   C   !   N   !   O   ! 
  0.9894  , 0.      , 0.     , 0.     , 0.    , 0.    , 0.    , 0.    ,&  ! 1 - H
  0.      , 0.      , 0.     , 0.     , 0.    , 0.    , 0.    , 0.    ,&  ! 2 - He
  0.      , 0.      , 0.     , 0.     , 0.    , 0.    , 0.    , 0.    ,&  ! 3 - Li
  0.      , 0.      , 0.     , 0.     , 0.    , 0.    , 0.    , 0.    ,&  ! 4 - Be
  0.      , 0.      , 0.     , 0.     , 0.    , 0.    , 0.    , 0.    ,&  ! 5 - B
  4.1688  , 0.      , 0.     , 0.     , 0.    , 33.107, 0.    , 0.    ,&  ! 6 - C
  3.6682  , 0.      , 0.     , 0.     , 0.    , 31.758, 30.170, 0.    ,&  ! 7 - N
  3.1914  , 0.      , 0.     , 0.     , 0.    , 31.191, 35.364, 41.648 /),&  ! 8 - O
                                                                        (/8,8/))


!Quantum of oscillation, [eV]

real, parameter, dimension(8,8) :: cOscQuantum_II = reshape(  (/   &
  !   H   !    He   !    Li  !    Be  !   B   !   C   !   N   !   O   ! 
  0.5456  , 0.      , 0.     , 0.     , 0.    , 0.    , 0.    , 0.    ,&  ! 1 - H
  0.      , 0.      , 0.     , 0.     , 0.    , 0.    , 0.    , 0.    ,&  ! 2 - He
  0.      , 0.      , 0.     , 0.     , 0.    , 0.    , 0.    , 0.    ,&  ! 3 - Li
  0.      , 0.      , 0.     , 0.     , 0.    , 0.    , 0.    , 0.    ,&  ! 4 - Be
  0.      , 0.      , 0.     , 0.     , 0.    , 0.    , 0.    , 0.    ,&  ! 5 - B
  0.3544  , 0.      , 0.     , 0.     , 0.    , 0.2299, 0.    , 0.    ,&  ! 6 - C
  0.3970  , 0.      , 0.     , 0.     , 0.    , 0.2565, 0.2924, 0.    ,&  ! 7 - N
  0.1528  , 0.      , 0.     , 0.     , 0.    , 0.2690, 0.2360, 0.1959 /),&  ! 8 - O
                                                                        (/8,8/))


 !Concentrations of diatomic molecules, [m-3]
 real, dimension(nMixMax,nMixMax) :: DiAtomicConcentration_II=-1.0

 logical:: UseDiatomicMolecules = .false.
 public:: get_chemical_equilibrium !Calculate concentrations of diatomic molecules
  
 ! logical::DoInit = .true.
contains
  subroutine get_chemical_equilibrium( Te, Na, NeutralFraction_I)
    !-------------------------------------------
    !For given \sigma_I (see the comment below) the procedure
    !calculates the array of molecular concentrations
    
    !Input parameters:
    real, intent(in):: Te    !Temperature in [eV]
    real, intent(in):: Na    !Total atomic concentrations

    !\
    !For each component, the fraction of non-ionized atom (tends to 0ne as the temperature drops. 
    !Therefore, the total atomic concentration of "i" component, which is denoted as \sigma_i in the
    !write-up equals: \sigma_{iMix} = Na * Concentration_I(iMix) *  NeutralFraction_I(iMix):
    !/
    real, intent(in):: NeutralFraction_I(1:nMix)

    real :: Sigma_I(nMix)

    !---------------------
    Sigma_I(1:nMix) =  Na * Concentration_I(1:nMix) *  NeutralFraction_I(1:nMix)
    !.............to be continued............
   
   


  end subroutine get_chemical_equilibrium
  !======================================
!  subroutine get_moment_of_inertia
!    integer :: iZ, jZ
!    if(.not.DoInit)return
!    cMomentOfInertia_II = 0.0
!
!    do jZ=1,8
!       do iZ=1,jZ
!           if(iZ == jZ) then
!           cMomentOfInertia_II(1:iZ,jZ) = 0.5/cRyToEv*cMomentOfInertia_II( 1:iZ, jZ); else
!           cMomentOfInertia_II(1:iZ,jZ) = 1/cRyToEv*cMomentOfInertia_II( 1:iZ, jZ); end if;
!       end do;
!    end do
!    DoInit = .false.
!  end subroutine get_moment_of_inertia;


end module CRASH_ModMolecularData


