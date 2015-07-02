!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CRASH_ModAtomicDataMix
  use CRASH_ModIonization
  use CRASH_ModAtomicMass,ONLY : nZMax,nMixMax
  use CRASH_ModExcitationData,ONLY: nExcitation, UseExcitation,&
       UsePressureIonization
  implicit none
  SAVE
  !Interface to databases

  integer :: nMix = -1    !Number of components in the mixture

  integer,dimension(nMixMax) :: nZ_I = -1  !Atomic numbers of elements in the mixture 


  ! relative concentrations of the elements in the mixture (part of whole comprised by the element)
  real,dimension(nMixMax) :: Concentration_I=0.0 

  ! Array of ionization potentials - energy needed to create i-level ion from (i-1)-level ion
  real,dimension(1:nZMax,nMixMax) :: IonizPotential_II

  !The degeneracy ('multiplicity') of the given quantum state:
  integer,dimension(0:nExcitation-1,nExcitation,0:nZMax-1,nMixMax):: Degeneracy_IIII = 0


  !Excitation energy in eV, with respect to the ground state energy as a function of the qunatum numbers l (first index) 
  !and n (second index for a given charge state (third index) for the given component of the mixture (fourth index)
  real,dimension(0:nExcitation-1,nExcitation,0:nZMax-1,nMixMax)::ExcitationEnergy_IIII = 0.0

  !The same, for the compression-dependent part of the state energy, assuming that the energy growth may be approximated as
  !\Delta E_{l,n,i)=VirialCoeff(l,n,i)*(a_0/r_{iono})**iPressureEffectIndex, 
  !where a_0 is the Bohr radius, r_{iono} is the iono-sphere radius
  real,dimension(0:nExcitation-1,nExcitation,0:nZMax-1,nMixMax)::VirialCoeff4Energy_IIII = 0.0

  integer,parameter :: iPressureEffectIndex = 2
contains
  !========================================================================
  subroutine set_data_for_mixture(nMixIn, nZIn_I, ConcentrationIn_I)
    use ModConst,ONLY : cTiny
    use CRASH_ModExcitationData, ONLY:cExcitationEnergy_IIII,&
         cDegeneracy_IIII,cVirialCoeff4Energy_IIII
    
    integer,intent(in)::nMixIn
    integer,dimension(nMixIn),intent(in) :: nZIn_I
    real,dimension(nMixIn),intent(in) :: ConcentrationIn_I

    integer            :: iZ, iMix  ! for loops
    !--------------------------!
    Concentration_I( 1:nMixIn ) =  ConcentrationIn_I( 1:nMixIn )
    if(abs (sum(Concentration_I( 1:nMixIn ) ) - 1.0 ) > cTiny)then
       write(*,*)&
            'Wrong input - the total of relative concentrations differs from 1.0: ',&
            Concentration_I( 1:nMixIn ) 
       call CON_stop('Stop')
    end if

    if(nMixIn==nMix)then
       if( all( nZIn_I( 1:nMix ) == nZ_I( 1:nMix ) ) )return
    end if

    nMix = nMixIn
    nZ_I( 1:nMix ) = nZIn_I( 1:nMix )
    do iMix=1,nMix
       call get_ioniz_potential(nZ_I(iMix),IonizPotential_II(1:nZ_I(iMix),iMix))
      
       ExcitationEnergy_IIII(:,:,0:nZ_I(iMix)-1,iMix) = cExcitationEnergy_IIII(:,:,0:nZ_I(iMix)-1,nZ_I(iMix))
       

       Degeneracy_IIII(:,:,0:nZ_I(iMix)-1,iMix) = cDegeneracy_IIII(:,:,0:nZ_I(iMix)-1,nZ_I(iMix))
    
       VirialCoeff4Energy_IIII(:,:,0:nZ_I(iMix)-1,iMix) = cVirialCoeff4Energy_IIII(:,:,0:nZ_I(iMix)-1,nZ_I(iMix))

    end do
   
  end subroutine set_data_for_mixture
end module CRASH_ModAtomicDataMix
