!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
 Module ModChemistry

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModInputs, only: iDebugLevel, UseIonChemistry, UseNeutralChemistry
  use ModConstants
  use ModSources, only: ChemicalHeatingS, IonPrecipIonRates,AuroralIonRates
  
  implicit None
    
  real :: Ions(nIons),Neutrals(nSpeciesTotal), eDensity
  logical :: useNeutralConstituent(nSpeciesTotal),useIonConstituent(nIons)

contains
!-----------------------------------
subroutine calc_reaction_rates(iLon,iLat,iAlt,iBlock)
  
  integer, intent(in) :: iLon,iLat,iAlt,iBlock

  return

  ! Don't do anything else....

end subroutine calc_reaction_rates

subroutine calc_chemical_sources(iLon,iLat,iAlt,iBlock,IonSources, &
       IonLosses,NeutralSources, &
       NeutralLosses,ChemicalHeatingSub,Emission)

  use ModSources

  integer, intent(in) :: ilon,ilat,ialt,iBlock
  real, intent(out) :: IonSources(nIons),IonLosses(nIons)
  real, intent(out) :: NeutralSources(nSpeciesTotal), NeutralLosses(nSpeciesTotal)
  real, intent(out) :: ChemicalHeatingSub,Emission(nEmissions)
  
  real :: Source, Reaction, rr
 
  integer :: iNeutral, iion
  
  real :: lon, lat, alt

  !---------------------------------------------------------------------------
 
  IonSources = 0.0
  NeutralSources = 0.0
  IonLosses  = 0.0
  NeutralLosses = 0.0
  
  ChemicalHeatingSub = 0.0
  Emission = 0.0

  rr=EuvDissRateS(iLon,iLat,iAlt,iH2_,iBlock)
  
  Reaction = rr * &
       Neutrals(iH2_)
  
  NeutralLosses(iH2_) = NeutralLosses(iH2_) + Reaction
  NeutralSources(iH_) = NeutralSources(iH_) + 2*Reaction

  ! Solar EUV
  
  Reaction = EuvIonRateS(iLon,iLat,iAlt,iH2P_,iBlock) * &
       Neutrals(iH2_)
  
  IonSources(iH2P_)   = IonSources(iH2P_)   + Reaction
  NeutralLosses(iH2_) = NeutralLosses(iH2_) + Reaction
  
    
end subroutine calc_chemical_sources



end Module ModChemistry

