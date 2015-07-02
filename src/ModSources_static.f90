!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModSources

  use ModSizeGitm
  use ModPlanet, only: nSpecies,nSpeciesTotal,nIons

  !\
  ! Sources for neutral temperature
  !/

  real, dimension(nLons, nLats, nAlts) :: &
       Conduction, NOCooling, OCooling, ElectronHeating, &
       AuroralHeating, JouleHeating, IonPrecipHeating, &
       EddyCond,EddyCondAdia,MoleConduction

  real, dimension(nLons, nLats, nAlts,nBlocksMax) :: &
       EuvHeating, eEuvHeating, PhotoElectronHeating, &
       RadCooling, RadCoolingRate, RadCoolingErgs, EuvHeatingErgs, &
       LowAtmosRadRate, UserHeatingRate

  real, dimension(nLons,nLats,nAlts,3) :: GWAccel = 0.0

  !\
  ! Reactions used in chemistry output
  ! i.e. in2p_e -->  n2+ + e
  !      ino_n  -->  no + n
  !/

  
  integer, parameter :: nReactions = 26
 real :: ChemicalHeatingSpecies(nLons, nLats, nAlts,nReactions)
  real :: ChemicalHeatingS(nReactions)
  real :: NeutralSourcesTotal(nAlts, nSpeciesTotal)
  real :: NeutralLossesTotal(nAlts, nSpeciesTotal)
  real :: ISourcesTotal(nLons,nLats,nAlts,nIons-1,nBlocksMax)
  real :: ILossesTotal(nLons,nLats,nAlts,nIons-1,nBlocksMax)

  integer, parameter ::    in2p_e = 1
  integer, parameter ::    io2p_e = 2
  integer, parameter ::    in2p_o = 3
  integer, parameter ::    inop_e = 4
  integer, parameter ::    inp_o2 = 5
  integer, parameter ::    ino_n  = 6
  integer, parameter ::    iop_o2 = 7
  integer, parameter ::    in_o2  = 8
  integer, parameter ::    io2p_n = 9
  integer, parameter ::    io2p_no= 10
  integer, parameter ::    io2p_n2= 11
  integer, parameter ::    in2p_o2= 12
  integer, parameter ::    inp_o  = 13
  integer, parameter ::    iop_n2 = 14
  integer, parameter ::    io1d_n2 = 15
  integer, parameter ::    io1d_o2 = 16
  integer, parameter ::    io1d_o  = 17
  integer, parameter ::    io1d_e  = 18
  integer, parameter ::    in2d_o2 = 19
  integer, parameter ::    iop2d_e = 20
  integer, parameter ::    in2d_o = 21
  integer, parameter ::    in2d_e = 22
  integer, parameter ::    iop2d_n2 =23
  integer, parameter ::    iop2p_e = 24
  integer, parameter ::    iop2p_o = 25
  integer, parameter ::    iop2p_n2 = 26




  !\
  ! Stuff for auroral energy deposition and ionization
  !/

  real, dimension(:), allocatable :: &
       ED_grid, ED_Energies, ED_Flux, ED_Ion, ED_Heating
  integer :: ED_N_Energies, ED_N_Alts
  real, dimension(nAlts) :: ED_Interpolation_Weight
  integer, dimension(nAlts) :: ED_Interpolation_Index

  real :: AuroralIonRateS(nLons, nLats, nAlts, nSpecies, nBlocksMax)
  real :: AuroralHeatingRate(nLons, nLats, nAlts, nBlocksMax)
  real :: IonPrecipIonRateS(nLons, nLats, nAlts, nSpecies, nBlocksMax)
  real :: IonPrecipHeatingRate(nLons, nLats, nAlts, nBlocksMax)
  real :: ChemicalHeatingRate(nLons, nLats, nAlts)

  real :: HorizontalTempSource(nLons, nLats, nAlts)

  real :: Diffusion(nLons, nLats, nAlts, nSpecies)
  real :: NeutralFriction(nLons, nLats, nAlts, nSpecies)
  real :: IonNeutralFriction(nLons, nLats, nAlts, nSpecies)

  real :: KappaEddyDiffusion(nLons, nLats, -1:nAlts+2, nBlocksMax)

contains
  !=========================================================================
  subroutine init_mod_sources
  end subroutine init_mod_sources
  !=========================================================================
  subroutine clean_mod_sources
  end subroutine clean_mod_sources
  !=========================================================================
end module ModSources
