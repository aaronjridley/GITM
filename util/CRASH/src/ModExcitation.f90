!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module CRASH_ModExcitation
  use CRASH_ModAtomicDataMix
  use CRASH_ModExcitationData,ONLY : n_ground
  implicit none
  PRIVATE !Except


  !The logical to handle whether energy states above continuum may be eliminated
  logical,public :: DoStateElimination = .false.


  !Excitation energy of ion of ionization state i, averaged by possible
  !excitation levels [eV]
  real,dimension(0:nZMax,nMixMax),public :: ExtraEnergyAv_II = 0.0

  !The average value of -V\frac{\partial \Delta E}{\partial V}
  !for given i and iMix [eV]
  real,dimension(0:nZMax,nMixMax),public :: VirialCoeffAv_II = 0.0

  !The average value of V^2\frac{\partial^2 \Delta E}{\partial V^2}
  !for given i and iMix [eV]
  real,dimension(0:nZMax,nMixMax),public :: SecondVirialCoeffAv_II = 0.0

  !The variance of the extra energy, E_x, for given i and iMix [eV^2]
  real,dimension(0:nZMax,nMixMax),public :: Cov2ExtraEnergy_II = 0.0

  !The covariance between the extra energy and -V\frac{\partial \Delta E}{\partial V}
  !for given i and iMix [eV^2]
  real,dimension(0:nZMax,nMixMax),public :: CovExtraEnergyVirial_II = 0.0

  !The variance of -V\frac{\partial \Delta E}{\partial V}
  !for given i and iMix [eV^2]
  real,dimension(0:nZMax,nMixMax),public :: Cov2VirialCoeff_II = 0.0


  !\logarithm of the statistical sub-sum over the excitated states, for
  ! a given sort of ions.
  real,dimension(0:nZMax,nMixMax),public :: LogGi_II = 0.0


  !\
  ! Partition function for the excited states
  !/
  real,dimension(0:nExcitation-1,nExcitation,0:nZMax,nMixMax) :: Partition_IIII = 0.0

  !Averages over l:
  real,dimension(0:nExcitation,0:nZMax  ,nMixMax),public :: Partition_III = 0.0
  real,dimension(nExcitation,0:nZMax-1,nMixMax),public :: ExcitationEnergy_III = 0.0

  integer,dimension(0:nZMax,1:nMixMax),public :: nExcitation_II = nExcitation


  public :: get_excitation_levels,&
            get_excitation_levels_zero
  public :: nMixMax
  real :: PowerOfRIono

  real,public :: IonizationPotentialLowering_I(0:nZMax) = 0.0
  logical,public :: UseGroundStatWeight = .false.
contains

  subroutine get_excitation_levels_zero(rIonoSphereInv)
    real,intent(in) :: rIonoSphereInv

    integer :: iMix, iZ, iN, iL
    integer :: nGround, iCount
    real :: Gin
    real,dimension(0:nExcitation-1,nExcitation,0:nZMax,nMixMax) :: ExtraEnergy_IIII = 0.0
    !------------

    Partition_IIII(     :,:,:,:) = 0.0
    ExtraEnergyAv_II(       :,:) = 0.0
    VirialCoeffAv_II(       :,:) = 0.0
    SecondVirialCoeffAv_II( :,:) = 0.0
    Cov2ExtraEnergy_II(     :,:) = 0.0
    CovExtraEnergyVirial_II(:,:) = 0.0
    Cov2VirialCoeff_II(     :,:) = 0.0
    Partition_III(        :,:,:) = 0.0
    ExcitationEnergy_III( :,:,:) = 0.0


    PowerOfRIono = rIonoSphereInv**iPressureEffectIndex

    do iMix = 1, nMix
       LogGi_II(nZ_I(iMix):nZMax,iMix) = 0.0 
       Partition_III(0,nZ_I(iMix),iMix) = 1.0
       do iZ = 0, nZ_I(iMix)-1

          nGround = n_ground(iZ, nZ_I(iMix))

          Partition_IIII(0,nGround,iZ,iMix) = 1.0
          Partition_III(nGround,iZ,iMix) = 1.0
          LogGi_II(iZ,iMix) = 0.0
          
          PRINCIPAL:do iN = nGround, nExcitation
             Gin = 0.0
             iCount = 0
             do iL = 0, iN-1

                if(DoStateElimination.and.(ExcitationEnergy_IIII(iL,iN,iZ,iMix) + &
                     VirialCoeff4Energy_IIII(iL,iN,iZ,iMix) * PowerOfRIono >= &
                     IonizPotential_II(iZ+1,iMix) - IonizationPotentialLowering_I(iZ)))CYCLE

                ExtraEnergy_IIII(iL,iN,iZ,iMix) = ExcitationEnergy_IIII(iL,iN,iZ,iMix) + &
                     VirialCoeff4Energy_IIII(iL,iN,iZ,iMix) * PowerOfRIono

                ExcitationEnergy_III(iN,iZ,iMix) = ExcitationEnergy_III(iN,iZ,iMix) + &
                     Degeneracy_IIII(iL,iN,iZ,iMix) * ExtraEnergy_IIII(iL,iN,iZ,iMix)

                Gin = Gin + Degeneracy_IIII(iL,iN,iZ,iMix)
                iCount = iCount + 1
             end do

             if (Gin /= 0) &
                ExcitationEnergy_III(iN,iZ,iMix) = &
                  ExcitationEnergy_III(iN,iZ,iMix) / Gin
   
             if(iCount == 0)then
                !Eliminate the whole shell:
                nExcitation_II(iZ, iMix) = iN - 1
                EXIT PRINCIPAL
             end if
          end do PRINCIPAL
          if(UseGroundStatWeight .and. GIn > 0.5 .and. iN==nGround)&
               LogGi_II(iZ,iMix) = log(GIn)      
       end do
    end do

  end subroutine get_excitation_levels_zero

  !Fill in arrays ExtraEnergyAv_II and LogGi_II
  subroutine get_excitation_levels(iMix, iZMin, iZMax, nZ, TeInv, rIonoSphereInv)
    integer,intent(in) :: iMix
    integer,intent(in) :: iZMin, iZMax
    integer,intent(in) :: nZ
    real,intent(in) :: TeInv
    real,intent(in) :: rIonoSphereInv

    integer :: iZ, iN, iL  !Loop variables

    real    :: DeltaEnergy
    real    :: Gi, GiInv
    integer :: nGround, iCount

    real,parameter :: IndexPerThree    = iPressureEffectIndex/3.0
    real,parameter :: IndexSecondDeriv = IndexPerThree * (1.0 + IndexPerThree)

    !Energy correction due to the pressure ionization effect, averaged by
    !possible excitation levels [eV]
    real,dimension(0:nZMax,nMixMax) :: DeltaEnergyAv_II = 0.0

    !Energies of excited levels relative to the ground level, E_x = E_{exc} + \Delta E,
    !where E_{exc} stands for the excitation energy, and \Delta E is the correction
    !accounting for the pressure ionization
    real,dimension(0:nExcitation-1,nExcitation,0:nZMax,nMixMax) :: ExtraEnergy_IIII = 0.0

    ! -V \partial E_x/\partial V
    real,dimension(0:nExcitation-1,nExcitation,0:nZMax,nMixMax) :: VirialCoeff_IIII = 0.0
    !------------

    if (.not.UseExcitation) return
 !      call get_excitation_levels_zero(rIonoSphereInv)
 !      return
 !   end if

    Partition_IIII(     :,:,:,iMix) = 0.0
    DeltaEnergyAv_II(       :,iMix) = 0.0
    ExtraEnergyAv_II(       :,iMix) = 0.0
    VirialCoeffAv_II(       :,iMix) = 0.0
    SecondVirialCoeffAv_II(      :,iMix) = 0.0
    Cov2ExtraEnergy_II(     :,iMix) = 0.0
    CovExtraEnergyVirial_II(:,iMix) = 0.0
    Cov2VirialCoeff_II(     :,iMix) = 0.0
    Partition_III(        :,:,iMix) = 0.0
    ExcitationEnergy_III( :,:,iMix) = 0.0


    PowerOfRIono = rIonoSphereInv**iPressureEffectIndex
    LogGi_II(nZ_I(iMix):nZMax, iMix) = 0.0 
    Partition_III(0,nZ_I(iMix),iMix) = 1.0
    nExcitation_II(nZ_I(iMix),iMix) = -1

    nExcitation_II(:, iMix) = nExcitation

    do iZ = iZMin, iZMax

       nGround = n_ground(iZ, nZ)
       Gi = 0.0
       PRINCIPAL:do iN = nGround, nExcitation
          iCount = 0 
          do iL = 0, iN-1
             
             if(DoStateElimination.and.(ExcitationEnergy_IIII(iL,iN,iZ,iMix) + &
                  VirialCoeff4Energy_IIII(iL,iN,iZ,iMix) * PowerOfRIono >= &
                  IonizPotential_II(iZ+1,iMix) - IonizationPotentialLowering_I(iZ)))CYCLE

             DeltaEnergy = VirialCoeff4Energy_IIII(iL,iN,iZ,iMix) * PowerOfRIono
             ExtraEnergy_IIII(iL,iN,iZ,iMix) = ExcitationEnergy_IIII(iL,iN,iZ,iMix) + &
                  DeltaEnergy

             Partition_IIII(iL,iN,iZ,iMix) = Degeneracy_IIII(iL,iN,iZ,iMix)*&
                  exp(-ExtraEnergy_IIII(iL,iN,iZ,iMix) * TeInv)

             iCount = iCount +1

             Gi = Gi + Partition_IIII(iL,iN,iZ,iMix)

             DeltaEnergyAv_II(iZ,iMix) = DeltaEnergyAv_II(iZ,iMix) + &
                  Partition_IIII(iL,iN,iZ,iMix) * DeltaEnergy
             ExtraEnergyAv_II(iZ,iMix) = ExtraEnergyAv_II(iZ,iMix) + &
                  Partition_IIII(iL,iN,iZ,iMix) * ExtraEnergy_IIII(iL,iN,iZ,iMix)

             VirialCoeff_IIII(iL,iN,iZ,iMix) = IndexPerThree * DeltaEnergy

          end do
          if(iCount == 0)then
             !Eliminate the whole shell:
             nExcitation_II(iZ, iMix) = iN - 1
             EXIT PRINCIPAL
          end if
       end do PRINCIPAL


       if(Gi==0.0)then
          LogGi_II(iZ,iMix) = -30.0
          GiInv = 0.0
          CYCLE
       else
          LogGi_II(iZ,iMix) = log(Gi)
          GiInv = 1.0/Gi
       end if
       
       !Normalize to obtain the partial partition function and average energies
       Partition_IIII(:,nGround:nExcitation_II(iZ,iMix), iZ, iMix) = &
            Partition_IIII(:,nGround:nExcitation_II(iZ,iMix), iZ, iMix) * GiInv

       ExtraEnergyAv_II(iZ,iMix) = ExtraEnergyAv_II(iZ,iMix) * GiInv
       DeltaEnergyAv_II(iZ,iMix) = DeltaEnergyAv_II(iZ,iMix) * GiInv
       !Take averages over l
       
       do iN = nGround, nExcitation_II(iZ,iMix)
          Partition_III(iN,iZ,iMix) = sum( Partition_IIII(0:iN-1,iN,iZ,iMix)) 
          ExcitationEnergy_III( iN,iZ,iMix) = sum(&
               Partition_IIII(0:iN-1,iN,iZ,iMix) * ExtraEnergy_IIII(0:iN-1,iN,iZ,iMix))

          if (Partition_III(iN,iZ,iMix) /= 0.0) ExcitationEnergy_III(iN,iZ,iMix) = &
               ExcitationEnergy_III(iN,iZ,iMix) / Partition_III(iN,iZ,iMix)
       end do

       VirialCoeffAv_II (iZ,iMix) = IndexPerThree    * DeltaEnergyAv_II(iZ,iMix)
       SecondVirialCoeffAv_II(iZ,iMix) = IndexSecondDeriv * DeltaEnergyAv_II(iZ,iMix)



       do iN = nGround, nExcitation_II(iZ,iMix)
          do iL = 0, iN-1
             if (Partition_IIII(iL,iN,iZ,iMix) == 0.0)CYCLE

             Cov2ExtraEnergy_II(iZ,iMix) = Cov2ExtraEnergy_II(iZ,iMix) + &
                  Partition_IIII(iL,iN,iZ,iMix) * &
                  (ExtraEnergy_IIII(iL,iN,iZ,iMix) - ExtraEnergyAv_II(iZ,iMix))**2

             CovExtraEnergyVirial_II(iZ,iMix) = CovExtraEnergyVirial_II(iZ,iMix) + &
                  Partition_IIII(iL,iN,iZ,iMix) * &
                  (ExtraEnergy_IIII(iL,iN,iZ,iMix) - ExtraEnergyAv_II(iZ,iMix)) * &
                  VirialCoeff_IIII(iL,iN,iZ,iMix)

             Cov2VirialCoeff_II(iZ,iMix) = Cov2VirialCoeff_II(iZ,iMix) + &
                  Partition_IIII(iL,iN,iZ,iMix) * &
                  (VirialCoeff_IIII(iL,iN,iZ,iMix) - VirialCoeffAv_II(iZ,iMix))**2

          end do
       end do

    end do
  end subroutine get_excitation_levels

end module CRASH_ModExcitation
