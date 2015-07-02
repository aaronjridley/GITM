!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=======================
module CRASH_ModThermoDynFunc
  use CRASH_ModPartition
  use ModConst
  implicit none
  save

  real,parameter:: cZMin = cTiny

contains
  !=======================================  

  ! Calculate the pressure in the plasma [Pa]
  ! Can only be called after set_ionization_equilibrium has executed

  real function pressure()
    use CRASH_ModFermiGas,ONLY:RPlus

    pressure = Na * cEV * Te * (1.0 + zAv*RPlus + VirialCoeffAv)

  end function pressure

  !============================================
  ! calculate the average internal energy per atomic unit [eV]
  ! Can only be called after set_ionization_equilibrium has executed 

  real function internal_energy()
    use CRASH_ModFermiGas,ONLY:RPlus

    internal_energy = 1.50 * Te *( 1.0 + zAv * RPlus) + EAv

  end function internal_energy

  !=======================================!
  ! Calculate the electron pressure in the plasma [Pa]
  ! Can only be called after set_ionization_equilibrium has executed

  real function pressure_e()
    use CRASH_ModFermiGas,ONLY:RPlus

    pressure_e = Na * cEV * Te * (zAv*RPlus + VirialCoeffAv)

  end function pressure_e

  !============================================
  ! calculate the average internal energy of electrons per atomic unit [eV]
  ! Can only be called after set_ionization_equilibrium has executed 

  real function internal_energy_e()
    use CRASH_ModFermiGas,ONLY:RPlus

    internal_energy_e = 1.50 * Te *zAv * RPlus + EAv

  end function internal_energy_e

  !=======================================!
  ! Calculating the Z average values from populations

  real function z_averaged()
    z_averaged = zAv
  end function z_averaged
  !======================
  ! Calculating the Z^2 average values from populations
  real function z2_averaged()
    !The average of square is the averaged squared plus dispersion squared
    z2_averaged = Z2
  end function z2_averaged
end module CRASH_ModThermoDynFunc
!=======================!
module CRASH_ModThermoDynDeriv
  use CRASH_ModPartition
  use CRASH_ModThermoDynFunc
  implicit none
  SAVE
Contains
  !==================================
  !Calculate the specific heat capacity at constant volume 
  !(derivative of internal energy wrt Te) from temperature:
  !Can only be called after set_ionization_equilibrium has executed
  real function heat_capacity()
    use CRASH_ModFermiGas
    !------------------

    if( zAv <= cZMin )then
       heat_capacity = 1.50; return
    end if

    ! calculate the heat capacity:
    heat_capacity = 1.50 * (1.0 +2.5*RPlus*zAv) + DeltaETeInv2Av &
         -( 1.5 * zAv - ETeInvDeltaZAv )**2 / (zAv*RMinus + DeltaZ2Av)

  end function heat_capacity 
 !==================================
  !Calculate the specific heat capacity of electrons at constant volume 
  !(derivative of internal energy wrt Te) from temperature:
  !Can only be called after set_ionization_equilibrium has executed
  real function heat_capacity_e()
    use CRASH_ModFermiGas
    !------------------

    ! calculate the heat capacity:
    heat_capacity_e = heat_capacity() - 1.50 

  end function heat_capacity_e 
  
  !==========================================================================!
  !Thermodynamic derivatives for pressure
  !Thermodynamic derivatives: use abbreviations:
  !Ov means over 
  !AtCrho means at constant  rho
  !=========================================================================!

  ! Calculate (1/Na)*(\partial P/\partial T)_{\rho=const}
  ! i.e the temperuture derivative of the specific pressure (P/Na) 
  real function d_pressure_over_d_te()
    use CRASH_ModFermiGas,ONLY:RPlus,RMinus
    !----------------------------------------------------------
    if(zAv<=cZMin)then
       d_pressure_over_d_te = 1.0
       return
    end if
    !\
    !Calculate derivatives at const Na:
    !/

    !For specific pressure:
    d_pressure_over_d_te = 1.0 + 2.5*zAv*RPlus - &
         (zAv + CovVirialZ) * ( 1.5 * zAv - ETeInvDeltaZAv ) / (zAv*RMinus + DeltaZ2Av) + &
         CovEnergyVirial

  end function d_pressure_over_d_te

  ! Calculate (1/Na)*(\partial P/\partial T)_{\rho=const}
  ! i.e the temperuture derivative of the specific pressure (P/Na)
  real function d_pressure_e_over_d_te()
    use CRASH_ModFermiGas,ONLY:RPlus,RMinus
    !----------------------------------------------------------
 
    !\
    !Calculate derivatives at const Na:
    !/

    !For specific pressure:
    d_pressure_e_over_d_te = d_pressure_over_d_te() - 1.0

  end function d_pressure_e_over_d_te
  !===============================================================
  ! Calculate (1/Na)*(\partial P/\partial T)_{\rho=const}
  ! i.e the temperuture derivative of the specific pressure (P/Na) 
  real function d_pressure_over_d_tz()
    use CRASH_ModFermiGas,ONLY:RPlus,RMinus
    !----------------------------------------------------------
    if(zAv<=cZMin)then
       d_pressure_over_d_tz = 1.0
       return
    end if
    !\
    !Calculate derivatives at const Na:
    !/

    !For specific pressure:
    d_pressure_over_d_tz = 2.5*zAv*RPlus - 1.5 * zAv / RMinus

  end function d_pressure_over_d_tz
  !========================================================================!
  !Calculate !(\rho/P)*(\partial P/\Partial \rho)_{T=const},
  !i.e the compressibility at constant temperature
  real function compressibility_at_const_te()
    use CRASH_ModFermiGas,ONLY:RPlus, RMinus


    !-------------------
    if(zAv<=cZMin)then
       compressibility_at_const_te = 1.0
       return
    end if
    !\
    !Derivatives at constant T:
    !/

    !For specific pressure:
    !(Na /P) \left(\partial P / \partial Na \right)_{T=const}

    compressibility_at_const_te = &
         (1.0 + (zAv + CovVirialZ)**2 / (zAv*RMinus + DeltaZ2Av) - Cov2Virial + SecondVirialAv) /&
         (1.0 + zAv*RPlus + VirialCoeffAv)

  end function compressibility_at_const_te

  !========================================================================!
  !Calculate !(\rho/P)*(\partial P/\Partial \rho)_{T=const},
  !i.e the compressibility at constant temperature
  real function compressibility_at_const_te_e()
    use CRASH_ModFermiGas,ONLY:RPlus, RMinus

    real :: Compr
    !-------------------
    if(zAv<=cZMin)then
       compressibility_at_const_te_e = 0.0
       return
    end if
    !\
    !Derivatives at constant T:
    !/

    !For specific pressure:
    !(Na /P) \left(\partial P / \partial Na \right)_{T=const}

    Compr = compressibility_at_const_te()
    compressibility_at_const_te_e = Compr + &
         (Compr - 1.0) / (zAv*RPlus + VirialCoeffAv)

  end function compressibility_at_const_te_e

  !===================================================================
  subroutine get_gamma(GammaOut,GammaSOut,GammaMaxOut,GammaeOut,GammaSeOut)
    real,optional,intent(out)::GammaOut    !1+P/UInternal
    real,optional,intent(out)::GammaSOut   !The speed of sound squared*Rho/P
    real,optional,intent(out)::GammaMaxOut !max(Gamma,GammaS)
    real,optional,intent(out)::GammaeOut   !1+P_e/UInternal_e
    real,optional,intent(out)::GammaSeOut  !The speed of sound in the electron gas squared*Rho/P
    real :: Gamma, GammaS, Gammae, GammaSe

    real,parameter::Gamma0=5.0/3.0
    !--------------------------------------!

    if(zAv<=cZMin)then
       if(present(GammaOut))   GammaOut=Gamma0
       if(present(GammaSOut))  GammaSOut=Gamma0
       if(present(GammaMaxOut))GammaMaxOut=Gamma0
       if(present(GammaeOut))   GammaeOut = Gamma0
       if(present(GammaSeOut)) GammaSeOut=Gamma0
       return
    end if
    Gamma  = 1.0 + pressure  ()/( Na * cEv * internal_energy  ())
    Gammae = 1.0 + pressure_e()/( Na * cEv * internal_energy_e())

    if(present(GammaOut ))GammaOut  = Gamma
    if(present(GammaeOut))GammaeOut = Gammae

    !Define GammaS in terms the speed of sound: GammaS*P/\rho=(\partial P/\partial \rho)_{s=const}
    !Use a formula 3.72 from R.P.Drake,{\it High-Energy Density Physics}, Springer, Berlin-Heidelberg-NY, 2006
    ! and apply the thermodinamic identity: 
    !(\partial \epsilon/\partial \rho)_{T=const}-P/\rho^2)=-T(\partial P/\partial T)_{\rho=const)/\rho^2

    GammaS =  compressibility_at_const_te  () + &
         d_pressure_over_d_te  ()**2 * (Na * Te * cEV) /(heat_capacity  () * pressure  ())

    if(present(GammaSOut ))GammaSOut  = GammaS
    if(present(GammaSeOut))then
       ! Note heat_capacity_e() is zero if zAv <= cZMin
       GammaSe = compressibility_at_const_te_e() + &
            d_pressure_e_over_d_te()**2 * (Na * Te * cEV) /(heat_capacity_e() * pressure_e())
       GammaSeOut = GammaSe
    end if

    if(present(GammaMaxOut))GammaMaxOut = max( GammaS,  Gamma)

  end subroutine get_gamma
end module CRASH_ModThermoDynDeriv
!=======================!
module CRASH_ModStatSum
  use CRASH_ModPartition
  use CRASH_ModThermoDynDeriv
  implicit none
  SAVE

  integer :: iIterTe !Temperature iterations counter
  logical :: UsePreviousTe = .false.

  !\
  ! WARNING: it is not recommended to change the tolerance parameters
  ! The algorithm convergence is not guaranteed with the reduced tolerance(s). 
  !/ 

  ! Accuracy of internal energy needed [(% deviation)/100]

  real, public:: Tolerance = 0.001

Contains
  !==========================================================================
  !Set the element and its Ionization Potentials
  subroutine set_element( nZIn)
    integer,intent(in) :: nZIn
    integer:: nZIn_I(1)
    real, parameter:: One_I(1) = (/ 1.0 /)
    !-----------------------------------------------------------------------
    nZIn_I(1) = nZIn
    call set_mixture(nMixIn=1, nZIn_I=nZIn_I, ConcentrationIn_I=One_I)
  end subroutine set_element
  !==========================================================================
  subroutine set_temperature(Uin, NaIn, iError)

    real,intent(in) :: Uin,& !Average internal energy per atomic unit [eV]
                       NaIn  !Density of heavy particles [# of particles/m^3]
    integer,intent(out),optional::iError

    real :: UDeviation,& ! The difference between the given internal energy and the calculated one
            ToleranceUeV ! The required accuracy of U in eV
    !-------------------------
    call init_madelung(NaIn)

    if( .not.UsePreviousTe ) then
       if(nMix>1) then 
          !To keep the backward compatibility, because the choice for mixture was different
          Te = UIn / 1.5
       else
          Te = max(UIn,IonizPotential_II(1,1)) * 0.1
       end if
    end if
    if( UIn <= 0.03 * minval(IonizPotential_II(1,1:nMix))) then

       Te=UIn/1.50 ;  call set_zero_ionization
       call check_applicability(iError)
       return
    end if

    ToleranceUeV = Tolerance * Uin
    iIterTe = 0

    !Use Newton-Rapson iterations to get a better approximation of Te:
    UDeviation = 2.*ToleranceUeV
    iterations: do while(abs(UDeviation) >ToleranceUeV .and. iIterTe<=20)

       !Find the populations for the trial Te

       call set_ionization_equilibrium(Te, Na) 
       UDeviation = internal_energy() - Uin

       !Calculate the improved value of Te, limiting the iterations so 
       !they can't jump too far out

       Te = min(2.0 * Te, max(0.5 * Te, Te - UDeviation/heat_capacity()))  

       iIterTe = iIterTe+1
    end do iterations

    call check_applicability(iError)
  end subroutine set_temperature
  !==========================================================================
  subroutine u_e_to_temperature(Uin, NaIn, iError)

    real,intent(in) :: Uin,& !Average internal energy per atomic unit [eV]
                       NaIn  !Density of heavy particles [# of particles/m^3]
    integer,intent(out),optional::iError

    real :: UDeviation,& ! The difference between the given internal energy and the calculated one
            ToleranceUeV ! The required accuracy of U in eV
    !-------------------------
    if (Uin < 0.0)&
         call CON_stop('Uin < 0.0 in u_e_to_temperature')

    call init_madelung(NaIn)

    if( .not.UsePreviousTe ) then
       
          Te = max(UIn,minval(IonizPotential_II(1,1:nMix))) * 0.1
   
    end if

    ToleranceUeV = Tolerance * Uin
    iIterTe = 0

    !Use Newton-Rapson iterations to get a better approximation of Te:
    UDeviation = 2.*ToleranceUeV
    iterations: do while(abs(UDeviation) >ToleranceUeV .and. iIterTe<=20)

       !Find the populations for the trial Te

       call set_ionization_equilibrium(Te, Na)
       UDeviation = internal_energy_e() - Uin

       !Calculate the improved value of Te, limiting the iterations so
       !they can't jump too far out

       Te = min(2.0 * Te, max(0.5 * Te, Te - UDeviation/heat_capacity_e()))

       iIterTe = iIterTe+1
    end do iterations

    call check_applicability(iError)
  end subroutine u_e_to_temperature
  !==========================================================================
  subroutine pressure_to_temperature(PToNaRatio, NaIn, iError)
    real,intent(in) :: PToNaRatio,& !Presure divided by Na [eV]
         NaIn !Density of heavy particles [# of particles/m^3]
    integer,intent(out),optional::iError

    ! The difference between the given pressure and the calculated one
    real :: PDeviation,& 
         TolerancePeV ! The required accuracy of P in eV*Na

    !Thermodinamical derivative (\partial P/\partial T)
    !at constant Na, divided by Na
    real :: NaInvDPOvDTAtCNa  
    !------------------------------------

    call init_madelung(NaIn)

    !It is difficult to make an initial guess about temperature
    !The suggested estimate optimizes the number of idle iterations:
    if(.not.UsePreviousTe) & 
         Te = max(1.50* PToNaRatio, IonizPotential_II(1,1) ) * 0.1

    if(PToNaRatio<=0.03*IonizPotential_II(1,1))then
       Te = PToNaRatio; call set_zero_ionization
       call check_applicability(iError)
       return
    end if

    TolerancePeV = Tolerance * PToNaRatio
    iIterTe = 0
    PDeviation = 2.*TolerancePeV

    !Use Newton-Rapson iterations to get a better approximation of Te:

    iterations: do while(abs(PDeviation) >TolerancePeV .and. iIterTe<=20)

       !Find the populations for the trial Te
       call set_ionization_equilibrium(Te, Na) 
       PDeviation = pressure()/(Na*cEV) - PToNaRatio

       !Calculate the improved value of Te, limiting the iterations so 
       !they can't jump too far out, if the initial guess for Te is bad.
       Te = min(2.0*Te, max(0.5*Te, Te - PDeviation/d_pressure_over_d_te() ) )  

       iIterTe = iIterTe+1  !Global variable, which is accessible from outside

    end do iterations
    call check_applicability(iError)
  end subroutine pressure_to_temperature
  !==========================================================================
  subroutine pressure_e_to_temperature(PeToNaRatio, NaIn, iError)
    real,intent(in) :: PeToNaRatio,& !Pressure divided by Na [eV]
         NaIn !Density of heavy particles [# of particles/m^3]
    integer,intent(out),optional::iError

    ! The difference between the given pressure and the calculated one
    real :: PDeviation,&
         TolerancePeV ! The required accuracy of P in eV*Na

    !Thermodinamical derivative (\partial P/\partial T)
    !at constant Na, divided by Na
    real :: NaInvDPOvDTAtCNa
    !------------------------------------
    if (PeToNaRatio < 0.0)&
         call CON_stop('PeToNaRatio < 0.0 in u_e_to_temperature')

    call init_madelung(NaIn)

    !It is difficult to make an initial guess about temperature
    !The suggested estimate optimizes the number of idle iterations:
    if(.not.UsePreviousTe) &
         Te = max(1.50* PeToNaRatio, IonizPotential_II(1,1) ) * 0.1

    TolerancePeV = Tolerance * PeToNaRatio
    iIterTe = 0
    PDeviation = 2.*TolerancePeV

    !Use Newton-Rapson iterations to get a better approximation of Te:

    iterations: do while(abs(PDeviation) >TolerancePeV .and. iIterTe<=20)

       !Find the populations for the trial Te
       call set_ionization_equilibrium(Te, Na)
       PDeviation = pressure_e()/(Na*cEV) - PeToNaRatio

       !Calculate the improved value of Te, limiting the iterations so
       !they can't jump too far out, if the initial guess for Te is bad.
       Te = min(2.0*Te, max(0.5*Te, Te - PDeviation/d_pressure_e_over_d_te() ) )

       iIterTe = iIterTe+1  !Global variable, which is accessible from outside

    end do iterations
    call check_applicability(iError)
  end subroutine pressure_e_to_temperature

end module CRASH_ModStatSum
