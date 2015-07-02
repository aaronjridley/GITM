!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module CRASH_ModFermiGas
  use ModConst,ONLY : cPi

  implicit none
  SAVE
  PRIVATE !Except

  !The module does nothing if this logical is false  
  logical,public :: UseFermiGas = .true.
  
  ! At LogGe >= LogGeMinBoltzmann the electrons are treated as a Boltzmann gas
  real, public :: LogGeMinBoltzmann = 4.0

  ! At LogGe >= LogGeMinFermi the effects of Fermi statistics are accounted for
  real, public :: LogGeMinFermi = -4.0

  ! If LogGeMinFermi < -3.0,
  ! init_fermi_function automatically resets it to -3.0 and
  ! sets UseAsymptLarge to .true., so that:
  !    At LogGe < -3.0 we use the first 3 sum terms of the asymptotic series
  !    At -3.0 <= LogGe <= 0.0 we use Taylor series
  !
  ! The value of -3.0 has been chosen because the Taylor series converges at |log(g_e)| < \pi
  logical :: UseAsymptLargeNeg = .false.

  ! At LogGeMinFermi <= LogGe <= 0.0 we use Taylor series to calculate Fermi Functions.
  !
  !
  !  At LogGe >= max(0.0,LogGeMinFermi) we use the series 
  !  \sum_{n=1}^\infty{\frac{(-1)^{n+1}}{n^{\nu+1 (g_e)^{n}}}}
  !

  ! Logarithm of the electron statistical weight = exp(-\mu/T)
  real, public :: LogGe

  !Correcting coefficients in thermodynamic functions  
  real, public :: rMinus = 1.0, rPlus = 1.0

  ! Named indices.
  integer, parameter :: NuEqMinus12_ = 1, NuEq12_ = 2, NuEq32_ = 3
  real, parameter :: NuPlus1_I(3) = (/ 0.5, 1.50, 2.50 /)



  !nStepSeries is a number of values between
  !LogGeMinSeries and LogGeMinBoltzmann, for which
  !FermiFunctionTableSeries_II is to be calculated.
  integer, parameter :: nStepSeries = 20
  real :: FermiFunctionTableSeries_II(0:nStepSeries, NuEqMinus12_ : NuEq32_) = 1.0

  integer, parameter :: nStepTaylorMax = 15
  integer :: nStepTaylor = nStepTaylorMax

  !The table for positive indices are only used in the test
  real :: FermiFunctionTableTaylor_II(&
       -nStepTaylorMax:nStepTaylorMax, NuEqMinus12_ : NuEq32_) = 1.0

  !Taylor coefficients = \frac{(1-2^{i-\nu}) \zeta(\nu+1-i)}{i!}
  !TaylorSeriesCoeff_II are Taylor coefficients multiplied by \pi^i
  !to maintain these numbers around 1
  integer, parameter :: nTaylor = 100
  real, dimension(0:nTaylor, NuEqMinus12_ : NuEq32_) :: TaylorSeriesCoeff_II = 1.0



  !The things related to the calculation using Tailor series.
  !The series over Semi-Integer arguments (....S_I arrays)
  !For i=1,2,3, etc in all arrays below the variable s passes 1/2,3/2,5/2 etc

  !2^{1-s)
  real, dimension(nTaylor+1):: TwoPoweredOneMinusS_I = 1.0

  !\sum_{i=1}^\infty{\frac{(-1)^{i+1}}{i^s}}
  real, dimension(nTaylor+1):: ZetaSeriesS_I = 1.0

  !Zeta function = ZetaSeries/(1-TwoPoweredOneMinusS_I)
  real, dimension(nTaylor+1):: ZetaFunctionS_I = 1.0


  public :: init_fermi_function, iterate_ge
  
  public :: test_fermi_function 

contains
  !=====================================
  !The Fermi functions are defined to be
  !Fe_{\nu}(g_e)=(1/\Gamma(\nu+1))\int_0^\infty{x^\nu dx/(g_e \exp(x)+1)}
  !At abs(g_e)>=1 the integral can be taken by developing into a convergent 
  !series in powers of 1/g_e. At g_e=1 the integral may be also expressed 
  !in terms of the Riemann zeta-function, which is used in testing
  !==============================
  !Initialaze the calculation of Fermi functions for \nu+1=1/2, 3/2, 5/2
  subroutine init_fermi_function
    if (LogGeMinFermi < -3.0) then
       UseAsymptLargeNeg = .true.
       LogGeMinFermi = -3.0
    end if

    call fill_lookup_table_series
    if (LogGeMinFermi < 0.0) then
       call init_taylor_series
       call fill_lookup_table_taylor
       FermiFunctionTableSeries_II(0,:) =  FermiFunctionTableTaylor_II(0,:)
    end if

  contains
    !========================================================
    subroutine init_taylor_series
      integer :: iNu, iTaylor
      
      integer :: i, n
      real,dimension(1000,nTaylor+1)::InvIPoweredS_II
      real :: SignOfTerm
      real :: CoreMultiplier
      real :: TwoPoweredIMinusNu
      !-----------------------------------
      InvIPoweredS_II = 0.0
      !Calculate TwoPoweredOneMinusS_I
      TwoPoweredOneMinusS_I(1)=sqrt(2.0)
      
      do i=2,nTaylor+1
         TwoPoweredOneMinusS_I(i) = TwoPoweredOneMinusS_I(i-1) * 0.5
      end do
      !\
      !Calculate series need for implement in Riemann's Zeta function
      
      InvIPoweredS_II(1,:) = 1.0
      SignOfTerm = 1.0
      do i=2,1000
         SignOfTerm = -SignOfTerm
         InvIPoweredS_II(i,1) = SignOfTerm/sqrt(real(i))
      end do
      
      do n = 2,nTaylor+1
         internal:do i = 2,1000
            InvIPoweredS_II(i,n) = InvIPoweredS_II(i,n-1)/real(i)
            if(abs(InvIPoweredS_II(i,n)) < 3.0e-8)exit internal
         end do internal
      end do
      do n=2,nTaylor+1
         ZetaSeriesS_I(n) = sum(InvIPoweredS_II(:,n))
      end do
      
      !Use explicit value for zeta(1/2) as the series converges too slowly
      ZetaSeriesS_I(1) = (-1.46035) * (1.0 - 2.0*sqrt(0.5))
      ZetaSeriesS_I(2) = (2.61237535) * (1.0 - sqrt(0.5))
      
      !Calculate Riemann's Zeta function
      ZetaFunctionS_I = ZetaSeriesS_I/(1.0 - TwoPoweredOneMinusS_I)
      !End of zeta function calculations
      !/
      
      !\
      !Calculate Taylor series coefficients
      do iNu = NuEqMinus12_ , NuEq32_
         TwoPoweredIMinusNu = 2 ** (1 - NuPlus1_I(iNu))
         
         !Add up sum terms with zeta function of a positive argument
         CoreMultiplier = 1.0
         !(iNu - 1) is the index of the last Taylor coefficient that doesn't
         !use zeta function of a negative argument
         do iTaylor = 0, min(nTaylor, iNu - 1)
            if (iTaylor /= 0)& 
                 CoreMultiplier = CoreMultiplier / iTaylor * cPi
            
            TaylorSeriesCoeff_II(iTaylor, iNu) = &
                 CoreMultiplier * (1 - TwoPoweredIMinusNu) * &
                 ZetaFunctionS_I(iNu - iTaylor)
            
            TwoPoweredIMinusNu = TwoPoweredIMinusNu * 2.0
         end do
         
         !Add up sum terms with zeta function of a negative argument using
         !this reflection formula:
         ! \zeta(x) = 2^x \pi^{x-1} \sin(\pi x / 2) \Gamma(1-x) \zeta(1-x)
         CoreMultiplier = cPi ** (iNu - 1) / 4.0
         do i = 2, iNu            !divide by (nPositiveZeta + 1)!
            CoreMultiplier = CoreMultiplier / i
         end do
         
         do iTaylor = iNu, nTaylor
            
            TaylorSeriesCoeff_II(iTaylor, iNu) = &
                 CoreMultiplier * (1 - TwoPoweredIMinusNu) * &
                 ZetaFunctionS_I(iTaylor - iNu + 2)

            if (mod(iTaylor - iNu, 4) < 2) then
               !set sign of \sin \frac{\pi (\nu + 1 - i)}{2}
               TaylorSeriesCoeff_II(iTaylor, iNu) = &
                    -TaylorSeriesCoeff_II(iTaylor, iNu)  
            end if

            CoreMultiplier = CoreMultiplier * (iTaylor - NuPlus1_I(iNu) + 1) / &
                 (2.0 * (iTaylor + 1))
            TwoPoweredIMinusNu = TwoPoweredIMinusNu * 2.0
         end do
      end do
    end subroutine init_taylor_series
  !========================================================
    subroutine fill_lookup_table_series
      integer :: iStep, iNu, iTaylor
      real :: GeInv, GeInvN, SumTerm, RealN
      
      real :: CoreMultiplier
      real :: X
      !-------------------------------------------
      
      !Fill in the lookup table
      do iStep=0,nStepSeries-1
         !We calculate \sum_{n=1}^\infty (-1/Ge)**(n-1) / n**(\mu+1)
         !THE RESULTING FERMI FUNCTIONS ARE ALL MULTIPLED BY Ge!!
         GeInv = (-1.0)*exp(-(max(LogGeMinFermi, 0.0)*(nStepSeries-iStep)&
              +iStep*LogGeMinBoltzmann)/nStepSeries)
         do iNu = NuEqMinus12_ ,  NuEq32_
            GeInvN = 1.0; SumTerm = 1.0 ; RealN = 1.0
            FermiFunctionTableSeries_II(iStep, iNu) = 1.0
            do while(abs(SumTerm) > 3.0e-4)
               RealN = RealN + 1.0
               GeInvN = GeInvN * GeInv
               SumTerm = GeInvN/RealN** NuPlus1_I(iNu)
               FermiFunctionTableSeries_II(iStep, iNu) = &
                    FermiFunctionTableSeries_II(iStep, iNu) + SumTerm
            end do
         end do
      end do
      
      !To achieve a continuity at G_e_Min_Boltzmann, put the functions  to be 1
      FermiFunctionTableSeries_II(nStepSeries, :) = 1.0
    end subroutine fill_lookup_table_series
    !========================================================
    subroutine fill_lookup_table_taylor
      real :: LogGeMin
      real :: LogGeMax
      
      integer :: iStep, iNu, iTaylor
      real :: GeInv, GeInvN, SumTerm, RealN
      
      real :: CoreMultiplier
      real :: X
      !-------------------------------------------
      LogGeMin = LogGeMinFermi
      LogGeMax = - LogGeMin
      !Fill in the lookup table
      do iStep = -nStepTaylor, nStepTaylor
         X = -(LogGeMin*(nStepTaylor-iStep)+&
              (iStep+nStepTaylor)*LogGeMax)/(2*nStepTaylor)
         
         !Calculate Taylor series
         do iNu = NuEqMinus12_, NuEq32_
            SumTerm = 100
            FermiFunctionTableTaylor_II(iStep,iNu) = 0.0
            CoreMultiplier = 1.0
            
            do iTaylor = 0, nTaylor
               if (abs(SumTerm) < 3.0e-8) exit
               
               SumTerm = TaylorSeriesCoeff_II(iTaylor,iNu) * CoreMultiplier
               FermiFunctionTableTaylor_II(iStep,iNu) = &
                    FermiFunctionTableTaylor_II(iStep,iNu) + SumTerm

               CoreMultiplier = CoreMultiplier / cPi * X  !CoreMultiplier = (x/\pi)^i
            end do
            
            !THE RESULTING FERMI FUNCTIONS ARE ALL MULTIPLED BY Ge!!
            FermiFunctionTableTaylor_II(iStep, iNu) = &
                 FermiFunctionTableTaylor_II(iStep, iNu) * exp(-X)
         end do
      end do
    end subroutine fill_lookup_table_taylor
  end subroutine init_fermi_function
  !========================================================
  subroutine iterate_ge(zAvr,Delta2I,LogGe1,Diff)
    real, intent(in) :: zAvr    !<i>
    real, intent(in) :: Delta2I !<i^2>-<i>^2
    real, intent(in) :: LogGe1  !log(G_{e1}

    real, intent(out):: Diff    !<i>-G_{e1}Fe_{1/2}(g_e)

    integer :: iStep, iStep1
    real :: Residual
    real :: FermiFunction_I(NuEqMinus12_:NuEq32_)
    !------------------------

    if (UseAsymptLargeNeg .and. LogGe < LogGeMinFermi) then
       FermiFunction_I = large_negative()
    else if (LogGe <= 0.0 .and. LogGeMinFermi < 0.0) then
       Residual = (LogGe - LogGeMinFermi)/&
            ( - LogGeMinFermi)*nStepTaylor

       Residual = max(Residual, 0.0)
       iStep    = int(Residual) - nStepTaylor
       iStep1   = min(iStep+1,0)

       Residual = Residual - real(iStep) - nStepTaylor
       FermiFunction_I = FermiFunctionTableTaylor_II(iStep,:)*(1.0-Residual)+&
            FermiFunctionTableTaylor_II(iStep1,:)*Residual
    else 
       Residual = (LogGe - max(LogGeMinFermi, 0.0))/&
            (LogGeMinBoltzmann - max(LogGeMinFermi, 0.0))*nStepSeries

       Residual = min(max(Residual,0.0),real(nStepSeries))
       iStep    = int(Residual)
       iStep1   = min(iStep+1,nStepSeries)
   
       Residual = Residual - real(iStep)
       FermiFunction_I = FermiFunctionTableSeries_II(iStep,:)*(1.0-Residual)+&
            FermiFunctionTableSeries_II(iStep1,:)*Residual
    end if
    


    rMinus = FermiFunction_I(NuEqMinus12_)/FermiFunction_I(NuEq12_)
    rPlus  = FermiFunction_I(NuEq32_) / FermiFunction_I(NuEq12_)

    Diff = (zAvr - exp(LogGe1-LogGe) * FermiFunction_I(NuEq12_))/&
         (Delta2I/zAvr + rMinus)
    LogGe = LogGe - Diff/zAvr
    !To improve the convergence at low zAvr uncomment this line
    ! Diff = Diff/min(zAvr,1.0)

  end subroutine iterate_ge
  !========================================================
  subroutine test_fermi_function
    integer :: iStep
    real :: FermiFunction_I(NuEqMinus12_:NuEq32_)
    !---------------

    write(*,*)'zeta-function(1/2) zeta-function(3/2) zeta-function(5/2):'
    write(*,*)'reference data: -1.46035 2.612 1.341'
    nStepTaylor = 13
    LogGeMinFermi = -1.3*log(10.)
    LogGeMinBoltzmann = 2.0*log(10.)
    call init_fermi_function

    write(*,*)'calculated values: ', ZetaFunctionS_I(1:3)
    write(*,*)'zeta-function multiplied by (1-2*0.5^\nu) = ',&
         -1.46035*(1.0 -2.0*sqrt(0.5)),&
         2.612*(1.0-sqrt(0.5)),1.341*(1.0-0.5*sqrt(0.5))
    write(*,*)'Calculated ZetaSeriesS_I: ',ZetaSeriesS_I(1:3)
    write(*,*)'Calculated TaylorSeriesCoeff_II for \nu = 1/2:', &
         TaylorSeriesCoeff_II(0:2, 2)
    write(*,*)'The computed values at g_e=1 are:', &
         FermiFunctionTableSeries_II(0 , :)


    write(*,*)'Standard Fermi function table'
    write(*,*)'At the boudary points the values are calculated twice,' 
    write(*,*)'using the different asymptotics.'    

    write(*,*)'log10g_e  g_e*Fe1/2   rMinus rPlus'

    do iStep = -20, -13
       LogGe = 0.1*iStep*log(10.0)
       FermiFunction_I = large_negative()
       write(*,*) LogGe/log(10.0), &
            FermiFunction_I(2), &
            FermiFunction_I(1)/FermiFunction_I(2), &
            FermiFunction_I(3)/FermiFunction_I(2)
    end do
    write(*,*) &
         'In the code the above line would be calculated via the Taylor series'
    do iStep = -nStepTaylor,0
       write(*,*) (LogGeMinFermi*(nStepTaylor-iStep)-&
            (iStep+nStepTaylor)*LogGeMinFermi)/(2*nStepTaylor)/&
            log(10.0), FermiFunctionTableTaylor_II(iStep,2),&
            FermiFunctionTableTaylor_II(iStep,1)/FermiFunctionTableTaylor_II(iStep,2),&
            FermiFunctionTableTaylor_II(iStep,3)/FermiFunctionTableTaylor_II(iStep,2)
    end do
   

    LogGeMinFermi = 0.0
    LogGeMinBoltzmann = 2.0*log(10.)
    call init_fermi_function

    write(*,*) 'Series: for LogGeMinFermi<0 the next line will be overwritten'
    do iStep = 0, nStepSeries
       write(*,*) (0.0 * (nStepSeries-iStep)+iStep*LogGeMinBoltzmann)/nStepSeries/&
            log(10.0), FermiFunctionTableSeries_II(iStep,2),&
            FermiFunctionTableSeries_II(iStep,1)/FermiFunctionTableSeries_II(iStep,2),&
            FermiFunctionTableSeries_II(iStep,3)/FermiFunctionTableSeries_II(iStep,2)
       if(iStep==13)write(*,*)'--------------------------------------------'
    end do

    LogGeMinFermi = -1.3*log(10.)
    call init_fermi_function

    !Test Taylor series for 0.0 <= \log(g_e) <= 3.0
   
    write(*,*)'Fermi function table for \log(g_e) >= 0 calculated via Taylor series:'
      
    do iStep = 0, nStepTaylor
       write(*,*) (0.0*(nStepTaylor-iStep)+iStep*abs(LogGeMinFermi))/nStepTaylor/&
            log(10.0), FermiFunctionTableTaylor_II(iStep,2),&
            FermiFunctionTableTaylor_II(iStep,1)/FermiFunctionTableTaylor_II(iStep,2),&
            FermiFunctionTableTaylor_II(iStep,3)/FermiFunctionTableTaylor_II(iStep,2)
    end do
 

  end subroutine test_fermi_function
  !========================================================
  function large_negative()
    real :: large_negative(NuEqMinus12_ : NuEq32_)
    integer :: iNu
    real :: X, X2Inv
    real, parameter :: ZetaFunction2 = 1.6449
    real, parameter :: TwoZetaFunctionSeries4 = 1.8941  ! 2*(1-2^{1-4})*\zeta(4)
    !-----------------
    X = -LogGe
    !Fe_{-1/2}(-log(Ge)\sim \sqrt{-Log(Ge)}/\Gamma(3/2)
    large_negative(NuEqMinus12_) = exp(LogGe) * sqrt(X/cPi) * 2.0

    !We calculate  \frac{x^{\nu+1}}{\Gamma(k+2)} (1+
    !     \sum_{i=1}^{2}{\zeta_series(2 i) \cdot
    !\frac{2 \Gamma(\nu+2)}{Gamma(\nu+2-2 i)} x^{-2 i}})
    do iNu = NuEqMinus12_ + 1, NuEq32_
       large_negative(iNu) = large_negative(iNu-1) * X / NuPlus1_I(iNu)
    end do
    X2Inv = 1.0/(X*X)

    do iNu = NuEqMinus12_ , NuEq32_
       large_negative(iNu) = large_negative(iNu)*(1.0 +&      ! 1st term
            NuPlus1_I(iNu) * (NuPlus1_I(iNu) - 1.0)*X2Inv*(&
            ZetaFunction2 +&                                  ! 2nd term
            (NuPlus1_I(iNu)-2.0)*(NuPlus1_I(iNu)-3.0)*&       ! 3rd term
            TwoZetaFunctionSeries4*X2Inv))
    end do
  end function large_negative
end module CRASH_ModFermiGas
