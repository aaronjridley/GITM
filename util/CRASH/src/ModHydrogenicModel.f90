!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModHydrogenicModel
  use ModConst
  use CRASH_ModAtomicNotation
  implicit none
  SAVE
  integer :: kMax = 100
  integer :: N_I(100),jPlusHalf_I(100)
  character(LEN=3) :: TypeTerm_I(100)=''
  integer :: nElectron_I(100)=0

  real,dimension(100,100)::Sigma_II
  character(LEN=*),parameter::NameFile = '../../../dataCRASH/AtomicData/sigma.dat'
  integer::nZ  !The nuclei charge.
contains
  function TypeConf()
    character(LEN=220)::TypeConf
    integer :: k,len
    len = 0
    do k=1,kMax
       if(nElectron_I(k)>0)then
          if(nElectron_I(k)<10)then
             write(TypeConf(len+1:220),'(a,i1,a)')TypeTerm_I(k),nElectron_I(k),' '
             len = len +5
          else
             write(TypeConf(len+1:220),'(a,i2,a)')TypeTerm_I(k),nElectron_I(k),' '
             len = len + 6
          end if
       end if
    end do
  end function TypeConf
  !==================================
  subroutine set_arrays(kMaxIn)
    integer, intent(in) :: kMaxIn
    integer:: N, L, K, nWrite
    !----------------
    kMax = min(kMaxIn,100)
    K=0
    !Loop over the principal quantum number
    do n=1,10
       if(n==10)then
          nWrite = 0
       else
          nWrite = n
       end if
       K = K +1

       N_I(K) = n
       !s- term:
       jPlusHalf_I(K) = 1
       write(TypeTerm_I(K),'(i1,a)') nWrite , 's+'
       if(k==kMax)return
       !Loop over the orbital quantum number
       do l=1, n-1
          K = K +1 
          N_I(K) = n
          jPlusHalf_I(K) = l
          write(TypeTerm_I(K),'(i1,a)') nWrite , TypeL_I(l)//'-'
          K = K +1 
          N_I(K) = n
          jPlusHalf_I(K) = l + 1
          write(TypeTerm_I(K),'(i1,a)') nWrite , TypeL_I(l)//'+'
          if(k==kMax)return
       end do
    end do
  end subroutine set_arrays
  !=====================================================
  real function hydrogen_level_energy(n,jPlusHalf,ZScreened)

    !Input parameters:
    !!Principal quantum number

    integer, intent(in) :: n 

    !j+1/2, also denoted as m, where j (=l+1/2 or l-1/2) is the total 
    !angular momentum, the latter being the sum of difference of the 
    !orbital angular momentum with the spin contribution. 
    !After all, j+1/2 = (l,l+1) for l/=0 or 1 for l=0

    integer, intent(in) :: jPlusHalf
    
    !Screened charge number. For a given relativistic subshell, k
    !ZSreened = nZ - \sum_{k^\prime}{ \sigma_{k k^\kprime}P_{\kprime}
    
    real, intent(in) :: ZScreened
    
    !The rest mass energy of the electron, in eV:

    real, parameter :: Emc2 = cElectronMass* cLightSpeed**2/cEV
    real, parameter :: cFine = 1.0/137.0360

    real :: AlphaZ !ZSreened multiplied by the ine structure constant
    !----------------------------!

    AlphaZ = ZScreened * cFine

    !The formula :
    !E =- m_e c^2 * ( sqrt(1 +(\alphaZ/(n-m+ sqrt(m^2- (alphaZ)^2) ) )^2) -1) 
    !where m = j+1/2

    hydrogen_level_energy = - &
          Emc2 * (&
          sqrt(1.0+ &
          (AlphaZ /&  !Enumerator
          (real(n - jPlusHalf) +sqrt(real(jPlusHalf)**2 - AlphaZ**2 ) )&
                         )**2  )  -1.0 )
  end function hydrogen_level_energy
  !====================================
  real function configuration_energy()
    integer:: k
    real    :: Charge_I(100),Energy_I(100)
    !----------------------------------!
    !Calculate screened charge, for each relativistic sub-shell

    Charge_I(1:kMax) = nZ - matmul(nElectron_I(1:kMax),Sigma_II(1:kMax,1:kMax))
    Energy_I = 0.0
    do k = 1,kMax

       !Balance self-screening (in the formula each of the electrtons screens itself.
       !The screening effect from k=k^\prime sub-shell shuld be reduced by 1

       Charge_I(k) = Charge_I(k) + Sigma_II(k,k)
       if(nElectron_I(k)==0)CYCLE
       Energy_I(k) = hydrogen_level_energy(N_I(k),jPlusHalf_I(k),Charge_I(k))
    end do
    
    configuration_energy = sum(Energy_I(1:kMax) * nElectron_I(1:kMax))
  end function configuration_energy
  !==================================
  subroutine read_sigma_coef
    use ModIoUnit,ONLY:io_unit_new
    integer::k,iUnit
    !---------
    iUnit = io_unit_new()
    open(iUnit,file=NameFile)
    read(iUnit,*)k
    call set_arrays(k)
  
    do k=1,kMax
       read(iUnit,*)Sigma_II(1:kMax,k)
    end do
    close(iUnit)
    
  end subroutine read_sigma_coef
  
end module ModHydrogenicModel
