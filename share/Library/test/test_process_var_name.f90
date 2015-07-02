!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program process_var_name_test

  ! This test is designed to be perfomed using, in turn, all the ModEquation
  ! files in BATSRUS.
  ! Automatic looping over all cases is done in the Makefile.

  ! The test converts the string array NameVar_V, defined in a specific
  ! ModEquation, to the same format it will take in a full SWMF run.

  use ModProcessVarName,  ONLY: process_var_name
  use ModVarIndexes,      ONLY: NameVar_V
  use ModUtilities,       ONLY: lower_case
  implicit none

  integer  :: nDensity, nSpeed, nP, nPpar, nVar, iVar
  integer  :: nWaveName, nMaterialName
  character(len=15),allocatable :: NameVarFixed_V(:)
  ! ----------------------------------------------------------
  nVar = size(NameVar_V,1)
  ! Replace '?' characters with numbers (for waves/materials)
  call set_namevar 
 
  ! Fix string array to comply with usual input to process_var_name.
  ! String length should be fixed(=15) regardless of length in ModEquation.
  ! All characters should be lower case.
  allocate(NameVarFixed_V(nVar))
  do iVar = 1, nVar
     NameVarFixed_V(iVar) = NameVar_V(iVar)
     call lower_case(NameVarFixed_V(iVar))
  end do

  call process_var_name(nVar, NameVarFixed_V, &
       nDensity, nSpeed, nP, nPpar, nWaveName, nMaterialName)

  write(*,*) 'Original  Standardized  nDensity=',nDensity,'nSpeed=',nSpeed
  write(*,*) '--------  ------------'
  do iVar = 1, nVar
     write(*,*) NameVar_V(iVar),'         ',trim(NameVarFixed_V(iVar))
  end do

  deallocate(NameVarFixed_V)

contains

 subroutine set_namevar

    integer :: iWave
    character(len=3):: NameWave
    integer :: iMaterial
    character(len=2):: NameMaterial
    !-------------------------------------------------------------------------           
    iWave = 1
    iMaterial = 1
    do iVar = 1, nVar
  
       ! Fix the NameVar_V string for waves                                              
       if(index(NameVar_V(iVar),'I?') >= 1) then
          write(NameWave,'(a,i2.2)') 'i',iWave
          NameVar_V(iVar) = NameWave
          iWave = iWave + 1
       end if

       ! Fix the NameVar_V string for material levels
       if(index(NameVar_V(iVar),'M?') >= 1) then
          write(NameMaterial,'(a,i1.1)') 'm',iMaterial
          NameVar_V(iVar) = NameMaterial
          iMaterial = iMaterial + 1
       end if

    end do

  end subroutine set_namevar

end program process_var_name_test
! =================================================================
subroutine CON_stop(String)
  implicit none
  character(len=*), intent(in):: String
  write(*,*)'ERROR: ',String
  stop
end subroutine CON_stop
! =================================================================
