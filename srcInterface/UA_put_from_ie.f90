!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==========================================================================
subroutine UA_put_from_ie(Buffer_II,iSize,jSize,NameVar)

  use CON_coupler, ONLY : IEt_Grid, check_allocate
  use ModInnerBc
  implicit none
  character(len=*), parameter :: NameSub='UA_put_from_ie'

  !\
  ! Buffer_II is a Latitude by Longitude array of potential.
  ! It starts at the pole and goes down to the equator for
  ! NameVar = 'PotNorth' and Equator to the pole for NameVar = 'PotSouth'.
  !
  ! iSize is the number of latitudes
  ! jSize is the number of longitudes
  !
  ! NameVar tells whether it is north or south, as described above.
  !/

  integer, intent(in) :: iSize,jSize
  real, intent(in) :: Buffer_II(iSize,jSize)
  character(len=*), intent(in) :: NameVar

  integer :: iError, i
  logical :: DoTest, DoTestMe

  !---------------------------------------------------------------------------

  call CON_set_do_test(NameSub,DoTest,DoTestMe)
  if(DoTest)write(*,*)NameSub,': NameVar,iSize,jSize=',NameVar,iSize,jSize

  !\
  ! Test to make sure that the code agrees with itself.
  !/

  if(  iSize /= IEt_Grid % Descriptor % nCells_D(1)+1 .or. &
       jSize /= IEt_Grid % Descriptor % nCells_D(2)+1 ) then

     write(*,*)NameSub//' grid sizes do not agree iSize,jSize,nCells=',&
          iSize,jSize,IEt_Grid % Descriptor % nCells_D(1:2)
     call CON_stop(NameSub//' SWMF_ERROR')

  end if

  !\
  ! Now lets start the interpolation into the UA coordinate system.
  !/

  if(IONO_nTheta < 0)then
     ! Allocate and calculate coordinates
     if(DoTest)write(*,*)NameSub,': allocating variables'
     IONO_nTheta = iSize
     IONO_nPsi   = jSize

     allocate(&
          IONO_NORTH_Theta(iSize,jSize), IONO_NORTH_Psi(iSize,jSize), &
          IONO_NORTH_Phi(iSize,jSize),                                &
          IONO_SOUTH_Theta(iSize,jSize), IONO_SOUTH_Psi(iSize,jSize), &
          IONO_SOUTH_Phi(iSize,jSize),                                &
          STAT = iError)
     call check_allocate(iError,NameSub//' IONO_NORTH/SOUTH_ThetaPsiPhi')

     if(DoTest)write(*,*)NameSub,': assigning coordinates'
     do i=1,iSize
        IONO_NORTH_Theta(i,:) = IEt_Grid % Coord1_I(i)
        IONO_SOUTH_Theta(i,:) = IEt_Grid % Coord1_I(i+iSize)
        IONO_NORTH_Psi(i,:)   = IEt_Grid % Coord2_I
        IONO_SOUTH_Psi(i,:)   = IEt_Grid % Coord2_I
     end do

  end if

  if(DoTest)write(*,*)NameSub,': putting potential'
  select case(NameVar)
  case('PotNorth')
     IONO_NORTH_Phi = Buffer_II
  case('PotSouth')
     IONO_SOUTH_Phi = Buffer_II
  case default
     call CON_stop(NameSub//' invalid NameVar='//NameVar)
  end select

  if(DoTest)write(*,*)NameSub,': done'

end subroutine UA_put_from_ie
