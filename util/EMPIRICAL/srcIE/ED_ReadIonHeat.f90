!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!------------------------------------------------------------------
!
!------------------------------------------------------------------

subroutine ReadIonHeat(filename, IsIonizationFile)

  use ModIons
  use ModIoUnit, ONLY: UnitTmp_

  implicit none

  character (len=100), intent(in) :: filename
  logical, intent(in)             :: IsIonizationFile

  integer ::            ilng,jlat,kalt,i
  character (len=120):: fileline

  logical :: IsFirstTime = .true.
  
!  write(*,*) "Reading Ion Precip file ",filename

  open(unit=UnitTmp_, file=TRIM(ADJUSTL(filename)))

  ! skip the first 3 lines
  do i=1, 3
     read(UnitTmp_, '(a120)') fileline
  end do
    
  ! number of altitude
  read(UnitTmp_, *) NAlt

  if (IsFirstTime) then
     ALLOCATE(glng(NLng))
     ALLOCATE(glat(NLat))
     ALLOCATE(alt(NAlt))
     IsFirstTime = .false.
  endif

  if (IsIonizationFile) then
     ALLOCATE(IonizationRate(NLng,NLat,NAlt))
  else
     ALLOCATE(HeatingRate(NLng,NLat,NAlt))
  endif

  ! altitude grid
  read(UnitTmp_, *)(alt(kalt), kalt=NAlt, 1, -1)
  read(UnitTmp_, '(a120)') fileline 

  do ilng=1, NLng
     do jlat=1, NLat

        if (IsIonizationFile) then
           read(UnitTmp_, *) glng(ilng), glat(jlat), &
                (IonizationRate(ilng,jlat,kalt),kalt=NAlt,1,-1)
        else
           read(UnitTmp_, *) glng(ilng), glat(jlat), &
                (HeatingRate(ilng,jlat,kalt),kalt=NAlt,-1)
        endif

        ! data(ilng,jlat,kalt) is on the geographic longitude of glng(ilng), 
        ! the geographic latitude of glat(jlat) and the altitude of alt(kalt)
   
     end do ! jlat
  end do ! ilng

  if (IsIonizationFile) then
     IonizationRate = IonizationRate*1.0e6
  else
     HeatingRate = HeatingRate*1.0e6*1.602e-19
  endif

  close(UnitTmp_)
 
end subroutine ReadIonHeat

!------------------------------------------------------------------
!
!------------------------------------------------------------------

subroutine interpolate_ions(nLons, nLats, nAlts, &
     GITMLons, GITMLats, GITMAlts, GITMIonRate, GITMHeatingRate)

  use ModIons
  use ModNumConst

  implicit none

  integer, intent(in) :: nLons, nLats, nAlts
  real, intent(in)    :: GITMLons(nLons)
  real, intent(in)    :: GITMLats(nLats)
  real, intent(in)    :: GITMAlts(nLons, nLats, nAlts)
  real, intent(out)   :: GITMIonRate(nLons, nLats, nAlts)
  real, intent(out)   :: GITMHeatingRate(nLons, nLats, nAlts)

  real :: rLons(nLons), LonsDeg(nLons)
  real :: rLats(nLats), LatsDeg(nLats)
  real :: rAlts(nAlts), AltsKms(nAlts)
  integer :: iLons(nLons)
  integer :: iLats(nLats)
  integer :: iAlts(nAlts)

  integer :: iLon, iLat, iAlt, i, j, k
  real    :: rlo, rla, ral
  !---------------------------------------------------------------------------
  GITMIonRate = 0.0
  GITMHeatingRate = 0.0

  LatsDeg = GITMLats*180/cPi
  LonsDeg = GITMLons*180/cPi

  if (minval(glat) > maxval(LatsDeg)) return

  rLons = 0.0
  rLats = 0.0
  rAlts = 0.0

  do iLon = 1, nLons
     rLons(iLon) = (LonsDeg(iLon)-glng(1)) / (glng(2) - glng(1)) + 1.0001
     iLons(iLon) = int(rLons(iLon))
     rLons(iLon) = rLons(iLon) - iLons(iLon)
     if (iLons(iLon) > NLng) iLons(iLon) = 0
  enddo

  do iLat = 1, nLats
     if (LatsDeg(iLat) <= glat(1)) iLats(iLat) = -1
     do i = 1, NLat
        if (LatsDeg(iLat) >= glat(i)) iLats(iLat) = i
     enddo
     i = iLats(iLat)
     if (i > 0) rLats(iLat) = (LatsDeg(iLat) - glat(i))/(glat(2) - glat(1))
  enddo

  do iLon = 1, nLons
     do iLat = 1, nLats
        if (iLats(iLat) < 1) CYCLE

        AltsKms = GITMAlts(iLon, iLat, :)/1000.0

        do iAlt = 1, nAlts
           if (AltsKms(iAlt) > alt(NAlt) .or. AltsKms(iAlt) < alt(1)) then
              iAlts(iAlt) = -1
           else
              do i = 1, NAlt
                 if (AltsKms(iAlt) >= alt(i)) iAlts(iAlt) = i
              enddo
           endif
           i = iAlts(iAlt)
           if (i >= 1) then 
              if (i == NAlt) then
                 rAlts(iAlt)= (AltsKms(iAlt)-alt(NAlt))/(alt(NAlt)-alt(NAlt-1))
              else
                 rAlts(iAlt)= (AltsKms(iAlt) - alt(i))/(alt(i+1) - alt(i))
              endif
           endif
        enddo

        do iAlt = 1, nAlts
           if (iAlts(iAlt) < 1 .or. iAlts(iAlt) >= nAlt) CYCLE

           i = iLons(iLon)
           j = iLats(iLat)
           k = iAlts(iAlt)
           rlo = rLons(iLon)
           rla = rLats(iLat)
           ral = rAlts(iAlt)

           GITMIonRate(iLon, iLat, iAlt) = &
                (1-rlo)*(1-rla)*(1-ral)*IonizationRate(i  ,j  ,k  ) + &
                (  rlo)*(1-rla)*(1-ral)*IonizationRate(i+1,j  ,k  ) + &
                (  rlo)*(  rla)*(1-ral)*IonizationRate(i+1,j+1,k  ) + &
                (1-rlo)*(  rla)*(1-ral)*IonizationRate(i  ,j+1,k  ) + &
                (1-rlo)*(1-rla)*(  ral)*IonizationRate(i  ,j  ,k+1) + &
                (  rlo)*(1-rla)*(  ral)*IonizationRate(i+1,j  ,k+1) + &
                (  rlo)*(  rla)*(  ral)*IonizationRate(i+1,j+1,k+1) + &
                (1-rlo)*(  rla)*(  ral)*IonizationRate(i  ,j+1,k+1)

           GITMHeatingRate(iLon, iLat, iAlt) = &
                (1-rlo)*(1-rla)*(1-ral)*HeatingRate(i  ,j  ,k  ) + &
                (  rlo)*(1-rla)*(1-ral)*HeatingRate(i+1,j  ,k  ) + &
                (  rlo)*(  rla)*(1-ral)*HeatingRate(i+1,j+1,k  ) + &
                (1-rlo)*(  rla)*(1-ral)*HeatingRate(i  ,j+1,k  ) + &
                (1-rlo)*(1-rla)*(  ral)*HeatingRate(i  ,j  ,k+1) + &
                (  rlo)*(1-rla)*(  ral)*HeatingRate(i+1,j  ,k+1) + &
                (  rlo)*(  rla)*(  ral)*HeatingRate(i+1,j+1,k+1) + &
                (1-rlo)*(  rla)*(  ral)*HeatingRate(i  ,j+1,k+1)

        enddo
     enddo
  enddo

end subroutine interpolate_ions
