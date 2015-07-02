!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! --------------------------------------------------------------------
!
!  This procedure feeds the previous electric potential pattern
!  into the background pattern for the current fit.
!
! --------------------------------------------------------------------

subroutine amiedt (ilat, ilon, ntime, etheta, ephi, potkv)

  use ModSize
  use ModAMIE
  use ModEField

  implicit none

  integer, intent(in) :: ilat, ilon, ntime
  real, intent(out) :: potkv, etheta, ephi

  if ((ntime.gt.1).and.(ilat.le.ithtrns)) then
     potkv  = Efield_Solution(ilon,ilat,potential_)/1000.0
     etheta = -1.0*Efield_Solution(ilon,ilat,efield_north_)
     ephi   = Efield_Solution(ilon,ilat,efield_east_)
  else
     potkv  = 0.0
     etheta = 0.0
     ephi   = 0.0
  endif

end subroutine amiedt


!!$! --------------------------------------------------------------------
!!$!
!!$!  This procedure is under development.
!!$!
!!$!  This takes the first potential pattern (if ntime > 1) and 
!!$!  changes it based on the statistical amie and the change in
!!$!  IMF from the first time period.  This is then used as the back ground
!!$!  pattern for the current inversion.
!!$!
!!$!  Possible Changes
!!$!  ----------------
!!$!  (1) Use current pattern with a better determination of how the change
!!$!      in potenial matches the change in IMF.
!!$!  (2) Use the current procedure and compute an average between
!!$!      that pattern and the previous pattern.
!!$!
!!$! --------------------------------------------------------------------
!!$
!!$subroutine amies_dt (iyear,imon,ida,ihr,imin,ilat,ilon,ntime,  &
!!$     ihsoln,ifrst,rmlat,rmlt,etheta,ephi,potkv,delta_ntime)
!!$
!!$      include 'param.h'
!!$      include 'qset.h'
!!$      include 'calc.h'
!!$      include 'amies.h'
!!$
!!$      character*100 filein
!!$      integer unitin
!!$      real*4 potkv, etheta, ephi
!!$      real*4 by,bz,sw
!!$
!!$      include_back = 0
!!$      if (ntime.eq.1) then
!!$         include_back = 1
!!$         by = bygsm(1)
!!$         bz = bzgsm(1)
!!$         sw = swv(1)
!!$         nsave = 1
!!$      else
!!$         if (delta_ntime.eq.-1) then
!!$            nsave = 1
!!$         else
!!$            if (delta_ntime.ge.ntime) then
!!$               nsave = 1
!!$            else
!!$               nsave = ntime-delta_ntime
!!$            endif
!!$         endif
!!$         by = bygsm(ntime)-bygsm(nsave)
!!$         bz = bzgsm(ntime)-bzgsm(nsave)
!!$         sw = swv(ntime)-swv(nsave)
!!$      endif
!!$
!!$      if (ntime.ge.2) then
!!$         save_efpot(ntime-1,ilon,ilat) = efpot(ilon,ilat)/1000.0
!!$         save_ee(ntime-1,ilon,ilat) = ee(ilon,ilat)
!!$         save_en(ntime-1,ilon,ilat) = en(ilon,ilat)
!!$      else
!!$         save_efpot(ntime,ilon,ilat) = 0.0
!!$         save_ee(ntime,ilon,ilat) = 0.0
!!$         save_en(ntime,ilon,ilat) = 0.0
!!$      endif
!!$
!!$      alte   = 6671.
!!$      alamn  = 58.
!!$      stepa  = 2.*alte*d2r
!!$
!!$      if (iamie_read_flag.ne.1) then
!!$         write(6,*) 'Reading AMIES parameters'
!!$         unitin = 99
!!$         call read_amies(unitin)
!!$         iamie_read_flag = 1
!!$      endif
!!$
!!$      if (ilat.le.ithtrns) then
!!$
!!$         potkv = 0.0
!!$         call amiemodel(rmlt, rmlat, by, bz, potkv, include_back)
!!$
!!$         potkv = potkv + save_efpot(nsave,ilon,ilat)
!!$
!!$         fmla = abs(rmlat)
!!$
!!$         kpol = 0
!!$         xmlt  = rmlt
!!$         xmlt1 = xmlt
!!$         fmla1 = fmla + 1.
!!$         if (fmla1 .gt. 90.) then
!!$            fmla1 = 180.  - fmla1
!!$            xmlt1 = xmlt1 + 12.
!!$         endif
!!$         p1 = 0.0
!!$         p2 = 0.0
!!$         call amiemodel(xmlt1, fmla1  , by, bz, p1, include_back)
!!$         call amiemodel(xmlt , fmla-1., by, bz, p2, include_back)
!!$
!!$         etheta = (p1 - p2) / stepa - save_en(nsave,ilon,ilat)
!!$
!!$!          calculate -(lon gradient).  For most latitudes, step 15 degrees
!!$!          in longitude (1 hr MLT = model resolution) along a great circle.
!!$!          However, limit minimum latitude to the model minimum, distorting
!!$!          the path onto a latitude line.  The step shrinks as latitude
!!$!          increases and would become zero at the pole, but a different
!!$!          scheme is used near the pole:  Assume lat=90 degrees and use
!!$!          Arts trick where Ephi(90,lon) = Etheta(90,lon+90.)
!!$
!!$         p1 = 0.0
!!$         p2 = 0.0
!!$
!!$         if (fmla .lt. 89.9) then
!!$            sl = sin (fmla*d2r)
!!$            cl = sqrt (1.-sl*sl)
!!$            sp = sin (15.*d2r)
!!$            sa = sqrt (1.-sl*sl*sp*sp)
!!$            fmla1 = acos(cl/sa)*r2d
!!$            call amiemodel(xmlt+1.,fmla1, by,bz, p1, include_back)
!!$            call amiemodel(xmlt-1.,fmla1, by,bz, p2, include_back)
!!$            step2 = 2.*alte*asin(cl*sp/sa)
!!$         else
!!$            step2 = stepa
!!$            xmlt = xmlt + 6.
!!$            if (xmlt.gt.24) xmlt = xmlt - 24.0
!!$            fmla = 90.
!!$            kpol = 1
!!$            xmlt1 = xmlt
!!$            fmla1 = fmla + 1.
!!$            if (fmla1 .gt. 90.) then
!!$               fmla1 = 180.  - fmla1
!!$               xmlt1 = xmlt1 + 12.
!!$               if (xmlt1.gt.24) xmlt1 = xmlt1 - 24.0
!!$            endif
!!$            call amiemodel(xmlt1, fmla1  , by,bz, p1, include_back)
!!$            call amiemodel(xmlt , fmla-1., by,bz, p2, include_back)
!!$         endif
!!$
!!$         ephi = (p2 - p1) / step2 + save_ee(nsave,ilon,ilat)
!!$         if (kpol .eq. 1) ephi = -ephi
!!$
!!$!          Below model minimum lat, the potential is value at min lat
!!$
!!$      else
!!$
!!$         potkv = 0.0
!!$         etheta = 0.0
!!$         ephi = 0.0
!!$
!!$      endif
!!$
!!$      end
!!$

! --------------------------------------------------------------------

subroutine amiespot (by,bz,rmlat,rmlt,etheta,ephi,potkv)

  use ModAMIE
  use ModEField
  use ModsAMIE
  use ModConstants

  implicit none

  real,intent(in) :: by, bz, rmlat, rmlt
  real, intent(out) :: potkv, etheta, ephi

  real :: alte, stepa, fmla, fmla1, xmlt, xmlt1, p1, p2
  real :: tmpsl, tmpcl, tmpsp, tmpsa, step2

  logical :: IncludeBackground, Pole

  IncludeBackground = .true.

  alte   = ri/1.0e3
  stepa  = 2.*alte*dtor

  call amiemodel(rmlt, rmlat, by, bz, potkv, IncludeBackground)

  fmla = abs(rmlat)

  Pole = .false.
  xmlt  = rmlt
  xmlt1 = xmlt
  fmla1 = fmla + 1.
  if (fmla1 > 90.) then
     fmla1 = 180.  - fmla1
     xmlt1 = xmlt1 + 12.
  endif
  call amiemodel(xmlt1, fmla1  , by,bz, p1, IncludeBackground)
  call amiemodel(xmlt , fmla-1., by,bz, p2, IncludeBackground)
  etheta = (p1 - p2) / stepa
  
  !          calculate -(lon gradient).  For most latitudes, step 15 degrees
  !          in longitude (1 hr MLT = model resolution) along a great circle.
  !          However, limit minimum latitude to the model minimum, distorting
  !          the path onto a latitude line.  The step shrinks as latitude
  !          increases and would become zero at the pole, but a different
  !          scheme is used near the pole:  Assume lat=90 degrees and use
  !          Arts trick where Ephi(90,lon) = Etheta(90,lon+90.)

  if (fmla < 89.9) then
     tmpsl = sin (fmla*dtor)
     tmpcl = sqrt (1.-tmpsl*tmpsl)
     tmpsp = sin (15.*dtor)
     tmpsa = sqrt (1.-tmpsl*tmpsl*tmpsp*tmpsp)
     fmla1 = acos(tmpcl/tmpsa)*rtod
     call amiemodel(xmlt+1.,fmla1, by,bz, p1, IncludeBackground)
     call amiemodel(xmlt-1.,fmla1, by,bz, p2, IncludeBackground)
     step2 = 2.*alte*asin(tmpcl*tmpsp/tmpsa)
  else
     step2 = stepa
     xmlt = xmlt + 6.
     fmla = 90.
     Pole = .true.
     xmlt1 = xmlt
     fmla1 = fmla + 1.
     if (fmla1 .gt. 90.) then
        fmla1 = 180.  - fmla1
        xmlt1 = xmlt1 + 12.
     endif
     call amiemodel(xmlt1, fmla1  , by,bz, p1, IncludeBackground)
     call amiemodel(xmlt , fmla-1., by,bz, p2, IncludeBackground)
  endif

  ephi = (p2 - p1) / step2

  if (Pole) ephi = -ephi

end subroutine amiespot

!--------------------------------------------------------------------------

subroutine amiemodel(inmlt, inmlat, by, bz, pot, IncludeBackground)

  use ModAMIE
  use ModEField
  use ModsAMIE

  implicit none

  real, intent(in) :: inmlt, inmlat, by, bz
  logical, intent(in) :: IncludeBackground
  real, intent(out) :: pot
  real :: rlat, rlon
  real :: lawt1, lawt2, lowt1, lowt2
  real :: poty_back_neg, poty_back_pos
  real :: potz_back_neg, potz_back_pos
  real :: back, potz, poty
  integer :: ilat1, ilat2, ilon1, ilon2

  if (DebugLevel > 3) then 
     write(*,*) "  ====> In subroutine amiemodel"
     write(*,*) "        inputs : ",inmlt,inmlat,by,bz,IncludeBackground
  endif

  ! turn latitude and mlt into grid point number
  ! with floating point, so if it is between a grid point, then we can
  ! interpolate

  rlat = (90.0-inmlat)/asdlat
  rlon = 15.0*inmlt/asdlon

  !     test for boundaries

  if (rlat < 0.0) rlat = 0.0
  if (rlat > iaslat-1) rlat = iaslat-1

  !     set points for interpolation

  ilat1 = rlat + 1
  ilat2 = ilat1 + 1
  if (ilat2 > iaslat) ilat2 = ilat1

  ilon1 = rlon + 1
  ilon2 = ilon1 + 1

  !     set weights for interpolation (all the wt2 values will be 0 if
  !     the lat and lon are on a grid point).

  lawt1 = 1.0 - ((rlat+1.0) - ilat1)
  lawt2 = 1.0 - (ilat2 - (rlat+1.0))

  lowt1 = 1.0 - ((rlon+1.0) - ilon1)
  lowt2 = 1.0 - (ilon2 - (rlon+1.0))

  if (ilon1 > iaslon) ilon1 = ilon1 - iaslon + 1
  if (ilon2 > iaslon) ilon2 = ilon2 - iaslon + 1

  if (IncludeBackground) then

     poty_back_neg =                                               &
          aspyni(ilon1,ilat1)*lowt1*lawt1+                         &
          aspyni(ilon1,ilat2)*lowt1*lawt2+                         &
          aspyni(ilon2,ilat1)*lowt2*lawt1+                         &
          aspyni(ilon2,ilat2)*lowt2*lawt2
     poty_back_pos =                                               &
          aspypi(ilon1,ilat1)*lowt1*lawt1+                         &
          aspypi(ilon1,ilat2)*lowt1*lawt2+                         &
          aspypi(ilon2,ilat1)*lowt2*lawt1+                         &
          aspypi(ilon2,ilat2)*lowt2*lawt2
         
     potz_back_neg =                                               &
          aspzni(ilon1,ilat1)*lowt1*lawt1+                         &
          aspzni(ilon1,ilat2)*lowt1*lawt2+                         &
          aspzni(ilon2,ilat1)*lowt2*lawt1+                         &
          aspzni(ilon2,ilat2)*lowt2*lawt2
     potz_back_pos =                                               &
          aspzpi(ilon1,ilat1)*lowt1*lawt1+                         &
          aspzpi(ilon1,ilat2)*lowt1*lawt2+                         &
          aspzpi(ilon2,ilat1)*lowt2*lawt1+                         &
          aspzpi(ilon2,ilat2)*lowt2*lawt2

     back = (poty_back_neg + poty_back_pos +                       &
          potz_back_neg + potz_back_pos)/                          &
          (4.0*(lowt1*lawt1+lowt1*lawt2+lowt2*lawt1+lowt2*lawt2))

  else

     back = 0.0

  endif

  if (by <= 0.0) then
     poty = by*aspyn(ilon1,ilat1)*lowt1*lawt1+                     &
          by*aspyn(ilon1,ilat2)*lowt1*lawt2+                       &
          by*aspyn(ilon2,ilat1)*lowt2*lawt1+                       &
          by*aspyn(ilon2,ilat2)*lowt2*lawt2
  else
     poty = by*aspyp(ilon1,ilat1)*lowt1*lawt1+                     &
          by*aspyp(ilon1,ilat2)*lowt1*lawt2+                       &
          by*aspyp(ilon2,ilat1)*lowt2*lawt1+                       &
          by*aspyp(ilon2,ilat2)*lowt2*lawt2
  endif

  poty = poty/(lowt1*lawt1+lowt1*lawt2+lowt2*lawt1+lowt2*lawt2)

  if (bz <= 0.0) then
     potz = bz*aspzn(ilon1,ilat1)*lowt1*lawt1+                     &
          bz*aspzn(ilon1,ilat2)*lowt1*lawt2+                       &
          bz*aspzn(ilon2,ilat1)*lowt2*lawt1+                       &
          bz*aspzn(ilon2,ilat2)*lowt2*lawt2
  else
     potz = bz*aspzp(ilon1,ilat1)*lowt1*lawt1+                     &
          bz*aspzp(ilon1,ilat2)*lowt1*lawt2+                       &
          bz*aspzp(ilon2,ilat1)*lowt2*lawt1+                       &
          bz*aspzp(ilon2,ilat2)*lowt2*lawt2
  endif

  potz = potz/(lowt1*lawt1+lowt1*lawt2+lowt2*lawt1+lowt2*lawt2)

  pot = (poty + potz + back)

end subroutine amiemodel

!--------------------------------------------------------------------------

subroutine read_amies(unitin)

  use ModAMIE
  use ModEField
  use ModsAMIE

  implicit none

  integer, intent(in) :: unitin
  character*100 ::line, fmt
  integer :: inlon, inlat, ierr, i, k

  read(unitin,'(2I8)',iostat = ierr) inlat,inlon

  if ((inlon /= iaslon).or.(inlat /= iaslat).or.(ierr < 0)) then

     write(6,*) 'Error in read_amies.f. Latitudes and longitude'
     write(6,*) 'do not match with code.'

     write(6,*) 'inlon : ',inlon,'  iaslon : ',iaslon
     write(6,*) 'inlat : ',inlat,'  iaslat : ',iaslat

     write(6,*) ierr

     stop

  else

     read(unitin,'(2F8.2)') asdlat, asdlon

     write(fmt,"(A1,I2,A5)") "(",inlon,"f8.4)"

     read(unitin,'(A100)') line
     do i=1,inlat
        read(unitin,fmt) (aspx(k,i),k=1,inlon)
     enddo

     read(unitin,'(A100)') line
     do i=1,inlat
        read(unitin,fmt) (aspxi(k,i),k=1,inlon)
     enddo

     read(unitin,'(A100)') line
     do i=1,inlat
        read(unitin,fmt) (aspyn(k,i),k=1,inlon)
     enddo

     read(unitin,'(A100)') line
     do i=1,inlat
        read(unitin,fmt) (aspyni(k,i),k=1,inlon)
     enddo

     read(unitin,'(A100)') line
     do i=1,inlat
        read(unitin,fmt) (aspyp(k,i),k=1,inlon)
     enddo

     read(unitin,'(A100)') line
     do i=1,inlat
        read(unitin,fmt) (aspypi(k,i),k=1,inlon)
     enddo

     read(unitin,'(A100)') line
     do i=1,inlat
        read(unitin,fmt) (aspzn(k,i),k=1,inlon)
     enddo

     read(unitin,'(A100)') line
     do i=1,inlat
        read(unitin,fmt) (aspzni(k,i),k=1,inlon)
     enddo

     read(unitin,'(A100)') line
     do i=1,inlat
        read(unitin,fmt) (aspzp(k,i),k=1,inlon)
     enddo

     read(unitin,'(A100)') line
     do i=1,inlat
        read(unitin,fmt) (aspzpi(k,i),k=1,inlon)
     enddo

     read(unitin,'(A100)') line
     do i=1,inlat
        read(unitin,fmt) (aspf(k,i),k=1,inlon)
     enddo

     read(unitin,'(A100)') line
     do i=1,inlat
        read(unitin,fmt) (aspfi(k,i),k=1,inlon)
     enddo

  endif

end subroutine read_amies
