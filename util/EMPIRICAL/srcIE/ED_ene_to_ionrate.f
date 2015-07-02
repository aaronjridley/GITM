      PROGRAM TEST1

      IMPLICIT NONE


      INCLUDE 'ED_R_elec_ed_lup_subs.inc'
      CHARACTER*80 DUMMY
      INTEGER*4 I,J,K
      INTEGER*4 ERR,ihold
      INTEGER*4 TEST_INDEX
      INTEGER*4 R_LOAD_EDEP_LOOKUP 
      INTEGER*4 R_EDEP_ALT_INDEX
      INTEGER*4 R_EDEP_ENG_INDEX
      integer*4 ind1,ind2
      REAL*4    R_EDEP_ALT_VALUE
      REAL*4    R_EDEP_ENG_VALUE
      real*4    R_LOG_INTERP
      real*4    R_SEMILOG_INTERP
      real*4    R_SEMILINEAR_INTERP
      real*4    R_LINEAR_INTERP
      REAL*4    alt_in_mbar
      REAL*4    alt_in_km
      REAL*4    alt_in_m
      REAL*4    alt_in_cm
      real*4    eng_in_erg
      real*4    eng_in_joule
      real*4    eng_in_ev
      real*4    eng_in_kev
      real*4    eng_in_mev
      real*4    test_alt
      real*4    test_eng
      real*4    x,x1,x2,y,y1,y2
      real*4    old_enrg_vals(ENERGY_LEVELS)
      real*4    old_flux_vals(ENERGY_LEVELS)
      real*4    old_erro_vals(ALTITUDES)
      real*4    old_prof_vals(ALTITUDES)
      real*4    old_alti_vals(ALTITUDES)
      real*4    new_prof_vals(ALTITUDES)
      real*4    alts(ALTITUDES)
      real*4    user_pres(15),hold

      real*4    ion_integral,zero,good_vals(altitudes)
      character*100 infile, outfile
      integer nfiles, ifile

      call init_all

      ERR = R_LOAD_EDEP_LOOKUP ()

      IF (ERR .NE. 1) THEN
        PRINT *,'FaTaL: DID NOT LOAD ENERGY DEPOSITION LOOKUP TABLES'
        STOP
      ENDIF

      nfiles = 1

      do ifile = 1,nfiles

         write(6,*) 'Enter input file name : '
         read(5,'(A100)') infile

         write(6,*) 'Enter output file name : '
         read(5,'(A100)') outfile

         open(unit=15,name=infile,form='formatted',type='old')

         read(15,'(A80)') dummy

         do i=1,energy_levels
            read(15,*) old_enrg_vals(i),old_flux_vals(i)
            if (old_enrg_vals(i).lt.1e-20) old_flux_vals(i) = 0.0
         enddo

         close(15)

         call R_ELEC_EDEP (old_flux_vals,30,old_enrg_vals,3,
     1        new_prof_vals,13)

         open(unit=16,name=outfile,type='unknown')

         write(16,*) 'ionization height profile ionizations/(cm**3-s)'

         write(16,'(I2)') altitudes

         do i=1,altitudes
            alt_in_km = R_EDEP_ALT_VALUE (I,3)
            alts(i) = alt_in_km
            if (new_prof_vals(i).lt.1e6) then
               good_vals(i) = new_prof_vals(i)
            else
               write(6,*) 'NaN found!'
               good_vals(i) = 
     |              (new_prof_vals(i-1)+new_prof_vals(i+1))/2.0
            endif
            write(16,'(1p,2e12.4)') alts(i),good_vals(i)
         enddo

         ion_integral = good_vals(1)*(alts(2)-alts(1))*1.0e5

         do i=2,altitudes
            ion_integral = ion_integral +
     |           good_vals(i)*(alts(i)-alts(i-1))*1.0e5
         enddo

         write(16,*) 'Integral'
         write(16,'(1p,e10.4)') ion_integral

         close(16)

      enddo

c  generate a user defined pressure grid
c      do i= 15,1,-1
c        user_pres(I) = 10.0**(float(4 - i))
c      enddo
c      print *
c      print *,'     New ionization rate mapped to user specified ',
c     1 'pressure grid -- '
c      print *
C  loop through each user defined pressure
c      do j= 1,15
c  set the first two low values of pressure
c        ind1 = 1
c        x1 = alog10(prof_pres(1)) - alog10(user_pres(j))
c        ind2 = 2
c        x2 = alog10(prof_pres(2)) - alog10(user_pres(j))
cc  switch if necessary so that x1<x2
c        if (abs(x1) .gt. abs(x2)) then
c          hold = x1
c          x1 = x2
c          x2 = hold
c          ihold = ind1
c          ind1 = ind2
c          ind2 = ihold
c        endif
cc  go through all profile values looking for closer pressures
c        do i=3,altitudes
c          x = alog10(prof_pres(i)) - alog10(user_pres(j))
cc  save lowest value in x1 and put the old x1 value into x2
c          if (abs(x) .le. abs(x1)) then
c            x2 = x1
c            ind2 = ind1
c            x1 = x
c            ind1 = i
c          else if (abs(x) .le. abs(x2)) then
cc  these are not the lowest, are they the second lowest?
c            x2 = x
c            ind2 = i
c          endif
c        enddo
cc  found lowest two, now interpolate/extrapolate as required
c        y = R_LOG_INTERP (user_pres(j),prof_pres(ind2),
c     1                  prof_pres(ind1),new_prof_vals(ind2),
c     2                  new_prof_vals(ind1))
cc  ionization rates can't be less than zero in this case
c        if (y .lt. 0.0) y = 0.0
cD       write(*,860) ind1,prof_pres(ind1),new_prof_vals(ind1)
c        write(*,850) user_pres(j),y
cD       write(*,860) ind2,prof_pres(ind2),new_prof_vals(ind2)
c      enddo
c  800 FORMAT(1x,i2,': ',f6.2,' km, ',f9.2,' m, ',f12.2,' cm, ',1pe10.3,
c     1     ' mbar')
c  810 FORMAT(1x,i2,': ',1pe9.2,' erg, ',e9.2,' Joule, ',e9.2,' eV, ',
c     1     e9.2,' keV, ',e9.2,' MeV')
c  830 format(1x,f6.2,' km: old =',1pe10.3,' ergs/(cm**3-s), new =,',
c     1     e10.3,' ergs/(cm**3-s)')
c  840 format(1x,1p,e10.3,' km: rate =',e10.3,' ionizations/(cm**3-s)')
c  850 format(1x,1pe10.3,' mbar: remapped =',
c     1     e10.3,' ionizations/(cm**3-s)')
c  860 format(1x,'from -> ',I2,': ',1pe10.3,' mbar AT ',
c     1     e10.3,' ionizations/(cm**3-s)')
c  900 format(a80)
      STOP
      END

c      print *,'Altitude tests --'
c      print *,'     check altitude matrix loaded properly'
c      print *,'     check requested units are proper'
c      print *,'     check index points to proper altitude'
c      DO I=1,ALTITUDES
c        alt_in_cm = R_EDEP_ALT_VALUE (I,1)
c        alt_in_m = R_EDEP_ALT_VALUE (I,2)
c        alt_in_km = R_EDEP_ALT_VALUE (I,3)
c        alt_in_mbar = R_EDEP_ALT_VALUE (I,4)
c        write (*,800) i,alt_in_km,alt_in_m,alt_in_cm,alt_in_mbar
c      ENDDO

c      print *,'     check altitude estimates index'

c      DO I=1,5
c        test_alt = 20.25/float(10**(i-1))
c        test_index = R_EDEP_ALT_INDEX (test_alt,4)
c        print *,'for a test altitude of ',test_alt,' mbar, ',
c     1          'the index found was ',test_index
c        test_alt = float(i)*20.25
c        test_index = R_EDEP_ALT_INDEX (test_alt,3)
c        print *,'for a test altitude of ',test_alt,' km, ',
c     1          'the index found was ',test_index
c        test_alt = test_alt*1000.0
c        test_index = R_EDEP_ALT_INDEX (test_alt,2)
c        print *,'for a test altitude of ',test_alt,' m, ',
c     1          'the index found was ',test_index
c        test_alt = test_alt*100.0
c        test_index = R_EDEP_ALT_INDEX (test_alt,1)
c        print *,'for a test altitude of ',test_alt,' cm, ',
c     1          'the index found was ',test_index
c      enddo

c      print *
c      print *
c      print *,'Energy tests --'
c      print *,'     check energy matrix loaded properly'
c      print *,'     check requested units are proper'
c      print *,'     check index points to proper energy'
c      DO I=1,ENERGY_LEVELS
c        eng_in_erg = R_EDEP_ENG_VALUE (I,1)
c        eng_in_joule = R_EDEP_ENG_VALUE (I,2)
c        eng_in_ev = R_EDEP_ENG_VALUE (I,3)
c        eng_in_kev = R_EDEP_ENG_VALUE (I,4)
c        eng_in_mev = R_EDEP_ENG_VALUE (I,5)
c        write (*,810) i,eng_in_erg,eng_in_joule,eng_in_ev,eng_in_kev,
c     1                eng_in_mev 
c      ENDDO

c      print *,'     check energy estimates index'
c      DO I=10,50,10
c        test_eng = float(i)/202.5
c        test_index = R_EDEP_ENG_INDEX (test_eng,5)
c        print *,'for a test energy of ',test_eng,' MeV, ',
c     1          'the index found was ',test_index
c        test_eng = test_eng*1000.0
c        test_index = R_EDEP_ENG_INDEX (test_eng,4)
c        print *,'for a test energy of ',test_eng,' keV, ',
c     1          'the index found was ',test_index
c        test_eng = test_eng*1000.0
c        test_index = R_EDEP_ENG_INDEX (test_eng,3)
c        print *,'for a test energy of ',test_eng,' eV, ',
c     1          'the index found was ',test_index
c        test_eng = test_eng*1.602E-19
c        test_index = R_EDEP_ENG_INDEX (test_eng,2)
c        print *,'for a test energy of ',test_eng,' Joule, ',
c     1          'the index found was ',test_index
c        test_eng = test_eng*1.0E+7
c        test_index = R_EDEP_ENG_INDEX (test_eng,1)
c        print *,'for a test energy of ',test_eng,' erg, ',
c     1          'the index found was ',test_index
c      enddo

c      print *
c      print *
c      print *,'Interpolation tests --'
c      print *,'     for x1,y1 = (4.0,1000.0) and x2,y2 = (2.0,100.0)'
c      x2 = 2.0
c      x1 = 4.0
c      y2 = 100.0
c      y1 = 1000.0
c      x = 3.0
c      y = R_LOG_INTERP (x,x2,x1,y2,y1)
c      print *,'     log interpolation at x=',x,' is ',y
c      y = R_SEMILOG_INTERP (x,x2,x1,y2,y1)
c      print *,'     semilog interpolation at x=',x,' is ',y
c      y = R_SEMILINEAR_INTERP (x,x2,x1,y2,y1)
c      print *,'     semilinear interpolation at x=',x,' is ',y
c      y = R_LINEAR_INTERP (x,x2,x1,y2,y1)
c      print *,'     linear interpolation at x=',x,' is ',y
c      print *,'     for x1,y1 = (2.0,1.0) and x2,y2 = (4.0,1000.0)'
c      print *,'     for x1,y1 = (2.0,100.0) and x2,y2 = (4.0,1000.0)'
c      x1 = 2.0
c      x2 = 4.0
c      y1 = 100.0
c      y2 = 1000.0
c      x = 3.0
c      y = R_LOG_INTERP (x,x2,x1,y2,y1)
c      print *,'     log interpolation at x=',x,' is ',y
c      y = R_SEMILOG_INTERP (x,x2,x1,y2,y1)
c      print *,'     semilog interpolation at x=',x,' is ',y
c      y = R_SEMILINEAR_INTERP (x,x2,x1,y2,y1)
c      print *,'     semilinear interpolation at x=',x,' is ',y
c      y = R_LINEAR_INTERP (x,x2,x1,y2,y1)
c      print *,'     linear interpolation at x=',x,' is ',y
c      print *,'     for x1,y1 = (2.0,1.0) and x2,y2 = (4.0,1000.0)'
c      y1 = 1.0
c      y = R_LOG_INTERP (x,x2,x1,y2,y1)
c      print *,'     log interpolation at x=',x,' is ',y
c      y = R_SEMILOG_INTERP (x,x2,x1,y2,y1)
c      print *,'     semilog interpolation at x=',x,' is ',y
c      y = R_SEMILINEAR_INTERP (x,x2,x1,y2,y1)
c      print *,'     semilinear interpolation at x=',x,' is ',y
c      y = R_LINEAR_INTERP (x,x2,x1,y2,y1)
c      print *,'     linear interpolation at x=',x,' is ',y
c      print *
c      print *
c      print *,'Extrapolation tests --'
c      print *,'     for x1,y1 = (4.0,150.0) and x2,y2 = (1.0,100.0)'
c      x2 = 1.0
c      x1 = 4.0
c      y2 = 100.0
c      y1 = 150.0
c      x = 0.2
c      y = R_LOG_INTERP (x,x2,x1,y2,y1)
c      print *,'     log extrapolation at x=',x,' is ',y
c      y = R_SEMILOG_INTERP (x,x2,x1,y2,y1)
c      print *,'     semilog extrapolation at x=',x,' is ',y
c      y = R_SEMILINEAR_INTERP (x,x2,x1,y2,y1)
c      print *,'     semilinear extrapolation at x=',x,' is ',y
c      y = R_LINEAR_INTERP (x,x2,x1,y2,y1)
c      print *,'     linear extrapolation at x=',x,' is ',y
c      print *,'     for x1,y1 = (2.0,1.0) and x2,y2 = (4.0,150.0)'
c      print *,'     for x1,y1 = (1.0,100.0) and x2,y2 = (4.0,150.0)'
c      x1 = 1.0
c      x2 = 4.0
c      y1 = 100.0
c      y2 = 150.0
c      x = 0.2
c      y = R_LOG_INTERP (x,x2,x1,y2,y1)
c      print *,'     log extrapolation at x=',x,' is ',y
c      y = R_SEMILOG_INTERP (x,x2,x1,y2,y1)
c      print *,'     semilog extrapolation at x=',x,' is ',y
c      y = R_SEMILINEAR_INTERP (x,x2,x1,y2,y1)
c      print *,'     semilinear extrapolation at x=',x,' is ',y
c      y = R_LINEAR_INTERP (x,x2,x1,y2,y1)
c      print *,'     linear extrapolation at x=',x,' is ',y
c      print *,'     for x1,y1 = (2.0,1.0) and x2,y2 = (4.0,150.0)'
c      y1 = 1.0
c      x = 7.1
c      y = R_LOG_INTERP (x,x2,x1,y2,y1)
c      print *,'     log extrapolation at x=',x,' is ',y
c      y = R_SEMILOG_INTERP (x,x2,x1,y2,y1)
c      print *,'     semilog extrapolation at x=',x,' is ',y
c      y = R_SEMILINEAR_INTERP (x,x2,x1,y2,y1)
c      print *,'     semilinear extrapolation at x=',x,' is ',y
c      y = R_LINEAR_INTERP (x,x2,x1,y2,y1)
c      print *,'     linear extrapolation at x=',x,' is ',y
c      print *
c      print *

c      print *,'Energy Deposition tests --'
