C Subroutine ROUT writes model atmosphere and excitation rates to an
C output file in order to transfer them to radiative transfer program.
C
C This software is part of the GLOW model.  Use is governed by the Open Source
C Academic Research License Agreement contained in the file glowlicense.txt.
C For more information see the file glow.txt.
C
C SMB/SCS, 9/94
C Replaced 834 with LBH, SCS, 2/03
C Reduced cascade contribution to 1356, SCS, 9/03
C Included radiative recombination in 1356 but commented out, SCS, 9/03
C
      SUBROUTINE ROUT(ROFILE,LUN,EF,EZ,ITAIL,FRACO,FRACO2,FRACN2)

      use Mod_GLOW

      dimension zout(jmax), zhe(jmax), e1356(jmax), e1304(jmax),
     >          e1027(jmax), e989(jmax), elbh(jmax)
      character*40 rofile
C
      do 50 j=1,jmax
      zout(j)=zz(j)/1.e5
      zhe(j)=0.
C Temporary fix-up to 1356:
      e1356(j)=aglw(3,1,j)+aglw(5,1,j)*0.5 !+ecalc(j)*zxden(3,j)*4.9e-13
      e1304(j)=aglw(4,1,j)+aglw(6,1,j)+aglw(7,1,j)+aglw(8,1,j)*0.1
      e1027(j)=aglw(7,1,j)
      e989(j)=aglw(8,1,j)
      elbh(j)=aglw(4,3,j)
  50  continue
C
      open(unit=lun,file=rofile,status='unknown')
C  
      write(lun,100)
 100  format('   JMAX ','   SZA  ','   UT   ','   IDATE',
     >             '   LAT  ','   LONG ','    DIP ')
      write(lun,200) jmax,sza*180./3.14159,ut,idate,glat,glong,dip
 200  format(i8,f8.2,f8.1,i8,3f8.2)
      write(lun,250)
 250  format('  F107  ','  F107p ','  F107a ','  HLyBr ',
     >       ' FeXVIr ','   HLYa ','  HeIew ',' XUVfac ')
      write(lun,300) f107,f107p,f107a,hlybr,fexvir,hlya,heiew,xuvfac
 300  format (5f8.2,1p,e10.2,0p,2f8.2)
      write(lun,350)
 350  format(' Eflux  ',' Ezero  ',' Itail  ',
     >       ' FracO  ',' FracO2  ',' FracN2  ')
      write(lun,400) ef, ez, itail, fraco, fraco2, fracn2
 400  format(f8.2,f8.1,i8,3f8.2)
      write(lun,500)
 500  format(' Alt    Tn    Ti    Te  ',
     >       '    O        O2       N2       He   ',
     >       '    N        Ne       O+      1356  ',
     >       '   1304     1027      989     LBH')
C
      do 700,j=1,jmax
        write(lun,600) zout(j),ztn(j),zti(j),zte(j),
     >               zo(j),zo2(j),zn2(j),zhe(j),
     >               zns(j),ecalc(j),zxden(3,j),e1356(j),
     >               e1304(j),e1027(j),e989(j),elbh(j)
 600    format(0p,f6.1,3f6.0,1p,12e9.2)
 700  continue
C
      close(lun)
      return
      end
