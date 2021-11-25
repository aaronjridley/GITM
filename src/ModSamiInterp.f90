! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!******************************************
!******************************************

!             find_point_test

!******************************************
!******************************************

module ModSamiInterp

  use parameter_mod, only: nz,nf,nlt
    
    real,parameter :: re = 6370.0
    real :: baltst_g(nz,nf,nlt),blonst_g(nz,nf,nlt),blatst_g(nz,nf,nlt)
    real :: SAMIVars_g(10,nz,nf,nlt)
    real :: gitm_mlon,gitm_mlat,gitm_alt
    real, allocatable :: mlon0(:),mlat0(:)
    integer :: i0,j0,k0,l0,i1,j1,k1,l1 
    real :: sami_lon0,sami_lon1,sami_lat0,sami_lat1
    
    real :: SAMIPhi_g(nf,nlt)
    real :: phi_gitm_mlon,phi_gitm_mlat
    integer :: phi_i0,phi_j0,phi_i1,phi_j1
    integer :: phi_flag
    real :: phi_sami_lon0,phi_sami_lon1, &
        phi_sami_lat0,phi_sami_lat1
    integer :: flag_eq

contains

subroutine process2(i_init,j_init,j,k,l)
    
    implicit none

    integer,intent(out) :: j,k,l
    integer,intent(in) :: i_init, j_init
    
    real, allocatable :: alt1(:)
    integer :: iAlt,i,ifl
    logical :: jflag 
    logical :: jflag2 
    real :: lat_k0,lat_l0,lat_l1

    jflag = .true.
    jflag2 = .true.


    allocate(alt1(1:nz))

    i = i_init
    j = j_init+1

    do while ( jflag )

        do iAlt = 1, nz

        alt1(iAlt) = baltst_g(iAlt,j,i)-re
        !alt1(iAlt) = baltst_g(iAlt,j,i)

        enddo

        if(gitm_alt <= maxval(alt1)) then
            
            alt1 = baltst_g(:,j-1,i)-re

            if (gitm_alt <= maxval(alt1)) then
                
                !print*,'check0', gitm_alt,maxval(alt1),minval(alt1)
                jflag = .false.
            endif
        endif
        
        j = j+1
    enddo

    j = j-1
    
    !print*,'check1: alt_j0,alt_j1 max',&
    !    gitm_alt,maxval(baltst_g(:,j-1,i)-re),maxval(baltst_g(:,j,i)-re)
    
    !j1 = j
    !j0 = j-1

    
    !j = j-1 !j0

 
    ifl =j
    
    k = -1
    l = -1

    do while (jflag2)

        call FindAltIndex(i,ifl,l)

        if (l>0) then
            ! gitm_alt le max(alt1_j0) but gitm_alt le max(alt_j1)
            lat_l0 = blatst_g(l,ifl,i)

            if (gitm_mlat < 0) then
                lat_l1 = blatst_g(l+1,ifl,i)
            else
                lat_l1 = blatst_g(l-1,ifl,i)
            endif

            if(abs(gitm_mlat) <= abs(lat_l0)) then
                
                !if (abs(gitm_mlat) > abs(lat_l1)) &
                !    print*,'GITM Point go outside SAMI3 grid( GT lat_l1)',&
                !         gitm_mlon,gitm_mlat,gitm_alt

                call FindAltIndex(i,ifl-1,k)
                lat_k0 = blatst_g(k,ifl-1,i)

                if ((abs(gitm_mlat) <  abs(lat_k0)) &
                    .or.(k<0)) stop
                jflag2 = .false.
            endif
            
            !l = l0
        endif


        ifl = ifl + 1

        if (ifl>= nf) then
            jflag2 = .false.
            print*,'GITM nan',gitm_mlon,gitm_mlat,gitm_alt
        endif
    enddo

    j = ifl - 1

    j = j-1
      
    deallocate(alt1)

endsubroutine process2


subroutine get_value_from_8points(gitm_x,gitm_y,gitm_z,data_inp)

    implicit none

    real :: gitm_x,gitm_y,gitm_z
    !integer :: i0,j0,k0,l0,i1,j1,k1,l1
    real  :: xijk     ,yijk     ,zijk, &
             xijk_1   ,yijk_1   ,zijk_1,&
             xij_1l   ,yij_1l   ,zij_1l , &
             xij_1l_1 ,yij_1l_1 ,zij_1l_1, &
             xi_1jk     ,yi_1jk     ,zi_1jk, &
             xi_1jk_1   ,yi_1jk_1   ,zi_1jk_1, &
             xi_1j_1l   ,yi_1j_1l   ,zi_1j_1l, &
             xi_1j_1l_1 ,yi_1j_1l_1 ,zi_1j_1l_1
    
    real :: xd00,xd01,xd10,xd11
    real :: x00,x01,x10,x11
    real :: y00,y01,y10,y11
    real :: z00,z01,z10,z11
    real :: data00,data01,data10,data11
    
    real :: yd0,yd1
    real :: x0,x1
    real :: y0,y1
    real :: z0,z1
    real :: data0,data1

    real :: zd
    real :: x_inp,y_inp,z_inp,data_inp
    !real  :: baltst_g(nz,nf,nlt),blonst_g(nz,nf,nlt),blatst_g(nz,nf,nlt)
   
    !blonst_g = SAMIVars_g(1,:,:,:)
    !blatst_g = SAMIVars_g(2,:,:,:)
    !baltst_g = SAMIVars_g(3,:,:,:)

    call sph_to_car(blonst_g(k0,j0,i0),blatst_g(k0,j0,i0),baltst_g(k0,j0,i0),&
                    xijk     ,yijk     ,zijk)
    call sph_to_car(blonst_g(k1,j0,i0),blatst_g(k1,j0,i0),baltst_g(k1,j0,i0),&
                    xijk_1   ,yijk_1   ,zijk_1)
    call sph_to_car(blonst_g(l0,j1,i0),blatst_g(l0,j1,i0),baltst_g(l0,j1,i0),&
                    xij_1l   ,yij_1l   ,zij_1l)
    call sph_to_car(blonst_g(l1,j1,i0),blatst_g(l1,j1,i0),baltst_g(l1,j1,i0),&
                    xij_1l_1 ,yij_1l_1 ,zij_1l_1)
    
    call sph_to_car(blonst_g(k0,j0,i1),blatst_g(k0,j0,i1),baltst_g(k0,j0,i1),&
                    xi_1jk     ,yi_1jk     ,zi_1jk)
    call sph_to_car(blonst_g(k1,j0,i1),blatst_g(k1,j0,i1),baltst_g(k1,j0,i1), &
                    xi_1jk_1   ,yi_1jk_1   ,zi_1jk_1 )
    call sph_to_car(blonst_g(l0,j1,i1),blatst_g(l0,j1,i1),baltst_g(l0,j1,i1), &
                    xi_1j_1l   ,yi_1j_1l   ,zi_1j_1l)
    call sph_to_car(blonst_g(l1,j1,i1),blatst_g(l1,j1,i1),baltst_g(l1,j1,i1), &
                    xi_1j_1l_1 ,yi_1j_1l_1 ,zi_1j_1l_1)
    
    !interp along x
    xd00 = (gitm_x-xijk)    /(xi_1jk    -xijk)
    xd01 = (gitm_x-xijk_1)  /(xi_1jk_1  -xijk_1)
    xd10 = (gitm_x-xij_1l)  /(xi_1j_1l  -xij_1l)
    xd11 = (gitm_x-xij_1l_1)/(xi_1j_1l_1-xij_1l_1)


    !data00 = data[i0,j0,k0] *(1.-xd00)+data[i1,j0,k0] *xd00
    !data01 = data[i0,j0,k1] *(1.-xd01)+data[i1,j0,k1] *xd01
    !data10 = data[i0,j1,l0] *(1.-xd10)+data[i1,j1,l0] *xd10
    !data11 = data[i0,j1,l1] *(1.-xd11)+data[i1,j1,l1] *xd11

    data00 = 1. *(1.- xd00) + 1. * xd00
    data01 = 1. *(1.- xd01) + 1. * xd01
    data10 = 1. *(1.- xd10) + 1. * xd10
    data11 = 1. *(1.- xd11) + 1. * xd11
    
    x00 = xijk      *(1.-xd00) + xi_1jk     *xd00
    x01 = xijk_1    *(1.-xd01) + xi_1jk_1   *xd01
    x10 = xij_1l    *(1.-xd10) + xi_1j_1l   *xd10
    x11 = xij_1l_1  *(1.-xd11) + xi_1j_1l_1 *xd11


    y00 = yijk     *(1.-xd00) + yi_1jk     *xd00
    y01 = yijk_1   *(1.-xd01) + yi_1jk_1   *xd01
    y10 = yij_1l   *(1.-xd10) + yi_1j_1l   *xd10
    y11 = yij_1l_1 *(1.-xd11) + yi_1j_1l_1 *xd11


    z00 = zijk     *(1.-xd00) + zi_1jk     *xd00
    z01 = zijk_1   *(1.-xd01) + zi_1jk_1   *xd01
    z10 = zij_1l   *(1.-xd10) + zi_1j_1l   *xd10
    z11 = zij_1l_1 *(1.-xd11) + zi_1j_1l_1 *xd11
    
    !interp along y

    yd0 = (gitm_y - y00)/(y10 - y00)
    yd1 = (gitm_y - y01)/(y11 - y01)

    data0 = data00*(1.-yd0)+data10*yd0
    data1 = data01*(1.-yd1)+data11*yd1

    x0 = x00*(1.-yd0)+x10*yd0
    x1 = x01*(1.-yd1)+x11*yd1

    y0 = y00*(1.-yd0)+y10*yd0
    y1 = y01*(1.-yd1)+y11*yd1

    z0 = z00*(1.-yd0)+z10*yd0
    z1 = z01*(1.-yd1)+z11*yd1

    !interp along z
    
    zd = (gitm_z - z0)/(z1-z0)
    data_inp = data0*(1.-zd)+data1*zd

    x_inp = x0*(1.-zd) + x1*zd
    y_inp = y0*(1.-zd) + y1*zd
    z_inp = z0*(1.-zd) + z1*zd
    
    !print*,x_inp,y_inp,z_inp

    return

end subroutine get_value_from_8points

subroutine interp_with_4points(xd00,xd10,yd0,data_inp)

    use ModGITM,only: iproc

    implicit none

    real, intent(out) :: data_inp
    
    real :: xd00,xd10
    real :: data00,data10
    
    real :: yd0


    !interp along x
       
       data00 = SAMIPhi_g(j0,i0) *(1.-xd00)+SAMIPhi_g(j0,i1) *xd00
       data10 = SAMIPhi_g(j1,i0) *(1.-xd10)+SAMIPhi_g(j1,i1) *xd10
   
    !interp along y  
       data_inp = data00*(1.-yd0)+data10*yd0

    return

end subroutine interp_with_4points


subroutine interp_with_8points(xd00,xd01,xd10,xd11,yd0,yd1,zd, &
               gitm_x,gitm_y,gitm_z,data_inp)

    use ModGITM,only: iproc

    implicit none

    !real, intent(in) :: var(nz,nf,nlt)
    real, intent(out) :: data_inp(10)
    real :: gitm_x,gitm_y,gitm_z
    
    real :: xd00,xd01,xd10,xd11
    real :: data00,data01,data10,data11
    
    real :: yd0,yd1
    real :: data0,data1

    real :: zd
    integer :: ivar

    !interp along x

    do ivar = 1, 10 

       !data00 = denit1_g(k0,j0,i0) *(1.-xd00)+denit1_g(k0,j0,i1) *xd00
       !data01 = denit1_g(k1,j0,i0) *(1.-xd01)+denit1_g(k1,j0,i1) *xd01
       !data10 = denit1_g(l0,j1,i0) *(1.-xd10)+denit1_g(l0,j1,i1) *xd10
       !data11 = denit1_g(l1,j1,i0) *(1.-xd11)+denit1_g(l1,j1,i1) *xd11
    
       data00 = SAMIVars_g(ivar,k0,j0,i0) *(1.-xd00)+SAMIVars_g(ivar,k0,j0,i1) *xd00
       data01 = SAMIVars_g(ivar,k1,j0,i0) *(1.-xd01)+SAMIVars_g(ivar,k1,j0,i1) *xd01
       data10 = SAMIVars_g(ivar,l0,j1,i0) *(1.-xd10)+SAMIVars_g(ivar,l0,j1,i1) *xd10
       data11 = SAMIVars_g(ivar,l1,j1,i0) *(1.-xd11)+SAMIVars_g(ivar,l1,j1,i1) *xd11
   
       !if (ivar ==1 ) then
       !    print*,'along x',iproc,data00,data01,data10,data11
       !    print*,SAMIVars_g(ivar,k0,j0,i0),SAMIVars_g(ivar,k0,j0,i1)
       !    print*,SAMIVars_g(ivar,k1,j0,i0),SAMIVars_g(ivar,k1,j0,i1)
       !    print*,SAMIVars_g(ivar,l0,j1,i0),SAMIVars_g(ivar,l0,j1,i1)
       !    print*,SAMIVars_g(ivar,l1,j1,i0),SAMIVars_g(ivar,l1,j1,i1)
       !   
       !    print*,'--xd00,xd01,xd10,xd11 b',xd00,xd01,xd10,xd11
       !
       !endif
       !interp along y

    
       data0 = data00*(1.-yd0)+data10*yd0
    
       data1 = data01*(1.-yd1)+data11*yd1

       !if (ivar ==1 ) print*,'along y',iproc,data0,data1
    
       !interp along z
    
    
       data_inp(ivar) = data0*(1.-zd)+data1*zd
       !if (ivar ==1 ) print*,'along z',iproc,data_inp(ivar)
    enddo


    return

end subroutine interp_with_8points

subroutine get_factors_from_4points_sph(gitm_x,gitm_y,xd00,xd10,yd0)

    use ModGITM, only: iproc

    implicit none

    real  :: gitm_x,gitm_y
    real  :: xijk     ,yijk   , &
             xij_1l   ,yij_1l   , &
             xi_1jk     ,yi_1jk     , &
             xi_1j_1l   ,yi_1j_1l  
    
    real,intent(out) :: xd00,xd10,yd0
    real :: x00,x10
    real :: y00,y10

    real :: x_inp,y_inp
    
    if (phi_gitm_mlat > 0) then

      yijk     = blatst_g(nz,phi_j0,phi_i0)
      yij_1l   = blatst_g(nz,phi_j1,phi_i0)
      
      yi_1jk     = blatst_g(nz,phi_j0,phi_i1)
      yi_1j_1l   = blatst_g(nz,phi_j1,phi_i1)
  else

      yijk     = blatst_g(1,phi_j0,phi_i0)
      yij_1l   = blatst_g(1,phi_j1,phi_i0)
      
      yi_1jk     = blatst_g(1,phi_j0,phi_i1)
      yi_1j_1l   = blatst_g(1,phi_j1,phi_i1)
  endif

  !print*,'before',iproc,phi_i0,phi_i1,nlt,gitm_x,gitm_y,phi_sami_lon0,phi_sami_lon1
  if ((phi_i0 == 1).and.(phi_i1==nlt)) then

        if (gitm_x < phi_sami_lon0) then
            phi_sami_lon0 = phi_sami_lon0 + 360.0
            gitm_x = gitm_x + 360.0
        else if (gitm_x > phi_sami_lon1) then
            phi_sami_lon0 = phi_sami_lon0+360.0
        endif
        
    endif
    
  !print*,'after',iproc,phi_i0,phi_i1,nlt,gitm_x,gitm_y,phi_sami_lon0,phi_sami_lon1
    if ((phi_j0 == 1).and.(phi_j1==1)) then
        if (phi_gitm_mlat <0 ) then
            
            yij_1l = -yijk
            yi_1j_1l = -yi_1jk
        else
            yijk = -yijk
            yi_1jk = -yi_1jk
        endif
    endif

    xijk         = phi_sami_lon0
    xij_1l       = phi_sami_lon0

    xi_1jk         = phi_sami_lon1
    xi_1j_1l       = phi_sami_lon1
   
    !print*,'==>>> 4 points', gitm_x,gitm_y
    !print*, '-- 1 ',         xijk     ,yijk     
    !print*, '-- 2 ',         xij_1l   ,yij_1l   
    !print*, '-- 3 ',         xi_1jk     ,yi_1jk     
    !print*, '-- 4 ',         xi_1j_1l   ,yi_1j_1l   


    !print*,'--> a sph_to_car',iproc
    !interp along x
    xd00 = (gitm_x-xijk)    /(xi_1jk    -xijk)

    xd10 = xd00
    
    x00 = xijk      *(1.-xd00) + xi_1jk     *xd00
    x10 = xij_1l    *(1.-xd10) + xi_1j_1l   *xd10

    y00 = yijk     *(1.-xd00) + yi_1jk     *xd00
    y10 = yij_1l   *(1.-xd10) + yi_1j_1l   *xd10
    
    !print*,'--xd00,xd01,xd10,xd11 a',xd00,xd01,xd10,xd11
    !print*,'-- x00,x01,x10,x11',x00,x01,x10,x11
    !print*,'-- y00,y01,y10,y11',y00,y01,y10,y11
    !print*,'-- z00,z01,z10,z11',z00,z01,z10,z11
    !interp along y

    yd0 = (gitm_y - y00)/(y10 - y00)

    x_inp = x00*(1.-yd0)+x10*yd0
    y_inp = y00*(1.-yd0)+y10*yd0

    !print*,'-- yd0,yd1',yd0,yd1
    !print*,'-- x0,x1',x0,x1
    !print*,'-- y0,y1',y0,y1
    !print*,'-- z0,z1',z0,z1
    !          if ((abs(mod(gitm_x,360.0)-(mod(x_inp,360.0)))>2).or. &
    !          (abs(gitm_y-y_inp)>4).or.&
    !          (abs(gitm_z-(z_inp))>10)) then
    !              print*,'===>> Magnetic Lon/Lat find',mod(gitm_x,360.0),gitm_y,gitm_z-re 
    !         
    !print*,'-- inp point',mod(x_inp,360.0),y_inp,z_inp-re,iproc
    !          stop
    !          endif
    !print*,'-- inp point',mod(x_inp,360.0),y_inp,iproc

end subroutine get_factors_from_4points_sph


subroutine get_factors_from_8points_sph_3(gitm_x,gitm_y,gitm_z,xd00,xd01,xd10,xd11,&
        !yd0,yd1,zd)
        yd0,yd1,zd,x_inp,y_inp,z_inp)

    use ModGITM, only: iproc

    implicit none

    real :: gitm_x,gitm_y,gitm_z
    !integer :: i0,j0,k0,l0,i1,j1,k1,l1
    real  :: xijk     ,yijk     ,zijk, &
             xijk_1   ,yijk_1   ,zijk_1,&
             xij_1l   ,yij_1l   ,zij_1l , &
             xij_1l_1 ,yij_1l_1 ,zij_1l_1, &
             xi_1jk     ,yi_1jk     ,zi_1jk, &
             xi_1jk_1   ,yi_1jk_1   ,zi_1jk_1, &
             xi_1j_1l   ,yi_1j_1l   ,zi_1j_1l, &
             xi_1j_1l_1 ,yi_1j_1l_1 ,zi_1j_1l_1
    
    real :: xd00,xd01,xd10,xd11
    real :: x00,x01,x10,x11
    real :: y00,y01,y10,y11
    real :: z00,z01,z10,z11
    real :: data00,data01,data10,data11
    
    real :: yd0,yd1
    real :: x0,x1
    real :: y0,y1
    real :: z0,z1
    real :: data0,data1

    real :: zd
    real :: x_inp,y_inp,z_inp,data_inp
 
    integer :: inz_ijk,inz_ijk_1,inz_ij_1l,inz_ij_1l_1, &
               inz_i_1jk,inz_i_1jk_1,inz_i_1j_1l,inz_i_1j_1l_1

    if (flag_eq == 1) then

        inz_ijk     = k0
        inz_ijk_1   = k1
        inz_ij_1l   = l0
        inz_ij_1l_1 = l1

        inz_i_1jk     = k0
        inz_i_1jk_1   = k1
        inz_i_1j_1l   = l0
        inz_i_1j_1l_1 = l1


    else
        if (gitm_mlat > 0) then
            inz_ijk     = nz
            inz_ijk_1   = nz
            inz_ij_1l   = nz
            inz_ij_1l_1 = nz

            inz_i_1jk     = nz
            inz_i_1jk_1   = nz
            inz_i_1j_1l   = nz
            inz_i_1j_1l_1 = nz

        else
            inz_ijk     = 1
            inz_ijk_1   = 1
            inz_ij_1l   = 1
            inz_ij_1l_1 = 1
    
            inz_i_1jk     = 1
            inz_i_1jk_1   = 1
            inz_i_1j_1l   = 1
            inz_i_1j_1l_1 = 1

        endif
    endif

      yijk     = blatst_g(inz_ijk     ,j0,i0)
      yijk_1   = blatst_g(inz_ijk_1   ,j0,i0)
      yij_1l   = blatst_g(inz_ij_1l   ,j1,i0)
      yij_1l_1 = blatst_g(inz_ij_1l_1 ,j1,i0)
      
      yi_1jk     = blatst_g(inz_i_1jk  ,   j0,i1)
      yi_1jk_1   = blatst_g(inz_i_1jk_1,   j0,i1)
      yi_1j_1l   = blatst_g(inz_i_1j_1l,   j1,i1)
      yi_1j_1l_1 = blatst_g(inz_i_1j_1l_1, j1,i1)
   
    if ((i0 == 1).and.(i1==nlt)) then
        
        !DoPrint =.true.

        !print*,'-- gitm_x,sami_lon0,a ',gitm_x,sami_lon0
        
        if (gitm_x < sami_lon0) then
            sami_lon0 = sami_lon0 + 360.0
            gitm_x = gitm_x + 360.0
        elseif (gitm_x > sami_lon1) then
            sami_lon0 = sami_lon0+360.0
        endif

        !print*,'-- gitm_x,sami_lon0,b ',gitm_x,sami_lon0
        
    endif

    xijk         = sami_lon0
    xijk_1       = sami_lon0
    xij_1l       = sami_lon0
    xij_1l_1     = sami_lon0

    xi_1jk         = sami_lon1
    xi_1jk_1       = sami_lon1
    xi_1j_1l       = sami_lon1
    xi_1j_1l_1     = sami_lon1
    
    zijk       = baltst_g(k0,j0,i0)
    zijk_1     = baltst_g(k1,j0,i0)
    zij_1l     = baltst_g(l0,j1,i0)
    zij_1l_1   = baltst_g(l1,j1,i0)

    zi_1jk     = baltst_g(k0,j0,i1)
    zi_1jk_1   = baltst_g(k1,j0,i1)
    zi_1j_1l   = baltst_g(l0,j1,i1)
    zi_1j_1l_1 = baltst_g(l1,j1,i1)
  
    !if (flag_eq == 1) then
    !print*,'==>>> 8 points', gitm_x,gitm_y,gitm_z-re
    !print*, '-- 1 ',         xijk     ,yijk     ,zijk-re
    !print*, '-- 2 ',         xijk_1   ,yijk_1   ,zijk_1-re
    !print*, '-- 3 ',         xij_1l   ,yij_1l   ,zij_1l-re
    !print*, '-- 4 ',         xij_1l_1 ,yij_1l_1 ,zij_1l_1-re
    !print*, '-- 5 ',         xi_1jk     ,yi_1jk     ,zi_1jk-re
    !print*, '-- 6 ',         xi_1jk_1   ,yi_1jk_1   ,zi_1jk_1-re
    !print*, '-- 7 ',         xi_1j_1l   ,yi_1j_1l   ,zi_1j_1l-re
    !print*, '-- 8 ',         xi_1j_1l_1 ,yi_1j_1l_1 ,zi_1j_1l_1-re

   !endif
    !print*,yijk,yi_1jk

    !print*,'--> a sph_to_car',iproc
    !interp along x
    xd00 = (gitm_x-xijk)    /(xi_1jk    -xijk) 
    !xd01 = (gitm_x-xijk_1)  /(xi_1jk_1  -xijk_1)
    !xd10 = (gitm_x-xij_1l)  /(xi_1j_1l  -xij_1l)
    !xd11 = (gitm_x-xij_1l_1)/(xi_1j_1l_1-xij_1l_1)
    
    !if (DoPrint) then
    !    print*,'-- xd00 ', iproc,(gitm_x-xijk)    /(xi_1jk    -xijk) 
    !
    !    DoPrint = .false.
    !endif
    xd01 = xd00
    xd10 = xd00
    xd11 = xd00

    x00 = xijk      *(1.-xd00) + xi_1jk     *xd00
    x01 = xijk_1    *(1.-xd01) + xi_1jk_1   *xd01
    x10 = xij_1l    *(1.-xd10) + xi_1j_1l   *xd10
    x11 = xij_1l_1  *(1.-xd11) + xi_1j_1l_1 *xd11

    y00 = yijk     *(1.-xd00) + yi_1jk     *xd00
    y01 = yijk_1   *(1.-xd01) + yi_1jk_1   *xd01
    y10 = yij_1l   *(1.-xd10) + yi_1j_1l   *xd10
    y11 = yij_1l_1 *(1.-xd11) + yi_1j_1l_1 *xd11


    z00 = zijk     *(1.-xd00) + zi_1jk     *xd00
    z01 = zijk_1   *(1.-xd01) + zi_1jk_1   *xd01
    z10 = zij_1l   *(1.-xd10) + zi_1j_1l   *xd10
    z11 = zij_1l_1 *(1.-xd11) + zi_1j_1l_1 *xd11
    
    !print*,'--xd00,xd01,xd10,xd11 a',xd00,xd01,xd10,xd11
    !print*,'-- x00,x01,x10,x11',x00,x01,x10,x11
    !print*,'-- y00,y01,y10,y11',y00,y01,y10,y11
    !print*,'-- z00,z01,z10,z11',z00,z01,z10,z11
    !interp along y

    yd0 = (gitm_y - y00)/(y10 - y00)
    yd1 = (gitm_y - y01)/(y11 - y01)

    !data0 = data00*(1.-yd0)+data10*yd0
    !data1 = data01*(1.-yd1)+data11*yd1

    x0 = x00*(1.-yd0)+x10*yd0
    x1 = x01*(1.-yd1)+x11*yd1

    y0 = y00*(1.-yd0)+y10*yd0
    y1 = y01*(1.-yd1)+y11*yd1

    z0 = z00*(1.-yd0)+z10*yd0
    z1 = z01*(1.-yd1)+z11*yd1

    !print*,'-- yd0,yd1',yd0,yd1
    !print*,'-- x0,x1',x0,x1
    !print*,'-- y0,y1',y0,y1
    !print*,'-- z0,z1',z0,z1
    !interp along z
    
    zd = (gitm_z - z0)/(z1-z0)
    !data_inp = data0*(1.-zd)+data1*zd

    !print*,'-- zd',zd
    x_inp = x0*(1.-zd) + x1*zd
    y_inp = y0*(1.-zd) + y1*zd
    z_inp = z0*(1.-zd) + z1*zd

end subroutine get_factors_from_8points_sph_3


subroutine get_factors_from_8points_sph_2(gitm_x,gitm_y,gitm_z,xd00,xd01,xd10,xd11,&
        yd0,yd1,zd,x_inp,y_inp,z_inp)

    use ModGITM, only: iproc

    implicit none

    real :: gitm_x,gitm_y,gitm_z
    !integer :: i0,j0,k0,l0,i1,j1,k1,l1
    real  :: xijk     ,yijk     ,zijk, &
             xijk_1   ,yijk_1   ,zijk_1,&
             xij_1l   ,yij_1l   ,zij_1l , &
             xij_1l_1 ,yij_1l_1 ,zij_1l_1, &
             xi_1jk     ,yi_1jk     ,zi_1jk, &
             xi_1jk_1   ,yi_1jk_1   ,zi_1jk_1, &
             xi_1j_1l   ,yi_1j_1l   ,zi_1j_1l, &
             xi_1j_1l_1 ,yi_1j_1l_1 ,zi_1j_1l_1
    
    real :: xd00,xd01,xd10,xd11
    real :: x00,x01,x10,x11
    real :: y00,y01,y10,y11
    real :: z00,z01,z10,z11
    real :: data00,data01,data10,data11
    
    real :: yd0,yd1
    real :: x0,x1
    real :: y0,y1
    real :: z0,z1
    real :: data0,data1

    real :: zd
    real :: x_inp,y_inp,z_inp,data_inp
  
    
    !xijk     = blonst_g(k0,j0,i0)
    !xijk_1   = blonst_g(k1,j0,i0)
    !xij_1l   = blonst_g(l0,j1,i0)
    !xij_1l_1 = blonst_g(l1,j1,i0)

    !xi_1jk     = blonst_g(k0,j0,i1)
    !xi_1jk_1   = blonst_g(k1,j0,i1)
    !xi_1j_1l   = blonst_g(l0,j1,i1)
    !xi_1j_1l_1 = blonst_g(l1,j1,i1)

    
    if (gitm_mlat > 0) then

      yijk     = blatst_g(nz,j0,i0)
      yijk_1   = blatst_g(nz,j0,i0)
      yij_1l   = blatst_g(nz,j1,i0)
      yij_1l_1 = blatst_g(nz,j1,i0)
      
      yi_1jk     = blatst_g(nz,j0,i1)
      yi_1jk_1   = blatst_g(nz,j0,i1)
      yi_1j_1l   = blatst_g(nz,j1,i1)
      yi_1j_1l_1 = blatst_g(nz,j1,i1)
  else

      yijk     = blatst_g(1,j0,i0)
      yijk_1   = blatst_g(1,j0,i0)
      yij_1l   = blatst_g(1,j1,i0)
      yij_1l_1 = blatst_g(1,j1,i0)
      
      yi_1jk     = blatst_g(1,j0,i1)
      yi_1jk_1   = blatst_g(1,j0,i1)
      yi_1j_1l   = blatst_g(1,j1,i1)
      yi_1j_1l_1 = blatst_g(1,j1,i1)
  endif
      
    !yijk     = blatst_g(k0,j0,i0)
    !yijk_1   = blatst_g(k1,j0,i0)
    !yij_1l   = blatst_g(l0,j1,i0)
    !yij_1l_1 = blatst_g(l1,j1,i0)

    !yi_1jk     = blatst_g(k0,j0,i1)
    !yi_1jk_1   = blatst_g(k1,j0,i1)
    !yi_1j_1l   = blatst_g(l0,j1,i1)
    !yi_1j_1l_1 = blatst_g(l1,j1,i1)

    if ((i0 == 1).and.(i1==nlt)) then

        if (gitm_x < sami_lon0) then
            sami_lon0 = sami_lon0 + 360.0
            gitm_x = gitm_x + 360.0
        elseif (gitm_x > sami_lon1) then
            sami_lon0 = sami_lon0+360.0
        endif
        
    endif

    !if ((j0 == 1).and.(j1==1)) then
    !    if (gitm_mlat <0 ) then
    !        sami_lat1 = -sami_lat0
    !    else
    !        sami_lat0 = -sami_lat0

    !    endif
    !endif

    if ((j0 == 1).and.(j1==1)) then
        if (gitm_mlat <0 ) then
            yij_1l = -yijk
            yij_1l_1 = -yijk_1

            yi_1j_1l = -yi_1jk
            yi_1j_1l_1 = -yi_1jk_1
        else
            yijk = -yijk
            yijk_1 = -yijk_1

            yi_1jk = -yi_1jk
            yi_1jk_1 = -yi_1jk_1
        endif
    endif

    !print*,sami_lon0,gitm_x,sami_lon1,sami_lat0,gitm_y,sami_lat1
    
   ! print*,yijk,yijk_1,yij_1l,yij_1l_1,yi_1jk,yi_1jk_1,yi_1j_1l,yi_1j_1l_1
    !elseif ((j0 == nlt).and.(j1==nlt)) then
    !    if (gitm_mlat <0 ) then 
    !        sami_lat1 = -(sami_lat0 + 180.0)
    !    else
    !        sami_lat1 = -sami_lat0 + 180.0
    !    endif
    !endif

    xijk         = sami_lon0
    xijk_1       = sami_lon0
    xij_1l       = sami_lon0
    xij_1l_1     = sami_lon0

    xi_1jk         = sami_lon1
    xi_1jk_1       = sami_lon1
    xi_1j_1l       = sami_lon1
    xi_1j_1l_1     = sami_lon1

    !yijk         = sami_lat0
    !yijk_1       = sami_lat0
    !yij_1l       = sami_lat0
    !yij_1l_1     = sami_lat0

    !yi_1jk         = sami_lat1
    !yi_1jk_1       = sami_lat1
    !yi_1j_1l       = sami_lat1
    !yi_1j_1l_1     = sami_lat1
    
    zijk       = baltst_g(k0,j0,i0)
    zijk_1     = baltst_g(k1,j0,i0)
    zij_1l     = baltst_g(l0,j1,i0)
    zij_1l_1   = baltst_g(l1,j1,i0)

    zi_1jk     = baltst_g(k0,j0,i1)
    zi_1jk_1   = baltst_g(k1,j0,i1)
    zi_1j_1l   = baltst_g(l0,j1,i1)
    zi_1j_1l_1 = baltst_g(l1,j1,i1)
   
    !print*,'==>>> 8 points', gitm_x,gitm_y,gitm_z-re
    !print*, '-- 1 ',         xijk     ,yijk     ,zijk-re
    !print*, '-- 2 ',         xijk_1   ,yijk_1   ,zijk_1-re
    !print*, '-- 3 ',         xij_1l   ,yij_1l   ,zij_1l-re
    !print*, '-- 4 ',         xij_1l_1 ,yij_1l_1 ,zij_1l_1-re
    !print*, '-- 5 ',         xi_1jk     ,yi_1jk     ,zi_1jk-re
    !print*, '-- 6 ',         xi_1jk_1   ,yi_1jk_1   ,zi_1jk_1-re
    !print*, '-- 7 ',         xi_1j_1l   ,yi_1j_1l   ,zi_1j_1l-re
    !print*, '-- 8 ',         xi_1j_1l_1 ,yi_1j_1l_1 ,zi_1j_1l_1-re

    !print*,yijk,yi_1jk

    !print*,'--> a sph_to_car',iproc
    !interp along x
    xd00 = (gitm_x-xijk)    /(xi_1jk    -xijk)
    !xd01 = (gitm_x-xijk_1)  /(xi_1jk_1  -xijk_1)
    !xd10 = (gitm_x-xij_1l)  /(xi_1j_1l  -xij_1l)
    !xd11 = (gitm_x-xij_1l_1)/(xi_1j_1l_1-xij_1l_1)

    xd01 = xd00
    xd10 = xd00
    xd11 = xd00

    !data00 = data[i0,j0,k0] *(1.-xd00)+data[i1,j0,k0] *xd00
    !data01 = data[i0,j0,k1] *(1.-xd01)+data[i1,j0,k1] *xd01
    !data10 = data[i0,j1,l0] *(1.-xd10)+data[i1,j1,l0] *xd10
    !data11 = data[i0,j1,l1] *(1.-xd11)+data[i1,j1,l1] *xd11

    
    x00 = xijk      *(1.-xd00) + xi_1jk     *xd00
    x01 = xijk_1    *(1.-xd01) + xi_1jk_1   *xd01
    x10 = xij_1l    *(1.-xd10) + xi_1j_1l   *xd10
    x11 = xij_1l_1  *(1.-xd11) + xi_1j_1l_1 *xd11

    y00 = yijk     *(1.-xd00) + yi_1jk     *xd00
    y01 = yijk_1   *(1.-xd01) + yi_1jk_1   *xd01
    y10 = yij_1l   *(1.-xd10) + yi_1j_1l   *xd10
    y11 = yij_1l_1 *(1.-xd11) + yi_1j_1l_1 *xd11


    z00 = zijk     *(1.-xd00) + zi_1jk     *xd00
    z01 = zijk_1   *(1.-xd01) + zi_1jk_1   *xd01
    z10 = zij_1l   *(1.-xd10) + zi_1j_1l   *xd10
    z11 = zij_1l_1 *(1.-xd11) + zi_1j_1l_1 *xd11
    
    !print*,'--xd00,xd01,xd10,xd11 a',xd00,xd01,xd10,xd11
    !print*,'-- x00,x01,x10,x11',x00,x01,x10,x11
    !print*,'-- y00,y01,y10,y11',y00,y01,y10,y11
    !print*,'-- z00,z01,z10,z11',z00,z01,z10,z11
    !interp along y

    yd0 = (gitm_y - y00)/(y10 - y00)
    yd1 = (gitm_y - y01)/(y11 - y01)

    !data0 = data00*(1.-yd0)+data10*yd0
    !data1 = data01*(1.-yd1)+data11*yd1

    x0 = x00*(1.-yd0)+x10*yd0
    x1 = x01*(1.-yd1)+x11*yd1

    y0 = y00*(1.-yd0)+y10*yd0
    y1 = y01*(1.-yd1)+y11*yd1

    z0 = z00*(1.-yd0)+z10*yd0
    z1 = z01*(1.-yd1)+z11*yd1

    !print*,'-- yd0,yd1',yd0,yd1
    !print*,'-- x0,x1',x0,x1
    !print*,'-- y0,y1',y0,y1
    !print*,'-- z0,z1',z0,z1
    !interp along z
    
    zd = (gitm_z - z0)/(z1-z0)
    !data_inp = data0*(1.-zd)+data1*zd

    !print*,'-- zd',zd
    x_inp = x0*(1.-zd) + x1*zd
    y_inp = y0*(1.-zd) + y1*zd
    z_inp = z0*(1.-zd) + z1*zd
    
    !          if ((abs(mod(gitm_x,360.0)-(mod(x_inp,360.0)))>2).or. &
    !          (abs(gitm_y-y_inp)>4).or.&
    !          (abs(gitm_z-(z_inp))>10)) then
    !              print*,'===>> Magnetic Lon/Lat find',mod(gitm_x,360.0),gitm_y,gitm_z-re 
    !         
    !print*,'-- inp point',mod(x_inp,360.0),y_inp,z_inp-re,iproc
    !          stop
    !          endif
    !print*,'-- inp point',mod(x_inp,360.0),y_inp,z_inp-re,iproc
    
    !print*,'end get_factors_from_8points',iproc

end subroutine get_factors_from_8points_sph_2


subroutine get_factors_from_8points_sph(gitm_x,gitm_y,gitm_z,xd00,xd01,xd10,xd11,&
        yd0,yd1,zd)

    use ModGITM, only: iproc

    implicit none

    real :: gitm_x,gitm_y,gitm_z
    !integer :: i0,j0,k0,l0,i1,j1,k1,l1
    real  :: xijk     ,yijk     ,zijk, &
             xijk_1   ,yijk_1   ,zijk_1,&
             xij_1l   ,yij_1l   ,zij_1l , &
             xij_1l_1 ,yij_1l_1 ,zij_1l_1, &
             xi_1jk     ,yi_1jk     ,zi_1jk, &
             xi_1jk_1   ,yi_1jk_1   ,zi_1jk_1, &
             xi_1j_1l   ,yi_1j_1l   ,zi_1j_1l, &
             xi_1j_1l_1 ,yi_1j_1l_1 ,zi_1j_1l_1
    
    real :: xd00,xd01,xd10,xd11
    real :: x00,x01,x10,x11
    real :: y00,y01,y10,y11
    real :: z00,z01,z10,z11
    real :: data00,data01,data10,data11
    
    real :: yd0,yd1
    real :: x0,x1
    real :: y0,y1
    real :: z0,z1
    real :: data0,data1

    real :: zd
    real :: x_inp,y_inp,z_inp,data_inp
   

    !call sph_to_car(blonst_g(k0,j0,i0),blatst_g(k0,j0,i0),baltst_g(k0,j0,i0),&
    !                xijk     ,yijk     ,zijk)
    !call sph_to_car(blonst_g(k1,j0,i0),blatst_g(k1,j0,i0),baltst_g(k1,j0,i0),&
    !                xijk_1   ,yijk_1   ,zijk_1)
    !call sph_to_car(blonst_g(l0,j1,i0),blatst_g(l0,j1,i0),baltst_g(l0,j1,i0),&
    !                xij_1l   ,yij_1l   ,zij_1l)
    !call sph_to_car(blonst_g(l1,j1,i0),blatst_g(l1,j1,i0),baltst_g(l1,j1,i0),&
    !                xij_1l_1 ,yij_1l_1 ,zij_1l_1)
    !
    !call sph_to_car(blonst_g(k0,j0,i1),blatst_g(k0,j0,i1),baltst_g(k0,j0,i1),&
    !                xi_1jk     ,yi_1jk     ,zi_1jk)
    !call sph_to_car(blonst_g(k1,j0,i1),blatst_g(k1,j0,i1),baltst_g(k1,j0,i1), &
    !                xi_1jk_1   ,yi_1jk_1   ,zi_1jk_1 )
    !call sph_to_car(blonst_g(l0,j1,i1),blatst_g(l0,j1,i1),baltst_g(l0,j1,i1), &
    !                xi_1j_1l   ,yi_1j_1l   ,zi_1j_1l)
    !call sph_to_car(blonst_g(l1,j1,i1),blatst_g(l1,j1,i1),baltst_g(l1,j1,i1), &
    !                xi_1j_1l_1 ,yi_1j_1l_1 ,zi_1j_1l_1)
  
    
    xijk     = blonst_g(k0,j0,i0)
    xijk_1   = blonst_g(k1,j0,i0)
    xij_1l   = blonst_g(l0,j1,i0)
    xij_1l_1 = blonst_g(l1,j1,i0)

    xi_1jk     = blonst_g(k0,j0,i1)
    xi_1jk_1   = blonst_g(k1,j0,i1)
    xi_1j_1l   = blonst_g(l0,j1,i1)
    xi_1j_1l_1 = blonst_g(l1,j1,i1)

    yijk     = blatst_g(k0,j0,i0)
    yijk_1   = blatst_g(k1,j0,i0)
    yij_1l   = blatst_g(l0,j1,i0)
    yij_1l_1 = blatst_g(l1,j1,i0)

    yi_1jk     = blatst_g(k0,j0,i1)
    yi_1jk_1   = blatst_g(k1,j0,i1)
    yi_1j_1l   = blatst_g(l0,j1,i1)
    yi_1j_1l_1 = blatst_g(l1,j1,i1)
 

    zijk       = baltst_g(k0,j0,i0)
    zijk_1     = baltst_g(k1,j0,i0)
    zij_1l     = baltst_g(l0,j1,i0)
    zij_1l_1   = baltst_g(l1,j1,i0)

    zi_1jk     = baltst_g(k0,j0,i1)
    zi_1jk_1   = baltst_g(k1,j0,i1)
    zi_1j_1l   = baltst_g(l0,j1,i1)
    zi_1j_1l_1 = baltst_g(l1,j1,i1)

    !xij_1l = blonst_g(l0,j1,i0)
    !xij_1l_1 = blonst_g(l1,j1,i0)




    !print*,'--- gitm_xyz',gitm_x,gitm_y,gitm_z
    !print*,'1',blonst_g(k0,j0,i0),xijk 
    !print*,'2',blonst_g(k1,j0,i0),xijk_1
    !print*,'3',blonst_g(l0,j1,i0),xij_1l 
    !print*,'4',blonst_g(l1,j1,i0),xij_1l_1
    !print*,'5',blonst_g(k0,j0,i1),xi_1jk 
    !print*,'6',blonst_g(k1,j0,i1),xi_1jk_1 
    !print*,'7',blonst_g(l0,j1,i1),xi_1j_1l 
    !print*,'8',blonst_g(l1,j1,i1),xi_1j_1l_1 


    !print*,'--> a sph_to_car',iproc
    !interp along x
    xd00 = (gitm_x-xijk)    /(xi_1jk    -xijk)
    !xd01 = (gitm_x-xijk_1)  /(xi_1jk_1  -xijk_1)
    !xd10 = (gitm_x-xij_1l)  /(xi_1j_1l  -xij_1l)
    !xd11 = (gitm_x-xij_1l_1)/(xi_1j_1l_1-xij_1l_1)

    xd01 = xd00
    xd10 = xd00
    xd11 = xd00

    !data00 = data[i0,j0,k0] *(1.-xd00)+data[i1,j0,k0] *xd00
    !data01 = data[i0,j0,k1] *(1.-xd01)+data[i1,j0,k1] *xd01
    !data10 = data[i0,j1,l0] *(1.-xd10)+data[i1,j1,l0] *xd10
    !data11 = data[i0,j1,l1] *(1.-xd11)+data[i1,j1,l1] *xd11

    
    x00 = xijk      *(1.-xd00) + xi_1jk     *xd00
    x01 = xijk_1    *(1.-xd01) + xi_1jk_1   *xd01
    x10 = xij_1l    *(1.-xd10) + xi_1j_1l   *xd10
    x11 = xij_1l_1  *(1.-xd11) + xi_1j_1l_1 *xd11

    !print*,'--xd00,xd01,xd10,xd11 a',xd00,xd01,xd10,xd11
    !print*,'-- x00,x01,x10,x11',x00,x01,x10,x11

    y00 = yijk     *(1.-xd00) + yi_1jk     *xd00
    y01 = yijk_1   *(1.-xd01) + yi_1jk_1   *xd01
    y10 = yij_1l   *(1.-xd10) + yi_1j_1l   *xd10
    y11 = yij_1l_1 *(1.-xd11) + yi_1j_1l_1 *xd11


    z00 = zijk     *(1.-xd00) + zi_1jk     *xd00
    z01 = zijk_1   *(1.-xd01) + zi_1jk_1   *xd01
    z10 = zij_1l   *(1.-xd10) + zi_1j_1l   *xd10
    z11 = zij_1l_1 *(1.-xd11) + zi_1j_1l_1 *xd11
    
    !interp along y

    yd0 = (gitm_y - y00)/(y10 - y00)
    yd1 = (gitm_y - y01)/(y11 - y01)

    !data0 = data00*(1.-yd0)+data10*yd0
    !data1 = data01*(1.-yd1)+data11*yd1

    x0 = x00*(1.-yd0)+x10*yd0
    x1 = x01*(1.-yd1)+x11*yd1

    y0 = y00*(1.-yd0)+y10*yd0
    y1 = y01*(1.-yd1)+y11*yd1

    z0 = z00*(1.-yd0)+z10*yd0
    z1 = z01*(1.-yd1)+z11*yd1

    !interp along z
    
    zd = (gitm_z - z0)/(z1-z0)
    !data_inp = data0*(1.-zd)+data1*zd

    x_inp = x0*(1.-zd) + x1*zd
    y_inp = y0*(1.-zd) + y1*zd
    z_inp = z0*(1.-zd) + z1*zd
    
   !print*,'-- inp point',x_inp,y_inp,z_inp,iproc
    
    !print*,'end get_factors_from_8points',iproc

end subroutine get_factors_from_8points_sph

subroutine get_factors_from_8points(gitm_x,gitm_y,gitm_z,xd00,xd01,xd10,xd11,&
        yd0,yd1,zd)

    use ModGITM, only: iproc

    implicit none

    real :: gitm_x,gitm_y,gitm_z
    !integer :: i0,j0,k0,l0,i1,j1,k1,l1
    real  :: xijk     ,yijk     ,zijk, &
             xijk_1   ,yijk_1   ,zijk_1,&
             xij_1l   ,yij_1l   ,zij_1l , &
             xij_1l_1 ,yij_1l_1 ,zij_1l_1, &
             xi_1jk     ,yi_1jk     ,zi_1jk, &
             xi_1jk_1   ,yi_1jk_1   ,zi_1jk_1, &
             xi_1j_1l   ,yi_1j_1l   ,zi_1j_1l, &
             xi_1j_1l_1 ,yi_1j_1l_1 ,zi_1j_1l_1
    
    real :: xd00,xd01,xd10,xd11
    real :: x00,x01,x10,x11
    real :: y00,y01,y10,y11
    real :: z00,z01,z10,z11
    real :: data00,data01,data10,data11
    
    real :: yd0,yd1
    real :: x0,x1
    real :: y0,y1
    real :: z0,z1
    real :: data0,data1

    real :: zd
    real :: x_inp,y_inp,z_inp,data_inp
   
   ! real  :: baltst_g(nz,nf,nlt),blonst_g(nz,nf,nlt),blatst_g(nz,nf,nlt)
   !
   ! blonst_g = SAMIVars_g(1,:,:,:)
   ! blatst_g = SAMIVars_g(2,:,:,:)
   ! baltst_g = SAMIVars_g(3,:,:,:)

    call sph_to_car(blonst_g(k0,j0,i0),blatst_g(k0,j0,i0),baltst_g(k0,j0,i0),&
                    xijk     ,yijk     ,zijk)
    call sph_to_car(blonst_g(k1,j0,i0),blatst_g(k1,j0,i0),baltst_g(k1,j0,i0),&
                    xijk_1   ,yijk_1   ,zijk_1)
    call sph_to_car(blonst_g(l0,j1,i0),blatst_g(l0,j1,i0),baltst_g(l0,j1,i0),&
                    xij_1l   ,yij_1l   ,zij_1l)
    call sph_to_car(blonst_g(l1,j1,i0),blatst_g(l1,j1,i0),baltst_g(l1,j1,i0),&
                    xij_1l_1 ,yij_1l_1 ,zij_1l_1)
    
    call sph_to_car(blonst_g(k0,j0,i1),blatst_g(k0,j0,i1),baltst_g(k0,j0,i1),&
                    xi_1jk     ,yi_1jk     ,zi_1jk)
    call sph_to_car(blonst_g(k1,j0,i1),blatst_g(k1,j0,i1),baltst_g(k1,j0,i1), &
                    xi_1jk_1   ,yi_1jk_1   ,zi_1jk_1 )
    call sph_to_car(blonst_g(l0,j1,i1),blatst_g(l0,j1,i1),baltst_g(l0,j1,i1), &
                    xi_1j_1l   ,yi_1j_1l   ,zi_1j_1l)
    call sph_to_car(blonst_g(l1,j1,i1),blatst_g(l1,j1,i1),baltst_g(l1,j1,i1), &
                    xi_1j_1l_1 ,yi_1j_1l_1 ,zi_1j_1l_1)
  
    !print*,'--- gitm_xyz',gitm_x,gitm_y,gitm_z
    !print*,'1',blonst_g(k0,j0,i0),xijk 
    !print*,'2',blonst_g(k1,j0,i0),xijk_1
    !print*,'3',blonst_g(l0,j1,i0),xij_1l 
    !print*,'4',blonst_g(l1,j1,i0),xij_1l_1
    !print*,'5',blonst_g(k0,j0,i1),xi_1jk 
    !print*,'6',blonst_g(k1,j0,i1),xi_1jk_1 
    !print*,'7',blonst_g(l0,j1,i1),xi_1j_1l 
    !print*,'8',blonst_g(l1,j1,i1),xi_1j_1l_1 


    !print*,'--> a sph_to_car',iproc
    !interp along x
    xd00 = (gitm_x-xijk)    /(xi_1jk    -xijk)
    xd01 = (gitm_x-xijk_1)  /(xi_1jk_1  -xijk_1)
    xd10 = (gitm_x-xij_1l)  /(xi_1j_1l  -xij_1l)
    xd11 = (gitm_x-xij_1l_1)/(xi_1j_1l_1-xij_1l_1)


    !data00 = data[i0,j0,k0] *(1.-xd00)+data[i1,j0,k0] *xd00
    !data01 = data[i0,j0,k1] *(1.-xd01)+data[i1,j0,k1] *xd01
    !data10 = data[i0,j1,l0] *(1.-xd10)+data[i1,j1,l0] *xd10
    !data11 = data[i0,j1,l1] *(1.-xd11)+data[i1,j1,l1] *xd11

    
    x00 = xijk      *(1.-xd00) + xi_1jk     *xd00
    x01 = xijk_1    *(1.-xd01) + xi_1jk_1   *xd01
    x10 = xij_1l    *(1.-xd10) + xi_1j_1l   *xd10
    x11 = xij_1l_1  *(1.-xd11) + xi_1j_1l_1 *xd11

    !print*,'--xd00,xd01,xd10,xd11 a',xd00,xd01,xd10,xd11
    !print*,'-- x00,x01,x10,x11',x00,x01,x10,x11

    y00 = yijk     *(1.-xd00) + yi_1jk     *xd00
    y01 = yijk_1   *(1.-xd01) + yi_1jk_1   *xd01
    y10 = yij_1l   *(1.-xd10) + yi_1j_1l   *xd10
    y11 = yij_1l_1 *(1.-xd11) + yi_1j_1l_1 *xd11


    z00 = zijk     *(1.-xd00) + zi_1jk     *xd00
    z01 = zijk_1   *(1.-xd01) + zi_1jk_1   *xd01
    z10 = zij_1l   *(1.-xd10) + zi_1j_1l   *xd10
    z11 = zij_1l_1 *(1.-xd11) + zi_1j_1l_1 *xd11
    
    !interp along y

    yd0 = (gitm_y - y00)/(y10 - y00)
    yd1 = (gitm_y - y01)/(y11 - y01)

    !data0 = data00*(1.-yd0)+data10*yd0
    !data1 = data01*(1.-yd1)+data11*yd1

    x0 = x00*(1.-yd0)+x10*yd0
    x1 = x01*(1.-yd1)+x11*yd1

    y0 = y00*(1.-yd0)+y10*yd0
    y1 = y01*(1.-yd1)+y11*yd1

    z0 = z00*(1.-yd0)+z10*yd0
    z1 = z01*(1.-yd1)+z11*yd1

    !interp along z
    
    zd = (gitm_z - z0)/(z1-z0)
    !data_inp = data0*(1.-zd)+data1*zd

    x_inp = x0*(1.-zd) + x1*zd
    y_inp = y0*(1.-zd) + y1*zd
    z_inp = z0*(1.-zd) + z1*zd
    
   !print*,'-- inp point',x_inp,y_inp,z_inp,iproc
    
    !print*,'end get_factors_from_8points',iproc

end subroutine get_factors_from_8points

subroutine sph_to_car(dlon,dlat,r,x,y,z)
    
    !use ModConstants,only: pi

    implicit none

    real,intent(in) :: dlon,dlat,r  !dgree,dgree,km
    real,intent(out) :: x,y,z       ! km
    real :: consPi
    real :: rlon,rcolat
    
    consPi = 3.141592653589793 / 180.0

    rlon =  dlon*consPi
    rcolat = (90.0-dlat)*consPi

    x = r*sin(rcolat)*cos(rlon)
    y = r*sin(rcolat)*sin(rlon)
    z = r*cos(rcolat)

    return
end subroutine sph_to_car

! find AltIndex
subroutine FindAltIndex_phi(i,j,k,phi_gitm_alt)
    
    implicit none
    
    integer,intent(in) :: i, j
    integer :: k
    
    real, allocatable :: alt1(:)
    integer :: iAlt
    integer :: mm
    real :: phi_gitm_alt
   

    allocate(alt1(1:nz))


    do iAlt = 1, nz

       alt1(iAlt) = baltst_g(iAlt,j,i)-re
       !alt1(iAlt) = baltst_g(iAlt,j,i)

    enddo

    !print*, maxval(alt1),minval(alt1)
 
  if( phi_gitm_alt < alt1(1)) then
      k = -1
      !print*, '--GITM point outside SAMI3 grid, smaller then alt1(1)'
  else if (  phi_gitm_alt  > maxval(alt1)) then
      k = -1
      !print*, '--GITM point outside SAMI3 grid, larger then maxval(alt1)' 
     !     gitm_mlon,gitm_mlat,gitm_alt
      !print*,alt1
  else
      if (gitm_mlat < 0) then
         mm = 1
         do while ( phi_gitm_alt  >= alt1(mm))
            mm = mm+1
         enddo
         k = mm-1
     else
         mm = nz
         do while( phi_gitm_alt  >= alt1(mm))
            mm = mm-1
         enddo
         k = mm+1
     endif
  endif
 
  !print*,'--k 100km',k
  !if (gitm_mlat < 0) then
  !    print*, '=alt1/k,k+1',alt1(k),alt1(k+1)
  !else
  !    print*, '=alt1/k,k-1',alt1(k),alt1(k-1)
  !endif

  deallocate(alt1)

endsubroutine FindAltIndex_phi



subroutine FindAltIndex(i,j,k)
    
    implicit none
    
    integer,intent(in) :: i, j
    integer :: k
    
    real, allocatable :: alt1(:)
    integer :: iAlt
    integer :: mm
    !real  :: baltst_g(nz,nf,nlt)
   
    !baltst_g = SAMIVars_g(3,:,:,:)

    allocate(alt1(1:nz))


    do iAlt = 1, nz

       alt1(iAlt) = baltst_g(iAlt,j,i)-re
       !alt1(iAlt) = baltst_g(iAlt,j,i)

    enddo

    !print*, maxval(alt1),minval(alt1)
 
  if(gitm_alt < alt1(1)) then
      k = -1
      !print*, '--GITM point outside SAMI3 grid, smaller then alt1(1)'
  else if (gitm_alt > maxval(alt1)) then
      k = -1
      !print*, '--GITM point outside SAMI3 grid, larger then maxval(alt1)' 
     !     gitm_mlon,gitm_mlat,gitm_alt
      !print*,alt1
  else
      if (gitm_mlat < 0) then
         mm = 1
         do while (gitm_alt >= alt1(mm))
            mm = mm+1
         enddo
         k = mm-1
     else
         mm = nz
         do while(gitm_alt >= alt1(mm))
            mm = mm-1
         enddo
         k = mm+1
     endif
  endif
 
  
  !print*,'--k',k
  !if (gitm_mlat < 0) then
  !    print*, '=alt1/k,k+1',alt1(k),alt1(k+1)
  !else
  !    print*, '=alt1/k,k-1',alt1(k),alt1(k-1)
  !endif

  deallocate(alt1)

endsubroutine FindAltIndex

subroutine get_index_4points_2(phi_gitm_alt)
   
   implicit none
  
  integer :: i,j,k,l
  logical :: IsFirstTime = .true.
  integer :: inlt,inf
  integer :: nn,nn2,inn,sizen
  real,intent(in) :: phi_gitm_alt
  integer :: mm


  phi_flag = 0

  if (IsFirstTime) then
     
      !allocate(mlon0(1:nlt))
      !allocate(mlat0(1:nf))

      do inlt = 1, nlt
         mlon0(inlt) = blonst_g(1,1,inlt)
      enddo

      do inf = 1, nf
         mlat0(inf)= blatst_g(nz,inf,1)
      enddo

      IsFirstTime = .false.
  endif
  
  k = -1
  l = -1

  ! --find lon_index : i
    if ((phi_gitm_mlon <= mlon0(1)).or. &
        (phi_gitm_mlon >= mlon0(nlt))) then
        !i = -1
        phi_i0 = 1
        phi_i1 = nlt
        !print*, '--GITM point outside SAMI3 grid in LONGITUDE'
    else
        mm = 1
        do while (phi_gitm_mlon > mlon0(mm))
           mm = mm+1
        enddo
        i = mm-1
        
        phi_i0 = i
        phi_i1 = i+1

         
        !print*, '=blon/i,i+1',i,mlon0(i),mlon0(i+1)
    endif
    !print*, '=blon/i,i+1',phi_i0,phi_i1,mlon0(phi_i0),mlon0(phi_i1)


  ! --find field_line_index : j
    if (abs(phi_gitm_mlat) <= mlat0(1)) then
        !j = -1
        phi_j0 = 1
        phi_j1 = 1
        !print*, '--GITM point outside SAMI3 grid in LATITUDE'
    elseif (abs(phi_gitm_mlat) >= mlat0(nf)) then
        !j = -1
        phi_j0 = -1
        phi_j1 = -1
        !print*, '--GITM point outside SAMI3 grid in LATITUDE'
    else
        mm = 1
        do while (abs(phi_gitm_mlat) > mlat0(mm))
           mm = mm+1
        enddo
        j = mm-1
        phi_j0 = j
        phi_j1 = j+1
        !print*, '=blat/j,j+1',j,mlat0(j),mlat0(j+1)
    endif

    !if (phi_j0 > 0) &
    !      print*, '=blat/j,j+1',phi_j0,phi_j1,mlat0(phi_j0),mlat0(phi_j1)
  ! get potential index


  if ((phi_j0>0).and.(phi_j0 /= phi_j1)) then
     
      call FindAltIndex_phi(phi_i0,phi_j0,k,phi_gitm_alt)
      call FindAltIndex_phi(phi_i0,phi_j1,l,phi_gitm_alt)

      if ((k>0).and.(l>0) ) then
          phi_flag = 1
      endif

  else if ((phi_j0 == 1).and.(phi_j1 == 1)) then
      
     call FindAltIndex_phi(phi_i0,phi_j0,k,phi_gitm_alt)
      if ((k>0) ) then
          phi_flag = 1
      endif
  
  endif

  if ((phi_j0>0).and.(phi_flag == 1)) then
  
      phi_sami_lon0 = mlon0(phi_i0)
      phi_sami_lon1 = mlon0(phi_i1)

      if (phi_gitm_mlat > 0) then
          phi_sami_lat0 = blatst_g(nz,phi_j0,1)
          phi_sami_lat1 = blatst_g(nz,phi_j1,1)
      else
          phi_sami_lat0 = blatst_g(1,phi_j0,1)
          phi_sami_lat1 = blatst_g(1,phi_j1,1)
      endif

      !print*,phi_flag,phi_sami_lon0,phi_sami_lat0,phi_sami_lon1,phi_sami_lat1
  endif

end subroutine get_index_4points_2

subroutine get_index_8points_3
   
  implicit none
  
  integer :: i,j,k,l

  integer, allocatable :: AltIndex0(:),LatIndex0(:)
  real, allocatable :: LatsSel0(:)
  integer :: ifl,mm,ifl_sel
  logical :: IsFirstTime = .true.
  integer :: inlt,inf
  integer :: nn,nn2,inn,sizen
  integer :: flagk,flagl
  integer :: jtmp

  flag_eq = 0

  if (IsFirstTime) then
     
      allocate(mlon0(1:nlt))
      allocate(mlat0(1:nf))

      do inlt = 1, nlt
         mlon0(inlt) = blonst_g(1,1,inlt)
      enddo

      do inf = 1, nf
         mlat0(inf)= blatst_g(nz,inf,1)
      enddo

      !print*,'-->> sami mlon0'
      !print*,mlon0
      !print*,'-->> sami mlat0'
      !
      !print*,'>>>>',blatst_g(:,nf,1)

      IsFirstTime = .false.
  endif
  
  k = -1
  l = -1

  ! --find lon_index : i
    if ((gitm_mlon <= mlon0(1)).or.(gitm_mlon >= mlon0(nlt))) then
        !i = -1
        i0 = 1
        i1 = nlt
        !print*, '--GITM point outside SAMI3 grid in LONGITUDE',gitm_mlon
    else
        mm = 1
        do while (gitm_mlon > mlon0(mm))
           mm = mm+1
        enddo
        i = mm-1
        
        i0 = i
        i1 = i+1
         
        !print*, '=blon/i,i+1',i,mlon0(i),mlon0(i+1)
    endif

  ! --find field_line_index : j
    if (abs(gitm_mlat) <= mlat0(1)) then
        !j = -1
        j0 = -1
        j1 = -1
        !print*, '--GITM point outside SAMI3 grid in LATITUDE'
    elseif (abs(gitm_mlat) >= mlat0(nf)) then
        !j = -1
        j0 = -1
        j1 = -1
        !print*, '--GITM point outside SAMI3 grid in LATITUDE'
    else
        mm = 1
        do while (abs(gitm_mlat) > mlat0(mm))
           mm = mm+1
        enddo
        j = mm-1
        j0 = j
        j1 = j+1
        !print*, '=blat/j,j+1',j,mlat0(j),mlat0(j+1)
    endif

    !if (j0>0) &
    !print*, '=blat/j,j+1',j0,mlat0(j0),mlat0(j1)
  
    if (j0>0) then

          call FindAltIndex(i0,j0,k)
          call FindAltIndex(i0,j1,l)
      
              
          !print*,'Get 8points_3',k,l,j0,j1
          if (k<0) then
             
              ! print*,'EQ points filling',i0,i1,j0,j1,k,l
              call process2(i0,j1,jtmp,k,l)
              
              j0 = jtmp
              j1 = j0+1

              flag_eq = 1
          endif

    endif


  if (gitm_mlat>0) then
     k1 = k-1
     l1 = l-1
  else
     k1 = k+1
     l1 = l+1
  endif

  k0 = k
  l0 = l


  if ((j0>0).and.(k0>0).and.(l0>0)) then
  
      sami_lon0 = mlon0(i0)
      sami_lon1 = mlon0(i1)

      if (gitm_mlat > 0) then
          sami_lat0 = blatst_g(nz,j0,1)
          sami_lat1 = blatst_g(nz,j1,1)
      else
          sami_lat0 = blatst_g(1,j0,1)
          sami_lat1 = blatst_g(1,j1,1)
      endif
  !print*, '--',sami_lon0,sami_lon1,sami_lat0,sami_lat1
  !print*,i0,j0,k0,l0,i1,j1,k1,l1
  endif

  !print*, '--',sami_lat0,sami_lat1
  
  !if ((i0>-1) .and. (j0>-1)) &
  !print*,i0,j0,k0,l0,i1,j1,k1,l1

end subroutine get_index_8points_3

subroutine get_index_8points_2
   
   implicit none

  !use parameter_mod, only: nz,nf,nlt,re
  
  !real :: gitm_mlat,gitm_mlon,gitm_alt
  integer :: i,j,k,l
  !integer :: i0,j0,k0,l0,i1,j1,k1,l1

  integer, allocatable :: AltIndex(:),LatIndex(:)
  real, allocatable :: LatsSel(:)
  integer, allocatable :: AltIndex0(:),LatIndex0(:)
  real, allocatable :: LatsSel0(:)
  integer :: ifl,ktmp,mm,ifl_sel
  logical :: IsFirstTime = .true.
  integer :: inlt,inf
  integer :: nn,nn2,inn,sizen
  !real :: sami_lat0,sami_lat1

  !real  :: baltst_g(nz,nf,nlt),blonst_g(nz,nf,nlt),blatst_g(nz,nf,nlt)
  ! 
  !  blonst_g = SAMIVars_g(1,:,:,:)
  !  blatst_g = SAMIVars_g(2,:,:,:)
  !  baltst_g = SAMIVars_g(3,:,:,:)

  if (IsFirstTime) then
     
      allocate(mlon0(1:nlt))
      allocate(mlat0(1:nf))

      do inlt = 1, nlt
         mlon0(inlt) = blonst_g(1,1,inlt)
      enddo

      do inf = 1, nf
         mlat0(inf)= blatst_g(nz,inf,1)
      enddo

      !print*,'-->> sami mlon0'
      !print*,mlon0
      !print*,'-->> sami mlat0'
      !
      !print*,'>>>>',blatst_g(:,nf,1)

      IsFirstTime = .false.
  endif
  
  k = -1
  l = -1

  ! --find lon_index : i
    if ((gitm_mlon <= mlon0(1)).or.(gitm_mlon >= mlon0(nlt))) then
        !i = -1
        i0 = 1
        i1 = nlt
        !print*, '--GITM point outside SAMI3 grid in LONGITUDE'
    else
        mm = 1
        do while (gitm_mlon > mlon0(mm))
           mm = mm+1
        enddo
        i = mm-1
        
        i0 = i
        i1 = i+1

         
        !print*, '=blon/i,i+1',i,mlon0(i),mlon0(i+1)
    endif
    !print*, '=blon/i,i+1',i0,mlon0(i0),mlon0(i1)


  ! --find field_line_index : j
    if (abs(gitm_mlat) <= mlat0(1)) then
        !j = -1
        j0 = -1
        j1 = -1
        !print*, '--GITM point outside SAMI3 grid in LATITUDE'
    elseif (abs(gitm_mlat) >= mlat0(nf)) then
        !j = -1
        j0 = -1
        j1 = -1
        !print*, '--GITM point outside SAMI3 grid in LATITUDE'
    else
        mm = 1
        do while (abs(gitm_mlat) > mlat0(mm))
           mm = mm+1
        enddo
        j = mm-1
        j0 = j
        j1 = j+1
        !print*, '=blat/j,j+1',j,mlat0(j),mlat0(j+1)
    endif

    !if (j0>0) &
    !print*, '=blat/j,j+1',j0,mlat0(j0),mlat0(j1)
  
  if ((i0>0).and.(j0>0).and.(j0 /= j1)) then
     
      call FindAltIndex(i0,j0,k)
      call FindAltIndex(i0,j1,l)
  

  else if ((j0 == 1).and.(j1 == 1)) then
      
     call FindAltIndex(i0,j0,k)
     if (k>0) l = nz+1-k
  
  endif

  if (gitm_mlat>0) then
     k1 = k-1
     l1 = l-1
  else
     k1 = k+1
     l1 = l+1
  endif

  k0 = k
  l0 = l

  !! get potential index

  !!phi_i0 = -1
  !!phi_i1 = -1
  !!phi_j0 = -1
  !!phi_j1 = -1
  !phi_flag = 0

  !if ((i0>0).and.(j0>0).and.(j0 /= j1)) then
  !   
  !    call FindAltIndex_100km(i0,j0,k)
  !    call FindAltIndex_100km(i0,j1,l)

  !    if ((k>0).and.(l>0) ) then
  !        phi_flag = 1
  !    endif

  !else if ((j0 == 1).and.(j1 == 1)) then
  !    
  !   call FindAltIndex_100km(i0,j0,k)
  !    if ((k>0) ) then
  !        phi_flag = 1
  !    endif
  !
  !endif



  !!print*,'i/j',i,j

  !if ((i>0).and.(j>0)) then
  !    ! --find alt index: k,l

  !   allocate(AltIndex0(1:nf-j+1))
  !   allocate(LatIndex0(1:nf-j+1))
  !   allocate(LatsSel0(1:nf-j+1))

  !   nn = 0
  !   do ifl = j,nf
  !  
  !      call FindAltIndex(i,ifl,ktmp)
  !   
  !      !print*,'ifl',ifl,ktmp
  !      
  !      AltIndex0(ifl-j+1) = ktmp
  !      LatsSel0(ifl-j+1) = blatst_g(ktmp,ifl,i)
  !      LatIndex0(ifl-j+1) = ifl

  !      if (ktmp>0) nn = nn+1

  !   enddo

  !   allocate(AltIndex(1:nn))
  !   allocate(LatIndex(1:nn))
  !   allocate(LatsSel(1:nn))

  !  ! Block0: where( AltIndex0 /= -1)
  !  !     AltIndex = AltIndex0
  !  !     LatIndex = LatIndex0
  !  !     LatsSel = LatsSel0
  !  ! end  where Block0
  !  
  !  sizen = size(AltIndex0)
  !  nn2 = 1
  !  do inn = 1,sizen
  !      if (AltIndex0(inn)>0) then
  !          AltIndex(nn2) = AltIndex0(inn)
  !          LatIndex(nn2) = LatIndex0(inn)
  !          LatsSel(nn2) = LatsSel0(inn)
  !          nn2= nn2+1
  !      endif
  !  enddo



  !   !print*,LatsSel
  !   !print*,LatIndex
  !   !print*, AltIndex
  !   !print*,'--',nn,size(AltIndex0),size(AltIndex)
  !  
  !   if ((abs(gitm_mlat)<minval(abs(LatsSel))).or. &
  !       (abs(gitm_mlat)>maxval(abs(LatsSel)))) then

  !     ! print*,'beyond LATITUDE RANGE'
  !      j = -1
  !      k = -1
  !      l = -1
  !   else
  !      mm=1
  !      do while(abs(gitm_mlat)>=abs(LatsSel(mm)))
  !         mm = mm+1
  !      enddo
 
  !      ifl_sel = mm-1

  !      !print*,'ifl_sel',ifl_sel

  !      j = LatIndex(ifl_sel)
  !      k = AltIndex(ifl_sel)
  !      l = AltIndex(ifl_sel+1)
  !   endif
  !  
  !   deallocate(AltIndex,LatIndex,LatsSel)
  !   
  !   !print*,'j b',j
  !endif

  !if (gitm_mlat>0) then
  !   k1 = k-1
  !   l1 = l-1
  !else
  !   k1 = k+1
  !   l1 = l+1
  !endif

  !i1 = i+1
  !j1 = j+1

  !i0 = i
  !j0 = j
  !k0 = k
  !l0 = l

  if ((j0>0).and.(k0>0).and.(l0>0)) then
  
      sami_lon0 = mlon0(i0)
      sami_lon1 = mlon0(i1)

      if (gitm_mlat > 0) then
          sami_lat0 = blatst_g(nz,j0,1)
          sami_lat1 = blatst_g(nz,j1,1)
      else
          sami_lat0 = blatst_g(1,j0,1)
          sami_lat1 = blatst_g(1,j1,1)
      endif
  !print*, '--',sami_lon0,sami_lon1,sami_lat0,sami_lat1
  !print*,i0,j0,k0,l0,i1,j1,k1,l1
  endif

  !print*, '--',sami_lat0,sami_lat1
  
  !if ((i0>-1) .and. (j0>-1)) &
  !print*,i0,j0,k0,l0,i1,j1,k1,l1

end subroutine get_index_8points_2

subroutine get_index_8points
   
   implicit none

  !use parameter_mod, only: nz,nf,nlt,re
  
  !real :: gitm_mlat,gitm_mlon,gitm_alt
  integer :: i,j,k,l
  !integer :: i0,j0,k0,l0,i1,j1,k1,l1

  integer, allocatable :: AltIndex(:),LatIndex(:)
  real, allocatable :: LatsSel(:)
  integer, allocatable :: AltIndex0(:),LatIndex0(:)
  real, allocatable :: LatsSel0(:)
  integer :: ifl,ktmp,mm,ifl_sel
  logical :: IsFirstTime = .true.
  integer :: inlt,inf
  integer :: nn,nn2,inn,sizen
  !real  :: baltst_g(nz,nf,nlt),blonst_g(nz,nf,nlt),blatst_g(nz,nf,nlt)
  ! 
  !  blonst_g = SAMIVars_g(1,:,:,:)
  !  blatst_g = SAMIVars_g(2,:,:,:)
  !  baltst_g = SAMIVars_g(3,:,:,:)

  if (IsFirstTime) then
     
      allocate(mlon0(1:nlt))
      allocate(mlat0(1:nf))

      do inlt = 1, nlt
         mlon0(inlt) = blonst_g(1,1,inlt)
      enddo

      do inf = 1, nf
         mlat0(inf)= blatst_g(nz,inf,1)
      enddo

      !print*,'-->> sami mlon0'
      !print*,mlon0
      !print*,'-->> sami mlat0'
      !
      !print*,'>>>>',blatst_g(:,nf,1)

      IsFirstTime = .false.
  endif
  
  k = -1
  l = -1

  ! --find lon_index : i
    if ((gitm_mlon < mlon0(1)).or.(gitm_mlon > mlon0(nlt))) then
        i = -1
        !print*, '--GITM point outside SAMI3 grid in LONGITUDE'
    else
        mm = 1
        do while (gitm_mlon >= mlon0(mm))
           mm = mm+1
        enddo
        i = mm-1
        !print*, '=blon/i,i+1',i,mlon0(i),mlon0(i+1)
    endif


  ! --find field_line_index : j
    if ((abs(gitm_mlat) < mlat0(1)).or.(abs(gitm_mlat) > mlat0(nf))) then
        j = -1
        !print*, '--GITM point outside SAMI3 grid in LATITUDE'
    else
        mm = 1
        do while (abs(gitm_mlat) >= mlat0(mm))
           mm = mm+1
        enddo
        j = mm-1
        !print*, '=blat/j,j+1',j,mlat0(j),mlat0(j+1)
    endif


  !print*,'i/j',i,j

  if ((i>0).and.(j>0)) then
      ! --find alt index: k,l

     allocate(AltIndex0(1:nf-j+1))
     allocate(LatIndex0(1:nf-j+1))
     allocate(LatsSel0(1:nf-j+1))

     nn = 0
     do ifl = j,nf
    
        call FindAltIndex(i,ifl,ktmp)
     
        !print*,'ifl',ifl,ktmp
        
        AltIndex0(ifl-j+1) = ktmp
        LatsSel0(ifl-j+1) = blatst_g(ktmp,ifl,i)
        LatIndex0(ifl-j+1) = ifl

        if (ktmp>0) nn = nn+1

     enddo

     allocate(AltIndex(1:nn))
     allocate(LatIndex(1:nn))
     allocate(LatsSel(1:nn))

    ! Block0: where( AltIndex0 /= -1)
    !     AltIndex = AltIndex0
    !     LatIndex = LatIndex0
    !     LatsSel = LatsSel0
    ! end  where Block0
    
    sizen = size(AltIndex0)
    nn2 = 1
    do inn = 1,sizen
        if (AltIndex0(inn)>0) then
            AltIndex(nn2) = AltIndex0(inn)
            LatIndex(nn2) = LatIndex0(inn)
            LatsSel(nn2) = LatsSel0(inn)
            nn2= nn2+1
        endif
    enddo



     !print*,LatsSel
     !print*,LatIndex
     !print*, AltIndex
     !print*,'--',nn,size(AltIndex0),size(AltIndex)
    
     if ((abs(gitm_mlat)<minval(abs(LatsSel))).or. &
         (abs(gitm_mlat)>maxval(abs(LatsSel)))) then

       ! print*,'beyond LATITUDE RANGE'
        j = -1
        k = -1
        l = -1
     else
        mm=1
        do while(abs(gitm_mlat)>=abs(LatsSel(mm)))
           mm = mm+1
        enddo
 
        ifl_sel = mm-1

        !print*,'ifl_sel',ifl_sel

        j = LatIndex(ifl_sel)
        k = AltIndex(ifl_sel)
        l = AltIndex(ifl_sel+1)
     endif
    
     deallocate(AltIndex,LatIndex,LatsSel)
     
     !print*,'j b',j
  endif

  if (gitm_mlat>0) then
     k1 = k-1
     l1 = l-1
  else
     k1 = k+1
     l1 = l+1
  endif

  i1 = i+1
  j1 = j+1

  i0 = i
  j0 = j
  k0 = k
  l0 = l
  
  !if ((i0>-1) .and. (j0>-1)) &
  !print*,i0,j0,k0,l0,i1,j1,k1,l1

end subroutine get_index_8points

! find indexes of 8 points in SAMI3
subroutine get_index_8points1
   
   implicit none

  !use parameter_mod, only: nz,nf,nlt,re
  
  !real :: gitm_mlat,gitm_mlon,gitm_alt
  integer :: i,j,k,l
  !integer :: i0,j0,k0,l0,i1,j1,k1,l1

  integer, allocatable :: AltIndex(:),LatIndex(:)
  real, allocatable :: LatsSel(:)
  integer :: ifl,ktmp,mm,ifl_sel
  !real  :: baltst_g(nz,nf,nlt),blonst_g(nz,nf,nlt),blatst_g(nz,nf,nlt)
  ! 
  !  blonst_g = SAMIVars_g(1,:,:,:)
  !  blatst_g = SAMIVars_g(2,:,:,:)
  !  baltst_g = SAMIVars_g(3,:,:,:)
  
  ! --find lon_index : i
    
  mm = 2
  do while (gitm_mlon >= mlon0(mm))
     mm = mm+1
  enddo
     i = mm-1
  !print*, '=blon/i,i+1',mlon0(i),mlon0(i+1)

  
  ! --find field_line_index : j

  mm = 2
  do while (abs(gitm_mlat) >= mlat0(mm))
     mm = mm+1
  enddo
  j = mm-1

  !print*, '=blat/j,j+1',j,mlat0(j),mlat0(j+1)

  ! --find alt index: k,l

  allocate(AltIndex(1:nf-j+1))
  allocate(LatIndex(1:nf-j+1))
  allocate(LatsSel(1:nf-j+1))

  do ifl = j,nf
    
     call FindAltIndex(i,ifl,ktmp)
     
     !print*,'ifl',ifl,ktmp
     AltIndex(ifl-j+1) = ktmp
     LatsSel(ifl-j+1) = blatst_g(ktmp,ifl,i)
     LatIndex(ifl-j+1) = ifl

  enddo

  ! print*,LatsSel
  
  mm=1
  do while(abs(gitm_mlat)>=abs(LatsSel(mm)))
     mm = mm+1
  enddo
 
  ifl_sel = mm-1

  j = LatIndex(ifl_sel)
  k = AltIndex(ifl_sel)
  l = AltIndex(ifl_sel+1)

  if (gitm_mlat>0) then
     k1 = k-1
     l1 = l-1
  else
     k1 = k+1
     l1 = l+1
  endif

  i1 = i+1
  j1 = j+1

  i0 = i
  j0 = j
  k0 = k
  l0 = l
   
  deallocate(AltIndex,LatIndex,LatsSel)
  return
end subroutine get_index_8points1

end module ModSamiInterp
