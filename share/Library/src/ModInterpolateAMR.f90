!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModInterpolateSimpleShape
  !\
  ! Construct interpolation weights for a point inside
  ! bodies, interolation being in terms of the values 
  ! in the vertexes of the body. Examples of the implemented
  ! bodies: tetrahedron, pyramids (rectangular and trapezoidal)
  ! and their different compositions. Interpolation in parallel 
  ! rays implements the following procedure; a ray of given 
  ! direction is passed trough a given point inside the body, 
  ! in the point of the ray intersects the body boundary, the value
  ! is interpolated over planar subface (traingle, trapezoid, 
  ! ractangle), then the point value is weighted in accordance with
  ! the distances from the given point to its intersections with the 
  ! body surface. Resolution_corner interpolates within the resolution 
  ! corner, which has "fine vertexes"  - grid points of finer 
  ! rectangular grid and "coarse vertexes" - grid points of coarser
  ! rectangular grid with twice larger grid size. Finer grid cube and 
  ! coarser grid cube have common central point. Such grid coniguration
  ! is typical for the resolution corner of a block adaptive grid, 
  ! so that the provided routine may be used for interpolating values
  ! in the proximity of such points.  
  !/
  implicit none
  PRIVATE
  SAVE
  integer, parameter:: nDim = 3, nGrid = 8
  integer, parameter:: Rectangular_=1, Trapezoidal_=2
  public :: interpolate_tetrahedron
  public :: interpolate_pyramid
  public:: interpolate_pyramids, interpolate_on_parallel_rays
contains
  !========================
  function cross_product(a_D, B_d)
    real, dimension(nDim), intent(in) :: a_D, b_D
    real, dimension(nDim) :: cross_product
    integer, parameter:: x_ = 1, y_ = 2, z_ = 3
    !-----------------------------------------------
    cross_product(x_) = a_D(y_)*b_D(z_) - a_D(z_)*b_D(y_)
    cross_product(y_) = a_D(z_)*b_D(x_) - a_D(x_)*b_D(z_)
    cross_product(z_) = a_D(x_)*b_D(y_) - a_D(y_)*b_D(x_)
  end function cross_product
  !=========================
  real function triple_product(a_D, b_D, c_D)
    real, dimension(nDim), intent(in) :: a_D, b_D, c_D
    !-----------------------------------------------
    triple_product = sum(a_D * cross_product(b_D, c_D))
  end function triple_product
  !===================================================
  subroutine interpolate_tetrahedron(X1_D, X2_D, X3_D, X4_D, Xyz_D, Weight_I)
    !\
    ! Interpolate in the tetrahedron
    !/
    real, dimension(nDim), intent(in) :: X1_D, X2_D, X3_D, X4_D, Xyz_D
    real, intent(out):: Weight_I(1:4)
    real:: Aux
    !-----------------------------------------------
    !\
    ! Need to solve an equation:
    ! Xyz = Weight_I(1)*X1+Weight_I(2)*X2+Weight_I(3)*X3+Weight_I(4)*X4
    ! Or, which is equivalent:
    ! Xyz -X1 = Weight_I(2)*(X2-X1)+Weight_I(3)*(X3-X1)+Weight_I(4)*(X4-X1)
    !/
    Aux= 1/triple_product( X4_D - X1_D,  X3_D - X1_D,  X2_D - X1_D)
    Weight_I(4) = &
         triple_product(Xyz_D - X1_D,  X3_D - X1_D,  X2_D - X1_D)*Aux
    Weight_I(3) = &
         triple_product( X4_D - X1_D, Xyz_D - X1_D,  X2_D - X1_D)*Aux
    Weight_I(2) = &
         triple_product( X4_D - X1_D,  X3_D - X1_D, Xyz_D - X1_D)*Aux
    Weight_I(1) = 1 - sum(Weight_I(2:4))
  end subroutine interpolate_tetrahedron
  !=======================================================================
  subroutine interpolate_pyramid(&
       iBase,X1_D, X2_D, X3_D, X4_D, X5_D, Xyz_D, Weight_I)
    !\
    ! Interpolate in the pyramid with base X1X2X3X4 and apex X5
    ! valid for case of rectangular or trapezoid base
    !/
    integer, intent(in) :: iBase
    real, dimension(nDim), intent(in) :: X1_D, X2_D, X3_D, X4_D, X5_D, Xyz_D
    real, intent(out):: Weight_I(1:5)
    !\
    !Projection of Xyz point on the base of the pyramid X1X2X3X4 
    !along the line X5Xyz
    !/
    real, dimension(nDim) :: XyzP_D 
    real ::Alpha5
    !-----------------------------------------------
    !\
    ! Solve equation: X5 + (Xyz - X5)/Alpha5 = XyzP
    ! where XyzP belongs to the pydamid base.
    ! As long as (XyzP-X1_D)\cdot[(X2-X1)\times(X3-X1)]=0,
    ! we have Alpha5 =(&
    ! X5 -Xyz)\cdot[(X2-X1)\times(X3-X1)]/(X5-X1)\cdot[(X2-X1)\times(X3-X1)]
    !/
    Alpha5 = triple_product(X5_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
         triple_product(X5_D - X1_D , X2_D - X1_D, X3_D - X1_D)
    if(Alpha5==0.0)then
       if(all(Xyz_D(:)==X5_D(:)))then
          Weight_I(5) = 1;   Weight_I(1:4) = 0
          RETURN
       else
          Weight_I(1:5) = -1
          RETURN
       end if
    elseif(Alpha5 < 0.0 .or. Alpha5 > 1.0)then
       Weight_I(1:5) = -1
       RETURN
    end if
    XyzP_D = X5_D + (Xyz_D-X5_D)/Alpha5
    !\
    ! Now, Xyz = Alpha5*XyzP + (1-Alpha5)*X5
    ! Find weight of the apex point X5_D
    !/
    Weight_I(5) = 1 - Alpha5
    !\
    ! Find weights of the base points X1, X2, X3, X4
    !/
    select case(iBase)
    case(Rectangular_)
       call rectangle(X1_D, X2_D, X3_D, X4_D, XyzP_D, Weight_I(1:4))
    case(Trapezoidal_)
       call trapezoid(X1_D, X2_D, X3_D, X4_D, XyzP_D, Weight_I(1:4))
    end select
    if(any(Weight_I(1:4) < 0.0))RETURN
    !\
    ! Correct weights due to pyramid geometry
    !/
    Weight_I(1:4) = Weight_I(1:4) * (1 - Weight_I(5))
  end subroutine interpolate_pyramid
  !============
  subroutine rectangle(X1_D, X2_D, X3_D, X4_D, XyzP_D, Weight_I)
    real, dimension(nDim), intent(in):: X1_D, X2_D, X3_D, X4_D, XyzP_D
    real, intent(out):: Weight_I(1:4)
    real:: x,y
    !---------
    !Calculate dimensionless coordinates with respect to vertex 1 
    y =  sum((X3_D - X1_D)*(XyzP_D - X1_D))/&
         sum((X3_D - X1_D)*(X3_D - X1_D))
    x =  sum((XyzP_D - X1_D)*(X2_D - X1_D))/&
         sum((X2_D - X1_D)*(X2_D - X1_D) )
    Weight_I(3)  =     y *(1 - x) ; Weight_I(4) =      y *x
    Weight_I(1) = (1 - y)*(1 - x) ; Weight_I(2) = (1 - y)*x
  end subroutine rectangle
  !==============================
  subroutine trapezoid(X1_D, X2_D, X3_D, X4_D, XyzP_D, Weight_I)
    real, dimension(nDim), intent(in):: X1_D, X2_D, X3_D, X4_D, XyzP_D
    real, intent(out)::Weight_I(1:4)
    real:: x,y, XyzMin_D(nDim), XyzMax_D(nDim)
    !We require that the lager base is X1X2
    !---------
    !Calculate dimensionless coordinates with respect to vertex 1
    XyzMin_D = 0.50*(X1_D + X2_D);  XyzMax_D = 0.50*(X3_D + X4_D)
    y =  sum((XyzMax_D - XyzMin_D)*(XyzP_D - XyzMin_D))/&
         sum((XyzMax_D - XyzMin_D)*(XyzMax_D - XyzMin_D))
    x =  sum((XyzP_D   - X1_D)*(X2_D - X1_D))/&
         sum( (X2_D    - X1_D)*(X2_D - X1_D))
    if( y < 0.0 .or. y > 1.0 .or. x < 0.0 .or. x > 1.0)then
       Weight_I(1:4) = -1
       RETURN
    end if
    if( x <= 0.250) then
       !Interpolation in triangle 132, 4th weight is zero
       Weight_I(3) = y                   ; Weight_I(4) = 0
       Weight_I(1) = 1 - 0.750*y - x     ; Weight_I(2) = x - 0.250*y
    elseif(x <= 0.750) then
       !Bilinear interpolation
       Weight_I(3) = y*(1 - 2*(x - 0.25)); Weight_I(4) = y*2*(x - 0.25)
       Weight_I(1) = (1 - y)*(1 - x)     ; Weight_I(2) = (1 - y)*x
    else
       !Interpolation in triangle 142, 3rd weight is zero
       Weight_I(3) = 0                   ; Weight_I(4) = y
       Weight_I(1) = 1 - 0.250*y - x     ; Weight_I(2) = x - 0.750*y
    end if
  end subroutine trapezoid
  !=====================================
  subroutine triangle(X1_D, X2_D, X3_D,  XyzP_D, Weight_I)
    real, dimension(nDim), intent(in):: X1_D, X2_D, X3_D, XyzP_D
    real, intent(out):: Weight_I(1:3)
    real:: CrossProductInv_D(nDim)
    !\
    ! In 2D case
    ! Weight_I(3) = cross_product(XyzP_D-X1_D,X2_D-X1_D)/&
    !               cross_product(X3_D-X1_D,X2_D-X1_D)
    ! Weight_I(2) = cross_product(X3_D-X1_D,Xyz_D-X1_D)/&
    !               cross_product(X3_D - X1_D, X2_D-X1_D)
    ! iN 3d case instead of /cross_product(X3_D - X1_D, X2_D-X1_D) we use 
    ! \cdot cross_product(X3_D - X1_D, X2_D-X1_D)/&
    !       sum(cross_product(X3_D - X1_D, X2_D-X1_D)**2)
    !/
    CrossProductInv_D = cross_product(X3_D - X1_D, X2_D - X1_D)
    CrossProductInv_D = CrossProductInv_D/sum(CrossProductInv_D**2)

    Weight_I(3) = sum(cross_product(XyzP_D - X1_D, X2_D  - X1_D)*&
         CrossProductInv_D)
    Weight_I(2) = sum(cross_product(X3_D  - X1_D, XyzP_D - X1_D)*&
         CrossProductInv_D)
    Weight_I(1) = 1 - sum(Weight_I(2:3))
  end subroutine triangle
  !=====================================
  subroutine interpolate_pyramids(&
       XyzGrid_DI,Xyz_D,iOrder_I, Weight_I, nGridOut, &
       iTetrahedron1_I,iTetrahedron2_I,iTetrahedron3_I,&
       iTetrahedron4_I,iTetrahedron5_I,iTetrahedron6_I,&
       iRectangular1_I,iRectangular2_I,iRectangular3_I,&
       iTrapezoidal1_I,iTrapezoidal2_I)
    integer,intent(in),optional,dimension(4)::&
         iTetrahedron1_I, iTetrahedron2_I, iTetrahedron3_I,&
         iTetrahedron4_I, iTetrahedron5_I, iTetrahedron6_I
    integer,intent(in),optional,dimension(5)::&
         iRectangular1_I, iRectangular2_I, iRectangular3_I,&
         iTrapezoidal1_I, iTrapezoidal2_I
    real, intent(in):: XyzGrid_DI(nDim,nGrid), Xyz_D(nDim)
    real, intent(out):: Weight_I(nGrid)
    integer, intent(out)::nGridOut, iOrder_I(1:5)
    !---------------------------
    nGridOut = -1; Weight_I = -1
    if(present(iTetrahedron1_I))then
       call interpolate_tetrahedron(&
            XyzGrid_DI(:,iTetrahedron1_I(1)),&
            XyzGrid_DI(:,iTetrahedron1_I(2)),&
            XyzGrid_DI(:,iTetrahedron1_I(3)),&
            XyzGrid_DI(:,iTetrahedron1_I(4)),Xyz_D, Weight_I(1:4))
       if(all(Weight_I(1:4).ge.0.0))then
          nGridOut = 4
          iOrder_I(1:4) = iTetrahedron1_I
          RETURN
       elseif(present(iTetrahedron2_I))then
          call interpolate_tetrahedron(&
               XyzGrid_DI(:,iTetrahedron2_I(1)),&
               XyzGrid_DI(:,iTetrahedron2_I(2)),&
               XyzGrid_DI(:,iTetrahedron2_I(3)),&
               XyzGrid_DI(:,iTetrahedron2_I(4)),Xyz_D, Weight_I(1:4))
          if(all(Weight_I(1:4).ge.0.0))then
             nGridOut = 4
             iOrder_I(1:4) = iTetrahedron2_I
             RETURN
          elseif(present(iTetrahedron3_I))then
             call interpolate_tetrahedron(&
                  XyzGrid_DI(:,iTetrahedron3_I(1)),&
                  XyzGrid_DI(:,iTetrahedron3_I(2)),&
                  XyzGrid_DI(:,iTetrahedron3_I(3)),&
                  XyzGrid_DI(:,iTetrahedron3_I(4)),&
                  Xyz_D, Weight_I(1:4))
             if(all(Weight_I(1:4).ge.0.0))then
                nGridOut = 4
                iOrder_I(1:4) = iTetrahedron3_I
                RETURN
             elseif(present(iTetrahedron4_I))then
                call interpolate_tetrahedron(&
                     XyzGrid_DI(:,iTetrahedron4_I(1)),&
                     XyzGrid_DI(:,iTetrahedron4_I(2)),&
                     XyzGrid_DI(:,iTetrahedron4_I(3)),&
                     XyzGrid_DI(:,iTetrahedron4_I(4)),&
                     Xyz_D, Weight_I(1:4))
                if(all(Weight_I(1:4).ge.0.0))then
                   nGridOut = 4
                   iOrder_I(1:4) = iTetrahedron4_I
                   RETURN
                elseif(present(iTetrahedron5_I))then
                   call interpolate_tetrahedron(&
                        XyzGrid_DI(:,iTetrahedron5_I(1)),&
                        XyzGrid_DI(:,iTetrahedron5_I(2)),&
                        XyzGrid_DI(:,iTetrahedron5_I(3)),&
                        XyzGrid_DI(:,iTetrahedron5_I(4)),&
                        Xyz_D, Weight_I(1:4))
                   if(all(Weight_I(1:4).ge.0.0))then
                      nGridOut = 4
                      iOrder_I(1:4) = iTetrahedron5_I
                      RETURN
                   elseif(present(iTetrahedron6_I))then
                      call interpolate_tetrahedron(&
                           XyzGrid_DI(:,iTetrahedron6_I(1)),&
                           XyzGrid_DI(:,iTetrahedron6_I(2)),&
                           XyzGrid_DI(:,iTetrahedron6_I(3)),&
                           XyzGrid_DI(:,iTetrahedron6_I(4)),&
                           Xyz_D, Weight_I(1:4))
                      if(all(Weight_I(1:4).ge.0.0))then
                         nGridOut = 4
                         iOrder_I(1:4) = iTetrahedron6_I
                         RETURN
                      end if ! 6
                   end if    ! 5
                end if       ! 4
             end if          ! 3
          end if             ! 2
       end if                ! 1
    end if                   ! no tetrahedron
    if(present(iRectangular1_I))then
       call interpolate_pyramid(Rectangular_,&
            XyzGrid_DI(:,iRectangular1_I(1)),&
            XyzGrid_DI(:,iRectangular1_I(2)),&
            XyzGrid_DI(:,iRectangular1_I(3)),&
            XyzGrid_DI(:,iRectangular1_I(4)),&
            XyzGrid_DI(:,iRectangular1_I(5)),Xyz_D, Weight_I(1:5) )
       if(all(Weight_I(1:5).ge.0.0))then
          nGridOut = 5
          iOrder_I(1:5) = iRectangular1_I
          RETURN
       elseif(present(iRectangular2_I))then
          call interpolate_pyramid(Rectangular_,&
               XyzGrid_DI(:,iRectangular2_I(1)),&
               XyzGrid_DI(:,iRectangular2_I(2)),&
               XyzGrid_DI(:,iRectangular2_I(3)),&
               XyzGrid_DI(:,iRectangular2_I(4)),&
               XyzGrid_DI(:,iRectangular2_I(5)),Xyz_D, Weight_I(1:5) )
          if(all(Weight_I(1:5).ge.0.0))then
             nGridOut = 5
             iOrder_I(1:5) = iRectangular2_I
             RETURN
          elseif(present(iRectangular3_I))then
             call interpolate_pyramid(Rectangular_,&
                  XyzGrid_DI(:,iRectangular3_I(1)),&
                  XyzGrid_DI(:,iRectangular3_I(2)),&
                  XyzGrid_DI(:,iRectangular3_I(3)),&
                  XyzGrid_DI(:,iRectangular3_I(4)),&
                  XyzGrid_DI(:,iRectangular3_I(5)),Xyz_D, Weight_I(1:5) )
             if(all(Weight_I(1:5).ge.0.0))then
                nGridOut = 5
                iOrder_I(1:5) = iRectangular3_I
                RETURN
             end if          ! 3
          end if             ! 2
       end if                ! 1
    end if                   ! no rectangular
    if(present(iTrapezoidal1_I))then
       call interpolate_pyramid(Trapezoidal_,&
            XyzGrid_DI(:,iTrapezoidal1_I(1)),&
            XyzGrid_DI(:,iTrapezoidal1_I(2)),&
            XyzGrid_DI(:,iTrapezoidal1_I(3)),&
            XyzGrid_DI(:,iTrapezoidal1_I(4)),&
            XyzGrid_DI(:,iTrapezoidal1_I(5)),Xyz_D, Weight_I(1:5) )
       if(all(Weight_I(1:5).ge.0.0))then
          nGridOut = 5
          iOrder_I(1:5) = iTrapezoidal1_I
          RETURN
       elseif(present(iTrapezoidal2_I))then
          call interpolate_pyramid(Trapezoidal_,&
               XyzGrid_DI(:,iTrapezoidal2_I(1)),&
               XyzGrid_DI(:,iTrapezoidal2_I(2)),&
               XyzGrid_DI(:,iTrapezoidal2_I(3)),&
               XyzGrid_DI(:,iTrapezoidal2_I(4)),&
               XyzGrid_DI(:,iTrapezoidal2_I(5)),Xyz_D, Weight_I(1:5) )
          if(all(Weight_I(1:5).ge.0.0))then
             nGridOut = 5
             iOrder_I(1:5) = iTrapezoidal2_I
             RETURN
          end if             ! 2
       end if                ! 1
    end if                   ! no trapezoid
  end subroutine interpolate_pyramids
  !======================
  subroutine interpolate_on_parallel_rays(&
       XyzGrid_DI, Xyz_D, iOrder_I, Weight_I, nGridOut, Dir_D, &
       iDRectangle1_I, iDTrapezoid1_I, iDTrapezoid2_I,   &
       iDTriangle1_I, iDTriangle2_I, iDTriangle3_I,&
       iDTriangle4_I, iDTriangle5_I, iDTriangle6_I,&
       iURectangle1_I, iURectangle2_I, iURectangle3_I, iUTrapezoid1_I,&
       iUTriangle1_I, iUTriangle2_I, iUTriangle3_I, iUTriangle4_I)
    real, intent(in):: XyzGrid_DI(nDim,nGrid), Xyz_D(nDim)
    real, intent(out):: Weight_I(nGrid)
    integer, intent(out)::nGridOut, iOrder_I(8)
    !Direction of parallel rays
    real, intent(in) :: Dir_D(nDim)
    !\
    !Up subfaces ure in the direction of  Dir_D from Xyz point
    !Down faces are in the direction of -Dir_D
    !/
    integer, intent(in), optional, dimension(3)::&
         iDTriangle1_I, iDTriangle2_I, iDTriangle3_I,&
         iDTriangle4_I, iDTriangle5_I, iDTriangle6_I,&
         iUTriangle1_I, iUTriangle2_I, iUTriangle3_I, iUTriangle4_I
    integer, intent(in), optional, dimension(4)::&
         iDRectangle1_I, iDTrapezoid1_I, iDTrapezoid2_I, iUTrapezoid1_I,&
         iURectangle1_I, iURectangle2_I, iURectangle3_I
    !\
    !The point belonging to one of the up and down subfaces
    !/ 
    real, dimension(nDim) :: XyzUp_D, XyzDown_D, X1_D, X2_D, X3_D, X4_D
    !Misc
    real:: AlphaUp, AlphaDown
    integer:: nGridOutUp, nGridOutDown
    !\
    ! Consider a ray passing through Xyz
    ! Solve equation:
    ! Xyz + Dir*AlphaUp = XyzUp 
    ! its solution is 
    ! AlphaUp = (UX1-Xyz)\cdot[(UX2-UX1)\times(UX3-UX1)]/&
    !          Dir\cdot[(UX2-UX1)\times(UX3-UX1)]
    ! XyzUp belongs to the up subface
    !/
    !\
    ! Solve equation:
    ! Xyz - Dir*AlphaDown = XyzDown 
    ! its solution is 
    ! AlphaDown = (Xyz-DX1)\cdot[(DX2-DX1)\times(DX3-DX1)]/&
    ! Dir\cdot[(DX2-DX1)\times(DX3-DX1)]
    ! XyzUp belongs to the Down subface
    !/
    !\
    ! As long as Xyz=(XyzDown*AlphaUp + XyzUp*AlphaDown)/(AlphaUp+AlphaDown)
    ! the weights for the upper face should be multiplied by 
    ! AlphaDown/(AlphaUp+AlphaDown)
    ! for the down face - by AlphaUp/(AlphaUp+AlphaDown)
    !/
    nGridOut = -1; nGridOutUp= -1; nGridOutDown = -1
    Weight_I = -1; iOrder_I = 0
    if(present(iUTriangle1_I))then
       X1_D = XyzGrid_DI(:,iUTriangle1_I(1))
       X2_D = XyzGrid_DI(:,iUTriangle1_I(2))
       X3_D = XyzGrid_DI(:,iUTriangle1_I(3))
       AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
            triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
       XyzUp_D = Xyz_D + AlphaUp * Dir_D
       call triangle(X1_D, X2_D, X3_D, XyzUp_D, Weight_I(1:3))
       if(all(Weight_I(1:3)>=0.0))then
          !\
          !The ray for point Xyz is projected into this triangle
          !/
          if(AlphaUp==0.0)then
             !\
             !Point Xyz belongs to this traingle
             !/
             iOrder_I(1:3) = iUTriangle1_I
             nGridOut = 3
             RETURN
          elseif(AlphaUp < 0)then
             !\
             ! Point is above the upper triangle
             ! return with negative weights and nGridOut = -1
             RETURN 
          else
             iOrder_I(1:3) = iUTriangle1_I
             nGridOutUp = 3
             !\
             ! Exit if for triangles with positive AlphaUp and nGridOut = 3
             !/
          end if
       elseif(present(iUTriangle2_I))then
          X1_D = XyzGrid_DI(:,iUTriangle2_I(1))
          X2_D = XyzGrid_DI(:,iUTriangle2_I(2))
          X3_D = XyzGrid_DI(:,iUTriangle2_I(3))
          AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
               triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
          XyzUp_D = Xyz_D + AlphaUp * Dir_D
          call triangle(X1_D, X2_D, X3_D, XyzUp_D, Weight_I(1:3))
          if(all(Weight_I(1:3)>=0.0))then
             !\
             !The ray for point Xyz is projected into this triangle
             !/
             if(AlphaUp==0.0)then
                !\
                !Point Xyz belongs to this traingle
                !/
                iOrder_I(1:3) = iUTriangle2_I
                nGridOut = 3
                RETURN
             elseif(AlphaUp < 0)then
                !\
                ! Point is above the upper triangle
                ! return with negative weights and nGridOut = -1
                RETURN 
             else             
                iOrder_I(1:3) = iUTriangle2_I
                nGridOutUp = 3
                !\
                ! Exit if for triangles with positive AlphaUp
                !/
             end if
          elseif(present(iUTriangle3_I))then
             X1_D = XyzGrid_DI(:,iUTriangle3_I(1))
             X2_D = XyzGrid_DI(:,iUTriangle3_I(2))
             X3_D = XyzGrid_DI(:,iUTriangle3_I(3))
             AlphaUp = &
                  triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                  triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
             XyzUp_D = Xyz_D + AlphaUp * Dir_D
             call triangle(X1_D, X2_D, X3_D, XyzUp_D, Weight_I(1:3))
             if(all(Weight_I(1:3)>=0.0))then
                !\
                !The ray for point Xyz is projected into this triangle
                !/
                if(AlphaUp==0.0)then
                   !\
                   !Point Xyz belongs to this traingle
                   !/
                   iOrder_I(1:3) = iUTriangle3_I
                   nGridOut = 3
                   RETURN
                elseif(AlphaUp < 0)then
                   !\
                   ! Point is above the upper triangle
                   ! return with negative weights and nGridOut = -1
                   RETURN 
                else        
                   iOrder_I(1:3) = iUTriangle3_I
                   nGridOutUp = 3
                   !\
                   ! Exit if for triangles with positive AlphaUp 
                   !/

                end if
             elseif(present(iUTriangle4_I))then
                X1_D = XyzGrid_DI(:,iUTriangle4_I(1))
                X2_D = XyzGrid_DI(:,iUTriangle4_I(2))
                X3_D = XyzGrid_DI(:,iUTriangle4_I(3))
                AlphaUp = &
                     triple_product(X1_D - Xyz_D,X2_D - X1_D,X3_D - X1_D)/&
                     triple_product(Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                XyzUp_D = Xyz_D + AlphaUp * Dir_D
                call triangle(X1_D, X2_D, X3_D, XyzUp_D, Weight_I(1:3))
                if(all(Weight_I(1:3)>=0.0))then
                   !\
                   !The ray for point Xyz is projected into this triangle
                   !/
                   if(AlphaUp==0.0)then
                      !\
                      !Point Xyz belongs to this traingle
                      !/
                      iOrder_I(1:3) = iUTriangle4_I
                      nGridOut = 3
                      RETURN
                   elseif(AlphaUp < 0.0)then
                      !\
                      ! Point is above the upper triangle
                      ! return with negative weights and nGridOut = -1
                      RETURN 
                   else   
                      iOrder_I(1:3) = iUTriangle4_I
                      nGridOutUp = 3
                      !\
                      ! Exit if for triangles with positive AlphaUp 
                      !/
                   end if
                end if
             end if       !4
          end if          !3
       end if             !2
    end if                !1
    if(nGridOutUp < 1)then
       if(present(iURectangle1_I))then
          X1_D = XyzGrid_DI(:,iURectangle1_I(1))
          X2_D = XyzGrid_DI(:,iURectangle1_I(2))
          X3_D = XyzGrid_DI(:,iURectangle1_I(3))
          X4_D = XyzGrid_DI(:,iURectangle1_I(4))
          AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
               triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
          XyzUp_D = Xyz_D + AlphaUp * Dir_D
          call rectangle(X1_D, X2_D, X3_D, X4_D, XyzUp_D, Weight_I(1:4))
          if(all(Weight_I(1:4)>=0.0))then
             !\
             !The ray for point Xyz is projected into this rectangle
             !/
             if(AlphaUp==0.0)then
                !\
                !Point Xyz belongs to this rectangle
                !/
                iOrder_I(1:4) = iURectangle1_I
                nGridOut = 4
                RETURN
             elseif(AlphaUp < 0.0)then
                !\
                ! Point is above the upper rectangle
                ! return with negative weights and nGridOut = -1
                RETURN 
             else 
                iOrder_I(1:4) = iURectangle1_I
                nGridOutUp = 4
                !\
                ! Exit if for rectangles with positive AlphaUp 
                !/
             end if
          elseif(present(iURectangle2_I))then
             X1_D = XyzGrid_DI(:,iURectangle2_I(1))
             X2_D = XyzGrid_DI(:,iURectangle2_I(2))
             X3_D = XyzGrid_DI(:,iURectangle2_I(3))
             X4_D = XyzGrid_DI(:,iURectangle2_I(4))
             AlphaUp = triple_product(X1_D - Xyz_D,X2_D - X1_D,X3_D - X1_D)/&
                  triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
             XyzUp_D = Xyz_D + AlphaUp * Dir_D
             call rectangle(X1_D, X2_D, X3_D, X4_D, XyzUp_D, Weight_I(1:4))
             if(all(Weight_I(1:4)>=0.0))then
                !\
                !The ray for point Xyz is projected into this rectangle
                !/
                if(AlphaUp==0.0)then
                   !\
                   !Point Xyz belongs to this rectangle
                   !/
                   iOrder_I(1:4) = iURectangle2_I
                   nGridOut = 4
                   RETURN
                elseif(AlphaUp < 0.0)then
                   !\
                   ! Point is above the upper rectangle
                   ! return with negative weights and nGridOut = -1
                   RETURN 
                else 
                   !\
                   ! Exit if for rectangles with positive AlphaUp 
                   !/
                   iOrder_I(1:4) = iURectangle2_I
                   nGridOutUp = 4
                end if
             elseif(present(iURectangle3_I))then
                X1_D = XyzGrid_DI(:,iURectangle3_I(1))
                X2_D = XyzGrid_DI(:,iURectangle3_I(2))
                X3_D = XyzGrid_DI(:,iURectangle3_I(3))
                X4_D = XyzGrid_DI(:,iURectangle3_I(4))
                AlphaUp = &
                     triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                     triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
                XyzUp_D = Xyz_D + AlphaUp * Dir_D
                call rectangle(X1_D, X2_D, X3_D, X4_D, XyzUp_D, Weight_I(1:4))
                if(all(Weight_I(1:4)>=0.0))then
                   !\
                   !The ray for point Xyz is projected into this rectangle
                   !/
                   if(AlphaUp==0.0)then
                      !\
                      !Point Xyz belongs to this rectangle
                      !/
                      iOrder_I(1:4) = iURectangle3_I
                      nGridOut = 4
                      RETURN
                   elseif(AlphaUp < 0.0)then
                      !\
                      ! Point is above the upper rectangle
                      ! return with negative weights and nGridOut = -1
                      RETURN 
                   else 
                      !\
                      ! Exit if for rectangles with positive AlphaUp 
                      !/
                      iOrder_I(1:4) = iURectangle3_I
                      nGridOutUp = 4
                   end if
                end if
             end if     !3
          end if        !2
       end if           !1
    end if
    if(nGridOutUp < 1)then
       if(present(iUTrapezoid1_I))then
          X1_D = XyzGrid_DI(:,iUTrapezoid1_I(1))
          X2_D = XyzGrid_DI(:,iUTrapezoid1_I(2))
          X3_D = XyzGrid_DI(:,iUTrapezoid1_I(3))
          X4_D = XyzGrid_DI(:,iUTrapezoid1_I(4))
          AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
               triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
          XyzUp_D = Xyz_D + AlphaUp * Dir_D
          call Trapezoid(X1_D, X2_D, X3_D, X4_D, XyzUp_D, Weight_I(1:4))
          if(all(Weight_I(1:4)>=0.0))then
             !\
             !The ray for point Xyz is projected into this trapezoid
             !/
             if(AlphaUp==0.0)then
                !\
                !Point Xyz belongs to this rectangle
                !/
                iOrder_I(1:4) = iUTrapezoid1_I
                nGridOut = 4
                RETURN
             elseif(AlphaUp < 0.0)then
                !\
                ! Point is above the upper trapezoid
                ! return with negative weights and nGridOut = -1
                RETURN 
             else 
                !\
                ! Exit if for rectangles with positive AlphaUp 
                !/
                iOrder_I(1:4) = iUTrapezoid1_I
                nGridOutUp = 4
             end if
          end if
       end if
    end if
    !\
    !If no intersection point with upper boundary is found
    !/
    if(nGridOutUp == -1)RETURN
    !\
    !Calculate low face
    if(present(iDTriangle1_I))then
       Weight_I(4:3+nGridOutUp) = Weight_I(1:nGridOutUp)
       Weight_I(1:3) = 0
       iOrder_I(4:3+nGridOutUp) = iOrder_I(1:nGridOutUp)
       iOrder_I(1:3) = 0
       X1_D = XyzGrid_DI(:,iDTriangle1_I(1))
       X2_D = XyzGrid_DI(:,iDTriangle1_I(2))
       X3_D = XyzGrid_DI(:,iDTriangle1_I(3))
       AlphaDown = &
            triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
            triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
       XyzDown_D = Xyz_D - AlphaDown*Dir_D
       call triangle(X1_D, X2_D, X3_D, XyzDown_D, Weight_I(1:3))
       if(all(Weight_I(1:3)>=0.0))then
          !\
          !The ray for point Xyz is projected into this triangle
          !/
          if(AlphaDown==0.0)then
             !\
             !Point Xyz belongs to this triangle
             !/
             iOrder_I(1:3) = iDTriangle1_I
             nGridOut = 3
             RETURN
          elseif(AlphaDown < 0.0)then
             !\
             ! Point is below the lower triangle
             ! return with negative weights and nGridOut = -1
             !/
             RETURN 
          else
             !\
             ! Exit if for triangles with positive AlphaDown 
             !/
             iOrder_I(1:3) = iDTriangle1_I
             nGridOutDown = 3
          end if
       elseif(present(iDTriangle2_I))then
          X1_D = XyzGrid_DI(:,iDTriangle2_I(1))
          X2_D = XyzGrid_DI(:,iDTriangle2_I(2))
          X3_D = XyzGrid_DI(:,iDTriangle2_I(3))
          AlphaDown = &
               triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
               triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
          XyzDown_D = Xyz_D - AlphaDown * Dir_D
          call triangle(X1_D, X2_D, X3_D, XyzDown_D, Weight_I(1:3))
          if(all(Weight_I(1:3)>=0.0))then
             !\
             !The ray for point Xyz is projected into this triangle
             !/
             if(AlphaDown==0.0)then
                !\
                !Point Xyz belongs to this triangle
                !/
                iOrder_I(1:3) = iDTriangle2_I
                nGridOut = 3
                RETURN
             elseif(AlphaDown < 0.0)then
                !\
                ! Point is below the lower triangle
                ! return with negative weights and nGridOut = -1
                !/
                RETURN
             else
                !\
                ! Exit if for triangles with positive AlphaDown 
                !/
                iOrder_I(1:3) = iDTriangle2_I
                nGridOutDown = 3
             end if
          elseif(present(iDTriangle3_I))then
             X1_D = XyzGrid_DI(:,iDTriangle3_I(1))
             X2_D = XyzGrid_DI(:,iDTriangle3_I(2))
             X3_D = XyzGrid_DI(:,iDTriangle3_I(3))
             AlphaDown = &
                  triple_product(Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                  triple_product(Dir_D       ,X2_D - X1_D,X3_D - X1_D)
             XyzDown_D = Xyz_D - AlphaDown * Dir_D
             call triangle(X1_D, X2_D, X3_D, XyzDown_D, Weight_I(1:3))
             if(all(Weight_I(1:3)>=0.0))then
                !\
                !The ray for point Xyz is projected into this triangle
                !/
                if(AlphaDown==0.0)then
                   !\
                   !Point Xyz belongs to this triangle
                   !/
                   iOrder_I(1:3) = iDTriangle3_I
                   nGridOut = 3
                   RETURN
                elseif(AlphaDown < 0.0)then
                   !\
                   ! Point is below the lower triangle
                   ! return with negative weights and nGridOut = -1
                   !/
                   RETURN
                else
                   !\
                   ! Exit if for triangles with positive AlphaDown 
                   !/
                   iOrder_I(1:3) = iDTriangle3_I
                   nGridOutDown = 3
                end if
             elseif(present(iDTriangle4_I))then
                X1_D = XyzGrid_DI(:,iDTriangle4_I(1))
                X2_D = XyzGrid_DI(:,iDTriangle4_I(2))
                X3_D = XyzGrid_DI(:,iDTriangle4_I(3))
                AlphaDown = &
                     triple_product(Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                     triple_product(Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                XyzDown_D = Xyz_D - AlphaDown * Dir_D
                call triangle(X1_D, X2_D, X3_D, XyzDown_D, Weight_I(1:3))
                if(all(Weight_I(1:3)>=0.0))then
                   !\
                   !The ray for point Xyz is projected into this triangle
                   !/
                   if(AlphaDown==0.0)then
                      !\
                      !Point Xyz belongs to this triangle
                      !/
                      iOrder_I(1:3) = iDTriangle4_I
                      nGridOut = 3
                      RETURN
                   elseif(AlphaDown < 0.0)then
                      !\
                      ! Point is below the lower triangle
                      ! return with negative weights and nGridOut = -1
                      !/
                      RETURN
                   else
                      !\
                      ! Exit if for triangles with positive AlphaDown 
                      !/
                      iOrder_I(1:3) = iDTriangle4_I
                      nGridOutDown = 3
                   end if
                elseif(present(iDTriangle5_I))then
                   X1_D = XyzGrid_DI(:,iDTriangle5_I(1))
                   X2_D = XyzGrid_DI(:,iDTriangle5_I(2))
                   X3_D = XyzGrid_DI(:,iDTriangle5_I(3))
                   AlphaDown = &
                        triple_product(&
                        Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                        triple_product(&
                        Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                   XyzDown_D = Xyz_D - AlphaDown * Dir_D
                   call triangle(X1_D, X2_D, X3_D, XyzDown_D, Weight_I(1:3))
                   if(all(Weight_I(1:3)>=0.0))then
                      !\
                      !The ray for point Xyz is projected into this triangle
                      !/
                      if(AlphaDown==0.0)then
                         !\
                         !Point Xyz belongs to this triangle
                         !/
                         iOrder_I(1:3) = iDTriangle5_I
                         nGridOut = 3
                         RETURN
                      elseif(AlphaDown < 0.0)then
                         !\
                         ! Point is below the lower triangle
                         ! return with negative weights and nGridOut = -1
                         !/
                         RETURN
                      else
                         !\
                         ! Exit if for triangles with positive AlphaDown 
                         !/
                         iOrder_I(1:3) = iDTriangle5_I
                         nGridOutDown = 3
                      end if
                   elseif(present(iDTriangle6_I))then
                      X1_D = XyzGrid_DI(:,iDTriangle6_I(1))
                      X2_D = XyzGrid_DI(:,iDTriangle6_I(2))
                      X3_D = XyzGrid_DI(:,iDTriangle6_I(3))
                      AlphaDown = &
                           triple_product(&
                           Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                           triple_product(&
                           Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                      XyzDown_D = Xyz_D - AlphaDown * Dir_D
                      call triangle(X1_D, X2_D, X3_D, XyzDown_D, Weight_I(1:3))
                      if(all(Weight_I(1:3)>=0.0))then
                         !\
                         !The ray for point Xyz is projected into this 
                         !triangle
                         !/
                         if(AlphaDown==0.0)then
                            !\
                            !Point Xyz belongs to this triangle
                            !/
                            iOrder_I(1:3) = iDTriangle6_I
                            nGridOut = 3
                            RETURN
                         elseif(AlphaDown < 0.0)then
                            !\
                            ! Point is below the lower triangle
                            ! return with negative weights and nGridOut = -1
                            !/
                            RETURN
                         else
                            !\
                            ! Exit if for triangles with positive AlphaDown 
                            !/
                            iOrder_I(1:3) = iDTriangle6_I
                            nGridOutDown = 3
                         end if
                      end if
                   end if !6
                end if    !5
             end if       !4
          end if          !3
       end if             !2
    end if                !1
    if(nGridOutDown <1)then
       if(present(iDRectangle1_I))then
          if(present(iDTriangle1_I))then
             Weight_I(5:4+nGridOutUp) = Weight_I(4:3+nGridOutUp)
             Weight_I(1:4) = 0
             iOrder_I(5:4+nGridOutUp) = iOrder_I(4:3+nGridOutUp)
             iOrder_I(1:4)=0
          else
             Weight_I(5:4+nGridOutUp) = Weight_I(1:nGridOutUp)
             Weight_I(1:4) = 0
             iOrder_I(5:4+nGridOutUp) = iOrder_I(1:nGridOutUp)
             iOrder_I(1:4)=0
          end if
          X1_D = XyzGrid_DI(:,iDRectangle1_I(1))
          X2_D = XyzGrid_DI(:,iDRectangle1_I(2))
          X3_D = XyzGrid_DI(:,iDRectangle1_I(3))
          X4_D = XyzGrid_DI(:,iDRectangle1_I(4))
          AlphaDown = triple_product(Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
               triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
          XyzDown_D = Xyz_D - AlphaDown * Dir_D
          call rectangle(X1_D, X2_D, X3_D, X4_D, XyzDown_D, Weight_I(1:4))
          if(all(Weight_I(1:4)>=0.0))then
             !\
             !The ray for point Xyz is projected into this rectangle
             !/
             if(AlphaDown==0.0)then
                !\
                !Point Xyz belongs to this rectangle
                !/
                iOrder_I(1:4) = iDRectangle1_I
                nGridOut = 4
                RETURN
             elseif(AlphaDown < 0.0)then
                !\
                ! Point is below the lower rectangle
                ! return with negative weights and nGridOut = -1
                !/
                RETURN
             else
                !\
                ! Exit if for rectangles with positive AlphaDown 
                !/
                iOrder_I(1:4) = iDRectangle1_I
                nGridOutDown = 4
             end if
          end if
       end if
    end if           !1
    if(nGridOutDown < 1)then
       if(present(iDTrapezoid1_I))then
          if(.not.present(iDRectangle1_I))then
             if(present(iDTriangle1_I))then
                Weight_I(5:4+nGridOutUp) = Weight_I(4:3+nGridOutUp)
                Weight_I(1:4)=0
                iOrder_I(5:4+nGridOutUp) = iOrder_I(4:3+nGridOutUp)
                iOrder_I(1:4) = 0
             else
                Weight_I(5:4+nGridOutUp) = Weight_I(1:nGridOutUp)
                Weight_I(1:4)=0
                iOrder_I(5:4+nGridOutUp) = iOrder_I(1:nGridOutUp)
                iOrder_I(1:4) = 0
             end if
          end if
          X1_D = XyzGrid_DI(:,iDTrapezoid1_I(1))
          X2_D = XyzGrid_DI(:,iDTrapezoid1_I(2))
          X3_D = XyzGrid_DI(:,iDTrapezoid1_I(3))
          X4_D = XyzGrid_DI(:,iDTrapezoid1_I(4))
          AlphaDown = triple_product(Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
               triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
          XyzDown_D = Xyz_D - AlphaDown * Dir_D
          call trapezoid(X1_D, X2_D, X3_D, X4_D, XyzDown_D, Weight_I(1:4))
          if(all(Weight_I(1:4)>=0.0))then
             !\
             !The ray for point Xyz is projected into this trapezoid
             !/
             if(AlphaDown==0.0)then
                !\
                !Point Xyz belongs to this rectangle
                !/
                iOrder_I(1:4) = iDTrapezoid1_I
                nGridOut = 4
                RETURN
             elseif(AlphaDown < 0.0)then
                !\
                ! Point is below the lower rectangle
                ! return with negative weights and nGridOut = -1
                !/
                RETURN
             else
                !\
                ! Exit if for rectangles with positive AlphaDown 
                !/
                iOrder_I(1:4) = iDTrapezoid1_I
                nGridOutDown = 4
             end if
          elseif(present(iDTrapezoid2_I))then
             X1_D = XyzGrid_DI(:,iDTrapezoid2_I(1))
             X2_D = XyzGrid_DI(:,iDTrapezoid2_I(2))
             X3_D = XyzGrid_DI(:,iDTrapezoid2_I(3))
             X4_D = XyzGrid_DI(:,iDTrapezoid2_I(4))
             AlphaDown = triple_product(Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                  triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
             XyzDown_D = Xyz_D - AlphaDown * Dir_D
             call trapezoid(X1_D, X2_D, X3_D, X4_D, XyzDown_D, Weight_I(1:4))
             if(all(Weight_I(1:4)>=0.0))then
                !\
                !The ray for point Xyz is projected into this trapezoid
                !/
                if(AlphaDown==0.0)then
                   !\
                   !Point Xyz belongs to this rectangle
                   !/
                   iOrder_I(1:4) = iDTrapezoid2_I
                   nGridOut = 4
                   RETURN
                elseif(AlphaDown < 0.0)then
                   !\
                   ! Point is below the lower rectangle
                   ! return with negative weights and nGridOut = -1
                   !/
                   RETURN
                else
                   !\
                   ! Exit if for rectangles with positive AlphaDown 
                   !/
                   iOrder_I(1:4) = iDTrapezoid2_I
                   nGridOutDown = 4
                end if
             end if
          end if
       end if           !1
    end if
    !\
    !If no intersection point with the down subface is found
    !/
    if(nGridOutDown == -1)RETURN 
    !\
    !Apply the weight for interpolation along the rate
    !/
    nGridOut = nGridOutUp + nGridOutDown
    Weight_I(1:nGridOutDown) = &
         Weight_I(1:nGridOutDown)*AlphaUp/(AlphaUp + AlphaDown)
    Weight_I(nGridOutDown+1:nGridOut) = &
         Weight_I(nGridOutDown+1:nGridOut)*AlphaDown/(AlphaUp + AlphaDown)
  end subroutine interpolate_on_parallel_rays
end module ModInterpolateSimpleShape
!============================
module ModCubeGeometry
  implicit none
  SAVE
  !=================ARRAYS FOR A CUBE====================================!
  !For a cubic stencil enumerated as follows:
  !
  !       ^
  ! z-axis|
  !       |  7----------8
  !       | /|         /|
  !       |/ |        / |
  !       5----------6  |
  !       |  |   _   |  |
  !       |  |   /|y-axi$
  !       |  |  /    |  |
  !       |  | /     |  |
  !       |  |/      |  |
  !       |  3----------4
  !       | /        | /
  !       |/         |/
  !       1----------2------> x-axis
  !       
  ! we provide several functions characterizing its geometry in terms of 
  ! the grid point numbers. The functions used in more than one place,
  ! therefore we delegate them here.   

  !Vertexes (enumerated by the first undex), which form 
  !the face of direction iDir (second index) including the 
  !given vertex (the third index). 
  !When the first index equals 1,2,3,4, the vertex, accordingly: 
  !(1) coincides with the given one
  !(2) is connected with the given one by the edge of direction iDir+1
  !(3) is connected with the given one by the edge of direction iDir+2
  !(4) is connected to the given one by the face diagonal of direction iDir
  integer, dimension(4,3,8), parameter:: iFace_IDI = &
       reshape((/&   ! yz face   ! xz face      ! xy face ! 
       1, 3, 5, 7,   1, 5, 2, 6,   1, 2, 3, 4, & 
       2, 4, 6, 8,   2, 6, 1, 5,   2, 1, 4, 3, &
       3, 1, 7, 5,   3, 7, 4, 8,   3, 4, 1, 2, &
       4, 2, 8, 6,   4, 8, 3, 7,   4, 3, 2, 1, &
       5, 7, 1, 3,   5, 1, 6, 2,   5, 6, 7, 8, &
       6, 8, 2, 4,   6, 2, 5, 1,   6, 5, 8, 7, &
       7, 5, 3, 1,   7, 3, 8, 4,   7, 8, 5, 6, &
       8, 6, 4, 2,   8, 4, 7, 3,   8, 7, 6, 5  /), (/4,3,8/))

  !Vertexes (enumerated by the first undex), which form 
  !the face of direction iDir (second index) and does not 
  !include the given vertex (the third index). 
  !When the first index equals 1,2,3,4, the vertex, accordingly: 
  !(1) is connected to the given one by the edge of direction iDir
  !(2) is connected to (1) by the edge of direction iDir+1
  !(3) is connected to (1) by the edge of direction iDir+2
  !(4) is connected to  by the face diagonal of direction iDir
  integer, dimension(4,3,8), parameter:: iOppositeFace_IDI = &
       reshape((/&   ! yz face   ! xz face      ! xy face ! 
       2, 4, 6, 8,   3, 7, 4, 8,   5, 6, 7, 8, & 
       1, 3, 5, 7,   4, 8, 3, 7,   6, 5, 8, 7, &
       4, 2, 8, 6,   1, 5, 2, 6,   7, 8, 5, 6, &
       3, 1, 7, 5,   2, 6, 1, 5,   8, 7, 6, 5, &
       6, 8, 2, 4,   7, 3, 8, 4,   1, 2, 3, 4, &
       5, 7, 1, 3,   8, 4, 7, 3,   2, 1, 4, 3, &
       8, 6, 4, 2,   5, 1, 6, 2,   3, 4, 1, 2, &
       7, 5, 3, 1,   6, 2, 5, 1,   4, 3, 2, 1  /), (/4,3,8/))
  !Number of the vertex connected by 
  !the edge of direction iDir (second index) 
  !with the given vertex (first index)
  integer, dimension(8,3), parameter:: iEdge_ID = &
       reshape((/&   !Number of the connected vertex
       2 ,1, 4, 3, 6, 5, 8, 7,        & !Edge x
       3, 4, 1, 2, 7, 8, 5, 6,        & !Edge y
       5, 6, 7, 8, 1, 2, 3, 4         & !Edge z
       /),(/8,3/))
  !=================ARRAYS FOR A RECTANGLE==============!
  !For a cubic stencil enumerated as follows:
  !
  !       ^
  ! y-axis|
  !       |   
  !       3----------4  
  !       |          |  
  !       |          !
  !       |          |  
  !       |          |  
  !       |          |  
  !       |          !
  !       |          |
  !       |          |
  !       1----------2------> x-axis
  !       
  ! we provide several functions characterizing its geometry in terms of 
  ! the grid point numbers. The functions used in more than one place,
  ! therefore we delegate them here. 


  !Vertexes (enumerated by the first undex), which form 
  !the side of direction iDir (second index) including the 
  !given vertex (the third index). 
  !When the first index equals 1,2,3,4, the vertex, accordingly: 
  !(1) coincides with the given one
  !(2) is connected with the given one by the edge of direction iDir
  integer, dimension(2,2,4), parameter:: iSide_IDI = &
       reshape((/&   ! x side   ! y side ! 
       1, 2,    1, 3, & !Vertex 1
       2, 1,    2, 4, & !Vertex 2
       3, 4,    3, 1, & !Vertex 3
       4, 3,    4, 2  & !Vertex 4
       /), (/2,2,4/))

  !Vertexes (enumerated by the first undex), which form 
  !the side of direction iDir (second index) and does not include the
  !given vertex (the third index). 
  !When the first index equals 1,2 the vertex, accordingly: 
  !(1) is connected to the given one by the side of direction 1+mod(iDir,2)
  !(2) is connected to (1) by the edge of direction iDir
  integer, dimension(2,2,4), parameter:: iOppositeSide_IDI = &
       reshape((/&   ! x side   ! y side ! 
       3, 4,    2, 4, & !Vertex 1
       4, 3,    1, 3, & !Vertex 2
       1, 2,    4, 2, & !Vertex 3
       2, 1,    3, 1  & !Vertex 4
       /), (/2,2,4/))
  !\
  !Array used to sort corners: 
  !(1) assign a type to it (Case_) 
  !(2) assign a grid point of the stencil chosen as the first one (Grid_)
  !(3) assign a direction, if the configuration can be characterized (Dir_)
  !by direction

  logical:: DoInit = .true.
  integer, parameter:: Grid_=1, Dir_=2, Case_=3
  integer:: iSortStencil3_II(Grid_:Case_,0:258) = 0
  integer:: iSortStencil2_II(Grid_:Case_, 0:15) = 0

  !\
  ! Different cases of stencil (will be stored in Case_ column
  !/
  integer, parameter:: Uniform_ = 0, Face_ = 1, Edge_ = 2
  integer, parameter:: FiveTetrahedra_ = 3
  integer, parameter:: OneFine_        = 4, OneCoarse_      = 5
  integer, parameter:: FineMainDiag_   = 6, FineFaceDiag_   = 7
  integer, parameter:: CoarseMainDiag_ = 8, CoarseFaceDiag_ = 9
  integer, parameter:: FineEdgePlusOne_ = 10, ThreeFineOnFace_ =11
  integer, parameter:: CoarseEdgePlusOne_ = 12, ThreeCoarseOnFace_=13
  integer, parameter:: ThreeCoarseOnFacePlusOne_ = 14, CoarseChain_   = 15
  integer, parameter:: Transition2Edge_ = 256, Transition2Corner_ = 257
  integer, parameter:: TransitionJunction_ = 258
  !\
  ! Analogous for 2D
  !/
  integer, parameter:: Trapezoid_ = 1
  integer, parameter:: Rhombus_ = 3
  integer, parameter:: &
       Coarse_  = 0,               &
       Fine_    = 1

  !Parameters to enumerate faces or face diagonal
  integer, parameter:: Xy_ = 3, Xz_ = 2, Yz_ = 1

  !Number of the vertex connected by 
  !the face diagonal across the face of direction iDir (second index) 
  !with the given vertex (first index)
  integer, dimension(8,3), parameter:: iFaceDiag_ID = &
       reshape((/&   !Number of the connected vertex
       7 ,8, 5, 6, 3, 4, 1, 2,        & !Face yz
       6, 5, 8, 7, 2, 1, 4, 3,        & !Face xz
       4, 3, 2, 1, 8, 7, 6, 5         & !Face xy
       /),(/8,3/))

  !Number of the vertex connected by 
  !the main diagonal with the given vertex (index)
  integer, dimension(8), parameter:: iMainDiag_I = &
       (/ 8, 7, 6, 5, 4, 3, 2, 1 /)
contains
  !==========================
  subroutine init_sort_stencil
    integer:: iCase, iGrid, iLoc, iDir, iMisc
    integer:: iLevel_I(8)
    integer :: nDim = 3, nGrid = 8
    integer, parameter:: x_ = 1, y_ = 2, z_ = 3
    !------------------------
    nDim = 3; nGrid = 8;  iSortStencil3_II = 0
    DoInit = .false. !Do this only once
    !Zero and 255 cases are both uniform. !2 cases, 2 totally, 254 left
    iSortStencil3_II(Grid_:Dir_,0:255:255) = 1
    iSortStencil3_II(Case_,0:255:255) = Uniform_
    CASE:do iCase = 1, 254
       !\
       ! Generate 'binary' iLevel_I from iCase
       !/
       iLevel_I=0
       iGrid = 1
       iMisc = iCase
       do while(iMisc > 0)
          iLevel_I(iGrid) = mod(iMisc,2)
          iMisc = (iMisc -  iLevel_I(iGrid))/2
          iGrid = iGrid +1
       end do
       !\
       ! Sort out faces 6 cases, 8 totally, 248 left
       !/
       do iDir = 1, 3
          if(all(iLevel_I(iFace_IDI(:,iDir,1))==Fine_).and.&
               all(iLevel_I(iOppositeFace_IDI(:,iDir,1))==Coarse_))then
             iSortStencil3_II(Grid_,iCase) = 1
             iSortStencil3_II(Dir_, iCase) = iDir !  2*iDir -1
             iSortStencil3_II(Case_,iCase) = Face_
             CYCLE CASE
          elseif(all(iLevel_I(iFace_IDI(:,iDir,1))==Coarse_).and.&
               all(iLevel_I(iOppositeFace_IDI(:,iDir,1))==Fine_))then
             iSortStencil3_II(Grid_,iCase) = 1
             iSortStencil3_II(Dir_, iCase) = iDir !2*iDir 
             iSortStencil3_II(Case_,iCase) = Face_
             CYCLE CASE
          end if
       end do
       !\
       ! Sort out edges 30 cases, 38 totally, 218 left
       !/
       do iDir = 1, 3
          if(all(iLevel_I(iFace_IDI(:,iDir,1))==&
               iLevel_I(iOppositeFace_IDI(:,iDir,1)) ) )then
             iSortStencil3_II(Grid_,iCase) = 1
             iSortStencil3_II(Dir_, iCase) = iDir
             iSortStencil3_II(Case_,iCase) = Edge_
             CYCLE CASE
          end if
       end do
       !\
       ! Find configurations which split into five tetrahedra
       !   26 cases 64 totally 192 left
       !/
       do iGrid = 1,nGrid
          if(  all(iLevel_I(iEdge_ID(iGrid,:))==Coarse_).and.&
               all(iLevel_I(iFaceDiag_ID(iGrid,:))==Fine_))then
             iSortStencil3_II(Dir_,iCase)  = 1
             iSortStencil3_II(Grid_,iCase) = iGrid
             iSortStencil3_II(Case_,iCase) = FiveTetrahedra_
             CYCLE CASE
          end if
       end do
       select case( count(iLevel_I==Fine_))
       case(1)                          ! 8 cases 72 totally 184 left 
          iLoc = maxloc(iLevel_I,DIM=1)
          iSortStencil3_II(Dir_,iCase) = 1
          iSortStencil3_II(Grid_,iCase) = iLoc
          iSortStencil3_II(Case_,iCase) = OneFine_
          CYCLE CASE 
       case(7)                          ! 8 cases 80 totally 176 left
          iLoc = minloc(iLevel_I,DIM=1)
          iSortStencil3_II(Dir_,iCase)  = 1
          iSortStencil3_II(Grid_,iCase) = iLoc
          iSortStencil3_II(Case_,iCase) = OneCoarse_
          CYCLE CASE
       case(2)
          iLoc = maxloc(iLevel_I,DIM=1)
          if(iLevel_I(iMainDiag_I(iLoc))==Fine_)then 
             ! 4 cases   84 totally 172 left
             iSortStencil3_II(Dir_,iCase)  = 1
             iSortStencil3_II(Grid_,iCase) = iLoc
             iSortStencil3_II(Case_,iCase) = FineMainDiag_
             CYCLE CASE
          else                      ! 12 cases  96 totally 160 left 
             do iDir = Yz_,Xy_
                if(iLevel_I(iFaceDiag_ID(iLoc, iDir))==Fine_)then
                   iSortStencil3_II(Dir_,iCase)  = iDir
                   iSortStencil3_II(Grid_,iCase) = iLoc
                   iSortStencil3_II(Case_,iCase) = FineFaceDiag_
                   CYCLE CASE
                end if
             end do
          end if
       case(6)                       ! 4 cases 100 totally 156 left
          iLoc = minloc(iLevel_I,DIM=1)
          if(iLevel_I(iMainDiag_I(iLoc))==Coarse_)then  
             iSortStencil3_II(Dir_,iCase)  = 1
             iSortStencil3_II(Grid_,iCase) = iLoc
             iSortStencil3_II(Case_,iCase) = CoarseMainDiag_
             CYCLE CASE
          else                      ! 12 cases 112 totally 144 left
             do iDir = Yz_,Xy_
                if(iLevel_I(iFaceDiag_ID(iLoc, iDir))==Coarse_)then
                   iSortStencil3_II(Dir_,iCase) = iDir
                   iSortStencil3_II(Grid_,iCase) = iLoc
                   iSortStencil3_II(Case_,iCase) = CoarseFaceDiag_
                   CYCLE CASE
                end if
             end do
          end if
       case(3)          !         24 cases  136 totally 120 left
          do iGrid = 1, nGrid
             if(iLevel_I(iGrid)==Fine_)then
                if(iLevel_I(iMainDiag_I(iGrid))==Fine_)then
                   do iDir = Yz_,Xy_
                      if(iLevel_I(iFaceDiag_ID(iGrid, iDir))==Fine_)then
                         iSortStencil3_II(Dir_,iCase)  = iDir
                         iSortStencil3_II(Grid_,iCase) = iGrid
                         iSortStencil3_II(Case_,iCase) = FineEdgePlusOne_
                         CYCLE CASE
                      end if
                   end do
                else            ! 24 cases 160 totally 96 left
                   do iDir = 1,nDim
                      if(all(iLevel_I(iFace_IDI(1:3,iDir,iGrid))==Fine_))then 
                         iSortStencil3_II(Dir_,iCase)  = iDir
                         iSortStencil3_II(Grid_,iCase) = iGrid
                         iSortStencil3_II(Case_,iCase) = ThreeFineOnFace_
                         CYCLE CASE
                      end if
                   end do
                end if
             end if
          end do
       case(5)          !         24 cases 184 totally 72 left 
          do iGrid = 1, nGrid
             if(iLevel_I(iGrid)/=Coarse_)CYCLE
             if(iLevel_I(iMainDiag_I(iGrid))==Coarse_)then
                do iDir = Yz_,Xy_
                   if(iLevel_I(iFaceDiag_ID(iGrid, iDir))==Coarse_)then 
                      iSortStencil3_II(Dir_,iCase)  = iDir
                      iSortStencil3_II(Grid_,iCase) = iGrid
                      iSortStencil3_II(Case_,iCase) = CoarseEdgePlusOne_
                      CYCLE CASE          
                   end if
                end do
             end if
             !                     24 cases 208 totally 48 left
             do iDir = 1,nDim
                if(all(iLevel_I(iFace_IDI(1:3, iDir, iGrid))==Coarse_))then
                   iSortStencil3_II(Dir_,iCase)  = iDir
                   iSortStencil3_II(Grid_,iCase) = iGrid
                   iSortStencil3_II(Case_,iCase) = ThreeCoarseOnFace_
                   CYCLE CASE
                end if
             end do
          end do
       case(4)                 
          do iGrid = 1,nGrid
             if(iLevel_I(iGrid)==Coarse_)then
                do iDir = 1,nDim    !   24 cases 232 totally 24 left 
                   if(all(&
                        iLevel_I(iOppositeFace_IDI(2:4,iDir,iGrid))==Coarse_&
                        ))then
                      iSortStencil3_II(Dir_,iCase)  = iDir
                      iSortStencil3_II(Grid_,iCase) = iGrid
                      iSortStencil3_II(Case_,iCase) = &
                           ThreeCoarseOnFacePlusOne_
                      CYCLE CASE
                   end if
                   !
                   !\
                   ! 24 cases 256 totally 0 left
                   !/
                   if(iLevel_I(iEdge_ID(iGrid,iDir))==Coarse_.and.&
                        iLevel_I(iEdge_ID(iGrid,1 + mod(iDir,3)))==Coarse_&
                        .and.iLevel_I(iEdge_ID(iEdge_ID(iGrid,iDir),&
                        1 + mod(1 + iDir,3)))==Coarse_)then
                      iSortStencil3_II(Dir_,iCase)  = iDir
                      iSortStencil3_II(Grid_,iCase) = iGrid
                      iSortStencil3_II(Case_,iCase) = CoarseChain_
                      CYCLE CASE
                   end if
                end do
             end if
          end do
       end select
    end do CASE
    iSortStencil3_II(Grid_:Dir_,Transition2Edge_:TransitionJunction_)  = 1
    iSortStencil3_II(Case_,Transition2Edge_:TransitionJunction_)  = &
         (/Transition2Edge_, Transition2Corner_, TransitionJunction_/)
    !=================2 dimensional case=======
    nDim = 2; nGrid = 4
    !Zero and 15 cases are both uniform. !2 cases, 2 totally, 14 left
    iSortStencil2_II(Grid_:Dir_,0:15:15) = 1
    iSortStencil2_II(Case_,0:15:15) = Uniform_
    CASE2:do iCase = 1, 14
       !\
       ! Generate 'binary' iLevel_I from iCase
       !/
       iLevel_I=0
       iGrid = 1
       iMisc = iCase
       do while(iMisc > 0)
          iLevel_I(iGrid) = mod(iMisc,2)
          iMisc = (iMisc -  iLevel_I(iGrid))/2
          iGrid = iGrid +1
       end do
       !\
       ! Sort out faces 4 cases, 6 totally, 10 left
       !/
       do iDir = 1, nDim
          if(all(iLevel_I(iSide_IDI(:,iDir,1))==Fine_).and.&
               all(iLevel_I(iOppositeSide_IDI(:,iDir,1))==Coarse_))then
             iSortStencil2_II(Grid_,iCase) = 1
             iSortStencil2_II(Dir_, iCase) = iDir !  2*iDir -1
             iSortStencil2_II(Case_,iCase) = Trapezoid_
             CYCLE CASE2
          elseif(all(iLevel_I(iSide_IDI(:,iDir,1))==Coarse_).and.&
               all(iLevel_I(iOppositeSide_IDI(:,iDir,1))==Fine_))then
             iSortStencil2_II(Grid_,iCase) = 1
             iSortStencil2_II(Dir_, iCase) = iDir !2*iDir 
             iSortStencil2_II(Case_,iCase) = Trapezoid_
             CYCLE CASE2
          end if
       end do
       select case( count(iLevel_I(1:4)==Fine_))
       case(1)                          ! 4 cases 10 totally 6 left 
          iLoc = maxloc(iLevel_I,DIM=1)
          iSortStencil2_II(Dir_,iCase) = 1
          iSortStencil2_II(Grid_,iCase) = iLoc
          iSortStencil2_II(Case_,iCase) = OneFine_
          CYCLE CASE2
       case(3)                          ! 4 cases 14 totally 2 left
          iLoc = minloc(iLevel_I(1:4),DIM=1)
          iSortStencil2_II(Dir_,iCase)  = 1
          iSortStencil2_II(Grid_,iCase) = iLoc
          iSortStencil2_II(Case_,iCase) = OneCoarse_
          CYCLE CASE2
       case(2)                          !          
          iLoc = maxloc(iLevel_I,DIM=1)
          iSortStencil2_II(Dir_,iCase) = 1
          iSortStencil2_II(Grid_,iCase) = iLoc
          iSortStencil2_II(Case_,iCase) = Rhombus_
       end select
    end do CASE2
  end subroutine init_sort_stencil
  !============================
  integer function i_case(iLevel_I)
    integer, intent(in):: iLevel_I(:)
    integer:: iGrid, nGrid
    !--------------------
    nGrid = size(iLevel_I,DIM=1)
    !\
    ! Generate iCase from 'binary' iLevel_I
    !/
    i_case = iLevel_I(nGrid)
    do iGrid = nGrid - 1, 1, -1
       i_case = i_case + i_case + iLevel_I(iGrid)
    end do
  end function i_case
  !============================
end module ModCubeGeometry
!=========================
module ModResolutionCorner
  use ModCubeGeometry, ONLY: iSortStencil3_II, iFace_IDI, iOppositeFace_IDI,&
       Case_, Grid_, Dir_, FiveTetrahedra_, & 
       OneFine_, OneCoarse_, FineMainDiag_, FineFaceDiag_,            &
       CoarseMainDiag_, CoarseFaceDiag_, FineEdgePlusOne_,            &
       ThreeFineOnFace_, CoarseEdgePlusOne_, ThreeCoarseOnFace_,      &
       ThreeCoarseOnFacePlusOne_, CoarseChain_, i_case
  implicit none
  PRIVATE
  SAVE
  public:: resolution_corner
contains
  !==================================
  subroutine  resolution_corner(&
       Xyz_D, XyzGridIn_DI, iLevel_I, nGridOut, Weight_I, iOrder_I)
    character(LEN=*),parameter:: NameSub='resolution_corner'
    integer, parameter:: nDim = 3, nGrid = 8
    !\
    !Input parameters
    !/
    !\
    !The location at which to interpolate the data
    !/
    real   , intent(in) :: Xyz_D(nDim) 
    !\
    !Grid point coordinates !3 coordinate, 8 points
    !/
    real  ,   intent(in) :: XyzGridIn_DI(nDim,nGrid)
    !\
    ! The same, but may be reordered, if needed 
    real:: XyzGrid_DI(nDim,nGrid)
    !The refinement level at each grid point. By one higher level
    ! of refinement assumes the cell size reduced by a factor of 0.5
    integer,  intent(in) :: iLevel_I(nGrid)
    !\
    !Output parameters
    !/
    !The number of grid points to be ultimately included into 
    !the interpolation stencil If nGridOut < nGrid only the
    !first nGridOut lines in the output arrays are meaningful.
    !/
    integer, intent(out) :: nGridOut
    !\
    !The weight coefficients array.
    real   , intent(out) :: Weight_I(nGrid)
    !\
    !Order(numbers) of grid points used for the interpolation
    integer, intent(out) :: iOrder_I(nGrid)    
    !\
    ! Minimum and maximum values of coordinates
    !/
    real    :: XyzMin_D(nDim), XyzMax_D(nDim)   
    !\
    ! To find the sort of corner stencil
    !/
    integer:: iCase, iDir, iGrid
    !-------------
    Weight_I = 0;  XyzGrid_DI = XyzGridIn_DI
    !\
    ! We store in advance the 'basic' grid point
    ! and orientation of all possible stencil configurations
    ! and now we extracted this information
    !/ 
    iCase = i_case(iLevel_I)
    iGrid = iSortStencil3_II(Grid_,iCase)
    iDir  = iSortStencil3_II(Dir_,iCase)
    iCase = iSortStencil3_II(Case_,iCase)
    iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
    iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
    XyzGrid_DI = XyzGrid_DI(:,iOrder_I)
    select case(iCase)
    case(FiveTetrahedra_)
       !\
       ! One configuration resulting in the corner split for five tetrahedra
       ! is occured in many cases (26), however, it is treated in the same  
       ! way in all these cases. If the given grid point is connected by 
       ! three edges with three coarse points, while the across the main  
       ! diagonal point is connected by three edges (or, equivalently, the  
       ! given point is connected by three face diagonals) with three coarse  
       ! points), then all faces are trangulated   
       !/
       !\
       ! View from  point iOrder
       !            xy face diag
       !   x-axis      |      y-axis
       !         \     V     /
       !          \  / F    /           
       !           C   | \C           
       !           |\  | / |            
       !           | \ |/  |            
       !           |  \/   |
       !           |   1   |
       !  xz facediag/ | \yz face diag
       !           /   |  \|
       !           F   |   F
       !             \ C /
       !               |
       !               V
       !               z-axis

       call interpolate_pyramids(&
            iTetrahedron1_I=(/1, 2, 6, 4/),&
            iTetrahedron2_I=(/1, 3, 7, 4/),&
            iTetrahedron3_I=(/1, 5, 7, 6/),&
            iTetrahedron4_I=(/1, 7, 6, 4/),&
            iTetrahedron5_I=(/8, 7, 6, 4/))
    case(OneFine_)    
       !\
       ! View from the fine point
       !            xy face diag
       !   x-axis      |      y-axis
       !         \     V     /
       !          \   C4    /      Connections to opposite faces: main diag+
       !           C2/ | \C3       3478 - y edge, yz face diag, xy face daig
       !           |\  | / |       2468 - x edge, xy face diag, xz face diag
       !           | \ |/  |       5687 - z face. yz face diag, xy face diag
       !           |  \/   |
       !           |   F1  |
       !  xz facediag/ | \yz face diag
       !           /   |  \|
       !          C6   |   C7
       !             \ C5/
       !               |
       !               V
       !               z-axis
       call interpolate_pyramids(&
            iRectangular1_I=(/2, 4, 6, 8, 1/) ,&
            iRectangular2_I=(/3, 4, 7, 8, 1/)  ,&
            iRectangular3_I=(/5, 6, 7, 8, 1/))
    case(OneCoarse_)                  ! 8 cases , totally 79 left 176
       !\
       ! View from the coarse point
       !            xy face diag
       !   x-axis            y-axis
       !         \          /
       !          \    F4  /           
       !           F2----- F3            
       !           |\    / |            
       !           | \  / /|            
       !           |\ \/   |            
       !           |   C1  |
       !           | \ | / |
       !           |   |   |
       !          F6  \ / F7
       !             \ F5/
       !               |
       !               V
       !               z-axis
       call interpolate_pyramids(iTetrahedron1_I=(/1,2,3,5/))
       if(nGridOut>1)RETURN
       call interpolate_on_parallel_rays(&
            Dir_D=XyzGrid_DI(:,8) - XyzGrid_DI(:,1), &
            iURectangle1_I=(/2, 4, 6, 8/),&
            iURectangle2_I=(/3, 4, 7, 8/),&
            iURectangle3_I=(/5, 6, 7, 8/),&
            iDTriangle1_I=(/2, 3, 5/),    &
            iDTriangle2_I=(/2, 3, 4/),    &
            iDTriangle3_I=(/2, 6, 5/),    &
            iDTriangle4_I=(/3, 5, 7/))
    case(FineMainDiag_)
       !Full traingulation of all faces, each coarse vertex is connected  
       !by edge with one fine point and with face diagonal to the other. 
       !\
       ! View from the fine point
       !            xy face diag
       !   x-axis      |      y-axis
       !         \     V     /
       !          \   C4    /           
       !           C2/ | \C3           The point F8 is also fine 
       !           |\  | / |           The line of sight goes along the
       !           | \ |/  |           main diagonal F1F8 which is the 
       !           |  \/   |           common side  of 6 tetrahedron
       !           |   F1  |
       !  xz facediag/ | \yz face diag
       !           /   |  \|
       !          C6   |   C7
       !             \ C5/
       !               |
       !               V
       !               z-axis
       call interpolate_pyramids(&
            iTetrahedron1_I=(/1, 2, 4, 8/), &
            iTetrahedron2_I=(/1, 2, 6, 8/), &
            iTetrahedron3_I=(/1, 3, 4, 8/), &
            iTetrahedron4_I=(/1, 3, 7, 8/), &
            iTetrahedron5_I=(/1, 5, 6, 8/), &
            iTetrahedron6_I=(/1, 5, 7, 8/))
    case(FineFaceDiag_)
       !
       !       ^
       !iDir-axis
       !       |  7C---------8C
       !       | /|         /|
       !       |/ |        / |
       !       5C---------6C |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3C---------4F
       !       | /        | /
       !       |/         |/
       !       1F---------2C-----> 1+mod(iDir,3)-axis
       !
       !
       ! 1 Remove tetrahedra 1F 4F 2C 6C and 1F 4F 3C 7C
       call interpolate_pyramids(&
            iTetrahedron1_I=(/1, 4, 2, 6/), &
            iTetrahedron2_I=(/1, 4, 3, 7/))
       if(nGridOut > 1)RETURN
       !
       !    7C------ 8C 
       !  5C--\---6C/|     
       !   |  _ \_  \|
       !   | /    --4F
       !  1F----/     
       !\
       ! Center of the lower face
       !/
       XyzMin_D = 0.50*(XyzGrid_DI(:,1) + XyzGrid_DI(:,4))
       !\
       ! Center of the upper face
       !/
       XyzMax_D = 0.50*&
            (XyzGrid_DI(:,5)+XyzGrid_DI(:,8))
       call interpolate_on_parallel_rays(Dir_D=XyzMax_D - XyzMin_D,&
            iURectangle1_I=(/5, 6, 7, 8/),&
            iDTriangle1_I=(/1, 4, 6/),&
            iDTriangle2_I=(/1, 4, 7/),&
            iDTriangle3_I=(/4, 6, 8/),&
            iDTriangle4_I=(/4, 7, 8/),&
            iDTriangle5_I=(/1, 5, 7/),&
            iDTriangle6_I=(/1, 5, 6/))
    case(CoarseMainDiag_)  
       !Two coarse pints across the main diagonal
       !\
       ! View from the coarse point
       !            xy face diag
       !   x-axis            y-axis
       !         \          /
       !          \    F4  /           
       !           F2----- F3            
       !           |\    / |            
       !           | \  / /|            
       !           |\ \/   |            
       !           |   C1  |
       !           | \ | / |
       !           |   |   |
       !          F6  \ / F7
       !             \ F5/
       !               |
       !               V
       !               z-axis
       call interpolate_pyramids(&
            iTetrahedron1_I=(/1, 2, 3, 5/),&
            iTetrahedron2_I=(/8, 7, 6, 4/))
       if(nGridOut>1)RETURN
       call interpolate_on_parallel_rays(&
            Dir_D=XyzGrid_DI(:,8) - XyzGrid_DI(:,1), &
            iUTriangle1_I=(/7, 6, 4/),&
            iUTriangle2_I=(/7, 6, 5/),&
            iUTriangle3_I=(/7, 4, 3/),&
            iUTriangle4_I=(/6, 4, 2/),&
            iDTriangle1_I=(/2, 3, 5/),&
            iDTriangle2_I=(/2, 3, 4/),&
            iDTriangle3_I=(/2, 5, 6/),&
            iDTriangle4_I=(/3, 5, 7/))
    case(CoarseFaceDiag_)
       !       ^
       !iDir-axis
       !       |  7F---------8F
       !       | /|         /|
       !       |/ |        / |
       !       5F---------6F |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3F---------4C
       !       | /        | /
       !       |/         |/
       !       1C---------2F-----> 1+mod(iDir,3)-axis
       !
       !
       ! 1 Remove tetrahedra 1C 2F 3F 5F and 4C 2F 3F 8F
       call interpolate_pyramids(&
            iTetrahedron1_I=(/1, 2, 3, 5/),&
            iTetrahedron2_I=(/4, 2, 3, 8/))
       if(nGridOut < 1)then
          !
          !    7F------ 8F 
          !  5F------6F//     
          !   \|      | |
          !    3F_    |/
          !        \_ 2F
          !\
          ! Center of the lower face
          !/
          XyzMin_D = 0.50*&
               (XyzGrid_DI(:,1)+XyzGrid_DI(:,4))
          !\
          ! Center of the upper face
          !/
          XyzMax_D = 0.50*&
               (XyzGrid_DI(:,5) + XyzGrid_DI(:,8))
          call interpolate_on_parallel_rays(Dir_D=XyzMax_D - XyzMin_D,&
               iURectangle1_I=(/5, 6, 7, 8/),&
               iDTriangle1_I=(/2, 3, 5/),&
               iDTriangle2_I=(/2, 3, 8/))
       end if
    case(FineEdgePlusOne_)
       !       ^
       !iDir-axis
       !       |  7C---------8F
       !       | /|         /|
       !       |/ |        / |
       !       5C---------6C |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3C---------4F
       !       | /        | /
       !       |/         |/
       !       1F---------2C-----> 1+mod(iDir,3)-axis
       !
       !                    F4 F8
       !         Trapezoid: C3 C7
       !                
       !         Trapezoid: F4 F8
       !                    C2 C6
       !     F1 is apex of trapezoidal pyramids 
       !     Tetrahedra left F1C5F8C6 + F1C5F8C7
       call interpolate_pyramids(&
            iTetrahedron1_I=(/1, 5, 8, 6/),&
            iTetrahedron2_I=(/1, 5, 8, 7/),&
            iTrapezoidal1_I=(/2, 6, 4, 8, 1/),&
            iTrapezoidal2_I=(/3, 7, 4, 8, 1/))
    case(ThreeFineOnFace_)
       !
       !       ^
       !iDir-axis
       !       |  7C---------8C
       !       | /|         /|
       !       |/ |        / |
       !       5C---------6C |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3F---------4C
       !       | /        | /
       !       |/         |/
       !       1F---------2F-----> 1+mod(iDir,3)-axis
       !
       !                  F1 F3
       !       Trapezoid: C5 C7       
       !               
       !       Trapezoid: F1 F2
       !                  C5 C6

       call interpolate_pyramids(&
            iTetrahedron1_I=(/2, 3, 8, 4/))
       if(nGridOut> 1)RETURN
       !\
       !Interpolation in parralel rays
       !/
       !The view for the lower subfaces is as follows
       !  C5----------C6
       !  | \         |
       !  |   F1---F2 |
       !  |   |   /   |
       !  !   |  /  \ | 
       !  |   F3      |
       !  | /     \   |
       !  C7----------C8      
       call interpolate_on_parallel_rays(Dir_D=&
            XyzGrid_DI(:,8) - &
            XyzGrid_DI(:,4),&
            iDTriangle1_I=(/1, 2, 3/),&
            iDTriangle2_I=(/8, 2, 3/),&
            iDTriangle3_I=(/8, 7, 3/),&
            iDTriangle4_I=(/8, 2, 6/),&
            iDTrapezoid1_I=(/5, 7, 1, 3/),&
            iDTrapezoid2_I=(/5, 6, 1, 2/),&
            iURectangle1_I=(/5, 6, 7, 8/))    
    case(CoarseEdgePlusOne_)
       !
       !       ^
       !iDir-axis
       !       |  7F---------8C
       !       | /|         /|
       !       |/ |        / |
       !       5F---------6F |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3F---------4C
       !       | /        | /
       !       |/         |/
       !       1C---------2F-----> 1+mod(iDir,3)-axis
       !
       !
       !                  F3 F7
       !       Trapezoid: C4 C8       
       !         
       !       Trapezoid: F2 F6
       !                  C4 C8
       !                
       call interpolate_pyramids(&
            iTetrahedron1_I=(/1, 2, 3, 5/),&
            iRectangular1_I=(/2, 3, 6, 7, 5/))
       if(nGridOut> 1)RETURN
       call interpolate_on_parallel_rays(Dir_D=&
            XyzGrid_DI(:,3) - &
            XyzGrid_DI(:,2),&
            iDTrapezoid1_I=(/4, 8, 2, 6/),&
            iUTrapezoid1_I=(/4, 8, 3, 7/)) 
    case(ThreeCoarseOnFace_)
       !
       !       ^
       !iDir-axis
       !       |  7F---------8F
       !       | /|         /|
       !       |/ |        / |
       !       5F---------6F |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3C---------4F
       !       | /        | /
       !       |/         |/
       !       1C---------2C-----> 1+mod(iDir,3)-axis
       !
       !
       !            F5 F7
       ! Trapezoid: C1 C3       
       !           
       ! Trapezoid: F5 F6
       !            C1 C2
       ! Rectangle: F5 F6 F7 F8
       !/               
       !    Common apex F4
       call interpolate_pyramids(&
            iRectangular1_I=(/7, 8, 5, 6, 4/),&
            iTrapezoidal1_I=(/1, 2, 5, 6, 4/),&
            iTrapezoidal2_I=(/1, 3, 5, 7, 4/))
    case(ThreeCoarseOnFacePlusOne_)
       !\
       ! face 1:4 has only one Coarse point 1
       ! face 5:8 has three Coarse points iOrder((/6,7,8/))
       !/
       !-----------------------------------------------
       !
       !       ^
       !iDir-axis
       !       |  7C---------8C
       !       | /|         /|
       !       |/ |        / |
       !       5F---------6C |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3F---------4F
       !       | /        | /
       !       |/         |/
       !       1C---------2F-----> 1+mod(iDir,3)-axis
       !
       !                F2 F4
       !      Trapezoid C6 C8
       !                                           
       !                F3 F4
       !      Trapezoid C7 C8
       !                
       !
       !Common apex F5
       call interpolate_pyramids(&
            iTrapezoidal1_I=(/6, 8, 2, 4, 5/),&
            iTrapezoidal2_I=(/7, 8, 3, 4, 5/),&
            iTetrahedron1_I=(/2, 3, 4, 5/),   &
            iTetrahedron2_I=(/1, 2, 3, 5/) )
    case(CoarseChain_)
       !\
       !Chain of four coarse points connected with mutually 
       !orthogonal edges
       !/
       !-----------------------------------------------
       !
       !       ^
       !iDir-axis
       !       |  C7---------8F
       !       | /|         /|
       !       |/ |        / |
       !       5C---------6F |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3F---------4F
       !       | /        | /
       !       |/         |/
       !       1C---------2C-----> 1+mod(iDir,3)-axis
       !
       !                  F6 F8
       !        Trapezoid C5 C7
       !                                           
       !      
       !                  F3 F4
       !        Trapezoid C1 C2                
       !
       call interpolate_pyramids(&
            iTrapezoidal1_I=(/1, 2, 3, 4, 6/),& 
            iTetrahedron1_I=(/8, 3, 4, 6/),&  !see #1 below 
            iTetrahedron2_I=(/1, 6, 3, 5/),&  !see #2 below
            iTrapezoidal2_I=(/5, 7, 6, 8, 3/))!see #3 below
       ! #1 The leftover is above subfaces 136 and 364
       ! #2 The leftover is above subfaces 136 and 368
       ! #3 The leftover is above subfaces 365 and 368
       !==========================
    end select
  contains
    subroutine interpolate_pyramids(&
         iTetrahedron1_I,iTetrahedron2_I,iTetrahedron3_I,&
         iTetrahedron4_I,iTetrahedron5_I,iTetrahedron6_I,&
         iRectangular1_I,iRectangular2_I,iRectangular3_I,&
         iTrapezoidal1_I,iTrapezoidal2_I)
      use ModInterpolateSimpleShape, ONLY:&
           gen_interpolate_pyramids=>interpolate_pyramids
      integer,intent(in),optional,dimension(4)::&
           iTetrahedron1_I, iTetrahedron2_I, iTetrahedron3_I,&
           iTetrahedron4_I, iTetrahedron5_I, iTetrahedron6_I
      integer,intent(in),optional,dimension(5)::&
           iRectangular1_I, iRectangular2_I, iRectangular3_I,&
           iTrapezoidal1_I, iTrapezoidal2_I
      integer:: iOrderHere_I(1:5)
      !---------------------------
      call gen_interpolate_pyramids(&
           XyzGrid_DI, Xyz_D,iOrderHere_I, Weight_I, nGridOut,&
           iTetrahedron1_I,iTetrahedron2_I,iTetrahedron3_I,&
           iTetrahedron4_I,iTetrahedron5_I,iTetrahedron6_I,&
           iRectangular1_I,iRectangular2_I,iRectangular3_I,&
           iTrapezoidal1_I,iTrapezoidal2_I)
      if(nGridOut > -1)iOrder_I(1:nGridOut) = &
           iOrder_I(iOrderHere_I(1:nGridOut))
    end subroutine interpolate_pyramids
    !======================
    subroutine interpolate_on_parallel_rays(Dir_D, &
         iDRectangle1_I, iDTrapezoid1_I, iDTrapezoid2_I, &
         iDTriangle1_I, iDTriangle2_I, iDTriangle3_I,&
         iDTriangle4_I, iDTriangle5_I, iDTriangle6_I,&
         iURectangle1_I, iURectangle2_I, iURectangle3_I,&
         iUTrapezoid1_I,                             &
         iUTriangle1_I, iUTriangle2_I, iUTriangle3_I,&
         iUTriangle4_I)
      use ModInterpolateSimpleShape, ONLY:&
           parallel_rays=>interpolate_on_parallel_rays
      !Direction of parallel rays
      real, intent(in) :: Dir_D(nDim)
      !\
      !Up subfaces ure in the direction of  Dir_D from Xyz point
      !Down faces are in the direction of -Dir_D
      !/
      integer, intent(in), optional, dimension(3)::&
           iDTriangle1_I, iDTriangle2_I, iDTriangle3_I,&
           iDTriangle4_I, iDTriangle5_I, iDTriangle6_I,&
           iUTriangle1_I, iUTriangle2_I, iUTriangle3_I,&
           iUTriangle4_I
      integer, intent(in), optional, dimension(4)::&
           iDRectangle1_I, iDTrapezoid1_I, iDTrapezoid2_I, iUTrapezoid1_I,&
           iURectangle1_I, iURectangle2_I, iURectangle3_I
      !\
      !Misc
      !/
      integer:: iOrderHere_I(nGrid)
      call parallel_rays(&
           XyzGrid_DI, Xyz_D, iOrderHere_I, Weight_I, nGridOut, Dir_D, &
           iDRectangle1_I, iDTrapezoid1_I, iDTrapezoid2_I,  &
           iDTriangle1_I, iDTriangle2_I, iDTriangle3_I,&
           iDTriangle4_I, iDTriangle5_I, iDTriangle6_I,&
           iURectangle1_I, iURectangle2_I, iURectangle3_I,&
           iUTrapezoid1_I,                             &
           iUTriangle1_I, iUTriangle2_I, iUTriangle3_I,&
           iUTriangle4_I)
      if(nGridOut > 0)iOrder_I(1:nGridOut) =  &
           iOrder_I(iOrderHere_I(1:nGridOut))
    end subroutine interpolate_on_parallel_rays
    !===========================
  end subroutine resolution_corner
end module ModResolutionCorner
!=========================
module ModInterpolateAMR
  !\
  !Generalize bilinear and trilinear interpolation for AMR grids
  !The data are given at the cell-centered grid which consists of AMR blocks.
  !The interpolation is free of any jumps at the resolution interfaces
  !including edges and corners
  !/
  !\
  !USE
  !/
  !=================ARRAYS FOR A CUBE====================================!
  !For a cubic stencil enumerated as follows:
  !
  !       ^
  ! z-axis|
  !       |  7----------8
  !       | /|         /|
  !       |/ |        / |
  !       5----------6  |
  !       |  |   _   |  |
  !       |  |   /|y-axi$
  !       |  |  /    |  |
  !       |  | /     |  |
  !       |  |/      |  |
  !       |  3----------4
  !       | /        | /
  !       |/         |/
  !       1----------2------> x-axis
  !       
  ! we provide several functions characterizing its geometry in terms of 
  ! the grid point numbers. 
  use ModCubeGeometry, ONLY: iFace_IDI
  !Vertexes (enumerated by the first undex), which form 
  !the face of direction iDir (second index) including the 
  !given vertex (the third index). 
  !When the first index equals 1,2,3,4, the vertex, accordingly: 
  !(1) coincides with the given one
  !(2) is connected with the given one by the edge of direction iDir+1
  !(3) is connected with the given one by the edge of direction iDir+2
  !(4) is connected to the given one by the face diagonal of direction iDir

  use ModCubeGeometry, ONLY: iOppositeFace_IDI
  !Vertexes (enumerated by the first undex), which form 
  !the face of direction iDir (second index) and does not 
  !include the given vertex (the third index). 
  !When the first index equals 1,2,3,4, the vertex, accordingly: 
  !(1) is connected to the given one by the edge of direction iDir
  !(2) is connected to (1) by the edge of direction iDir+1
  !(3) is connected to (1) by the edge of direction iDir+2
  !(4) is connected to  by the face diagonal of direction iDir

  use ModCubeGeometry, ONLY: iEdge_ID  
  !Number of the vertex connected by 
  !the edge of direction iDir (second index) 
  !with the given vertex (first index)

  !=================ARRAYS FOR A RECTANGLE==============!
  !For a rectangular stencil enumerated as follows:
  !
  !       ^
  ! y-axis|
  !       |   
  !       3----------4  
  !       |          |  
  !       |          !
  !       |          |  
  !       |          |  
  !       |          |  
  !       |          !
  !       |          |
  !       |          |
  !       1----------2------> x-axis
  !       
  ! we provide several functions characterizing its geometry in terms of 
  ! the grid point numbers. 

  use ModCubeGeometry, ONLY: iSide_IDI
  !Vertexes (enumerated by the first undex), which form 
  !the side of direction iDir (second index) including the 
  !given vertex (the third index). 
  !When the first index equals 1,2,3,4, the vertex, accordingly: 
  !(1) coincides with the given one
  !(2) is connected with the given one by the edge of direction iDir

  use ModCubeGeometry, ONLY: iOppositeSide_IDI
  !Vertexes (enumerated by the first undex), which form 
  !the side of direction iDir (second index) and does not include the
  !given vertex (the third index). 
  !When the first index equals 1,2 the vertex, accordingly: 
  !(1) is connected to the given one by the side of direction 1+mod(iDir,2)
  !(2) is connected to (1) by the edge of direction iDir
  !---------------------REFINEMENT CHARACTERSTICS------------
  use ModCubeGeometry, ONLY: Coarse_, Fine_

  !\
  ! Different cases of 3D stencil 
  !/
  use ModCubeGeometry, ONLY: Uniform_, Face_, Edge_, OneFine_, OneCoarse_,& 
       Transition2Edge_, Transition2Corner_, TransitionJunction_
  !\
  ! Analogous for 2D
  !/
  use ModCubeGeometry, ONLY: Trapezoid_, Rhombus_ 
  !\
  ! Array to sort stencil
  !/
  use ModCubeGeometry, ONLY: iSortStencil3_II, Case_, Dir_, Grid_, i_case
  !\
  !Array used to sort stencil: 
  !(1) assign a type to it (Case_) 
  !(2) assign a grid point of the stencil chosen as the basic one.
  !(Grid_) example, if a a single point in the stencil is coarse and all
  !other are fine, the 'basic point' is the coarse point.
  !(3) assign a direction, if the orientation if this sort of stencil
  !has a characteristic direction (Dir_). For example, for the 'resolution 
  !edge' of y direction iDir equals 1.
  !/
  implicit none
  PRIVATE !Except
  SAVE
  integer, parameter:: x_ = 1, y_ = 2, z_ = 3
  !\
  ! To improve the algorithm stability against roundoff errors
  !/
  real, parameter:: cTol = 0.00000010
  !============================================================================
  !Interpolation on the block AMR grid
  !\
  !Calculates interpolation weights
  !/
  !\
  !Example of application for SERIAL calculation of the
  !interpolated value of the state vector sampled in the 
  !grid points as
  !State_VGB(nVar,nI,nJ,nK,nBlock) array, where nI=nJ=nK=4
  !in point Xyz_D looks as follows:
  !
  !call interpolate_amr(&
  !  nDim=3,              &!number of dimensions
  !  XyzIn_D=,Xyz_D,      &!Point in which to interpolate
  !  nIndexes = 4,        &!Three cell indexes plus one block index
  !  find = find_subroutine , &! Search in grid
  !  nCell_D  = (/4,4,4/) &!
  !  nGridOut=nGridOut,   &! Number of points in the output stencil
  !  Weight_I=Weight_I,   &! Weight coefficients
  !  iIndexes_II=iIndexes_II) !Cell+block indexes to be used in interpolation
  !
  !  if(nGridOut < 1) call CON_stop('Interpolation failed')
  !  Value_V(1:nVar) = 0
  !  do iGrid = 1, nGridOut
  !     Value_V = Value_V + &
  !       State_VGB(:,iIndexes_II(1,iGrid), iIndexes_II(2,iGrid), &
  !             iIndexes_II(3,iGrid), iIndexes_II(4,iGrid))*&
  !                               Weight_I(iGrid)
  !  end do
  ! For PARALLEL code the processor number at which the corresponding part
  ! of  State_VGB is allocated is provided in iIndex_II(0,:) components
  ! of the index output array
  !/
  public interpolate_amr
contains
  !=================================
  subroutine interpolate_amr(nDim, XyzIn_D, nIndexes, find, nCell_D,  &
       nGridOut, Weight_I, iIndexes_II, IsSecondOrder)
    use ModResolutionCorner, ONLY: resolution_corner
    use ModKind, ONLY: nByteReal
    !\
    ! INPUT PARAMETERS
    !/
    !\
    ! Number of dimensions
    !/
    integer, intent(in) :: nDim 
    !\
    ! Number of indexes. Usually, nIndexes = nDim + 1, for three cell indexes 
    ! and a block number. If the latter information is not needed, 
    ! nIndexes=nDim should be added
    !/
    integer, intent(in) :: nIndexes
    !\
    ! Point coordinates
    !/
    real,    intent(in) :: XyzIn_D(nDim)
    !\
    ! Block AMR Grid characteristic: number of cells
    !/ 
    integer, intent(in) :: nCell_D(nDim)
    !\
    ! Yet another Block AMR Grid characteristic:
    ! the search routine, which returns, for a given point,
    ! the block to which this points belong and the processor, at which
    ! this block is allocated, as well as the block parameters:
    ! coordinates of the left corner (the point in the block with the 
    ! minvalue of coordinates) and (/Dx, Dy, Dz/)
    !/
    interface 
       subroutine find(nDim, Xyz_D, &
            iProc, iBlock, XyzCorner_D, Dxyz_D, IsOut)
         implicit none
         integer, intent(in) :: nDim
         !\
         ! "In"- the coordinates of the point, "out" the coordinates of the
         ! point with respect to the block corner. In the most cases
         ! XyzOut_D = XyzIn_D - XyzCorner_D, the important distinction,
         ! however, is the periodic boundary, near which the jump in the
         ! stencil coordinates might occur. To handle the latter problem,
         ! we added the "out" intent. The coordinates for the stencil
         ! and input point are calculated and recalculated below with
         ! respect to the block corner. 
         !/
         real,  intent(inout):: Xyz_D(nDim)
         integer, intent(out):: iProc, iBlock !processor and block number
         !\
         ! Block left corner coordinates and the grid size:
         !/
         real,    intent(out):: XyzCorner_D(nDim), Dxyz_D(nDim)
         logical, intent(out):: IsOut !Point is out of the domain.
       end subroutine find
    end interface
    !\
    !OUTPUT PARAMETERS
    !/
    !\
    !Number of grid points involved into interpolation
    !/
    integer, intent(out):: nGridOut      
    !\
    ! Interpolation weights (only the first nGridOut values are meaningful
    !/
    real,    intent(out):: Weight_I(2**nDim) 
    !\
    ! Cell(+block) indexes and processor number for grid points to be
    ! invilved into interpolation. iProc numbers are stored in 
    ! iIndexes_II(0,:) 
    !/
    integer, intent(out):: iIndexes_II(0:nIndexes,2**nDim)

    !\
    !The following is true if stencil does not employ 
    !the out-of-grid points
    !/
    logical, intent(out), optional:: IsSecondOrder  
  
    !\
    ! Coordinates of the input point may be recalculated within
    ! routine; inverse of (/Dx, Dy, Dz/) is reused
    !/
    real, dimension(nDim)      :: Xyz_D, DxyzInv_D, XyzStencil_D
    !\
    ! The extended stencil in a structured form:
    ! A cubic 2*2*2 grid with 2*2*2 subgrids covering each vertex
    !/
    real       :: XyzGrid_DII(nDim,0:2**nDim,2**nDim)
    integer, dimension(2**nDim):: iBlock_I  , iProc_I 
    integer                    :: iCellIndexes_DII(nDim,2**nDim,2**nDim)
    integer, dimension(2**nDim):: iLevelSubgrid_I
    logical, dimension(2**nDim):: IsOut_I

    real, dimension(nDim) :: XyzBasicBlock_D, DxyzBasicBlock_D
    logical               :: IsOutOfDomain
    !\
    !Just 2**nDim
    integer:: nGrid
    !/
    !\
    ! To improve the algorithm stability against roundoff errors
    !/
    real           :: cTol2
    !\
    ! Shift of the iGrid point in the stencil with respect to the
    ! first one
    !/
    integer, dimension(3,8),parameter :: iShift_DI = reshape((/&
         0, 0, 0,   1, 0, 0,   0, 1, 0,   1, 1, 0, &
         0, 0, 1,   1, 0, 1,   0, 1, 1,   1, 1, 1/),(/3,8/)) 

    integer:: iGridOutOfBlock
    !------------------------
    cTol2 = cTol**(nByteReal/4)

    nGrid = 2**nDim !Number of points in a basic stencil

    !\
    ! Initialize 
    !/
    iIndexes_II = 0; Weight_I    = 0
    nGridOut = -1; Xyz_D = XyzIn_D ; IsOut_I = .false.
    if(present(IsSecondOrder))IsSecondOrder = .false.
    !\
    ! Implemented version:
    ! only physical points are involved into interpolation,
    ! but, in general case, from many blocks. An alternative 
    ! version, when the ghost cells may be involved but the 
    ! stencil points should be all from one block is under
    ! development.
    !/
    !\
    ! Find block to which the point belong
    !/ 
    call find(nDim, Xyz_D, iProc_I(1), iBlock_I(1), &
         XyzBasicBlock_D, DxyzBasicBlock_D, IsOutOfDomain)
    if(IsOutOfDomain)then
       !\
       ! The algorithm does not work for a point out of the computation
       ! domain. It could, but in this case too much information
       ! about grid should be brought to the table - the domain size,
       ! periodicity etc
       !/
       nGridOut = -1
       RETURN
    end if
    !\
    ! Now Xyz_D is given  with respect to the main block corner
    !/
    call get_main_block(iGridOutOfBlock, IsFirstCall = .true.)
    !\
    !The interpolation is done, if the stencil is within a single block
    !/
    if(nGridOut > 0) then
       if(present(IsSecondOrder))IsSecondOrder = .true.
       RETURN
    end if
    call get_other_blocks(iGridOutOfBlock)
    call interpolate_extended_stencil(nDim, Xyz_D, nIndexes, &
         XyzGrid_DII, iCellIndexes_DII, iBlock_I, iProc_I,   &
         iLevelSubgrid_I, IsOut_I,DxyzInv_D,                 &
         nGridOut, Weight_I, iIndexes_II, IsSecondOrder)
contains
    subroutine get_main_block(iGridOutOfBlock,IsFirstCall)
      integer, intent(out):: iGridOutOfBlock
      logical,optional,intent(in):: IsFirstCall
      !\
      ! Fills in XyzGrid_DII(:,0,:) - coarse grid
      ! Fills in indexes for the points of the stencil 
      ! belonging to this block. Returns the maximum
      ! number of the grid poit which is out of the main
      ! block
      !/
      !\
      ! Loop variable
      !/
      integer:: iGrid
      !\
      !Displacement measured in grid sizes or in their halfs
      !/ 
      integer, dimension(nDim) :: iShift_D
      !\
      ! Misc
      !/
      real :: XyzMisc_D(nDim)
      !------------------------------------
      iLevelSubgrid_I =  0
      iCellIndexes_DII = 0; XyzGrid_DII     = 0
      DxyzInv_D = 1/DxyzBasicBlock_D
      XyzMisc_D = Xyz_D*DxyzInv_D + 0.50
      iCellIndexes_DII(:,1,1) = floor(XyzMisc_D)
      !\
      !Calculate coordinates of the left corner of a stencil
      !/
      XyzMisc_D = XyzMisc_D - iCellIndexes_DII(:,1,1)
      if(present(IsFirstCall))then
         if(all(iCellIndexes_DII(:,1,1) > 0.and.&
              iCellIndexes_DII(:,1,1) < nCell_D))then
            !\
            ! The whole interpolation stencil is within this block
            !/
            call interpolate_uniform(nDim, XyzMisc_D, Weight_I)
            !\
            ! Form index array and sort out zero weights
            !/
            nGridOut = 0 ; cTol2 = 2*cTol2
            do iGrid = 1, nGrid
               if(Weight_I(iGrid) < cTol2)CYCLE
               nGridOut = nGridOut + 1
               iIndexes_II(0,       nGridOut) = iProc_I(1)
               iIndexes_II(nIndexes,nGridOut) = iBlock_I(1)
               iIndexes_II(1:nDim,  nGridOut) = iCellIndexes_DII(:,1,1) + &
                    iShift_DI(1:nDim,iGrid)
               Weight_I(nGridOut) = Weight_I(iGrid)
            end do
            RETURN  !All interpolation is done, ready to exit
            ! elseif(all(abs(XyzMisc_D) < cTol2))then
            !    !\
            !    ! Xyz coincides with the grid point
            !    ! Commented out to satisfy iFort
            !    !/
            !    nGridOut = 1; Weight_I =0; Weight_I(1) = 1
            !    iIndexes_II(0,       1) = iProc_I(1)
            !    iIndexes_II(nIndexes,1) = iBlock_I(1)
            !    iIndexes_II(1:nDim,  1) = iCellIndexes_DII(:,1,1)
            !    RETURN
         end if
      end if
      XyzGrid_DII(:,0,1) = DxyzBasicBlock_D*(iCellIndexes_DII(:,1,1) - 0.50)
      !\ 
      ! now XyzMisc_D = (Xyz_D-XyzGrid_DII(:,0,1))/Dxyz satisfies 
      ! inequalities: XyzMisc_D >= 0 and XyzMisc_D < 1. Strengthen 
      ! these inequalities 
      !/
      XyzMisc_D = min(1 - cTol2,&
           max(XyzMisc_D, cTol2 ))
      Xyz_D = XyzGrid_DII(:,0,1) + XyzMisc_D*DxyzBasicBlock_D
      !\
      !Calculate other grid points, check if all points belong to 
      !the found block
      !/
      iGridOutOfBlock = -1
      do iGrid = nGrid, 1, -1
         iShift_D = iShift_DI(1:nDim,iGrid)
         iBlock_I(iGrid) = iBlock_I(1)
         iCellIndexes_DII(:,1,iGrid) = &
              iCellIndexes_DII(:,1,1) + iShift_D
         XyzGrid_DII(:,0,iGrid) = &
              XyzGrid_DII(:,0,1) + iShift_D*DxyzBasicBlock_D
         if(any(iCellIndexes_DII(:,1,iGrid) < 1).or.&
              any(iCellIndexes_DII(:,1,iGrid) > nCell_D))then
            !\
            !This grid point is out of block, mark it
            !/
            iProc_I(iGrid) = -1 
            iGridOutOfBlock = max(iGridOutOfBlock,iGrid)
         else
            XyzGrid_DII(:,1,iGrid) = XyzGrid_DII(:,0,iGrid)
            iProc_I(iGrid) = iProc_I(1)
         end if
      end do
    end subroutine get_main_block
    !=====================
    recursive subroutine get_other_blocks(iGridOutOfBlock)
      integer, intent(inout)::iGridOutOfBlock
      !\
      ! Misc
      !/
      integer:: iGrid, iGridStored
      !\
      ! For using find routine with inout argument
      !/
      real  :: XyzMisc_D(nDim)
      !\
      ! Output parameters of find routine
      !/
      real    :: XyzCorner_D(nDim), Dxyz_D(nDim)
      logical :: IsOut
      !----------------------------
      if(iGridOutOfBlock == -1)RETURN
      iGridStored = iGridOutOfBlock
      !\
      ! For the grid point not belonging to the block
      ! find the block they belong to
      !/ 
      !\
      !Recalculate absolute coordinates for
      !the grid point which is out of the block
      !/
      XyzMisc_D = XyzBasicBlock_D + XyzGrid_DII(:,0,iGridStored)
      !\
      ! Find neighboring block
      !/
      call find(nDim, XyzMisc_D, &
           iProc_I(iGridStored), iBlock_I(iGridStored), &
           XyzCorner_D, Dxyz_D, IsOut)
      if(IsOut)then
         iProc_I(iGridStored) = 0 !For not processing this point again
         IsOut_I(iGridStored) = .true.
         XyzGrid_DII(:,1,iGridStored) = XyzGrid_DII(:,0,iGridStored)
         !\
         ! Find the next out-of-block point
         !/
         do iGrid = iGridStored - 1, 1, -1
            if(iProc_I(iGrid)==-1)then
               iGridOutOfBlock = iGrid
               call get_other_blocks(iGridOutOfBlock)
               RETURN
            end if
         end do
         RETURN
      end if
      iLevelSubgrid_I(iGridStored) = 1 - &
           floor(Dxyz_D(1)*DXyzInv_D(1)+ cTol)
      !\                     ^
      ! For expression above | equal to 2 , 1, 0.5 correspondingly
      ! iLevel = -1, 0, Fine_, meaning that the neighboring block
      ! is coarser, at the same resolution or finer than the basic 
      ! one.
      !/
      select case(iLevelSubgrid_I(iGridStored))
      case(-1)
         !The neighboring block is coarser and should be used as the 
         !stencil base. Now Xyz_D and XyzGrid_DII(:,0,iGridStored) are
         !defined with respect to the original block, the latter point
         !in the coarser neighboring block has a coordinates XyzMisc_D
         !So that
         Xyz_D = Xyz_D - XyzGrid_DII(:,0,iGridStored) + XyzMisc_D
         !\
         ! Return to the beginning of algorithm. Before the COARSEN
         ! loop the block and processor number for the basic block
         ! were stored in iProc_I(1), iBlock_I(1)
         !/
         iProc_I( 1) = iProc_I(iGridStored)
         iBlock_I(1) = iBlock_I(iGridStored)
         XyzBasicBlock_D = XyzCorner_D
         DxyzBasicBlock_D = Dxyz_D
         call get_main_block(iGridOutOfBlock)
      case(0)  ! (New Dxyz_D)*Stored DXyzInv =1
         call get_block(iGridOutOfBlock, XyzMisc_D, Dxyz_D)
      case(Fine_  )  !1, (New Dxyz_D)*Stored DXyzInv =0.5 
         call get_fine_block(iGridOutOfBlock, XyzMisc_D, Dxyz_D)
      end select
      call get_other_blocks(iGridOutOfBlock)
    end subroutine get_other_blocks
    !========================
    subroutine get_block(iGridOutOfBlock, Xyz_D, Dxyz_D)
      integer, intent(inout)::iGridOutOfBlock
      real, dimension(nDim), intent(in):: Xyz_D, Dxyz_D
      !\
      ! Fills in the indexes for the grid boints belonging
      ! to the block, which it at the same resolution as the
      ! main block Returns the maximum number of the grid poit
      ! which is out of all blocks  found so far
      !/
      !\
      ! Misc
      !/
      integer:: iGrid, iGridStored
      !\
      !Displacement measured in grid sizes or in their halfs
      !/ 
      integer, dimension(nDim) :: iShift_D
      !---------------------
      iGridStored = iGridOutOfBlock
      XyzGrid_DII(:,1,iGridStored) = XyzGrid_DII(:,0,iGridStored)
      !\
      ! Calculate cell indexes as we did before. Use nint
      ! instead of int as long as XyzMisc_D is very close to 
      ! the grid point 
      !/
      iCellIndexes_DII(:,1,iGridStored) = nint(Xyz_D*DxyzInv_D + 0.50)
      !\
      ! Check if there are more grid points belonging to the
      ! newly found block
      !/
      iShift_D = iShift_DI(1:nDim,iGridStored)
      iGridOutOfBlock = -1
      do iGrid = iGridStored - 1, 1, -1
         if(iProc_I(iGrid)/=-1)CYCLE !This point is done earlier
         iCellIndexes_DII(:,1,iGrid) = &
              iCellIndexes_DII(:,1,iGridStored) +&
              iShift_DI(1:nDim,iGrid) - iShift_D 
         !\
         ! Check if the point is in the newly found block
         !/ 
         if(any(iCellIndexes_DII(:,1,iGrid) < 1).or.&
              any(iCellIndexes_DII(:,1,iGrid) > nCell_D))then
            !\
            !This grid point is out of block, mark it for further work
            !/
            iGridOutOfBlock = max(iGridOutOfBlock, iGrid)
         else
            iProc_I( iGrid) = iProc_I( iGridStored)
            iBlock_I(iGrid) = iBlock_I(iGridStored)
            XyzGrid_DII(:,1,iGrid) = XyzGrid_DII(:,0,iGrid)
            iLevelSubGrid_I(iGrid) = 0
         end if
      end do
    end subroutine get_block
    !=====================
    subroutine get_fine_block(iGridOutOfBlock, Xyz_D, Dxyz_D)
      integer, intent(inout)::iGridOutOfBlock
      real, dimension(nDim), intent(in):: Xyz_D, Dxyz_D
      !\
      ! Fills in the indexes for the grid boints belonging
      ! to the block, which it at the same resolution as the
      ! main block Returns the maximum number of the grid poit
      ! which is out of all blocks  found so far
      !/
      !\
      ! Loop variables
      !/
      integer:: iGrid, iSubGrid
      !\
      !Displacement measured in grid sizes or in their halfs
      !/ 
      integer, dimension(nDim) :: iShift_D
      !/
      integer :: iGridStored
      !------------------
      iGridStored = iGridOutOfBlock
      !\
      ! Fine subgrid is displaced by half of finer subgrid size
      !/
      XyzGrid_DII(:,1,iGridStored) = XyzGrid_DII(:,0,iGridStored) -&
           0.50*Dxyz_D
      !\
      ! Calculate cell indexes as we did above. Note that DxyzInv_D
      ! is twice less than needed, because it is calculated for the
      ! whole stencil, not for the finer subgrid
      !/
      iCellIndexes_DII(:,1,iGridStored) = &
           floor(2*Xyz_D*DxyzInv_D + 0.50)
      !\
      ! All points in the 2*2*2 finer subgrid are involved
      !/
      do iSubGrid = 2, nGrid
         iShift_D = iShift_DI(1:nDim,iSubGrid)
         XyzGrid_DII(:,iSubGrid,iGridStored) = &
              XyzGrid_DII(:,1,iGridStored) + Dxyz_D*iShift_D
         iCellIndexes_DII(:,iSubGrid,iGridStored) = &
              iCellIndexes_DII(:,1,iGridStored) + iShift_D
      end do
      !\
      ! Check if there are more grid points belonging to the
      ! newly found block
      !/
      iGridOutOfBlock = -1
      do iGrid = iGridStored -1, 1, -1
         if(iProc_I(iGrid)/=-1)CYCLE !This point is done earlier
         iCellIndexes_DII(:,1,iGrid) = &
              iCellIndexes_DII(:,1,iGridStored) + 2*(&
              iShift_DI(1:nDim,iGrid) - iShift_DI(1:nDim,iGridStored))
         if(any(iCellIndexes_DII(:,1,iGrid) < 1).or.&
              any(iCellIndexes_DII(:,1,iGrid) > nCell_D))then
            !\
            !This grid point is out of block
            !/
            iGridOutOfBlock = max(iGridOutOfBlock, iGrid)
         else
            iProc_I(iGrid ) = iProc_I( iGridStored)
            iBlock_I(iGrid) = iBlock_I(iGridStored)
            XyzGrid_DII(:,1,iGrid) = XyzGrid_DII(:,0,iGrid) -&
                 0.50*Dxyz_D
            iLevelSubGrid_I(iGrid) = Fine_
            do iSubGrid = 2, nGrid
               iShift_D = iShift_DI(1:nDim,iSubGrid)
               XyzGrid_DII(:,iSubGrid,iGrid) = &
                    XyzGrid_DII(:,1,iGrid) + Dxyz_D*iShift_D
               iCellIndexes_DII(:,iSubGrid,iGrid) = &
                    iCellIndexes_DII(:,1,iGrid) + iShift_D
            end do
         end if
      end do
    end subroutine get_fine_block
  end subroutine interpolate_amr
  !============================
  subroutine interpolate_uniform(nDim, Dimless_D, Weight_I, IsOut_I)
    integer,intent(in)::nDim
    !\
    !The displacement of the point from the stencil left corner
    !normalized by the grid size. Not exceeding zero and is less than 1.
    !/
    real, intent(in)::Dimless_D(nDim)
    real, intent(out):: Weight_I(2**nDim)
    logical, optional, intent(in):: IsOut_I(2**nDim)
    !\
    !Loop variables
    !/
    !------------------
    Weight_I(1) = (1 - Dimless_D(1))*(1 - Dimless_D(2))
    Weight_I(2) =      Dimless_D(1) *(1 - Dimless_D(2))
    Weight_I(3) = (1 - Dimless_D(1))*     Dimless_D(2)
    Weight_I(4) =      Dimless_D(1) *     Dimless_D(2)
    if(nDim==3)then
       Weight_I(5:8) = Weight_I(1:4)*     Dimless_D(3)
       Weight_I(1:4) = Weight_I(1:4)*(1 - Dimless_D(3))
    end if
    if(present(IsOut_I))then
       if(any(IsOut_I))then
          where(IsOut_I)Weight_I = 0
          Weight_I = Weight_I/sum(Weight_I)
       end if
    end if
  end subroutine interpolate_uniform
  !=================================
  subroutine interpolate_extended_stencil(nDim, Xyz_D, nIndexes, &
    XyzGrid_DII, iCellIndexes_DII, iBlock_I, iProc_I,   &
    iLevelSubgrid_I, IsOut_I, DxyzInv_D,               &
    nGridOut, Weight_I, iIndexes_II, IsSecondOrder)
    use ModResolutionCorner, ONLY: resolution_corner
    !\
    !USE: The tools to fill in sort stencil array
    !/
    use ModCubeGeometry, ONLY: DoInit, init_sort_stencil
    use ModKind, ONLY: nByteReal
    !\
    ! INPUT PARAMETERS
    !/
    !\
    ! Number of dimensions
    !/
    integer, intent(in) :: nDim 
    !\
    ! Number of indexes. Usually, nIndexes = nDim + 1, for three cell indexes 
    ! and a block number. If the latter information is not needed, 
    ! nIndexes=nDim should be added
    !/
    integer, intent(in) :: nIndexes
    !\
    ! Point coordinates
    !/
    !\
    ! Beyond the boundary the grid is prolonged and near the boundary
    ! pooint Xyz may be shifted to physical domain for interpolation.
    ! Therefore, Xyz_D, XyzGrid_DII and iLevelSubGrid_I have intent inout
    ! and their actual values used may be found, if desired. 
    !/
    real,    intent(inout) :: Xyz_D(nDim)
    !\
    ! The extended stencil in a structured form:
    ! A cubic 2*2*2 grid with 2*2*2 subgrids covering fine vertexes
    !/
    real,    intent(inout):: XyzGrid_DII(nDim,0:2**nDim,2**nDim) 
    integer, intent(in):: iCellIndexes_DII(nDim,2**nDim,2**nDim)
    integer, intent(in), dimension(2**nDim):: iBlock_I, iProc_I
    integer, intent(inout):: iLevelSubgrid_I(2**nDim) 
    logical, intent(inout):: IsOut_I(2**nDim)
    real,    intent(in):: dXyzInv_D(nDim)  !The inverse of grid size

    !\
    !OUTPUT PARAMETERS
    !/
    !\
    !Number of grid points involved into interpolation
    !/
    integer, intent(out):: nGridOut      
    !\
    ! Interpolation weights (only the first nGridOut values are meaningful
    !/
    real,    intent(out):: Weight_I(2**nDim) 
    !\
    ! Cell(+block) indexes and processor number for grid points to be
    ! invilved into interpolation. iProc numbers are stored in 
    ! iIndexes_II(0,:) 
    !/
    integer, intent(out):: iIndexes_II(0:nIndexes,2**nDim)
    
    !\
    !The following is true if stencil does not employ 
    !the out-of-grid points
    !/
    logical, intent(out), optional:: IsSecondOrder  
    !\
    ! Local variables
    !/
    !\
    ! The output array in interpolate_amr2,3  routine
    ! This is the set of numbers of the elements of the basic stencil
    ! to be involved in interpolation
    !/
    integer, dimension(2**nDim):: iOrder_I
    !\
    ! Output parameter of interpolate_amr3
    ! If .true., the basic stencil should be re-evaluated
    !/
    logical                    :: IsCorner
    !\
    ! Basic stencil; points:
    !/
    real                       :: XyzGrid_DI(nDim,2**nDim)
    !\
    ! Basic stencil; refinement level in the grid points:
    !/
    integer                    :: iLevel_I(2**nDim)
    
    real, dimension(nDim)      :: XyzStencil_D
   
    logical:: IsNearBoundary
    !\
    ! The sort of stencil as derived from EXTENDED stencil
    !/
    integer:: iCaseExtended
    !\
    !Just 2**nDim and a loop variable
    !/
    integer:: nGrid, iGrid
    !/
    !\
    ! To improve the algorithm stability against roundoff errors
    !/
    real           :: cTol2
    !------------------------  
    !\
    ! Initialize 
    !/
    iIndexes_II = 0; Weight_I    = 0; iOrder_I = 0
    nGridOut = -1
    cTol2 = cTol**(nByteReal/4)
    nGrid = 2**nDim !Number of points in a basic stencil
    IsNearBoundary = any(IsOut_I)
    if(maxval(iLevelSubgrid_I,MASK=.not.IsOut_I)==Coarse_) then
       !\
       ! Calculate nDim-linear interpolation on uniform grid
       !/
       call interpolate_uniform(nDim=nDim, &
            Dimless_D=(Xyz_D - XyzGrid_DII(:,1,1))*DxyzInv_D, &
            Weight_I=Weight_I, IsOut_I=IsOut_I)
       !\
       ! Form index array and sort out zero weights
       !/
       nGridOut = 0 ; cTol2 = 2*cTol2
       do iGrid = 1, nGrid
          if(Weight_I(iGrid) < cTol2)CYCLE
          nGridOut = nGridOut + 1
          iIndexes_II(0,       nGridOut) = iProc_I( iGrid)
          iIndexes_II(nIndexes,nGridOut) = iBlock_I(iGrid)
          iIndexes_II(1:nDim,  nGridOut) = iCellIndexes_DII(:,1,iGrid)
          Weight_I(nGridOut) = Weight_I(iGrid)
       end do
       if(present(IsSecondOrder))IsSecondOrder = .not.IsNearBoundary
       RETURN
    end if
    !\
    !Start non-umifirm grid
    !/
    if(IsNearBoundary)call prolong_beyond_boundary
    if(DoInit)call init_sort_stencil
    select case(nDim)
    case(2)
       call check_transition2(&
            Xyz_D, DxyzInv_D, XyzGrid_DII(:,0,:), iLevelSubgrid_I, &
            XyzStencil_D, iCaseExtended)
       call generate_basic_stencil(XyzStencil_D) 
       call interpolate_amr2(&
            Xyz_D , XyzGrid_DI, iLevel_I, IsOut_I, iCaseExtended, &
            nGridOut, Weight_I, iOrder_I)
    case(3)
       !\
       ! For edges and corners we need to find, if point Xyz is
       ! really within the corner, or it falls into some transition 
       ! region (from edge to corner, from resolution interface to
       ! edge. iCaseExtended for these cases takes the values 
       ! Transition2Edge_, Transition2Corner_, TransitionJunction_
       !/ 
       call check_transition3(&
            Xyz_D, DxyzInv_D, XyzGrid_DII(:,0,:), iLevelSubgrid_I, &
            XyzStencil_D, iCaseExtended, IsCorner)
       if(.not.IsCorner)then
          call generate_basic_stencil(XyzStencil_D)
          call interpolate_amr3(&
               Xyz_D , XyzGrid_DI, iLevel_I, IsOut_I, iCaseExtended,&
               nGridOut, Weight_I, iOrder_I, IsCorner)
       end if
       if(IsCorner)then
          !\
          ! This is not necessarily esleif,
          ! in case of sophisticated transition from edge to corner 
          ! amr3 routins finally decides if the point should be 
          ! interpolated with the corner stencil (if the interpolation 
          ! algorithm for tramsition region fails)
          !/
          call generate_corner_stencil
          call resolution_corner(Xyz_D , XyzGrid_DI, iLevel_I,&
               nGridOut, Weight_I, iOrder_I)
       end if
    end select
    !\
    ! Eliminate repeating grid points and points with zero weight
    ! which may be behind the domain boundary
    !/
    call sort_out
    iIndexes_II(:, 1:nGridOut) = iIndexes_II(:,iOrder_I(1:nGridOut))
    if(present(IsSecondOrder))IsSecondOrder = .not.IsNearBoundary
  contains
    !====================
    subroutine prolong_beyond_boundary
      !\
      ! Handle points behind the boundary
      !/
      !\
      ! Misc
      !/
      integer:: iGrid, iSubgrid, iDir, iLoc
      !---------------- 
      select case(count(IsOut_I))
      case(4)
         !\
         ! nDim = 3, one of the faces is fully out of the domain
         !/
         !\
         !Find point in the domain
         !/
         iLoc = maxloc(iLevelSubGrid_I,MASK=.not.IsOut_I,DIM=1)
         do iDir = 1, nDim
            if(all(IsOut_I(iOppositeFace_IDI(:,iDir,iLoc))))then
                  !\
                  ! Prolong grid from the physical face to the ghost one
                  ! accounting for the difference in resolution
                  !/
               do iGrid = 1,4
                  call prolong(iFace_IDI(iGrid,iDir,iLoc), &
                       iOppositeFace_IDI(iGrid,iDir,iLoc))
               end do
               RETURN
            end if
         end do
      case(2)
         ! Analogous to the previous case, but nDim = 2 
         !\
         !Find point in the domain
         !/
         iLoc = maxloc(iLevelSubGrid_I,MASK=.not.IsOut_I, DIM=1)
         do iDir = 1, 2
            if(all(&
                 IsOut_I(iOppositeSide_IDI(:,iDir,iLoc))))then
               !\
               ! Prolong grid from the physical face to the ghost 
               ! accounting for the difference in resolution
               !/
               do iGrid = 1,2
                  call prolong(iSide_IDI(iGrid,iDir,iLoc), &
                       iOppositeSide_IDI(iGrid,iDir,iLoc))
               end do
               RETURN
            end if
         end do
      case(6)
         !\
         ! Three-dimensional case, only two grid vertexes are
         ! inside the domain. From each of the two physical
         ! points the grid should be prolonged to three
         ! ghost points along the plane
         !/
         iLoc = maxloc(iLevelSubGrid_I,MASK=.not.IsOut_I, DIM=1)
         do iDir  = 1, nDim
            if(.not.IsOut_I(iEdge_ID(iLoc,iDir)))then
               do iGrid = 2,4
                  call prolong(iLoc,&
                       iFace_IDI(iGrid,iDir,iLoc))
                  call prolong(iEdge_ID(iLoc,iDir),&
                       iOppositeFace_IDI(iGrid,iDir,iLoc))
               end do
               RETURN
            end if
         end do
      end select
    end subroutine prolong_beyond_boundary
    !======================
    !\
    ! Used to prolong grid behind the boundary
    !/
    subroutine prolong(iGridPhys, iGridGhost)          
      integer, intent(in) :: iGridPhys, iGridGhost
      !Loop variable
      integer :: iSubGrid
      !--------------------
      iLevelSubgrid_I(iGridGhost) = iLevelSubgrid_I(iGridPhys) 
      !\
      !nGrid for refined subgrid 
      !/
      do iSubGrid = 1, (nGrid - 1)*iLevelSubgrid_I(iGridPhys) + 1 
         XyzGrid_DII(:,iSubGrid,iGridGhost) = &
              XyzGrid_DII(:,iSubGrid,iGridPhys) +&
              XyzGrid_DII(:,0,iGridGhost) - XyzGrid_DII(:,0,iGridPhys )
      end do
    end subroutine prolong
    !==================
    subroutine generate_corner_stencil
      integer:: iGrid
      iLevel_I = iLevelSubgrid_I
      ! IsOut_I = IsOut_I
      do iGrid = 1, nGrid
         iIndexes_II(0,        iGrid) =  iProc_I(iGrid)
         iIndexes_II(nIndexes, iGrid) = iBlock_I(iGrid)
         if(iLevel_I(iGrid)==Coarse_)then
            XyzGrid_DI(:,iGrid) = XyzGrid_DII(:,1,iGrid)
            iIndexes_II(1:nDim,iGrid) = iCellIndexes_DII(:,1,iGrid)
         else
            XyzGrid_DI(:,iGrid) = XyzGrid_DII(:,9 - iGrid,iGrid)
            iIndexes_II(1:nDim,iGrid) = &
                 iCellIndexes_DII(:,9 - iGrid, iGrid)
         end if
      end do
    end subroutine generate_corner_stencil
    !==================
    subroutine generate_basic_stencil(XyzStencil_D)  
      real,    intent(in) :: XyzStencil_D(nDim)

      integer:: nExtendedStencil, iGrid_I(64), iSubgrid_I(64)
      integer:: iGrid, iPoint, iSubgrid !Loop variables
      logical:: IsMask_I(64), IsOutSaved_I(nGrid) 
      logical:: IsBelowExtended_DI(nDim,64)
      real   :: Distance_I(64)
      !\
      ! Arrays of logical values Xyz_D < XyzGrid_DI(:,iGrid) for iGrid'th grid
      ! point of the basic stencil for point Xyz
      !/
      logical, parameter, dimension(3,8):: IsBelow_DI=reshape((/&
           .false., .false., .false., &
           .true. , .false., .false., &
           .false., .true. , .false., &
           .true. , .true. , .false., &
           .false., .false., .true. , &
           .true. , .false., .true. , &
           .false., .true. , .true. , &
           .true. , .true. , .true. /)&
           , (/3,8/) )
      !-------------------------------------

      IsOutSaved_I = IsOut_I
      iPoint = 0
      do iGrid = 1, nGrid
         do iSubgrid = 1, 1 + (nGrid - 1)*iLevelSubgrid_I(iGrid)
            iPoint = iPoint + 1
            iGrid_I(iPoint) = iGrid
            iSubGrid_I(iPoint) = iSubGrid
            IsBelowExtended_DI(:,iPoint) = &
                 XyzStencil_D < XyzGrid_DII(:,iSubGrid,iGrid) 
            Distance_I(iPoint) = sum( &
                 ((XyzStencil_D - XyzGrid_DII(:,iSubGrid,iGrid))*dXyzInv_D)**2)
         end do
      end do
      nExtendedStencil = iPoint
      do iGrid = 1, nGrid
         do iPoint = 1, nExtendedStencil
            !\
            ! Mark all the candidates for the role of iGrid point of  
            ! the basic stencil
            !/
            IsMask_I(iPoint) = all(&
                 IsBelowExtended_DI(:,iPoint).eqv.IsBelow_DI(1:nDim,iGrid))
         end do
         !\
         !Among the candidates chose the closest one
         !/
         iPoint = minloc(Distance_I(1:nExtendedStencil),&
              MASK = IsMask_I(1:nExtendedStencil), DIM=1)
         XyzGrid_DI(:,iGrid) = &
              XyzGrid_DII(:,iSubGrid_I(iPoint),iGrid_I(iPoint)) 
         iIndexes_II(0,iGrid) = &
              iProc_I(iGrid_I(iPoint))
         iIndexes_II(nIndexes,iGrid) = &
              iBlock_I(iGrid_I(iPoint))
         iIndexes_II(1:nDim,iGrid)   = &
              iCellIndexes_DII(:,iSubGrid_I(iPoint),iGrid_I(iPoint))
         iLevel_I(iGrid) = iLevelSubgrid_I(iGrid_I(iPoint))
         IsOut_I(iGrid) = IsOutSaved_I(iGrid_I(iPoint))
      end do
    end subroutine generate_basic_stencil
    !==================
    subroutine sort_out
      !\
      ! Sorts out zero weights and repeating points
      !/
      integer:: nStored !To store starting nGridOut
      integer:: iLoc    !To find location of repeating index
      integer:: iGrid   !Loop variable
      !\
      ! Form index array and sort out zero weights
      !/
      nStored =  nGridOut 
      nGridOut = 0 
      cTol2 = 2*cTol2
      iLoc = 0
      ALL:do iGrid = 1, nStored
         if(Weight_I(iGrid) < cTol2)CYCLE
         do iLoc = 1, nGridOut
            if(iOrder_I(iLoc)==iOrder_I(iGrid))then
               Weight_I(iLoc) = Weight_I(iLoc) + Weight_I(iGrid)
               CYCLE ALL
            end if
         end do
         nGridOut = nGridOut + 1
         iOrder_I(nGridOut) = iOrder_I(iGrid)
         Weight_I(nGridOut) = Weight_I(iGrid)
      end do ALL
    end subroutine sort_out
  end subroutine interpolate_extended_stencil
  !=====================================================================
  subroutine check_transition2( Xyz_D, dXyzInv_D, XyzGrid_DI, iLevel_I,&
       XyzStencil_D, iCase)
    integer, parameter:: nDim = 2, nGrid = 4
    !\
    ! Point where to interpolate; inverse of (/Dx,Dy,Dz/)
    !/
    real,    intent(in) :: Xyz_D(nDim), dXyzInv_D(nDim)
    !\
    ! The rectangular grid combining centers of the refined subgrids and
    ! coarse vertexes
    !/
    real, intent(in):: XyzGrid_DI(nDim,nGrid)
    !\
    ! The refinement level pattern of the extended stencil
    !/
    integer, intent(in)::iLevel_I(nGrid)
    !\
    ! Point about which to construct basic stencil.
    !/
    real, intent(out) :: XyzStencil_D(nDim)
    !\
    ! For three-dimensional corner within this routine
    ! it is more convenient to figure out if Xyz point is
    ! inside the domain of transition from an edge to the central
    ! corner part or it belongs to a junction of such two domains
    !/
    integer, intent(out):: iCase
    !\
    !Dimensionless displacement from the first grid point to Xyz_D
    !/
    real :: Dimless_D(nDim)
    !\
    !Three components of Discr_D equal to -1,0 or 1 each. The first component 
    !equals -1, if Xyz point is close to refined face x=0, +1 if it is close
    !to refined face x=1, 0 otherwise
    !/
    integer:: iDiscr_D(nDim), iDir, iDim, iGrid, jGrid
    !\ 
    !Dir = 0 - corner, positive iDirs - face, negative iDir - edges, which
    !occur at the intersection of the refined faces, which are close to
    !Xyz
    !/
    integer, parameter, dimension(2,-1:1,-1:1):: &
         iGridDir_III = reshape((/&
         1,0,     1,2,  2,0,  1,1,  0,0, 2,1,  3,0,   3,2, 4,0/),(/2,3,3/))
    !x=0,y=0!y=0 !x=1,y=0!x=0 !    !x=1 !x=0,y=1!y=1 !x=1,y=1!
    real:: XyMin, z
    !\
    ! Average coordinates for the center of the resolution corner,
    ! for centers of faces or edges and distance from them to Xyz
    !/
    real   :: XyzAvr_D(nDim), XyzAvr_DD(nDim, nDim), Distance_D(nDim) 
    !----------------------
    XyzStencil_D = Xyz_D
    iCase = i_case(iLevel_I)
    !\                                    F F
    !We do not care about trapezoids like  C   C
    !/
    if(iSortStencil3_II(Case_,iCase) <= Face_)RETURN
    Dimless_D = (Xyz_D - XyzGrid_DI(:,1))*DxyzInv_D 
    iDiscr_D = 0
    do iDim = 1, nDim
       if(Dimless_D(iDim) <  0.250.and.any(iLevel_I(&
            iSide_IDI(:,3 - iDim,1))==Fine_))iDiscr_D(iDim) = -1
       if(Dimless_D(iDim) >=  0.750.and.any(iLevel_I(&
            iOppositeSide_IDI(:,3 - iDim,1))==Fine_))iDiscr_D(iDim) = 1
    end do
    iGrid = iGridDir_III(1, iDiscr_D(1), iDiscr_D(2))
    XyzAvr_D = 0.50*(XyzGrid_DI(:,1) + XyzGrid_DI(:,4))
    if(iGrid==0)then 
       !\
       ! The point is in the close proximity of the resolution edge
       !/
       XyzStencil_D = XyzAvr_D 
       RETURN
    end if
    !\
    !Xyz point is in the transition domain
    !/   
    iDir = iGridDir_III(2, iDiscr_D(1), iDiscr_D(2))
    if(iDir== 0)then
       !\
       !The point is close to two faces
       !/
       do iDim = 1,nDim
          XyzAvr_DD(:,iDim) = 0.50*(&
               XyzGrid_DI(:,iSide_IDI(1,3 - iDim,iGrid)) +&
               XyzGrid_DI(:,iSide_IDI(2,3 - iDim,iGrid)) )
          Distance_D(iDim) = &
               sum(((Xyz_D - XyzAvr_DD(:,iDim))*DxyzInv_D)**2) 
       end do
       iDir = minloc(Distance_D(1:2),DIM=1)
       !\
       !else Xyz point is close to a single refined face 
       !/
    end if
    iCase = Transition2Edge_
    !\
    ! Displace the stencil center toward the face center
    !/
    XyzStencil_D(3 - iDir) = XyzAvr_D(3 - iDir)
  end subroutine check_transition2
  !=========================================================
  subroutine  interpolate_amr2(&
       Xyz_D, XyzGridIn_DI, iLevel_I, IsOut_I, iCaseExtended,&
       nGridOut, Weight_I, iOrder_I)
    use ModCubeGeometry, ONLY: iSortStencil2_II
    integer,parameter :: nGrid = 4, nDim = 2

    character(LEN=*),parameter:: NameSub='interpolate_amr2'
    !\
    !Input parameters
    !/
    !\
    !The location at which to interpolate the data
    !/
    real, intent(in):: Xyz_D(nDim) 

    !Grid point coordinates !2 coordinate, 4 points
    real  ,   intent(in) :: XyzGridIn_DI(nDim,nGrid) 
    !The same, but reordered
    real                 :: XyzGrid_DI(nDim,nGrid) 

    !The refinement level at each grid point. By one higher level of refinement
    !assumes the cell size reduced by a factor of 0.5
    integer,  intent(in) :: iLevel_I(nGrid)

    !\
    ! Logical which marks "ghost" points, not belonging to the computational 
    ! domain
    !/
    logical, intent(in) :: IsOut_I(nGrid)
    integer, intent(in) :: iCaseExtended
    !\
    !Output parameters
    !/
    !The number of grid points to be ultimately included into 
    !the interpolation stencil. If nGridOut < nGridIn, only the first 
    !nGridOut lines are meaningful in the output
    integer, intent(out) :: nGridOut

    !The weight coefficients array.
    real   , intent(out) :: Weight_I(nGrid)

    !Order(numbers) of grid points used for the interpolation
    integer, intent(out) :: iOrder_I(nGrid)
    !\
    !Parameters fully characterizing a pattern of refinement for a stencil
    !Calculated using iSortStencil_III array
    !/
    integer :: iCase, iGrid, iDir  
    !\
    ! 3 - iDir: 
    ! for coordinate x_ the perpendicular cordinate is y_ and vice versa 
    !/
    integer ::  iDirPerp              
    !\
    ! Misc
    !/ 
    real:: Aux               
    !-------------
    !\
    ! For input iLevel_I find a type of a stencil
    !/
    iCase = i_case(iLevel_I)
    iGrid = iSortStencil2_II(Grid_,iCase)
    iDir  = iSortStencil2_II(Dir_, iCase)
    iCase = iSortStencil2_II(Case_,iCase)
    !\
    ! Check if this is a transion to the edge, not an edge.
    !/
    if(iCaseExtended == Transition2Edge_.and.iCase==OneCoarse_)&
         iCase = Transition2Edge_
    iOrder_I(1:2) = iSide_IDI(:,iDir,iGrid)
    iOrder_I(3:4) = iOppositeSide_IDI(:,iDir,iGrid)
    XyzGrid_DI = XyzGridIn_DI(:,iOrder_I)
    select case(iCase)
    case(Face_)
       !\                                        C   C
       ! Resolution line going along iDir-axis    F F  or  F F   ---> iDir
       !/                                                C  C           
       !\
       ! Interpolate along edge X1X2 (iGrid = 1) and X3X4 (iGrid = 3)
       !/
       do iGrid = 1,3,2
          !\
          ! Check if the grid points are out of the grid boundary
          !/
          if(IsOut_I(iOrder_I(iGrid)))then
             Weight_I(iGrid) = 0; Weight_I(iGrid + 1) = 1 
          elseif(IsOut_I(iOrder_I(iGrid + 1)))then
             Weight_I(iGrid + 1) = 0; Weight_I(iGrid) = 1
          else
             Weight_I(iGrid + 1) = (Xyz_D(iDir) -  XyzGrid_DI(iDir,iGrid) ) /&
                  ( XyzGrid_DI(iDir,iGrid + 1) -  XyzGrid_DI(iDir,iGrid) )
             Weight_I(iGrid) = 1 - Weight_I(iGrid + 1) 
          end if
       end do
       !\
       ! Apply weight for interpolation along another axis
       !/
       iDirPerp = 3 - iDir
       Aux = (Xyz_D(iDirPerp) - XyzGrid_DI(iDirPerp,1 )  )/&
            (XyzGrid_DI(iDirPerp,4) - XyzGrid_DI(iDirPerp,1 )) 
       Weight_I(1:2) =  Weight_I(1:2)*(1 - Aux)
       Weight_I(3:4) =  Weight_I(3:4)*Aux
       nGridOut = 4
    case(OneCoarse_)
       !     F-F
       !    / \|
       !   C---F
       call triangles(iTriangle1_I=(/2, 3, 1/), iTriangle2_I=(/2, 3, 4/))
    case(OneFine_,Rhombus_,Transition2Edge_)
       !   C---C
       !    \ /|
       !     F-C
       !   C-F
       !   |/\
       !   F--C
       !       |   |  
       !     3F--4F|   
       !------------    2F here
       !      \ /  !
       !       1C  | 
       call triangles(iTriangle1_I=(/1, 4, 2/), iTriangle2_I=(/1, 4, 3/))
    end select
  contains
    !=========================
    real function cross_product(a_D, b_D)
      real,dimension(nDim),intent(in) :: a_D, b_D
      !----------
      cross_product = a_D(x_)* b_D(y_) - a_D(y_)*b_D(x_)
    end function cross_product
    !=======
    subroutine triangles(iTriangle1_I,iTriangle2_I)
      integer,          intent(in) :: iTriangle1_I(3)
      integer,optional, intent(in) :: iTriangle2_I(3)
      real, dimension(nDim) :: X1_D, X2_D, X3_D
      !-------
      X1_D=XyzGrid_DI(:,iTriangle1_I(1))
      X2_D=XyzGrid_DI(:,iTriangle1_I(2))
      X3_D=XyzGrid_DI(:,iTriangle1_I(3))
      call triangle(X1_D, X2_D, X3_D)
      if(all(Weight_I(1:3)>=0.0))then
         iOrder_I = (/iOrder_I(iTriangle1_I(1)),iOrder_I(iTriangle1_I(2)),&
              iOrder_I(iTriangle1_I(3)), iOrder_I(10 - sum(iTriangle1_I))/) 
      else
         if(.not.present(iTriangle2_I))RETURN
         X1_D=XyzGrid_DI(:,iTriangle2_I(1))
         X2_D=XyzGrid_DI(:,iTriangle2_I(2))
         X3_D=XyzGrid_DI(:,iTriangle2_I(3))
         call triangle(X1_D, X2_D, X3_D)
         iOrder_I = (/iOrder_I(iTriangle2_I(1)),iOrder_I(iTriangle2_I(2)),&
              iOrder_I(iTriangle2_I(3)), iOrder_I(10 - sum(iTriangle2_I))/) 
      end if
    end subroutine triangles
    !============
    subroutine triangle(X1_D, X2_D, X3_D)
      real, dimension(nDim), intent(in):: X1_D, X2_D, X3_D
      !-------------
      Weight_I = 0; nGridOut = -1
      Weight_I(3) = cross_product(Xyz_D - X1_D,X2_D - X1_D)/&
           cross_product(X3_D - X1_D,X2_D - X1_D)
      if(Weight_I(3) >= 0.0)then
         nGridOut = 3
         Weight_I(2) = cross_product(X3_D-X1_D,Xyz_D - X1_D)/&
              cross_product(X3_D - X1_D,X2_D - X1_D)
         Weight_I(1) = 1 - Weight_I(2) - Weight_I(3)
      end if
    end subroutine triangle
  end subroutine interpolate_amr2
  !==============================
  subroutine check_transition3( Xyz_D, dXyzInv_D, XyzGrid_DI, iLevel_I,&
       XyzStencil_D, iCase, IsCorner)
    integer, parameter:: nDim = 3, nGrid = 8
    !\
    ! Point where to interpolate; inverse of (/Dx,Dy,Dz/)
    !/
    real,    intent(in) :: Xyz_D(nDim), dXyzInv_D(nDim)
    !\
    ! The rectangular grid combining centers of the refined subgrids and
    ! coarse vertexes
    !/
    real, intent(in):: XyzGrid_DI(nDim,nGrid)
    !\
    ! The refinement level pattern of the extended stencil
    !/
    integer, intent(in)::iLevel_I(nGrid)
    !\
    ! Point about which to construct basic stencil.
    !/
    real, intent(out) :: XyzStencil_D(nDim)
    !\
    ! For three-dimensional corner within this routine
    ! it is more convenient to figure out if Xyz point is
    ! inside the domain of transition from an edge to the central
    ! corner part or it belongs to a junction of such two domains
    !/
    integer, intent(out):: iCase
    !\
    ! Is true if the point is inside the resolution corner
    !/
    logical, intent(out):: IsCorner
    !\
    !Dimensionless displacement from the first grid point to Xyz_D
    !/
    real :: Dimless_D(nDim)
    !\
    !Three components of Discr_D equal to -1,0 or 1 each. The first component 
    !equals -1, if Xyz point is close to refined face x=0, +1 if it is close
    !to refined face x=1, 0 otherwise
    !/
    integer:: iDiscr_D(nDim), iDir, iDim, iGrid, jGrid
    !\ 
    !Dir = 0 - corner, positive iDirs - face, negative iDir - edges, which
    !occur at the intersection of the refined faces, which are close to
    !Xyz
    !/
    integer, parameter, dimension(2,-1:1,-1:1,-1:1):: &
         iGridDir_IIII = reshape((/&
         1, 0,    1,-1,  2, 0, 1,-2,1,3,2,-2,  3, 0, 3,-1, 4, 0,   &!z=0!
         1,-3,    1, 2,  2,-3, 1, 1,0,0,2, 1,  3,-3, 3, 2, 4,-3,   &!   !
         5, 0,    5,-1,  6, 0, 5,-2,5,3,6,-2,  7, 0, 7,-1, 8, 0 /),&!z=1!
         (/2,3,3,3/))
    !x=0,y=0!y=0 !x=1,y=0!x=0 !    !x=1 !x=0,y=1!y=1 !x=1,y=1!
    !\
    !Transitions junction's signatures
    !/
    real:: XyMin, z
    !\
    ! Average coordinates for the center of the resolution corner,
    ! for centers of faces or edges and distance from them to Xyz
    !/
    real   :: XyzAvr_D(nDim), XyzAvr_DD(nDim, nDim), Distance_D(nDim) 
    !\
    ! To treat Edges
    !/
    logical:: IsEdge
    integer:: iDirEdge
    !----------------------
    XyzStencil_D = Xyz_D
    IsCorner = .false.
    iCase = i_case(iLevel_I)
    !\
    !There is also no need to look for transitions in faces
    !/
    if(iSortStencil3_II(Case_,iCase) <= Face_)RETURN
    IsEdge = iSortStencil3_II(Case_,iCase) == Edge_
    if(IsEdge) iDirEdge = iSortStencil3_II(Dir_,iCase)
    Dimless_D = (Xyz_D - XyzGrid_DI(:,1))*DxyzInv_D 
    iDiscr_D = 0
    do iDim = 1, nDim
       if(Dimless_D(iDim) <  0.250.and.any(iLevel_I(&
            iFace_IDI(:,iDim,1))==Fine_))iDiscr_D(iDim) = -1
       if(Dimless_D(iDim) >=  0.750.and.any(iLevel_I(&
            iOppositeFace_IDI(:,iDim,1))==Fine_))iDiscr_D(iDim) = 1
    end do
    if(IsEdge) iDiscr_D(iDirEdge) = 0
    iGrid = iGridDir_IIII(1, iDiscr_D(1), iDiscr_D(2), iDiscr_D(nDim))
    if(iGrid==0)then
       IsCorner = .not.IsEdge
       RETURN
    end if
    !\
    !Xyz point is in the transition or transition junction domain
    !/
    XyzAvr_D = 0.50*(XyzGrid_DI(:,1) + XyzGrid_DI(:,8))
    iDir = iGridDir_IIII(2, iDiscr_D(1), iDiscr_D(2), iDiscr_D(nDim))
    if(iDir > 0)then
       !Xyz point os close to a single refined face - transition!
       iCase = Transition2Corner_
       iSortStencil3_II(Grid_,Transition2Corner_) = iGrid
       iSortStencil3_II(Dir_ ,Transition2Corner_) = iDir
       !\
       ! Displace the stencil center toward the face center
       !/
       XyzStencil_D(1 + mod(iDir,3)) = XyzAvr_D(1 + mod(iDir,3))
       XyzStencil_D(1 + mod(iDir + 1,3)) = &
            XyzAvr_D(1 + mod(iDir + 1,3))
    elseif(iDir < 0)then
       !\
       !The point is close to two faces intersescting at the
       !edge of direction -iDir
       !/
       iDir = -iDir
       jGrid = iEdge_ID(iGrid,iDir)
       !\
       ! Edge consists of the coarse points or finer subgrid 
       ! iGrid,jGrid
       !/
       if(iLevel_I(iGrid)==Coarse_.and.iLevel_I(jGrid)==Fine_)then
          !\
          !Jucntion of two transition regions which goes along
          !the edge of direction iDir and closed with fine
          !subgrid at jGrid end
          iCase = TransitionJunction_
          iSortStencil3_II(Grid_,TransitionJunction_) = iGrid
          iSortStencil3_II(Dir_ ,TransitionJunction_) = iDir    
          XyzStencil_D(iDir) = XyzAvr_D(iDir)
       elseif(iLevel_I(jGrid)==Coarse_.and.iLevel_I(iGrid)==Fine_)then
          !\
          !Jucntion of two transition regions which goes along
          !the edge of direction iDir and closed with fine
          !subgrid at iGrid end
          iCase = TransitionJunction_
          iSortStencil3_II(Grid_,TransitionJunction_) = jGrid
          iSortStencil3_II(Dir_ ,TransitionJunction_) = iDir   
          XyzStencil_D(iDir) = XyzAvr_D(iDir)
       else
          !\
          ! We need to judge to which of the two transition regions
          ! intersecting along the coarse edge, that is without 
          ! forming a specific junction region as closed with the  
          ! fine subface, point Xyz belongs 
          !/
          do iDim = 0,1
             XyzAvr_DD(:,1+iDim) = 0.50*(&
                  XyzGrid_DI(:,&
                  iFace_IDI(1,1 + mod(iDir + iDim,3),iGrid)) +&
                  XyzGrid_DI(:,&
                  iFace_IDI(4,1 + mod(iDir + iDim,3),iGrid)) )
             Distance_D(1 + iDim) = &
                  sum(((Xyz_D - XyzAvr_DD(:,1+iDim))*DxyzInv_D)**2) 
          end do
          iDim = minloc(Distance_D(1:2),DIM=1)
          iDir = 1 + mod(iDir + iDim - 1,3)
          iCase = Transition2Corner_
          iSortStencil3_II(Grid_,Transition2Corner_) = iGrid
          iSortStencil3_II(Dir_ ,Transition2Corner_) = iDir
          !\
          ! Displace the stencil center toward the face center
          !/
          XyzStencil_D(1 + mod(iDir,3)) = XyzAvr_D(1 + mod(iDir,3))
          XyzStencil_D(1 + mod(iDir + 1,3)) = &
               XyzAvr_D(1 + mod(iDir + 1,3))
       end if
    else
       !\
       ! Xyz point is close to three refined planes
       !/
       select case(count(iLevel_I(iEdge_ID(iGrid,:))==Fine_))
       case(0)
          !\
          ! We need to judge to which of the three transition regions
          ! intersecting at the corner without forming a specific
          ! junction region, point Xyz belongs 
          !/
          do iDim = 1,nDim
             XyzAvr_DD(:,iDim) = 0.50*(&
                  XyzGrid_DI(:,iFace_IDI(1,iDim,iGrid)) +&
                  XyzGrid_DI(:,iFace_IDI(4,iDim,iGrid)) )
             Distance_D(iDim) = &
                  sum(((Xyz_D - XyzAvr_DD(:,iDim))*DxyzInv_D)**2) 
          end do
          iDir = minloc(Distance_D,DIM=1)
          iCase = Transition2Corner_
          iSortStencil3_II(Grid_,Transition2Corner_) = iGrid
          iSortStencil3_II(Dir_ ,Transition2Corner_) = iDir
          !\
          ! Displace the stencil center toward the face center
          !/
          XyzStencil_D(1 + mod(iDir,3)) = XyzAvr_D(1 + mod(iDir,3))
          XyzStencil_D(1 + mod(iDir + 1,3)) = &
               XyzAvr_D(1 + mod(iDir + 1,3))
       case(2,3)
          !\
          ! We need to judge to which of the two or three transition 
          ! junctions point Xyz belongs 
          !/
          do iDim = 1,nDim
             !\
             !If there is no junction region along this direction,
             !that is iLevel_I(iEdge_ID(iGrid,iDim))=0, move the
             !center of this edge far apart, otherwise take the 
             !center of the edge  
             !/
             XyzAvr_DD(:,iDim) = (&
                  XyzGrid_DI(:,iGrid)*iLevel_I(iEdge_ID(iGrid,iDim))&
                  + XyzGrid_DI(:,iEdge_ID(iGrid,iDim)) )/&
                  (iLevel_I(iEdge_ID(iGrid,iDim)) + 1)
             Distance_D(iDim) = &
                  sum(((Xyz_D - XyzAvr_DD(:,iDim))*DxyzInv_D)**2) 
          end do
          iDir = minloc(Distance_D,DIM=1)
          iCase = TransitionJunction_
          iSortStencil3_II(Grid_,TransitionJunction_) = iGrid
          iSortStencil3_II(Dir_ ,TransitionJunction_) = iDir
          XyzStencil_D(iDir) = XyzAvr_D(iDir)
       case(1)
          !\
          ! We need to find the only transition junction and then  
          ! judge if point Xyz belongs to this transition junction  
          ! or to the transition region near the lower face of the
          ! same direction. 
          !/
          iDir = minloc(iLevel_I(iEdge_ID(iGrid,:)), MASK=&
               iLevel_I(iEdge_ID(iGrid,:))==Fine_, DIM=1)
          !
          !        ^     7F----8F
          !   iDir |     /    /
          !            5F----6F
          !                  /3C
          !                 /     4F
          !                /________
          !               1C       2C
          !\
          ! As a separator between the transition junction 1C5F6F7F8F
          ! and a transition region near face 1C2C3C4F we use a
          !  surface, z = min(x,y)*Alpha, which passes through
          !  point 4F at Alpha = 1/3 and through point 8F at Alpha=3. 
          !/
          z = (Xyz_D(iDir) - XyzGrid_DI(iDir,iGrid))/  &
               (XyzGrid_DI(iDir,iEdge_ID(iGrid,iDir)) -&
               XyzGrid_DI(iDir,iGrid))
          XyMin = min( (Xyz_D(1 + mod(iDir,3))        -&
               XyzGrid_DI(1 + mod(iDir,3),iGrid))/     &
               (XyzGrid_DI(1 + mod(iDir,3),            &
               iEdge_ID(iGrid,1 + mod(iDir,3)))       -&
               XyzGrid_DI(1 + mod(iDir,3),iGrid)),     &
               (Xyz_D(1 + mod(iDir + 1,3))            -&
               XyzGrid_DI(1 + mod(iDir + 1,3),iGrid))/ &
               (XyzGrid_DI(1 + mod(iDir + 1,3),        &
               iEdge_ID(iGrid,1 + mod(iDir + 1,3)))   -&
               XyzGrid_DI(1 + mod(iDir + 1,3),iGrid)) )
          if(z  >= 3*XyMin)then
             iCase = TransitionJunction_
             iSortStencil3_II(Grid_,TransitionJunction_) = iGrid
             iSortStencil3_II(Dir_ ,TransitionJunction_) = iDir
             XyzStencil_D(iDir) = XyzAvr_D(iDir)
          elseif(z <= XyMin/3)then
             iCase = Transition2Corner_
             iSortStencil3_II(Grid_,Transition2Corner_) = iGrid
             iSortStencil3_II(Dir_ ,Transition2Corner_) = iDir
             !\
             ! Displace the stencil center toward the face center
             !/
             XyzStencil_D(1 + mod(iDir,3)) = XyzAvr_D(1 + mod(iDir,3))
             XyzStencil_D(1 + mod(iDir + 1,3)) = &
                  XyzAvr_D(1 + mod(iDir + 1,3))
          else
             IsCorner = .true.
             XyzStencil_D = XyzAvr_D 
          end if
       end select
    end if
    if(IsEdge)then
       iCase = Transition2Edge_
       iSortStencil3_II(Dir_, iCase) = iDirEdge
       XyzStencil_D(iDirEdge) = Xyz_D(iDirEdge)
    end if
  end subroutine check_transition3
  !============================
  subroutine  interpolate_amr3(&
       XyzIn_D, XyzGridIn_DI, iLevel_I, IsOut_I, iCaseExtended,&
       nGridOut, Weight_I, iOrder_I, IsCorner)
    integer,parameter :: nGrid = 8, nDim = 3
    character(LEN=*),parameter:: NameSub='interpolate_amr3'
    !\
    !Input parameters
    !/
    !\
    !The location at which to interpolate the data
    !/
    real   , intent(in) :: XyzIn_D(nDim)
    !\
    ! Same, but reordered
    !/
    real                :: Xyz_D(nDim)
    !\
    !Grid point coordinates !3 coordinate, 8 points
    !/
    real  ,   intent(in) :: XyzGridIn_DI(nDim,nGrid)
    real :: XyzGrid_DI(nDim,nGrid) !The same with changed order
    !\
    !The refinement level at each grid point. By one higher level
    ! of refinement assumes the cell size reduced by a factor of 0.5
    !/
    integer,  intent(in) :: iLevel_I(nGrid)
    !\
    ! Logical which marks "ghost" points, not belonging to the 
    ! computational domain
    !/
    logical, intent(in) :: IsOut_I(nGrid)
    integer, intent(in) :: iCaseExtended
    !\
    !Output parameters
    !/
    !The number of grid points to be ultimately included into 
    !the interpolation stencil If nGridOut < nGrid only the
    !first nGridOut lines in the output arrays are meaningful.
    !/
    integer, intent(out) :: nGridOut
    !\
    !The weight coefficients array.
    real   , intent(out) :: Weight_I(nGrid)
    !\
    !Order(numbers) of grid points used for the interpolation
    integer, intent(out) :: iOrder_I(nGrid)
    !\
    ! The logical is true, if (for some resolution combinations) the 
    ! basic stencil for the point at which to inpertpolate is not 
    ! applicable
    !/
    logical, intent(out):: IsCorner
    !\
    ! Misc
    !/
    real    :: zMin, zMax                
    real    :: Aux, AuxCoarse, AuxFine   ! Misc
    !\
    ! To call two-dimensional interpolation
    !/
    integer :: nGridOut2, iOrder2_I(4),iLevel2_I(4)
    logical :: IsOut2_I(4)

    !\
    ! Auxilary variables to determine vertical weight
    ! for transition
    !/
    real :: dXyzUp, dXyzDown
    integer :: nFine, nCoarse

    !\
    ! To find the sort of corner stencil
    !/
    integer:: iCase, iDir, iGrid
    !\
    ! Switches in the algorithm
    !/
    logical,parameter:: UseFaceInsideCorner =.true.,&
         UseEdgeInsideCorner = .true.
    !-------------
    Weight_I = 0
    !\
    ! We store in advance the 'basic' grid point
    ! and orientation of all possible stencil configurations
    ! and now we extracted this information
    !/ 
    iCase = i_case(iLevel_I)
    !\
    ! Check junction and transitions
    !/
    if(.not.UseFaceInsideCorner)then
       !\
       !The most general and easy-to-understand algorithm
       !Face algorithm is only applied to transitions to edges
       !/
       if(.not.(iSortStencil3_II(Case_,iCase)==Face_.and.&
            iCaseExtended==Transition2Edge_))&
            iCase = iCaseExtended
    elseif(UseEdgeInsideCorner)then
       !\
       ! Both edge and faces basic stencils are used if
       ! occur in the transition or junction regions
       !/
       if(iSortStencil3_II(Case_,iCase) > Edge_&
            .and.iCaseExtended >= Transition2Corner_)&
            iCase = iCaseExtended
       if(iSortStencil3_II(Case_,iCase) == Edge_&
            .and.iCaseExtended == Transition2Edge_)&
            iCase = iCaseExtended 
    else
       !\
       ! face interpolation only is applied whenever possible
       !/
       if(.not.iSortStencil3_II(Case_,iCase)==Face_)&
            iCase = iCaseExtended
    end if
    iGrid = iSortStencil3_II(Grid_,iCase)
    iDir  = iSortStencil3_II(Dir_,iCase)
    iCase = iSortStencil3_II(Case_,iCase)
    IsCorner = .false. 
    iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
    iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
    !\
    ! Reorder the coordinate arrays
    !/
    Xyz_D = XyzIn_D((/1 + mod(iDir,3), 1 + mod(iDir + 1,3), iDir/))
    XyzGrid_DI = XyzGridIn_DI((/1 + mod(iDir,3), 1 + mod(iDir + 1,3), iDir/),&
         iOrder_I)
    !\
    ! Reassign iCase for edges appearing from transitions, which require 
    ! stencil fix
    !/
    if(UseEdgeInsideCorner.and.iCase==Edge_.and.&
         (iCaseExtended == TransitionJunction_.or.&
         (iCaseExtended == Transition2Corner_.and.&
         iDir/=iSortStencil3_II(Dir_,Transition2Corner_))))&
         iCase = Transition2Edge_
    select case(iCase)
    case(Face_)
       !\
       ! Faces going along xy - plane 
       !/
       !\
       ! Interpolate along lower face (iGrid = 1) and upper face (iGrid = 5)
       !/
       do iGrid = 1,5,4
          call interpolate_uniform(nDim=2,&
               Dimless_D = (Xyz_D(x_:y_) - XyzGrid_DI(x_:y_,iGrid)) /&
               (XyzGrid_DI(x_:y_,iGrid + 3) - XyzGrid_DI(x_:y_,iGrid)), &
               Weight_I=Weight_I(iGrid:iGrid + 3), &
               IsOut_I=IsOut_I(iOrder_I(iGrid:iGrid + 3)))
       end do
       !\
       ! Apply weight for interolation along z
       !/
       Aux = (Xyz_D(z_) - XyzGrid_DI(z_,1))/&
            (XyzGrid_DI(z_,5) - XyzGrid_DI(z_,1))
       Weight_I(1:4) = Weight_I(1:4)*(1 - Aux)
       Weight_I(5:8) = Weight_I(5:8)*Aux
       nGridOut = nGrid
    case(Edge_,Transition2Edge_)
       !\
       ! Edges going along resolution edge
       ! Opposite "faces" have the same geometry in projection
       ! term "face" is used to refer to 4 upper-, left-, 
       ! front-most etc. grid points they do not belong to plane
       !/      
       !\
       ! "Faces" have the same geometry in projection on xy-plane
       !/
       !\
       ! Interpolation weights along z_ axis are calculated separately 
       ! for coarse and fine points
       !/
       zMin   = minval(XyzGrid_DI(z_,1:4),MASK=iLevel_I(iOrder_I(1:4))/=Fine_)
       AuxCoarse  = (Xyz_D(z_) - zMin)/(maxval(&
            XyzGrid_DI(z_,5:8),MASK=iLevel_I(iOrder_I(5:8))/=Fine_) - zMin)

       zMin = minval(XyzGrid_DI(z_,1:4),MASK=iLevel_I(iOrder_I(1:4))==Fine_)
       AuxFine  = (Xyz_D(z_) - zMin)/(maxval(&
            XyzGrid_DI(z_,5:8),MASK=iLevel_I(iOrder_I(5:8))==Fine_) - zMin) 
       !\
       ! Interpolate along lower face
       !/
       iLevel2_I(1:4) = iLevel_I( iOrder_I(1:4)); IsOut2_I = .false.
       call interpolate_amr2(&
            Xyz_D(x_:y_), XyzGrid_DI(x_:y_,1:4), iLevel2_I, IsOut2_I, iCase, &
            nGridOut2, Weight_I(1:4), iOrder2_I)
       iOrder_I(1:3) = iOrder_I(iOrder2_I(1:3))
       iOrder_I(4:6) = iOrder_I(iOrder2_I(1:3)+4)
       nGridOut = 6
       !\
       ! Apply weight for interpolation along iDir
       !/
       do iGrid = 1,3
          if(IsOut_I(iOrder_I(iGrid)))then
             Weight_I(3 + iGrid) = Weight_I(iGrid);  Weight_I(iGrid) = 0
          elseif(&
               IsOut_I(iOrder_I(3 + iGrid)))then
             Weight_I(3 + iGrid) = 0
          elseif(iLevel_I(iOrder_I(iGrid))==Fine_)then
             Weight_I(3 + iGrid) = Weight_I(iGrid)*     AuxFine
             Weight_I(    iGrid) = Weight_I(iGrid)*(1-  AuxFine)
          else
                Weight_I(3 + iGrid) = Weight_I(iGrid)*     AuxCoarse
                Weight_I(    iGrid) = Weight_I(iGrid)*(1 - AuxCoarse)
          end if
       end do
    case(TransitionJunction_)
       !\
       ! Check corner transition junction
       !/
       !\
       ! Junction around direction z_
       !/
       ! Face with all Fine points is numbered 5:8
       ! Face with one to three Coarse points is numbered 1:4
       !
       !
       !  F5--F7
       !  |\ /|
       !  | C1|-----3
       !  |/|\|
       !  F6--F8
       !    |      4
       !    |
       !    2
       !
       !----------------------------------------------------------------
       !\
       ! Stencil part C1F5F6F7F8 is symmetric with respect to the plane 158
       ! Find if point Xyz is from the same side of the plane as point 2 is
       ! or as point 3 is
       !/
       if(sum((Xyz_D - XyzGrid_DI(:,1))*&
            (XyzGrid_DI(:,6) - XyzGrid_DI(:,7))) > 0)then
          iGrid = 2
       else
          iGrid = 3
       end if
       if(iLevel_I(iOrder_I(iGrid))==Coarse_)then
          call interpolate_on_parallel_rays(Dir_D=0.50*(&
               XyzGrid_DI(:,8) + XyzGrid_DI(:,5)) - XyzGrid_DI(:,1),&
               iURectangle1_I=(/5, 6, 7, 8/),&
               iDTriangle1_I=(/1, iGrid, 8/))
       else
          call interpolate_pyramids(iRectangular1_I = (/5, 6, 7, 8, 1/),&
               iRectangular2_I= (/iGrid, 4, iGrid +4, 8, 1/))
       end if
       IsCorner = nGridOut <= 3
    case(Transition2Corner_)
       !\
       ! Subroutine intepolates in transitional near a resolution corner
       ! Points facing a corner are numbered 5:8
       !/
       !-------------------------------------------------------------------
       !\
       ! Interpolate over face 1:4
       !/
       iLevel2_I = iLevel_I(iOrder_I(1:4)); IsOut2_I = IsOut_I(iOrder_I(1:4))
       call interpolate_amr2(&
            Xyz_D(x_:y_), XyzGrid_DI(x_:y_,1:4), iLevel2_I, IsOut2_I, Edge_,&
            nGridOut2, Weight_I(1:4), iOrder2_I)
       !\
       ! Rearrange points according to 2D interpolation output
       !/
       iOrder_I(1:4) = iOrder_I(iOrder2_I  )
       iOrder_I(5:8) = iOrder_I(iOrder2_I + 4)
       XyzGrid_DI(:,1:4) = XyzGrid_DI(:,iOrder2_I)
       XyzGrid_DI(:,5:8) = XyzGrid_DI(:,4 + iOrder2_I)
       !\
       ! Find number of Fine points in the output
       ! Move them to the begining of iOrder_I
       !/
       nFine   = 0
       do iGrid = 1, nGridOut2
          if(iLevel_I(iOrder_I(iGrid))==Fine_)then
             nFine   = nFine   + 1
             iOrder_I(    (/nFine,iGrid/)) = iOrder_I(    (/iGrid,nFine/))
             XyzGrid_DI(:,(/nFine,iGrid/)) = XyzGrid_DI(:,(/iGrid,nFine/)) 
             Weight_I(    (/nFine,iGrid/)) = Weight_I(    (/iGrid,nFine/))
             iOrder_I(     (/nFine + 4,iGrid + 4/)) = &
                  iOrder_I(     (/iGrid + 4,nFine + 4/))
             XyzGrid_DI(:,(/nFine + 4,iGrid + 4 /)) =  &
                  XyzGrid_DI(:,(/iGrid + 4, nFine + 4/)) 
          end if
       end do
       nCoarse = nGridOut2 - nFine
       !\
       ! Two exceptions result from our desire to use whenever possible
       ! interolation for edges and faces, even if the point is in the 
       ! transition region of the corner. To maintain continuity at the 
       ! boundaries with the domains in which these fake faces and edges
       ! are employed we need to go farther these way and may need either 
       ! to prolong a triangulated piece with THREE (not 4 as usually 
       ! takes place forresolution edges) parallel edges in the stencil..
       !/
       if(UseEdgeInsideCorner.and.&
            nCoarse==1.and.iLevel_I(iOrder_I(nFine+5))==Coarse_)then
          !\                            View along iDir:
          !     F-----F                    F        Both up from C point
          !    /|    /|                   /  \      and to the right from
          !   C-----C | ---->iDir         |    \    it there are junctions
          !    \|    \|                  /  ---F    treated as 'edges':
          !     F-----F                 C---    there is a ray from C point
          !/                            and one more ray from fine neighbor
          AuxCoarse = (XyzGrid_DI(z_,nFine + 5) - Xyz_D(z_))/&
               (XyzGrid_DI(z_,nFine + 5) - XyzGrid_DI(z_,nFine +1))
          AuxFine = (XyzGrid_DI(z_,5) - Xyz_D(z_))/&
               (XyzGrid_DI(z_,5) - XyzGrid_DI(z_,1))
          nGridOut = 2*nGridOut2
          iOrder_I(nGridOut2 + 1:nGridOut    ) = iOrder_I(5:4 + nGridOut2)
          Weight_I(nGridOut2  +1:nGridOut2+nFine) = &
               Weight_I(1:nFine)*(1 - AuxFine)
          Weight_I(nGridOut2 + nFine+1:nGridOut) = &
               Weight_I(nFine + 1:nGridOut2)*(1 - AuxCoarse)
          Weight_I(1:nFine) = Weight_I(1:nFine)*AuxFine
          Weight_I(nFine + 1:nGridOut2) = &
               Weight_I(nFine + 1:nGridOut2)*AuxCoarse
          IsCorner = .false.
          RETURN
       elseif(UseFaceInsideCorner.and.&
            nCoarse==2.and.nFine==1.and.iLevel_I(iOrder_I(6))==Coarse_.and.&
            iLevel_I(iOrder_I(7))==Coarse_)then
          !\
          !... or we need to choose some unnatural direction for parallel
          !rays interpolation, because the adjacent 'transition junction' was
          !treated as the resolution interface, so that with the 'natural' 
          !choice the continuity would break.
          !/                    View along iDir (point order is set in amr2):
          !\                                   
          !     C---C           C4         C2  Points more F, F1, C3, C2 and 
          !    /|\  |                     /|   their 'out-of-plane' neighbors 
          !   F---F |-->iDir            /  |   form a resolution interface.
          !    \|/  |                 /    |   To keep continuity, interpolate
          !     C---C               F1--\  |   on parallel rays, the direction
          !                              --C3  of rays being perpendicular to 
          !                    more F          that of fake resolution 
          !/                                   interface nearby
          call interpolate_on_parallel_rays(Dir_D = &
               XyzGrid_DI(:,4)-XyzGrid_DI(:,2),&
               iDRectangle1_I = (/2,3,6,7/),&
               iUTriangle1_I = (/1,2, 5/))
          IsCorner = nGridOut < 1
          RETURN
       end if
       !\
       ! Otherwise use parallel-rays procedure with rays parallel to
       ! iDir, that is in the direction toward the corner. We employ the 
       ! property of interpolation to exactly reproduce linear functions 
       ! of coordinates. This allows us to easily find the coordinate of 
       ! the point in which the ray intersects the corner boundary. 
       !/
       dXyzDown = Xyz_D(z_) - &
            sum( XyzGrid_DI(z_,1:nGridOut2)*Weight_I(1:nGridOut2) )
       dXyzUp   = Xyz_D(z_) - (sum(XyzGrid_DI(z_,5:4 + nFine)*&
            Weight_I(1:nFine)) + sum(&
            XyzGrid_DI(z_,nFine + 1:nGridOut2)*Weight_I(nFine + 1:nGridOut2)))

       IsCorner = dXyzUp * dXyzDown > 0.0
       if(IsCorner)RETURN

       if(dXyzUp==0)then
          AuxFine = 0.0
       else
          AuxFine = abs(dXyzUp)/(abs(dXyzUp) + abs(dXyzDown))
       end if
       !\
       ! Remove points, apply weights for fine points. 
       !/
       nGridOut = nGridOut2 + nFine
       iOrder_I(nGridOut2 + 1:nGridOut) = iOrder_I(5:4+nFine)
       Weight_I(nGridOut2 + 1:nGridOut) = Weight_I(1:nFine)*(1 - AuxFine)
       Weight_I(     1:nFine) = Weight_I(1:  nFine)*AuxFine
    case default
       call CON_stop('How could this happen?')
    end select
  contains
    !==========================
    subroutine interpolate_pyramids(&
         iRectangular1_I,iRectangular2_I)
      use ModInterpolateSimpleShape, ONLY:&
           gen_interpolate_pyramids=>interpolate_pyramids
      integer,intent(in),optional,dimension(5)::&
           iRectangular1_I, iRectangular2_I
      integer:: iOrderHere_I(1:5)
      !---------------------------
      call gen_interpolate_pyramids(&
           XyzGrid_DI, Xyz_D,iOrderHere_I, Weight_I, nGridOut,&
           iRectangular1_I=iRectangular1_I,iRectangular2_I=iRectangular2_I)
      if(nGridOut > -1)iOrder_I(1:nGridOut) = &
           iOrder_I(iOrderHere_I(1:nGridOut))
    end subroutine interpolate_pyramids
    !======================
    subroutine interpolate_on_parallel_rays(Dir_D, &
         iDRectangle1_I,  iDTriangle1_I, &
         iURectangle1_I,  iUTriangle1_I)
      use ModInterpolateSimpleShape,ONLY: &
           parallel_rays=>interpolate_on_parallel_rays
      !Direction of parallel rays
      real, intent(in) :: Dir_D(nDim)
      !\
      !Up subfaces ure in the direction of  Dir_D from Xyz point
      !Down faces are in the direction of -Dir_D
      !/
      integer, intent(in), optional ::&
           iDTriangle1_I(3), iUTriangle1_I(3), &
           iDRectangle1_I(4), iURectangle1_I(4)
      integer:: iOrderHere_I(nGrid)
      !----------------------------
      call parallel_rays(&
           XyzGrid_DI, Xyz_D, iOrderHere_I, Weight_I, nGridOut, Dir_D, &
           iDRectangle1_I=iDRectangle1_I, iDTriangle1_I =iDTriangle1_I,  &
           iURectangle1_I=iURectangle1_I, iUTriangle1_I= iUTriangle1_I)
      if(nGridOut > 0)iOrder_I(1:nGridOut) =  &
           iOrder_I(iOrderHere_I(1:nGridOut))
    end subroutine interpolate_on_parallel_rays
    !===========================
  end subroutine interpolate_amr3
end module ModInterpolateAMR
