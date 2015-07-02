!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModTriangulate

  implicit none

  private ! except
  public mesh_triangulation  ! Triangulation of a distorted 2D mesh
  public calc_triangulation  ! Delaunay triangulation of a set of points in 2D
  public find_triangle       ! Find the triangle containing a point
  public triangle_area       ! Returns the area of a triangle
  public test_triangulation  ! Unit test

contains

  !============================================================================
  subroutine mesh_triangulation(n1, n2, CoordXy_DI, &
       iNodeTriangle_II, nTriangle)

    integer, intent(in) :: n1, n2
    real,    intent(in) :: CoordXy_DI(2,n1*n2)
    integer, intent(out):: iNodeTriangle_II(3,2*n1*n2)
    integer, intent(out):: nTriangle

    ! Triangulate a distorted n1*n2 logically Cartesian mesh by splitting
    ! each quadrange along the shorter diagonal. The coordinates are
    ! assumed to be in a natural ordering, with first index changing faster.
    ! Triangle nodes are returned in counter clockwise order assuming that
    ! the distorted mesh is based on a right handed coordinate system.

    integer:: i1, i2, iCell, jCell, kCell, lCell
    !-------------------------------------------------------------------------
    iCell = 0
    nTriangle = 0
    do i2 = 1, n2-1; do i1 = 1, n1
       iCell = iCell + 1

       ! Skip max boundary in the first dimension (unless periodic)
       if(i1 == n1) CYCLE

       ! The vertices around this quadrangle are iCell and
       jCell = iCell + 1
       kCell = iCell + n1
       lCell = iCell + n1 + 1

       if(  sum((CoordXy_DI(:,iCell) - CoordXy_DI(:,lCell))**2) <= &
            sum((CoordXy_DI(:,jCell) - CoordXy_DI(:,kCell))**2))then

          ! Split along iCell--lCell diagonal
          iNodeTriangle_II(:,nTriangle+1) = (/iCell, jCell, lCell/)
          iNodeTriangle_II(:,nTriangle+2) = (/iCell, lCell, kCell/)
       else
          ! Split along jCell--kCell diagonal
          iNodeTriangle_II(:,nTriangle+1) = (/jCell, kCell, iCell/)
          iNodeTriangle_II(:,nTriangle+2) = (/jCell, lCell, kCell/)
       end if
       nTriangle = nTriangle + 2

    end do; end do

  end subroutine mesh_triangulation
  !============================================================================
  subroutine calc_triangulation(nPoint, CoordXy_DI, &
       iNodeTriangle_II, nTriangle)

    ! calc_triangulation: origionally TABLE_DELAUNAY.
    !
    ! Computes the Delaunay triangulation for a set of points in the plane.
    !
    !    Delaunay triangulation. The Delaunay triangulation is an organization
    !    of the data into triples, forming a triangulation of the data, with
    !    the property that the circumcircle of each triangle never contains
    !    another data point.  
    !
    !  Author: John Burkardt, Modified by Alex Glocer:2/2007

    implicit none

    integer, intent(in) :: nPoint
    real, intent(in)    :: CoordXy_DI(2,nPoint)
    integer, intent(out):: iNodeTriangle_II(3,2*nPoint)

    integer, allocatable:: WorkArray_II(:,:)
    integer,intent(out) :: nTriangle
    !--------------------------------------------------------------------------
    allocate(WorkArray_II(3,3*nPoint))

    call dtris2(nPoint, CoordXy_DI, nTriangle, iNodeTriangle_II, WorkArray_II)
    
    deallocate(WorkArray_II)

  end subroutine calc_triangulation

  !============================================================================

  subroutine find_triangle(&
       nPoint, nTriangle, CoordIn_D, XyNode_DI, iNodeTriangle_II,&
       iNode1, iNode2, iNode3, Weight1, Weight2, Weight3, IsTriangleFound, &
       nTriangle_C, iTriangle_IC)

    ! Determine the triangle containing a given point. Calculate the
    ! weights for the 3 nodes if the weights are present.
    ! If the point is not contained in any triangle, then either
    ! return IsTriangleFound false if the argument is present,
    ! or the closest node as iNode1 (and Weight1 = 1.0 if present)

    ! created by Alex Glocer, 01/2007
    ! extended with weights by G. Toth, 10/2008

    integer,intent(in)   :: nPoint,nTriangle
    real,   intent(in)   :: CoordIN_D(2)
    real,   intent(in)   :: XyNode_DI(2,nPoint)
    integer,intent(in)   :: iNodeTriangle_II(3,nTriangle)

    integer,intent(out)  :: iNode1, iNode2, iNode3

    real,    optional, intent(out):: Weight1, Weight2, Weight3
    logical, optional, intent(out):: IsTriangleFound
    
    ! Number and indexes of triangles intersecting cells in a uniform mesh
    integer, optional, intent(inout):: nTriangle_C(:,:)
    integer, optional, pointer      :: iTriangle_IC(:,:,:)

    integer, parameter:: x_=1, y_=2

    integer :: iTriangle

    ! Size of optional lookup mesh
    integer:: nX, nY
    ! Limits of domain and resolution of mesh
    real   :: xMin, xMax, yMin, yMax, DxInv, DyInv

    integer:: iNode_I(3)
    real   :: x_I(3), y_I(3)
    integer:: i, j, n, iLoop, iMin, iMax, jMin, jMax
    !--------------------------------------------------------------------------

    if(present(nTriangle_C) .and. present(iTriangle_IC))then

       nX    = size(nTriangle_C, DIM=1)
       nY    = size(nTriangle_C, DIM=2)

       xMin  = minval(XyNode_DI(1,:))
       xMax  = maxval(XyNode_DI(1,:))
       yMin  = minval(XyNode_DI(2,:))
       yMax  = maxval(XyNode_DI(2,:))
       DxInv = nX/(xMax - xMin)
       DyInv = nY/(yMax - yMin)

       if(nTriangle_C(1,1)<0)then
          do iLoop = 1, 2
             nTriangle_C = 0
             do iTriangle = 1, nTriangle
                iNode_I = iNodeTriangle_II(:,iTriangle)
                x_I   = XyNode_DI(1,iNode_I)
                y_I   = XyNode_DI(2,iNode_I)
                iMin = max(1,  floor(DxInv*(minval(x_I) - xMin)) + 1)
                iMax = min(nX, floor(DxInv*(maxval(x_I) - xMin)) + 1)
                jMin = max(1,  floor(DyInv*(minval(y_I) - yMin)) + 1)
                jMax = min(nY, floor(DyInv*(maxval(y_I) - yMin)) + 1)

                nTriangle_C(iMin:iMax,jMin:jMax) = &
                     nTriangle_C(iMin:iMax,jMin:jMax) + 1

                if(iLoop == 1)CYCLE
                
                do j = jMin, jMax; do i = iMin, iMax
                   n = nTriangle_C(i,j)
                   iTriangle_IC(n,i,j) = iTriangle
                end do; end do
             end do
             if(iLoop == 1) then
                n = maxval(nTriangle_C)
                allocate(iTriangle_IC(n,nX,nY))
             end if
          end do
       end if

       i = max(1, min(nX, floor(DxInv*(CoordIn_D(1) - xMin)) + 1))
       j = max(1, min(nY, floor(DyInv*(CoordIn_D(2) - yMin)) + 1))

       do n = 1, nTriangle_C(i,j)
          iTriangle = iTriangle_IC(n,i,j)
          if(is_inside_triangle()) RETURN
       end do

    else
       do iTriangle = 1, nTriangle 
          if(is_inside_triangle()) RETURN
       end do
    end if

    ! If no triangle contains the CoordIn_D point and IsTriangleFound is
    ! present then return a false value
    if(present(IsTriangleFound))then
       IsTriangleFound = .false.
       RETURN
    end if

    ! If IsTriangleFound is not present then set iNode1 to the closest node
    iNode1 = minloc( (XyNode_DI(x_,:)-CoordIn_D(x_))**2 &
         +           (XyNode_DI(y_,:)-CoordIn_D(y_))**2 , DIM=1)

    ! Set iNode2 and iNode3 to valid values so that general interpolation
    ! formulas can work but assign same index to signal that there is no
    ! triangle containing the point
    iNode2 = 1
    iNode3 = 1

    if(present(Weight1))then
       ! Assign full weight to the closest node and zero to the others
       Weight1 = 1.0
       Weight2 = 0.0
       Weight3 = 0.0
    end if

  contains
    !=======================================================================
    logical function is_inside_triangle()

      real    :: x1, x2, x3, y1, y2, y3
      real    :: Area1, Area2, Area3, AreaSum
      
      iNode1 = iNodeTriangle_II(1,iTriangle)
      iNode2 = iNodeTriangle_II(2,iTriangle)
      iNode3 = iNodeTriangle_II(3,iTriangle)
      !--------------------------------------------------------------------
      ! fill in x,y pairs for current triangle nodes
      x1 = XyNode_DI(1,iNode1)
      y1 = XyNode_DI(2,iNode1)

      x2 = XyNode_DI(1,iNode2)
      y2 = XyNode_DI(2,iNode2)

      x3 = XyNode_DI(1,iNode3)
      y3 = XyNode_DI(2,iNode3)

      ! Check if triangle contains point or if point lies on edge.
      ! Area1 is twice the area of the triangle formed by Node2, Node3 
      ! and the point. Area1 is zero if the point lies on the Node2-Node3
      ! edge and positive if point lies to the left of the Node2-Node3 vector.

      Area1 = (CoordIn_D(y_)-y2)*(x3-x2) - (CoordIn_D(x_)-x2)*(y3-y2)
      Area2 = (CoordIn_D(y_)-y3)*(x1-x3) - (CoordIn_D(x_)-x3)*(y1-y3)
      Area3 = (CoordIn_D(y_)-y1)*(x2-x1) - (CoordIn_D(x_)-x1)*(y2-y1)

      ! The point is inside the triangle if all areas are positive
      if (Area1 >= 0.0 .and. Area2 >= 0.0 .and. Area3 >=0.0) then
         if(present(IsTriangleFound)) IsTriangleFound = .true.
         if(present(Weight1)) then
            ! The node weights are proportional to the area of 
            ! the opposite triangle.
            AreaSum = max(Area1 + Area2 + Area3, tiny(1.0))
            Weight1 = Area1/AreaSum
            Weight2 = Area2/AreaSum
            Weight3 = Area3/AreaSum
          end if
          is_inside_triangle = .true.
       else
          is_inside_triangle = .false.
       end if

     end function is_inside_triangle

  end subroutine find_triangle

  !============================================================================

  real function triangle_area(Node_DI)

    real,intent(in)       :: Node_DI(2,3)
    real :: a_D(2), b_D(2)
    !-------------------------------------------------------------------------
    a_D = Node_DI(:,1)-Node_DI(:,2)
    b_D = Node_DI(:,1)-Node_DI(:,3)
    triangle_area = 0.5*abs(a_D(1)*b_D(2)-a_D(2)*b_D(1))

  end function triangle_area

  !============================================================================

  integer function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )

    !**************************************************************************
    !
    !! DIAEDG chooses a diagonal edge.
    !
    !  Discussion:
    !
    !    The routine determines whether 0--2 or 1--3 is the diagonal edge
    !    that should be chosen, based on the circumcircle criterion, where
    !    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
    !    quadrilateral in counterclockwise order.
    !
    !  Modified:
    !
    !    19 February 2001
    !
    !  Author:
    !
    !    Barry Joe,
    !    Department of Computing Science,
    !    University of Alberta,
    !    Edmonton, Alberta, Canada  T6G 2H1
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !      using geometric algorithms, 
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, real  X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
    !    coordinates of the vertices of a quadrilateral, given in
    !    counter clockwise order.
    !
    !    Output, integer DIAEDG, chooses a diagonal:
    !    +1, if diagonal edge 02 is chosen;
    !    -1, if diagonal edge 13 is chosen;
    !     0, if the four vertices are cocircular.
    !
    implicit none

    real  ca
    real  cb
    real  dx10
    real  dx12
    real  dx30
    real  dx32
    real  dy10
    real  dy12
    real  dy30
    real  dy32
    real  s
    real  tol
    real  tola
    real  tolb
    real  x0
    real  x1
    real  x2
    real  x3
    real  y0
    real  y1
    real  y2
    real  y3

    tol = 100.0 * epsilon ( tol )

    dx10 = x1 - x0
    dy10 = y1 - y0
    dx12 = x1 - x2
    dy12 = y1 - y2
    dx30 = x3 - x0
    dy30 = y3 - y0
    dx32 = x3 - x2
    dy32 = y3 - y2

    tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
    tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

    ca = dx10 * dx30 + dy10 * dy30
    cb = dx12 * dx32 + dy12 * dy32

    if ( tola < ca .and. tolb < cb ) then

       diaedg = -1

    else if ( ca < -tola .and. cb < -tolb ) then

       diaedg = 1

    else

       tola = max ( tola, tolb )
       s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

       if ( tola < s ) then
          diaedg = -1
       else if ( s < -tola ) then
          diaedg = 1
       else
          diaedg = 0
       end if

    end if

    return
  end function diaedg

  subroutine dtris2 ( point_num, point_xy, tri_num, tri_vert, tri_nabe )

    !**************************************************************************
    !
    !! DTRIS2 constructs a Delaunay triangulation of 2D vertices.
    !
    !  Discussion:
    !
    ! The routine constructs the Delaunay triangulation of a set of 2D vertices
    ! using an incremental approach and diagonal edge swaps.  Vertices are
    ! first sorted in lexicographically increasing (X,Y) order, and
    ! then are inserted one at a time from outside the convex hull.
    !
    !  Modified:
    !
    !    25 August 2001
    !
    !  Author:
    !
    !    Barry Joe,
    !    Department of Computing Science,
    !    University of Alberta,
    !    Edmonton, Alberta, Canada  T6G 2H1
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !      using geometric algorithms, 
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, integer POINT_NUM, the number of vertices.
    !
    !    Input/output, real  POINT_XY(2,POINT_NUM), the coordinates 
    !    of the vertices.  On output, the vertices have been sorted into 
    !    dictionary order.
    !
    !    Output, integer TRI_NUM, the number of triangles in the triangulation;
    !    TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the number
    !    of boundary vertices.
    !
    !Output, integer TRI_VERT(3,TRI_NUM), the nodes that make up each triangle.
    !The elements are indices of POINT_XY.  The vertices of the triangles are
    !in counter clockwise order.
    !
    !Output, integer TRI_NABE(3,TRI_NUM), the triangle neighbor list.
    !Positive elements are indices of TIL; negative elements are used for links
    !of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
    !where I, J = triangle, edge index; TRI_NABE(J,I) refers to
    !the neighbor along edge from vertex J to J+1 (mod 3).
    !
    implicit none

    integer point_num

    real  cmax
    integer e
    integer i
    integer ierr
    integer indx(point_num)
    integer j
    integer k
    integer l
    integer ledg
    integer lr
    integer ltri
    integer m
    integer m1
    integer m2
    integer n
    real  point_xy(2,point_num)
    integer redg
    integer rtri
    integer stack(point_num)
    integer t
    real  tol
    integer top
    integer tri_nabe(3,point_num*2)
    integer tri_num
    integer tri_vert(3,point_num*2)

    tol = 100.0 * epsilon ( tol )

    ierr = 0
    !
    !  Sort the vertices by increasing (x,y).
    !
    call r82vec_sort_heap_index_a ( point_num, point_xy, indx )

    call r82vec_permute ( point_num, point_xy, indx )
    !
    !  Make sure that the data points are "reasonably" distinct.
    !
    m1 = 1

    do i = 2, point_num

       m = m1
       m1 = i

       k = 0

       do j = 1, 2

          cmax = max ( abs ( point_xy(j,m) ), abs ( point_xy(j,m1) ) )

          if ( tol * ( cmax + 1.0 ) &
               < abs ( point_xy(j,m) - point_xy(j,m1) ) ) then
             k = j
             exit
          end if

       end do

       if ( k == 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
          write ( *, '(a,i8)' ) '  Fails for point number I = ', i
          write ( *, '(a,i8)' ) '  M = ', m
          write ( *, '(a,i8)' ) '  M1 = ', m1
          write ( *, '(a,2g14.6)' ) '  X,Y(M)  = ', point_xy(1,m), point_xy(2,m)
          write ( *, '(a,2g14.6)' ) '  X,Y(M1) = ', point_xy(1,m1), point_xy(2,m1)
          ierr = 224
          return
       end if

    end do
    !
    !  Starting from points M1 and M2, search for a third point M that
    !  makes a "healthy" triangle (M1,M2,M)
    !
    m1 = 1
    m2 = 2
    j = 3

    do

       if ( point_num < j ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
          ierr = 225
          return
       end if

       m = j

       lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
            point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0 )

       if ( lr /= 0 ) then
          exit
       end if

       j = j + 1

    end do
    !
    !  Set up the triangle information for (M1,M2,M), and for any other
    !  triangles you created because points were collinear with M1, M2.
    !
    tri_num = j - 2

    if ( lr == -1 ) then

       tri_vert(1,1) = m1
       tri_vert(2,1) = m2
       tri_vert(3,1) = m
       tri_nabe(3,1) = -3

       do i = 2, tri_num

          m1 = m2
          m2 = i+1
          tri_vert(1,i) = m1
          tri_vert(2,i) = m2
          tri_vert(3,i) = m
          tri_nabe(1,i-1) = -3 * i
          tri_nabe(2,i-1) = i
          tri_nabe(3,i) = i - 1

       end do

       tri_nabe(1,tri_num) = -3 * tri_num - 1
       tri_nabe(2,tri_num) = -5
       ledg = 2
       ltri = tri_num

    else

       tri_vert(1,1) = m2
       tri_vert(2,1) = m1
       tri_vert(3,1) = m
       tri_nabe(1,1) = -4

       do i = 2, tri_num
          m1 = m2
          m2 = i+1
          tri_vert(1,i) = m2
          tri_vert(2,i) = m1
          tri_vert(3,i) = m
          tri_nabe(3,i-1) = i
          tri_nabe(1,i) = -3 * i - 3
          tri_nabe(2,i) = i - 1
       end do

       tri_nabe(3,tri_num) = -3 * tri_num
       tri_nabe(2,1) = -3 * tri_num - 2
       ledg = 2
       ltri = 1

    end if
    !
    !  Insert the vertices one at a time from outside the convex hull,
    !  determine visible boundary edges, and apply diagonal edge swaps until
    !  Delaunay triangulation of vertices (so far) is obtained.
    !
    top = 0

    do i = j+1, point_num

       m = i
       m1 = tri_vert(ledg,ltri)

       if ( ledg <= 2 ) then
          m2 = tri_vert(ledg+1,ltri)
       else
          m2 = tri_vert(1,ltri)
       end if

       lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
            point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0 )

       if ( 0 < lr ) then
          rtri = ltri
          redg = ledg
          ltri = 0
       else
          l = -tri_nabe(ledg,ltri)
          rtri = l / 3
          redg = mod(l,3) + 1
       end if

       call vbedg ( point_xy(1,m), point_xy(2,m), point_num, point_xy, tri_num, &
            tri_vert, tri_nabe, ltri, ledg, rtri, redg )

       n = tri_num + 1
       l = -tri_nabe(ledg,ltri)

       do

          t = l / 3
          e = mod ( l, 3 ) + 1
          l = -tri_nabe(e,t)
          m2 = tri_vert(e,t)

          if ( e <= 2 ) then
             m1 = tri_vert(e+1,t)
          else
             m1 = tri_vert(1,t)
          end if

          tri_num = tri_num + 1
          tri_nabe(e,t) = tri_num
          tri_vert(1,tri_num) = m1
          tri_vert(2,tri_num) = m2
          tri_vert(3,tri_num) = m
          tri_nabe(1,tri_num) = t
          tri_nabe(2,tri_num) = tri_num - 1
          tri_nabe(3,tri_num) = tri_num + 1
          top = top + 1

          if ( point_num < top ) then
             ierr = 8
             write ( *, '(a)' ) ' '
             write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
             write ( *, '(a)' ) '  Stack overflow.'
             return
          end if

          stack(top) = tri_num

          if ( t == rtri .and. e == redg ) then
             exit
          end if

       end do

       tri_nabe(ledg,ltri) = -3 * n - 1
       tri_nabe(2,n) = -3 * tri_num - 2
       tri_nabe(3,tri_num) = -l
       ltri = n
       ledg = 2

       call swapec ( m, top, ltri, ledg, point_num, point_xy, tri_num, &
            tri_vert, tri_nabe, stack, ierr )

       if ( ierr /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
          write ( *, '(a)' ) '  Error return from SWAPEC.'
          return
       end if

    end do
    !
    !  Now account for the sorting that we did.
    !
    do i = 1, 3
       do j = 1, tri_num
          tri_vert(i,j) = indx ( tri_vert(i,j) )
       end do
    end do

    call perm_inv ( point_num, indx )

    call r82vec_permute ( point_num, point_xy, indx )

    return
  end subroutine dtris2





  integer function i4_modp ( i, j )

    !*******************************************************************************
    !
    !! I4_MODP returns the nonnegative remainder of integer division.
    !
    !  Formula:
    !
    !    If
    !      NREM = I4_MODP ( I, J )
    !      NMULT = ( I - NREM ) / J
    !    then
    !      I = J * NMULT + NREM
    !    where NREM is always nonnegative.
    !
    !  Comments:
    !
    !    The MOD function computes a result with the same sign as the
    !    quantity being divided.  Thus, suppose you had an angle A,
    !    and you wanted to ensure that it was between 0 and 360.
    !    Then mod(A,360) would do, if A was positive, but if A
    !    was negative, your result would be between -360 and 0.
    !
    !    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
    !
    !  Examples:
    !
    !        I     J     MOD  I4_MODP    Factorization
    !
    !      107    50       7       7    107 =  2 *  50 + 7
    !      107   -50       7       7    107 = -2 * -50 + 7
    !     -107    50      -7      43   -107 = -3 *  50 + 43
    !     -107   -50      -7      43   -107 =  3 * -50 + 43
    !
    !  Modified:
    !
    !    02 March 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer I, the number to be divided.
    !
    !    Input, integer J, the number that divides I.
    !
    !    Output, integer I4_MODP, the nonnegative remainder when I is
    !    divided by J.
    !
    implicit none

    integer i
    integer j

    if ( j == 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'I4_MODP - Fatal error!'
       write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
       call CON_stop('share/Library/src/ModTriangulate ERROR')
    end if

    i4_modp = mod ( i, j )

    if ( i4_modp < 0 ) then
       i4_modp = i4_modp + abs ( j )
    end if

    return
  end function i4_modp
  integer function i4_wrap ( ival, ilo, ihi )

    !*******************************************************************************
    !
    !! I4_WRAP forces an integer to lie between given limits by wrapping.
    !
    !  Example:
    !
    !    ILO = 4, IHI = 8
    !
    !    I  I4_WRAP
    !
    !    -2     8
    !    -1     4
    !     0     5
    !     1     6
    !     2     7
    !     3     8
    !     4     4
    !     5     5
    !     6     6
    !     7     7
    !     8     8
    !     9     4
    !    10     5
    !    11     6
    !    12     7
    !    13     8
    !    14     4
    !
    !  Modified:
    !
    !    19 August 2003
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer IVAL, an integer value.
    !
    !    Input, integer ILO, IHI, the desired bounds for the integer value.
    !
    !    Output, integer I4_WRAP, a "wrapped" version of IVAL.
    !
    implicit none

    integer ihi
    integer ilo
    integer ival
    integer jhi
    integer jlo
    integer wide

    jlo = min ( ilo, ihi )
    jhi = max ( ilo, ihi )

    wide = jhi - jlo + 1

    if ( wide == 1 ) then
       i4_wrap = jlo
    else
       i4_wrap = jlo + i4_modp ( ival - jlo, wide )
    end if

    return
  end function i4_wrap



  subroutine i4vec_indicator ( n, a )

    !*******************************************************************************
    !
    !! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
    !
    !  Modified:
    !
    !    09 November 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer N, the number of elements of A.
    !
    !    Output, integer A(N), the array to be initialized.
    !
    implicit none

    integer n

    integer a(n)
    integer i

    do i = 1, n
       a(i) = i
    end do

    return
  end subroutine i4vec_indicator
  subroutine i4vec_sort_heap_index_a ( n, a, indx )

    !*******************************************************************************
    !
    !! I4VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an I4VEC.
    !
    !  Discussion:
    !
    !    The sorting is not actually carried out.  Rather an index array is
    !    created which defines the sorting.  This array may be used to sort
    !    or index the array, or to sort or index related arrays keyed on the
    !    original array.
    !
    !    Once the index array is computed, the sorting can be carried out
    !    "implicitly:
    !
    !      A(INDX(I)), I = 1 to N is sorted,
    !
    !    or explicitly, by the call
    !
    !      call I4VEC_PERMUTE ( N, A, INDX )
    !
    !    after which A(I), I = 1 to N is sorted.
    !
    !  Modified:
    !
    !    25 September 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer N, the number of entries in the array.
    !
    !    Input, integer A(N), an array to be index-sorted.
    !
    !    Output, integer INDX(N), the sort index.  The
    !    I-th element of the sorted array is A(INDX(I)).
    !
    implicit none

    integer n

    integer a(n)
    integer aval
    integer i
    integer indx(n)
    integer indxt
    integer ir
    integer j
    integer l

    if ( n <= 1 ) then
       return
    end if

    do i = 1, n
       indx(i) = i
    end do

    l = n / 2 + 1
    ir = n

    do

       if ( 1 < l ) then

          l = l - 1
          indxt = indx(l)
          aval = a(indxt)

       else

          indxt = indx(ir)
          aval = a(indxt)
          indx(ir) = indx(1)
          ir = ir - 1

          if ( ir == 1 ) then
             indx(1) = indxt
             exit
          end if

       end if

       i = l
       j = l + l

       do while ( j <= ir )

          if ( j < ir ) then
             if ( a(indx(j)) < a(indx(j+1)) ) then
                j = j + 1
             end if
          end if

          if ( aval < a(indx(j)) ) then
             indx(i) = indx(j)
             i = j
             j = j + j
          else
             j = ir + 1
          end if

       end do

       indx(i) = indxt

    end do

    return
  end subroutine i4vec_sort_heap_index_a
  integer function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

    !*******************************************************************************
    !
    !! LRLINE determines if a point is left of, right or, or on a directed line.
    !
    !  Discussion:
    !
    !    The directed line is parallel to, and at a signed distance DV from
    !    a directed base line from (XV1,YV1) to (XV2,YV2).
    !
    !  Modified:
    !
    !    14 July 2001
    !
    !  Author:
    !
    !    Barry Joe,
    !    Department of Computing Science,
    !    University of Alberta,
    !    Edmonton, Alberta, Canada  T6G 2H1
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, real  XU, YU, the coordinates of the point whose
    !    position relative to the directed line is to be determined.
    !
    !    Input, real  XV1, YV1, XV2, YV2, the coordinates of two points
    !    that determine the directed base line.
    !
    !    Input, real  DV, the signed distance of the directed line
    !    from the directed base line through the points (XV1,YV1) and (XV2,YV2).
    !    DV is positive for a line to the left of the base line.
    !
    !    Output, integer LRLINE, the result:
    !    +1, the point is to the right of the directed line;
    !     0, the point is on the directed line;
    !    -1, the point is to the left of the directed line.
    !
    implicit none

    real  dv
    real  dx
    real  dxu
    real  dy
    real  dyu
    real  t
    real  tol
    real  tolabs
    real  xu
    real  xv1
    real  xv2
    real  yu
    real  yv1
    real  yv2

    tol = 100.0 * epsilon ( tol )

    dx = xv2 - xv1
    dy = yv2 - yv1
    dxu = xu - xv1
    dyu = yu - yv1

    tolabs = tol * max ( abs ( dx ), abs ( dy ), abs ( dxu ), &
         abs ( dyu ), abs ( dv ) )

    t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy )

    if ( tolabs < t ) then
       lrline = 1
    else if ( -tolabs <= t ) then
       lrline = 0
    else
       lrline = -1
    end if

    return
  end function lrline
  subroutine perm_inv ( n, p )

    !*******************************************************************************
    !
    !! PERM_INV inverts a permutation "in place".
    !
    !  Modified:
    !
    !    25 July 2000
    !
    !  Parameters:
    !
    !    Input, integer N, the number of objects being permuted.
    !
    !    Input/output, integer P(N), the permutation, in standard index form.
    !    On output, P describes the inverse permutation
    !
    implicit none

    integer n

    integer i
    integer i0
    integer i1
    integer i2
    integer is
    integer p(n)

    if ( n <= 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'PERM_INV - Fatal error!'
       write ( *, '(a,i8)' ) '  Input value of N = ', n
       call CON_stop('share/Library/src/ModTriangulate ERROR')
    end if

    is = 1

    do i = 1, n

       i1 = p(i)

       do while ( i < i1 )
          i2 = p(i1)
          p(i1) = -i2
          i1 = i2
       end do

       is = -sign ( 1, p(i) )
       p(i) = sign ( p(i), is )

    end do

    do i = 1, n

       i1 = -p(i)

       if ( 0 <= i1 ) then

          i0 = i

          do

             i2 = p(i1)
             p(i1) = i0

             if ( i2 < 0 ) then
                exit
             end if

             i0 = i1
             i1 = i2

          end do

       end if

    end do

    return
  end subroutine perm_inv
  subroutine r82vec_permute ( n, a, p )

    !*******************************************************************************
    !
    !! R82VEC_PERMUTE permutes an R82VEC in place.
    !
    !  Discussion:
    !
    !    This routine permutes an array of real "objects", but the same
    !    logic can be used to permute an array of objects of any arithmetic
    !    type, or an array of objects of any complexity.  The only temporary
    !    storage required is enough to store a single object.  The number
    !    of data movements made is N + the number of cycles of order 2 or more,
    !    which is never more than N + N/2.
    !
    !  Example:
    !
    !    Input:
    !
    !      N = 5
    !      P = (   2,    4,    5,    1,    3 )
    !      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
    !          (11.0, 22.0, 33.0, 44.0, 55.0 )
    !
    !    Output:
    !
    !      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
    !             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
    !
    !  Modified:
    !
    !    11 January 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer N, the number of objects.
    !
    !    Input/output, real  A(2,N), the array to be permuted.
    !
    !    Input, integer P(N), the permutation.  P(I) = J means
    !    that the I-th element of the output array should be the J-th
    !    element of the input array.  P must be a legal permutation
    !    of the integers from 1 to N, otherwise the algorithm will
    !    fail catastrophically.
    !
    implicit none

    integer n

    real  a(2,n)
    real  a_temp(2)
    integer iget
    integer iput
    integer istart
    integer p(n)
    !
    !  Search for the next element of the permutation that has not been used.
    !
    do istart = 1, n

       if ( p(istart) < 0 ) then

          cycle

       else if ( p(istart) == istart ) then

          p(istart) = - p(istart)
          cycle

       else

          a_temp(1:2) = a(1:2,istart)
          iget = istart
          !
          !  Copy the new value into the vacated entry.
          !
          do

             iput = iget
             iget = p(iget)

             p(iput) = - p(iput)

             if ( iget < 1 .or. n < iget ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
                call CON_stop('share/Library/src/ModTriangulate ERROR')
             end if

             if ( iget == istart ) then
                a(1:2,iput) = a_temp(1:2)
                exit
             end if

             a(1:2,iput) = a(1:2,iget)

          end do

       end if

    end do
    !
    !  Restore the signs of the entries.
    !
    p(1:n) = -p(1:n)

    return
  end subroutine r82vec_permute
  subroutine r82vec_sort_heap_index_a ( n, a, indx )

    !*******************************************************************************
    !
    !! R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.
    !
    !  Discussion:
    !
    !    The sorting is not actually carried out.  Rather an index array is
    !    created which defines the sorting.  This array may be used to sort
    !    or index the array, or to sort or index related arrays keyed on the
    !    original array.
    !
    !    Once the index array is computed, the sorting can be carried out
    !    "implicitly:
    !
    !      A(1:2,INDX(I)), I = 1 to N is sorted,
    !
    !    or explicitly, by the call
    !
    !      call R82VEC_PERMUTE ( N, A, INDX )
    !
    !    after which A(1:2,I), I = 1 to N is sorted.
    !
    !  Modified:
    !
    !    11 January 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer N, the number of entries in the array.
    !
    !    Input, real  A(2,N), an array to be index-sorted.
    !
    !    Output, integer INDX(N), the sort index.  The
    !    I-th element of the sorted array is A(1:2,INDX(I)).
    !
    implicit none

    integer n

    real  a(2,n)
    real  aval(2)
    integer i
    integer indx(n)
    integer indxt
    integer ir
    integer j
    integer l

    if ( n < 1 ) then
       return
    end if

    if ( n == 1 ) then
       indx(1) = 1
       return
    end if

    call i4vec_indicator ( n, indx )

    l = n / 2 + 1
    ir = n

    do

       if ( 1 < l ) then

          l = l - 1
          indxt = indx(l)
          aval(1:2) = a(1:2,indxt)

       else

          indxt = indx(ir)
          aval(1:2) = a(1:2,indxt)
          indx(ir) = indx(1)
          ir = ir - 1

          if ( ir == 1 ) then
             indx(1) = indxt
             exit
          end if

       end if

       i = l
       j = l + l

       do while ( j <= ir )

          if ( j < ir ) then
             if (   a(1,indx(j)) <  a(1,indx(j+1)) .or. &
                  ( a(1,indx(j)) == a(1,indx(j+1)) .and. &
                  a(2,indx(j)) <  a(2,indx(j+1)) ) ) then
                j = j + 1
             end if
          end if

          if (   aval(1) <  a(1,indx(j)) .or. &
               ( aval(1) == a(1,indx(j)) .and. &
               aval(2) <  a(2,indx(j)) ) ) then
             indx(i) = indx(j)
             i = j
             j = j + j
          else
             j = ir + 1
          end if

       end do

       indx(i) = indxt

    end do

    return
  end subroutine r82vec_sort_heap_index_a



  subroutine swapec ( i, top, btri, bedg, point_num, point_xy, tri_num, &
       tri_vert, tri_nabe, stack, ierr )

    !*******************************************************************************
    !
    !! SWAPEC swaps diagonal edges until all triangles are Delaunay.
    !
    !  Discussion:
    !
    !    The routine swaps diagonal edges in a 2D triangulation, based on
    !    the empty circumcircle criterion, until all triangles are Delaunay,
    !    given that I is the index of the new vertex added to the triangulation.
    !
    !  Modified:
    !
    !    14 July 2001
    !
    !  Author:
    !
    !    Barry Joe,
    !    Department of Computing Science,
    !    University of Alberta,
    !    Edmonton, Alberta, Canada  T6G 2H1
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !      using geometric algorithms, 
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, integer I, the index of the new vertex.
    !
    !    Input/output, integer TOP, the index of the top of the stack.
    !    On output, TOP is zero.
    !
    !    Input/output, integer BTRI, BEDG; on input, if positive, are the
    !    triangle and edge indices of a boundary edge whose updated indices
    !    must be recorded.  On output, these may be updated because of swaps.
    !
    !    Input, intger POINT_NUM, the number of points.
    !
    !    Input, real  POINT_XY(2,POINT_NUM), the coordinates
    !    of the points.
    !
    !    Input, integer TRI_NUM, the number of triangles.
    !
    !    Input/output, integer TRI_VERT(3,TRI_NUM), the triangle incidence list.
    !    May be updated on output because of swaps.
    !
    !    Input/output, integer TRI_NABE(3,TRI_NUM), the triangle neighbor list;
    !    negative values are used for links of the counter-clockwise linked
    !    list of boundary edges;  May be updated on output because of swaps.
    !
    !      LINK = -(3*I + J-1) where I, J = triangle, edge index.
    !
    !    Workspace, integer STACK(MAXST); on input, entries 1 through TOP
    !    contain the indices of initial triangles (involving vertex I)
    !    put in stack; the edges opposite I should be in interior;  entries
    !    TOP+1 through MAXST are used as a stack.
    !
    !    Output, integer IERR is set to 8 for abnormal return.
    !
    implicit none

    integer point_num
    integer tri_num

    integer a
    integer b
    integer bedg
    integer btri
    integer c
    integer e
    integer ee
    integer em1
    integer ep1
    integer f
    integer fm1
    integer fp1
    integer i
    integer ierr
    integer l
    integer r
    integer s
    integer stack(point_num)
    integer swap
    integer t
    integer top
    integer tri_nabe(3,tri_num)
    integer tri_vert(3,tri_num)
    integer tt
    integer u
    real  point_xy(2,point_num)
    real  x
    real  y
    !
    !  Determine whether triangles in stack are Delaunay, and swap
    !  diagonal edge of convex quadrilateral if not.
    !
    x = point_xy(1,i)
    y = point_xy(2,i)

    do

       if ( top <= 0 ) then
          exit
       end if

       t = stack(top)
       top = top - 1

       if ( tri_vert(1,t) == i ) then
          e = 2
          b = tri_vert(3,t)
       else if ( tri_vert(2,t) == i ) then
          e = 3
          b = tri_vert(1,t)
       else
          e = 1
          b = tri_vert(2,t)
       end if

       a = tri_vert(e,t)
       u = tri_nabe(e,t)

       if ( tri_nabe(1,u) == t ) then
          f = 1
          c = tri_vert(3,u)
       else if ( tri_nabe(2,u) == t ) then
          f = 2
          c = tri_vert(1,u)
       else
          f = 3
          c = tri_vert(2,u)
       end if

       swap = diaedg ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,c), &
            point_xy(2,c), point_xy(1,b), point_xy(2,b) )

       if ( swap == 1 ) then

          em1 = i4_wrap ( e - 1, 1, 3 )
          ep1 = i4_wrap ( e + 1, 1, 3 )
          fm1 = i4_wrap ( f - 1, 1, 3 )
          fp1 = i4_wrap ( f + 1, 1, 3 )

          tri_vert(ep1,t) = c
          tri_vert(fp1,u) = i
          r = tri_nabe(ep1,t)
          s = tri_nabe(fp1,u)
          tri_nabe(ep1,t) = u
          tri_nabe(fp1,u) = t
          tri_nabe(e,t) = s
          tri_nabe(f,u) = r

          if ( 0 < tri_nabe(fm1,u) ) then
             top = top + 1
             stack(top) = u
          end if

          if ( 0 < s ) then

             if ( tri_nabe(1,s) == u ) then
                tri_nabe(1,s) = t
             else if ( tri_nabe(2,s) == u ) then
                tri_nabe(2,s) = t
             else
                tri_nabe(3,s) = t
             end if

             top = top + 1

             if ( point_num < top ) then
                ierr = 8
                return
             end if

             stack(top) = t

          else

             if ( u == btri .and. fp1 == bedg ) then
                btri = t
                bedg = e
             end if

             l = - ( 3 * t + e - 1 )
             tt = t
             ee = em1

             do while ( 0 < tri_nabe(ee,tt) )

                tt = tri_nabe(ee,tt)

                if ( tri_vert(1,tt) == a ) then
                   ee = 3
                else if ( tri_vert(2,tt) == a ) then
                   ee = 1
                else
                   ee = 2
                end if

             end do

             tri_nabe(ee,tt) = l

          end if

          if ( 0 < r ) then

             if ( tri_nabe(1,r) == t ) then
                tri_nabe(1,r) = u
             else if ( tri_nabe(2,r) == t ) then
                tri_nabe(2,r) = u
             else
                tri_nabe(3,r) = u
             end if

          else

             if ( t == btri .and. ep1 == bedg ) then
                btri = u
                bedg = f
             end if

             l = - ( 3 * u + f - 1 )
             tt = u
             ee = fm1

             do while ( 0 < tri_nabe(ee,tt) )

                tt = tri_nabe(ee,tt)

                if ( tri_vert(1,tt) == b ) then
                   ee = 3
                else if ( tri_vert(2,tt) == b ) then
                   ee = 1
                else
                   ee = 2
                end if

             end do

             tri_nabe(ee,tt) = l

          end if

       end if

    end do

    return
  end subroutine swapec


  subroutine vbedg ( x, y, point_num, point_xy, tri_num, tri_vert, tri_nabe, &
       ltri, ledg, rtri, redg )

    !*******************************************************************************
    !
    !! VBEDG determines which boundary edges are visible to a point.
    !
    !  Discussion:
    !
    !    The point (X,Y) is assumed to be outside the convex hull of the
    !    region covered by the 2D triangulation.
    !
    !  Modified:
    !
    !    25 August 2001
    !
    !  Author:
    !
    !    Barry Joe,
    !    Department of Computing Science,
    !    University of Alberta,
    !    Edmonton, Alberta, Canada  T6G 2H1
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !      using geometric algorithms, 
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, real  X, Y, the coordinates of a point outside
    !    the convex hull of the current triangulation.
    !
    !    Input, integer POINT_NUM, the number of points.
    !
    !    Input, real  POINT_XY(2,POINT_NUM), the coordinates 
    !    of the vertices.
    !
    !    Input, integer TRI_NUM, the number of triangles.
    !
    !    Input, integer TRI_VERT(3,TRI_NUM), the triangle incidence list.
    !
    !    Input, integer TRI_NABE(3,TRI_NUM), the triangle neighbor list; negative
    !    values are used for links of a counter clockwise linked list of boundary
    !    edges;
    !      LINK = -(3*I + J-1) where I, J = triangle, edge index.
    !
    !    Input/output, integer LTRI, LEDG.  If LTRI /= 0 then these values are
    !    assumed to be already computed and are not changed, else they are updated.
    !    On output, LTRI is the index of boundary triangle to the left of the
    !    leftmost boundary triangle visible from (X,Y), and LEDG is the boundary
    !    edge of triangle LTRI to the left of the leftmost boundary edge visible
    !    from (X,Y).  1 <= LEDG <= 3.
    !
    !    Input/output, integer RTRI.  On input, the index of the boundary triangle
    !    to begin the search at.  On output, the index of the rightmost boundary
    !    triangle visible from (X,Y).
    !
    !    Input/output, integer REDG, the edge of triangle RTRI that is visible
    !    from (X,Y).  1 <= REDG <= 3.
    !
    implicit none

    integer point_num
    integer tri_num

    integer a
    integer b
    integer e
    integer l
    logical ldone
    integer ledg
    integer lr
    integer ltri
    real  point_xy(2,point_num)
    integer redg
    integer rtri
    integer t
    integer tri_nabe(3,tri_num)
    integer tri_vert(3,tri_num)
    real  x
    real  y
    !
    !  Find the rightmost visible boundary edge using links, then possibly
    !  leftmost visible boundary edge using triangle neighbor information.
    !
    if ( ltri == 0 ) then
       ldone = .false.
       ltri = rtri
       ledg = redg
    else
       ldone = .true.
    end if

    do

       l = -tri_nabe(redg,rtri)
       t = l / 3
       e = mod ( l, 3 ) + 1
       a = tri_vert(e,t)

       if ( e <= 2 ) then
          b = tri_vert(e+1,t)
       else
          b = tri_vert(1,t)
       end if

       lr = lrline ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), &
            point_xy(2,b), 0.0 )

       if ( lr <= 0 ) then
          exit
       end if

       rtri = t
       redg = e

    end do

    if ( ldone ) then
       return
    end if

    t = ltri
    e = ledg

    do

       b = tri_vert(e,t)
       e = i4_wrap ( e-1, 1, 3 )

       do while ( 0 < tri_nabe(e,t) )

          t = tri_nabe(e,t)

          if ( tri_vert(1,t) == b ) then
             e = 3
          else if ( tri_vert(2,t) == b ) then
             e = 1
          else
             e = 2
          end if

       end do

       a = tri_vert(e,t)

       lr = lrline ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), &
            point_xy(2,b), 0.0 )

       if ( lr <= 0 ) then
          exit
       end if

    end do

    ltri = t
    ledg = e

  end subroutine vbedg

  !============================================================================
  subroutine test_triangulation

    character(len=*), parameter:: NameSub = 'test_triangulation'

    integer, parameter:: n1 = 5, n2 = 3, nPoint = n1*n2
    real   :: CoordXy_DI(2,n1*n2)
    integer:: iNodeTriangle_II(3,2*n1*n2)
    logical:: IsTriangleFound
    integer:: nTriangle, iNode1, iNode2, iNode3
    real   :: Area, Weight1, Weight2, Weight3
    integer:: i1, i2, n
    integer, parameter:: nX=100, nY=100
    integer:: nTriangle_C(nX,nY)
    integer, pointer, save:: iTriangle_IC(:,:,:)
    !------------------------------------------------------------------------
    
    ! Create a slanted mesh
    n = 0
    do i2 = 1, n2
       do i1 = 1, n1
          n = n + 1
          CoordXy_DI(:,n) = (/ 10.0*i1 + i2, 3.0*i2 /)
       end do
    end do

    write(*,'(a)')'Testing mesh_triangulation'
    call mesh_triangulation(n1, n2, CoordXy_DI, &
         iNodeTriangle_II, nTriangle)

    if(nTriangle /= 2*(n1-1)*(n2-1))write(*,*) NameSub, &
         ' error: nTriangle=', nTriangle, ' should be', 2*(n1-1)*(n2-1)

    if(any( iNodeTriangle_II(:,1) /= (/ 2, n1+1, 1 /))) write(*,*) NameSub, &
         ' error: iNodeTriangle_II(:,1) =', iNodeTriangle_II(:,1) , &
         ' should be ', 2, n1+1, 1

    if(any( iNodeTriangle_II(:,nTriangle) /= (/ n1*(n2-1), n1*n2, n1*n2-1 /)))&
         write(*,*) NameSub, &
         ' error: iNodeTriangle_II(:,1) =', iNodeTriangle_II(:,nTriangle) , &
         ' should be ', n1*(n2-1), n1*n2, n1*n2-1
 
    write(*,'(a)')'Testing calc_triangulation'
    
    call calc_triangulation(nPoint, CoordXy_DI, &
         iNodeTriangle_II, nTriangle)

    if(nTriangle /= 2*(n1-1)*(n2-1))write(*,*) NameSub, &
         ' error: nTriangle=', nTriangle, ' should be', 2*(n1-1)*(n2-1)

    if(any( iNodeTriangle_II(:,1) /= (/ n1+1, 1, 2 /))) write(*,*) NameSub, &
         ' error: iNodeTriangle_II(:,1) =', iNodeTriangle_II(:,1) , &
         ' should be ', n1+1, 1, 2

    if(any( iNodeTriangle_II(:,nTriangle) /= (/ n1*n2-1, n1*(n2-1), n1*n2 /)))&
         write(*,*) NameSub, &
         ' error: iNodeTriangle_II(:,1) =', iNodeTriangle_II(:,nTriangle) , &
         ' should be ', n1*n2-1, n1*(n2-1), n1*n2

    write(*,'(a)')'Testing triangle_area'

    Area = triangle_area(CoordXy_DI(:,iNodeTriangle_II(:,1)))
    if(abs(Area - 15.0) > 1e-6) write(*,*) NameSub, &
         ' error: Area = ', Area, ' should be 15'
 
    write(*,'(a)')'Testing find_triangle'
    ! test point outside of region with IsTriangleFound logical
    call find_triangle(&
         nPoint, nTriangle, (/10.9, 2.9/), CoordXy_DI, iNodeTriangle_II,&
         iNode1, iNode2, iNode3, IsTriangleFound=IsTriangleFound)

    if(IsTriangleFound) write(*,*) NameSub, &
         ' error: IsTriangleFound should be false'

    ! test point outside of region without IsTriangleFound logical
    call find_triangle(&
         nPoint, nTriangle, (/10.9, 2.9/), CoordXy_DI, iNodeTriangle_II,&
         iNode1, iNode2, iNode3, Weight1, Weight2, Weight3)

    if(any( (/iNode1, iNode2, iNode3/) /= (/ 1, 1, 1 /))) write(*,*) NameSub, &
          ' error: iNode1..3=', iNode1, iNode2, iNode3, &
          ' should be ', 1, 1, 1

    if(maxval( abs( (/ Weight1 - 1.0, Weight2, Weight3 /)) ) > 1e-6) &
         write(*,*) NameSub, &
         ' error: Weight1..3=, ', Weight1, Weight2, Weight3, &
         ' should be 1, 0, 0'
    
    ! test point inside the first triangle
    call find_triangle(&
         nPoint, nTriangle, (/11.1, 3.1/), CoordXy_DI, iNodeTriangle_II,&
         iNode1, iNode2, iNode3, Weight1, Weight2, Weight3, IsTriangleFound)

    if(.not. IsTriangleFound) write(*,*) NameSub, &
         ' error: IsTriangleFound 1st should be true'

    if(any( (/iNode1, iNode2, iNode3/) /= iNodeTriangle_II(:,1) )) &
          write(*,*) NameSub, 'error: iNode1..3=', iNode1, iNode2, iNode3, &
          ' should be ', iNodeTriangle_II(:,1)

    if( 1e-6 < maxval(abs( &
         (/Weight1-0.5/Area, Weight2-(Area-0.6)/Area, Weight3-0.1/Area/)))) &
         write(*,*) NameSub, &
         ' error: Weight1..3=, ', Weight1, Weight2, Weight3, &
         ' should be ', 0.5/Area, (Area-0.6)/Area, 0.1/Area

    nTriangle_C(1,1) = -1
    call find_triangle(                                                     &
         nPoint, nTriangle, (/11.1, 3.1/), CoordXy_DI, iNodeTriangle_II,    &
         iNode1, iNode2, iNode3, Weight1, Weight2, Weight3, IsTriangleFound,&
         nTriangle_C, iTriangle_IC)

    if(.not. IsTriangleFound) write(*,*) NameSub, &
         ' error: IsTriangleFound 2nd should be true'

    if(nTriangle_C(1,1) /= 1)write(*,*) NameSub, &
         'error nTriangle_C(1,1)=',nTriangle_C(1,1),' should be 1'

    if(iTriangle_IC(1,1,1) /= 1)write(*,*) NameSub, &
         'error iTriangle_IC(1,1,1)=',iTriangle_IC(1,1,1),' should be 1'

    if(nTriangle_C(nX,nY) /= 1)write(*,*) NameSub, &
         'error nTriangle_C(nX,nY)=',nTriangle_C(nX,nY),' should be 1'

    if(iTriangle_IC(1,nX,nY) /= nTriangle)write(*,*) NameSub, &
         'error iTriangle_IC(1,nX,nY)=',iTriangle_IC(1,nX,nY),&
         ' should be ', nTriangle

    if(any( (/iNode1, iNode2, iNode3/) /= iNodeTriangle_II(:,1) )) &
          write(*,*) NameSub, 'error: iNode1..3=', iNode1, iNode2, iNode3, &
          ' should be ', iNodeTriangle_II(:,1)

    if( 1e-6 < maxval(abs( &
         (/Weight1-0.5/Area, Weight2-(Area-0.6)/Area, Weight3-0.1/Area/)))) &
         write(*,*) NameSub, &
         ' error: Weight1..3=, ', Weight1, Weight2, Weight3, &
         ' should be ', 0.5/Area, (Area-0.6)/Area, 0.1/Area

    call find_triangle(                                                     &
         nPoint, nTriangle, (/53.0, 9.0/), CoordXy_DI, iNodeTriangle_II,    &
         iNode1, iNode2, iNode3, Weight1, Weight2, Weight3, IsTriangleFound,&
         nTriangle_C, iTriangle_IC)

    if(.not. IsTriangleFound) write(*,*) NameSub, &
         ' error: IsTriangleFound 3rd should be true'

    if(any( (/iNode1, iNode2, iNode3/) /= iNodeTriangle_II(:,nTriangle) )) &
          write(*,*) NameSub, 'error: iNode1..3=', iNode1, iNode2, iNode3, &
          ' should be ', iNodeTriangle_II(:,nTriangle)

    if( 1e-6 < maxval(abs( &
         (/Weight1, Weight2, Weight3-1.0/)))) &
         write(*,*) NameSub, &
         ' error: Weight1..3=, ', Weight1, Weight2, Weight3, &
         ' should be 0, 0, 1'


  end subroutine test_triangulation

end module ModTriangulate
