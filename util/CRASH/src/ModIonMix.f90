!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CRASH_ModIonMix
  !This is an interface to some methods of the code iomix.f
  ! IONMIX, A CODE FOR COMPUTING THE EQUATION OF STATE AND   
  ! RADIATIVE PROPERTIES OF LTE AND NON-LTE PLASMAS.  J.J. MACFARLANE.  
  ! REF. IN COMP. PHYS. COMMUN. 56 (1989) 259 
  use ModNumConst
  implicit none
  SAVE
  real,parameter:: con(9)=(/&
       0.0,&
       1.0e-10,& !min. species concentration to compute bb and bf transitions
       1.0e-10,& !min. ionization concentration to compute bb and bf transitions
       1.0e-10,& !min. concentration of an atomic state 4 bb and bf transitions
       1.0e10, & !range, in # of line widths (fwhm)to compute contribution from bb
       10.0,   & !width, in # of line widths (fwhm),of line core           
       10.0,   & ! abs. cfs. are weighted by Planck fn. when 1/con7 < hv/kT < con7
       0.0,    &
       1.0     /)! multiplier for mesh pt spacing about lines

  !\
  ! Normalization parameters
  !/                    \approx  1/ 6.4939     \approx 1/ 25.976 
  real,parameter :: cNormG5 = 15.0/(cPi**4), cNormG6 = 0.250 * cNormG5

  !common variables used by {\it gint}s
  real,parameter :: xMin =  .005 
  real ::  xmnlog =  -5.298317366548036  !This is log(xMin)

  !This is log(30/xMin)/(50-1), where 30 was the original value of xMax and 50 was the original number 
  real:: DxLog = .177541117310412
  
  !This is 30/exp(4*DxLog), where 4 is the number of truncated points
  real, parameter :: xMax =  14.74690031285409   
  
  !Misc
  real:: xLog, XII, EX
  integer :: II 
  

contains
 
  !===============================================
  real function gint ( igint,xStart,xEnd )
    real,intent(in) :: xStart, xEnd
    integer,intent(in) :: iGInt
    !                                                  
    ! ... this function returns the integral <from xmin to xmax> of  
    !     GINTn (where n = 5 thru 6)                                                                                    

   
    !----------------                                               

    select case(igint)                             
    case(5)
       if(xStart > xMax) then
          gint = gint2(xEnd) - gint2(xStart)
       else
          gint = gint5(xEnd) - gint5(xStart)
       end if
    case(6)   
       if(xStart > xMax) then
          gint = gint1(xEnd) - gint1(xStart)
       else
          gint = gint6(xEnd) - gint6(xStart)   
       end if
    case default
       call CON_stop('Error in g_int')                                                        
    end select
  end function gint

  !==================
  real function gint1 ( x )                                         
    real,intent(in) :: x   
    
    ! ... this routine returns the integral (from infinity to x) of the            
    !     function:                                                         
    !                                                                       
    !                     x_                                                
    !           gint1  = _/  dy y**4 exp(-y)                                
    !                   \infty                                                   
    !                                                                       
    gint1 = -exp( -x ) * (x**4 + 4.0 * x**3 + 12.0 * x**2 + 24.0 * x + 24.0 )               
  end function gint1
  
  !==================
  real function gint2 ( x )                                         
    real,intent(in) :: x   
    
    ! ... this routine returns the integral (from infinity to x) of the            
    !     function:                                                         
    !                                                                       
    !                     x_                                                
    !           gint2  = _/  dy y**3 exp(-y)                                
    !                   \infty                                                   
    !                                                                       
    gint2 = -exp( -x ) * ( x**3 + 3.*x**2 + 6.*x + 6. )               
  end function gint2
  
  !==================
  real  function gint5 ( x )                                          
    real, intent(in):: x                                                 
    ! ... this routine returns the integral (from 0 to x) of the           
    !     function:                                                         
    !                                                                       
    !                     x_                                                
    !           gint5  = _/  dy y**3 / ( exp(y)-1 )                         
    !                   0                                                   
    !                                                                       
    ! ... the table values are the logarithms of the integrals for          
    !     equally spaced values of  log(x).                                 
    !                                                                                                                  
    real,parameter:: gtable(46) = (/&                                                     
         -1.6995E+01,-1.6463E+01,-1.5931E+01,-1.5399E+01,-1.4867E+01,&    
         -1.4335E+01,-1.3803E+01,-1.3272E+01,-1.2740E+01,-1.2209E+01,&     
         -1.1678E+01,-1.1148E+01,-1.0618E+01,-1.0088E+01,-9.5594E+00,&     
         -9.0312E+00,-8.5039E+00,-7.9775E+00,-7.4524E+00,-6.9289E+00,&     
         -6.4070E+00,-5.8874E+00,-5.3703E+00,-4.8563E+00,-4.3460E+00,&     
         -3.8402E+00,-3.3399E+00,-2.8462E+00,-2.3605E+00,-1.8846E+00,&     
         -1.4204E+00,-9.7072E-01,-5.3858E-01,-1.2780E-01, 2.5710E-01,&     
         6.1092E-01, 9.2788E-01, 1.2021E+00, 1.4283E+00, 1.6031E+00,&     
         1.7268E+00, 1.8043E+00, 1.8456E+00, 1.8634E+00, 1.8693E+00,&     
         1.8706E+00/) ! truncated: give zero integrals: 1.8708E+00, 1.8709E+00, 1.8709E+00, 1.8709E+00 
    !========================                                                                  
    if ( x .lt. xmin ) then                                           
       
       ! ...    analytic form using expansion of exp(x)                        
       gint5 = x**3 / 3.                                             
       
    else if ( x .gt. xmax*0.9999 ) then                               
       
       ! ...    the integral as x -> infinity + gint2(x)                                  
       gint5 = 1.0/cNormG5 + gint2(x)                                      
       
    else                                                              
       
       ! ...    otherwise, interpolate using table                             
       xlog = log( x )                                                
       ii = 1 + ( xlog-xmnlog ) / dxlog                              
       xii = xmnlog + ( ii-1 )*dxlog                                 
       ex = gtable(ii) + ( xlog-xii ) * &                             
            ( gtable(ii+1) - gtable(ii) ) / dxlog              
       gint5 = exp( ex )                                             
       
    endif
    
  end function gint5
  !-----------------------------------------------------------------
  real function gint6 ( x )
    real,intent(in) :: x                                             
    
    ! ... this routine returns the integral (from 0 to x) of the            
    !     function:                                                         
    !                                                                       
    !                     x_                                                
    !           gint6  = _/  dy y**4 exp(-y) / ( 1-exp(-y) )**2             
    !                   0                                                   
    !                                                                       
    ! ... the table values are the logarithms of the integrals for          
    !     equally spaced values of  log(x).                                 
    !                                                                       
    
    
    real,parameter:: gtable(46)= (/   &                                                 
         -1.6994E+01,-1.6461E+01,-1.5928E+01,-1.5396E+01,-1.4863E+01,&
         -1.4330E+01,-1.3798E+01,-1.3265E+01,-1.2733E+01,-1.2200E+01,&
         -1.1667E+01,-1.1135E+01,-1.0602E+01,-1.0070E+01,-9.5370E+00,&
         -9.0045E+00,-8.4720E+00,-7.9395E+00,-7.4071E+00,-6.8748E+00,&
         -6.3426E+00,-5.8106E+00,-5.2789E+00,-4.7476E+00,-4.2169E+00,&
         -3.6869E+00,-3.1581E+00,-2.6309E+00,-2.1060E+00,-1.5843E+00,&
         -1.0671E+00,-5.5646E-01,-5.4768E-02, 4.3443E-01, 9.0650E-01,&
         1.3554E+00, 1.7736E+00, 2.1524E+00, 2.4822E+00, 2.7542E+00,&
         2.9625E+00, 3.1063E+00, 3.1926E+00, 3.2353E+00, 3.2516E+00,&
         3.2562E+00/) ! truncated: give zero integrals: 3.2571E+00, 3.2572E+00, 3.2572E+00, 3.2572E+00 /)    
    !---------------------------
    
    if ( x .lt. xmin ) then                     
       ! ...    analytic form using expansion of exp(x)               
       gint6 = x**3 / 3.                                              
       
    else if ( x .gt. xmax*0.9999 ) then                               
       
       ! ...    the integral as x -> infinity                                  
       gint6 = 1.0/cNormG6 + gint1(x)                               
       
    else                                                              
       
       ! ...    otherwise, interpolate using table                             
       xlog = log( x )                                                
       ii = 1 + ( xlog-xmnlog ) / dxlog                               
       xii = xmnlog + ( ii-1 )*dxlog                                  
       ex = gtable(ii) + ( xlog-xii ) * &                             
            ( gtable(ii+1) - gtable(ii) ) / dxlog               
       gint6 = exp( ex )
    end if
  end function gint6
end module CRASH_ModIonMix
