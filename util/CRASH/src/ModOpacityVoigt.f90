!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!A part of the package ionmix
module CRASH_ModOpacityVoigt
  implicit none
  SAVE
  PRIVATE !Except

  real,parameter:: aTab_I(11)=(/  1.00E-01,1.58E-01,2.51E-01,3.98E-01,6.31E-01,&       
       1.00E+00,1.58E+00,2.51E+00,3.98E+00,6.31E+00,1.00E+01 /) 

  real,parameter :: vTab_I(11)=(/  0.00E+00,5.00E-01,1.00E+00,1.50E+00,2.00E+00, &       
       2.50E+00,3.00E+00,3.50E+00,4.00E+00,4.50E+00,5.00E+00 /)  

  real,parameter:: hTable_II(11,11)=reshape( (/&
       8.96E-01,8.44E-01,7.69E-01,6.72E-01,5.54E-01,4.28E-01,3.08E-01, &  
       2.10E-01,1.38E-01,8.83E-02,5.61E-02,7.18E-01,6.85E-01,6.38E-01, & 
       5.72E-01,4.89E-01,3.91E-01,2.92E-01,2.04E-01,1.36E-01,8.78E-02, &  
       5.60E-02,3.73E-01,3.74E-01,3.72E-01,3.63E-01,3.43E-01,3.05E-01, & 
       2.50E-01,1.88E-01,1.30E-01,8.63E-02,5.56E-02,1.34E-01,1.48E-01, & 
       1.66E-01,1.87E-01,2.05E-01,2.12E-01,1.98E-01,1.65E-01,1.22E-01, &
       8.39E-02,5.49E-02,4.02E-02,5.18E-02,6.85E-02,9.07E-02,1.17E-01, &
       1.40E-01,1.51E-01,1.40E-01,1.12E-01,8.07E-02,5.40E-02,1.47E-02, &
       2.19E-02,3.28E-02,4.86E-02,6.97E-02,9.38E-02,1.13E-01,1.17E-01, &
       1.02E-01,7.69E-02,5.29E-02,7.94E-03,1.25E-02,1.95E-02,3.01E-02, &
       4.55E-02,6.53E-02,8.53E-02,9.64E-02,9.11E-02,7.27E-02,5.16E-02, &
       5.34E-03,8.44E-03,1.33E-02,2.08E-02,3.21E-02,4.77E-02,6.58E-02, &
       7.97E-02,8.09E-02,6.83E-02,5.01E-02,3.92E-03,6.21E-03,9.81E-03, &
       1.54E-02,2.40E-02,3.63E-02,5.18E-02,6.62E-02,7.16E-02,6.39E-02, &
       4.85E-02,3.02E-03,4.79E-03,7.57E-03,1.19E-02,1.86E-02,2.85E-02, &
       4.17E-02,5.55E-02,6.33E-02,5.94E-02,4.69E-02,2.41E-03,3.81E-03, &
       6.03E-03,9.52E-03,1.49E-02,2.30E-02,3.42E-02,4.69E-02,5.59E-02, &
       5.51E-02,4.51E-02 /),(/11,11/))

  !PUBLIC MEMBERS:

  !Locigal to determine, if we use the Voigt profile.
  !Otherwise, only the Lorenz broadenning is accounted for
  !with the Gaussin profile for the line.

  logical,public:: UseVoigt = .true.

 
  
  ! ... this subroutine computes the Voigt function "h" as a function     
  !     of "a" and "v".  See Mihalas, "Stellar Atmospheres", 1978, for    
  !     definitions and details.
                                          
  public :: voigt_profile 

  ! ... compute the total line width of a bound-bound transition.  this   
  !     sums the contributions from natural, Doppler, and collisional     
  !     broadening.   
  public :: line_width

contains
  !================================================
  subroutine line_width( tp,densnn,ennp,AtomicWeight, &      
                           gamma,avoigt,dnudop )
    real,intent(in) :: tp,densnn,ennp
    real,intent(in):: AtomicWeight
    real,intent(out):: gamma,avoigt,dnudop 
    ! ... compute the total line width of a bound-bound transition.  this   
    !     sums the contributions from natural, Doppler, and collisional     
    !     broadening.                                                       
    !                                                                       
    ! ... input variables:                                                  
    !       tp      -  plasma temperature (eV)                              
    !       densnn  -  number density of all nuclei (cm^{-3})                
    !       ennp    -  energy of the transition from "n" to "np" (eV)       
    !       AtomicWeight  -  atomic weight of the ion (amu)                       
    !                                                                       
    ! ... output variables:                                                 
    !       gamma   -  total line width (sec**-1)                           
    !       avoigt  -  "a" used to compute the Voigt profile                
    !       dnudop  -  doppler width (sec**-1)                              
    !                                                                       
    !Local variables:
    real :: Vel,    &!Mean squared velocity
            WidNat, &!Natural width
            WidDop, &!Doppler width
            WidCol   !Collisional width
    
    vel = sqrt( tp / AtomicWeight )                                         
    
    widnat = 2.29e6 * ennp**2                                         
    widdop = 1.41e11 * ennp * vel                                     
    widcol = 4.58e6 * densnn**0.3333 * vel                           
                                                                        
    gamma = widnat + widdop + widcol  
    if(widdop<1.e-30)then
       write(*,*)ennp,densnn,tp,vel
       call CON_stop('Zero Doppler width of the line')
    end if
    avoigt = ( widnat + widcol ) / widdop                             
    dnudop = widdop / 12.566  
  end subroutine line_width
  !====================================================
  real function voigt_profile ( a,vv )                                           
    real,intent(in)::a,vv       
    ! ... this subroutine computes the Voigt function "h" as a function 
    !     of "a" and "v".  See Mihalas, "Stellar Atmospheres", 1978, for
    !     definitions and details.                                          


    real,parameter:: cSqrtPi = 1.7725                                                                                                
    real :: v,LogA, fa, fv, g1, g2
    integer:: iA, iV
    !---------------------------

    v = abs( vv )                                                     

    if ( v.gt.5. .or. a.ge.10. ) then                                 

       ! ...    Lorentzian profile in limit of large a,v                       
       voigt_profile = a / cSqrtPi / ( a*a + v*v )                             

    else if ( a.le.0.1 ) then                                         

       ! ...    use Dawson's integral (see Mihalas) for small a     
       voigt_profile = exp( -min(v*v,50.0) ) + 2.0*a/cSqrtPi * &                 
            ( 2.0*v*dawson(v) - 1.0 )                                

    else                                                              

       ! ...    use a table lookup                                   
       LogA = log10( a )                                              
       ia = int(LogA/0.2 + 6)                                              
       iv = int(v/0.5 + 1)                                                 
       fa = LogA/0.2 - (ia-6)                                         
       fv = v/0.5 - (iv-1)                                            
       g1 = (1.-fa)*hTable_II(ia,iv) + fa*hTable_II(ia+1,iv)                
       g2 = (1.-fa)*hTable_II(ia,iv+1) + fa*hTable_II(ia+1,iv+1)            
       voigt_profile = (1.-fv)*g1 + fv*g2                                     

    endif
  contains
   real function dawson ( v )                                             
     real,intent(in) :: v
           
     ! ... computes dawson's integral (which is sometimes used to compute    
     !     the Voigt function):                                              
     !                               _v                                      
     !       D.I.  =  exp(-v**2) * _/ dt exp(t**2)                           
     !                            0                                         
     !                                                                       
                                       
                                                                        
      real,parameter:: Table1_I(21)=(/ .00000, .09933, .19475, .28263, .35994, &   
            .42443, .47476, .51050, .53210, .54072, .53807, &           
            .52620, .50727, .48339, .45650, .42824, .39993, &           
            .37255, .34677, .32297, .30134 /)                          
      real,parameter:: Table2_I(11) = &
           (/ .50000, .50650, .51358, .52142, .53037, &           
              .54079, .55265, .56547, .57852, .59110, .60268 /)  
      real:: f,y,x
      integer::iX
      !------------------------------------------------------

      if ( v .le. 2. ) then                                             
         x = v                                                          
         iX = 1 + x / 0.1                                               
         f = x/0.1 - (iX-1)                                             
         y = (1.-f) * Table1_I(iX) + f * Table1_I(iX+1)                     
         dawson = y                                                     
      else                                                              
         x = 1. / (v*v)                                                 
         ix = 1 + x / 0.025                                             
         f = x/0.025 - (iX-1)                                           
         y = (1.-f) * Table2_I(iX) + f * Table2_I(iX+1)                     
         dawson = y / v                                                 
      endif                                                             
                                                                                                                          
    end function dawson

  end function voigt_profile
end module CRASH_ModOpacityVoigt
