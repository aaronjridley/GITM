!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CRASH_ModExcitationData
  use CRASH_ModAtomicMass,ONLY : nZMax
  use CRASH_ModIonization,ONLY : get_ioniz_potential
  implicit none
  PRIVATE !except
  !\
  !  Maximum principal quantum number for the excitation states.
  !/
  integer,parameter,public:: nExcitation = 10

  !public methods
  !Principal quantum number as a function of the charge state and the element number
  public:: n_ground  

  !The number of electrons at inner shells
  public:: n_screened

  public:: init_excitation

  !For test of formula for excitation energies
  public :: cExcitationN_III,cExcitationBe_III

  real,   public, allocatable:: &
       cExcitationEnergy_IIII(:,:,:,:), cVirialCoeff4Energy_IIII(:,:,:,:)

  integer, public, allocatable:: cDegeneracy_IIII(:,:,:,:)

  !\
  ! The logical to handle whether excitation levels should
  ! be accounted for
  !/
  logical,public :: UseExcitation = .false.

  logical,public :: UsePressureIonization = .false. ! UseExcitation should be also .true.

  !Determine if we account for the optical transitions with no change 
  !in the principal quantum number

  logical,parameter,public :: UseDeltaNEq0Transition = .false.

  !Determine if we account for bound-free transition, in which the
  !bound electron is photo-ionized from inner shells ("core electron")

  logical,public,parameter :: UseCoreElectron = .false.


  !+++++++++++++++++++++++INTERNAL DATABASE==================================
  ! Tables of excitation energies
  ! He, nZ = 2
  ! From SPECTR-w3
  !          S     |      P     |      D     |      F     |      G     |
  real, parameter, dimension(0:4, 2:5, 0:1) :: cExcitationHe_III = reshape((/&
                                !Ionization potential = 24.58741, nGround = 1
       19.820112,   20.964486,         0.0,         0.0,         0.0,&  !n = 2
       22.718862,   23.006505,   23.070977,         0.0,         0.0,&  !n = 3
       23.594191,   23.708256,   23.735533,         0.0,         0.0,&  !n = 4
       23.972342,   24.026895,   24.038054,         0.0,         0.0,&  !n = 5
                                !Ionization potential = 54.41778, nGround = 1
       40.813085,   40.813027,         0.0,         0.0,         0.0,&  !n = 2
       48.371310,   48.371292,   48.371507,         0.0,         0.0,&  !n = 3
       51.016663,   51.016656,   51.016747,   51.016777,         0.0,&  !n = 4
       52.241073,   52.241069,   52.241115,   52.241131,   52.241138 &  !n = 5
       /), (/5, 4, 2/))

  ! Be, nZ = 4
  ! From SPECTR-w3
  !          S     |      P     |      D     |      F     |      G     |
  real, parameter, dimension(0:4, 2:5, 0:3) :: cExcitationBe_III = reshape((/&
                                !Ionization potential = 9.32263, nGround = 2
       0.0,    2.724652,         0.0,         0.0,         0.0,&  !n = 2
       6.455745,    7.288175,    7.693182,         0.0,         0.0,&  !n = 3
       7.949928,    8.220561,    8.366949,         0.0,         0.0,&  !n = 4
       8.484647,    8.611297,    8.680418,         0.0,         0.0,&  !n = 5
                                !Ionization potential = 18.21116, nGround = 2
       0.0,    3.956807,         0.0,         0.0,         0.0,&  !n = 2
       10.938703,   11.963904,   12.156649,         0.0,         0.0,&  !n = 3
       14.315214,   14.724362,   14.806191,         0.0,         0.0,&  !n = 4
       15.787526,   15.990240,   16.031031,         0.0,         0.0,&  !n = 5
                                !Ionization potential = 153.897, nGround = 1
       118.590874,  121.923569,         0.0,         0.0,         0.0,&  !n = 2
       139.008589,  139.891357,  140.273228,         0.0,         0.0,&  !n = 3
       145.718613,  146.074448,  146.233148,         0.0,         0.0,&  !n = 4
       148.731429,  148.912446,  148.991796,  148.996755,         0.0,&  !n = 5
                                !Ionization potential = 217.713, nGround = 1
       163.285338,  163.284594,         0.0,         0.0,         0.0,&  !n = 2
       193.527201,  193.526977,  193.530424,         0.0,         0.0,&  !n = 3
       204.111433,  204.111346,  204.112797,  204.113280,         0.0,&  !n = 4
       209.010234,  209.010185,  209.010928,  209.011176,  209.011300 &  !n = 5
       /), (/5, 4, 4/))

  ! C, nZ = 6
  ! From SPECTR-w3 and
  ! http://physics.nist.gov/PhysRefData/Handbook/Tables/carbontable5.htm
  !          S     |      P     |      D     |      F     |      G     |
  real, parameter, dimension(0:4, 2:5, 0:5) :: cExcitationC_III = reshape((/&
                                !Ionization potential = 11.26030, nGround = 2
       0.0,    0.002480,         0.0,         0.0,         0.0,&  !n = 2
       7.479966,    8.537097,    9.631092,         0.0,         0.0,&  !n = 3
       9.683165,         0.0,   10.352680,         0.0,         0.0,&  !n = 4
       10.382436,         0.0,   10.686197,         0.0,         0.0,&  !n = 5
                                !Ionization potential = 24.38332, nGround = 2
       0.0,    0.007439,         0.0,         0.0,         0.0,&  !n = 2
       14.448869,   16.331197,   18.045898,         0.0,         0.0,&  !n = 3
       19.494034,   20.149786,   20.844221,   20.950848,         0.0,&  !n = 4
       21.492411,   21.732940,   22.130557,         0.0,         0.0,&  !n = 5
                                !Ionization potential = 47.888, nGround = 2
       0.0,    6.488266,         0.0,         0.0,         0.0,&  !n = 2
       29.538612,   32.096530,   33.476226,         0.0,         0.0,&  !n = 3
       38.409805,   39.376633,   39.849757,   39.922908,         0.0,&  !n = 4
       41.969887,   42.494836,   42.734249,   43.042350,         0.0,&  !n = 5
                                !Ionization potential = 64.492, nGround = 2
       0.0,    7.994252,         0.0,         0.0,         0.0,&  !n = 2
       37.548611,   39.679981,   40.189981,         0.0,         0.0,&  !n = 3
       49.736504,   50.623983,   50.875175,   50.886829,         0.0,&  !n = 4
       55.200363,   55.651294,   55.779245,   55.785445,         0.0,&  !n = 5
                                !Ionization potential = 392.08, nGround = 1
       298.958108,  304.399774,         0.0,         0.0,         0.0,&  !n = 2
       352.060534,  353.530987,  354.260014,         0.0,         0.0,&  !n = 3
       369.911777,  370.509381,  370.809423,  370.826781,         0.0,&  !n = 4
       378.019103,  378.316665,  378.472885,  378.480324,  378.481564,&  !n = 5
                                !Ionization potential = 489.98, nGround = 1
       367.477261,  367.474025,         0.0,         0.0,         0.0,&  !n = 2
       435.547716,  435.546749,  435.564169,         0.0,         0.0,&  !n = 3
       459.370211,  459.369802,  459.377154,  459.379596,         0.0,&  !n = 4
       470.395666,  470.395455,  470.399212,  470.400464,  470.401096 &  !n = 5
       /), (/5, 4, 6/))

  ! N, nZ = 7
  ! From SPECTR-w3
  !          S     |      P     |      D     |      F     |      G     |
  real, parameter, dimension(0:4, 2:5, 0:6) :: cExcitationN_III = reshape((/&
                                !Ionization potential = 14.53414, nGround = 2
       0.0,    2.382976,         0.0,         0.0,         0.0,&  !n = 2
       10.325403,   11.844209,   12.971226,         0.0,         0.0,&  !n = 3
       12.847241,         0.0,   13.661817,         0.0,         0.0,&  !n = 4
       13.614703,         0.0,   13.980457,         0.0,         0.0,&  !n = 5
                                !Ionization potential = 29.6013, nGround = 2
       0.0,    0.006199,         0.0,         0.0,         0.0,&  !n = 2
       18.462485,   20.409037,   23.196201,         0.0,         0.0,&  !n = 3
       24.367852,   25.065883,   26.028000,         0.0,         0.0,&  !n = 4
       26.558652,         0.0,   27.338513,         0.0,         0.0,&  !n = 5
                                !Ionization potential = 47.449, nGround = 2
       0.0,    0.021077,         0.0,         0.0,         0.0,&  !n = 2
       27.437700,   30.458699,   33.133410,         0.0,         0.0,&  !n = 3
       37.329159,   38.644755,   39.395975,   39.710895,         0.0,&  !n = 4
       41.374763,   41.898472,   42.395896,   42.495580,         0.0,&  !n = 5
                                !Ionization potential = 77.472, nGround = 2
       0.0,    8.338842,         0.0,         0.0,         0.0,&  !n = 2
       46.776754,   50.142428,   52.079557,         0.0,         0.0,&  !n = 3
       60.455929,   62.448355,   63.343149,   64.045271,         0.0,&  !n = 4
       67.437230,   68.081080,   68.411498,   68.811223,         0.0,&  !n = 5
                                !Ionization potential = 97.89, nGround = 2
       0.0,    9.975768,         0.0,         0.0,         0.0,&  !n = 2
       56.552287,   59.235304,   60.057940,         0.0,         0.0,&  !n = 3
       75.149667,   76.267632,   76.611440,   76.629542,         0.0,&  !n = 4
       83.530874,   84.098473,   84.273291,   84.283210,   84.284449,&  !n = 5
                                !Ionization potential = 552.06, nGround = 1
       419.796815,  426.294826,         0.0,         0.0,         0.0,&  !n = 2
       494.927512,  496.684367,  497.600611,         0.0,         0.0,&  !n = 3
       520.336831,  521.052219,  521.432851,  521.453928,         0.0,&  !n = 4
       531.913234,  532.272788,  532.458764,  532.474882,         0.0,&  !n = 5
                                !Ionization potential = 667.03, nGround = 1
       500.252260,  500.246681,         0.0,         0.0,         0.0,&  !n = 2
       592.926844,  592.925108,  592.957344,         0.0,         0.0,&  !n = 3
       625.358875,  625.358131,  625.371769,  625.376233,         0.0,&  !n = 4
       640.368400,  640.368028,  640.375096,  640.377327,  640.378567 &  !n = 5
       /), (/5, 4, 7/))

  ! O, nZ = 8
  ! Partly from SPECTR-w3, partly from
  ! http://physics.nist.gov/PhysRefData/Handbook/Tables/oxygentable5.htm and
  ! http://physics.nist.gov/PhysRefData/Handbook/Tables/oxygentable6.htm
  !          S     |      P     |      D     |      F     |      G     |
  real, parameter, dimension(0:4, 2:5, 0:7) :: cExcitationO_III = reshape((/&
                                !Ionization potential = 13.61806, nGround = 2
       0.0,    0.019837,         0.0,         0.0,         0.0,&  !n = 2
       9.146313,   10.740224,   12.078539,         0.0,         0.0,&  !n = 3
       11.838010,         0.0,   12.754253,         0.0,         0.0,&  !n = 4
       12.661265,         0.0,   13.069173,         0.0,         0.0,&  !n = 5
                                !Ionization potential = 35.11730, nGround = 2
       0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 2
       22.966245,   25.285619,   28.677062,         0.0,         0.0,&  !n = 3
       29.586031,         0.0,         0.0,         0.0,         0.0,&  !n = 4
       0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 5
                                !Ionization potential = 54.936, nGround = 2
       0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 2
       0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 3
       0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 4
       0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 5
                                !Ionization potential = 77.413, nGround = 2
       0.0,    0.047920,         0.0,         0.0,         0.0,&  !n = 2
       44.337985,   48.373670,   52.015085,         0.0,         0.0,&  !n = 3
       60.233997,   61.927001,   63.301366,   63.627444,         0.0,&  !n = 4
       66.873102,   67.775459,   68.442990,   68.501263,         0.0,&  !n = 5
                                !Ionization potential = 113.90, nGround = 2
       0.0,   10.174961,         0.0,         0.0,         0.0,&  !n = 2
       67.818358,   72.004808,   74.498006,         0.0,         0.0,&  !n = 3
       89.425082,   91.082006,   91.951755,         0.0,         0.0,&  !n = 4
       98.605987,   99.423290,   99.847936,         0.0,         0.0,&  !n = 5
                                !Ionization potential = 138.12, nGround = 2
       0.0,   11.948852,         0.0,         0.0,         0.0,&  !n = 2
       79.353970,   82.587230,   83.642955,         0.0,         0.0,&  !n = 3
       105.688583,  107.038647,  107.478667,  107.503092,         0.0,&  !n = 4
       117.601604,  118.287732,  118.518963,         0.0,         0.0,&  !n = 5
                                !Ionization potential = 739.29, nGround = 1
       560.983806,  568.544362,         0.0,         0.0,         0.0,&  !n = 2
       661.929250,  664.018384,  665.104485,         0.0,         0.0,&  !n = 3
       696.337342,  697.113483,  697.597021,  697.561065,         0.0,&  !n = 4
       711.991585,  712.343700,  712.613985,  712.561912,         0.0,&  !n = 5
                                !Ionization potential = 871.41, nGround = 1
       653.502789,  653.493738,         0.0,         0.0,         0.0,&  !n = 2
       774.581653,  774.578926,  774.633975,         0.0,         0.0,&  !n = 3
       816.952257,  816.951141,  816.974326,  816.982137,         0.0,&  !n = 4
       836.560728,  836.560232,  836.572134,  836.575978,  836.577962 &  !n = 5
       /), (/5, 4, 8/))

  ! Al, nZ = 13
  ! From SPECTR-w3 and W.C.Martin and Romuald Zalubas,
  ! Energy Levels of Aluminium, Al I through Al XIII

  !Ionization potential = 5.98577, nGround = 3
  !          S     |      P     |      D     |      F     |      G     |
  real, parameter, dimension(20) :: cExcitationAl0_I = (/&
       0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 2
       0.0,     3.59807,     4.02148,         0.0,         0.0,&  !n = 3
       3.14272,     4.08525,     4.82663,     5.12295,         0.0,&  !n = 4
       4.67289,     4.99309,     5.23631,     5.43436,         0.0 &  !n = 5
       /)

  !Ionization potential = 18.82856, nGround = 3
  real, parameter, dimension(20) :: cExcitationAl1_I = (/&
       0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 2
       0.0,     4.63614,    11.84662,         0.0,         0.0,&  !n = 3
       11.31659,    13.07135,    15.06203,    15.30194,         0.0,&  !n = 4
       14.88963,    15.58520,    16.46793,    16.54417,    16.63666 &  !n = 5
       /)

  !Ionization potential = 28.448, nGround = 3
  real, parameter, dimension(20) :: cExcitationAl2_I = (/&
       0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 2
       71.757087,    6.552254,   14.149199,         0.0,         0.0,&  !n = 3
       15.64235 ,   17.80827 ,   20.55488 ,   20.78133 ,         0.0,&  !n = 4
       21.15633 ,   22.12292 ,   23.41783 ,   23.54163 ,    23.54813 &  !n = 5
       /)

  !Ionization potential = 119.99, nGround = 2
  real, parameter, dimension(20) :: cExcitationAl3_I = (/&
       0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 2
       75.475869,   83.27181 ,   92.820265,         0.0,         0.0,&  !n = 3
       99.42072 ,  101.84171 ,  105.60007 ,  106.30523 ,         0.0,&  !n = 4
       107.98922 ,         0.0,  110.91749 ,  111.23825 ,  111.27212  &  !n = 5
       /)

  !Ionization potential = 153.83, nGround = 2
  real, parameter, dimension(20) :: cExcitationAl4_I = (/&
       0.0,    0.422492,         0.0,         0.0,         0.0,&  !n = 2
       91.969981,  100.240842,  111.727853,         0.0,         0.0,&  !n = 3
       124.119328,         0.0,  131.734437,  132.418830,         0.0,&  !n = 4
       0.0,         0.0,  139.798369,  140.131886,  140.203797 &  !n = 5
       /)

  !Ionization potential = 190.49, nGround = 2
  real, parameter, dimension(20) :: cExcitationAl5_I = (/&
       0.0,    0.319992,         0.0,         0.0,         0.0,&  !n = 2
       109.492418,  118.861779,  130.740084,         0.0,         0.0,&  !n = 3
       151.048694,         0.0,  159.066751,         0.0,         0.0,&  !n = 4
       174.168025,         0.0,  170.509251,         0.0,         0.0 &  !n = 5
       /)

  !Ionization potential = 241.76, nGround = 2
  real, parameter, dimension(20) :: cExcitationAl6_I = (/&
       0.0,   34.74037 ,         0.0,         0.0,         0.0,&  !n = 2
       142.23466 ,         0.0,  163.31445 ,         0.0,         0.0,&  !n = 3
       0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 4
       0.0,         0.0,         0.0,         0.0,         0.0 &  !n = 5
       /)

  !Ionization potential = 284.66, nGround = 2
  real, parameter, dimension(20) :: cExcitationAl7_I = (/&
       0.0,   16.55313 ,         0.0,         0.0,         0.0,&  !n = 2
       163.57730 ,  173.84567 ,  182.09557 ,         0.0,         0.0,&  !n = 3
       0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 4
       0.0,         0.0,         0.0,         0.0,         0.0 &  !n = 5
       /)

  !Ionization potential = 330.1, nGround = 2
  real, parameter, dimension(20) :: cExcitationAl8_I = (/&
       0.0,    0.606283,         0.0,         0.0,         0.0,&  !n = 2
       186.102742,  194.920498,  203.599391,         0.0,         0.0,&  !n = 3
       253.411277,  256.195962,  259.625365,         0.0,         0.0,&  !n = 4
       281.842091,  283.508438,  285.300010,         0.0,         0.0 &  !n = 5
       /)

  !Ionization potential = 398.8, nGround = 2
  real, parameter, dimension(20) :: cExcitationAl9_I = (/&
       0.0,   19.239742,         0.0,         0.0,         0.0,&  !n = 2
       230.011742,  238.659639,  243.827300,         0.0,         0.0,&  !n = 3
       306.910453,  310.267945,  312.126468,  313.49897 ,         0.0,&  !n = 4
       340.782933,  342.466638,  343.390320,         0.0,         0.0 &  !n = 5
       /)

  !Ionization potential = 442.0, nGround = 2
  real, parameter, dimension(20) :: cExcitationAl10_I = (/&
       0.0,   21.821713,         0.0,         0.0,         0.0,&  !n = 2
       250.516246,  256.504683,  258.898817,         0.0,         0.0,&  !n = 3
       335.508646,  338.025525,  339.022357,  339.095508,         0.0,&  !n = 4
       374.316936,  375.602652,  376.118426,         0.0,         0.0 &  !n = 5
       /)

  !Ionization potential = 2086., nGround = 1
  real, parameter, dimension(20) :: cExcitationAl11_I = (/&
       1575.033103, 1588.175427,         0.0,         0.0,         0.0,&  !n = 2
       1862.304461, 1865.813214, 1868.044929,         0.0,         0.0,&  !n = 3
       1961.045467, 1962.458887, 1963.401166, 1964.393040,         0.0,&  !n = 4
       2006.349288, 2007.055998, 2007.539536, 2008.481816,         0.0 &  !n = 5
       /)

  !Ionization potential = 2304., nGround = 1
  real, parameter, dimension(20) :: cExcitationAl12_I = (/&
       1727.735746, 1727.684912,         0.0,         0.0,         0.0,&  !n = 2
       2048.097243, 2048.082365, 2048.467956,         0.0,         0.0,&  !n = 3
       2160.175228, 2160.167789, 2160.331448, 2160.384761,         0.0,&  !n = 4
       2212.029134, 2212.025414, 2212.108483, 2212.137000, 2212.150638 &  !n = 5
       /)

  real, parameter, dimension(0:4, 2:5, 0:12) :: cExcitationAl_III = reshape((/&
       cExcitationAl0_I,  cExcitationAl1_I,  cExcitationAl2_I,&
       cExcitationAl3_I,  cExcitationAl4_I,  cExcitationAl5_I,&
       cExcitationAl6_I,  cExcitationAl7_I,  cExcitationAl8_I,&
       cExcitationAl9_I, cExcitationAl10_I, cExcitationAl11_I,&
       cExcitationAl12_I &
       /), (/5, 4, 13/))

  ! Xe, nZ = 54
  ! E.B.Saloman, Energy Levels and Observed Spectral Lines of Xenon,
  ! Xe I through Xe LIV

  !Ionization potential = 12.129874, nGround = 5
  real, parameter, dimension(0:4, 4:10) :: cExcitationXe0_II = reshape((/&
       0.0,       0.0,       0.0,  11.26270,       0.0,&  !n = 4
       0.0,       0.0,   9.89038,  11.57549,  11.58275,&  !n = 5
       8.31532,   9.58015,  10.97149,  11.74563,       0.0,&  !n = 6
       10.56206,  10.90157,  11.43877,  11.84806,       0.0,&  !n = 7
       11.25833,  11.42555,  11.66993,  11.91442,       0.0,&  !n = 8
       11.57991,  11.66281,  11.80076,  11.95984,       0.0,&  !n = 9
       11.74873,  11.79764,  11.88912,  11.99231,       0.0 &  !n = 10
       /), (/5, 7/))

  !Ionization potential = 21.20979, nGround = 5
  real, parameter, dimension(0:4, 4:8) :: cExcitationXe1_II = reshape((/&
       0.0,       0.0,       0.0,  17.22982,       0.0,&  !n = 4
       0.0,  11.26692,  11.82769,  18.60964,  18.77645,&  !n = 5
       11.53901,  13.86046,  16.80076,  19.32030,  19.44976,&  !n = 6
       16.43024,  17.24977,       0.0,       0.0,  19.85531,&  !n = 7
       18.26418,       0.0,       0.0,       0.0,       0.0 &  !n = 8
       /), (/5, 5/))

  !Ionization potential = 32.12295, nGround = 5
  real, parameter, dimension(0:3, 4:7) :: cExcitationXe2_II = reshape((/&
       0.0,       0.0,       0.0,  20.62542,&  !n = 4
       0.0,  12.18299,  13.83731,  24.46339,&  !n = 5
       15.06110,  18.19858,  22.60701,       0.0,&  !n = 6
       22.62497,       0.0,       0.0,       0.0 &  !n = 7
       /), (/4, 4/))

  !======================CORRECTION!!!!
  !The statistical weight from the Allen's
  !Astrophysical quantities book is incorrect (thanks to Rafael Rodrigez)
  !For terms p^2-p^4 and d^2-d^8 the used stat.weight is off by a factor of 2 to 50
  logical,public:: UseFullGroundState = .false. 
  !====================THE END OF INTERNAL DATABASE============================
  logical,public::UseDataBase = .false.
contains
  !============================================================================
  subroutine init_excitation

    use ModConst, ONLY: cRyToEV
    integer :: nZ

    real,dimension(1:nZMax) :: IonizPotential_I

    integer :: iL, iN, iZ, nGround
    integer :: InnerDegeneracy

    !Degeneracies of ground states for positive ions of the first 10 elements.
    !The data for the first 3 columns are taken from the book 
    !"Allen' Astrophysical Quantities, Edn iV" p. 33
    !by Allen, Clabon W. Editor: Cox, Arthur N.
    !Publisher: Springer (2000) (g_0 for Y=1, Y=2, Y=3

    !First 10 elements - full ionizations 
    integer, dimension(10,10)           :: cDegeneracy10_II
    integer, dimension(10,10),parameter :: cDegeneracy10Bad_II = reshape(  (/   &
                                !   I  !  II  !  III !  IV  !   V  !  VI  !  VII ! VIII !  IX  !   X  !
         2 ,     1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1     ,&  ! 1 - H
         1 ,     2 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1     ,&  ! 2 - He
         2 ,     1 ,    2 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1     ,&  ! 3 - Li
         1 ,     2 ,    1 ,    2 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1     ,&  ! 4 - Be
         6 ,     1 ,    2 ,    1 ,    2 ,    1 ,    1 ,    1 ,    1 ,    1     ,&  ! 5 - B
         9 ,     6 ,    1 ,    2 ,    1 ,    2 ,    1 ,    1 ,    1 ,    1     ,&  ! 6 - C
         4 ,     9 ,    6 ,    1 ,    2 ,    1 ,    2 ,    1 ,    1 ,    1     ,&  ! 7 - N
         9 ,     4 ,    9 ,    6 ,    1 ,    2 ,    1 ,    2 ,    1 ,    1     ,&  ! 8 - O
         6 ,     9 ,    4 ,    9 ,    6 ,    1 ,    2 ,    1 ,    2 ,    1     ,&  ! 9 - F
         1 ,     6 ,    9 ,    4 ,    9 ,    6 ,    1 ,    2 ,    1 ,    2   /),&  !10 - Ne
         (/10,10/))
    integer, dimension(10,10),parameter :: cDegeneracy10Good_II = reshape(  (/   &
                                !   I  !  II  !  III !  IV  !   V  !  VI  !  VII ! VIII !  IX  !   X  !
         2 ,     1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1     ,&  ! 1 - H
         1 ,     2 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1     ,&  ! 2 - He
         2 ,     1 ,    2 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1     ,&  ! 3 - Li
         1 ,     2 ,    1 ,    2 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1     ,&  ! 4 - Be
         6 ,     1 ,    2 ,    1 ,    2 ,    1 ,    1 ,    1 ,    1 ,    1     ,&  ! 5 - B
        15 ,     6 ,    1 ,    2 ,    1 ,    2 ,    1 ,    1 ,    1 ,    1     ,&  ! 6 - C
        20 ,    15 ,    6 ,    1 ,    2 ,    1 ,    2 ,    1 ,    1 ,    1     ,&  ! 7 - N
        15 ,    20 ,   15 ,    6 ,    1 ,    2 ,    1 ,    2 ,    1 ,    1     ,&  ! 8 - O
         6 ,    15 ,   20 ,    9 ,    6 ,    1 ,    2 ,    1 ,    2 ,    1     ,&  ! 9 - F
         1 ,     6 ,   15 ,    4 ,    9 ,    6 ,    1 ,    2 ,    1 ,    2   /),&  !10 - Ne
         (/10,10/))
    
    ! Degeneracies of ground states for positive ions of aluminium
    integer, dimension(13)              :: cDegeneracyAl_I
    integer, dimension(13),parameter :: cDegeneracyAlBad_I = (/ &
         6, 1, 2, 1, 6, 9, 4, 9, 6, 1, 2, 1, 2 /)
    integer, dimension(13),parameter :: cDegeneracyAlGood_I = (/ &
         6, 1, 2, 1, 6,15,20,15, 6, 1, 2, 1, 2 /)

    ! Degeneracies of ground states for positive ions of xenon
    ! Electron configurations used are confirmed in E.B.Saloman,
    ! Energy Levels and Observed Spectral Lines of Xenon, Xe I through Xe LIV
    integer, dimension(54)              :: cDegeneracyXe_I
    integer, parameter, dimension(54) :: cDegeneracyXeBad_I = (/ &
         1, 6,   9,  4, 9, 6, 1, 2, &
         1, 10, 21, 28, 25, 6, 25, 28, 21, 10, &
         1, 6,   9,  4, 9, 6, 1, 2,&
         1, 10, 21, 28, 25, 6, 25, 28, 21, 10, &
         1, 6, 9, 4, 9, 6, 1, 2, &
         1, 6, 9, 4, 9, 6, 1, 2,&
         1, 2 /)
    integer, dimension(54),parameter :: cDegeneracyXeGood_I = (/ &
         1, 6,  15, 20, 15,  6,  1,  2, &
         1, 10, 45,120,210,252,210,120,45, 10, &
         1, 6,  15, 20, 15,  6,  1,  2,&
         1, 10, 45,120,210,252,210,120,45, 10, &
         1, 6,  15, 20, 15, 6, 1, 2, &
         1, 6,  15, 20, 15, 6, 1, 2,&
         1, 2 /)
    
    
    !-------------------------------------------------------------------------
    if(.not.allocated(cExcitationEnergy_IIII)) allocate( &
         cExcitationEnergy_IIII(0:nExcitation-1,nExcitation,0:nZMax-1,nZMax),&
         cVirialCoeff4Energy_IIII(0:nExcitation-1,nExcitation,0:nZMax-1,nZMax)&
         , cDegeneracy_IIII(0:nExcitation-1,nExcitation,0:nZMax-1,nZMax))

    if(UseFullGroundState)then
       cDegeneracy10_II = cDegeneracy10Good_II
       cDegeneracyAl_I  = cDegeneracyAlGood_I
       cDegeneracyXe_I  = cDegeneracyXeGood_I
    else
       cDegeneracy10_II = cDegeneracy10Bad_II
       cDegeneracyAl_I  = cDegeneracyAlBad_I
       cDegeneracyXe_I  = cDegeneracyXeBad_I
    end if

    !Fill in the array ExcitationEnergy_III with tabulated or calculted excitation energies
    !for quantum states defined by n, l of atoms or ions of a particular element defined
    !by their charge state
    cExcitationEnergy_IIII   = 0.0
    cVirialCoeff4Energy_IIII = 0.0
    cDegeneracy_IIII = 0
    do nZ=1,10
       call for_particular_nz
    end do
    nZ = 13; call for_particular_nz
    nZ = 54; call for_particular_nz
  contains
    subroutine for_particular_nz
      call get_ioniz_potential(nZ, IonizPotential_I(1:nZ))

      do iZ = 0, nZ-1
         nGround = n_ground(iZ, nZ)

         do iN = nGround+1, nExcitation
            cExcitationEnergy_IIII(0:iN-1,iN,iZ,nZ) = IonizPotential_I(iZ+1) - &
                 cRyToEV * (real(iZ+1)/iN)**2
         end do
         InnerDegeneracy = 1
         select case (nZ)
         case (1:10)
            cDegeneracy_IIII(0,nGround,iZ,nZ)   = cDegeneracy10_II(iZ+1,nZ)
            if (iZ < nZ-1) InnerDegeneracy = cDegeneracy10_II(iZ+2,nZ)
         case (13)
            cDegeneracy_IIII(0,nGround,iZ,nZ)   = cDegeneracyAl_I(iZ+1)
            if (iZ < nZ-1) InnerDegeneracy = cDegeneracyAl_I(iZ+2)
         case (54)
            cDegeneracy_IIII(0,nGround,iZ,nZ)   = cDegeneracyXe_I(iZ+1)
            if (iZ < nZ-1) InnerDegeneracy = cDegeneracyXe_I(iZ+2)
         case default
            write(*,*)nZ
            call CON_stop('No such element found in the database'//&
                 ' of ground state degeneracies')
         end select

         do iN = nGround+1, nExcitation
            do iL = 0, iN-1
               cDegeneracy_IIII(iL,iN,iZ,nZ) = InnerDegeneracy * 2*(2*iL + 1)
            end do
         end do

      end do

      select case (nZ)
      case (2)
         where (cExcitationHe_III(0:4,2:5,0:nZ-1) /= 0.0)&
              cExcitationEnergy_IIII(0:4,2:5,0:nZ-1,nZ) = cExcitationHe_III(0:4,2:5,0:nZ-1)

      case (4)
         where (cExcitationBe_III(0:4,2:5,0:nZ-1) /= 0.0)&
              cExcitationEnergy_IIII(0:4,2:5,0:nZ-1,nZ) = cExcitationBe_III(0:4,2:5,0:nZ-1)
                  
         if(UseDataBase)call read_atomic_data(&
              4,3, cExcitationEnergy_IIII(:,:,0:3,4),cDegeneracy_IIII(:,:,0:3,4))
         call get_ioniz_potential(4,IonizPotential_I(1:4))
      case (6)
         where (cExcitationC_III(0:4,2:5,0:nZ-1) /= 0.0)&
              cExcitationEnergy_IIII(0:4,2:5,0:nZ-1,nZ) = cExcitationC_III(0:4,2:5,0:nZ-1)

      case (7)
         where (cExcitationN_III(0:4,2:5,0:nZ-1) /= 0.0)&
              cExcitationEnergy_IIII(0:4,2:5,0:nZ-1,nZ) = cExcitationN_III(0:4,2:5,0:nZ-1)

      case (8)
         where (cExcitationO_III(0:4,2:5,0:nZ-1) /= 0.0)&
              cExcitationEnergy_IIII(0:4,2:5,0:nZ-1,nZ) = cExcitationO_III(0:4,2:5,0:nZ-1)

      case (13)
         where (cExcitationAl_III(0:4,2:5,0:nZ-1) /= 0.0)&
              cExcitationEnergy_IIII(0:4,2:5,0:nZ-1,nZ) = cExcitationAl_III(0:4,2:5,0:nZ-1)

      case (54)
         where (cExcitationXe0_II(0:4,4:10) /= 0.0)&
              cExcitationEnergy_IIII(0:4,4:10,0,nZ) = cExcitationXe0_II(0:4,4:10)
         where (cExcitationXe1_II(0:4,4:8) /= 0.0)&
              cExcitationEnergy_IIII(0:4,4:8,1,nZ) = cExcitationXe1_II(0:4,4:8)
         where (cExcitationXe2_II(0:3,4:7) /= 0.0)&
              cExcitationEnergy_IIII(0:3,4:7,2,nZ) = cExcitationXe2_II(0:3,4:7)
         
         if(UseDataBase)call read_atomic_data(&
              54,19, cExcitationEnergy_IIII(:,:,0:19,54),cDegeneracy_IIII(:,:,0:19,54))
      end select

      if(UsePressureIonization)then    

         do iZ = 0, nZ-1
            nGround = n_ground(iZ, nZ)

            do iN = nGround, nExcitation
               do iL = 0, iN-1
                  cVirialCoeff4Energy_IIII(iL,iN,iZ,nZ) = IonizPotential_I(iZ+1) * &
                       real(nGround)**2 * real(iN - iL)**2 / real(iZ+1)**2
               end do
            end do
         end do
      end if
    end subroutine for_particular_nz
  end subroutine init_excitation
  !===============================
  integer function n_ground(iZ,nZ)
    integer,intent(in)::iZ,nZ
    !The principal quantum number of the outermost electron in bounded with an ion
    !with I electrons

    ! Works for H through Ar, Xe and Au (to work with gold uncomment the line with if)
    integer,parameter :: nGround_I(0:79) = (/ &
         0, &                                                                 ! For fully stripped ion
         1, 1, &                                                              !  2
         2, 2, 2, 2, 2, 2, 2, 2, &                                            !  8
         3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &              ! 18
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &              ! 18
         5, &           !\
         5, &           ! \
         5, &           !  \
         5, &           !   \
         5, &           !in Au:
         5, &           !14*(4f) -electrons 
         5, &           !
         5, &           !However in Xe
         4, &           !the first 8 of them are
         4, &           !2*(5s) and 6*(5p) electrons
         4, &           !   /
         4, &           !  /
         4, &           ! /
         4, &           !/
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &              ! 18
         6 /)                                                                 ! 1
    !1-6-s--- !---2-5--p-------!----3-5------d----------------!--------4f-----!
    !------------------------------------------------------
    n_ground = nGround_I(nZ - iZ)
    if(nZ==79.and.iZ<=54.and.iZ>=47)n_ground=4
  end function n_ground
  !======================================================================================
  integer function n_screened(iZ,nZ)
    !The number of electrons at inner shells
    integer,intent(in) :: iZ, nZ

    integer :: nGround
    integer,parameter :: nScreenedShell_I(0:6)  = (/0,0,2,10,28,46,78 /)
    !---------------------

    nGround = n_ground(iZ, nZ)
    n_screened = nScreenedShell_I(nGround)
    if(nZ==79.and.nGround==5)n_screened=60

  end function n_screened
  !=======================
  !\
  !Read atomic data from FAC output files. Collapses the unused quantum numbers
  !/
  subroutine read_atomic_data(&
       nZ,    &!The number of chemical element
       iZMax, &!The maximal charge state to be involved
       ExcitationEnergy_III, &
       Degeneracy_III) 
    use ModIoUnit, ONLY: io_unit_new
    use CRASH_ModAtomicNotation
    use ModUtilities,ONLY: split_string
    use CRASH_ModIonization,ONLY:put_ioniz_potential
    implicit none

    integer,intent(in)::nZ,iZMax
    real,intent(out),dimension(0:nExcitation-1,1:nExcitation,0:iZMax)::&
         ExcitationEnergy_III
    integer,intent(out),dimension(0:nExcitation-1,1:nExcitation,0:iZMax)::Degeneracy_III

    integer::iUnit

    integer::  iZ !Charge state

    character(len=30)::NameFile
    character(len=120)::NameRead
    integer,parameter::nHeader=3
    integer::iPosition,nLine,iLine,iAux,Degeneracy,iAux1,iAux2,iAux3

    !Principal and orbital quantum numbers
    integer::n,l
    real::Energy

    integer,parameter::MaxString = 20
    character(LEN=20),dimension(MaxString)::NameSplit_I

    integer::nString,iString
    integer::nInner,lInner,nGround,lGround
    character(len=1)::TypeL
    real::Energy0,Energy0Old
    !================
    Degeneracy_III       = 0.0
    ExcitationEnergy_III = 0.0
    do iZ=0,min(iZMax,nZ-1)

       NameFile = '../../../dataCRASH/AtomicData/' &
            // NameElement_I(nZ) // '/' // NameElement_I(nZ-iZ) // '.lev'
       write(*,'(a)')NameFile
       iUnit = io_unit_new()
       open(iUnit,file=NameFile,status='old')

       read(iUnit,'(a)')NameRead
       iPosition=index(NameRead,',')

       read( NameRead(iPosition+1:len(NameRead)),*)Energy0
       if(iZ>0)call put_ioniz_potential(nZ,iZ,Energy0 - Energy0Old)
       Energy0Old = Energy0

       do iLine=2,nHeader
          read(iUnit,*)
       end do
       read(iUnit,'(a)')NameRead
       iPosition=index(NameRead,'=')
       read( NameRead(iPosition+1:len(NameRead)),*)nLine
       read(iUnit,*)

       LINES: do iLine=1,nLine
          read(iUnit,'(a)')NameRead
          !First, second and fourth column are unusable.
          read(NameRead,*)iAux1,iAux2,Energy,iAux3,iAux,Degeneracy 
          Degeneracy = Degeneracy + 1
          n = iAux/100
          l = iAux - 100 * n
          if(iLine==1)then
             nGround = n
             lGround = l
          else
             !Skip the shell distrubution (like 1*2 )
             call split_string(NameRead,MaxString,NameSplit_I,nString)
             
             STRING:   do iString = 7,nString
                if(index(NameSplit_I(iString),'*')>0)CYCLE STRING
                NameRead = NameSplit_I(iString)
                !Now the NameRead variable starts from the term of the first unclosed shell
                if(NameRead(2:2)=='0')then
                   
                   !To avoid the problem with n=10
                   
                   read( NameRead(1:2),'(i2)')nInner
                   lInner = l_orbital(NameRead(3:3))
                else
                   
                   !Read the principle quantum number for the inner shell with the vacancy
                   
                   read( NameRead(1:1),'(i1)')nInner
                   lInner = l_orbital(NameRead(2:2))
                end if
                if(nInner < nGround .or. (nInner==nGround .and. lInner <lGround))then
                   n = nInner
                   l = lInner
                   if(.not.UseCoreElectron)CYCLE LINES
                end if
                EXIT STRING
             end do STRING
             if((.not.UseDeltaNEq0Transition)&
                  .and.(n == nGround.and.l>lGround))CYCLE LINES
          end if
          ExcitationEnergy_III(l,n,iZ) = (ExcitationEnergy_III(l,n,iZ) * Degeneracy_III(l,n,iZ) +&
               Energy * Degeneracy)/(Degeneracy_III(l,n,iZ) + Degeneracy)
          Degeneracy_III(l,n,iZ) = Degeneracy_III(l,n,iZ) + Degeneracy
          
       end do LINES
       close(iUnit)
    end do
  end subroutine read_atomic_data

end module CRASH_ModExcitationData
