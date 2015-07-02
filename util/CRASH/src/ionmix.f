ABJTIONMIX.  IONMIX, A CODE FOR COMPUTING THE EQUATION OF STATE AND     ABJT0000
1   RADIATIVE PROPERTIES OF LTE AND NON-LTE PLASMAS.  J.J. MACFARLANE.  ABJT0000
REF. IN COMP. PHYS. COMMUN. 56 (1989) 259                               ABJT0000
c ============== BEGIN TEST DECK ====================================   ABJT0001
    $data                                                               ABJT0002
         ngases = 1,                                                    ABJT0003
         izgas(1) =  2,                                                 ABJT0004
         atomwt(1) = 4.,                                                ABJT0005
         fracsp(1) = 1.,                                                ABJT0006
         ntemp = 10,                                                    ABJT0007
         ndens = 1,                                                     ABJT0008
         dlgtmp = 0.25,                                                 ABJT0009
         dlgden = 1.,                                                   ABJT0010
         tplsma(1) = 1.,                                                ABJT0011
         densnn(1) = 3.5e16,                                            ABJT0012
         ngrups = 20,                                                   ABJT0013
         nptspg = 10,                                                   ABJT0014
         nfrqbb = 9,                                                    ABJT0015
    $end                                                                ABJT0016
c ============== END TEST DECK ======================================   ABJT0017
c ============== BEGIN PROGRAM ======================================   ABJT0018
      program ionmix                                                    ABJT0019
                                                                        ABJT0020
c ... driver program                                                    ABJT0021
                                                                        ABJT0022
c ... set the maximum number of:                                        ABJT0023
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT0024
c       temperatures (mxtemp), densities (mxdens),                      ABJT0025
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT0026
                                                                        ABJT0027
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT0028
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT0029
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT0030
                                                                        ABJT0031
                                                                        ABJT0032
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT0033
c ................................................................      ABJT0034
                                                                        ABJT0035
      common / opacs / nptspg, ngrups, engrup(mxgrps+1)                 ABJT0036
c ................................................................      ABJT0037
                                                                        ABJT0038
      common / params / npqmax(mxgass), npmaxp, nfrqbb,                 ABJT0039
     2                  npring(100), neopen(5,100), bbnton(6,55),       ABJT0040
     3                  defpot(27,55), noccdf(27,55)                    ABJT0041
c ................................................................      ABJT0042
                                                                        ABJT0043
      common / tdmesh / ntemp,ndens,ntrad,                              ABJT0044
     2                  tplsma(mxtemp),densnn(mxdens),trad(mxtemp),     ABJT0045
     2                  dlgtmp,dlgden,dlgtrd                            ABJT0046
c ................................................................      ABJT0047
                                                                        ABJT0048
c ... type declaration statements                                       ABJT0049
                                                                        ABJT0050
      dimension enrgy(mxtemp,mxdens),heatcp(mxtemp,mxdens)              ABJT0051
      dimension dzdt(mxtemp,mxdens),densne(mxtemp,mxdens)               ABJT0052
      dimension dedden(mxtemp,mxdens)                                   ABJT0053
      dimension opgpa(mxtemp,mxdens,mxgrps),opma(mxtemp,mxdens)         ABJT0054
      dimension opgpe(mxtemp,mxdens,mxgrps),opme(mxtemp,mxdens)         ABJT0055
      dimension orgp(mxtemp,mxdens,mxgrps),orm(mxtemp,mxdens)           ABJT0056
      dimension oppat(mxgrps),oppet(mxgrps),oprt(mxgrps)                ABJT0057
      dimension op2tp(mxtemp,mxdens,mxtemp),op2tpt(mxtemp)              ABJT0058
      dimension op2tr(mxtemp,mxdens,mxtemp),op2trt(mxtemp)              ABJT0059
      dimension photen(mxphot),abscfs(mxphot),sctcfs(mxphot)            ABJT0060
      dimension emscfs(mxphot)                                          ABJT0061
                                                                        ABJT0062
c *****************************************************************     ABJT0063
c                                                                       ABJT0064
c                           begin execution                             ABJT0065
c                                                                       ABJT0066
c *****************************************************************     ABJT0067
                                                                        ABJT0068
                                                                        ABJT0069
c ... open files                                                        ABJT0070
      open ( unit=5,file='ionmxinp',status='old',form='formatted' )     ABJT0071
      open ( unit=6,file='ionmxout',status='new',form='formatted' )     ABJT0072
                                                                        ABJT0073
c ... read input                                                        ABJT0074
      call input                                                        ABJT0075
                                                                        ABJT0076
c ... loop over densities                                               ABJT0077
                                                                        ABJT0078
      do 200 idens=1,ndens                                              ABJT0079
                                                                        ABJT0080
c ...    loop over temperature                                          ABJT0081
                                                                        ABJT0082
         do 100 itemp=1,ntemp                                           ABJT0083
                                                                        ABJT0084
            call eos ( tplsma(itemp),densnn(idens),                     ABJT0085
     &                         enrgy(itemp,idens),heatcp(itemp,idens),  ABJT0086
     &                         dzdt(itemp,idens),densne(itemp,idens),   ABJT0087
     &                         dedden(itemp,idens) )                    ABJT0088
                                                                        ABJT0089
            if ( isw(2) .eq. 0 ) then                                   ABJT0090
                                                                        ABJT0091
c ...          set up the group structure for the opacities unless      ABJT0092
c              it was specified by the user                             ABJT0093
               if ( isw(13) .eq. 2 ) then                               ABJT0094
                 engrup(1) = 0.01 * tplsma(itemp)                       ABJT0095
                 engrup(ngrups+1) = 100. * tplsma(itemp)                ABJT0096
                 if ( ngrups .gt. 1 ) then                              ABJT0097
                   elnmin = log( engrup(1) )                            ABJT0098
                   elnmax = log( engrup(ngrups+1) )                     ABJT0099
                   delog = ( elnmax-elnmin ) / ngrups                   ABJT0100
                   elog = elnmin                                        ABJT0101
                   do 50 igrp=2,ngrups                                  ABJT0102
                     elog = elog + delog                                ABJT0103
                     engrup(igrp) = exp( elog )                         ABJT0104
   50              continue                                             ABJT0105
                 endif                                                  ABJT0106
               endif                                                    ABJT0107
                                                                        ABJT0108
c ...          set up the mesh points at which the absorption           ABJT0109
c              coefficients are to be evaluated                         ABJT0110
               call meshhv ( tplsma(itemp),densnn(idens),densne,        ABJT0111
     &                       engrup,ngrups,nptspg,                      ABJT0112
     &                                              photen,nphot )      ABJT0113
                                                                        ABJT0114
c ...          calculate the absorption coefficients                    ABJT0115
               call abscon ( tplsma(itemp),densnn(idens),               ABJT0116
     &                       densne(itemp,idens),photen,nphot,          ABJT0117
     &                                          abscfs,emscfs,sctcfs )  ABJT0118
                                                                        ABJT0119
c ...          compute the group opacities                              ABJT0120
               call opacys ( tplsma(itemp),densnn(idens),               ABJT0121
     &                      densne(itemp,idens),photen,nphot,abscfs,    ABJT0122
     &                      emscfs,sctcfs,engrup,ngrups,ntrad,trad,     ABJT0123
     &                                 oppat,oppet,oprt,                ABJT0124
     &                                 oppma,oppme,oprm,culrat,         ABJT0125
     &                                 op2tpt,op2trt )                  ABJT0126
                                                                        ABJT0127
c ...          save opacities in arrays                                 ABJT0128
               opma(itemp,idens) = oppma                                ABJT0129
               opme(itemp,idens) = oppme                                ABJT0130
               orm(itemp,idens)  = oprm                                 ABJT0131
               do 60 ig=1,ngrups                                        ABJT0132
                  opgpa(itemp,idens,ig) = oppat(ig)                     ABJT0133
                  opgpe(itemp,idens,ig) = oppet(ig)                     ABJT0134
                  orgp(itemp,idens,ig)  = oprt(ig)                      ABJT0135
  60           continue                                                 ABJT0136
               if ( ntrad.ne.0 ) then                                   ABJT0137
                  do 70 itrad=1,ntrad                                   ABJT0138
                     op2tp(itemp,idens,itrad) = op2tpt(itrad)           ABJT0139
                     op2tr(itemp,idens,itrad) = op2trt(itrad)           ABJT0140
   70             continue                                              ABJT0141
               endif                                                    ABJT0142
                                                                        ABJT0143
            endif                                                       ABJT0144
                                                                        ABJT0145
c ...       write results to output files                               ABJT0146
            call owt1 ( tplsma(itemp),densnn(idens),                    ABJT0147
     &                  densne(itemp,idens),enrgy(itemp,idens),         ABJT0148
     &                  heatcp(itemp,idens),dzdt(itemp,idens),          ABJT0149
     &                  dedden(itemp,idens),oppat,oppet,oprt,           ABJT0150
     &                  ngrups,engrup,oppma,oppme,oprm,culrat )         ABJT0151
                                                                        ABJT0152
  100    continue                                                       ABJT0153
                                                                        ABJT0154
  200 continue                                                          ABJT0155
                                                                        ABJT0156
c ... print out final results                                           ABJT0157
                                                                        ABJT0158
      call owtf ( densne,enrgy,heatcp,dzdt,dedden,                      ABJT0159
     &            oppma,oppme,oprm,opgpa,opgpe,orgp,op2tp,op2tr )       ABJT0160
                                                                        ABJT0161
                                                                        ABJT0162
                                                                        ABJT0163
      end                                                               ABJT0164
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   ABJT0165
                                                                        ABJT0166
      block data menu                                                   ABJT0167
                                                                        ABJT0168
c ... this block data routine is used to initialize constants at        ABJT0169
c     the beginning of the computation, and to provide miscellaneous    ABJT0170
c     information, such as common block variable definitions.           ABJT0171
c                                                                       ABJT0172
                                                                        ABJT0173
c ... set the maximum number of:                                        ABJT0174
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT0175
c       temperatures (mxtemp), densities (mxdens),                      ABJT0176
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT0177
                                                                        ABJT0178
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT0179
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT0180
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT0181
                                                                        ABJT0182
c                                                                       ABJT0183
                                                                        ABJT0184
      common / consts / pi, avgdro, sbcon, hplank                       ABJT0185
c ...............................................................       ABJT0186
c                                                                       ABJT0187
c       pi      -  3.14159                                              ABJT0188
c       avgdro  -  6.022e23    Avogadro's number                        ABJT0189
c       sbcon   -  5.667e-5    Stefan-Boltzmann constant (erg/cm**2/s)  ABJT0190
c       hplank  -  4.136e-15   Planck's constant (eV sec)               ABJT0191
c ....................................................................  ABJT0192
                                                                        ABJT0193
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT0194
c ................................................................      ABJT0195
c                                                                       ABJT0196
c       isw     -  control switches                                     ABJT0197
c       con     -  real constants used throughout IONMIX                ABJT0198
c       iplot   -  control switch for plot data files                   ABJT0199
c       dtheat  -  fractional temperature increment for heat capacity   ABJT0200
c                  calculation                                          ABJT0201
c       critsc  -  critical value for switching between Saha and        ABJT0202
c                  coronal ionization models                            ABJT0203
c ....................................................................  ABJT0204
                                                                        ABJT0205
      common / dbugcm / nobug,ncycld(10),icycld(10)                     ABJT0206
                                                                        ABJT0207
      logical nobug                                                     ABJT0208
c ................................................................      ABJT0209
c                                                                       ABJT0210
c       nobug   -  logical flag; "true" if no debug output is requested ABJT0211
c       ncycld  -  debug output is printed every "ncycld"th cycle       ABJT0212
c       icycld  -  cycle counter for each subroutine name               ABJT0213
c                                                                       ABJT0214
                                                                        ABJT0215
      common / gases / ngases,izgas(mxgass),atomwt(mxgass),             ABJT0216
     2                 avgatw,avgatn,                                   ABJT0217
     3                 fracsp(mxgass),fraciz(mxgass,mxatom),            ABJT0218
     4                 fraclv(mxgass,mxatom,20),pot(mxgass,mxatom),     ABJT0219
     5                 numocc(mxgass,mxatom)                            ABJT0220
c ................................................................      ABJT0221
c                                                                       ABJT0222
c       ngases  -  number of gas species (max. = 5)                     ABJT0223
c       izgas   -  atomic number of each gas species                    ABJT0224
c       atomwt  -  atomic weight of each gas species                    ABJT0225
c       avgatw  -  average atomic weight                                ABJT0226
c       avgatn  -  average atomic number                                ABJT0227
c       fracsp  -  number fraction of each gas species, relative to     ABJT0228
c                  the total number of nuclei                           ABJT0229
c       fraciz  -  number fraction of nuclei of each ionization state,  ABJT0230
c                  relative to the number of nuclei of the particular   ABJT0231
c                  gas species                                          ABJT0232
c       fraclv  -  number fraction of ions in each atomic level,relativeABJT0233
c                  to the number of ions in the particular ionization   ABJT0234
c                  state                                                ABJT0235
c       pot     -  ionization potentials for each ionization state and  ABJT0236
c                  each gas species (eV)                                ABJT0237
c       numocc  -  electron shell occupation number for each ionization ABJT0238
c                  state and each gas species                           ABJT0239
c                                                                       ABJT0240
                                                                        ABJT0241
      common / opacs / nptspg, ngrups, engrup(mxgrps+1)                 ABJT0242
c ................................................................      ABJT0243
c                                                                       ABJT0244
c       nptspg  -  number of mesh points per group for which absorption ABJT0245
c                  coefficients are to be computed                      ABJT0246
c       ngrups  -  number of photon energy groups                       ABJT0247
c       engrup  -  photon energy group boundaries (eV)                  ABJT0248
c                                                                       ABJT0249
                                                                        ABJT0250
      common / params / npqmax(mxgass), npmaxp, nfrqbb,                 ABJT0251
     2                  npring(100), neopen(5,100), bbnton(6,55),       ABJT0252
     3                  defpot(27,55), noccdf(27,55)                    ABJT0253
c ................................................................      ABJT0254
c                                                                       ABJT0255
c       npqmax(l)- maximum principal quantum number used in computing   ABJT0256
c                  absorption coefficients of the l'th gas              ABJT0257
c       npmaxp  -  maximum principal quantum number used in computing   ABJT0258
c                  populations                                          ABJT0259
c       nfrqbb  -  number of photon energy points near a line center    ABJT0260
c                  at which absorption coefficients will be computed    ABJT0261
c       npring  -  principal quantum number of the valence electrons    ABJT0262
c                  in their ground state, as a function of the number   ABJT0263
c                  of bound electrons                                   ABJT0264
c       neopen  -  number of open holes available in the valence shell  ABJT0265
c       bbnton  -  constants for bound-bound transitions in which       ABJT0266
c                  the principal quantum number does not change         ABJT0267
c       defpot  -  default ionization potentials (eV)                   ABJT0268
c       noccdf  -  default electron shell occupation numbers            ABJT0269
c                                                                       ABJT0270
                                                                        ABJT0271
      common / strngs / header, namedb(10)                              ABJT0272
                                                                        ABJT0273
      character header*80, namedb*6                                     ABJT0274
c ................................................................      ABJT0275
c                                                                       ABJT0276
c       namedb  -  array of subroutine names for which debug output     ABJT0277
c                  is requested                                         ABJT0278
c       header  -  header record for output files; contains time and    ABJT0279
c                  date of calculation                                  ABJT0280
c                                                                       ABJT0281
                                                                        ABJT0282
      common / tdmesh / ntemp,ndens,ntrad,                              ABJT0283
     2                  tplsma(mxtemp),densnn(mxdens),trad(mxtemp),     ABJT0284
     2                  dlgtmp,dlgden,dlgtrd                            ABJT0285
c ................................................................      ABJT0286
c                                                                       ABJT0287
c       ntemp   -  number of plasma temperature points for              ABJT0288
c                  calculation                                          ABJT0289
c       ndens   -  number of density points for calculation             ABJT0290
c       ntrad   -  number of radiation temperature points for           ABJT0291
c                  calculation                                          ABJT0292
c       tplsma  -  plasma temperatures (eV)                             ABJT0293
c       densnn  -  number density of all nuclei (cm**-3)                ABJT0294
c       trad    -  radiation temperatures (eV)                          ABJT0295
c       dlgtmp  -  logarithmic increment in plasma temperature          ABJT0296
c       dlgden  -  logarithmic increment in number density              ABJT0297
c       dlgtrd  -  logarithmic increment in radiation temperature       ABJT0298
c                                                                       ABJT0299
c *******************************************************************   ABJT0300
c                                                                       ABJT0301
c                        begin data initialization                      ABJT0302
c                                                                       ABJT0303
c *******************************************************************   ABJT0304
                                                                        ABJT0305
      data  pi / 3.14159 /, avgdro / 6.0222e23 /, sbcon / 5.6697e-5 /   ABJT0306
      data  hplank / 4.136e-15 /                                        ABJT0307
                                                                        ABJT0308
c ... principal quantum number of valence electrons in their            ABJT0309
c     ground state, as a function of the number of electrons attached   ABJT0310
c     to the ion                                                        ABJT0311
      data  npring / 2*1, 8*2, 18*3, 32*4, 40*5 /                       ABJT0312
                                                                        ABJT0313
c ... neopen(i,j) =  number of available spots in a principal quantum   ABJT0314
c     level open to an electron after it is captured, as a function     ABJT0315
c     of the principal quantum number "i", and the number of electrons, ABJT0316
c     "j", attached to the ion before capture.                          ABJT0317
      data (neopen(1,j),j=1,100) /  2,1, 98*0 /                         ABJT0318
      data (neopen(2,j),j=1,100) /  2*8, 8,7,6,5,4,3,2,1, 90*0 /        ABJT0319
      data (neopen(3,j),j=1,100) /  10*18, 18,17,16,15,14,13,12,11,10,  ABJT0320
     2      9, 8, 7, 6, 5, 4, 3, 2, 1, 72*0 /                           ABJT0321
      data (neopen(4,j),j=1,100) /  28*32, 32,31,30,                    ABJT0322
     3 29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,     ABJT0323
     4  9, 8, 7, 6, 5 , 4, 3, 2, 1, 40*0 /                              ABJT0324
      data (neopen(5,j),j=1,100) /    60*50,                            ABJT0325
     2 50,49,48,47,46,45,44,43,42,41,40,39,38,37,36,35,34,33,32,31,30,  ABJT0326
     3 29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11 /       ABJT0327
                                                                        ABJT0328
c ... parameters for bound-bound transitions when "delta n" = 0.        ABJT0329
c     these values are taken from Post, et.al., Atomic and Nuclear      ABJT0330
c     Data Tables, Vol. 20, p.397 (1977).                               ABJT0331
                                                                        ABJT0332
      data (bbnton(1,n),n=1,55) /                                       ABJT0333
     1      0.,  0.,-.02,-.01,-.04, .04,-.01,-.01, .01,  0.,            ABJT0334
     2     .01, .10, .25, .17, .08,-.02,-.14,-.27,-.29,-.30,            ABJT0335
     3    -.30,-.29,-.27,-.24,-.20,-.14,-.08,  0., .97,1.96,            ABJT0336
     4    1.92,1.89,1.86,1.83,1.78,1.73,1.41,1.05, .67, .26,            ABJT0337
     5    -.17,-.64,-1.14,-1.67,-2.26,-2.88,-2.90,-2.83,-2.72,-2.61,    ABJT0338
     6    -2.45,-2.27,-2.05,-1.81,-1.55 /                               ABJT0339
      data (bbnton(2,n),n=1,55) /                                       ABJT0340
     1      0.,  0.,2.00,4.33,4.68,3.60,2.85,2.40,1.30,  0.,            ABJT0341
     2    10.5,20.0,16.4,24.7,34.1,44.8,56.8,70.3,68.6,65.8,            ABJT0342
     3    61.9,56.9,50.7,45.5,34.4,24.2,12.8,  0.,-25.8,-53.1,          ABJT0343
     4    -28.1,-2.97,23.0,50.3,79.1,109.,142.,176.,214.,257.,          ABJT0344
     5    300.,348.,401.,457.,518.,584.,572.,551.,525.,497.,            ABJT0345
     6    463.,427.,385.,340.,291. /                                    ABJT0346
      data (bbnton(3,n),n=1,55) /                                       ABJT0347
     1      0.,  0.,2.04,4.49,6.80,8.16,6.80,10.6,13.1,  0.,            ABJT0348
     2    2.30,3.88,5.71,5.44,8.16,6.80,8.30,4.32,5.11,6.04,            ABJT0349
     3    7.08,8.12,9.69,11.3,13.7,14.9,17.1,  0.,.015,.019,            ABJT0350
     4    .056,.120,.213,.334,.466,.666,.910,1.25,1.69,1.98,            ABJT0351
     5    3.03,3.92,5.19,6.40,8.05,9.80,10.8,11.9,13.2,14.7,            ABJT0352
     6    15.9,18.1,19.2,21.0,23.3 /                                    ABJT0353
      data (bbnton(4,n),n=1,55) /                                       ABJT0354
     1      1.,  1., 1.0, .93, .86, .80, .87, .77, .77,  1.,            ABJT0355
     2     .99, .90, .83, .87, .77, .82, .79,1.06,1.03,1.00,            ABJT0356
     3     .07, .94, .91, .88, .85, .83, .80,  1.,2.46,2.40,            ABJT0357
     4    2.14,1.96,1.83,1.72,1.64,1.57,1.49,1.41,1.34,1.30,            ABJT0358
     5    1.19,1.13,1.06,1.01, .95, .91, .89, .87, .85, .83,            ABJT0359
     6     .82, .80, .78, .76, .74 /                                    ABJT0360
      data (bbnton(5,n),n=1,55) /                                       ABJT0361
     1      0.,  0., 1.2, .91, .89, .92, .86, .87, .85,  0.,            ABJT0362
     2     .72, .78, .76, .78, .75, .70, .68, .65, .65, .70,            ABJT0363
     3      .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,            ABJT0364
     4      .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,            ABJT0365
     5      .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,  .7,            ABJT0366
     6      .7,  .7,  .7,  .7,  .7 /                                    ABJT0367
      data (bbnton(6,n),n=1,55) /                                       ABJT0368
     1      0.,  0., .54, .77, .58, .57, .41, .63, .66,  0.,            ABJT0369
     2     .97, .96, .88, .86, .87, .85, .83, .81, .80, .80,            ABJT0370
     3      .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,            ABJT0371
     4      .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,            ABJT0372
     5      .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,  .8,            ABJT0373
     6      .8,  .8,  .8,  .8,  .8 /                                    ABJT0374
                                                                        ABJT0375
c ... default ionization potentials for ground state ions,              ABJT0376
c     taken from Carlson, et.al., Atomic Data, 2, p. 63 (1970).         ABJT0377
                                                                        ABJT0378
c ... H, Xe                                                             ABJT0379
      data (defpot(1,j),j=1,55) /                                       ABJT0380
     1   13.6,12.0,23.5,35.1,46.7,59.7,71.8,98.1,112.3,170.8,           ABJT0381
     2   201.7,232.6,263.5,294.4,325.3,358.3,389.6,420.9,452.2,572.5,   ABJT0382
     3   607.7,642.9,678.1,726.0,762.4,852.7,890.6,1394.,1491.,1587.,   ABJT0383
     4   1684.,1781.,1877.,1987.,2085.,2183.,2281.,2548.,2637.,2726.,   ABJT0384
     5   2814.,3001.,3093.,3296.,3386.,7224.,7491.,7758.,8024.,8617.,   ABJT0385
     6   8899.,9330.,9569.,39250.,40270. /                              ABJT0386
c ... He, I                                                             ABJT0387
      data (defpot(2,j),j=1,55) /                                       ABJT0388
     1   25.0,54.3,                                                     ABJT0389
     2   53*0. /                                                        ABJT0390
c ... Li, Te                                                            ABJT0391
      data (defpot(3,j),j=1,55) /                                       ABJT0392
     1   5.34,74.5,121.9,                                               ABJT0393
     2   52*0. /                                                        ABJT0394
c ... Be, Sb                                                            ABJT0395
      data (defpot(4,j),j=1,55) /                                       ABJT0396
     1   8.42,18.7,149.3,214.9,                                         ABJT0397
     2   51*0. /                                                        ABJT0398
c ... B, Sn                                                             ABJT0399
      data (defpot(5,j),j=1,55) / 55*0. /                               ABJT0400
c ... C, In                                                             ABJT0401
      data (defpot(6,j),j=1,55) /                                       ABJT0402
     1   10.6,26.1,50.4,67.6,374.0,475.6,                               ABJT0403
     2   49*0. /                                                        ABJT0404
c ... N, Cd                                                             ABJT0405
      data (defpot(7,j),j=1,55) /                                       ABJT0406
     1   13.4,32.0,51.6,82.6,103.1,524.0,643.3,                         ABJT0407
     2   48*0. /                                                        ABJT0408
c ... O, Ag                                                             ABJT0409
      data (defpot(8,j),j=1,55) /                                       ABJT0410
     1   16.4,38.2,60.7,82.8,121.9,145.8,698.8,836.0,                   ABJT0411
     2   47*0. /                                                        ABJT0412
c ... F, Pd                                                             ABJT0413
      data (defpot(9,j),j=1,55) / 55*0. /                               ABJT0414
c ... Ne, Rh                                                            ABJT0415
      data (defpot(10,j),j=1,55) /                                      ABJT0416
     1   23.1,51.3,79.2,107.6,135.9,164.1,221.8,252.4,1123.,1296.,      ABJT0417
     2   45*0. /                                                        ABJT0418
c ... Mg, Tc                                                            ABJT0419
      data (defpot(12,j),j=1,55) /                                      ABJT0420
     1   6.90,15.3,78.7,118.4,158.1,197.8,237.8,277.7,356.7,396.2,      ABJT0421
     2   1671.,1880.,                                                   ABJT0422
     3   43*0. /                                                        ABJT0423
c ... Si, Nb                                                            ABJT0424
      data (defpot(14,j),j=1,55) /                                      ABJT0425
     1   7.26,17.0,34.3,46.7,159.8,210.5,261.3,312.0,364.0,415.1,       ABJT0426
     2   503.7,552.2,2324.,2569.,                                       ABJT0427
     3   41*0. /                                                        ABJT0428
c ... S, Y                                                              ABJT0429
      data (defpot(16,j),j=1,55) /                                      ABJT0430
     1   11.4,24.5,38.0,51.1,76.7,92.6,265.8,327.4,389.1,450.7,         ABJT0431
     2   514.1,576.2,675.6,733.0,3082.,3364.,                           ABJT0432
     3   39*0. /                                                        ABJT0433
c ... Ar, Rb                                                            ABJT0434
      data (defpot(18,j),j=1,55) /                                      ABJT0435
     1   16.0,32.3,48.7,65.0,81.6,98.0,133.,153.,396.,469.,             ABJT0436
     2   541.,614.,689.,762.,873.,939.,3947.,4264.,                     ABJT0437
     3   37*0. /                                                        ABJT0438
c ... K, Kr                                                             ABJT0439
      data (defpot(19,j),j=1,55) /                                      ABJT0440
     1   19*0.,                                                         ABJT0441
     2   14.0,27.9,41.8,55.7,70.3,84.5,116.,133.,219.,269.,             ABJT0442
     3   318.,367.,416.,465.,515.,565.,614.,664.,837.,887.,             ABJT0443
     4   937.,987.,1046.,1097.,1219.,1271.,2728.,2897.,3065.,3234.,     ABJT0444
     5   3457.,3630.,3872.,4021.,16750.,17400. /                        ABJT0445
c ... Fe, Cu                                                            ABJT0446
      data (defpot(26,j),j=1,55) /                                      ABJT0447
     1   7.17,15.7,32.7,57.4,83.0,108.1,133.2,158.3,240.2,271.7,        ABJT0448
     2   303.1,334.6,372.2,404.0,472.4,506.0,1168.,1283.,1399.,1514.,   ABJT0449
     3   1644.,1760.,1923.,2025.,8503.,8969.,                           ABJT0450
     4   29*0. /                                                        ABJT0451
                                                                        ABJT0452
                                                                        ABJT0453
c ... The corresponding occupation numbers (again, from Carlson)        ABJT0454
                                                                        ABJT0455
      data (noccdf(1,j),j=1,55) /                                       ABJT0456
     1   1,4,3,2,1,2,1,2,1,6,5,4,3,2,1,4,3,2,1,4,3,2,1,2,1,2,1,6,       ABJT0457
     2   5,4,3,2,1,4,3,2,1,4,3,2,1,2,1,2,1,4,3,2,1,2,1,2,1,2,1 /        ABJT0458
      data (noccdf(2,j),j=1,55) /                                       ABJT0459
     1   2,1, 53*0 /                                                    ABJT0460
      data (noccdf(3,j),j=1,55) /                                       ABJT0461
     1   1,2,1, 52*0 /                                                  ABJT0462
      data (noccdf(4,j),j=1,55) /                                       ABJT0463
     1   2,1,2,1, 51*0 /                                                ABJT0464
      data (noccdf(5,j),j=1,55) / 55*0 /                                ABJT0465
      data (noccdf(6,j),j=1,55) /                                       ABJT0466
     1   2,1,2,1,2,1, 49*0 /                                            ABJT0467
      data (noccdf(7,j),j=1,55) /                                       ABJT0468
     1   2,1,1,2,1,2,1, 48*0 /                                          ABJT0469
      data (noccdf(8,j),j=1,55) /                                       ABJT0470
     1   2,1,2,1,2,1,2,1, 47*0 /                                        ABJT0471
      data (noccdf(9,j),j=1,55) / 55*0 /                                ABJT0472
      data (noccdf(10,j),j=1,55) /                                      ABJT0473
     1   4,3,2,1,2,1,2,1,2,1, 45*0 /                                    ABJT0474
      data (noccdf(12,j),j=1,55) /                                      ABJT0475
     1   2,1,4,3,2,1,2,1,2,1,2,1, 43*0 /                                ABJT0476
      data (noccdf(14,j),j=1,55) /                                      ABJT0477
     1   2,1,2,1,4,3,2,1,2,1,2,1,2,1, 41*0 /                            ABJT0478
      data (noccdf(16,j),j=1,55) /                                      ABJT0479
     1   2,1,2,1,2,1,4,3,2,1,2,1,2,1,2,1, 39*0 /                        ABJT0480
      data (noccdf(18,j),j=1,55) /                                      ABJT0481
     1   4,3,2,1,2,1,2,1,4,3,2,1,2,1,2,1,2,1, 37*0 /                    ABJT0482
      data (noccdf(19,j),j=1,55) /                                      ABJT0483
     1   19*0, 4,3,2,1,2,1,2,1,6,5,4,3,2,1,4,3,2,1,4,3,2,1,2,1,2,1,     ABJT0484
     2   4,3,2,1,2,1,2,1,2,1 /                                          ABJT0485
      data (noccdf(26,j),j=1,55) /                                      ABJT0486
     1   2,1,2,1,4,3,2,1,4,3,2,1,2,1,2,1,4,3,2,1,2,1,2,1,2,1, 29*0 /    ABJT0487
                                                                        ABJT0488
                                                                        ABJT0489
      end                                                               ABJT0490
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   ABJT0491
                                                                        ABJT0492
      subroutine abscon ( tp,densnn,densne,photen,nphot,                ABJT0493
     &                                         abscfs,emscfs,sctcfs )   ABJT0494
                                                                        ABJT0495
c ... this routine calculates the absorption, emission, and scattering  ABJT0496
c     coefficients for an array of photon energies                      ABJT0497
c                                                                       ABJT0498
c ... input variables:                                                  ABJT0499
c       tp      -  plasma temperature (eV)                              ABJT0500
c       densnn  -  number density of all nuclei (cm**-3)                ABJT0501
c       densne  -  electron density (cm**-3)                            ABJT0502
c       photen  -  array of photon energies at which the absorption     ABJT0503
c                  coefficient will be evaluated (eV)                   ABJT0504
c       nphot   -  number of elements in the "photen" array             ABJT0505
c                                                                       ABJT0506
c ... output variables:                                                 ABJT0507
c       abscfs  -  array of absorption coefficients (cm**-1)            ABJT0508
c       emscfs  -  array of emission coefficients (cm**-1)              ABJT0509
c       sctcfs  -  array of scattering coefficients (cm**-1)            ABJT0510
c                                                                       ABJT0511
                                                                        ABJT0512
c ... set the maximum number of:                                        ABJT0513
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT0514
c       temperatures (mxtemp), densities (mxdens),                      ABJT0515
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT0516
                                                                        ABJT0517
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT0518
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT0519
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT0520
                                                                        ABJT0521
                                                                        ABJT0522
      common / consts / pi, avgdro, sbcon, hplank                       ABJT0523
c ...............................................................       ABJT0524
                                                                        ABJT0525
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT0526
c ................................................................      ABJT0527
                                                                        ABJT0528
      common / gases / ngases,izgas(mxgass),atomwt(mxgass),             ABJT0529
     2                 avgatw,avgatn,                                   ABJT0530
     3                 fracsp(mxgass),fraciz(mxgass,mxatom),            ABJT0531
     4                 fraclv(mxgass,mxatom,20),pot(mxgass,mxatom),     ABJT0532
     5                 numocc(mxgass,mxatom)                            ABJT0533
c ................................................................      ABJT0534
                                                                        ABJT0535
      common / params / npqmax(mxgass), npmaxp, nfrqbb,                 ABJT0536
     2                  npring(100), neopen(5,100), bbnton(6,55),       ABJT0537
     3                  defpot(27,55), noccdf(27,55)                    ABJT0538
c ................................................................      ABJT0539
                                                                        ABJT0540
      dimension photen(nphot),abscfs(nphot),sctcfs(nphot)               ABJT0541
      dimension brems(mxgass),fotiza(mxgass),fotize(mxgass)             ABJT0542
      dimension emscfs(nphot)                                           ABJT0543
      dimension nscrsh(6)                                               ABJT0544
                                                                        ABJT0545
      logical lbug                                                      ABJT0546
                                                                        ABJT0547
      data nscrsh / 0,2,10,28,60,110 /                                  ABJT0548
                                                                        ABJT0549
c ********************************************************************  ABJT0550
c                                                                       ABJT0551
c                            begin execution                            ABJT0552
c                                                                       ABJT0553
c ********************************************************************  ABJT0554
                                                                        ABJT0555
c ... check debug option                                                ABJT0556
      call debug ( 'abscon',                                            ABJT0557
     &                       lbug )                                     ABJT0558
                                                                        ABJT0559
      if ( lbug ) write (6,900) tp,densnn,densne,nphot                  ABJT0560
                                                                        ABJT0561
                                                                        ABJT0562
c ... loop over photon energies                                         ABJT0563
                                                                        ABJT0564
      do 500 iphot=1,nphot                                              ABJT0565
                                                                        ABJT0566
         exhvot = exp( -photen(iphot) / tp )                            ABJT0567
         photn3 = photen(iphot)**3                                      ABJT0568
                                                                        ABJT0569
c ...    loop over gas species                                          ABJT0570
                                                                        ABJT0571
         do 300 lgas=1,ngases                                           ABJT0572
                                                                        ABJT0573
            brems(lgas) = 0.                                            ABJT0574
            fotiza(lgas) = 0.                                           ABJT0575
            fotize(lgas) = 0.                                           ABJT0576
                                                                        ABJT0577
c ...       Bremsstrahlung                                              ABJT0578
c           --------------                                              ABJT0579
                                                                        ABJT0580
c ...       loop over ionization levels                                 ABJT0581
            sum = 0.                                                    ABJT0582
            do 100 izlevl=1,izgas(lgas)                                 ABJT0583
               izp1 = izlevl + 1                                        ABJT0584
               sum = sum + izlevl**2 * fraciz(lgas,izp1)                ABJT0585
  100       continue                                                    ABJT0586
                                                                        ABJT0587
            densi2 = sum * fracsp(lgas) * densnn * 1.e-16               ABJT0588
                                                                        ABJT0589
c ...       the free-free gaunt factor is a simple fit to the results   ABJT0590
c           of Karzas and Latter (Ap. J. Suppl., 6, 167 (1961))         ABJT0591
            gam2lg = log10( 13.6 * sum / tp )                           ABJT0592
            gbrem  = 1. + 0.44 * exp( -0.25*(gam2lg+0.25)**2 )          ABJT0593
                                                                        ABJT0594
            brems(lgas) = 2.4e-21 * densi2 * gbrem * densne *           ABJT0595
     &          (1.-exhvot) / ( sqrt( tp ) * photn3 )                   ABJT0596
                                                                        ABJT0597
            if ( lbug ) write (6,902) lgas,brems(lgas),densi2,sum,      ABJT0598
     &                  gbrem                                           ABJT0599
                                                                        ABJT0600
c ...       photoionization                                             ABJT0601
c           ----------------                                            ABJT0602
                                                                        ABJT0603
c ...       sum over ionization levels for the transition from          ABJT0604
c           "izlevl" to "izlevl+1"                                      ABJT0605
                                                                        ABJT0606
            sum1a = 0.                                                  ABJT0607
            sum1e = 0.                                                  ABJT0608
            conpi = 1.66e-22 * densne / tp**1.5                         ABJT0609
                                                                        ABJT0610
            do 200 izlevl=0,izgas(lgas)-1                               ABJT0611
                                                                        ABJT0612
               izp1 = izlevl + 1                                        ABJT0613
                                                                        ABJT0614
c ...          find the principal quantum number of the valence electronABJT0615
c              in the ground state for the ion before ("nprin0");       ABJT0616
c              "nbound" is the number of electrons bound to the ion     ABJT0617
c              after another is captured (or, before ionization).       ABJT0618
                                                                        ABJT0619
               nbound = izgas(lgas) - izlevl                            ABJT0620
               nprin0 = npring(nbound)                                  ABJT0621
                                                                        ABJT0622
c ...          first, consider the contibution from valence shell       ABJT0623
c              electrons                                                ABJT0624
                                                                        ABJT0625
               sum2a = 0.                                               ABJT0626
               sum2e = 0.                                               ABJT0627
               if ( fracsp(lgas)*fraciz(lgas,izp1).ge.con(3) ) then     ABJT0628
                                                                        ABJT0629
c ...            sum over quantum states                                ABJT0630
                 do 150 nprin=nprin0,npqmax(lgas)                       ABJT0631
                                                                        ABJT0632
c ...              calculate the energy to excite the electron into     ABJT0633
c                  the continuum                                        ABJT0634
                   ennp = nprin0**2 * pot(lgas,izp1) / nprin**2         ABJT0635
                                                                        ABJT0636
c ...              the photon energy must exceed the binding energy     ABJT0637
                   if ( photen(iphot) .lt. ennp ) go to 150             ABJT0638
                                                                        ABJT0639
c ...              find the number of "screening" electrons "nscren"    ABJT0640
c                  and the number of electrons in the outermost shell   ABJT0641
c                  "nvalen"                                             ABJT0642
                   if ( nprin.eq.nprin0 ) then                          ABJT0643
c ...                ground state ion                                   ABJT0644
                     nscren = nscrsh( nprin )                           ABJT0645
                     nvalen = nbound - nscren                           ABJT0646
                   else                                                 ABJT0647
c ...                ion is excited                                     ABJT0648
                     nscren = nbound - 1                                ABJT0649
                     nvalen = 1                                         ABJT0650
                   endif                                                ABJT0651
                                                                        ABJT0652
                   ennpot = ennp / tp                                   ABJT0653
                   n5 = nprin**5                                        ABJT0654
                                                                        ABJT0655
c ...              use an "effective charge" seen by the electron,      ABJT0656
c                  corrected for screening                              ABJT0657
                   izeff = izgas(lgas) - nscren                         ABJT0658
                                                                        ABJT0659
                   deni = fraciz(lgas,izp1)*fraclv(lgas,izp1,nprin)     ABJT0660
                   if ( ennpot .lt. 50. ) then                          ABJT0661
                      eqdeni = nprin**2 * conpi * exp( ennpot ) *       ABJT0662
     &                  fraciz(lgas,izp1+1)*fraclv(lgas,izp1+1,nprin0)  ABJT0663
                   else                                                 ABJT0664
                      eqdeni = 0.                                       ABJT0665
                   endif                                                ABJT0666
                                                                        ABJT0667
c ...              the degeneracy level of the fully stripped ion       ABJT0668
c                  is 1, while a value of 2 is used for other ions.     ABJT0669
                   if ( izp1 .eq. izgas(lgas) ) eqdeni = eqdeni * 2.    ABJT0670
                                                                        ABJT0671
c ...              correction for stimulated emission; do not allow     ABJT0672
c                  this to be < 0 for the abs. coef.                    ABJT0673
                   if ( isw(6) .ne. 1 ) then                            ABJT0674
c ...                 non-LTE correction                                ABJT0675
                      dumden = deni                                     ABJT0676
                   else                                                 ABJT0677
c ...                 LTE correction                                    ABJT0678
                      dumden = eqdeni                                   ABJT0679
                   endif                                                ABJT0680
                   stcorr = max ( 0., dumden - eqdeni * exhvot )        ABJT0681
                                                                        ABJT0682
                   xsec = nvalen * izeff**4 / n5                        ABJT0683
                   sum2a = sum2a + stcorr * xsec                        ABJT0684
                   sum2e = sum2e + eqdeni*(1.-exhvot) * xsec            ABJT0685
                                                                        ABJT0686
  150            continue                                               ABJT0687
                                                                        ABJT0688
               endif                                                    ABJT0689
                                                                        ABJT0690
c ...          now, add core electron photoionization cross-sections    ABJT0691
c              to the absorption term                                   ABJT0692
                                                                        ABJT0693
c ...          loop over inner shells (K,L,M,...); each inner shell     ABJT0694
c              is assume to be full                                     ABJT0695
                                                                        ABJT0696
               sum1ac = 0.                                              ABJT0697
               nshels = npring(nbound)-1                                ABJT0698
               if ( nshels .gt. 0 ) then                                ABJT0699
                 do 175 ishell=1,nshels                                 ABJT0700
                                                                        ABJT0701
c ...              determine the photoionization cutoff energy (eV); useABJT0702
c                  the ionization potential of the outermost bound      ABJT0703
c                  electron; "nocc" is the number of electrons occupyingABJT0704
c                  shell "ishell", "nscren" is the number of electrons  ABJT0705
c                  screening shell "ishell"                             ABJT0706
                                                                        ABJT0707
                   nocc = 2*ishell*ishell                               ABJT0708
                   nscren = nscrsh( ishell )                            ABJT0709
                   izeff = izgas(lgas) - nscren                         ABJT0710
                   enpi = pot(lgas,izgas(lgas)+1-nscrsh(ishell+1))      ABJT0711
                   if ( photen(iphot) .ge. enpi ) then                  ABJT0712
                     sum1ac = sum1ac + nocc * izeff**4 / ishell**5      ABJT0713
                   endif                                                ABJT0714
                                                                        ABJT0715
  175            continue                                               ABJT0716
                                                                        ABJT0717
                 sum1a = sum1a + sum1ac * fraciz(lgas,izp1)             ABJT0718
                 if ( lbug ) write (6,915) izp1,nbound,sum1ac,sum2a     ABJT0719
               endif                                                    ABJT0720
                                                                        ABJT0721
               sum1a = sum1a + sum2a                                    ABJT0722
               sum1e = sum1e + sum2e                                    ABJT0723
                                                                        ABJT0724
  200       continue                                                    ABJT0725
                                                                        ABJT0726
            const = (1.99e-14*densnn) * fracsp(lgas) / photn3           ABJT0727
            fotiza(lgas) = const * sum1a                                ABJT0728
            fotize(lgas) = const * sum1e                                ABJT0729
                                                                        ABJT0730
            if ( lbug ) then                                            ABJT0731
               write (6,903) const,fotiza(lgas),sum1a,fotize(lgas),sum1eABJT0732
            endif                                                       ABJT0733
                                                                        ABJT0734
  300    continue                                                       ABJT0735
                                                                        ABJT0736
c ...    Scattering contributions                                       ABJT0737
c        ------------------------                                       ABJT0738
                                                                        ABJT0739
c ...    Thomson scattering contribution                                ABJT0740
c        first, find the "effective" electron density for scattering;   ABJT0741
c        i.e., if the photon energy is greater than the binding energy  ABJT0742
c        of bound electrons, include bound electrons to the density.    ABJT0743
         scatne = 0.                                                    ABJT0744
         do 330 lgas=1,ngases                                           ABJT0745
            iq = 0                                                      ABJT0746
            do 325 izlevl=0,izgas(lgas)-1                               ABJT0747
               izp1 = izlevl + 1                                        ABJT0748
               if ( photen(iphot) .gt. pot(lgas,izp1) ) then            ABJT0749
                  iq = iq + 1                                           ABJT0750
                  scatne = scatne + fracsp(lgas)                        ABJT0751
               endif                                                    ABJT0752
  325       continue                                                    ABJT0753
                                                                        ABJT0754
            do 326 izlevl=1,izgas(lgas)                                 ABJT0755
               izp1 = izlevl + 1                                        ABJT0756
               if ( photen(iphot) .le. pot(lgas,izlevl) ) then          ABJT0757
                  scatne = scatne + (izlevl-iq) *                       ABJT0758
     &                     fraciz(lgas,izp1) * fracsp(lgas)             ABJT0759
               endif                                                    ABJT0760
  326       continue                                                    ABJT0761
  330    continue                                                       ABJT0762
         scatne = scatne * densnn                                       ABJT0763
         tscatt = 6.66e-25 * scatne                                     ABJT0764
                                                                        ABJT0765
c ...    contribution from plasma oscillations                          ABJT0766
         hnucut = sqrt( densne / 7.25e20 )                              ABJT0767
         if ( photen(iphot) .lt. hnucut ) then                          ABJT0768
            pscatt = 5.05e4 * sqrt( hnucut**2 - photen(iphot)**2 )      ABJT0769
         else                                                           ABJT0770
            pscatt = 0.                                                 ABJT0771
         endif                                                          ABJT0772
                                                                        ABJT0773
         if ( isw(20).eq.0 ) sctcfs(iphot) = tscatt + pscatt            ABJT0774
                                                                        ABJT0775
         if ( lbug ) write (6,931) densne,scatne,tscatt,hnucut,pscatt   ABJT0776
                                                                        ABJT0777
         abscfs(iphot) = 0.                                             ABJT0778
         emscfs(iphot) = 0.                                             ABJT0779
         do 400 lgas=1,ngases                                           ABJT0780
           if (isw(17).eq.0) abscfs(iphot) = abscfs(iphot)+brems(lgas)  ABJT0781
           if (isw(18).eq.0) abscfs(iphot) = abscfs(iphot)+fotiza(lgas) ABJT0782
           if (isw(17).eq.0) emscfs(iphot) = emscfs(iphot)+brems(lgas)  ABJT0783
           if (isw(18).eq.0) emscfs(iphot) = emscfs(iphot)+fotize(lgas) ABJT0784
  400    continue                                                       ABJT0785
                                                                        ABJT0786
c ...    line contributions                                             ABJT0787
c        ------------------                                             ABJT0788
                                                                        ABJT0789
c ...    add in the contribution from bound-bound transitions           ABJT0790
         if ( isw(19).eq.0 .or. isw(19).eq.2 ) then                     ABJT0791
            call lines ( tp,densnn,densne,photen(iphot),                ABJT0792
     &                                             abslns,emslns )      ABJT0793
         else                                                           ABJT0794
            abslns = 0.                                                 ABJT0795
            emslns = 0.                                                 ABJT0796
         endif                                                          ABJT0797
                                                                        ABJT0798
         abscfs(iphot) = abscfs(iphot) + abslns                         ABJT0799
         emscfs(iphot) = emscfs(iphot) + emslns                         ABJT0800
                                                                        ABJT0801
         if ( lbug ) then                                               ABJT0802
            write (6,906) iphot,photen(iphot),abscfs(iphot),            ABJT0803
     &                    sctcfs(iphot),abslns                          ABJT0804
            write (6,908) emscfs(iphot),emslns                          ABJT0805
         endif                                                          ABJT0806
                                                                        ABJT0807
c ...    write to plot file                                             ABJT0808
         rho = densnn * avgatw / avgdro                                 ABJT0809
         if ( iplot(1).eq.1 ) then                                      ABJT0810
            brmtot = 0.                                                 ABJT0811
            piztot = 0.                                                 ABJT0812
            do 460 lgas=1,ngases                                        ABJT0813
               brmtot = brmtot + brems(lgas)                            ABJT0814
               piztot = piztot + fotiza(lgas)                           ABJT0815
  460       continue                                                    ABJT0816
            write (11,911) photen(iphot),abscfs(iphot),                 ABJT0817
     &         brmtot,max(1.e-30,piztot),max(1.e-30,abslns),            ABJT0818
     &         sctcfs(iphot),abscfs(iphot)+sctcfs(iphot)                ABJT0819
         endif                                                          ABJT0820
                                                                        ABJT0821
         if ( iplot(3).eq.1 ) then                                      ABJT0822
            brmtot = 0.                                                 ABJT0823
            piztot = 0.                                                 ABJT0824
            do 470 lgas=1,ngases                                        ABJT0825
               brmtot = brmtot + brems(lgas)                            ABJT0826
               piztot = piztot + fotize(lgas)                           ABJT0827
  470       continue                                                    ABJT0828
            write (13,911) photen(iphot),emscfs(iphot),                 ABJT0829
     &         brmtot,max(1.e-30,piztot),max(1.e-10,emslns)             ABJT0830
         endif                                                          ABJT0831
                                                                        ABJT0832
  500 continue                                                          ABJT0833
                                                                        ABJT0834
                                                                        ABJT0835
      return                                                            ABJT0836
                                                                        ABJT0837
c ... format statements                                                 ABJT0838
                                                                        ABJT0839
  900 format (' debug output from -abscon-:'/t4,'tp',t18,'densnn',      ABJT0840
     2  t32,'densne',t46,'nphot'/t4,1p3e14.4,0p,i4)                     ABJT0841
  902 format (t4,'lgas',t18,'brems',t32,'densi2',t46,'sum',             ABJT0842
     2  t60,'gbrem'/t4,i4,t18,1p6e14.4)                                 ABJT0843
  903 format (t4,'const',t18,'fotiza',t32,'sum1a',t46,'fotize',         ABJT0844
     2  t60,'sum1e'/t4,1p5e14.4)                                        ABJT0845
  906 format (t4,'iphot',t18,'photen',t32,'abscfs',t46,'sctcfs',        ABJT0846
     2  t60,'abslns'/t4,i4,t18,1p4e14.4)                                ABJT0847
  908 format (t4,'emscfs',t18,'emslns'/                                 ABJT0848
     2  t4,1p2e14.4)                                                    ABJT0849
  911 format (1p7e14.5)                                                 ABJT0850
  915 format (t4,'izp1,nbound,sum1ac,sum2a =',2i5,1p2e14.3)             ABJT0851
  931 format (t4,'densne',t18,'scatne',t32,'tscatt',t46,'hnucut',       ABJT0852
     2  t60,'pscatt'/t4,1p5e14.4)                                       ABJT0853
                                                                        ABJT0854
      end                                                               ABJT0855
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT0856
                                                                        ABJT0857
      subroutine abslin ( tp,densnn,densne,denlqn,dnlqnp,ephot,potiz,   ABJT0858
     &                    atomwt,n,np,nprin0,nbound,izgas,              ABJT0859
     &                                                 abscof,emscof )  ABJT0860
                                                                        ABJT0861
c ... computes the absorption coefficient from a particular             ABJT0862
c     bound-bound transition from quantum state "n" to "np".            ABJT0863
c                                                                       ABJT0864
c ... input variables                                                   ABJT0865
c       tp      =  plasma temperature (eV)                              ABJT0866
c       densnn  =  number density of all nuclei (cm**-3)                ABJT0867
c       densne  =  electron density (cm**-3)                            ABJT0868
c       denlqn  =  number density of nuclei of the "l"th gas, in        ABJT0869
c                  the "q"th ionization state, with an electron in      ABJT0870
c                  the "n"th quantum state                              ABJT0871
c       dnlqnp  =  number density of nuclei of the "l"th gas, in        ABJT0872
c                  the "q"th ionization state, with an electron in      ABJT0873
c                  the "np"th quantum state                             ABJT0874
c       ephot   =  photon energy (eV)                                   ABJT0875
c       potiz   =  ionization potential for the ion in its ground       ABJT0876
c                  state (eV)                                           ABJT0877
c       atomwt  =  atomic weight of the ion (amu)                       ABJT0878
c       n       =  initial principal quantum number for transition      ABJT0879
c       np      =  final principal quantum number                       ABJT0880
c       nprin0  =  principal quantum number of the valence electrons    ABJT0881
c                  in their ground state                                ABJT0882
c       nbound  =  number of electrons bound to the ion                 ABJT0883
c       izgas   =  atomic number of the ion                             ABJT0884
c                                                                       ABJT0885
c ... output variable                                                   ABJT0886
c       abscof  =  absorption coefficient (cm**-1)                      ABJT0887
c       emscof  =  emission coefficient (cm**-1)                        ABJT0888
c                                                                       ABJT0889
                                                                        ABJT0890
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT0891
c ................................................................      ABJT0892
                                                                        ABJT0893
      common / consts / pi, avgdro, sbcon, hplank                       ABJT0894
c ...............................................................       ABJT0895
                                                                        ABJT0896
c ... declaration statements                                            ABJT0897
      logical lbug                                                      ABJT0898
                                                                        ABJT0899
c ... statement function (oscillator strength taken from Zeldovich      ABJT0900
c     & Raizor for the upward transition from "ni" to "nf")             ABJT0901
      oscstr(ni,nf) = 1.96 / ni**5 / nf**3 / (1./ni**2-1./nf**2)**3     ABJT0902
                                                                        ABJT0903
c ********************************************************************  ABJT0904
c                                                                       ABJT0905
c                             begin execution                           ABJT0906
c                                                                       ABJT0907
c ********************************************************************  ABJT0908
                                                                        ABJT0909
c ... check debug option                                                ABJT0910
      call debug ( 'abslin',                                            ABJT0911
     &                        lbug )                                    ABJT0912
                                                                        ABJT0913
      if ( lbug ) write (6,900) n,np,denlqn,potiz,ephot,tp,densnn       ABJT0914
                                                                        ABJT0915
c ... compute the transition energy                                     ABJT0916
      if ( n.eq.np ) then                                               ABJT0917
         call bbneq0 ( tp,izgas,nbound,                                 ABJT0918
     &                                   ennp,fnn,gnn )                 ABJT0919
c ...    some transitions are not allowed                               ABJT0920
         if ( ennp.le.0. .or. isw(11).ne.0 ) return                     ABJT0921
      else                                                              ABJT0922
         ennp = nprin0**2 * potiz * (1./n**2-1./np**2)                  ABJT0923
      endif                                                             ABJT0924
                                                                        ABJT0925
c ... compute the line widths for natural, Doppler, and pressure        ABJT0926
c     broadening to be used with Lorentzian line shape                  ABJT0927
      call linwid ( tp,densnn,ennp,atomwt,                              ABJT0928
     &                                     gamma,avoigt,dnudop )        ABJT0929
                                                                        ABJT0930
c ... compute this contribution if the photon energy is not far         ABJT0931
c     from the line center                                              ABJT0932
      deltav = ( ennp-ephot ) / hplank                                  ABJT0933
      if ( abs( deltav ) .le. con(5)*gamma ) then                       ABJT0934
                                                                        ABJT0935
c ...    if the contribution from line "cores" will be added            ABJT0936
c        elsewhere (-opacbb-), then use the value at the                ABJT0937
c        line core boundary.                                            ABJT0938
         if ( isw(19).eq.2 .and. abs(deltav).lt.con(6)*gamma ) then     ABJT0939
            deltav = con(6)*gamma                                       ABJT0940
         endif                                                          ABJT0941
                                                                        ABJT0942
         if ( isw(14).eq.0 ) then                                       ABJT0943
c ...       use Voigt profile                                           ABJT0944
            vvoigt = deltav / dnudop                                    ABJT0945
            shape = voigt ( avoigt,vvoigt ) / dnudop / 1.7725           ABJT0946
         else                                                           ABJT0947
c ...       use Lorentzian profile                                      ABJT0948
            shape = (gamma/39.48) / (deltav**2 + (gamma/12.57)**2)      ABJT0949
         endif                                                          ABJT0950
                                                                        ABJT0951
         ex1 = exp( -ephot/tp )                                         ABJT0952
         ex2 = exp( - ennp/tp )                                         ABJT0953
                                                                        ABJT0954
         if ( n.eq.np ) then                                            ABJT0955
                                                                        ABJT0956
            alpha = 2.65e-2 * fnn * shape * denlqn * ( 1.-ex1 )         ABJT0957
                                                                        ABJT0958
c ...       correct for the "effective" stimulated emission to          ABJT0959
c           give the proper form for the cooling rate                   ABJT0960
            if ( isw(6) .ne. 1 ) then                                   ABJT0961
c ...          general form                                             ABJT0962
               xnjoni = gnn * densne * ex2 /                            ABJT0963
     &                ( gnn * densne + ennp**3 * sqrt(tp) * 2.74e12 )   ABJT0964
               corsee = gnn * densne /                                  ABJT0965
     &                ( gnn * densne + ennp**3 * sqrt(tp) * 2.74e12 )   ABJT0966
            else                                                        ABJT0967
c ...          LTE assumed                                              ABJT0968
               xnjoni = ex2                                             ABJT0969
               corsee = 1.                                              ABJT0970
            endif                                                       ABJT0971
                                                                        ABJT0972
            corsea = ( 1.-xnjoni ) / ( 1.-ex2 )                         ABJT0973
            abscof = alpha * corsea                                     ABJT0974
            emscof = alpha * corsee                                     ABJT0975
                                                                        ABJT0976
         else                                                           ABJT0977
                                                                        ABJT0978
            fnnp = oscstr ( n,np )                                      ABJT0979
            dum = ( n*n*dnlqnp ) / ( np*np*denlqn )                     ABJT0980
                                                                        ABJT0981
c ...       correct for stimulated emission                             ABJT0982
            corrse = ( 1.-dum ) * ( 1.-ex1 ) / ( 1.-ex2 )               ABJT0983
            alpha  = 2.65e-2 * fnnp * shape * denlqn                    ABJT0984
            abscof = alpha * corrse                                     ABJT0985
            if ( ex2 .gt. 0. ) then                                     ABJT0986
              emscof = alpha * dum * ( 1.-ex1 ) / ex2                   ABJT0987
            else                                                        ABJT0988
              emscof = abscof                                           ABJT0989
            endif                                                       ABJT0990
                                                                        ABJT0991
         endif                                                          ABJT0992
                                                                        ABJT0993
      else                                                              ABJT0994
                                                                        ABJT0995
         abscof = 0.                                                    ABJT0996
         emscof = 0.                                                    ABJT0997
                                                                        ABJT0998
      endif                                                             ABJT0999
                                                                        ABJT1000
      if ( lbug ) then                                                  ABJT1001
         write (6,906) abscof,fnnp,deltav,shape,gamma,ennp              ABJT1002
         write (6,907) emscof,corrse,dnlqnp                             ABJT1003
      endif                                                             ABJT1004
                                                                        ABJT1005
      return                                                            ABJT1006
                                                                        ABJT1007
c ... format statements                                                 ABJT1008
                                                                        ABJT1009
  900 format (t2,'debug output from -abslin-'/t4,'n',t18,'np',          ABJT1010
     2  t32,'denlqn',t46,'potiz',t60,'ephot',t74,'tp',t88,'densnn'/     ABJT1011
     3  t4,i5,t18,i5,t32,1p5e12.3)                                      ABJT1012
  906 format (t4,'abscof',t18,'fnnp',t32,'deltav',t46,'shape',          ABJT1013
     2  t60,'gamma',t74,'ennp'/t4,1p6e14.4)                             ABJT1014
  907 format (t4,'emscof',t18,'corrse',t32,'dnlqnp'/t4,1p3e14.4)        ABJT1015
                                                                        ABJT1016
      end                                                               ABJT1017
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT1018
                                                                        ABJT1019
      subroutine atmlv ( lgas,izlevl,tp,izgas,pot,densne,               ABJT1020
     &                                                     fraclv )     ABJT1021
                                                                        ABJT1022
c ... this routine calculates the fractional densities of the           ABJT1023
c     various atomic levels -- normalized to the fraction in            ABJT1024
c     ionization state "izlevl" -- for ionic species "lgas" at          ABJT1025
c     a temperature(gas) "tp".                                          ABJT1026
c                                                                       ABJT1027
c ... input variables:                                                  ABJT1028
c       lgas    -  gas species index                                    ABJT1029
c       izlevl  -  ionization level (0 => neutral)                      ABJT1030
c       tp      -  plasma temperature (eV)                              ABJT1031
c       izgas   -  array of atomic numbers of all gas species           ABJT1032
c       pot     -  ionization potentials for all gas species and        ABJT1033
c                  all ionization states (eV)                           ABJT1034
c       densne  -  electron density (cm**-3)                            ABJT1035
c                                                                       ABJT1036
c ... output variables:                                                 ABJT1037
c       fraclv  -  fraction of a given ionization level in a particular ABJT1038
c                  principal quantum state                              ABJT1039
c                                                                       ABJT1040
                                                                        ABJT1041
c ... set the maximum number of:                                        ABJT1042
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT1043
c       temperatures (mxtemp), densities (mxdens),                      ABJT1044
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT1045
                                                                        ABJT1046
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT1047
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT1048
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT1049
                                                                        ABJT1050
                                                                        ABJT1051
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT1052
c ................................................................      ABJT1053
                                                                        ABJT1054
      common / params / npqmax(mxgass), npmaxp, nfrqbb,                 ABJT1055
     2                  npring(100), neopen(5,100), bbnton(6,55),       ABJT1056
     3                  defpot(27,55), noccdf(27,55)                    ABJT1057
c ................................................................      ABJT1058
                                                                        ABJT1059
c ... type declaration statements                                       ABJT1060
      dimension pot(mxgass,mxatom),fraclv(mxgass,mxatom,20)             ABJT1061
      dimension izgas(mxgass)                                           ABJT1062
                                                                        ABJT1063
      logical lbug                                                      ABJT1064
                                                                        ABJT1065
c ********************************************************************  ABJT1066
c                                                                       ABJT1067
c                         begin execution                               ABJT1068
c                                                                       ABJT1069
c ********************************************************************  ABJT1070
                                                                        ABJT1071
c ... check debug option                                                ABJT1072
      call debug ( 'atmlv ',                                            ABJT1073
     &                       lbug )                                     ABJT1074
                                                                        ABJT1075
      if ( lbug ) write (6,901) lgas,izlevl,npmaxp,tp,densne            ABJT1076
                                                                        ABJT1077
      izp1 = izlevl + 1                                                 ABJT1078
                                                                        ABJT1079
                                                                        ABJT1080
c ... determine the principal quantum # of the valence electrons        ABJT1081
c     in their ground state                                             ABJT1082
      nbound = izgas(lgas) - izlevl                                     ABJT1083
      if ( nbound .gt. 0 ) then                                         ABJT1084
         nprin0 = npring(nbound)                                        ABJT1085
      else                                                              ABJT1086
         nprin0 = 1                                                     ABJT1087
      endif                                                             ABJT1088
      fraclv(lgas,izp1,nprin0) = 1.                                     ABJT1089
                                                                        ABJT1090
c ... set the populations with lower principal quantum numbers to zero  ABJT1091
      if ( nprin0 .gt. 1 ) then                                         ABJT1092
         do 3 n=1,nprin0-1                                              ABJT1093
           fraclv(lgas,izp1,n) = 0.                                     ABJT1094
    3    continue                                                       ABJT1095
      endif                                                             ABJT1096
                                                                        ABJT1097
c ... if the ionization level "izlevl" corresponds to full ionization,  ABJT1098
c     the electrons can be thought of as all being in the ground state  ABJT1099
c     within this ionization level; also, inputting a nonzero value     ABJT1100
c     for switch 12 will automatically restrict ions to their ground    ABJT1101
c     states.                                                           ABJT1102
                                                                        ABJT1103
      if ( izlevl .eq. izgas(lgas) .or. isw(12).ne.0 ) then             ABJT1104
         do 5 n=nprin0+1,npmaxp                                         ABJT1105
            fraclv(lgas,izp1,n) = 0.                                    ABJT1106
    5    continue                                                       ABJT1107
         go to 90                                                       ABJT1108
      endif                                                             ABJT1109
                                                                        ABJT1110
                                                                        ABJT1111
c ... now, for the ionization stages below full ionization              ABJT1112
                                                                        ABJT1113
      sum = 1.                                                          ABJT1114
                                                                        ABJT1115
c ... set a lower limit for the electron density used to avoid          ABJT1116
c     overflows                                                         ABJT1117
      elden = max( densne , 1.e0 )                                      ABJT1118
                                                                        ABJT1119
c ... the relative densities of atoms with valence electrons            ABJT1120
c     in each of the first "npmaxp" atomic states are equal to their    ABJT1121
c     Boltzmann factors times the degeneracies of the levels.           ABJT1122
c     the second term in "denom" includes the radiative decay rate,     ABJT1123
c     which dominates when the electron density is small.               ABJT1124
                                                                        ABJT1125
      do 10 nprin=nprin0+1,npmaxp                                       ABJT1126
         xnm12 = (nprin-1)**2                                           ABJT1127
         xn2   = nprin * nprin                                          ABJT1128
         ennp  = nprin0**2 * pot(lgas,izp1) * (1./xnm12-1./xn2 )        ABJT1129
         xnnp = ennp / tp                                               ABJT1130
                                                                        ABJT1131
c ...    compute the Gaunt factor; this is based on Van Regimorter      ABJT1132
c        (1962); (see also Tucker & Gould (1967).                       ABJT1133
         if ( xnnp .ge. 1. ) then                                       ABJT1134
            gauntf = 0.2                                                ABJT1135
         else if ( xnnp .le. 0.01 ) then                                ABJT1136
            gauntf = 2.76 * ( -.577-log(xnnp) )                         ABJT1137
         else                                                           ABJT1138
            yy = 1. / sqrt( xnnp )                                      ABJT1139
            gauntf = 0.2 + 0.15 * (yy-1)                                ABJT1140
         endif                                                          ABJT1141
                                                                        ABJT1142
c ...    if "isw6".eq.1, then LTE populations are calculated            ABJT1143
         if ( isw(6) .ne. 1 ) then                                      ABJT1144
            denom = xnm12/xn2 * ( 1. + ennp**3 * sqrt(tp) *             ABJT1145
     &                          ( 2.74e12/elden/gauntf ) )              ABJT1146
         else                                                           ABJT1147
            denom = xnm12 / xn2                                         ABJT1148
         endif                                                          ABJT1149
                                                                        ABJT1150
         fraclv(lgas,izp1,nprin) = fraclv(lgas,izp1,nprin-1) *          ABJT1151
     &          exp( -xnnp ) / denom                                    ABJT1152
         sum = sum + fraclv(lgas,izp1,nprin)                            ABJT1153
   10 continue                                                          ABJT1154
                                                                        ABJT1155
c ... the density of atoms must by normalized to conserve particles     ABJT1156
                                                                        ABJT1157
      do 20 nprin=nprin0,npmaxp                                         ABJT1158
        fraclv(lgas,izp1,nprin) = fraclv(lgas,izp1,nprin) / sum         ABJT1159
   20 continue                                                          ABJT1160
                                                                        ABJT1161
   90 continue                                                          ABJT1162
      if ( lbug ) write (6,906) (fraclv(lgas,izp1,n),n=1,npmaxp)        ABJT1163
                                                                        ABJT1164
      return                                                            ABJT1165
                                                                        ABJT1166
c ... format statements                                                 ABJT1167
                                                                        ABJT1168
  901 format (' debug output from -atmlv-:'/t4,'lgas',t18,              ABJT1169
     2  'izlevl',t32,'npmaxp',t46,'tp',t60,'densne'/                    ABJT1170
     3  t4,3(i4,10x),1p2e14.4)                                          ABJT1171
  906 format ('  fraclv:'/3(t4,1p7e11.3/))                              ABJT1172
                                                                        ABJT1173
      end                                                               ABJT1174
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT1175
c MOVED TO ModMultiGroup.f90                                            ABJT1176
      subroutine bbneq0 ( tp,izgas,nbound,                              ABJT1177
     &                                      enn,fnn,gnn )               ABJT1178
                                                                        ABJT1179
c ... calculates the effective energy, oscillator strength, and         ABJT1180
c     gaunt factor for bound-bound transitions in which the             ABJT1181
c     principal quantum number does not change.  This approximation     ABJT1182
c     was taken from Post, et.al., Atomic & Nuclear Data Tables,        ABJT1183
c     vol. 20, p.397 (1977).                                            ABJT1184
c                                                                       ABJT1185
c ... input variables                                                   ABJT1186
c       tp      -  plasma temperature (eV)                              ABJT1187
c       izgas   -  nuclear charge of the ion                            ABJT1188
c       nbound  -  number of electrons bound to the ion                 ABJT1189
c                                                                       ABJT1190
c ... output variable                                                   ABJT1191
c       enn  -  effective energy of the transition (eV)                 ABJT1192
c       fnn  -  effective oscillator strength                           ABJT1193
c       gnn  -  effective gaunt factor                                  ABJT1194
c                                                                       ABJT1195
                                                                        ABJT1196
c ... set the maximum number of:                                        ABJT1197
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT1198
c       temperatures (mxtemp), densities (mxdens),                      ABJT1199
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT1200
                                                                        ABJT1201
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT1202
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT1203
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT1204
                                                                        ABJT1205
                                                                        ABJT1206
      common / params / npqmax(mxgass), npmaxp, nfrqbb,                 ABJT1207
     2                  npring(100), neopen(5,100), bbnton(6,55),       ABJT1208
     3                  defpot(27,55), noccdf(27,55)                    ABJT1209
c ................................................................      ABJT1210
                                                                        ABJT1211
c ... declaration statements                                            ABJT1212
      logical lbug                                                      ABJT1213
                                                                        ABJT1214
c ********************************************************************  ABJT1215
c                                                                       ABJT1216
c                             begin execution                           ABJT1217
c                                                                       ABJT1218
c ********************************************************************  ABJT1219
                                                                        ABJT1220
c ... check debug option                                                ABJT1221
      call debug ( 'bbneq0',                                            ABJT1222
     &                        lbug )                                    ABJT1223
                                                                        ABJT1224
      iq =izgas + 1 - nbound                                            ABJT1225
      nb = min( nbound , 55 )                                           ABJT1226
                                                                        ABJT1227
      fnn = bbnton(1,nb) + bbnton(2,nb) / izgas                         ABJT1228
      enn = bbnton(3,nb) * iq**bbnton(4,nb)                             ABJT1229
                                                                        ABJT1230
      if ( enn .le. 0. ) return                                         ABJT1231
                                                                        ABJT1232
      xnn = enn / tp                                                    ABJT1233
      if ( iq.eq.1 ) then                                               ABJT1234
         c = 0.06 * (sqrt( xnn )-2.) / (1.+xnn)                         ABJT1235
      else                                                              ABJT1236
         c = bbnton(5,nb) * ( 1. - bbnton(6,nb)/(iq-1) )                ABJT1237
      endif                                                             ABJT1238
                                                                        ABJT1239
c ... compute the exponential integral (see Abramowitz & Stegun), and   ABJT1240
c     multiply by exp(xnn)                                              ABJT1241
      if ( xnn .ge. 1. ) then                                           ABJT1242
         expint = (xnn**2+2.3347*xnn+0.25062) /                         ABJT1243
     &            (xnn**2+3.3306*xnn+1.6815) / xnn                      ABJT1244
      else                                                              ABJT1245
         expint = -log( xnn ) - 0.57722 + 0.99999*xnn - 0.24991*xnn**2  ABJT1246
     &            + 0.05520*xnn**3 - 0.00976*xnn**4 + 0.00108*xnn**5    ABJT1247
         expint = expint * exp( xnn )                                   ABJT1248
      endif                                                             ABJT1249
                                                                        ABJT1250
      gnn = c + 0.276 * expint                                          ABJT1251
                                                                        ABJT1252
      if ( lbug ) write (6,901) izgas,iq,nbound,enn,fnn,gnn             ABJT1253
                                                                        ABJT1254
      return                                                            ABJT1255
                                                                        ABJT1256
c ... format statements                                                 ABJT1257
                                                                        ABJT1258
  901 format (t2,'debug output from -bbneq0-'/t4,'izgas',t12,           ABJT1259
     2  'iq',t20,'nbound',t32,'enn',t46,'fnn',t60,'gnn'/                ABJT1260
     3  t4,i4,t12,i4,t20,i4,t32,1p3e14.4)                               ABJT1261
c END OF THE MOVED FRAGMENT                                             ABJT1262
      end                                                               ABJT1263
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT1264
                                                                        ABJT1265
      function colliz ( q,n,nprin0,tempel,atomnm,potiz )                ABJT1266
c NON-LTE PART NOT NEEDED                                               ABJT1267
c ... collisional ionization rate function.  calculates the             ABJT1268
c     rate (per ion per free electron => cm**3/sec) at which            ABJT1269
c     electrons leave the "n"th level of an ion with charge state       ABJT1270
c     "q".  "q+1" is the resulting ionization state of the atom.        ABJT1271
c                                                                       ABJT1272
c ... input variables                                                   ABJT1273
c       q       =  charge state of the ion before ionization            ABJT1274
c                  (0 < q < atomnm-1)                                   ABJT1275
c       n       =  principal quantum # of the electron being ionized    ABJT1276
c       nprin0  =  principal quantum # of the valence electrons in      ABJT1277
c                  their ground state of this ion                       ABJT1278
c       tempel  =  electron temperature (eV)                            ABJT1279
c       atomnm  =  atomic number of the ion                             ABJT1280
c       potiz   =  ionization potential for an unexcited ion with       ABJT1281
c                  charge state q (eV)                                  ABJT1282
c                                                                       ABJT1283
c ... output variable                                                   ABJT1284
c       colliz  =  collisional ionization rate (cm**3/sec)              ABJT1285
c                                                                       ABJT1286
c ... declaration statements                                            ABJT1287
      logical lbug                                                      ABJT1288
                                                                        ABJT1289
c ********************************************************************  ABJT1290
c                                                                       ABJT1291
c                             begin execution                           ABJT1292
c                                                                       ABJT1293
c ********************************************************************  ABJT1294
                                                                        ABJT1295
c ... check debug option                                                ABJT1296
      call debug ( 'colliz',                                            ABJT1297
     &                        lbug )                                    ABJT1298
                                                                        ABJT1299
      if ( lbug ) write (6,900) n,q,tempel,potniz                       ABJT1300
                                                                        ABJT1301
      tkev = tempel / 1000.                                             ABJT1302
                                                                        ABJT1303
c ... ionization potential based on Bohr model                          ABJT1304
      potniz = nprin0**2 * potiz / n**2                                 ABJT1305
                                                                        ABJT1306
      xn = potniz / tempel                                              ABJT1307
      emxn = exp( -xn )                                                 ABJT1308
                                                                        ABJT1309
c ... compute Gaunt factor                                              ABJT1310
      if ( q .lt. 1. ) then                                             ABJT1311
         gauntf = 1.98                                                  ABJT1312
      else                                                              ABJT1313
         y = log10( xn ) * 0.25                                         ABJT1314
         fofy = 0.23 + 0.046*y + 0.1074*y**2 - 0.0459*y**3              ABJT1315
     &          - 0.01505*y**4                                          ABJT1316
         arg = n / (n+5)                                                ABJT1317
         gauntf = 12.18 * exp( -arg ) *                                 ABJT1318
     &        (1.-0.622/q-0.0745/q/q) * (1.+0.0335*n) * abs( fofy )     ABJT1319
      endif                                                             ABJT1320
                                                                        ABJT1321
c ... compute ionization rate                                           ABJT1322
      colliz = 1.86e-7 * sqrt( tkev ) * emxn * (1.-emxn) * gauntf       ABJT1323
     &         * (13.6/potniz)**2                                       ABJT1324
                                                                        ABJT1325
      if ( lbug ) write (6,906) xn,emxn,gauntf,colliz                   ABJT1326
                                                                        ABJT1327
      return                                                            ABJT1328
                                                                        ABJT1329
c ... format statements                                                 ABJT1330
                                                                        ABJT1331
  900 format (t2,'debug output from -colliz-'/t4,'n',t18,'q',           ABJT1332
     2  t32,'tempel',t46,'potiz'/t4,i5,t18,1p3e12.3)                    ABJT1333
  906 format (t4,'xn',t18,'emxn',t32,'gauntf',t46,'colliz'/t4,          ABJT1334
     2  1p4e12.3)                                                       ABJT1335
                                                                        ABJT1336
      end                                                               ABJT1337
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT1338
                                                                        ABJT1339
      function collrc ( q,n,nprin0,tempel,atomnm,potiz,edens0,pfrat )   ABJT1340
c NON-LTE PART -NOT NEEDED                                              ABJT1341
c ... collisional recombination rate function.  calculates the          ABJT1342
c     rate (per ion per free electron => cm**3/sec) at which            ABJT1343
c     electrons enter the "n"th level of an ion with charge state       ABJT1344
c     "q".  "q+1" is the initial ionization state of the atom.          ABJT1345
c                                                                       ABJT1346
c ... input variables                                                   ABJT1347
c       q       =  charge state of the ion before ionization            ABJT1348
c                  (0 < q < atomnm-1)                                   ABJT1349
c       n       =  principal quantum # of the electron after            ABJT1350
c                  recombining                                          ABJT1351
c       nprin0  =  principal quantum # of the valence electrons in      ABJT1352
c                  their ground state of this ion                       ABJT1353
c       tempel  =  electron temperature (eV)                            ABJT1354
c       atomnm  =  atomic number of the ion                             ABJT1355
c       potiz   =  ionization potential for an unexcited ion with       ABJT1356
c                  charge state q (eV)                                  ABJT1357
c       edens0  =  initial guess at electron density (cm**-3)           ABJT1358
c       pfrat   =  ratio of the electronic partition functions (q/q+1)  ABJT1359
c                                                                       ABJT1360
c ... output variable                                                   ABJT1361
c       collrc  =  collisional recombination rate (cm**3/sec)           ABJT1362
c                                                                       ABJT1363
c ... declaration statements                                            ABJT1364
      logical lbug                                                      ABJT1365
                                                                        ABJT1366
c ********************************************************************  ABJT1367
c                                                                       ABJT1368
c                             begin execution                           ABJT1369
c                                                                       ABJT1370
c ********************************************************************  ABJT1371
                                                                        ABJT1372
c ... check debug option                                                ABJT1373
      call debug ( 'collrc',                                            ABJT1374
     &                        lbug )                                    ABJT1375
                                                                        ABJT1376
      if ( lbug ) write (6,900) n,q,tempel,potniz                       ABJT1377
                                                                        ABJT1378
      tkev = tempel / 1000.                                             ABJT1379
                                                                        ABJT1380
c ... ionization potential based on Bohr model                          ABJT1381
      potniz = nprin0**2 * potiz / n**2                                 ABJT1382
                                                                        ABJT1383
      xn = potniz / tempel                                              ABJT1384
      emxn = exp( -xn )                                                 ABJT1385
                                                                        ABJT1386
c ... compute Gaunt factor                                              ABJT1387
      if ( q .lt. 1. ) then                                             ABJT1388
         gauntf = 1.98                                                  ABJT1389
      else                                                              ABJT1390
         y = log10( xn ) * 0.25                                         ABJT1391
         fofy = 0.23 + 0.046*y + 0.1074*y**2 - 0.0459*y**3              ABJT1392
     &          - 0.01505*y**4                                          ABJT1393
         arg = n / (n+5)                                                ABJT1394
         gauntf = 12.18 * exp( -arg ) *                                 ABJT1395
     &        (1.-0.622/q-0.0745/q/q) * (1.+0.0335*n) * abs( fofy )     ABJT1396
      endif                                                             ABJT1397
                                                                        ABJT1398
c ... compute recombination rate                                        ABJT1399
      collrc = 1.86e-7  * (1.-emxn) * gauntf * (13.6/potniz)**2         ABJT1400
     &       * 5.25e-27 * edens0 / tkev * pfrat                         ABJT1401
                                                                        ABJT1402
      if ( lbug ) write (6,906) xn,emxn,gauntf,collrc,pfrat             ABJT1403
                                                                        ABJT1404
      return                                                            ABJT1405
                                                                        ABJT1406
c ... format statements                                                 ABJT1407
                                                                        ABJT1408
  900 format (t2,'debug output from -collrc-'/t4,'n',t18,'q',           ABJT1409
     2  t32,'tempel',t46,'potiz'/t4,i5,t18,1p3e12.3)                    ABJT1410
  906 format (t4,'xn',t18,'emxn',t32,'gauntf',t46,'collrc',             ABJT1411
     2  t60,'pfrat'/t4,1p5e12.3)                                        ABJT1412
                                                                        ABJT1413
      end                                                               ABJT1414
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT1415
                                                                        ABJT1416
      subroutine corona ( tp,densnn,ngases,izgas,pot,numocc,fracsp,     ABJT1417
     &                    edens0,                                       ABJT1418
     &                                                 densne,fraciz )  ABJT1419
c NON-LTE  PART - NOT NEEDED                                            ABJT1420
c ... calculate the ionization populations according to the             ABJT1421
c     coronal model                                                     ABJT1422
c                                                                       ABJT1423
c ... input variables                                                   ABJT1424
c       tp      -  plasma temperature (eV)                              ABJT1425
c       densnn  -  ion density (cm**-3)                                 ABJT1426
c       ngases  -  number of gas species                                ABJT1427
c       izgas   -  atomic number of each gas                            ABJT1428
c       pot(l,j)-  ionization potential for the ground state of "j-1"   ABJT1429
c                  charge state of the "l"th gas (eV)                   ABJT1430
c       numocc  -  number of electrons in the valence shell of an ion   ABJT1431
c                  in its ground state                                  ABJT1432
c       fracsp  -  gas species number fraction                          ABJT1433
c       edens0  -  initial guess at the electron density (cm**-3)       ABJT1434
c                                                                       ABJT1435
c ... output variables                                                  ABJT1436
c       densne      -  electron density (cm**-3)                        ABJT1437
c       fraciz(l,j) -  fraction of ions of gas "l" in the "j-1" charge  ABJT1438
c                      state; i.e., fraciz(l,1) => neutral atom.        ABJT1439
c                                                                       ABJT1440
                                                                        ABJT1441
c ... set the maximum number of:                                        ABJT1442
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT1443
c       temperatures (mxtemp), densities (mxdens),                      ABJT1444
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT1445
                                                                        ABJT1446
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT1447
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT1448
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT1449
                                                                        ABJT1450
                                                                        ABJT1451
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT1452
c ................................................................      ABJT1453
                                                                        ABJT1454
      common / params / npqmax(mxgass), npmaxp, nfrqbb,                 ABJT1455
     2                  npring(100), neopen(5,100), bbnton(6,55),       ABJT1456
     3                  defpot(27,55), noccdf(27,55)                    ABJT1457
c ................................................................      ABJT1458
c                                                                       ABJT1459
c ... type declaration statements                                       ABJT1460
      dimension izgas(mxgass),dumm(mxatom),fracsp(mxgass)               ABJT1461
      dimension pot(mxgass,mxatom),numocc(mxgass,mxatom)                ABJT1462
      dimension fraciz(mxgass,mxatom), partfn(mxgass,mxatom)            ABJT1463
                                                                        ABJT1464
      logical lbug                                                      ABJT1465
                                                                        ABJT1466
c ********************************************************************  ABJT1467
c                                                                       ABJT1468
c                          begin execution                              ABJT1469
c                                                                       ABJT1470
c ********************************************************************  ABJT1471
                                                                        ABJT1472
c ... check debug option                                                ABJT1473
      call debug ( 'corona',                                            ABJT1474
     &                       lbug )                                     ABJT1475
                                                                        ABJT1476
      if ( lbug ) then                                                  ABJT1477
         write (6,901) ngases,edens0,tp,izgas(1),izgas(2)               ABJT1478
      endif                                                             ABJT1479
                                                                        ABJT1480
                                                                        ABJT1481
c ... first, calculate the electronic partition functions if 3-body     ABJT1482
c     recombination is to be included                                   ABJT1483
                                                                        ABJT1484
      if ( isw(6) .eq. 3 ) then                                         ABJT1485
                                                                        ABJT1486
         do 5 lgas=1,ngases                                             ABJT1487
                                                                        ABJT1488
            izgm1 = izgas(lgas)-1                                       ABJT1489
            do 4 izlevl=0,izgm1                                         ABJT1490
                                                                        ABJT1491
c ...          find the principal quantum number of the valence         ABJT1492
c              electrons in their ground state                          ABJT1493
               nbound = izgas(lgas) - izlevl                            ABJT1494
               nprin0 = npring(nbound)                                  ABJT1495
                                                                        ABJT1496
               izp1 = izlevl + 1                                        ABJT1497
               eizokt = pot(lgas,izp1) / tp                             ABJT1498
                                                                        ABJT1499
c ...          compute the electronic partition function                ABJT1500
                                                                        ABJT1501
               qsum = 0.                                                ABJT1502
               do 3 n=nprin0,npmaxp                                     ABJT1503
c ...             calculate the energy relative to the ground state     ABJT1504
                  ennpot = nprin0**2 * eizokt * (1./nprin0**2-1./n**2)  ABJT1505
                  qsum = qsum + 2.*n**2 * exp( -ennpot )                ABJT1506
    3          continue                                                 ABJT1507
               partfn(lgas,izp1) = qsum                                 ABJT1508
                                                                        ABJT1509
    4       continue                                                    ABJT1510
                                                                        ABJT1511
c ...       the last partition fn corresponds to the fully ionized stateABJT1512
            partfn(lgas,izgas(lgas)+1) = 1.                             ABJT1513
                                                                        ABJT1514
    5    continue                                                       ABJT1515
                                                                        ABJT1516
      endif                                                             ABJT1517
                                                                        ABJT1518
c ... initialize variables                                              ABJT1519
      densne = 0.                                                       ABJT1520
                                                                        ABJT1521
c ... sum over gas species                                              ABJT1522
                                                                        ABJT1523
      do 100 lgas = 1,ngases                                            ABJT1524
                                                                        ABJT1525
c ...    find relative populations over ionization states               ABJT1526
                                                                        ABJT1527
         dumm(1) = 0.                                                   ABJT1528
         dummax = 0.                                                    ABJT1529
         do 10 izlevl=1,izgas(lgas)                                     ABJT1530
            izp1 = izlevl + 1                                           ABJT1531
                                                                        ABJT1532
c ...       find the ratio of the ionization rate to recombination      ABJT1533
c           rate between charge states "izlevl-1" and "izlevl".         ABJT1534
c ...       assume the electron and plasma temperatures are equil       ABJT1535
                                                                        ABJT1536
            q = izlevl - 1                                              ABJT1537
            atomnm = izgas(lgas)                                        ABJT1538
            if ( isw(6) .eq. 3 ) then                                   ABJT1539
               pfrat = partfn(lgas,izlevl) / partfn(lgas,izp1)          ABJT1540
            else                                                        ABJT1541
               pfrat = 1.                                               ABJT1542
            endif                                                       ABJT1543
                                                                        ABJT1544
            call rizrec ( tp,densnn,q,atomnm,edens0,pfrat,              ABJT1545
     &                    pot(lgas,izlevl),numocc(lgas,izlevl),         ABJT1546
     &                                                         ratio )  ABJT1547
                                                                        ABJT1548
            dumm(izp1) = dumm(izlevl) + log( max( ratio,1.e-30 ) )      ABJT1549
            dummax = max( dummax , dumm(izp1) )                         ABJT1550
   10    continue                                                       ABJT1551
                                                                        ABJT1552
c ...    this "re-normalization" procedure for the relative ionization  ABJT1553
c        populations was inserted to stop overflow problems             ABJT1554
                                                                        ABJT1555
         sumfrc = 0.                                                    ABJT1556
         do 20 izlevl=0,izgas(lgas)                                     ABJT1557
            izp1 = izlevl + 1                                           ABJT1558
            fraciz(lgas,izp1) = exp( dumm(izp1) - dummax )              ABJT1559
            sumfrc = sumfrc + fraciz(lgas,izp1)                         ABJT1560
 20      continue                                                       ABJT1561
                                                                        ABJT1562
c ...    normalize the populations (relative to the total               ABJT1563
c        number of that species)                                        ABJT1564
c        sum the electron density to find average ionization            ABJT1565
c        state of species "lgas"                                        ABJT1566
                                                                        ABJT1567
         sumne = 0.                                                     ABJT1568
         do 30 izlevl=0,izgas(lgas)                                     ABJT1569
            izp1 = izlevl + 1                                           ABJT1570
            fraciz(lgas,izp1) = fraciz(lgas,izp1) / sumfrc              ABJT1571
            sumne = sumne + fraciz(lgas,izp1) * izlevl                  ABJT1572
   30    continue                                                       ABJT1573
                                                                        ABJT1574
c ...    sum for the total electron density                             ABJT1575
         densne = densne + sumne*fracsp(lgas)                           ABJT1576
                                                                        ABJT1577
  100 continue                                                          ABJT1578
                                                                        ABJT1579
      densne = densne * densnn                                          ABJT1580
                                                                        ABJT1581
      if ( lbug ) then                                                  ABJT1582
         write (6,906)  densne                                          ABJT1583
         do 510 lgas=1,ngases                                           ABJT1584
            write (6,907)  ( fraciz(lgas,izp1),izp1=1,izgas(lgas)+1 )   ABJT1585
  510    continue                                                       ABJT1586
      endif                                                             ABJT1587
                                                                        ABJT1588
      return                                                            ABJT1589
                                                                        ABJT1590
c ... format statements                                                 ABJT1591
                                                                        ABJT1592
  901 format (' debug output from -corona-:'/t4,'ngases',t18,           ABJT1593
     2  'edens0',t32,'tp',t46,'izgas(1)',t60,'izgas(2)'/                ABJT1594
     &  t4,i4,t18,1p2e14.4,0p,t46,i4,t60,i4)                            ABJT1595
  906 format (t4,'densne'/t4,1pe14.4)                                   ABJT1596
  907 format (t2,'fraciz:'/8(t4,1p7e11.3/))                             ABJT1597
                                                                        ABJT1598
      end                                                               ABJT1599
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT1600
c MOVED TO ModVoigt                                                     ABJT1601
      function dawson ( v )                                             ABJT1602
                                                                        ABJT1603
c ... computes dawson's integral (which is sometimes used to compute    ABJT1604
c     the Voigt function):                                              ABJT1605
c                               _v                                      ABJT1606
c       D.I.  =  exp(-v**2) * _/ dt exp(t**2)                           ABJT1607
c                            0                                          ABJT1608
c                                                                       ABJT1609
      dimension table1(21),table2(11)                                   ABJT1610
                                                                        ABJT1611
      data table1 / .00000, .09933, .19475, .28263, .35994,             ABJT1612
     2      .42443, .47476, .51050, .53210, .54072, .53807,             ABJT1613
     3      .52620, .50727, .48339, .45650, .42824, .39993,             ABJT1614
     4      .37255, .34677, .32297, .30134 /                            ABJT1615
      data table2 / .50000, .50650, .51358, .52142, .53037,             ABJT1616
     2      .54079, .55265, .56547, .57852, .59110, .60268 /            ABJT1617
                                                                        ABJT1618
c ********************************************************************  ABJT1619
c                                                                       ABJT1620
c                           begin execution                             ABJT1621
c                                                                       ABJT1622
c ********************************************************************  ABJT1623
                                                                        ABJT1624
      if ( v .le. 2. ) then                                             ABJT1625
         x = v                                                          ABJT1626
         ix = 1 + x / 0.1                                               ABJT1627
         f = x/0.1 - (ix-1)                                             ABJT1628
         y = (1.-f) * table1(ix) + f * table1(ix+1)                     ABJT1629
         dawson = y                                                     ABJT1630
      else                                                              ABJT1631
         x = 1. / (v*v)                                                 ABJT1632
         ix = 1 + x / 0.025                                             ABJT1633
         f = x/0.025 - (ix-1)                                           ABJT1634
         y = (1.-f) * table2(ix) + f * table2(ix+1)                     ABJT1635
         dawson = y / v                                                 ABJT1636
      endif                                                             ABJT1637
                                                                        ABJT1638
      return                                                            ABJT1639
      end
c END OF MOVED DAWSON                                                   ABJT1640
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT1641
c NOT NEEDED                                                            ABJT1642
      subroutine debug ( name,                                          ABJT1643
     &                          lbug )                                  ABJT1644
                                                                        ABJT1645
c ... this routine sets the logical flag "lbug" to true if debug        ABJT1646
c     output has been requested for subroutine "name".                  ABJT1647
c                                                                       ABJT1648
c ... input:                                                            ABJT1649
c        name  =  name of the subroutine (character string)             ABJT1650
c                                                                       ABJT1651
c ... output:                                                           ABJT1652
c        lbug  =  logical flag; set to .true. if debug output is        ABJT1653
c                 requested                                             ABJT1654
c                                                                       ABJT1655
                                                                        ABJT1656
      common / dbugcm / nobug,ncycld(10),icycld(10)                     ABJT1657
                                                                        ABJT1658
      logical nobug                                                     ABJT1659
c ................................................................      ABJT1660
                                                                        ABJT1661
      common / strngs / header, namedb(10)                              ABJT1662
                                                                        ABJT1663
      character header*80, namedb*6                                     ABJT1664
c ................................................................      ABJT1665
c                                                                       ABJT1666
c ... declaration statements:                                           ABJT1667
      character*6 name                                                  ABJT1668
                                                                        ABJT1669
      logical lbug                                                      ABJT1670
                                                                        ABJT1671
c ********************************************************************  ABJT1672
c                                                                       ABJT1673
c                            begin execution                            ABJT1674
c                                                                       ABJT1675
c ********************************************************************  ABJT1676
                                                                        ABJT1677
      if ( nobug ) return                                               ABJT1678
                                                                        ABJT1679
      lbug = .false.                                                    ABJT1680
                                                                        ABJT1681
      do 100 i=1,10                                                     ABJT1682
                                                                        ABJT1683
         if ( namedb(i).eq.'      ' ) return                            ABJT1684
                                                                        ABJT1685
         if ( namedb(i).eq.name ) then                                  ABJT1686
                                                                        ABJT1687
            icycld(i) = icycld(i) + 1                                   ABJT1688
                                                                        ABJT1689
c ...       set the debug flag to .true. only if we are on the          ABJT1690
c           right cycle number                                          ABJT1691
            iovern = icycld(i) / ncycld(i)                              ABJT1692
            if ( iovern*ncycld(i) .eq. icycld(i) ) then                 ABJT1693
c ...          do not allow over 50 debug prints per subroutine         ABJT1694
               if ( iovern.le.50 ) lbug = .true.                        ABJT1695
               return                                                   ABJT1696
            else                                                        ABJT1697
               return                                                   ABJT1698
            endif                                                       ABJT1699
                                                                        ABJT1700
         endif                                                          ABJT1701
                                                                        ABJT1702
  100 continue                                                          ABJT1703
                                                                        ABJT1704
      return                                                            ABJT1705
                                                                        ABJT1706
      end                                                               ABJT1707
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT1708
c END NOT NEEDED                                                        ABJT1709
      function dielrc ( q,n,nprin0,nbound,izgas,tempel,edens0,potiz )   ABJT1710
c NON-LTE PART 
                                                                        ABJT1711
c ... dielectronic recombination rate function.  calculates the         ABJT1712
c     rate (per ion per free electron => cm**3/sec) at which            ABJT1713
c     electrons enter the "n"th level of an ion with charge state       ABJT1714
c     "q".  "q+1" is the original ionization state of the atom.         ABJT1715
c                                                                       ABJT1716
c ... input variables                                                   ABJT1717
c       q       =  charge state of the ion after recombination          ABJT1718
c                  (0 < q < atomnm-1)                                   ABJT1719
c       n       =  principal quantum # of the electron after being      ABJT1720
c                  recombined                                           ABJT1721
c       nprin0  =  principal quantum # of the valence electrons in      ABJT1722
c                  the ion's ground state                               ABJT1723
c       nbound  =  number of bound electrons of the ion before          ABJT1724
c                  recombination                                        ABJT1725
c       izgas   =  atomic number of the ion                             ABJT1726
c       tempel  =  electron temperature (eV)                            ABJT1727
c       edens0  =  electron density (cm**-3)                            ABJT1728
c       potiz   =  ionization potential (eV)                            ABJT1729
c                                                                       ABJT1730
c ... output variable                                                   ABJT1731
c       dielrc  =  dielectronic recombination rate (cm**3/sec)          ABJT1732
c                                                                       ABJT1733
c ... declaration statements                                            ABJT1734
      logical lbug                                                      ABJT1735
                                                                        ABJT1736
c ... statement function (oscillator strength taken from Zeldovich      ABJT1737
c     & Raizor for the upward transition from "ni" to "nf")             ABJT1738
      oscstr(ni,nf) = 1.96 / ni**5 / nf**3 / (1./ni**2-1./nf**2)**3     ABJT1739
                                                                        ABJT1740
c ********************************************************************  ABJT1741
c                                                                       ABJT1742
c                             begin execution                           ABJT1743
c                                                                       ABJT1744
c ********************************************************************  ABJT1745
                                                                        ABJT1746
c ... check debug option                                                ABJT1747
      call debug ( 'dielrc',                                            ABJT1748
     &                        lbug )                                    ABJT1749
                                                                        ABJT1750
      tkev = tempel / 1000.                                             ABJT1751
      qp1 = q + 1.                                                      ABJT1752
      qp2 = q + 2.                                                      ABJT1753
                                                                        ABJT1754
      bq = sqrt( qp1 * qp2**5 / (qp1**2+13.4) )                         ABJT1755
      if ( edens0 .gt. 1.e0 ) then                                      ABJT1756
        xnt = ( 4.77e18 * qp1**6 * sqrt( tkev ) / edens0 ) ** (1./7.)   ABJT1757
      else                                                              ABJT1758
        xnt = 20.                                                       ABJT1759
      endif                                                             ABJT1760
      nt = min( xnt , 20. )                                             ABJT1761
      a = 1. + 0.015 * qp1**3 / qp2**2                                  ABJT1762
      ebarcn = qp2**2 / ( 73.5*tkev*a )                                 ABJT1763
      dqtnn = 0.005*xnt / ( 1.+0.005*xnt )                              ABJT1764
      dum = 0.0015 * (qp2*xnt)**2                                       ABJT1765
      dqtnnp = dum / (1.+dum)                                           ABJT1766
                                                                        ABJT1767
      if ( lbug ) write (6,900) n,q,tempel,edens0,xnt,nt                ABJT1768
                                                                        ABJT1769
      summ = 0.                                                         ABJT1770
      do 100 np=n,nt                                                    ABJT1771
                                                                        ABJT1772
         if ( np.eq.n ) then                                            ABJT1773
                                                                        ABJT1774
            call bbneq0 ( tempel,izgas,nbound,                          ABJT1775
     &                                         enn,fnnp,dum1 )          ABJT1776
c ...       convert to keV                                              ABJT1777
            eij = enn / 1000.                                           ABJT1778
c ...       convert to units of 13.6 eV (changed 11/87)                 ABJT1779
            eij = eij / 0.0136                                          ABJT1780
                                                                        ABJT1781
            y = qp2 * eij                                               ABJT1782
            dqt = dqtnn                                                 ABJT1783
            ay = sqrt( y ) / ( 1.+0.105*y+0.015*y**2 )                  ABJT1784
                                                                        ABJT1785
         else                                                           ABJT1786
                                                                        ABJT1787
            fnnp = oscstr(n,np)                                         ABJT1788
            eij = nprin0**2 * ( 1./n**2 - 1./np**2 )                    ABJT1789
            y = qp2 * eij                                               ABJT1790
            dqt = dqtnnp                                                ABJT1791
            ay = 0.5 * sqrt( y ) / (1.+0.210*y+0.030*y**2)              ABJT1792
                                                                        ABJT1793
         endif                                                          ABJT1794
                                                                        ABJT1795
         ebarot = eij * ebarcn                                          ABJT1796
         expon = exp( -ebarot )                                         ABJT1797
         summ = summ + dqt * fnnp * ay * expon                          ABJT1798
         if ( n.eq.np ) summ1 = dqt*fnnp*ay*expon                       ABJT1799
                                                                        ABJT1800
  100 continue                                                          ABJT1801
                                                                        ABJT1802
                                                                        ABJT1803
c ... compute recombination rate                                        ABJT1804
      dielrc = 7.59e-14 * bq * summ / tkev**1.5                         ABJT1805
                                                                        ABJT1806
      if ( lbug ) write (6,905) bq,summ,dielrc,enn,summ1                ABJT1807
                                                                        ABJT1808
      return                                                            ABJT1809
                                                                        ABJT1810
c ... format statements                                                 ABJT1811
                                                                        ABJT1812
  900 format (t2,'debug output from -dielrc-'/t4,'n',t18,'q',           ABJT1813
     2  t32,'tempel',t46,'edens0',t60,'xnt',t74,'nt'/t4,i5,             ABJT1814
     3  t18,1p4e12.3,0p,t74,i5)                                         ABJT1815
  905 format (t4,'bq',t18,'summ',t32,'dielrc',t46,'enn',t60,            ABJT1816
     2  'summ1'/t4,1p5e12.3)                                            ABJT1817
                                                                        ABJT1818
      end
c END OF NON-LTE PART                                                   ABJT1819
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT1820
                                                                        ABJT1821
      subroutine energy ( tp,densnn,                                    ABJT1822
     &                                enrgy,densne )                    ABJT1823
c  NON-LTE PART                                                         ABJT1824
c ... this subroutine computes the ionization states of gases as a      ABJT1825
c     function of density and temperature using either the saha         ABJT1826
c     or coronal ionization models                                      ABJT1827
c     after the populations have been computed, the energy density      ABJT1828
c     and electron density are calculated                               ABJT1829
c                                                                       ABJT1830
c ... input variables:                                                  ABJT1831
c       tp      -  plasma temperature (eV)                              ABJT1832
c       densnn  -  total density of nuclei (cm**-3)                     ABJT1833
c                                                                       ABJT1834
c ... output variables:                                                 ABJT1835
c       enrgy   -  specific energy (J/g)                                ABJT1836
c       densne  -  electron density (cm**-3)                            ABJT1837
c                                                                       ABJT1838
                                                                        ABJT1839
c ... set the maximum number of:                                        ABJT1840
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT1841
c       temperatures (mxtemp), densities (mxdens),                      ABJT1842
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT1843
                                                                        ABJT1844
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT1845
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT1846
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT1847
                                                                        ABJT1848
                                                                        ABJT1849
      common / consts / pi, avgdro, sbcon, hplank                       ABJT1850
c ...............................................................       ABJT1851
                                                                        ABJT1852
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT1853
c ................................................................      ABJT1854
                                                                        ABJT1855
      common / gases / ngases,izgas(mxgass),atomwt(mxgass),             ABJT1856
     2                 avgatw,avgatn,                                   ABJT1857
     3                 fracsp(mxgass),fraciz(mxgass,mxatom),            ABJT1858
     4                 fraclv(mxgass,mxatom,20),pot(mxgass,mxatom),     ABJT1859
     5                 numocc(mxgass,mxatom)                            ABJT1860
c ................................................................      ABJT1861
                                                                        ABJT1862
      common / params / npqmax(mxgass), npmaxp, nfrqbb,                 ABJT1863
     2                  npring(100), neopen(5,100), bbnton(6,55),       ABJT1864
     3                  defpot(27,55), noccdf(27,55)                    ABJT1865
c ................................................................      ABJT1866
                                                                        ABJT1867
c ... type declaration statements                                       ABJT1868
      dimension fracc(mxgass,mxatom),fracs(mxgass,mxatom)               ABJT1869
                                                                        ABJT1870
      logical lbug                                                      ABJT1871
                                                                        ABJT1872
c ********************************************************************  ABJT1873
c                                                                       ABJT1874
c                           begin execution                             ABJT1875
c                                                                       ABJT1876
c ********************************************************************  ABJT1877
                                                                        ABJT1878
c ... check debug option                                                ABJT1879
      call debug ( 'energy',                                            ABJT1880
     &                        lbug )                                    ABJT1881
                                                                        ABJT1882
      if ( lbug ) write (6,901) tp,densnn,ngases,isw(6)                 ABJT1883
                                                                        ABJT1884
c ... set the maximum principal quantum number used in the summations.  ABJT1885
c     "npmaxp" is used for computing the populations; "npqmax" is       ABJT1886
c     used for calculation absorption coefficients and opacities        ABJT1887
                                                                        ABJT1888
c ... max. quantum # based on mean atomic separation                    ABJT1889
      npmaxp = sqrt( izgas(1) / .5292e-8 / densnn**.33 )                ABJT1890
      npmaxp = min( npmaxp , 20 )                                       ABJT1891
                                                                        ABJT1892
c ... if "isw(15)" is zero, use the default value, which is 1 plus the  ABJT1893
c     principal quantum number of the neutral atom's outer electrons,   ABJT1894
c     thereby insuring at least the first excited state of each ion is  ABJT1895
c     considered in the absorption coefficient calculations.            ABJT1896
                                                                        ABJT1897
      maxnpq = 0                                                        ABJT1898
      do 5 lgas=1,ngases                                                ABJT1899
         if ( isw(15) .ge. 0 ) then                                     ABJT1900
            npqmax(lgas) = npring( izgas(lgas) ) + max( 1, isw(15) )    ABJT1901
         else                                                           ABJT1902
            npqmax(lgas) = -isw(15)                                     ABJT1903
         endif                                                          ABJT1904
         npqmax(lgas) = min ( 20 , max( 2 , npqmax(lgas) ) )            ABJT1905
         maxnpq = max( maxnpq , npqmax(lgas) )                          ABJT1906
    5 continue                                                          ABJT1907
                                                                        ABJT1908
      npmaxp = max( npmaxp , maxnpq )                                   ABJT1909
      if ( isw(9).ne.0 ) npmaxp = npring( izgas(1) ) + isw(9)           ABJT1910
                                                                        ABJT1911
      if ( lbug ) write (6,903) tp,densnn,(npqmax(lgas),lgas=1,ngases)  ABJT1912
                                                                        ABJT1913
                                                                        ABJT1914
c ... find the populations of the ionization states of each gas         ABJT1915
c     ---------------------------------------------------------         ABJT1916
                                                                        ABJT1917
      if ( isw(6).eq.1 ) then                                           ABJT1918
                                                                        ABJT1919
c ...    Saha model only                                                ABJT1920
         call saha ( tp,densnn,fracsp,ngases,izgas,pot,                 ABJT1921
     &                                                 densne,fraciz )  ABJT1922
                                                                        ABJT1923
      else if ( isw(6).eq.2 ) then                                      ABJT1924
                                                                        ABJT1925
c ...    Coronal model only                                             ABJT1926
         edens0 = densnn * avgatn                                       ABJT1927
         call corona ( tp,densnn,ngases,izgas,pot,numocc,fracsp,        ABJT1928
     &                 edens0,                                          ABJT1929
     &                                                 densne,fraciz )  ABJT1930
                                                                        ABJT1931
      else if ( isw(6).eq.3 ) then                                      ABJT1932
                                                                        ABJT1933
c ...    Coronal model with 3-body recombination; iterate to find the   ABJT1934
c        appropriate electron density                                   ABJT1935
                                                                        ABJT1936
c ...    first guess at the electron density (based on hydrogenic Saha eABJT1937
                                                                        ABJT1938
         avgpot = 0.                                                    ABJT1939
         do 15 lgas=1,ngases                                            ABJT1940
            avgpot = avgpot + pot(lgas,1)*fracsp(lgas)                  ABJT1941
   15    continue                                                       ABJT1942
         beta = 3.e21 / densnn * tp**1.5 * exp( -avgpot/tp )            ABJT1943
         if ( beta .eq. 0. ) then                                       ABJT1944
            edens0 = densnn * 1.e-20                                    ABJT1945
         else if ( beta .lt. 100. ) then                                ABJT1946
            edens0 = avgatn * densnn * beta*0.5*( sqrt(1.+4./beta) -1. )ABJT1947
         else                                                           ABJT1948
            edens0 = avgatn * densnn                                    ABJT1949
         endif                                                          ABJT1950
                                                                        ABJT1951
c ...    iterate on the electron density until converged                ABJT1952
                                                                        ABJT1953
         iter = 0                                                       ABJT1954
   25    iter = iter + 1                                                ABJT1955
                                                                        ABJT1956
            edens = edens0                                              ABJT1957
            call corona ( tp,densnn,ngases,izgas,pot,numocc,fracsp,     ABJT1958
     &                    edens,                                        ABJT1959
     &                                             densne,fraciz )      ABJT1960
            q1 = densne - edens                                         ABJT1961
            edens = edens0 * 1.01                                       ABJT1962
            call corona ( tp,densnn,ngases,izgas,pot,numocc,fracsp,     ABJT1963
     &                    edens,                                        ABJT1964
     &                                             densne,fraciz )      ABJT1965
            q2 = densne - edens                                         ABJT1966
            dqdne = ( q2-q1 ) / ( 0.01*edens0 )                         ABJT1967
            if ( dqdne .ge. 0. ) then                                   ABJT1968
               dne = q1                                                 ABJT1969
            else                                                        ABJT1970
               dne = -q1 / dqdne                                        ABJT1971
            endif                                                       ABJT1972
            if ( abs( dne ) .lt. 1.e-4*edens0 ) then                    ABJT1973
               densne = edens0                                          ABJT1974
            else if ( iter .le. 20 ) then                               ABJT1975
               if ( edens0+dne .le. 0. ) then                           ABJT1976
                 edens0 = edens0 * 1.e-2                                ABJT1977
               else                                                     ABJT1978
                 edens0 = edens0 + dne                                  ABJT1979
               endif                                                    ABJT1980
               edens0 = min( edens0 , avgatn*densnn )                   ABJT1981
               if ( lbug ) write (6,990) iter,edens0,q1,q2,dqdne,       ABJT1982
     &                     dne                                          ABJT1983
               go to 25                                                 ABJT1984
            else                                                        ABJT1985
               write (6,991) iter,dne,edens0,densne                     ABJT1986
               stop ' loop did not converge in -energy-'                ABJT1987
            endif                                                       ABJT1988
                                                                        ABJT1989
      else                                                              ABJT1990
                                                                        ABJT1991
c ...    let the program decide which model to use based on the         ABJT1992
c        criterion of Mosher (NRL report)                               ABJT1993
                                                                        ABJT1994
         edens0 = densnn * avgatn                                       ABJT1995
         call corona ( tp,densnn,ngases,izgas,pot,numocc,fracsp,        ABJT1996
     &                 edens0,                                          ABJT1997
     &                                             densne,fracc )       ABJT1998
                                                                        ABJT1999
         critr1 = densne / ( 1.e16 * tp**3.5 )                          ABJT2000
         if ( critr1 .gt. 1. ) then                                     ABJT2001
c ...       use Saha model                                              ABJT2002
            call saha ( tp,densnn,fracsp,ngases,izgas,pot,              ABJT2003
     &                                                 densne,fracs )   ABJT2004
         endif                                                          ABJT2005
                                                                        ABJT2006
c ...    if the electron density is very low, use the Coronal model;    ABJT2007
c        if the electron density is very high, use the Saha model;      ABJT2008
c        otherwise interpolate between the two                          ABJT2009
                                                                        ABJT2010
         fcrit = ( critr1 - 1. ) / ( critsc - 1. )                      ABJT2011
         fcrit = max( 0. , min( fcrit , 1. ) )                          ABJT2012
                                                                        ABJT2013
c ...    note:  using a simple average like this, the states            ABJT2014
c        remain normalized                                              ABJT2015
                                                                        ABJT2016
         sumne = 0.                                                     ABJT2017
         do 100 lgas=1,ngases                                           ABJT2018
         do 100 izlevl=0,izgas(lgas)                                    ABJT2019
            izp1 = izlevl + 1                                           ABJT2020
            fraciz(lgas,izp1) = (1.-fcrit) * fracc(lgas,izp1) +         ABJT2021
     &                          fcrit * fracs(lgas,izp1)                ABJT2022
            sumne = sumne + fracsp(lgas)*fraciz(lgas,izp1)*izlevl       ABJT2023
  100    continue                                                       ABJT2024
                                                                        ABJT2025
         densne = sumne * densnn                                        ABJT2026
                                                                        ABJT2027
      endif                                                             ABJT2028
                                                                        ABJT2029
                                                                        ABJT2030
c ... find the population of levels within each ionization state        ABJT2031
c     ----------------------------------------------------------        ABJT2032
                                                                        ABJT2033
      do 200 lgas=1,ngases                                              ABJT2034
      do 200 izlevl=0,izgas(lgas)                                       ABJT2035
         call atmlv ( lgas,izlevl,tp,izgas,pot,densne,                  ABJT2036
     &                                                  fraclv )        ABJT2037
  200 continue                                                          ABJT2038
                                                                        ABJT2039
                                                                        ABJT2040
c ... calulate the energy density                                       ABJT2041
c     ---------------------------                                       ABJT2042
                                                                        ABJT2043
c ... calculate the average ionization state                            ABJT2044
      avizst = densne / densnn                                          ABJT2045
                                                                        ABJT2046
c ... calculate the thermal energy per nucleus                          ABJT2047
      etherm = 1.5 * ( 1.+avizst ) * tp                                 ABJT2048
                                                                        ABJT2049
c ... sum for the energy in ionization and excitation                   ABJT2050
                                                                        ABJT2051
      tsumex = 0.                                                       ABJT2052
      tsumiz = 0.                                                       ABJT2053
      do 320 lgas=1,ngases                                              ABJT2054
                                                                        ABJT2055
         sumex = 0.                                                     ABJT2056
         sumiz = 0.                                                     ABJT2057
         pottot = 0.                                                    ABJT2058
         do 310 izlevl=0,izgas(lgas)                                    ABJT2059
                                                                        ABJT2060
            izp1 = izlevl + 1                                           ABJT2061
            if ( izlevl .lt. izgas(lgas) ) then                         ABJT2062
                                                                        ABJT2063
c ...          find the principal quantum number of the valence electronABJT2064
c              in their ground state                                    ABJT2065
               nbound = izgas(lgas) - izlevl                            ABJT2066
               nprin0 = npring(nbound)                                  ABJT2067
                                                                        ABJT2068
               sum = 0.                                                 ABJT2069
               do 300 nprin=nprin0,npmaxp                               ABJT2070
c ...             calculate the energy relative to the ground state     ABJT2071
                  ennp = nprin0**2 * pot(lgas,izp1) *                   ABJT2072
     &                   (1./nprin0**2 - 1./nprin**2)                   ABJT2073
                  sum = sum + fraclv(lgas,izp1,nprin) * ennp            ABJT2074
  300          continue                                                 ABJT2075
               sumex = sumex + sum * fraciz(lgas,izp1)                  ABJT2076
            endif                                                       ABJT2077
            if ( izlevl .ne. 0 ) then                                   ABJT2078
               pottot = pottot + pot(lgas,izlevl)                       ABJT2079
               sumiz = sumiz + fraciz(lgas,izp1) * pottot               ABJT2080
            endif                                                       ABJT2081
                                                                        ABJT2082
  310    continue                                                       ABJT2083
                                                                        ABJT2084
         tsumex = tsumex + sumex * fracsp(lgas)                         ABJT2085
         tsumiz = tsumiz + sumiz * fracsp(lgas)                         ABJT2086
                                                                        ABJT2087
  320 continue                                                          ABJT2088
                                                                        ABJT2089
c ... calculate the energy density                                      ABJT2090
      enrgy = ( etherm + tsumex + tsumiz ) * densnn                     ABJT2091
                                                                        ABJT2092
c ... convert to J/g                                                    ABJT2093
      rho = densnn * avgatw / avgdro                                    ABJT2094
      enrgy = enrgy / ( rho * 6.242e18 )                                ABJT2095
                                                                        ABJT2096
      if ( lbug ) then                                                  ABJT2097
         write (6,906) enrgy,etherm,tsumex,tsumiz                       ABJT2098
      endif                                                             ABJT2099
                                                                        ABJT2100
c ... plot ionization states ?                                          ABJT2101
      if ( iplot(9) .eq. 1 ) then                                       ABJT2102
         do 505 izp1=1,izgas(1)+1                                       ABJT2103
            write (19,909) tp,log10( max( fraciz(1,izp1) , 1.e-30 ) )   ABJT2104
  505    continue                                                       ABJT2105
      endif                                                             ABJT2106
                                                                        ABJT2107
                                                                        ABJT2108
      return                                                            ABJT2109
                                                                        ABJT2110
c ... format statements                                                 ABJT2111
                                                                        ABJT2112
  901 format (' debug output from -energy-:'/t4,'tp',t18,'densnn',      ABJT2113
     2  t32,'ngases',t46,'isw(6)'/t4,1p2e14.4,0p,2(i4,10x))             ABJT2114
  903 format (/t2,'Temperature           =',1pe11.3/                    ABJT2115
     2         t2,'Number Density        =',1pe11.3/                    ABJT2116
     3         t2,'Max. Prin. Quantum #s =',0p,20i5/)                   ABJT2117
  906 format (t4,'enrgy',t18,'etherm',t32,'tsumex',t46,'tsumiz'/        ABJT2118
     2  t4,1p4e14.4)                                                    ABJT2119
  909 format (t2,1p5e11.3)                                              ABJT2120
  990 format (t4,'iter',t18,'edens0',t32,'q1',t46,'q2',t60,'dqdne',     ABJT2121
     2  t74,'dne'/t4,i4,t18,1p5e14.4)                                   ABJT2122
  991 format (///' loop did not converge in -energy-'/                  ABJT2123
     2  t4,'iter',t18,'dne',t32,'edens0',t46,'densne'/                  ABJT2124
     3  t4,i4,t18,1p3e14.4)                                             ABJT2125
                                                                        ABJT2126
      end                                                               ABJT2127
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT2128
c NOT NEEDED: NON-LTE                                                   ABJT2129
      subroutine eos ( tp,densnn,                                       ABJT2130
     &                             enrgy,heatcp,dzdt,densne,dedden )    ABJT2131
c                                                                       ABJT2132
c ... compute the energy density, heat capacity, the temperature        ABJT2133
c     derivative of the energy density, and the electron density        ABJT2134
c                                                                       ABJT2135
c ... input variables:                                                  ABJT2136
c       tp      -  plasma temperature (eV)                              ABJT2137
c       densnn  -  total density of nuclei (cm**-3)                     ABJT2138
c                                                                       ABJT2139
c ... output variables:                                                 ABJT2140
c       enrgy   -  energy density (J/g)                                 ABJT2141
c       heatcp  -  heat capacity (J/g/eV)                               ABJT2142
c       dzdt    -  temperature derivative of the average charge         ABJT2143
c                  state (eV**-1)                                       ABJT2144
c       densne  -  electron density (cm**-3)                            ABJT2145
c       dedden  -  ion density derivative of the specific energy        ABJT2146
c                  (J*cm**3/g)                                          ABJT2147
c                                                                       ABJT2148
                                                                        ABJT2149
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT2150
c ................................................................      ABJT2151
                                                                        ABJT2152
      logical lbug                                                      ABJT2153
                                                                        ABJT2154
c *******************************************************************   ABJT2155
c                                                                       ABJT2156
c                           begin execution                             ABJT2157
c                                                                       ABJT2158
c *******************************************************************   ABJT2159
                                                                        ABJT2160
c ... check debug option                                                ABJT2161
      call debug ( 'eos   ',                                            ABJT2162
     &                       lbug )                                     ABJT2163
                                                                        ABJT2164
c ... if isw4=0, calculate the heat capacity, etc.                      ABJT2165
                                                                        ABJT2166
      if( isw(4) .eq. 0 ) then                                          ABJT2167
                                                                        ABJT2168
c ...    find energy densities and electron densities at tp*(1+dtheat)  ABJT2169
c        for calculating temperature derivatives                        ABJT2170
         call energy ( tp*(1.+dtheat),densnn,                           ABJT2171
     &                                        energ2,eld2 )             ABJT2172
                                                                        ABJT2173
c ...    find the energy density at densnn*(1+dtheat) to calculate the  ABJT2174
c        density derivative                                             ABJT2175
         call energy ( tp,densnn*(1.+dtheat),                           ABJT2176
     &                                        energ3,eld3 )             ABJT2177
                                                                        ABJT2178
      endif                                                             ABJT2179
                                                                        ABJT2180
c ... calculate the energy density and electron density at T,rho;       ABJT2181
c     this must be calculated last because the ionization populations   ABJT2182
c     are used to compute the absorption coefficients                   ABJT2183
                                                                        ABJT2184
      call energy ( tp,densnn,                                          ABJT2185
     &                         enrgy,densne )                           ABJT2186
                                                                        ABJT2187
                                                                        ABJT2188
      if( isw(4) .eq. 0 ) then                                          ABJT2189
                                                                        ABJT2190
c ...    calculate the heat capacity                                    ABJT2191
         heatcp = ( energ2-enrgy ) / ( dtheat*tp )                      ABJT2192
                                                                        ABJT2193
c ...    calculate the temperature derivative of the average            ABJT2194
c        charge state                                                   ABJT2195
         dzdt = ( eld2-densne ) / ( dtheat*tp*densnn )                  ABJT2196
                                                                        ABJT2197
c ...    calculate the density derivative of the specific energy        ABJT2198
         dedden = ( energ3-enrgy ) / ( dtheat*densnn )                  ABJT2199
                                                                        ABJT2200
      else                                                              ABJT2201
                                                                        ABJT2202
         heatcp = 0.                                                    ABJT2203
         dzdt = 0.                                                      ABJT2204
         dedden = 0.                                                    ABJT2205
                                                                        ABJT2206
      endif                                                             ABJT2207
                                                                        ABJT2208
      return                                                            ABJT2209
                                                                        ABJT2210
c ... format statements                                                 ABJT2211
                                                                        ABJT2212
  901 format (' debug output from -eos-:'/)                             ABJT2213
                                                                        ABJT2214
      end                                                               ABJT2215
c NON-LTE PART
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ABJT2216
c MOVED INTO ModMultiGroup                                              ABJT2217
      function gint ( igint,xmin,xmax )                                 ABJT2218
                                                                        ABJT2219
c ... this function returns the integral <from xmin to xmax> of         ABJT2220
c     GINTn (where n = 1 thru 6)                                        ABJT2221
                                                                        ABJT2222
      logical lbug                                                      ABJT2223
                                                                        ABJT2224
c ... check debug option                                                ABJT2225
      call debug ( 'gint  ',                                            ABJT2226
     &                        lbug )                                    ABJT2227
                                                                        ABJT2228
c ... compute the integral                                              ABJT2229
                                                                        ABJT2230
      go to ( 10,20,30,40,50,60 ), igint                                ABJT2231
                                                                        ABJT2232
   10 gint = gint1(xmax) - gint1(xmin)                                  ABJT2233
      go to 100                                                         ABJT2234
   20 gint = gint2(xmax) - gint2(xmin)                                  ABJT2235
      go to 100                                                         ABJT2236
   30 gint = gint3(xmax) - gint3(xmin)                                  ABJT2237
      go to 100                                                         ABJT2238
   40 gint = gint4(xmax) - gint4(xmin)                                  ABJT2239
      go to 100                                                         ABJT2240
   50 gint = gint5(xmax) - gint5(xmin)                                  ABJT2241
      go to 100                                                         ABJT2242
   60 gint = gint6(xmax) - gint6(xmin)                                  ABJT2243
                                                                        ABJT2244
  100 continue                                                          ABJT2245
      if ( lbug ) write (6,900) igint,xmin,xmax,gint                    ABJT2246
                                                                        ABJT2247
                                                                        ABJT2248
      return                                                            ABJT2249
                                                                        ABJT2250
c ... format statements                                                 ABJT2251
                                                                        ABJT2252
  900 format (' debug output from -gint-:'/t4,'igint',t18,              ABJT2253
     &  'xmin',t32,'xmax',t46,'gint'/t4,i4,t18,1p3e14.4)                ABJT2254
                                                                        ABJT2255
      end                                                               ABJT2256
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ABJT2257
                                                                        ABJT2258
      function gint1 ( x )                                              ABJT2259
                                                                        ABJT2260
c ... this routine returns the integral (from 0 to x) of the            ABJT2261
c     function:                                                         ABJT2262
c                                                                       ABJT2263
c                     x_                                                ABJT2264
c           gint1  = _/  dy exp(-y)                                     ABJT2265
c                   0                                                   ABJT2266
c                                                                       ABJT2267
      logical lbug                                                      ABJT2268
                                                                        ABJT2269
      gint1 = -exp( -x )                                                ABJT2270
                                                                        ABJT2271
      return                                                            ABJT2272
      end                                                               ABJT2273
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ABJT2274
                                                                        ABJT2275
      function gint2 ( x )                                              ABJT2276
                                                                        ABJT2277
c ... this routine returns the integral (from 0 to x) of the            ABJT2278
c     function:                                                         ABJT2279
c                                                                       ABJT2280
c                     x_                                                ABJT2281
c           gint2  = _/  dy y**3 exp(-y)                                ABJT2282
c                   0                                                   ABJT2283
c                                                                       ABJT2284
      gint2 = -exp( -x ) * ( x**3 + 3.*x**2 + 6.*x + 6. )               ABJT2285
                                                                        ABJT2286
      return                                                            ABJT2287
      end                                                               ABJT2288
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ABJT2289
                                                                        ABJT2290
      function gint3 ( x )                                              ABJT2291
                                                                        ABJT2292
c ... this routine returns the integral (from 0 to x) of the            ABJT2293
c     function:                                                         ABJT2294
c                                                                       ABJT2295
c                     x_                                                ABJT2296
c           gint3  = _/  dy y**4 exp(-y) / ( 1-exp(-y) )**3             ABJT2297
c                   0                                                   ABJT2298
c                                                                       ABJT2299
c ... the table values are the logarithms of the integrals for          ABJT2300
c     equally spaced values of  log(x).                                 ABJT2301
c                                                                       ABJT2302
      dimension gtable(50)                                              ABJT2303
      logical lbug                                                      ABJT2304
                                                                        ABJT2305
      data xmin / .005 /, xmax / 30. /, dxlog / .17754 /, npts / 50 /   ABJT2306
      data xmnlog / -5.2983 /                                           ABJT2307
      data gtable /                                                     ABJT2308
     1 -1.1288E+01,-1.0933E+01,-1.0577E+01,-1.0222E+01,-9.8661E+00,     ABJT2309
     2 -9.5103E+00,-9.1545E+00,-8.7984E+00,-8.4422E+00,-8.0858E+00,     ABJT2310
     3 -7.7292E+00,-7.3722E+00,-7.0149E+00,-6.6571E+00,-6.2988E+00,     ABJT2311
     4 -5.9399E+00,-5.5803E+00,-5.2199E+00,-4.8584E+00,-4.4958E+00,     ABJT2312
     5 -4.1318E+00,-3.7661E+00,-3.3986E+00,-3.0290E+00,-2.6568E+00,     ABJT2313
     6 -2.2819E+00,-1.9038E+00,-1.5224E+00,-1.1374E+00,-7.4873E-01,     ABJT2314
     7 -3.5677E-01, 3.7811E-02, 4.3361E-01, 8.2825E-01, 1.2180E+00,     ABJT2315
     8  1.5975E+00, 1.9592E+00, 2.2940E+00, 2.5913E+00, 2.8406E+00,     ABJT2316
     9  3.0339E+00, 3.1686E+00, 3.2499E+00, 3.2903E+00, 3.3058E+00,     ABJT2317
     &  3.3101E+00, 3.3109E+00, 3.3110E+00, 3.3110E+00, 3.3110E+00/     ABJT2318
                                                                        ABJT2319
c ********************************************************************  ABJT2320
c                                                                       ABJT2321
c                         begin execution                               ABJT2322
c                                                                       ABJT2323
c ********************************************************************  ABJT2324
                                                                        ABJT2325
c ... check debug option                                                ABJT2326
      call debug ( 'gint3',                                             ABJT2327
     &                       lbug )                                     ABJT2328
                                                                        ABJT2329
      if ( lbug ) write (6,900) x                                       ABJT2330
                                                                        ABJT2331
      if ( x .lt. xmin ) then                                           ABJT2332
                                                                        ABJT2333
c ...    analytic form using expansion of exp(x)                        ABJT2334
         gint3 = x**2 / 2.                                              ABJT2335
                                                                        ABJT2336
      else if ( x .gt. xmax*0.9999 ) then                               ABJT2337
                                                                        ABJT2338
c ...    the integral as x -> infinity                                  ABJT2339
         gint3 = exp( 3.3110 )                                          ABJT2340
                                                                        ABJT2341
      else                                                              ABJT2342
                                                                        ABJT2343
c ...    otherwise, interpolate using table                             ABJT2344
         xlog = log( x )                                                ABJT2345
         ii = 1 + ( xlog-xmnlog ) / dxlog                               ABJT2346
         xii = xmnlog + ( ii-1 )*dxlog                                  ABJT2347
         ex = gtable(ii) + ( xlog-xii ) *                               ABJT2348
     2              ( gtable(ii+1) - gtable(ii) ) / dxlog               ABJT2349
         gint3 = exp( ex )                                              ABJT2350
                                                                        ABJT2351
      endif                                                             ABJT2352
                                                                        ABJT2353
      if ( lbug ) write (6,901) x,gint3,ii,xlog,xii,ex                  ABJT2354
                                                                        ABJT2355
      return                                                            ABJT2356
                                                                        ABJT2357
c ... format statements                                                 ABJT2358
                                                                        ABJT2359
  900 format (' debug output from -gint3-:    x=',1pe14.4)              ABJT2360
  901 format (t4,'x',t18,'gint3',t32,'ii',t46,'xlog',t60,               ABJT2361
     2  'xii',t74,'ex'/t4,1p2e14.4,0p,i4,t46,1p3e14.4)                  ABJT2362
                                                                        ABJT2363
      end                                                               ABJT2364
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ABJT2365
                                                                        ABJT2366
      function gint4 ( x )                                              ABJT2367
                                                                        ABJT2368
c ... this routine returns the integral (from 0 to x) of the            ABJT2369
c     function:                                                         ABJT2370
c                                                                       ABJT2371
c                     x_                                                ABJT2372
c           gint4  = _/  dy y**7 exp(-y) / ( 1-exp(-y) )**3             ABJT2373
c                   0                                                   ABJT2374
c                                                                       ABJT2375
c ... the table values are the logarithms of the integrals for          ABJT2376
c     equally spaced values of  log(x).                                 ABJT2377
c                                                                       ABJT2378
      dimension gtable(50)                                              ABJT2379
      logical lbug                                                      ABJT2380
                                                                        ABJT2381
      data xmin / .005 /, xmax / 30. /, dxlog / .17754 /, npts / 50 /   ABJT2382
      data xmnlog / -5.2983 /                                           ABJT2383
      data gtable /                                                     ABJT2384
     1 -2.8099E+01,-2.7211E+01,-2.6323E+01,-2.5434E+01,-2.4546E+01,     ABJT2385
     2 -2.3657E+01,-2.2769E+01,-2.1880E+01,-2.0991E+01,-2.0101E+01,     ABJT2386
     3 -1.9212E+01,-1.8322E+01,-1.7431E+01,-1.6540E+01,-1.5648E+01,     ABJT2387
     4 -1.4756E+01,-1.3863E+01,-1.2968E+01,-1.2073E+01,-1.1176E+01,     ABJT2388
     5 -1.0277E+01,-9.3763E+00,-8.4734E+00,-7.5678E+00,-6.6594E+00,     ABJT2389
     6 -5.7477E+00,-4.8324E+00,-3.9134E+00,-2.9904E+00,-2.0638E+00,     ABJT2390
     7 -1.1342E+00,-2.0272E-01, 7.2834E-01, 1.6554E+00, 2.5731E+00,     ABJT2391
     8  3.4736E+00, 4.3462E+00, 5.1772E+00, 5.9497E+00, 6.6448E+00,     ABJT2392
     9  7.2434E+00, 7.7287E+00, 8.0902E+00, 8.3287E+00, 8.4606E+00,     ABJT2393
     &  8.5172E+00, 8.5343E+00, 8.5375E+00, 8.5379E+00, 8.5379E+00/     ABJT2394
                                                                        ABJT2395
c ********************************************************************  ABJT2396
c                                                                       ABJT2397
c                         begin execution                               ABJT2398
c                                                                       ABJT2399
c ********************************************************************  ABJT2400
                                                                        ABJT2401
c ... check debug option                                                ABJT2402
      call debug ( 'gint4',                                             ABJT2403
     &                       lbug )                                     ABJT2404
                                                                        ABJT2405
      if ( lbug ) write (6,900) x                                       ABJT2406
                                                                        ABJT2407
      if ( x .lt. xmin ) then                                           ABJT2408
                                                                        ABJT2409
c ...    analytic form using expansion of exp(x)                        ABJT2410
         gint4 = x**5 / 5.                                              ABJT2411
                                                                        ABJT2412
      else if ( x .gt. xmax*0.9999 ) then                               ABJT2413
                                                                        ABJT2414
c ...    the integral as x -> infinity                                  ABJT2415
         gint4 = exp( 8.5379 )                                          ABJT2416
                                                                        ABJT2417
      else                                                              ABJT2418
                                                                        ABJT2419
c ...    otherwise, interpolate using table                             ABJT2420
         xlog = log( x )                                                ABJT2421
         ii = 1 + ( xlog-xmnlog ) / dxlog                               ABJT2422
         xii = xmnlog + ( ii-1 )*dxlog                                  ABJT2423
         ex = gtable(ii) + ( xlog-xii ) *                               ABJT2424
     2              ( gtable(ii+1) - gtable(ii) ) / dxlog               ABJT2425
         gint4 = exp( ex )                                              ABJT2426
                                                                        ABJT2427
      endif                                                             ABJT2428
                                                                        ABJT2429
      if ( lbug ) write (6,901) x,gint4,ii,xlog,xii,ex                  ABJT2430
                                                                        ABJT2431
      return                                                            ABJT2432
                                                                        ABJT2433
c ... format statements                                                 ABJT2434
                                                                        ABJT2435
  900 format (' debug output from -gint4-:    x=',1pe14.4)              ABJT2436
  901 format (t4,'x',t18,'gint4',t32,'ii',t46,'xlog',t60,               ABJT2437
     2  'xii',t74,'ex'/t4,1p2e14.4,0p,i4,t46,1p3e14.4)                  ABJT2438
                                                                        ABJT2439
      end                                                               ABJT2440
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ABJT2441
                                                                        ABJT2442
      function gint5 ( x )                                              ABJT2443
                                                                        ABJT2444
c ... this routine returns the integral (from 0 to x) of the            ABJT2445
c     function:                                                         ABJT2446
c                                                                       ABJT2447
c                     x_                                                ABJT2448
c           gint5  = _/  dy y**3 / ( exp(y)-1 )                         ABJT2449
c                   0                                                   ABJT2450
c                                                                       ABJT2451
c ... the table values are the logarithms of the integrals for          ABJT2452
c     equally spaced values of  log(x).                                 ABJT2453
c                                                                       ABJT2454
      dimension gtable(50)                                              ABJT2455
      logical lbug                                                      ABJT2456
                                                                        ABJT2457
      data xmin / .005 /, xmax / 30. /, dxlog / .17754 /, npts / 50 /   ABJT2458
      data xmnlog / -5.2983 /                                           ABJT2459
      data gtable /                                                     ABJT2460
     1 -1.6995E+01,-1.6463E+01,-1.5931E+01,-1.5399E+01,-1.4867E+01,     ABJT2461
     2 -1.4335E+01,-1.3803E+01,-1.3272E+01,-1.2740E+01,-1.2209E+01,     ABJT2462
     3 -1.1678E+01,-1.1148E+01,-1.0618E+01,-1.0088E+01,-9.5594E+00,     ABJT2463
     4 -9.0312E+00,-8.5039E+00,-7.9775E+00,-7.4524E+00,-6.9289E+00,     ABJT2464
     5 -6.4070E+00,-5.8874E+00,-5.3703E+00,-4.8563E+00,-4.3460E+00,     ABJT2465
     6 -3.8402E+00,-3.3399E+00,-2.8462E+00,-2.3605E+00,-1.8846E+00,     ABJT2466
     7 -1.4204E+00,-9.7072E-01,-5.3858E-01,-1.2780E-01, 2.5710E-01,     ABJT2467
     8  6.1092E-01, 9.2788E-01, 1.2021E+00, 1.4283E+00, 1.6031E+00,     ABJT2468
     9  1.7268E+00, 1.8043E+00, 1.8456E+00, 1.8634E+00, 1.8693E+00,     ABJT2469
     &  1.8706E+00, 1.8708E+00, 1.8709E+00, 1.8709E+00, 1.8709E+00 /    ABJT2470
                                                                        ABJT2471
c ********************************************************************  ABJT2472
c                                                                       ABJT2473
c                         begin execution                               ABJT2474
c                                                                       ABJT2475
c ********************************************************************  ABJT2476
                                                                        ABJT2477
c ... check debug option                                                ABJT2478
      call debug ( 'gint5',                                             ABJT2479
     &                       lbug )                                     ABJT2480
                                                                        ABJT2481
      if ( lbug ) write (6,900) x                                       ABJT2482
                                                                        ABJT2483
      if ( x .lt. xmin ) then                                           ABJT2484
                                                                        ABJT2485
c ...    analytic form using expansion of exp(x)                        ABJT2486
         gint5 = x**3 / 3.                                              ABJT2487
                                                                        ABJT2488
      else if ( x .gt. xmax*0.9999 ) then                               ABJT2489
                                                                        ABJT2490
c ...    the integral as x -> infinity                                  ABJT2491
         gint5 = exp( 1.8709 )                                          ABJT2492
                                                                        ABJT2493
      else                                                              ABJT2494
                                                                        ABJT2495
c ...    otherwise, interpolate using table                             ABJT2496
         xlog = log( x )                                                ABJT2497
         ii = 1 + ( xlog-xmnlog ) / dxlog                               ABJT2498
         xii = xmnlog + ( ii-1 )*dxlog                                  ABJT2499
         ex = gtable(ii) + ( xlog-xii ) *                               ABJT2500
     2              ( gtable(ii+1) - gtable(ii) ) / dxlog               ABJT2501
         gint5 = exp( ex )                                              ABJT2502
                                                                        ABJT2503
      endif                                                             ABJT2504
                                                                        ABJT2505
      if ( lbug ) write (6,901) x,gint5,ii,xlog,xii,ex                  ABJT2506
                                                                        ABJT2507
      return                                                            ABJT2508
                                                                        ABJT2509
c ... format statements                                                 ABJT2510
                                                                        ABJT2511
  900 format (' debug output from -gint5-:    x=',1pe14.4)              ABJT2512
  901 format (t4,'x',t18,'gint5',t32,'ii',t46,'xlog',t60,               ABJT2513
     2  'xii',t74,'ex'/t4,1p2e14.4,0p,i4,t46,1p3e14.4)                  ABJT2514
                                                                        ABJT2515
      end                                                               ABJT2516
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ABJT2517
                                                                        ABJT2518
      function gint6 ( x )                                              ABJT2519
                                                                        ABJT2520
c ... this routine returns the integral (from 0 to x) of the            ABJT2521
c     function:                                                         ABJT2522
c                                                                       ABJT2523
c                     x_                                                ABJT2524
c           gint6  = _/  dy y**4 exp(-y) / ( 1-exp(-y) )**2             ABJT2525
c                   0                                                   ABJT2526
c                                                                       ABJT2527
c ... the table values are the logarithms of the integrals for          ABJT2528
c     equally spaced values of  log(x).                                 ABJT2529
c                                                                       ABJT2530
      dimension gtable(50)                                              ABJT2531
      logical lbug                                                      ABJT2532
                                                                        ABJT2533
      data xmin / .005 /, xmax / 30. /, dxlog / .17754 /, npts / 50 /   ABJT2534
      data xmnlog / -5.2983 /                                           ABJT2535
      data gtable /                                                     ABJT2536
     1 -1.6994E+01,-1.6461E+01,-1.5928E+01,-1.5396E+01,-1.4863E+01,     ABJT2537
     2 -1.4330E+01,-1.3798E+01,-1.3265E+01,-1.2733E+01,-1.2200E+01,     ABJT2538
     3 -1.1667E+01,-1.1135E+01,-1.0602E+01,-1.0070E+01,-9.5370E+00,     ABJT2539
     4 -9.0045E+00,-8.4720E+00,-7.9395E+00,-7.4071E+00,-6.8748E+00,     ABJT2540
     5 -6.3426E+00,-5.8106E+00,-5.2789E+00,-4.7476E+00,-4.2169E+00,     ABJT2541
     6 -3.6869E+00,-3.1581E+00,-2.6309E+00,-2.1060E+00,-1.5843E+00,     ABJT2542
     7 -1.0671E+00,-5.5646E-01,-5.4768E-02, 4.3443E-01, 9.0650E-01,     ABJT2543
     8  1.3554E+00, 1.7736E+00, 2.1524E+00, 2.4822E+00, 2.7542E+00,     ABJT2544
     9  2.9625E+00, 3.1063E+00, 3.1926E+00, 3.2353E+00, 3.2516E+00,     ABJT2545
     &  3.2562E+00, 3.2571E+00, 3.2572E+00, 3.2572E+00, 3.2572E+00 /    ABJT2546
                                                                        ABJT2547
c ********************************************************************  ABJT2548
c                                                                       ABJT2549
c                         begin execution                               ABJT2550
c                                                                       ABJT2551
c ********************************************************************  ABJT2552
                                                                        ABJT2553
c ... check debug option                                                ABJT2554
      call debug ( 'gint6',                                             ABJT2555
     &                       lbug )                                     ABJT2556
                                                                        ABJT2557
      if ( lbug ) write (6,900) x                                       ABJT2558
                                                                        ABJT2559
      if ( x .lt. xmin ) then                                           ABJT2560
                                                                        ABJT2561
c ...    analytic form using expansion of exp(x)                        ABJT2562
         gint6 = x**3 / 3.                                              ABJT2563
                                                                        ABJT2564
      else if ( x .gt. xmax*0.9999 ) then                               ABJT2565
                                                                        ABJT2566
c ...    the integral as x -> infinity                                  ABJT2567
         gint6 = exp( 3.2572 )                                          ABJT2568
                                                                        ABJT2569
      else                                                              ABJT2570
                                                                        ABJT2571
c ...    otherwise, interpolate using table                             ABJT2572
         xlog = log( x )                                                ABJT2573
         ii = 1 + ( xlog-xmnlog ) / dxlog                               ABJT2574
         xii = xmnlog + ( ii-1 )*dxlog                                  ABJT2575
         ex = gtable(ii) + ( xlog-xii ) *                               ABJT2576
     2              ( gtable(ii+1) - gtable(ii) ) / dxlog               ABJT2577
         gint6 = exp( ex )                                              ABJT2578
                                                                        ABJT2579
      endif                                                             ABJT2580
                                                                        ABJT2581
      if ( lbug ) write (6,901) x,gint6,ii,xlog,xii,ex                  ABJT2582
                                                                        ABJT2583
      return                                                            ABJT2584
                                                                        ABJT2585
c ... format statements                                                 ABJT2586
                                                                        ABJT2587
  900 format (' debug output from -gint6-:    x=',1pe14.4)              ABJT2588
  901 format (t4,'x',t18,'gint6',t32,'ii',t46,'xlog',t60,               ABJT2589
     2  'xii',t74,'ex'/t4,1p2e14.4,0p,i4,t46,1p3e14.4)                  ABJT2590
                                                                        ABJT2591
      end                                                               ABJT2592
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT2593
                                                                        ABJT2594
      subroutine input                                                  ABJT2595
c                                                                       ABJT2596
c ... set default values for input data, read in namelist input and     ABJT2597
c     atomic ionization potential files, and write out various          ABJT2598
c     information.                                                      ABJT2599
c                                                                       ABJT2600
c ... data that are read in are loaded into common blocks               ABJT2601
c                                                                       ABJT2602
                                                                        ABJT2603
c ... set the maximum number of:                                        ABJT2604
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT2605
c       temperatures (mxtemp), densities (mxdens),                      ABJT2606
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT2607
                                                                        ABJT2608
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT2609
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT2610
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT2611
                                                                        ABJT2612
                                                                        ABJT2613
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT2614
c ................................................................      ABJT2615
                                                                        ABJT2616
      common / dbugcm / nobug,ncycld(10),icycld(10)                     ABJT2617
                                                                        ABJT2618
      logical nobug                                                     ABJT2619
c ................................................................      ABJT2620
                                                                        ABJT2621
      common / gases / ngases,izgas(mxgass),atomwt(mxgass),             ABJT2622
     2                 avgatw,avgatn,                                   ABJT2623
     3                 fracsp(mxgass),fraciz(mxgass,mxatom),            ABJT2624
     4                 fraclv(mxgass,mxatom,20),pot(mxgass,mxatom),     ABJT2625
     5                 numocc(mxgass,mxatom)                            ABJT2626
c ................................................................      ABJT2627
                                                                        ABJT2628
      common / opacs / nptspg, ngrups, engrup(mxgrps+1)                 ABJT2629
c ................................................................      ABJT2630
                                                                        ABJT2631
      common / params / npqmax(mxgass), npmaxp, nfrqbb,                 ABJT2632
     2                  npring(100), neopen(5,100), bbnton(6,55),       ABJT2633
     3                  defpot(27,55), noccdf(27,55)                    ABJT2634
c ................................................................      ABJT2635
                                                                        ABJT2636
      common / strngs / header, namedb(10)                              ABJT2637
                                                                        ABJT2638
      character header*80, namedb*6                                     ABJT2639
c ................................................................      ABJT2640
                                                                        ABJT2641
      common / tdmesh / ntemp,ndens,ntrad,                              ABJT2642
     2                  tplsma(mxtemp),densnn(mxdens),trad(mxtemp),     ABJT2643
     2                  dlgtmp,dlgden,dlgtrd                            ABJT2644
c ................................................................      ABJT2645
                                                                        ABJT2646
c ... declaration statements                                            ABJT2647
      dimension grupbd(mxgrps+1)                                        ABJT2648
                                                                        ABJT2649
      character potfil*10, pltfil*12                                    ABJT2650
      character ctime*8, cdate*8, cmonth(12)*4                          ABJT2651
      character cstrng*80                                               ABJT2652
                                                                        ABJT2653
      data cmonth / 'Jan.','Feb.','Mar.','Apr.',' May','June',          ABJT2654
     2              'July','Aug.','Sep.','Oct.','Nov.','Dec.' /         ABJT2655
                                                                        ABJT2656
      namelist / data / ngases,  izgas,   atomwt,  fracsp,              ABJT2657
     2         ntemp,   ndens,   tplsma,  densnn,  dlgtmp,  dlgden,     ABJT2658
     3         ntrad,   trad,    dlgtrd,                                ABJT2659
     4         nptspg,  ngrups,  nfrqbb,  grupbd,                       ABJT2660
     5         isw,     con,     iplot,   dtheat,  critsc               ABJT2661
                                                                        ABJT2662
c ********************************************************************  ABJT2663
c                                                                       ABJT2664
c                              begin execution                          ABJT2665
c                                                                       ABJT2666
c ********************************************************************  ABJT2667
                                                                        ABJT2668
c ... set default values for namelist input                             ABJT2669
                                                                        ABJT2670
      ngases = 1                                                        ABJT2671
      ntemp = 1                                                         ABJT2672
      ndens = 1                                                         ABJT2673
      ntrad = 0                                                         ABJT2674
      dlgtmp = 0.                                                       ABJT2675
      dlgden = 0.                                                       ABJT2676
      dlgtrd = 0.                                                       ABJT2677
                                                                        ABJT2678
      dtheat = 0.01                                                     ABJT2679
      critsc = 10.                                                      ABJT2680
                                                                        ABJT2681
      nptspg = 5                                                        ABJT2682
      ngrups = 20                                                       ABJT2683
      nfrqbb = 5                                                        ABJT2684
                                                                        ABJT2685
c ... default energy group boundaries (in eV)                           ABJT2686
                                                                        ABJT2687
      engrup( 1) = 1.00e-1                                              ABJT2688
      engrup( 2) = 3.16e-1                                              ABJT2689
      engrup( 3) = 1.00e0                                               ABJT2690
      engrup( 4) = 1.58e0                                               ABJT2691
      engrup( 5) = 2.51e0                                               ABJT2692
      engrup( 6) = 3.97e0                                               ABJT2693
      engrup( 7) = 6.29e0                                               ABJT2694
      engrup( 8) = 1.00e1                                               ABJT2695
      engrup( 9) = 1.58e1                                               ABJT2696
      engrup(10) = 2.51e1                                               ABJT2697
      engrup(11) = 3.97e1                                               ABJT2698
      engrup(12) = 6.29e1                                               ABJT2699
      engrup(13) = 1.00e2                                               ABJT2700
      engrup(14) = 1.78e2                                               ABJT2701
      engrup(15) = 3.16e2                                               ABJT2702
      engrup(16) = 5.62e2                                               ABJT2703
      engrup(17) = 1.00e3                                               ABJT2704
      engrup(18) = 3.16e3                                               ABJT2705
      engrup(19) = 1.00e4                                               ABJT2706
      engrup(20) = 1.00e5                                               ABJT2707
      engrup(21) = 1.00e6                                               ABJT2708
                                                                        ABJT2709
      do 20 i=1,mxgrps+1                                                ABJT2710
         grupbd(i) = 0.                                                 ABJT2711
   20 continue                                                          ABJT2712
      do 30 i=1,mxtemp                                                  ABJT2713
         tplsma(i) = 0.                                                 ABJT2714
         trad(i)   = 0.                                                 ABJT2715
   30 continue                                                          ABJT2716
      do 35 i=1,mxdens                                                  ABJT2717
         densnn(i) = 0.                                                 ABJT2718
   35 continue                                                          ABJT2719
                                                                        ABJT2720
      do 40 i=1,30                                                      ABJT2721
         isw(i) = 0                                                     ABJT2722
         con(i) = 0.                                                    ABJT2723
         iplot(i) = 0                                                   ABJT2724
   40 continue                                                          ABJT2725
                                                                        ABJT2726
      isw(6)  = 3                                                       ABJT2727
      isw(15) = 2                                                       ABJT2728
      isw(19) = 2                                                       ABJT2729
                                                                        ABJT2730
      con(2) = 1.e-10                                                   ABJT2731
      con(3) = 1.e-10                                                   ABJT2732
      con(4) = 1.e-10                                                   ABJT2733
      con(5) = 1.e10                                                    ABJT2734
      con(6) = 10.                                                      ABJT2735
      con(7) = 10.                                                      ABJT2736
      con(9) = 1.                                                       ABJT2737
                                                                        ABJT2738
      do 100 lgas=1,mxgass                                              ABJT2739
                                                                        ABJT2740
         izgas(lgas) = 0                                                ABJT2741
         atomwt(lgas) = 0.                                              ABJT2742
         fracsp(lgas) = 0.                                              ABJT2743
                                                                        ABJT2744
         do 80 iz=1,mxatom                                              ABJT2745
                                                                        ABJT2746
            pot(lgas,iz) = 0.                                           ABJT2747
            numocc(lgas,iz) = 0                                         ABJT2748
            fraciz(lgas,iz) = 0.                                        ABJT2749
            do 60 nprin=1,20                                            ABJT2750
               fraclv(lgas,iz,nprin) = 0.                               ABJT2751
   60       continue                                                    ABJT2752
                                                                        ABJT2753
   80    continue                                                       ABJT2754
                                                                        ABJT2755
  100 continue                                                          ABJT2756
                                                                        ABJT2757
c ... default values for debug variables                                ABJT2758
      nobug = .true.                                                    ABJT2759
      do 110 i=1,10                                                     ABJT2760
         namedb(i) = '      '                                           ABJT2761
         ncycld(i) = 1                                                  ABJT2762
         icycld(i) = 0                                                  ABJT2763
  110 continue                                                          ABJT2764
                                                                        ABJT2765
                                                                        ABJT2766
c ... read namelist input                                               ABJT2767
c     -------------------                                               ABJT2768
                                                                        ABJT2769
      read ( 5,data )                                                   ABJT2770
                                                                        ABJT2771
c ... read debug info                                                   ABJT2772
c     ---------------                                                   ABJT2773
                                                                        ABJT2774
      if ( isw(3) .ne. 0 ) then                                         ABJT2775
                                                                        ABJT2776
         read (5,811) (namedb(i),i=1,isw(3) )                           ABJT2777
         read (5,812) (ncycld(i),i=1,isw(3) )                           ABJT2778
                                                                        ABJT2779
      endif                                                             ABJT2780
                                                                        ABJT2781
c ... set up a header record listing the time and date of the run       ABJT2782
                                                                        ABJT2783
      call time ( ctime )                                               ABJT2784
      call idate ( imon,iday,iyear )                                    ABJT2785
      write (header,955) cmonth(imon),iday,iyear,ctime                  ABJT2786
      write (6,911) header                                              ABJT2787
                                                                        ABJT2788
c ... if this calculation requests data to be written for use by CONRAD.ABJT2789
c     check the input variables to make sure they are OK                ABJT2790
                                                                        ABJT2791
      if ( isw(8).eq.1 .or. isw(8).eq.12 .or. isw(8).eq.13 ) then       ABJT2792
                                                                        ABJT2793
c ...    write to file with old format used by CONRAD                   ABJT2794
                                                                        ABJT2795
         if ( (ntemp.ne.20)  .or. (ntrad.ne.20) .or.                    ABJT2796
     &        (ndens.ne.10)  .or. (dlgtmp.eq.0.) .or.                   ABJT2797
     &        (dlgtrd.eq.0.) .or. (dlgden.eq.0.) ) then                 ABJT2798
            write (6,941) ntemp,ndens,ntrad,dlgtmp,dlgden,dlgtrd        ABJT2799
            write (6,942) 20,10,20                                      ABJT2800
            stop                                                        ABJT2801
     &     ' improper input for output to be written for CONRAD'        ABJT2802
         endif                                                          ABJT2803
                                                                        ABJT2804
c ...    open a file for the CONRAD data                                ABJT2805
         open (unit=8,file='cnrdeos0',status='new',form='formatted')    ABJT2806
                                                                        ABJT2807
      endif                                                             ABJT2808
                                                                        ABJT2809
      if ( isw(8).eq.2 .or. isw(8).eq.12 ) then                         ABJT2810
                                                                        ABJT2811
         stop ' SESAME-like files can only be written with CRAY runs'   ABJT2812
                                                                        ABJT2813
      endif                                                             ABJT2814
                                                                        ABJT2815
      if ( isw(8).eq.3 .or. isw(8).eq.13 ) then                         ABJT2816
                                                                        ABJT2817
c ...    write to file with new (post-9/87) format used by CONRAD       ABJT2818
                                                                        ABJT2819
         if ( (ntemp.ne.20)  .or. (ndens.ne.10)  .or.                   ABJT2820
     &        (dlgtmp.eq.0.) .or. (dlgden.eq.0.) ) then                 ABJT2821
            write (6,943) ntemp,ndens,dlgtmp,dlgden                     ABJT2822
            write (6,944) 20,10                                         ABJT2823
            stop                                                        ABJT2824
     &     ' improper input for output to be written for CONRAD'        ABJT2825
         endif                                                          ABJT2826
                                                                        ABJT2827
c ...    open a file for the CONRAD data                                ABJT2828
         open (unit=10,file='cnrdeos',status='new',form='formatted')    ABJT2829
                                                                        ABJT2830
      endif                                                             ABJT2831
                                                                        ABJT2832
c ... set up the temperature and density grid                           ABJT2833
                                                                        ABJT2834
      if ( dlgtmp.ne.0. ) then                                          ABJT2835
         tlog = log10( tplsma(1) )                                      ABJT2836
         do 310 itemp=2,ntemp                                           ABJT2837
            tlog = tlog + dlgtmp                                        ABJT2838
            tplsma(itemp) = 10.**tlog                                   ABJT2839
  310    continue                                                       ABJT2840
      endif                                                             ABJT2841
                                                                        ABJT2842
      if ( dlgden.ne.0. ) then                                          ABJT2843
         dlog = log10( densnn(1) )                                      ABJT2844
         do 320 idens=2,ndens                                           ABJT2845
            dlog = dlog + dlgden                                        ABJT2846
            densnn(idens) = 10.**dlog                                   ABJT2847
  320    continue                                                       ABJT2848
      endif                                                             ABJT2849
                                                                        ABJT2850
      if ( ntrad.ne.0 ) then                                            ABJT2851
         if ( dlgtrd.ne.0. ) then                                       ABJT2852
            tlog = log10( trad(1) )                                     ABJT2853
            do 325 itrad=2,ntrad                                        ABJT2854
              tlog = tlog + dlgtrd                                      ABJT2855
              trad(itrad) = 10.**tlog                                   ABJT2856
  325       continue                                                    ABJT2857
         endif                                                          ABJT2858
      endif                                                             ABJT2859
                                                                        ABJT2860
                                                                        ABJT2861
c ... normalize the number fractions to 1; then calculate the           ABJT2862
c     average atomic weight of the gas                                  ABJT2863
                                                                        ABJT2864
      summ = 0.                                                         ABJT2865
      do 328 lgas=1,ngases                                              ABJT2866
         summ = summ + fracsp(lgas)                                     ABJT2867
  328 continue                                                          ABJT2868
                                                                        ABJT2869
      avgatw = 0.                                                       ABJT2870
      avgatn = 0.                                                       ABJT2871
      do 330 lgas=1,ngases                                              ABJT2872
         if ( izgas(lgas) .ge. mxatom ) then                            ABJT2873
            stop ' array size too small; increase -mxatom-'             ABJT2874
         endif                                                          ABJT2875
         fracsp(lgas) = fracsp(lgas) / summ                             ABJT2876
         avgatw = avgatw + fracsp(lgas) * atomwt(lgas)                  ABJT2877
         avgatn = avgatn + fracsp(lgas) * izgas(lgas)                   ABJT2878
  330 continue                                                          ABJT2879
                                                                        ABJT2880
                                                                        ABJT2881
c ... write out some info                                               ABJT2882
                                                                        ABJT2883
      write (6,906) ngases,(i,izgas(i),fracsp(i),atomwt(i),i=1,ngases)  ABJT2884
                                                                        ABJT2885
                                                                        ABJT2886
c ... set debug flag                                                    ABJT2887
      if ( namedb(1) .ne. '      ' ) nobug = .false.                    ABJT2888
                                                                        ABJT2889
c ... photon energy group boundaries                                    ABJT2890
                                                                        ABJT2891
      if ( isw(13) .eq. 1 ) then                                        ABJT2892
         do 350 i=1,ngrups+1                                            ABJT2893
            engrup(i) = grupbd(i)                                       ABJT2894
  350    continue                                                       ABJT2895
      endif                                                             ABJT2896
                                                                        ABJT2897
      if ( isw(1) .ne. 0 ) then                                         ABJT2898
                                                                        ABJT2899
c ...    read the ionization potentials and the number of               ABJT2900
c        bound electrons in the outermost shell                         ABJT2901
                                                                        ABJT2902
         do 250 lgas=1,ngases                                           ABJT2903
                                                                        ABJT2904
c ...       set the name of the input file                              ABJT2905
            write ( potfil,951 ) izgas(lgas)                            ABJT2906
                                                                        ABJT2907
c ...       open the file                                               ABJT2908
            open ( unit=7,file=potfil,status='old',form='formatted' )   ABJT2909
                                                                        ABJT2910
            write (6,901) lgas,izgas(lgas)                              ABJT2911
                                                                        ABJT2912
c ...       read in the data                                            ABJT2913
            do 240 iz=1,izgas(lgas)                                     ABJT2914
               read (7,*) pot(lgas,iz)                                  ABJT2915
c ...          since "numocc" is only used for the old Mosher formula,  ABJT2916
c              just set this equal to 1                                 ABJT2917
               numocc(lgas,iz) = 1                                      ABJT2918
               write (6,902) iz-1,pot(lgas,iz)                          ABJT2919
  240       continue                                                    ABJT2920
                                                                        ABJT2921
c ...       close the file                                              ABJT2922
            close ( unit=7 )                                            ABJT2923
                                                                        ABJT2924
  250    continue                                                       ABJT2925
                                                                        ABJT2926
      else                                                              ABJT2927
                                                                        ABJT2928
c ...    use default ionization potentials (from Carlson, et.al.,       ABJT2929
c        Atomic Data (1970)).                                           ABJT2930
                                                                        ABJT2931
         do 260 lgas=1,ngases                                           ABJT2932
                                                                        ABJT2933
            if ( izgas(lgas) .le. 27 ) then                             ABJT2934
               ii = izgas(lgas)                                         ABJT2935
               jj = 0                                                   ABJT2936
            else                                                        ABJT2937
               ii = 55 - izgas(lgas)                                    ABJT2938
               jj = ii                                                  ABJT2939
            endif                                                       ABJT2940
                                                                        ABJT2941
            write (6,901) lgas,izgas(lgas)                              ABJT2942
                                                                        ABJT2943
            do 255 iz=1,izgas(lgas)                                     ABJT2944
               if ( izgas(lgas) .le. 54 ) then                          ABJT2945
                  pot(lgas,iz) = defpot(ii,jj+iz)                       ABJT2946
                  numocc(lgas,iz) = noccdf(ii,jj+iz)                    ABJT2947
               endif                                                    ABJT2948
               write (6,902) iz-1,pot(lgas,iz)                          ABJT2949
               if ( pot(lgas,1).eq.0. ) then                            ABJT2950
                 stop ' No default ionization potentials for this gas'  ABJT2951
               endif                                                    ABJT2952
  255       continue                                                    ABJT2953
                                                                        ABJT2954
  260    continue                                                       ABJT2955
                                                                        ABJT2956
      endif                                                             ABJT2957
                                                                        ABJT2958
c ... open plot files                                                   ABJT2959
                                                                        ABJT2960
      do 450 i=1,30                                                     ABJT2961
         if ( iplot(i).eq.1 ) then                                      ABJT2962
            write (pltfil,952) i                                        ABJT2963
            ip10 = i + 10                                               ABJT2964
            open ( unit=ip10,file=pltfil,status='new',form='formatted' )ABJT2965
         endif                                                          ABJT2966
  450 continue                                                          ABJT2967
                                                                        ABJT2968
c ... print out some variables                                          ABJT2969
                                                                        ABJT2970
      write (6,929) (isw(i),i=1,10)                                     ABJT2971
      write (6,930) (isw(i),i=11,20)                                    ABJT2972
      write (6,932) (con(i),i=1,10)                                     ABJT2973
      write (6,933) (iplot(i),i=1,10)                                   ABJT2974
                                                                        ABJT2975
c ... copy input file to output file                                    ABJT2976
                                                                        ABJT2977
      if ( isw(5) .eq. 0 ) then                                         ABJT2978
         rewind ( unit=5 )                                              ABJT2979
         write (6,800)                                                  ABJT2980
         icard = 0                                                      ABJT2981
  510    icard = icard + 1                                              ABJT2982
           read (5,801,end=520) cstrng                                  ABJT2983
           write (6,802) cstrng                                         ABJT2984
           go to 510                                                    ABJT2985
  520    continue                                                       ABJT2986
         write (6,803)                                                  ABJT2987
      endif                                                             ABJT2988
                                                                        ABJT2989
      return                                                            ABJT2990
                                                                        ABJT2991
c ... format statements                                                 ABJT2992
                                                                        ABJT2993
  800 format (//t4,'Contents of input file -- IONMXINP'//)              ABJT2994
  801 format (a80)                                                      ABJT2995
  802 format (1x,a80)                                                   ABJT2996
  803 format (///)                                                      ABJT2997
  811 format (10(1x,a6))                                                ABJT2998
  812 format (10(1x,i6))                                                ABJT2999
  901 format (//t4,'For gas # ',i2,', with atomic # ',i3,',',/          ABJT3000
     2  t7,'the ionization potentials are:'/)                           ABJT3001
  902 format (t7,'ionization state =',i3,t37,'ionization potential =',  ABJT3002
     2  f10.3,t79)                                                      ABJT3003
  906 format (//t4,'number of gases  = ',i2//t6,'gas #',t26,'atomic #', ABJT3004
     2  t46,'number fraction',t66,'atomic weight'//                     ABJT3005
     3  20(t6,i3,t26,i4,t46,f15.8,t66,f11.3/))                          ABJT3006
  911 format (//t2,a80////)                                             ABJT3007
  929 format (////t30,'Switches used for current calculation'//         ABJT3008
     2  t4,'isw( 1) = ',i3,                                             ABJT3009
     &  '  User supplies ionization potentials (0=>no)'/                ABJT3010
     3  t4,'isw( 2) = ',i3,'  Compute opacities ?  (0=>yes)'/           ABJT3011
     4  t4,'isw( 3) = ',i3,'  Request debug output (# => # of subrts)'/ ABJT3012
     5  t4,'isw( 4) = ',i3,'  Compute heat cap. and dZ/dT ? (0=>yes)'/  ABJT3013
     6  t4,'isw( 5) = ',i3,'  Copy input file to output file (0=>yes)'/ ABJT3014
     7  t4,'isw( 6) = ',i3,'  Saha/Coronal model (0=>S/C, 1=>S, 2=>C,', ABJT3015
     &     ' 3=>C w/ 3-body rec.)'/                                     ABJT3016
     8  t4,'isw( 7) = ',i3,'  '/                                        ABJT3017
     9  t4,'isw( 8) = ',i3,'  Calculation to make a file for CONRAD ?'/ ABJT3018
     &  t17,'    (0=>no; 1,12=>CONRAD; 2,12=>SESAME)'/                  ABJT3019
     6  t4,'isw( 9) = ',i3,'  Maximum prin. quantum # in populations'   ABJT3020
     &  ' calculation (0=>code picks)'/                                 ABJT3021
     &  t17,'     (if>0,npmaxp=isw9+grnd.st.#)'/                        ABJT3022
     1  t4,'isw(10) = ',i3,'  ')                                        ABJT3023
  930 format(                                                           ABJT3024
     2  t4,'isw(11) = ',i3,'  Include delta-n=0 transitions ?(0=>yes)'/ ABJT3025
     3  t4,'isw(12) = ',i3,'  Restrict ions to ground state ? (0=>no)'/ ABJT3026
     4  t4,'isw(13) = ',i3,'  Specify group boundaries (0=>default)'/   ABJT3027
     &  t17,'     (1=>user specifies)  (2=>default T-dep. values)'/     ABJT3028
     5  t4,'isw(14) = ',i3,'  Voigt or Lorentzian line profile (0=>V)'/ ABJT3029
     6  t4,'isw(15) = ',i3,'  Maximum prin. quantum # in opacity'       ABJT3030
     &  ' calculations (0=>code picks)'/                                ABJT3031
     &  t17,'     (if<0,npqmax=-isw15; if>0,npqmax=isw15+grnd.st.#)'/   ABJT3032
     7  t4,'isw(16) = ',i3,'  Include dielectronic recomb. ? (0=>yes)'/ ABJT3033
     8  t4,'isw(17) = ',i3,'  Include bremsstrahlung ? (0=>yes)'/       ABJT3034
     9  t4,'isw(18) = ',i3,'  Include photoionization ? (0=>yes)'/      ABJT3035
     &  t4,'isw(19) = ',i3,'  Include line contributions ? (0=>yes)'/   ABJT3036
     &  t17,'     (1=>no)    (2=>core/wings computed separately)'/      ABJT3037
     1  t4,'isw(20) = ',i3,'  Include scattering contribs. ? (0=>yes)') ABJT3038
  931 format (                                                          ABJT3039
     2  t4,'isw(21) = ',i3,'  '/                                        ABJT3040
     3  t4,'isw(22) = ',i3,'  '/                                        ABJT3041
     4  t4,'isw(23) = ',i3,'  '/                                        ABJT3042
     5  t4,'isw(24) = ',i3,'  '/                                        ABJT3043
     6  t4,'isw(25) = ',i3,'  '/,                                       ABJT3044
     7  t4,'isw(26) = ',i3,'  '/                                        ABJT3045
     8  t4,'isw(27) = ',i3,'  '/                                        ABJT3046
     9  t4,'isw(28) = ',i3,'  '/                                        ABJT3047
     &  t4,'isw(29) = ',i3,'  '/                                        ABJT3048
     1  t4,'isw(30) = ',i3,'  '/)                                       ABJT3049
  932 format (////t30,'Constants used for current calculation'//        ABJT3050
     2  t4,'con( 1) = ',1pe10.2,'  '/                                   ABJT3051
     3  t4,'con( 2) = ',1pe10.2,                                        ABJT3052
     &  '  min. species concentration to compute bb and bf transitions'/ABJT3053
     4  t4,'con( 3) = ',1pe10.2,'  min. ionization concentration '      ABJT3054
     &  'to compute bb and bf transitions'/                             ABJT3055
     5  t4,'con( 4) = ',1pe10.2,'  min. concentration of ',             ABJT3056
     &  'an atomic state to consider bb and bf transitions'/            ABJT3057
     6  t4,'con( 5) = ',1pe10.2,'  range, in # of line widths (fwhm),'  ABJT3058
     &  ' to compute contribution from a bb transition'/                ABJT3059
     7  t4,'con( 6) = ',1pe10.2,'  width, in # of line widths (fwhm),'  ABJT3060
     &  ' of line core'/                                                ABJT3061
     8  t4,'con( 7) = ',1pe10.2,'  abs. cfs. are weighted by Planck '   ABJT3062
     &  'fn. when 1/con7 < hv/kT < con7'/                               ABJT3063
     9  t4,'con( 8) = ',1pe10.2,'  '/                                   ABJT3064
     &  t4,'con( 9) = ',1pe10.2,'  multiplier for mesh pt spacing'      ABJT3065
     &  ' about lines'/                                                 ABJT3066
     1  t4,'con(10) = ',1pe10.2,'  '/)                                  ABJT3067
  933 format (////t30,'Plot files opened for current calculation'//     ABJT3068
     2  t4,'iplot( 1) = ',i3,'  Absorption coefs. vs. photon energy'/   ABJT3069
     3  t4,'iplot( 2) = ',i3,'  Mean opacities vs. temperature'/        ABJT3070
     4  t4,'iplot( 3) = ',i3,'  Emission coefs. vs. photon energy'/     ABJT3071
     5  t4,'iplot( 4) = ',i3,'  Mean opacities vs. density'/            ABJT3072
     6  t4,'iplot( 5) = ',i3,'  Charge and cool. rate vs. temp.'/       ABJT3073
     7  t4,'iplot( 6) = ',i3,'  '/                                      ABJT3074
     8  t4,'iplot( 7) = ',i3,'  '/                                      ABJT3075
     9  t4,'iplot( 8) = ',i3,'  '/                                      ABJT3076
     &  t4,'iplot( 9) = ',i3,'  Ioniz. pop.s vs. temperature'/          ABJT3077
     1  t4,'iplot(10) = ',i3,'  '/)                                     ABJT3078
  941 format (///t2,'Calculation stopped in -input-.'/t2,               ABJT3079
     2  'This run is supposed to write a file in CONRAD format, but'/   ABJT3080
     3  t2,'the input data is not consistent with this format.'//       ABJT3081
     4  t2,'The input variables are:'/t4,'ntemp',t18,'ndens',           ABJT3082
     5  t32,'ntrad',t46,'dlgtmp',t60,'dlgden',t74,'dlgtrd'/             ABJT3083
     6  t4,3(i4,10x),1p3e14.4)                                          ABJT3084
  942 format (/t2,'They should be:'/t4,'ntemp',t18,'ndens',             ABJT3085
     2  t32,'ntrad',t46,'dlgtmp',t60,'dlgden',t74,'dlgtrd'/             ABJT3086
     3  t4,3(i4,10x),3('   .gt. 0.',4x))                                ABJT3087
  943 format (///t2,'Calculation stopped in -input-.'/t2,               ABJT3088
     2  'This run is supposed to write a file in the NEW CONRAD',       ABJT3089
     &  ' format, but'/                                                 ABJT3090
     3  t2,'the input data is not consistent with this format.'//       ABJT3091
     4  t2,'The input variables are:'/t4,'ntemp',t18,'ndens',           ABJT3092
     5  t32,'dlgtmp',t46,'dlgden'/                                      ABJT3093
     6  t4,2(i4,10x),1p2e14.4)                                          ABJT3094
  944 format (/t2,'They should be:'/t4,'ntemp',t18,'ndens',             ABJT3095
     2  t32,'dlgtmp',t46,'dlgden'/                                      ABJT3096
     3  t4,2(i4,10x),2('   .gt. 0.',4x))                                ABJT3097
  951 format ('atom',i2.2)                                              ABJT3098
  952 format ('implot',i2.2)                                            ABJT3099
  955 format (t10,' IONMIX calculation performed on  ',a4,i3,', 19',    ABJT3100
     2  i2,'  at  ',a8)                                                 ABJT3101
  956 format (t10,' IONMIX calculation performed on date  ',a8,         ABJT3102
     2  '  at time  ',a8)                                               ABJT3103
                                                                        ABJT3104
      end                                                               ABJT3105
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT3106
                                                                        ABJT3107
      subroutine lines ( tp,densnn,densne,ephot,                        ABJT3108
     &                                            abstot,emstot )       ABJT3109
c                                                                       ABJT3110
c ... compute the contribution to the absorption coefficient from       ABJT3111
c     all lines (bound-bound transitions)                               ABJT3112
c                                                                       ABJT3113
c ... input variables:                                                  ABJT3114
c       tp      -  plasma temperature (eV)                              ABJT3115
c       densnn  -  number density of all nuclei (cm**-3)                ABJT3116
c       densne  -  electron density (cm**-3)                            ABJT3117
c       ephot   -  photon energy (eV)                                   ABJT3118
c                                                                       ABJT3119
c ... output variables:                                                 ABJT3120
c       abstot  -  absorption coefficient due to lines (cm**-1)         ABJT3121
c       emstot  -  emission coefficient due to lines (cm**-1)           ABJT3122
c                                                                       ABJT3123
                                                                        ABJT3124
c ... set the maximum number of:                                        ABJT3125
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT3126
c       temperatures (mxtemp), densities (mxdens),                      ABJT3127
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT3128
                                                                        ABJT3129
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT3130
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT3131
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT3132
                                                                        ABJT3133
                                                                        ABJT3134
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT3135
c ................................................................      ABJT3136
                                                                        ABJT3137
      common / gases / ngases,izgas(mxgass),atomwt(mxgass),             ABJT3138
     2                 avgatw,avgatn,                                   ABJT3139
     3                 fracsp(mxgass),fraciz(mxgass,mxatom),            ABJT3140
     4                 fraclv(mxgass,mxatom,20),pot(mxgass,mxatom),     ABJT3141
     5                 numocc(mxgass,mxatom)                            ABJT3142
c ................................................................      ABJT3143
                                                                        ABJT3144
      common / params / npqmax(mxgass), npmaxp, nfrqbb,                 ABJT3145
     2                  npring(100), neopen(5,100), bbnton(6,55),       ABJT3146
     3                  defpot(27,55), noccdf(27,55)                    ABJT3147
c ................................................................      ABJT3148
                                                                        ABJT3149
      logical lbug                                                      ABJT3150
                                                                        ABJT3151
c ********************************************************************  ABJT3152
c                                                                       ABJT3153
c                         begin execution                               ABJT3154
c                                                                       ABJT3155
c ********************************************************************  ABJT3156
                                                                        ABJT3157
c ... check debug option                                                ABJT3158
      call debug ( 'lines',                                             ABJT3159
     &                        lbug )                                    ABJT3160
                                                                        ABJT3161
      if ( lbug ) write (6,900) ephot,tp,densnn                         ABJT3162
                                                                        ABJT3163
      isav = 0                                                          ABJT3164
      abstot = 0.                                                       ABJT3165
      emstot = 0.                                                       ABJT3166
                                                                        ABJT3167
c ... loop over gas species                                             ABJT3168
      do 400 lgas=1,ngases                                              ABJT3169
         if ( fracsp(lgas) .lt. con(2) ) go to 400                      ABJT3170
                                                                        ABJT3171
c ...    loop over ionization states                                    ABJT3172
         do 300 izlevl=0,izgas(lgas)-1                                  ABJT3173
            izp1 = izlevl + 1                                           ABJT3174
            if ( fracsp(lgas)*fraciz(lgas,izp1) .lt. con(3)) go to 300  ABJT3175
                                                                        ABJT3176
c ...       find the principal quantum number of the valence electrons  ABJT3177
c           in their ground state                                       ABJT3178
            nbound = izgas(lgas) - izlevl                               ABJT3179
            nprin0 = npring(nbound)                                     ABJT3180
                                                                        ABJT3181
            denlq  = densnn * fracsp(lgas) * fraciz(lgas,izp1)          ABJT3182
                                                                        ABJT3183
c ...       loop over initial quantum states                            ABJT3184
            do 200 n=nprin0,npqmax(lgas)-1                              ABJT3185
                                                                        ABJT3186
               if ( fracsp(lgas)*fraciz(lgas,izp1)*fraclv(lgas,izp1,n)  ABJT3187
     &              .lt. con(4) ) go to 200                             ABJT3188
                                                                        ABJT3189
c ...          compute the density of ions with an electron in          ABJT3190
c              in state "n"                                             ABJT3191
               denlqn = denlq * fraclv(lgas,izp1,n)                     ABJT3192
                                                                        ABJT3193
c ...          loop over final quantum states                           ABJT3194
               do 100 np=n,npqmax(lgas)                                 ABJT3195
                                                                        ABJT3196
                  dnlqnp = denlq * fraclv(lgas,izp1,np)                 ABJT3197
                  call abslin ( tp,densnn,densne,denlqn,dnlqnp,ephot,   ABJT3198
     &                          pot(lgas,izp1),atomwt(lgas),n,np,       ABJT3199
     &                          nprin0,nbound,izgas(lgas),              ABJT3200
     &                                                 abscof,emscof )  ABJT3201
                                                                        ABJT3202
                  if ( abscof .ne. 0. ) then                            ABJT3203
                     abstot = abstot + abscof                           ABJT3204
                     emstot = emstot + emscof                           ABJT3205
                     isav = isav + 1                                    ABJT3206
                     if ( isav.eq.200 ) then                            ABJT3207
                        stop ' you are keeping track of too many lines' ABJT3208
                     endif                                              ABJT3209
                  endif                                                 ABJT3210
                                                                        ABJT3211
  100          continue                                                 ABJT3212
  200       continue                                                    ABJT3213
  300    continue                                                       ABJT3214
  400 continue                                                          ABJT3215
                                                                        ABJT3216
      if ( lbug ) write (6,905) abstot,emstot,isav                      ABJT3217
                                                                        ABJT3218
      return                                                            ABJT3219
                                                                        ABJT3220
c ... format statements                                                 ABJT3221
                                                                        ABJT3222
  900 format (t2,'debug output from -lines-'/t4,'ephot',t18,'tp',       ABJT3223
     2  t32,'densnn'/t4,1p3e14.4)                                       ABJT3224
  905 format (t4,'abstot',t18,'emstot',t32,'# of lines contributing'    ABJT3225
     2  /t4,1p2e14.4,t40,0p,i4)                                         ABJT3226
                                                                        ABJT3227
      end                                                               ABJT3228
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT3229
                                                                        ABJT3230
      subroutine linwid ( tp,densnn,ennp,atomwt,                        ABJT3231
     &                                         gamma,avoigt,dnudop )    ABJT3232
c    MOVED TO ModOpacityVoigt                                           ABJT3233
c ... compute the total line width of a bound-bound transition.  this   ABJT3234
c     sums the contributions from natural, Doppler, and collisional     ABJT3235
c     broadening.                                                       ABJT3236
c                                                                       ABJT3237
c ... input variables:                                                  ABJT3238
c       tp      -  plasma temperature (eV)                              ABJT3239
c       densnn  -  number density of all nuclei (cm**-3)                ABJT3240
c       ennp    -  energy of the transition from "n" to "np" (eV)       ABJT3241
c       atomwt  -  atomic weight of the ion (amu)                       ABJT3242
c                                                                       ABJT3243
c ... output variables:                                                 ABJT3244
c       gamma   -  total line width (sec**-1)                           ABJT3245
c       avoigt  -  "a" used to compute the Voigt profile                ABJT3246
c       dnudop  -  doppler width (sec**-1)                              ABJT3247
c                                                                       ABJT3248
      logical lbug                                                      ABJT3249
                                                                        ABJT3250
c ********************************************************************  ABJT3251
c                                                                       ABJT3252
c                         begin execution                               ABJT3253
c                                                                       ABJT3254
c ********************************************************************  ABJT3255
                                                                        ABJT3256
c ... check debug option                                                ABJT3257
      call debug ( 'linwid',                                            ABJT3258
     &                        lbug )                                    ABJT3259
                                                                        ABJT3260
      vel = sqrt( tp / atomwt )                                         ABJT3261
                                                                        ABJT3262
      widnat = 2.29e6 * ennp**2                                         ABJT3263
      widdop = 1.41e11 * ennp * vel                                     ABJT3264
      widcol = 4.58e6 * densnn**0.3333 * vel                            ABJT3265
                                                                        ABJT3266
      gamma = widnat + widdop + widcol                                  ABJT3267
      avoigt = ( widnat + widcol ) / widdop                             ABJT3268
      dnudop = widdop / 12.566                                          ABJT3269
                                                                        ABJT3270
      if ( lbug ) write (6,901) tp,densnn,ennp,gamma,widnat,widdop,widcoABJT3271
                                                                        ABJT3272
      return                                                            ABJT3273
                                                                        ABJT3274
c ... format statements                                                 ABJT3275
                                                                        ABJT3276
  901 format (t2,'debug output from -linwid-'/t4,'tp',t18,'densnn',     ABJT3277
     &  t32,'ennp',t46,'gamma',t60,'widnat',t74,'widdop',t88,'widcol'/  ABJT3278
     &  t4,1p7e14.4)                                                    ABJT3279
                                                                        ABJT3280
      end                                                               ABJT3281
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT3282
                                                                        ABJT3283
      subroutine meshhv ( tp,densnn,densne,engrup,ngrups,nptspg,        ABJT3284
     &                                                  photen,nphot )  ABJT3285
c  COPIED TO ModMultiGroup.f90                                          ABJT3286
c ... set a mesh of photon energies at which we would like to evaluate  ABJT3287
c     absorption coefficients                                           ABJT3288
c                                                                       ABJT3289
c ... input variables:                                                  ABJT3290
c       tp      -  plasma temperature (eV)                              ABJT3291
c       densnn  -  number density of all nuclei (cm**-3)                ABJT3292
c       densne  -  electron density (cm**-3)                            ABJT3293
c       engrup  -  photon energy group boundaries (ngrups+1) (eV)       ABJT3294
c       ngrups  -  number of photon energy groups                       ABJT3295
c       nptspg  -  minimum number of mesh points per energy group       ABJT3296
c                                                                       ABJT3297
c ... output variables:                                                 ABJT3298
c       photen  -  array of photon energies (eV)                        ABJT3299
c       nphot   -  number of mesh points including contrbutions from    ABJT3300
c                  lines and photoionization edges                      ABJT3301
c                                                                       ABJT3302
                                                                        ABJT3303
c ... set the maximum number of:                                        ABJT3304
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT3305
c       temperatures (mxtemp), densities (mxdens),                      ABJT3306
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT3307
                                                                        ABJT3308
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT3309
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT3310
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT3311
                                                                        ABJT3312
                                                                        ABJT3313
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT3314
c ................................................................      ABJT3315
                                                                        ABJT3316
      common / gases / ngases,izgas(mxgass),atomwt(mxgass),             ABJT3317
     2                 avgatw,avgatn,                                   ABJT3318
     3                 fracsp(mxgass),fraciz(mxgass,mxatom),            ABJT3319
     4                 fraclv(mxgass,mxatom,20),pot(mxgass,mxatom),     ABJT3320
     5                 numocc(mxgass,mxatom)                            ABJT3321
c ................................................................      ABJT3322
                                                                        ABJT3323
      common / params / npqmax(mxgass), npmaxp, nfrqbb,                 ABJT3324
     2                  npring(100), neopen(5,100), bbnton(6,55),       ABJT3325
     3                  defpot(27,55), noccdf(27,55)                    ABJT3326
c ................................................................      ABJT3327
                                                                        ABJT3328
c ... type declaration statements                                       ABJT3329
      dimension photen(mxphot),engrup(ngrups+1)                         ABJT3330
                                                                        ABJT3331
      logical lbug                                                      ABJT3332
                                                                        ABJT3333
      data dplus / 1.001 /, dminus / 0.999 /                            ABJT3334
                                                                        ABJT3335
c ********************************************************************  ABJT3336
c                                                                       ABJT3337
c                         begin execution                               ABJT3338
c                                                                       ABJT3339
c ********************************************************************  ABJT3340
                                                                        ABJT3341
c ... check debug option                                                ABJT3342
      call debug ( 'meshhv',                                            ABJT3343
     &                        lbug )                                    ABJT3344
                                                                        ABJT3345
c ... initialize variables                                              ABJT3346
                                                                        ABJT3347
      nphot = 0                                                         ABJT3348
      do 5 i=1,mxphot                                                   ABJT3349
         photen(i) = 0.                                                 ABJT3350
    5 continue                                                          ABJT3351
                                                                        ABJT3352
      if ( lbug ) write (6,901) nptspg,ngrups,tp,densnn,densne          ABJT3353
                                                                        ABJT3354
c ... set up initial grid within each photon energy group               ABJT3355
                                                                        ABJT3356
      do 100 ig=1,ngrups                                                ABJT3357
                                                                        ABJT3358
         hvmin = engrup(ig)                                             ABJT3359
         hvmax = engrup(ig+1)                                           ABJT3360
         dloghv = log( hvmax/hvmin ) / nptspg                           ABJT3361
         hvlog = log( hvmin ) - dloghv                                  ABJT3362
         if ( ig.lt.ngrups ) then                                       ABJT3363
            nmax = nptspg                                               ABJT3364
         else                                                           ABJT3365
            nmax = nptspg + 1                                           ABJT3366
         endif                                                          ABJT3367
                                                                        ABJT3368
         do 60 i=1,nmax                                                 ABJT3369
            hvlog = hvlog + dloghv                                      ABJT3370
            photen(nphot+1) = exp( hvlog )                              ABJT3371
            nphot = nphot + 1                                           ABJT3372
   60    continue                                                       ABJT3373
                                                                        ABJT3374
  100 continue                                                          ABJT3375
                                                                        ABJT3376
c ... add 2 points near the plasma cutoff frequency                     ABJT3377
                                                                        ABJT3378
      hvcut = sqrt( densne / 7.25e20 )                                  ABJT3379
      if ( hvcut.gt.engrup(1) .and. hvcut.lt.engrup(ngrups+1) ) then    ABJT3380
         photen(nphot+1) = hvcut * dminus                               ABJT3381
         photen(nphot+2) = hvcut * dplus                                ABJT3382
         nphot = nphot + 2                                              ABJT3383
      endif                                                             ABJT3384
                                                                        ABJT3385
c ... add points near line centers and ionization edges                 ABJT3386
                                                                        ABJT3387
      do 300 lgas=1,ngases                                              ABJT3388
         if ( fracsp(lgas) .lt. con(2) ) go to 300                      ABJT3389
                                                                        ABJT3390
         do 290 izlevl=0,izgas(lgas)-1                                  ABJT3391
            izp1 = izlevl + 1                                           ABJT3392
            if ( fracsp(lgas)*fraciz(lgas,izp1) .lt. con(3) ) go to 290 ABJT3393
                                                                        ABJT3394
c ...       find the principal quantum number of the valence electrons  ABJT3395
c           in their ground state                                       ABJT3396
            nbound = izgas(lgas) - izlevl                               ABJT3397
            nprin0 = npring(nbound)                                     ABJT3398
                                                                        ABJT3399
            do 280 n=nprin0,npqmax(lgas)                                ABJT3400
               if ( fracsp(lgas)*fraciz(lgas,izp1)*fraclv(lgas,izp1,n)  ABJT3401
     &              .lt. con(4) ) go to 275                             ABJT3402
                                                                        ABJT3403
               if ( isw(19).ne.2 .and. n.lt.npqmax(lgas) ) then         ABJT3404
                                                                        ABJT3405
                  do 270 np=n,npqmax(lgas)                              ABJT3406
c ...               calculate the energy of the transition              ABJT3407
                    if ( n .eq. np ) then                               ABJT3408
                       call bbneq0 ( tp,izgas(lgas),nbound,             ABJT3409
     &                                                ennp,dum1,dum2 )  ABJT3410
                       if ( ennp.le.0. .or. isw(11).ne.0 ) go to 270    ABJT3411
                    else                                                ABJT3412
                       ennp = nprin0**2 * pot(lgas,izp1) *              ABJT3413
     &                        (1./n**2-1./np**2)                        ABJT3414
                    endif                                               ABJT3415
                    if ( ennp .gt. engrup(1) .and.                      ABJT3416
     &                   ennp .lt. engrup(ngrups+1) ) then              ABJT3417
                                                                        ABJT3418
                       call linwid ( tp,densnn,ennp,atomwt(lgas),       ABJT3419
     &                                            gamma,avoigt,dnudop ) ABJT3420
                                                                        ABJT3421
c ...                  set points at and near the line center           ABJT3422
                       dhnu = con(9) * 4.14e-15 * gamma / 12.57         ABJT3423
                       photen(nphot+1) = ennp                           ABJT3424
                       nphot = nphot + 1                                ABJT3425
                       ilast = nfrqbb / 2 + 1                           ABJT3426
                       do 260 i=1,ilast                                 ABJT3427
                         kk = 2**(i-1)                                  ABJT3428
                         photen(nphot+1) = ennp - kk * dhnu             ABJT3429
                         photen(nphot+2) = ennp + kk * dhnu             ABJT3430
                         nphot = nphot + 2                              ABJT3431
  260                  continue                                         ABJT3432
                                                                        ABJT3433
c ...                  stop if too many mesh pts are requested          ABJT3434
                       if ( nphot .gt. mxphot - 5 - 2*ilast )           ABJT3435
     &                    stop ' too many pts in -meshhv-'              ABJT3436
                     endif                                              ABJT3437
                                                                        ABJT3438
  270             continue                                              ABJT3439
                                                                        ABJT3440
               endif                                                    ABJT3441
                                                                        ABJT3442
c ...          add 2 points near the ionization edge                    ABJT3443
  275          edge = nprin0**2 * pot(lgas,izp1) / n**2                 ABJT3444
               if ( edge.gt.engrup(1) .and. edge.lt.engrup(ngrups+1) )  ABJT3445
     &            then                                                  ABJT3446
                  photen(nphot+1) = edge * dminus                       ABJT3447
                  photen(nphot+2) = edge * dplus                        ABJT3448
                  nphot = nphot + 2                                     ABJT3449
               endif                                                    ABJT3450
                                                                        ABJT3451
  280       continue                                                    ABJT3452
                                                                        ABJT3453
  290    continue                                                       ABJT3454
                                                                        ABJT3455
c ...    add 2 points near photoionization edges of core electrons      ABJT3456
         do 295 ncore=1,izgas(lgas)-1                                   ABJT3457
            edge = pot(lgas,izgas(lgas)+1-ncore)                        ABJT3458
            if ( edge.gt.engrup(1) .and.                                ABJT3459
     &           edge.lt.engrup(ngrups+1) ) then                        ABJT3460
               photen(nphot+1) = edge * dminus                          ABJT3461
               photen(nphot+2) = edge * dplus                           ABJT3462
               nphot = nphot + 2                                        ABJT3463
            endif                                                       ABJT3464
  295    continue                                                       ABJT3465
                                                                        ABJT3466
  300 continue                                                          ABJT3467
                                                                        ABJT3468
c ... now, arrange the photon energies in monotonically increasing      ABJT3469
c     order                                                             ABJT3470
      call sort ( nphot,                                                ABJT3471
     &                   photen )                                       ABJT3472
                                                                        ABJT3473
      if ( lbug ) write (6,902) nphot,(photen(i),i=1,nphot)             ABJT3474
                                                                        ABJT3475
      return                                                            ABJT3476
                                                                        ABJT3477
c ... format statements                                                 ABJT3478
                                                                        ABJT3479
  901 format ('  debug output from -meshhv-:'/t4,'nptspg',t18,          ABJT3480
     2  'ngrups',t32,'tp',t46,'densnn',t60,'densne'/t4,i4,t18,          ABJT3481
     3  i4,t32,1p3e11.3)                                                ABJT3482
  902 format ('  nphot =',i5,4x,'photen:'/50(1p10e11.3/))               ABJT3483
c END OF THE MOVED FRAGMENT                                             ABJT3484
      end                                                               ABJT3485
c COPIED TO ModMultiGroup
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT3486
                                                                        ABJT3487
      subroutine opacbb ( tp,densnn,densne,trad,xgmin,xgmax,            ABJT3488
     &                                            opabs,opems,opross )  ABJT3489
c                                                                       ABJT3490
c ... compute the contribution to the group opacity from                ABJT3491
c     all lines (bound-bound transitions)                               ABJT3492
c                                                                       ABJT3493
c ... input variables:                                                  ABJT3494
c       tp      -  plasma temperature (eV)                              ABJT3495
c       densnn  -  number density of all nuclei (cm**-3)                ABJT3496
c       densne  -  electron density (cm**-3)                            ABJT3497
c       trad    -  radiation temperature (eV)                           ABJT3498
c       xgmin   -  minimum photon energy of group (in units of kT)      ABJT3499
c       xgmax   -  maximum photon energy of group (in units of kT)      ABJT3500
c                                                                       ABJT3501
c ... output variables:                                                 ABJT3502
c       opabs   -  Planck absorption due to lines (cm**2/g)             ABJT3503
c       opems   -  Planck emission due to lines (cm**2/g)               ABJT3504
c       opross  -  Rosseland opacity due to lines (cm**2/g)             ABJT3505
c                                                                       ABJT3506
                                                                        ABJT3507
c ... set the maximum number of:                                        ABJT3508
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT3509
c       temperatures (mxtemp), densities (mxdens),                      ABJT3510
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT3511
                                                                        ABJT3512
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT3513
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT3514
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT3515
                                                                        ABJT3516
                                                                        ABJT3517
      common / consts / pi, avgdro, sbcon, hplank                       ABJT3518
c ...............................................................       ABJT3519
                                                                        ABJT3520
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT3521
c ................................................................      ABJT3522
                                                                        ABJT3523
      common / gases / ngases,izgas(mxgass),atomwt(mxgass),             ABJT3524
     2                 avgatw,avgatn,                                   ABJT3525
     3                 fracsp(mxgass),fraciz(mxgass,mxatom),            ABJT3526
     4                 fraclv(mxgass,mxatom,20),pot(mxgass,mxatom),     ABJT3527
     5                 numocc(mxgass,mxatom)                            ABJT3528
c ................................................................      ABJT3529
                                                                        ABJT3530
      common / params / npqmax(mxgass), npmaxp, nfrqbb,                 ABJT3531
     2                  npring(100), neopen(5,100), bbnton(6,55),       ABJT3532
     3                  defpot(27,55), noccdf(27,55)                    ABJT3533
c ................................................................      ABJT3534
                                                                        ABJT3535
c ... type declaration statements                                       ABJT3536
      logical lbug                                                      ABJT3537
                                                                        ABJT3538
c ... statement function (oscillator strength taken from Zeldovich      ABJT3539
c     & Raizor for the upward transition from "ni" to "nf")             ABJT3540
      oscstr(ni,nf) = 1.96 / ni**5 / nf**3 / (1./ni**2-1./nf**2)**3     ABJT3541
                                                                        ABJT3542
c ********************************************************************  ABJT3543
c                                                                       ABJT3544
c                         begin execution                               ABJT3545
c                                                                       ABJT3546
c ********************************************************************  ABJT3547
                                                                        ABJT3548
c ... check debug option                                                ABJT3549
      call debug ( 'opacbb',                                            ABJT3550
     &                        lbug )                                    ABJT3551
                                                                        ABJT3552
      if ( lbug ) write (6,900) tp,densnn,trad,xgmin,xgmax              ABJT3553
                                                                        ABJT3554
      opabs = 0.                                                        ABJT3555
      opems = 0.                                                        ABJT3556
      opross = 0.                                                       ABJT3557
                                                                        ABJT3558
      rho = densnn * avgatw / avgdro                                    ABJT3559
      g5  = gint ( 5,xgmin,xgmax )                                      ABJT3560
                                                                        ABJT3561
      if ( g5 .eq. 0. ) return                                          ABJT3562
                                                                        ABJT3563
      const = 1.10e-16 * densnn / ( tp * rho * g5 )                     ABJT3564
                                                                        ABJT3565
c ... loop over gas species                                             ABJT3566
      do 400 lgas=1,ngases                                              ABJT3567
         if ( fracsp(lgas) .lt. con(2) ) go to 400                      ABJT3568
                                                                        ABJT3569
c ...    loop over ionization states                                    ABJT3570
         do 300 izlevl=0,izgas(lgas)-1                                  ABJT3571
            izp1 = izlevl + 1                                           ABJT3572
            if ( fracsp(lgas)*fraciz(lgas,izp1) .lt. con(3)) go to 300  ABJT3573
                                                                        ABJT3574
c ...       find the principal quantum number of the valence electrons  ABJT3575
c           in their ground state                                       ABJT3576
            nbound = izgas(lgas) - izlevl                               ABJT3577
            nprin0 = npring(nbound)                                     ABJT3578
                                                                        ABJT3579
            const1 = const * fracsp(lgas) * fraciz(lgas,izp1)           ABJT3580
                                                                        ABJT3581
c ...       loop over initial quantum states                            ABJT3582
            do 200 n=nprin0,npqmax(lgas)-1                              ABJT3583
                                                                        ABJT3584
               if ( fracsp(lgas)*fraciz(lgas,izp1)*fraclv(lgas,izp1,n)  ABJT3585
     &              .lt. con(4) ) go to 200                             ABJT3586
                                                                        ABJT3587
c ...          calculate opacity due to "delta n" = 0 transitions       ABJT3588
                                                                        ABJT3589
               call bbneq0 ( tp,izgas(lgas),nbound,                     ABJT3590
     &                                             enn,fnn,gnn )        ABJT3591
               if ( isw(11).ne.0 ) enn = 0.                             ABJT3592
               x0 = enn / trad                                          ABJT3593
               if ( x0.gt.xgmin .and. x0.le.xgmax ) then                ABJT3594
                  ex0 = exp( -x0 )                                      ABJT3595
                  opacnn = const1 * fnn * fraclv(lgas,izp1,n) *         ABJT3596
     &                     ex0 * x0**3                                  ABJT3597
                                                                        ABJT3598
c ...             correct for the "effective" stimulated emission to    ABJT3599
c                 give the proper form for the cooling rate             ABJT3600
                  xnjoni = gnn * densne * ex0 /                         ABJT3601
     &                   ( gnn * densne + enn**3 * sqrt(tp) * 2.74e12 ) ABJT3602
                  corsee = xnjoni / ex0                                 ABJT3603
                  corsea = ( 1.-xnjoni ) / ( 1.-ex0 )                   ABJT3604
                  opabs = opabs + opacnn * corsea                       ABJT3605
                  opems = opems + opacnn * corsee                       ABJT3606
               endif                                                    ABJT3607
                                                                        ABJT3608
c ...          loop over final quantum states                           ABJT3609
               do 100 np=n+1,npqmax(lgas)                               ABJT3610
                                                                        ABJT3611
c ...            compute the transition energy                          ABJT3612
                 ennp = nprin0**2 * pot(lgas,izp1) * (1./n**2-1./np**2) ABJT3613
                 x0 = ennp / trad                                       ABJT3614
                 if ( x0 .le. xgmin .or. x0 .gt. xgmax ) go to 100      ABJT3615
                                                                        ABJT3616
                 ex0 = exp( -x0 )                                       ABJT3617
                 fnnp = oscstr(n,np)                                    ABJT3618
                 opacij = const1 * fnnp * fraclv(lgas,izp1,n) *         ABJT3619
     &                    ex0 * x0**3                                   ABJT3620
                                                                        ABJT3621
c ...            note that in LTE, "dum" = "ex0"                        ABJT3622
                 dum = ( n*n*fraclv(lgas,izp1,np) ) /                   ABJT3623
     &                 ( np*np*fraclv(lgas,izp1,n) )                    ABJT3624
                 opaca = opacij * ( 1.-dum ) / ( 1.-ex0 )               ABJT3625
                 opace = opacij * dum / ex0                             ABJT3626
                                                                        ABJT3627
                 opabs = opabs + opaca                                  ABJT3628
                 opems = opems + opace                                  ABJT3629
                                                                        ABJT3630
  100          continue                                                 ABJT3631
  200       continue                                                    ABJT3632
  300    continue                                                       ABJT3633
  400 continue                                                          ABJT3634
                                                                        ABJT3635
      if ( lbug ) write (6,905) opabs,opems                             ABJT3636
                                                                        ABJT3637
      return                                                            ABJT3638
                                                                        ABJT3639
c ... format statements                                                 ABJT3640
                                                                        ABJT3641
  900 format (t2,'debug output from -opacbb-'/t4,'tp',t18,'densnn',     ABJT3642
     2  t32,'trad',t46,'xgmin',t60,'xgmax'/t4,1p5e14.4)                 ABJT3643
  905 format (t4,'opabs',t18,'opems'/t4,1p2e14.4)                       ABJT3644
                                                                        ABJT3645
      end                                                               ABJT3646
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT3647
                                                                        ABJT3648
      subroutine opacgp ( tp,rho,trad,photen,nphot,abscfs,emscfs,       ABJT3649
     &                    sctcfs,xgmin,xgmax,                           ABJT3650
     &                                            opacpa,opacpe,opacr ) ABJT3651
                                                                        ABJT3652
c ... this routine computes the Planck and Rosseland opacities for      ABJT3653
c     the photon energy group from "xgmin" to "xgmax".                  ABJT3654
c                                                                       ABJT3655
c ... input variables:                                                  ABJT3656
c       tp     - plasma temperature (eV)                                ABJT3657
c       rho    - mass density (g/cm**3)                                 ABJT3658
c       trad   - radiation temperature (eV)                             ABJT3659
c       photen - photon energies at which the absorption coefficients   ABJT3660
c                were calculated                                        ABJT3661
c       nphot  - number of photon energy mesh points                    ABJT3662
c       abscfs - array of absorption coefficients (cm**-1) (w/o scatt)  ABJT3663
c       emscfs - array of emission coefficients (cm**-1) (w/o scatt)    ABJT3664
c       sctcfs - array of scattering coefficients (cm**-1)              ABJT3665
c       xgmin  - minimum photon energy (in units of kT)                 ABJT3666
c       xgmax  - maximum photon energy (in units of kT)                 ABJT3667
c                                                                       ABJT3668
c ... output variables:                                                 ABJT3669
c       opacpa - Planck opacity for absorption (cm**-1)                 ABJT3670
c       opacpe - Planck opacity for emission (cm**-1)                   ABJT3671
c       opacr  - Rosseland opacity (cm**-1)                             ABJT3672
c                                                                       ABJT3673
                                                                        ABJT3674
c ... set the maximum number of:                                        ABJT3675
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT3676
c       temperatures (mxtemp), densities (mxdens),                      ABJT3677
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT3678
                                                                        ABJT3679
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT3680
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT3681
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT3682
                                                                        ABJT3683
                                                                        ABJT3684
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT3685
c ................................................................      ABJT3686
                                                                        ABJT3687
      common / gases / ngases,izgas(mxgass),atomwt(mxgass),             ABJT3688
     2                 avgatw,avgatn,                                   ABJT3689
     3                 fracsp(mxgass),fraciz(mxgass,mxatom),            ABJT3690
     4                 fraclv(mxgass,mxatom,20),pot(mxgass,mxatom),     ABJT3691
     5                 numocc(mxgass,mxatom)                            ABJT3692
c ................................................................      ABJT3693
                                                                        ABJT3694
      common / params / npqmax(mxgass), npmaxp, nfrqbb,                 ABJT3695
     2                  npring(100), neopen(5,100), bbnton(6,55),       ABJT3696
     3                  defpot(27,55), noccdf(27,55)                    ABJT3697
c ................................................................      ABJT3698
                                                                        ABJT3699
c ... type declaration statements                                       ABJT3700
      dimension photen(nphot),abscfs(nphot),sctcfs(nphot)               ABJT3701
      dimension emscfs(nphot)                                           ABJT3702
                                                                        ABJT3703
      logical lbug                                                      ABJT3704
                                                                        ABJT3705
c *******************************************************************   ABJT3706
c                                                                       ABJT3707
c                            begin execution                            ABJT3708
c                                                                       ABJT3709
c *******************************************************************   ABJT3710
                                                                        ABJT3711
c ... check debug option                                                ABJT3712
      call debug ( 'opacgp',                                            ABJT3713
     &                       lbug )                                     ABJT3714
                                                                        ABJT3715
      if ( lbug ) write (6,900) tp,rho,trad,xgmin,xgmax                 ABJT3716
                                                                        ABJT3717
      hvmin = xgmin * trad                                              ABJT3718
      hvmax = xgmax * trad                                              ABJT3719
                                                                        ABJT3720
c ... find the index where the photon energy equals the lower group     ABJT3721
c     boundary (note: mesh pts should be placed at each boundary)       ABJT3722
      do 20 i=1,nphot                                                   ABJT3723
         if ( photen(i) .ge. hvmin*0.99999 ) go to 25                   ABJT3724
   20 continue                                                          ABJT3725
   25 iphot0 = i                                                        ABJT3726
                                                                        ABJT3727
c ... integrate to get the group opacities                              ABJT3728
c     ------------------------------------                              ABJT3729
                                                                        ABJT3730
c ... initialize values for the first point                             ABJT3731
                                                                        ABJT3732
      xg2 = photen(iphot0) / trad                                       ABJT3733
      emx = exp( -xg2 )                                                 ABJT3734
                                                                        ABJT3735
      fpa2 = xg2**3 * emx * abscfs(iphot0) / ( 1.-emx )                 ABJT3736
      fpa20 = abscfs(iphot0)                                            ABJT3737
      fpe2 = xg2**3 * emx * emscfs(iphot0) / ( 1.-emx )                 ABJT3738
      fr2 = xg2**4 * emx / (1.-emx)**2 /                                ABJT3739
     &      ( abscfs(iphot0) + sctcfs(iphot0) )                         ABJT3740
      fr20 = 1. / ( abscfs(iphot0) + sctcfs(iphot0) )                   ABJT3741
                                                                        ABJT3742
c ... loop over photon energy mesh points                               ABJT3743
                                                                        ABJT3744
      sumpa  = 0.                                                       ABJT3745
      sumpa0 = 0.                                                       ABJT3746
      sumpe  = 0.                                                       ABJT3747
      sumr   = 0.                                                       ABJT3748
      sumr0  = 0.                                                       ABJT3749
                                                                        ABJT3750
      iphot = iphot0                                                    ABJT3751
  100 iphot = iphot + 1                                                 ABJT3752
                                                                        ABJT3753
         if ( photen(iphot) .gt. hvmax*1.00001 ) go to 200              ABJT3754
                                                                        ABJT3755
         xg1 = xg2                                                      ABJT3756
         xg2 = photen(iphot) / trad                                     ABJT3757
                                                                        ABJT3758
         fpa1  = fpa2                                                   ABJT3759
         fpa10 = fpa20                                                  ABJT3760
         fpe1  = fpe2                                                   ABJT3761
         fr1   = fr2                                                    ABJT3762
         fr10  = fr20                                                   ABJT3763
                                                                        ABJT3764
         if ( xg2 .lt. 1.e20 ) then                                     ABJT3765
                                                                        ABJT3766
            emx = exp( -xg2 )                                           ABJT3767
                                                                        ABJT3768
            fpa2  = xg2**3 * emx * abscfs(iphot) / ( 1.-emx )           ABJT3769
            fpa20 = abscfs(iphot)                                       ABJT3770
            fpe2  = xg2**3 * emx * emscfs(iphot) / ( 1.-emx )           ABJT3771
            fr2   = xg2**4 * emx / (1.-emx)**2 /                        ABJT3772
     &              ( abscfs(iphot) + sctcfs(iphot) )                   ABJT3773
            fr20  = 1. / ( abscfs(iphot) + sctcfs(iphot) )              ABJT3774
                                                                        ABJT3775
            dxg = xg2 - xg1                                             ABJT3776
                                                                        ABJT3777
c ...       integrate using a logarithmic interpolation scheme          ABJT3778
                                                                        ABJT3779
c ...       Planck mean absorption opacities                            ABJT3780
            if ( abs( fpa1-fpa2 ) .gt. 1.e-3*fpa2 ) then                ABJT3781
              if ( fpa2 .ne. 0. ) then                                  ABJT3782
                dsumpa = dxg * ( fpa2-fpa1 ) / log( fpa2/fpa1 )         ABJT3783
              else                                                      ABJT3784
                dsumpa = 0.                                             ABJT3785
              endif                                                     ABJT3786
            else                                                        ABJT3787
              dsumpa = dxg * fpa1                                       ABJT3788
            endif                                                       ABJT3789
                                                                        ABJT3790
c ...       mean (non-weighted) absorption opacities                    ABJT3791
            if ( abs( fpa10-fpa20 ) .gt. 1.e-3*fpa20 ) then             ABJT3792
              if ( fpa20 .ne. 0. ) then                                 ABJT3793
                dsmpa0 = dxg * ( fpa20-fpa10 ) / log( fpa20/fpa10 )     ABJT3794
              else                                                      ABJT3795
                dsmpa0 = 0.                                             ABJT3796
              endif                                                     ABJT3797
            else                                                        ABJT3798
              dsmpa0 = dxg * fpa10                                      ABJT3799
            endif                                                       ABJT3800
                                                                        ABJT3801
c ...       Planck mean emission opacities                              ABJT3802
            if ( abs( fpe1-fpe2 ) .gt. 1.e-3*fpe2 ) then                ABJT3803
              if ( fpe2 .ne. 0. ) then                                  ABJT3804
                dsumpe = dxg * ( fpe2-fpe1 ) / log( fpe2/fpe1 )         ABJT3805
              else                                                      ABJT3806
                dsumpe = 0.                                             ABJT3807
              endif                                                     ABJT3808
            else                                                        ABJT3809
              dsumpe = dxg * fpe1                                       ABJT3810
            endif                                                       ABJT3811
                                                                        ABJT3812
c ...       Rosseland mean opacities                                    ABJT3813
            if ( abs( fr1-fr2 ) .gt. 1.e-3*fr2 ) then                   ABJT3814
              if ( fr2 .ne. 0. ) then                                   ABJT3815
                dsumr = dxg * ( fr2-fr1 ) / log( fr2/fr1 )              ABJT3816
              else                                                      ABJT3817
                dsumr = 0.                                              ABJT3818
              endif                                                     ABJT3819
            else                                                        ABJT3820
              dsumr = dxg * fr1                                         ABJT3821
            endif                                                       ABJT3822
                                                                        ABJT3823
c ...       mean (non-weighted) "transport" opacities                   ABJT3824
            if ( abs( fr10-fr20 ) .gt. 1.e-3*fr20 ) then                ABJT3825
              if ( fr20 .ne. 0. ) then                                  ABJT3826
                dsumr0 = dxg * ( fr20-fr10 ) / log( fr20/fr10 )         ABJT3827
              else                                                      ABJT3828
                dsumr0 = 0.                                             ABJT3829
              endif                                                     ABJT3830
            else                                                        ABJT3831
              dsumr0 = dxg * fr10                                       ABJT3832
            endif                                                       ABJT3833
                                                                        ABJT3834
            sumpa  = sumpa + dsumpa                                     ABJT3835
            sumpa0 = sumpa0 + dsmpa0                                    ABJT3836
            sumpe  = sumpe + dsumpe                                     ABJT3837
            sumr   = sumr + dsumr                                       ABJT3838
            sumr0  = sumr0 + dsumr0                                     ABJT3839
                                                                        ABJT3840
         endif                                                          ABJT3841
                                                                        ABJT3842
         if ( iphot .lt. nphot ) go to 100                              ABJT3843
                                                                        ABJT3844
  200 continue                                                          ABJT3845
                                                                        ABJT3846
c ... normalize to get the group opacities                              ABJT3847
                                                                        ABJT3848
      g5 = gint(5,xgmin,xgmax)                                          ABJT3849
      g6 = gint(6,xgmin,xgmax)                                          ABJT3850
      if ( xgmax.ge.1./con(7) .and. xgmin.le.con(7) .and.               ABJT3851
     &     g5.gt.0. ) then                                              ABJT3852
c ...    weight absorption coef. by Planck function                     ABJT3853
         opacpa = sumpa / rho / g5                                      ABJT3854
         opacpe = sumpe / rho / g5                                      ABJT3855
      else                                                              ABJT3856
c ...    use straight average of absorption coef.                       ABJT3857
         opacpa = sumpa0 / rho / ( xgmax-xgmin )                        ABJT3858
         opacpe = 0.                                                    ABJT3859
      endif                                                             ABJT3860
                                                                        ABJT3861
      if ( xgmax.ge.1./con(7) .and. xgmin.le.con(7) .and.               ABJT3862
     &     sumr.gt.0. .and. g6.gt.0. ) then                             ABJT3863
c ...    weight absorption and scattering coefs. by Planck function     ABJT3864
c        derivative wrt temperature                                     ABJT3865
         opacr = g6 / rho / sumr                                        ABJT3866
      else if ( sumr0.gt.0. ) then                                      ABJT3867
c ...    use straight average of absorption plus scattering coefs.      ABJT3868
         opacr = ( xgmax-xgmin ) / rho / sumr0                          ABJT3869
      else                                                              ABJT3870
c ...    in the very high photen energy limit, this just becomes equal  ABJT3871
c        to Thomson scattering contribution (assuming the absorption    ABJT3872
c        is much less)                                                  ABJT3873
         opacr = ( sctcfs(nphot) + abscfs(nphot) ) / rho                ABJT3874
      endif                                                             ABJT3875
                                                                        ABJT3876
      if ( lbug ) write (6,904) sumpa,sumpa0*g5/(xgmax-xgmin),          ABJT3877
     &                          g6/sumr,(xgmax-xgmin)/sumr0             ABJT3878
      if ( lbug ) write (6,905) iphot0,iphot-1,sumpa,sumpe,sumr,        ABJT3879
     &                          opacpa,opacpe,opacr                     ABJT3880
                                                                        ABJT3881
      return                                                            ABJT3882
                                                                        ABJT3883
c ... format statements                                                 ABJT3884
                                                                        ABJT3885
  900 format (' debug output from -opacgp-:'/t4,'tp',t18,'rho',         ABJT3886
     2  t32,'trad',t46,'xgmin',t60,'xgmax-1'/t4,1p5e14.4)               ABJT3887
  904 format (t4,'sumpa',t18,'sumpa0*g5/dxg',t32,'g6/sumr',             ABJT3888
     2  t46,'dxg/sumr0'/t4,1p4e14.4)                                    ABJT3889
  905 format (t4,'iphot0',t11,'iphot',t18,'sumpa',t32,'sumpe',          ABJT3890
     2  t46,'sumr',t60,'opacpa',t74,'opacpe',t88,'opacr'                ABJT3891
     3  /t4,i5,t11,i5,t18,1p6e14.4)                                     ABJT3892
                                                                        ABJT3893
      end                                                               ABJT3894
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT3895
                                                                        ABJT3896
      subroutine opacys ( tp,densnn,densne,photen,nphot,abscfs,emscfs,  ABJT3897
     &                    sctcfs,engrup,ngrups,ntrad,trad,              ABJT3898
     &                                 opacpa,opacpe,opacr,             ABJT3899
     &                                 opatot,opetot,ortot,culrat,      ABJT3900
     &                                 op2tp,op2tr )                    ABJT3901
                                                                        ABJT3902
c ... this routine computes the Planck and Rosseland opacities for      ABJT3903
c     the photon energy group "engrup".  "opacpm" and "opacrm" are the  ABJT3904
c     Planck and Rosseland mean opacities (integrated over frequency).  ABJT3905
c                                                                       ABJT3906
c ... input variables:                                                  ABJT3907
c       tp     - plasma temperature (eV)                                ABJT3908
c       densnn - number density of all nuclei (cm**-3)                  ABJT3909
c       densne - electron density (cm**-3)                              ABJT3910
c       photen - photon energies (eV)                                   ABJT3911
c       nphot  - number of elements in "photen" array                   ABJT3912
c       abscfs - array of absorption coefficients (cm**-1)              ABJT3913
c       emscfs - array of emission coefficients (cm**-1)                ABJT3914
c       sctcfs - array of scattering coefficients (cm**-1)              ABJT3915
c       engrup - photon energy group boundaries (eV) (ngrups+1)         ABJT3916
c       ngrups - number of energy bins                                  ABJT3917
c       ntrad  - number of radiation temperatures                       ABJT3918
c       trad   - array of radiation temperatures (eV)                   ABJT3919
c                                                                       ABJT3920
c ... output variables:                                                 ABJT3921
c       opacpa - Planck group opacities for absorption (cm**2/g)        ABJT3922
c       opacpe - Planck group opacities for emission (cm**2/g)          ABJT3923
c       opacr  - Rosseland group opacity (cm**2/g)                      ABJT3924
c       opatot - Planck mean opacity for absorption (cm**2/g)           ABJT3925
c       opetot - Planck mean opacity for emission (cm**2/g)             ABJT3926
c       ortot  - Rosseland mean opacity (cm**2/g)                       ABJT3927
c       culrat - plasma cooling rate (erg cm**3/sec)                    ABJT3928
c       op2tp  - Planck 2-temperature opacities (cm**2/g)               ABJT3929
c       op2tr  - Rosseland 2-temperature opacities (cm**2/g)            ABJT3930
c                                                                       ABJT3931
                                                                        ABJT3932
c ... set the maximum number of:                                        ABJT3933
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT3934
c       temperatures (mxtemp), densities (mxdens),                      ABJT3935
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT3936
                                                                        ABJT3937
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT3938
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT3939
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT3940
                                                                        ABJT3941
                                                                        ABJT3942
      common / consts / pi, avgdro, sbcon, hplank                       ABJT3943
c ...............................................................       ABJT3944
                                                                        ABJT3945
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT3946
c ................................................................      ABJT3947
                                                                        ABJT3948
      common / gases / ngases,izgas(mxgass),atomwt(mxgass),             ABJT3949
     2                 avgatw,avgatn,                                   ABJT3950
     3                 fracsp(mxgass),fraciz(mxgass,mxatom),            ABJT3951
     4                 fraclv(mxgass,mxatom,20),pot(mxgass,mxatom),     ABJT3952
     5                 numocc(mxgass,mxatom)                            ABJT3953
c ................................................................      ABJT3954
                                                                        ABJT3955
c ... declaration statements:                                           ABJT3956
      logical lbug                                                      ABJT3957
                                                                        ABJT3958
      dimension photen(nphot),sctcfs(nphot),abscfs(nphot)               ABJT3959
      dimension engrup(ngrups+1),opacpa(ngrups),opacr(ngrups)           ABJT3960
      dimension emscfs(nphot),opacpe(ngrups)                            ABJT3961
      dimension trad(mxtemp),op2tp(mxtemp),op2tr(mxtemp)                ABJT3962
                                                                        ABJT3963
c *******************************************************************   ABJT3964
c                                                                       ABJT3965
c                            begin execution                            ABJT3966
c                                                                       ABJT3967
c *******************************************************************   ABJT3968
                                                                        ABJT3969
c ... check debug option                                                ABJT3970
      call debug ( 'opacys',                                            ABJT3971
     &                       lbug )                                     ABJT3972
                                                                        ABJT3973
      if ( lbug ) write (6,900) (engrup(i),i=1,ngrups+1)                ABJT3974
                                                                        ABJT3975
      rho = densnn * avgatw / avgdro                                    ABJT3976
                                                                        ABJT3977
c ... compute the Planck and Rosseland opacities for each photon        ABJT3978
c     energy group                                                      ABJT3979
                                                                        ABJT3980
      opatot = 0.                                                       ABJT3981
      opetot = 0.                                                       ABJT3982
      ortot = 0.                                                        ABJT3983
      do 100 i=1,ngrups                                                 ABJT3984
         xgmin = engrup(i) / tp                                         ABJT3985
         xgmax = engrup(i+1) / tp                                       ABJT3986
         trad0 = tp                                                     ABJT3987
         call opacgp ( tp,rho,trad0,photen,nphot,abscfs,emscfs,sctcfs,  ABJT3988
     &                 xgmin,xgmax,                                     ABJT3989
     &                                  opacpa(i),opacpe(i),opacr(i) )  ABJT3990
         if ( isw(19) .eq. 2 ) then                                     ABJT3991
c ...       use analytic solution to bound-bound opacities              ABJT3992
            call opacbb ( tp,densnn,densne,trad0,xgmin,xgmax,           ABJT3993
     &                                          opabsl,opemsl,opross )  ABJT3994
            opacpa(i) = opacpa(i) + opabsl                              ABJT3995
            opacpe(i) = opacpe(i) + opemsl                              ABJT3996
            opacr(i)  = opacr(i)  + opross                              ABJT3997
         endif                                                          ABJT3998
                                                                        ABJT3999
         opatot = opatot + gint(5,xgmin,xgmax) * opacpa(i)              ABJT4000
         opetot = opetot + gint(5,xgmin,xgmax) * opacpe(i)              ABJT4001
         ortot  = ortot  + gint(6,xgmin,xgmax) / opacr(i)               ABJT4002
                                                                        ABJT4003
  100 continue                                                          ABJT4004
                                                                        ABJT4005
c ... compute the total opacities based on the group opacities          ABJT4006
      opatot = opatot / 6.4939                                          ABJT4007
      opetot = opetot / 6.4939                                          ABJT4008
      ortot = 25.976 / ortot                                            ABJT4009
                                                                        ABJT4010
      if ( lbug ) write (6,906) opatot,opetot,ortot                     ABJT4011
                                                                        ABJT4012
c ... use this temporary integration scheme to check Planck opacity     ABJT4013
c     and plasma emission rate                                          ABJT4014
                                                                        ABJT4015
      sumpa = 0.                                                        ABJT4016
      sumpe = 0.                                                        ABJT4017
      sumr = 0.                                                         ABJT4018
      do 200 i=2,nphot                                                  ABJT4019
         im1 = i - 1                                                    ABJT4020
         exui = exp( -photen(i)/tp )                                    ABJT4021
         exuim1 = exp( -photen(im1)/tp )                                ABJT4022
         sumpa = sumpa + ( photen(i)-photen(im1) ) * 0.5 *              ABJT4023
     2        ( abscfs(i) * photen(i)**3 * exui / ( 1.-exui )           ABJT4024
     3        + abscfs(im1)*photen(im1)**3*exuim1 / ( 1.-exuim1 ) )     ABJT4025
         sumpe = sumpe + ( photen(i)-photen(im1) ) * 0.5 *              ABJT4026
     2        ( emscfs(i) * photen(i)**3 * exui / ( 1.-exui )           ABJT4027
     3        + emscfs(im1)*photen(im1)**3*exuim1 / ( 1.-exuim1 ) )     ABJT4028
         sumr = sumr + ( photen(i)-photen(im1) ) * 0.5 *                ABJT4029
     2        ( photen(i)**4 * exui /                                   ABJT4030
     3        ((abscfs(i) + sctcfs(i) ) * (1.-exui)**2)                 ABJT4031
     4        + photen(im1)**4*exuim1/                                  ABJT4032
     5        ((abscfs(im1)+sctcfs(im1))*(1.-exuim1)**2) )              ABJT4033
  200 continue                                                          ABJT4034
                                                                        ABJT4035
      opmpa = sumpa * 15. / ( pi**4 * tp**4 * rho )                     ABJT4036
      opmpe = sumpe * 15. / ( pi**4 * tp**4 * rho )                     ABJT4037
      opmr  = 4. * pi**4 * tp**4 / ( 15. * rho * sumr )                 ABJT4038
                                                                        ABJT4039
      if ( isw(19) .eq. 2 ) then                                        ABJT4040
c ...    use analytic solution to the bound-bound opacities             ABJT4041
         xgmin = 0.01                                                   ABJT4042
         xgmax = 100.                                                   ABJT4043
         trad0 = tp                                                     ABJT4044
         call opacbb ( tp,densnn,densne,trad0,xgmin,xgmax,              ABJT4045
     &                                          opabsl,opemsl,opross )  ABJT4046
         opmpa = opmpa + opabsl                                         ABJT4047
         opmpe = opmpe + opemsl                                         ABJT4048
         opmr  = opmr  + opross                                         ABJT4049
      endif                                                             ABJT4050
                                                                        ABJT4051
c ... compute the plasma emission (cooling) rate                        ABJT4052
      if ( opetot.gt.0. .and. densne.gt.0. ) then                       ABJT4053
         culrat = log10( opetot*4.*sbcon*rho ) + 4.*log10( 11606.*tp )  ABJT4054
     &          - log10( densne ) - log10( densnn )                     ABJT4055
      else                                                              ABJT4056
         culrat = 0.                                                    ABJT4057
      endif                                                             ABJT4058
                                                                        ABJT4059
      if ( lbug ) write (6,907) opmpa,opmpe,opmr,culrat,rho,tp          ABJT4060
                                                                        ABJT4061
                                                                        ABJT4062
c ... compute the 2-temperature opacities                               ABJT4063
                                                                        ABJT4064
      if ( ntrad .ne. 0 ) then                                          ABJT4065
                                                                        ABJT4066
         xgmin = 0.01                                                   ABJT4067
         xgmax = 100.                                                   ABJT4068
         do 300 itrad=1,ntrad                                           ABJT4069
                                                                        ABJT4070
           call opacgp ( tp,rho,trad(itrad),photen,nphot,abscfs,emscfs, ABJT4071
     &                   sctcfs,xgmin,xgmax,                            ABJT4072
     &                                op2tp(itrad),dummy,op2tr(itrad) ) ABJT4073
           if ( isw(19) .eq. 2 ) then                                   ABJT4074
             call opacbb ( tp,densnn,densne,trad(itrad),xgmin,xgmax,    ABJT4075
     &                                           opabsl,opemsl,opross ) ABJT4076
             op2tp(itrad) = op2tp(itrad) + opabsl                       ABJT4077
             op2tr(itrad) = op2tr(itrad) + opross                       ABJT4078
           endif                                                        ABJT4079
                                                                        ABJT4080
  300    continue                                                       ABJT4081
                                                                        ABJT4082
      endif                                                             ABJT4083
                                                                        ABJT4084
                                                                        ABJT4085
      return                                                            ABJT4086
                                                                        ABJT4087
c ... format statements                                                 ABJT4088
                                                                        ABJT4089
  900 format (t2,'debug output from -opacys-:  engrup='/                ABJT4090
     2  10(t4,1p5e14.4/))                                               ABJT4091
  906 format (t4,'opatot',t18,'opetot',t32,'ortot'/t4,1p3e14.4)         ABJT4092
  907 format (t4,'opmpa',t18,'opmpe',t32,'opmr',t46,'culrat',t60,       ABJT4093
     2  'rho',t74,'tp'/t4,1p6e14.4)                                     ABJT4094
                                                                        ABJT4095
      end                                                               ABJT4096
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT4097
c END OF THE COPIED FRAGMENT
                                                                        ABJT4098
      subroutine owt1 ( tp,densnn,densne,enrgy,heatcp,dzdt,dedden,      ABJT4099
     &                  oppat,oppet,oprt,ngrups,                        ABJT4100
     &                  engrup,oppma,oppme,oprm,culrat )                ABJT4101
c NOT NEEDED                                                            ABJT4102
c ... subroutine to write out various results                           ABJT4103
c                                                                       ABJT4104
c ... input variables                                                   ABJT4105
c                                                                       ABJT4106
c ... output variable                                                   ABJT4107
c        none                                                           ABJT4108
c                                                                       ABJT4109
                                                                        ABJT4110
c ... set the maximum number of:                                        ABJT4111
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT4112
c       temperatures (mxtemp), densities (mxdens),                      ABJT4113
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT4114
                                                                        ABJT4115
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT4116
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT4117
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT4118
                                                                        ABJT4119
                                                                        ABJT4120
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT4121
c ................................................................      ABJT4122
                                                                        ABJT4123
      common / consts / pi, avgdro, sbcon, hplank                       ABJT4124
c ...............................................................       ABJT4125
                                                                        ABJT4126
      common / gases / ngases,izgas(mxgass),atomwt(mxgass),             ABJT4127
     2                 avgatw,avgatn,                                   ABJT4128
     3                 fracsp(mxgass),fraciz(mxgass,mxatom),            ABJT4129
     4                 fraclv(mxgass,mxatom,20),pot(mxgass,mxatom),     ABJT4130
     5                 numocc(mxgass,mxatom)                            ABJT4131
c ................................................................      ABJT4132
                                                                        ABJT4133
      common / params / npqmax(mxgass), npmaxp, nfrqbb,                 ABJT4134
     2                  npring(100), neopen(5,100), bbnton(6,55),       ABJT4135
     3                  defpot(27,55), noccdf(27,55)                    ABJT4136
c ................................................................      ABJT4137
                                                                        ABJT4138
c ... declaration statements                                            ABJT4139
      logical lbug                                                      ABJT4140
                                                                        ABJT4141
      dimension oppat(mxgrps),oppet(mxgrps),oprt(mxgrps)                ABJT4142
      dimension engrup(mxgrps+1)                                        ABJT4143
                                                                        ABJT4144
c ********************************************************************  ABJT4145
c                                                                       ABJT4146
c                             begin execution                           ABJT4147
c                                                                       ABJT4148
c ********************************************************************  ABJT4149
                                                                        ABJT4150
c ... check debug option                                                ABJT4151
      call debug ( 'owt1',                                              ABJT4152
     &                        lbug )                                    ABJT4153
                                                                        ABJT4154
c ... find the average charge state                                     ABJT4155
      zbar = densne / densnn                                            ABJT4156
      rho = densnn * avgatw / avgdro                                    ABJT4157
      pres = 8.314e7 * rho * 11605. * tp * (1.+zbar) / avgatw           ABJT4158
      tdegk = tp * 11605.                                               ABJT4159
                                                                        ABJT4160
c ... write results to output file                                      ABJT4161
                                                                        ABJT4162
      write (6,900)                                                     ABJT4163
      write (6,901) tp,densnn,rho,densne,zbar,enrgy,pres                ABJT4164
      if ( isw(4).eq.0 ) write (6,921) heatcp,dzdt,dedden               ABJT4165
      write (6,902) npmaxp,(npqmax(lgas),lgas=1,ngases)                 ABJT4166
                                                                        ABJT4167
      if ( isw(2) .eq. 0 ) then                                         ABJT4168
         write (6,903) oppma,oppme,oprm,culrat                          ABJT4169
                                                                        ABJT4170
c ...   if the temperature is outside the photon energy group boundariesABJT4171
c        print out a warning                                            ABJT4172
         if ( tp.lt.engrup(1) .or. tp.gt.engrup(ngrups+1) )             ABJT4173
     &      write (6,909)                                               ABJT4174
                                                                        ABJT4175
         write (6,904)                                                  ABJT4176
         do 100 ig=1,ngrups                                             ABJT4177
            write (6,905) ig,engrup(ig),engrup(ig+1),oppat(ig),         ABJT4178
     &                    oppet(ig),oprt(ig)                            ABJT4179
  100    continue                                                       ABJT4180
                                                                        ABJT4181
      endif                                                             ABJT4182
                                                                        ABJT4183
c ... write results to a plot file                                      ABJT4184
                                                                        ABJT4185
      if ( iplot(2) .eq. 1 ) then                                       ABJT4186
         write (12,911) tp,culrat,oppma,oppme,oprm                      ABJT4187
      endif                                                             ABJT4188
                                                                        ABJT4189
      if ( iplot(4) .eq. 1 ) then                                       ABJT4190
         write (14,911) densnn,zbar,enrgy,pres                          ABJT4191
      endif                                                             ABJT4192
                                                                        ABJT4193
      if ( iplot(5) .eq. 1 ) then                                       ABJT4194
         write (15,911) tp,culrat,zbar,enrgy,pres                       ABJT4195
      endif                                                             ABJT4196
                                                                        ABJT4197
      if ( iplot(6) .eq. 1 ) then                                       ABJT4198
         dpdt = densnn * 11605.*1.381e-16 * (1.+zbar+tp*dzdt)           ABJT4199
         if ( heatcp .ne. 0. ) then                                     ABJT4200
            ratio = dpdt / heatcp                                       ABJT4201
         else                                                           ABJT4202
            ratio = 0.                                                  ABJT4203
         endif                                                          ABJT4204
         write (16,911) tp,heatcp,dzdt,dpdt,ratio                       ABJT4205
      endif                                                             ABJT4206
                                                                        ABJT4207
      if ( iplot(7) .eq. 1 ) then                                       ABJT4208
         write (17,911) densnn,oppme/oppma,oppma,oppme,oprm             ABJT4209
      endif                                                             ABJT4210
                                                                        ABJT4211
      if ( iplot(8) .eq. 1 ) then                                       ABJT4212
         do 220 ig=1,ngrups                                             ABJT4213
            g5 = gint( 5, engrup(ig)/tp, engrup(ig+1)/tp )              ABJT4214
            g5odhv = g5 / ( engrup(ig+1) - engrup(ig) )                 ABJT4215
c ...       divide by the photon energy to get the "photon flux"        ABJT4216
            hvavg = 0.5 * (engrup(ig)+engrup(ig+1))                     ABJT4217
            dum1 = g5odhv * tp**4 / hvavg * hvavg**isw(10)              ABJT4218
            write (18,911) engrup(ig),  oppet(ig)*dum1,g5odhv           ABJT4219
            write (18,911) engrup(ig+1),oppet(ig)*dum1,g5odhv           ABJT4220
  220    continue                                                       ABJT4221
      endif                                                             ABJT4222
                                                                        ABJT4223
      return                                                            ABJT4224
                                                                        ABJT4225
c ... format statements                                                 ABJT4226
                                                                        ABJT4227
  900 format (///t20,'*******************************************'/     ABJT4228
     2           t20,'Results for next temperature, density point'/     ABJT4229
     3           t20,'*******************************************'//)   ABJT4230
  901 format (t2,'Temperature             =',1pe11.3,' eV'/             ABJT4231
     2        t2,'Number density          =',1pe11.3,' cm**-3'/         ABJT4232
     3        t2,'Mass density            =',1pe11.3,' grams/cm**3'//   ABJT4233
     4        t2,'Electron density        =',1pe11.3,' cm**-3'/         ABJT4234
     5        t2,'Average charge state    =',1pe11.3/                   ABJT4235
     6        t2,'Specific energy         =',1pe11.3,' J/gram'/         ABJT4236
     7        t2,'Pressure                =',1pe11.3,' dyne/cm**2')     ABJT4237
  902 format (/t2,'Max. prin. quantum # used to compute ',              ABJT4238
     2            'populations    =',i5/                                ABJT4239
     3         t2,'Max. prin. quantum # used to compute ',              ABJT4240
     4            'absorp. coefs. =',20i5/)                             ABJT4241
  903 format (t2,'Planck mean opacity (abs.) =',1pe11.3,' cm**2/gram'/  ABJT4242
     2        t2,'Planck mean opacity (ems.) =',1pe11.3,' cm**2/gram'/  ABJT4243
     3        t2,'Rosseland mean opacity     =',1pe11.3,' cm**2/gram'/  ABJT4244
     4        t2,'Plasma cooling rate        =',1pe11.3,                ABJT4245
     5           ' erg*cm**3/sec   (Log10 value)'/)                     ABJT4246
  904 format (/t20,'Group opacities:'//t2,'group',t10,                  ABJT4247
     2  'lower',t24,'upper',t38,'Planck (abs)',t52,'Planck (ems)',      ABJT4248
     3  t66,'Rosseland'/                                                ABJT4249
     4  t2,'number',t10,'boundary',t24,'boundary',t38,                  ABJT4250
     5  'opacity',t52,'opacity',t66,                                    ABJT4251
     6  'opacity'/t10,'(eV)',t24,'(eV)',t38,'(cm**2/g)',                ABJT4252
     7  t52,'(cm**2/g)',t66,'(cm**2/g)'//)                              ABJT4253
  905 format (t2,i3,t6,1p5e14.3)                                        ABJT4254
  909 format (/t2,'*** warning ***  the temperature is not within the', ABJT4255
     2  ' photon energy boundaries --'/                                 ABJT4256
     3  t4,'this means the mean opacities are suspect'/)                ABJT4257
  911 format (t2,1p5e14.4)                                              ABJT4258
  921 format (t2,'Heat capacity           =',1pe11.3,' J/gram/eV'/      ABJT4259
     2        t2,'d(Charge st.)/d(Temp.)  =',1pe11.3,' ev**-1'/         ABJT4260
     3        t2,'d(Spec. En.)/d(Dens.)   =',1pe11.3,' J*cm**3/g' )     ABJT4261
                                                                        ABJT4262
      end                                                               ABJT4263
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT4264
                                                                        ABJT4265
      subroutine owtf ( densne,enrgy,heatcp,dzdt,dedden,                ABJT4266
     &                  opma,opme,orm,opgpa,opgpe,orgp,op2tp,op2tr )    ABJT4267
                                                                        ABJT4268
c ... prints out final results                                          ABJT4269
c                                                                       ABJT4270
c ... input variables                                                   ABJT4271
c       densne  -  electron densities (cm**-3)                          ABJT4272
c       enrgy   -  specific energies (J/g)                              ABJT4273
c       heatcp  -  heat capacities (J/g/eV)                             ABJT4274
c       dzdt    -  temperature derivatives of the average charge        ABJT4275
c                  states (eV**-1)                                      ABJT4276
c       dedden  -  density derivative of the specific energy (J*cm**3/g)ABJT4277
c       opm     -  Planck mean opacities (cm**2/g)                      ABJT4278
c       orm     -  Rosseland mean opacities (cm**2/g)                   ABJT4279
c       opgp    -  Planck group opacities (cm**2/g)                     ABJT4280
c       orgp    -  Rosseland group opacities (cm**2/g)                  ABJT4281
c       op2tp   -  Planck 2-temperature opacities (cm**2/g)             ABJT4282
c       op2tr   -  Rosseland 2-temperature opacities (cm**2/g)          ABJT4283
c                                                                       ABJT4284
c ... output variable                                                   ABJT4285
c        none                                                           ABJT4286
c                                                                       ABJT4287
                                                                        ABJT4288
c ... set the maximum number of:                                        ABJT4289
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT4290
c       temperatures (mxtemp), densities (mxdens),                      ABJT4291
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT4292
                                                                        ABJT4293
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT4294
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT4295
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT4296
                                                                        ABJT4297
                                                                        ABJT4298
      common / consts / pi, avgdro, sbcon, hplank                       ABJT4299
c ...............................................................       ABJT4300
                                                                        ABJT4301
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT4302
c ................................................................      ABJT4303
                                                                        ABJT4304
      common / gases / ngases,izgas(mxgass),atomwt(mxgass),             ABJT4305
     2                 avgatw,avgatn,                                   ABJT4306
     3                 fracsp(mxgass),fraciz(mxgass,mxatom),            ABJT4307
     4                 fraclv(mxgass,mxatom,20),pot(mxgass,mxatom),     ABJT4308
     5                 numocc(mxgass,mxatom)                            ABJT4309
c ................................................................      ABJT4310
                                                                        ABJT4311
      common / opacs / nptspg, ngrups, engrup(mxgrps+1)                 ABJT4312
c ................................................................      ABJT4313
                                                                        ABJT4314
      common / strngs / header, namedb(10)                              ABJT4315
                                                                        ABJT4316
      character header*80, namedb*6                                     ABJT4317
c ................................................................      ABJT4318
                                                                        ABJT4319
      common / tdmesh / ntemp,ndens,ntrad,                              ABJT4320
     2                  tplsma(mxtemp),densnn(mxdens),trad(mxtemp),     ABJT4321
     2                  dlgtmp,dlgden,dlgtrd                            ABJT4322
c ................................................................      ABJT4323
                                                                        ABJT4324
c ... declaration statements                                            ABJT4325
                                                                        ABJT4326
      dimension enrgy(mxtemp,mxdens),heatcp(mxtemp,mxdens)              ABJT4327
      dimension dzdt(mxtemp,mxdens),densne(mxtemp,mxdens)               ABJT4328
      dimension dedden(mxtemp,mxdens)                                   ABJT4329
      dimension opgpa(mxtemp,mxdens,mxgrps),opma(mxtemp,mxdens)         ABJT4330
      dimension opgpe(mxtemp,mxdens,mxgrps),opme(mxtemp,mxdens)         ABJT4331
      dimension orgp(mxtemp,mxdens,mxgrps),orm(mxtemp,mxdens)           ABJT4332
      dimension op2tp(mxtemp,mxdens,mxtemp)                             ABJT4333
      dimension op2tr(mxtemp,mxdens,mxtemp)                             ABJT4334
                                                                        ABJT4335
      logical lbug                                                      ABJT4336
                                                                        ABJT4337
      character*80 headr2,headr3                                        ABJT4338
                                                                        ABJT4339
c ********************************************************************  ABJT4340
c                                                                       ABJT4341
c                             begin execution                           ABJT4342
c                                                                       ABJT4343
c ********************************************************************  ABJT4344
                                                                        ABJT4345
c ... check debug option                                                ABJT4346
      call debug ( 'owtf',                                              ABJT4347
     &                        lbug )                                    ABJT4348
                                                                        ABJT4349
      if ( isw(8).eq.1 .or. isw(8).eq.12 .or. isw(8).eq.13 ) then       ABJT4350
                                                                        ABJT4351
c ...    write output in "original" (i.e., late 1986) CONRAD-acceptable ABJT4352
c        format                                                         ABJT4353
                                                                        ABJT4354
         write (8,980) header                                           ABJT4355
         write (8,981) dlgden,log10( densnn(1) ),dlgtmp,                ABJT4356
     &                 log10( tplsma(1) ),ngrups                        ABJT4357
         write (8,991) ((densne(it,id)/densnn(id),it=1,ntemp),          ABJT4358
     &                  id=1,ndens)                                     ABJT4359
         write (8,991) ((enrgy(it,id),it=1,ntemp),id=1,ndens)           ABJT4360
         write (8,991) (((op2tr(itp,id,itr),itr=1,ntrad),               ABJT4361
     &                  itp=1,ntemp),id=1,ndens)                        ABJT4362
         write (8,991) (((op2tp(itp,id,itr),itr=1,ntrad),               ABJT4363
     &                  itp=1,ntemp),id=1,ndens)                        ABJT4364
         write (8,991) (engrup(ig),ig=1,ngrups+1)                       ABJT4365
         write (8,991) (((orgp(it,id,ig),it=1,ntemp),id=1,ndens),       ABJT4366
     &                  ig=1,ngrups)                                    ABJT4367
         write (8,991) (((opgpe(it,id,ig),it=1,ntemp),id=1,ndens),      ABJT4368
     &                  ig=1,ngrups)                                    ABJT4369
                                                                        ABJT4370
      endif                                                             ABJT4371
                                                                        ABJT4372
                                                                        ABJT4373
      if ( isw(8).eq.2 .or. isw(8).eq.12 ) then                         ABJT4374
                                                                        ABJT4375
                                                                        ABJT4376
      endif                                                             ABJT4377
                                                                        ABJT4378
      if ( isw(8).eq.3 .or. isw(8).eq.13 ) then                         ABJT4379
                                                                        ABJT4380
c ...    write output in "new" (i.e., post-September 1987)              ABJT4381
c        CONRAD-acceptable format                                       ABJT4382
                                                                        ABJT4383
c ...    set up additional header records first                         ABJT4384
         write (headr2,921) (izgas(l),l=1,ngases)                       ABJT4385
         write (headr3,922) (fracsp(l),l=1,ngases)                      ABJT4386
                                                                        ABJT4387
c ...    now, write data tables                                         ABJT4388
                                                                        ABJT4389
         write (10,980) header                                          ABJT4390
         write (10,980) headr2                                          ABJT4391
         write (10,980) headr3                                          ABJT4392
         write (10,981) dlgden,log10( densnn(1) ),dlgtmp,               ABJT4393
     &                  log10( tplsma(1) ),ngrups                       ABJT4394
         write (10,991) ((densne(it,id)/densnn(id),it=1,ntemp),         ABJT4395
     &                  id=1,ndens)                                     ABJT4396
         write (10,991) ((enrgy(it,id),it=1,ntemp),id=1,ndens)          ABJT4397
         write (10,991) ((heatcp(it,id),it=1,ntemp),id=1,ndens)         ABJT4398
         condd = -6.242e18 * avgatw / avgdro                            ABJT4399
         write (10,991) ((dedden(it,id)*condd*densnn(id)/tplsma(it)**2, ABJT4400
     &                  it=1,ntemp),id=1,ndens)                         ABJT4401
         write (10,991) (engrup(ig),ig=1,ngrups+1)                      ABJT4402
         write (10,991) (((orgp(it,id,ig),it=1,ntemp),id=1,ndens),      ABJT4403
     &                  ig=1,ngrups)                                    ABJT4404
         write (10,991) (((opgpa(it,id,ig),it=1,ntemp),id=1,ndens),     ABJT4405
     &                  ig=1,ngrups)                                    ABJT4406
         write (10,991) (((opgpe(it,id,ig),it=1,ntemp),id=1,ndens),     ABJT4407
     &                  ig=1,ngrups)                                    ABJT4408
                                                                        ABJT4409
      endif                                                             ABJT4410
                                                                        ABJT4411
      return                                                            ABJT4412
                                                                        ABJT4413
c ... format statements                                                 ABJT4414
                                                                        ABJT4415
  921 format (' atomic #s of gases: ',5i10)                             ABJT4416
  922 format (' relative fractions: ',1p5e10.2)                         ABJT4417
  980 format (a80)                                                      ABJT4418
  981 format (4e12.6,i12)                                               ABJT4419
  991 format (4e12.6)                                                   ABJT4420
c END OF NON-NEEDED FRAGMENT                                            ABJT4421
      end                                                               ABJT4422
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT4423
                                                                        ABJT4424
      function radrec ( q,n,nprin0,tempel,potiz )                       ABJT4425
c NON-LTE, NOT NEEDED                                                   ABJT4426
c ... radiative recombination rate function.  calculates the            ABJT4427
c     rate (per ion per free electron => cm**3/sec) at which            ABJT4428
c     electrons enter the "n"th level of an ion with charge state       ABJT4429
c     "q".  "q+1" is the original ionization state of the atom.         ABJT4430
c                                                                       ABJT4431
c ... input variables                                                   ABJT4432
c       q       =  charge state of the ion after recombination          ABJT4433
c                  (0 < q < atomic # - 1)                               ABJT4434
c       n       =  principal quantum # of the electron after being      ABJT4435
c                  recombined                                           ABJT4436
c       nprin0  =  principal quantum # of the valence electrons in      ABJT4437
c                  the ion's ground state                               ABJT4438
c       tempel  =  electron temperature (eV)                            ABJT4439
c       potiz   =  ionization potential for unexcited ion with charge   ABJT4440
c                  state q (eV)                                         ABJT4441
c                                                                       ABJT4442
c ... output variable                                                   ABJT4443
c       radrec  =  radiative recombination rate (cm**3/sec)             ABJT4444
c                                                                       ABJT4445
c ... declaration statements                                            ABJT4446
      logical lbug                                                      ABJT4447
                                                                        ABJT4448
c ********************************************************************  ABJT4449
c                                                                       ABJT4450
c                             begin execution                           ABJT4451
c                                                                       ABJT4452
c ********************************************************************  ABJT4453
                                                                        ABJT4454
c ... check debug option                                                ABJT4455
      call debug ( 'radrec',                                            ABJT4456
     &                        lbug )                                    ABJT4457
                                                                        ABJT4458
      tkev = tempel / 1000.                                             ABJT4459
                                                                        ABJT4460
c ... ionization potential based on Bohr model                          ABJT4461
      potniz = nprin0**2 * potiz / n**2                                 ABJT4462
                                                                        ABJT4463
      xn = potniz / tempel                                              ABJT4464
                                                                        ABJT4465
c ... compute the exponential integral (see Abramowitz & Stegun), and   ABJT4466
c     multiply by exp(xn)                                               ABJT4467
      if ( xn .ge. 1. ) then                                            ABJT4468
         expint = (xn**2+2.3347*xn+0.25062) / (xn**2+3.3306*xn+1.6815)  ABJT4469
     &            / xn                                                  ABJT4470
      else                                                              ABJT4471
         expint = -log( xn ) - 0.57722 + 0.99999*xn - 0.24991*xn**2     ABJT4472
     &            + 0.05520*xn**3 - 0.00976*xn**4 + 0.00108*xn**5       ABJT4473
         expint = expint * exp( xn )                                    ABJT4474
      endif                                                             ABJT4475
                                                                        ABJT4476
c ... compute the radiative recombination rate                          ABJT4477
      radrec = 5.20e-14 * (q+1.) * xn**1.5 * expint                     ABJT4478
                                                                        ABJT4479
      if ( lbug ) write (6,901) radrec,xn,q,expint,tempel               ABJT4480
                                                                        ABJT4481
      return                                                            ABJT4482
                                                                        ABJT4483
c ... format statements                                                 ABJT4484
                                                                        ABJT4485
  901 format (t2,'debug output from -radrec-'/t4,'radrec',t18,          ABJT4486
     2  'xn',t32,'q',t46,'expint',t60,'tempel'/t4,1p5e14.4)             ABJT4487
                                                                        ABJT4488
      end                                                               ABJT4489
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT4490
                                                                        ABJT4491
      subroutine rizrec ( tempel,densnn,q,atomnm,edens0,pfrat,          ABJT4492
     &                    potniz,nvalel,                                ABJT4493
     &                                                         ratio )  ABJT4494
c                                                                       ABJT4495
c ... computes the ratio of the total ionization rate to the            ABJT4496
c     to the total recombination rate for use in the Coronal            ABJT4497
c     Equilibrium model.  The ionization rate is computed for the       ABJT4498
c     transition from "q" -> "q+1", and the recombination rate          ABJT4499
c     for "q+1" -> "q".                                                 ABJT4500
c                                                                       ABJT4501
c ... input variables                                                   ABJT4502
c       tempel  =  electron temperature (eV)                            ABJT4503
c       densnn  =  atomic number density (cm**-3)                       ABJT4504
c       q       =  charge state of the ion before ionization (or        ABJT4505
c                  after recombination) (0 < q < atomnm-1)              ABJT4506
c       atomnm  =  atomic number of the ion                             ABJT4507
c       edens0  =  initial guess at electron density (cm**-3)           ABJT4508
c       pfrat   =  ratio of electronic partition functions (q/q+1)      ABJT4509
c       potniz  =  ionization potential for an unexcited ion with       ABJT4510
c                  charge state q (eV)                                  ABJT4511
c       nvalel  =  number of valence electrons for the ion in its       ABJT4512
c                  ground state                                         ABJT4513
c                                                                       ABJT4514
c ... output variable                                                   ABJT4515
c       ratio   =  ratio of the ionization rate to recombination rate   ABJT4516
c                                                                       ABJT4517
                                                                        ABJT4518
c ... set the maximum number of:                                        ABJT4519
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT4520
c       temperatures (mxtemp), densities (mxdens),                      ABJT4521
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT4522
                                                                        ABJT4523
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT4524
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT4525
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT4526
                                                                        ABJT4527
                                                                        ABJT4528
      common / contrl / isw(30),con(30),iplot(30),dtheat,critsc         ABJT4529
c ................................................................      ABJT4530
                                                                        ABJT4531
      common / params / npqmax(mxgass), npmaxp, nfrqbb,                 ABJT4532
     2                  npring(100), neopen(5,100), bbnton(6,55),       ABJT4533
     3                  defpot(27,55), noccdf(27,55)                    ABJT4534
c ................................................................      ABJT4535
                                                                        ABJT4536
c ... declaration statements                                            ABJT4537
      logical lbug                                                      ABJT4538
                                                                        ABJT4539
c ********************************************************************  ABJT4540
c                                                                       ABJT4541
c                             begin execution                           ABJT4542
c                                                                       ABJT4543
c ********************************************************************  ABJT4544
                                                                        ABJT4545
c ... check debug option                                                ABJT4546
      call debug ( 'rizrec',                                            ABJT4547
     &                        lbug )                                    ABJT4548
                                                                        ABJT4549
      if ( lbug ) write (6,900) q,tempel,potniz,atomnm,nvalel           ABJT4550
                                                                        ABJT4551
c ... Mosher formula (from NRL report, and originally used by MIXERG)   ABJT4552
c      potnot = potniz / tempel                                         ABJT4553
c      ratio1 = 1.27e8 * exp( -potnot ) * nvalel /                      ABJT4554
c     &        ( potnot**0.75 * potniz**2 )                             ABJT4555
                                                                        ABJT4556
c ... more exact calculation                                            ABJT4557
                                                                        ABJT4558
c ... find the principal quantum number of the valence electrons        ABJT4559
c     in the ground state (before recombination)                        ABJT4560
      nbound = atomnm - q - 1.                                          ABJT4561
      if ( nbound .gt. 0 ) then                                         ABJT4562
         nprin0 = npring(nbound)                                        ABJT4563
      else                                                              ABJT4564
         nprin0 = 1                                                     ABJT4565
      endif                                                             ABJT4566
      nprin = nprin0                                                    ABJT4567
      izgas = atomnm                                                    ABJT4568
                                                                        ABJT4569
      coll = colliz ( q,nprin,nprin0,tempel,atomnm,potniz )             ABJT4570
      radr = radrec ( q,nprin,nprin0,tempel,potniz )                    ABJT4571
      if ( nbound.ge.1 .and. isw(16).eq.0 ) then                        ABJT4572
         diel = dielrc ( q,nprin,nprin0,nbound,izgas,tempel,edens0,     ABJT4573
     &                   potniz )                                       ABJT4574
      else                                                              ABJT4575
         diel = 0.                                                      ABJT4576
      endif                                                             ABJT4577
                                                                        ABJT4578
c ... 3-body recombination calculation                                  ABJT4579
      if ( isw(6) .eq. 3 ) then                                         ABJT4580
         r3body = collrc ( q,nprin,nprin0,tempel,atomnm,potniz,         ABJT4581
     &                     edens0,pfrat )                               ABJT4582
      else                                                              ABJT4583
         r3body = 0.                                                    ABJT4584
      endif                                                             ABJT4585
                                                                        ABJT4586
      ratio = coll / ( radr + diel + r3body )                           ABJT4587
                                                                        ABJT4588
      if ( lbug ) write (6,906) coll,radr,diel,ratio,r3body             ABJT4589
                                                                        ABJT4590
      return                                                            ABJT4591
                                                                        ABJT4592
c ... format statements                                                 ABJT4593
                                                                        ABJT4594
  900 format (t2,'debug output from -rizrec-'/t4,'q',t18,'tempel',      ABJT4595
     2  t32,'potniz',t46,'atomnm',t60,'nvalel'/t4,1p4e12.3,             ABJT4596
     3  t60,0p,i5)                                                      ABJT4597
  906 format (t4,'coll',t18,'radr',t32,                                 ABJT4598
     2  'diel',t46,'ratio',t60,'r3body'/t4,1p7e12.3)                    ABJT4599
                                                                        ABJT4600
      end                                                               ABJT4601
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT4602
                                                                        ABJT4603
      subroutine saha ( tp,densnn,fracsp,ngases,izgas,pot,              ABJT4604
     &                                                 densne,fraciz )  ABJT4605
c                                                                       ABJT4606
c ... this routine computes the populations of charge states using      ABJT4607
c     the Saha model.                                                   ABJT4608
c                                                                       ABJT4609
c ... input variables                                                   ABJT4610
c       tp        -  plasma temperature (eV)                            ABJT4611
c       densnn    -  ion density (cm**-3)                               ABJT4612
c       fracsp(l) -  fraction of total gas atoms that is "l"th species  ABJT4613
c       ngases    -  number of gas species                              ABJT4614
c       izgas     -  atomic number of each gas                          ABJT4615
c       pot(l,j)  -  ionization potential for the ground state of "j-1" ABJT4616
c                    charge state of the "l"th gas (eV)                 ABJT4617
c                                                                       ABJT4618
c ... output variables                                                  ABJT4619
c       densne      -  electron density (cm**-3)                        ABJT4620
c       fraciz(l,j) -  fraction of ions of gas "l" in the "j-1" charge  ABJT4621
c                      state; i.e., fraciz(l,1) => neutral atom.        ABJT4622
c                                                                       ABJT4623
                                                                        ABJT4624
c ... set the maximum number of:                                        ABJT4625
c       gases (mxgass), photen energy mesh points (mxphot),             ABJT4626
c       temperatures (mxtemp), densities (mxdens),                      ABJT4627
c       opacity groups (mxgrps), and gas atomic number (= mxatom-1)     ABJT4628
                                                                        ABJT4629
      parameter ( mxtemp = 20, mxdens = 10 )                            ABJT4630
      parameter ( mxphot = 3000, mxgrps = 100 )                         ABJT4631
      parameter ( mxgass = 10, mxatom = 55 )                            ABJT4632
                                                                        ABJT4633
                                                                        ABJT4634
      common / params / npqmax(mxgass), npmaxp, nfrqbb,                 ABJT4635
     2                  npring(100), neopen(5,100), bbnton(6,55),       ABJT4636
     3                  defpot(27,55), noccdf(27,55)                    ABJT4637
c ................................................................      ABJT4638
                                                                        ABJT4639
c ... type declaration statements                                       ABJT4640
                                                                        ABJT4641
      dimension izgas(mxgass),s(mxgass)                                 ABJT4642
      dimension pot(mxgass,mxatom),fraciz(mxgass,mxatom)                ABJT4643
      dimension fracsp(mxgass),phi(mxgass,mxatom),pfsum(mxatom)         ABJT4644
      dimension sumphi(mxgass,mxatom),check(mxatom),xlnp(mxgass,mxatom) ABJT4645
      dimension dpdne(mxgass,mxatom),dsdne(mxgass)                      ABJT4646
                                                                        ABJT4647
      logical lbug, lincr                                               ABJT4648
                                                                        ABJT4649
c ... csaha = .5 * ( 2 pi me k / h**2 / 11605 )**-1.5                   ABJT4650
      data csaha / 1.66e-22 /                                           ABJT4651
                                                                        ABJT4652
c ********************************************************************  ABJT4653
c                                                                       ABJT4654
c                         begin execution                               ABJT4655
c                                                                       ABJT4656
c ********************************************************************  ABJT4657
                                                                        ABJT4658
c ... check debug option                                                ABJT4659
      call debug ( 'saha  ',                                            ABJT4660
     &                      lbug )                                      ABJT4661
                                                                        ABJT4662
      if ( lbug ) write (6,901) tp,densnn,ngases,npmaxp                 ABJT4663
                                                                        ABJT4664
c ... initialize variables                                              ABJT4665
      csahap = csaha * densnn / tp**1.5                                 ABJT4666
                                                                        ABJT4667
c ... set up arrays used in the Saha recursive relations                ABJT4668
c                                                                       ABJT4669
c ... first, calculate the electronic partition functions               ABJT4670
      sphimx = 0.                                                       ABJT4671
      do 250 lgas=1,ngases                                              ABJT4672
                                                                        ABJT4673
         izgm1 = izgas(lgas)-1                                          ABJT4674
         do 100 izlevl=0,izgm1                                          ABJT4675
                                                                        ABJT4676
c ...       find the principal quantum number of the valence electrons  ABJT4677
c           in their ground state                                       ABJT4678
            nbound = izgas(lgas) - izlevl                               ABJT4679
            nprin0 = npring(nbound)                                     ABJT4680
                                                                        ABJT4681
            izp1 = izlevl + 1                                           ABJT4682
            eizokt = pot(lgas,izp1) / tp                                ABJT4683
                                                                        ABJT4684
c ...       compute the electronic partition function                   ABJT4685
                                                                        ABJT4686
            qsum = 0.                                                   ABJT4687
            do 50 n=nprin0,npmaxp                                       ABJT4688
c ...          calculate the energy relative to the ground state        ABJT4689
               ennpot = nprin0**2 * eizokt * (1./nprin0**2-1./n**2)     ABJT4690
               qsum = qsum + 2.*n**2 * exp( -ennpot )                   ABJT4691
   50       continue                                                    ABJT4692
            pfsum(izp1) = qsum                                          ABJT4693
                                                                        ABJT4694
  100    continue                                                       ABJT4695
                                                                        ABJT4696
c ...    the last partition fn corresponds to the fully ionized nucleus ABJT4697
         pfsum(izgas(lgas)+1) = 1.                                      ABJT4698
                                                                        ABJT4699
c ...    next, calculate "sumphi", the constants (i.e., independent of  ABJT4700
c        the electron density) in the recursion relations               ABJT4701
c        Note: the phi(l,i) have been multiplied by the total number    ABJT4702
c        density to "normalize"                                         ABJT4703
                                                                        ABJT4704
         do 200 izlevl=0,izgm1                                          ABJT4705
            izp1 = izlevl + 1                                           ABJT4706
            phi(lgas,izp1) = log( csahap*pfsum(izp1)/pfsum(izp1+1) )    ABJT4707
     &                     + pot(lgas,izp1) / tp                        ABJT4708
  200    continue                                                       ABJT4709
                                                                        ABJT4710
c ...    "sumphi" is the sum of (phi)*(total # density)                 ABJT4711
                                                                        ABJT4712
         do 235 izlevl=0,izgm1                                          ABJT4713
            izp1 = izlevl + 1                                           ABJT4714
            sumphi(lgas,izp1) = 0.                                      ABJT4715
            do 230 k=izlevl,izgm1                                       ABJT4716
               sumphi(lgas,izp1) = sumphi(lgas,izp1) + phi(lgas,k+1)    ABJT4717
               sphimx = max ( sphimx , sumphi(lgas,izp1) )              ABJT4718
  230       continue                                                    ABJT4719
  235    continue                                                       ABJT4720
                                                                        ABJT4721
         sumphi(lgas,izgas(lgas)+1) = 0.                                ABJT4722
                                                                        ABJT4723
  250 continue                                                          ABJT4724
                                                                        ABJT4725
c ... normalize "sumphi"                                                ABJT4726
                                                                        ABJT4727
      do 270 lgas=1,ngases                                              ABJT4728
      do 270 izlevl=0,izgas(lgas)                                       ABJT4729
         izp1 = izlevl + 1                                              ABJT4730
         sumphi(lgas,izp1) = sumphi(lgas,izp1) - sphimx                 ABJT4731
  270 continue                                                          ABJT4732
                                                                        ABJT4733
      if ( lbug ) then                                                  ABJT4734
         write (6,902) (pfsum(i),i=1,izgas(1)+1)                        ABJT4735
         write (6,903) (phi(1,i),i=1,izgas(1)+1)                        ABJT4736
         write (6,904) (sumphi(1,i),i=1,izgas(1)+1)                     ABJT4737
      endif                                                             ABJT4738
                                                                        ABJT4739
                                                                        ABJT4740
c ... iterate to find the electron density and charge state populations ABJT4741
c     Note: the electron density is "normalized"; i.e., it is in units  ABJT4742
c     of the total number density                                       ABJT4743
c                                                                       ABJT4744
c ... first guess at the electron density (based on hydrogenic Saha eq.)ABJT4745
                                                                        ABJT4746
      avgne = 0.                                                        ABJT4747
      avgpot = 0.                                                       ABJT4748
      do 298 lgas=1,ngases                                              ABJT4749
         avgne = avgne + izgas(lgas)*fracsp(lgas)                       ABJT4750
         avgpot = avgpot + pot(lgas,1)*fracsp(lgas)                     ABJT4751
  298 continue                                                          ABJT4752
      beta = 3.e21 / densnn * tp**1.5 * exp( -avgpot/tp )               ABJT4753
      if ( beta .lt. 100. ) then                                        ABJT4754
         eldold = avgne * beta * 0.5 * ( sqrt(1.+4./beta) -1. )         ABJT4755
      else                                                              ABJT4756
         eldold = avgne                                                 ABJT4757
      endif                                                             ABJT4758
                                                                        ABJT4759
      lincr = .false.                                                   ABJT4760
      iter = 0                                                          ABJT4761
  300 iter = iter + 1                                                   ABJT4762
                                                                        ABJT4763
         sigma = 0.                                                     ABJT4764
         dsgdne = 0.                                                    ABJT4765
         eldlog = log( eldold )                                         ABJT4766
                                                                        ABJT4767
         do 400 lgas=1,ngases                                           ABJT4768
                                                                        ABJT4769
            xlnpmx = -1.e20                                             ABJT4770
                                                                        ABJT4771
c ...       find the relative populations of the ionization states      ABJT4772
                                                                        ABJT4773
            do 310 izlevl=0,izgas(lgas)                                 ABJT4774
               izp1 = izlevl + 1                                        ABJT4775
               xlnp(lgas,izp1) = (izgas(lgas)-izlevl) * eldlog          ABJT4776
     &                         + sumphi(lgas,izp1)                      ABJT4777
               xlnpmx = max( xlnpmx , xlnp(lgas,izp1) )                 ABJT4778
  310       continue                                                    ABJT4779
                                                                        ABJT4780
c ...       normalize relative populations to avoid overflows           ABJT4781
                                                                        ABJT4782
            do 315 izlevl=0,izgas(lgas)                                 ABJT4783
               izp1 = izlevl + 1                                        ABJT4784
               xlnp(lgas,izp1) = xlnp(lgas,izp1) - xlnpmx               ABJT4785
  315       continue                                                    ABJT4786
                                                                        ABJT4787
c ...       sum for the electron density                                ABJT4788
                                                                        ABJT4789
            s(lgas) = 0.                                                ABJT4790
            dsdne(lgas) = 0.                                            ABJT4791
                                                                        ABJT4792
            do 320 izlevl=0,izgas(lgas)                                 ABJT4793
               izp1 = izlevl + 1                                        ABJT4794
               dpdne(lgas,izp1) = (izgas(lgas)-izlevl) *                ABJT4795
     &                            exp( xlnp(lgas,izp1) ) / eldold       ABJT4796
               s(lgas) = s(lgas) + exp( xlnp(lgas,izp1) )               ABJT4797
               dsdne(lgas) = dsdne(lgas) + dpdne(lgas,izp1)             ABJT4798
  320       continue                                                    ABJT4799
                                                                        ABJT4800
            sum1 = 0.                                                   ABJT4801
            sum2 = 0.                                                   ABJT4802
            do 330 izlevl=0,izgas(lgas)                                 ABJT4803
               izp1 = izlevl + 1                                        ABJT4804
               fraciz(lgas,izp1) = exp( xlnp(lgas,izp1) ) / s(lgas)     ABJT4805
               sum1 = sum1 + izlevl * fraciz(lgas,izp1)                 ABJT4806
               sum2 = sum2 + izlevl * dpdne(lgas,izp1)                  ABJT4807
  330       continue                                                    ABJT4808
                                                                        ABJT4809
c ...       now, sum the electron density contribution from this speciesABJT4810
                                                                        ABJT4811
            sigma = sigma +  fracsp(lgas) * sum1                        ABJT4812
            dsgdne = dsgdne + fracsp(lgas) / s(lgas) *                  ABJT4813
     &               ( sum2 - sum1*dsdne(lgas) )                        ABJT4814
                                                                        ABJT4815
c ...       write debug output                                          ABJT4816
                                                                        ABJT4817
            if ( lbug ) then                                            ABJT4818
               write (6,906) iter,lgas,sigma,dsgdne,eldold              ABJT4819
            endif                                                       ABJT4820
                                                                        ABJT4821
  400    continue                                                       ABJT4822
                                                                        ABJT4823
c ...    calculate the adjustment to the new electron density           ABJT4824
         dne = ( sigma-eldold ) / ( 1.-dsgdne )                         ABJT4825
                                                                        ABJT4826
c ...    if the new value for the electron density does not "equal" the ABJT4827
c        old estimate, we need another iteration                        ABJT4828
                                                                        ABJT4829
         if ( abs( dne ) .lt. 1.e-4*eldold ) then                       ABJT4830
            densne = eldold + dne                                       ABJT4831
         else if ( iter .le. 20 ) then                                  ABJT4832
            if ( (dne.gt.0.) .neqv. lincr ) then                        ABJT4833
               eldnew = eldold + 0.5 * dne                              ABJT4834
            else                                                        ABJT4835
               eldnew = eldold + dne                                    ABJT4836
            endif                                                       ABJT4837
            if ( eldnew .le. 0. ) eldnew = eldold / 5.                  ABJT4838
            eldold = eldnew                                             ABJT4839
            lincr = dne .gt. 0.                                         ABJT4840
            if ( lbug ) write (6,907) dne,eldold                        ABJT4841
            go to 300                                                   ABJT4842
         else                                                           ABJT4843
            write (6,991) iter,dne,eldold                               ABJT4844
            stop ' loop did not converge in -saha-'                     ABJT4845
         endif                                                          ABJT4846
                                                                        ABJT4847
c ... end of iteration                                                  ABJT4848
                                                                        ABJT4849
c ... multiply by the total density of nuclei to get the electron       ABJT4850
c     density                                                           ABJT4851
      densne = densne * densnn                                          ABJT4852
                                                                        ABJT4853
      if ( lbug ) then                                                  ABJT4854
         write (6,910)  densne                                          ABJT4855
         do 510 lgas=1,ngases                                           ABJT4856
            write (6,911)  ( fraciz(lgas,izp1),izp1=1,izgas(lgas)+1 )   ABJT4857
            do 505 izlevl=0,izgas(lgas)-1                               ABJT4858
               izp1 = izlevl + 1                                        ABJT4859
               check(izp1+1) = fraciz(lgas,izp1)*(densnn/densne)        ABJT4860
     &                       * exp( -phi(lgas,izp1) )                   ABJT4861
  505       continue                                                    ABJT4862
            write (6,912) 0.,( check(izp1),izp1=2,izgas(lgas)+1)        ABJT4863
  510    continue                                                       ABJT4864
      endif                                                             ABJT4865
                                                                        ABJT4866
      return                                                            ABJT4867
                                                                        ABJT4868
c ... format statements                                                 ABJT4869
                                                                        ABJT4870
  901 format (' debug output from -saha-:'/t4,'tp',t18,'densnn',t32,    ABJT4871
     2  'ngases',t46,'npmaxp'/t4,1p2e14.4,0p,2(i4,10x))                 ABJT4872
  902 format ('  pfsum:'/6(t4,1p10e11.3/))                              ABJT4873
  903 format ('  phi:'/6(t4,1p10e11.3/))                                ABJT4874
  904 format ('  sumphi:'/6(t4,1p10e11.3/))                             ABJT4875
  906 format (t4,'iter',t18,'lgas',t32,'sigma',t46,'dsgdne',t60,        ABJT4876
     2  'eldold'/t4,i4,t18,i4,t32,1p3e11.3)                             ABJT4877
  907 format ( t4,'dne',t18,'eldold'/t4,1p2e11.3)                       ABJT4878
  910 format (t4,'densne'/t4,1pe14.4)                                   ABJT4879
  911 format (t2,'fraciz:'/8(t4,1p7e11.3/))                             ABJT4880
  912 format (t2,'check :'/8(t4,1p7e11.3/))                             ABJT4881
  991 format (///' loop did not converge in -saha-:'//t4,'iter',        ABJT4882
     2 t18,'eldnew',t32,'eldold'/t4,i4,t18,1p2e14.3)                    ABJT4883
c END OF NON-NEEDED FRAGMENT                                            ABJT4884
      end                                                               ABJT4885
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT4886
                                                                        ABJT4887
      subroutine sort ( npts,                                           ABJT4888
     &                            xin )                                 ABJT4889
c COPIED TO ModMultiGroup                                               ABJT4890
c ... puts the elements of "xin" into ascending order                   ABJT4891
                                                                        ABJT4892
      dimension xin(npts)                                               ABJT4893
      logical flag                                                      ABJT4894
                                                                        ABJT4895
c ********************************************************************  ABJT4896
c                                                                       ABJT4897
c                            begin execution                            ABJT4898
c                                                                       ABJT4899
c ********************************************************************  ABJT4900
                                                                        ABJT4901
      ncycle = 0                                                        ABJT4902
                                                                        ABJT4903
c ... first, throw out those values < or = zero                         ABJT4904
                                                                        ABJT4905
      j = 0                                                             ABJT4906
      do 5 i=1,npts                                                     ABJT4907
         if ( xin(i) .gt. 0. ) then                                     ABJT4908
            j = j + 1                                                   ABJT4909
            xin(j) = xin(i)                                             ABJT4910
         endif                                                          ABJT4911
    5 continue                                                          ABJT4912
      npts = j                                                          ABJT4913
                                                                        ABJT4914
   15 flag = .false.                                                    ABJT4915
                                                                        ABJT4916
         do 10 i=2,npts                                                 ABJT4917
            im1 = i - 1                                                 ABJT4918
            if ( xin(im1) .gt. xin(i) ) then                            ABJT4919
               save = xin(im1)                                          ABJT4920
               xin(im1) = xin(i)                                        ABJT4921
               xin(i) = save                                            ABJT4922
               flag = .true.                                            ABJT4923
             endif                                                      ABJT4924
   10    continue                                                       ABJT4925
         ncycle = ncycle + 1                                            ABJT4926
                                                                        ABJT4927
c ...    check for error                                                ABJT4928
                                                                        ABJT4929
         if ( ncycle .gt. npts ) then                                   ABJT4930
            write (6,901) ncycle,npts,(xin(i),i=1,npts)                 ABJT4931
            stop ' error in -sort-'                                     ABJT4932
         endif                                                          ABJT4933
                                                                        ABJT4934
c ...    if any of the elements have been swapped, try again            ABJT4935
         if ( flag ) goto 15                                            ABJT4936
                                                                        ABJT4937
      return                                                            ABJT4938
                                                                        ABJT4939
c ... format statements                                                 ABJT4940
                                                                        ABJT4941
  901 format (////'  error in -sort-:  ncycle =',i5,4x,'npts =',        ABJT4942
     2  i5,4x,'xin:'/50(t4,1p10e11.3/))                                 ABJT4943
                                                                        ABJT4944
      end                                                               ABJT4945
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  ABJT4946
                                                                        ABJT4947
      function voigt ( a,vv )                                           ABJT4948
c COPIED TO ModOpacityVoigt                                             ABJT4949
c ... this subroutine computes the Voigt function "h" as a function     ABJT4950
c     of "a" and "v".  See Mihalas, "Stellar Atmospheres", 1978, for    ABJT4951
c     definitions and details.                                          ABJT4952
                                                                        ABJT4953
      dimension atab(11),vtab(11),htable(11,11)                         ABJT4954
                                                                        ABJT4955
      data atab /  1.00E-01,1.58E-01,2.51E-01,3.98E-01,6.31E-01,        ABJT4956
     2    1.00E+00,1.58E+00,2.51E+00,3.98E+00,6.31E+00,1.00E+01 /       ABJT4957
      data vtab /  0.00E+00,5.00E-01,1.00E+00,1.50E+00,2.00E+00,        ABJT4958
     2    2.50E+00,3.00E+00,3.50E+00,4.00E+00,4.50E+00,5.00E+00 /       ABJT4959
      data htable /                                                     ABJT4960
     2 8.96E-01,8.44E-01,7.69E-01,6.72E-01,5.54E-01,4.28E-01,3.08E-01,  ABJT4961
     3 2.10E-01,1.38E-01,8.83E-02,5.61E-02,7.18E-01,6.85E-01,6.38E-01,  ABJT4962
     4 5.72E-01,4.89E-01,3.91E-01,2.92E-01,2.04E-01,1.36E-01,8.78E-02,  ABJT4963
     5 5.60E-02,3.73E-01,3.74E-01,3.72E-01,3.63E-01,3.43E-01,3.05E-01,  ABJT4964
     6 2.50E-01,1.88E-01,1.30E-01,8.63E-02,5.56E-02,1.34E-01,1.48E-01,  ABJT4965
     7 1.66E-01,1.87E-01,2.05E-01,2.12E-01,1.98E-01,1.65E-01,1.22E-01,  ABJT4966
     8 8.39E-02,5.49E-02,4.02E-02,5.18E-02,6.85E-02,9.07E-02,1.17E-01,  ABJT4967
     9 1.40E-01,1.51E-01,1.40E-01,1.12E-01,8.07E-02,5.40E-02,1.47E-02,  ABJT4968
     1 2.19E-02,3.28E-02,4.86E-02,6.97E-02,9.38E-02,1.13E-01,1.17E-01,  ABJT4969
     2 1.02E-01,7.69E-02,5.29E-02,7.94E-03,1.25E-02,1.95E-02,3.01E-02,  ABJT4970
     3 4.55E-02,6.53E-02,8.53E-02,9.64E-02,9.11E-02,7.27E-02,5.16E-02,  ABJT4971
     4 5.34E-03,8.44E-03,1.33E-02,2.08E-02,3.21E-02,4.77E-02,6.58E-02,  ABJT4972
     5 7.97E-02,8.09E-02,6.83E-02,5.01E-02,3.92E-03,6.21E-03,9.81E-03,  ABJT4973
     6 1.54E-02,2.40E-02,3.63E-02,5.18E-02,6.62E-02,7.16E-02,6.39E-02,  ABJT4974
     7 4.85E-02,3.02E-03,4.79E-03,7.57E-03,1.19E-02,1.86E-02,2.85E-02,  ABJT4975
     8 4.17E-02,5.55E-02,6.33E-02,5.94E-02,4.69E-02,2.41E-03,3.81E-03,  ABJT4976
     9 6.03E-03,9.52E-03,1.49E-02,2.30E-02,3.42E-02,4.69E-02,5.59E-02,  ABJT4977
     1 5.51E-02,4.51E-02 /                                              ABJT4978
                                                                        ABJT4979
      data sqrtpi / 1.7725 /                                            ABJT4980
      logical lbug                                                      ABJT4981
                                                                        ABJT4982
c ********************************************************************  ABJT4983
c                                                                       ABJT4984
c                           begin execution                             ABJT4985
c                                                                       ABJT4986
c ********************************************************************  ABJT4987
                                                                        ABJT4988
c ... check debug option                                                ABJT4989
      call debug ( 'voigt',                                             ABJT4990
     &                      lbug )                                      ABJT4991
                                                                        ABJT4992
      v = abs( vv )                                                     ABJT4993
                                                                        ABJT4994
      if ( v.gt.5. .or. a.ge.10. ) then                                 ABJT4995
                                                                        ABJT4996
c ...    Lorentzian profile in limit of large a,v                       ABJT4997
         voigt = a / sqrtpi / ( a*a + v*v )                             ABJT4998
                                                                        ABJT4999
      else if ( a.le.0.1 ) then                                         ABJT5000
                                                                        ABJT5001
c ...    use Dawson's integral (see Mihalas) for small a                ABJT5002
         voigt = exp( -min(v*v,50.) ) + 2.*a/sqrtpi *                   ABJT5003
     &           ( 2.*v*dawson(v) - 1. )                                ABJT5004
                                                                        ABJT5005
      else                                                              ABJT5006
                                                                        ABJT5007
c ...    use a table lookup                                             ABJT5008
         alog = log10( a )                                              ABJT5009
         ia = alog/0.2 + 6                                              ABJT5010
         iv = v/0.5 + 1                                                 ABJT5011
         fa = alog/0.2 - (ia-6)                                         ABJT5012
         fv = v/0.5 - (iv-1)                                            ABJT5013
         g1 = (1.-fa)*htable(ia,iv) + fa*htable(ia+1,iv)                ABJT5014
         g2 = (1.-fa)*htable(ia,iv+1) + fa*htable(ia+1,iv+1)            ABJT5015
         voigt = (1.-fv)*g1 + fv*g2                                     ABJT5016
                                                                        ABJT5017
      endif                                                             ABJT5018
                                                                        ABJT5019
      if ( lbug ) write (6,901) a,v,voigt,g1,g2,ia,iv                   ABJT5020
  901 format (t2,'debug output from -voigt-'/t4,'a',t18,'v',t32,        ABJT5021
     2  'voigt',t46,'g1',t60,'g2',t74,'ia',t80,'iv'/t4,1p5e14.4,        ABJT5022
     3  0p,t74,2i6)                                                     ABJT5023
                                                                        ABJT5024
      return                                                            ABJT5025
      end                                                               ABJT5026
c ============== END PROGRAM ========================================== ABJT5027
                                                                        ABJT****
