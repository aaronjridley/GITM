!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

Module ModEIEFiles

  use ModCharSize
  
  character (len=iCharLenIE_) :: ihp_file              = 'hpke.noaa'
  character (len=iCharLenIE_) :: pem_file              = 'hpke2.pem'
  character (len=iCharLenIE_) :: izmem_file            = 'iz94.cofcnts'
  character (len=iCharLenIE_) :: hepner_maynard_file   = 'hmr89.cofcnts'
  character (len=iCharLenIE_) :: millstone_hill_i_file = 'mhi.cofcnts'
  character (len=iCharLenIE_) :: millstone_hill_s_file = 'mhs.cofcnts'
  character (len=iCharLenIE_) :: stat_amie_file        = 'amie.ascii'
  character (len=iCharLenIE_) :: weimer96_file         = 'wei96.cofcnts'
  character (len=iCharLenIE_) :: weimer01_file         = 'wei01.cofcnts'
  character (len=iCharLenIE_) :: AMIEFileNorth, AMIEFileSouth

end Module ModEIEFiles
