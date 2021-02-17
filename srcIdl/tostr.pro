;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
function tostr,value
  return, strcompress(string(long(value)),/remove_all)
end
