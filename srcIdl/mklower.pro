;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
function mklower, string

  temp = byte(string)

  loc = where((temp ge 65) and (temp le 90), count)

  if count ne 0 then temp(loc) = temp(loc)+32

  return, string(temp)

end

