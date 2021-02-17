;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
pro c_s_to_a, timearray, strtime

  timearray = intarr(6)

  mon='JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC'

  timearray(2) = fix(strmid(strtime,0,2))
  smon = strmid(strtime,3,3)
  bs = byte(smon)
  loc = where((bs ge 97) and (bs lt 122), count)
  if count gt 0 then bs(loc) = bs(loc)-byte(32)
  smon = string(bs)

  for j=0,11 do 							      $
    if strmid(mon,j*3,3) eq smon then timearray(1)=j+1
  timearray(0)=fix(strmid(strtime,7,2))+1900
  if timearray(0) lt 1965 then timearray(0) = timearray(0) + 100

  if (strlen(strmid(strtime,10,2)) gt 0) and			$
     (strmid(strtime,10,2) ne '  ') then			$
    timearray(3)=fix(strmid(strtime,10,2))			$
  else timearray(3) = 0
  if (strlen(strmid(strtime,13,2)) gt 0) and			$
     (strmid(strtime,13,2) ne '  ') then			$
    timearray(4)=fix(strmid(strtime,13,2))			$
  else timearray(4) = 0
  if (strlen(strmid(strtime,16,6)) gt 0) and			$
     (strmid(strtime,16,2) ne '  ') then			$
    timearray(5)=fix(strmid(strtime,16,6))			$
  else timearray(5) = 0

  RETURN

END

