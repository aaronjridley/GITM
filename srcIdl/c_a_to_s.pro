;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
pro c_a_to_s, timearray, strtime

  mon='JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC' 

  sd = '0'+tostr(timearray(2))
  sd = strmid(sd,strlen(sd)-2,2)
  sm = strmid(mon,(timearray(1)-1)*3,3)
  if timearray(0) lt 1900 then year = timearray(0) 		$
  else year = timearray(0)-1900
  if (year ge 100) then year = year - 100
  sy = chopr('0'+tostr(year),2)
  sh = '0'+tostr(timearray(3))
  sh = strmid(sh,strlen(sh)-2,2)
  si = '0'+tostr(timearray(4))
  si = strmid(si,strlen(si)-2,2)
  ss = '0'+tostr(timearray(5))
  ss = strmid(ss,strlen(ss)-2,2)

  strtime = sd+'-'+sm+'-'+sy+' '+sh+':'+si+':'+ss+'.000'

  RETURN

END

