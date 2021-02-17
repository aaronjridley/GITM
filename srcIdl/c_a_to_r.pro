;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
PRO c_a_to_r, timearray, timereal

  dayofmon = [31,28,31,30,31,30,31,31,30,31,30,31]
  if ((timearray(0) mod 4) eq 0) then dayofmon(1) = dayofmon(1) + 1

  timereal = double(0.0)

  if timearray(0) lt 65 then timearray(0) = timearray(0) + 2000
  if timearray(0) gt 1900 then numofyears = timearray(0)-1965 		      $
  else numofyears = timearray(0)-65	
  numofleap = floor(float(numofyears)/4.0)
  numofmonths = timearray(1) - 1
  numofdays = 0

  for i = 0, numofmonths-1 do begin

    numofdays = numofdays + dayofmon(i)

  endfor

  numofdays = numofdays + timearray(2) - 1
  numofhours = timearray(3)
  numofminutes = timearray(4)
  numofseconds = timearray(5)

  timereal = double(numofseconds*1.0) +       $
	     double(numofminutes*60.0) +             $
	     double(numofhours*60.0*60.0) +          $
	     double(numofdays*24.0*60.0*60.0) +      $
	     double(numofleap*24.0*60.0*60.0) +      $
	     double(numofyears*365.0*24.0*60.0*60.0)

  RETURN

END

