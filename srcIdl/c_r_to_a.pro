;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
PRO c_r_to_a, timearray, timereal

  dayofmon = [31,28,31,30,31,30,31,31,30,31,30,31]

  timearray = intarr(6)

  speryear = double(31536000.0)
  sperday  = double(86400.0)
  sperhour = double(3600.0)
  spermin  = double(60.0)

  numofyears = floor(timereal/speryear)
  if (numofyears+65) mod 4 eq 0 then dayofmon(1) = dayofmon(1) + 1
  numofdays = floor((timereal mod speryear)/sperday)
  numofleap = floor(numofyears / 4)
  numofdays = numofdays - numofleap
  if numofdays lt 0 then begin
    if (numofyears+65) mod 4 eq 0 then dayofmon(1) = dayofmon(1) - 1
    numofyears = numofyears - 1
    numofdays = numofdays + numofleap + 365
    if (numofyears+65) mod 4 eq 0 then dayofmon(1) = dayofmon(1) + 1
    numofleap = floor(numofyears / 4)
    numofdays = numofdays - numofleap
  endif
  numofhours = floor((timereal mod sperday)/sperhour)
  numofmin = floor((timereal mod sperhour)/spermin)
  numofsec = floor(timereal mod spermin)

  numofmon = 0

  while numofdays ge dayofmon(numofmon) do begin

    numofdays = numofdays - dayofmon(numofmon)
    numofmon = numofmon + 1

  endwhile

  timearray(0) = numofyears + 1965
  timearray(1) = numofmon + 1
  timearray(2) = numofdays + 1
  timearray(3) = numofhours
  timearray(4) = numofmin
  timearray(5) = numofsec

  RETURN

END
