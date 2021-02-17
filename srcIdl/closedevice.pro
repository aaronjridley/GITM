;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
pro closedevice

  if !d.name eq 'PS' then begin
    device, /close
    set_plot, 'X'
  endif

  return

end
