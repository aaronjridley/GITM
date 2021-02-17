;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
;
; function mm
;
;  computes min and max of an array and returns them in an array
;

function mm, array
return, [min(array),max(array)]
end
