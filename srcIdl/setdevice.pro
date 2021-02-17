;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
pro setdevice, psfile, orient, psfont, percent, eps=eps, 	$
	       psname_inq = psname_inq

  if n_elements(psfile) eq 0 then begin

    psfile = ''
    if n_elements(psname_inq) gt 0 then begin
      read, 'Enter ps filename : ',psfile
    endif
    if strlen(psfile) eq 0 then psfile = 'idl.ps'

  endif

  if n_elements(percent) eq 0 then percent = 1.0		$
  else if percent gt 1.0 then percent = float(percent)/100.0
  if n_elements(orient) eq 0 then orient = 'landscape'
  if n_elements(psfont) eq 0 then psfont = 28
  if n_elements(eps) eq 0 then eps = 0 else eps = 1
  set_plot, 'ps', /copy, /interpolate

  !p.font = 0

  if (strmid(orient,0,1) eq 'p') or (strmid(orient,0,1) eq 'P') then begin

    changep = percent
    xs = 7.5
    ys = 9.5

    if eps eq 0 then begin

      case (psfont) of

	0  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Courier 
	1  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Courier, /Bold 
    	2  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Courier, /Oblique 
	3  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Courier, /Bold, /Oblique
       	4  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Helvetica
      	5  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Helvetica, /Bold
    	6  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Helvetica, /Oblique
       	8  : device, file = psfile, /color, bits=8,      	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Helvetica, /Bold, /Oblique 
    	12 : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Avantgarde, /Book 
     	13 : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Avantgarde, /Book, /Oblique
	 14 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Avantgarde, /Demi 
      	15 : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Avantgarde, /Demi, /Oblique
       	20 : device, file = psfile, /color, bits=8, $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Schoolbook
   	21 : device, file = psfile, /color, bits=8,$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Schoolbook, /Bold
      	22 : device, file = psfile, /color, bits=8,$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Schoolbook, /Italic
       	23 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Schoolbook, /Bold, /Italic 
	28 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Times
	29 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Times, /Bold
	 30 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Times, /Italic
	31 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Times, /Bold, /Italic

      endcase

    endif else begin

      case (psfont) of

	0  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Courier, /encapsulated 
	1  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Courier, /Bold, /encapsulated 
    	2  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Courier, /Oblique, /encapsulated 
	3  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Courier, /Bold, /Oblique, /encapsulated
       	4  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Helvetica, /encapsulated
      	5  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Helvetica, /Bold, /encapsulated
    	6  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Helvetica, /Oblique, /encapsulated
       	8  : device, file = psfile, /color, bits=8,      	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Helvetica, /Bold, /Oblique, /encapsulated 
    	12 : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Avantgarde, /Book, /encapsulated 
     	13 : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Avantgarde, /Book, /Oblique, /encapsulated
	 14 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Avantgarde, /Demi, /encapsulated 
      	15 : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Avantgarde, /Demi, /Oblique, /encapsulated
       	20 : device, file = psfile, /color, bits=8, $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Schoolbook, /encapsulated
   	21 : device, file = psfile, /color, bits=8,$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Schoolbook, /Bold, /encapsulated
      	22 : device, file = psfile, /color, bits=8,$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Schoolbook, /Italic, /encapsulated
       	23 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Schoolbook, /Bold, /Italic, /encapsulated 
	28 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Times, /encapsulated
	29 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Times, /Bold, /encapsulated
	 30 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Times, /Italic, /encapsulated
	31 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Times, /Bold, /Italic, /encapsulated 

      endcase

    endelse

  endif else begin

    xs = 10.0
    ys = 7.0
    change = percent

    if eps eq 0 then begin

      case (psfont) of

	0  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Courier 
	1  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Courier, /Bold 
    	2  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Courier, /Oblique 
	3  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Courier, /Bold, /Oblique
       	4  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Helvetica
      	5  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Helvetica, /Bold
    	6  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Helvetica, /Oblique
       	8  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Helvetica, /Bold, /Oblique 
    	12 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Avantgarde, /Book 
     	13 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Avantgarde, /Book, /Oblique
	14 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
	 	/Avantgarde, /Demi 
      	15 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Avantgarde, /Demi, /Oblique
       	20 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Schoolbook
   	21 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Schoolbook, /Bold
      	22 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Schoolbook, /Italic
       	23 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Schoolbook, /Bold, /Italic 
	28 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Times
	29 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Times, /Bold
	30 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Times, /Italic
	31 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Times, /Bold, /Italic

      endcase

    endif else begin

      case (psfont) of

	0  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Courier, /encapsulated  
	1  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Courier, /Bold, /encapsulated  
    	2  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Courier, /Oblique, /encapsulated  
	3  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Courier, /Bold, /Oblique, /encapsulated 
       	4  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Helvetica, /encapsulated 
      	5  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Helvetica, /Bold, /encapsulated 
    	6  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Helvetica, /Oblique, /encapsulated 
       	8  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Helvetica, /Bold, /Oblique, /encapsulated  
    	12 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Avantgarde, /Book, /encapsulated  
     	13 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Avantgarde, /Book, /Oblique, /encapsulated 
	14 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
	 	/Avantgarde, /Demi, /encapsulated  
      	15 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Avantgarde, /Demi, /Oblique, /encapsulated 
       	20 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Schoolbook, /encapsulated 
   	21 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Schoolbook, /Bold, /encapsulated 
      	22 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Schoolbook, /Italic, /encapsulated 
       	23 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Schoolbook, /Bold, /Italic, /encapsulated  
	28 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Times, /encapsulated 
	29 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Times, /Bold, /encapsulated 
	30 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Times, /Italic, /encapsulated 
	31 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Times, /Bold, /Italic, /encapsulated 

      endcase

    endelse

  endelse

  return

end

