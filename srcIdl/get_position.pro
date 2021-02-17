;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
;
; get_position
;
; used in conjunction with pos_space. Determines the position of the current
; plotting region, given the output parameters from pos_space.
;
; Input parameters:
; nb, space, bs, nbx, nby, xoff, yoff, xf, yf - Outputs from pos_space
; pos_num - the number of the plot, ranges from 0 : bs-1
;
; Output parameters:
;
; pos - the position of the plot, used in the plot command
;
; modified to make rectangles on Jan 2, 1998

pro get_position, nb, space, sizes, pos_num, pos, rect = rect,		$
		  xmargin = xmargin, ymargin = ymargin

  xipos = fix(pos_num) mod sizes.nbx
  yipos = fix(pos_num)/sizes.nbx

  yf2 = sizes.yf
  yf = sizes.yf*(1.0-space)
  xf2 = sizes.xf
  xf = sizes.xf*(1.0-space)

  if n_elements(rect) gt 0 then begin

    if n_elements(xmargin) gt 0 then xmar = xmargin(0) 			$
    else xmar = space/2.0

    if n_elements(ymargin) gt 0 then ymar = ymargin(0) 			$
    else ymar = space/2.0

    xtotal = 1.0 - (space*float(sizes.nbx-1) + xmar + xf2*space/2.0)
    xbs = xtotal/(float(sizes.nbx)*xf)

    xoff = xmar - xf2*space/2.0

    ytotal = 1.0 - (space*float(sizes.nby-1) + ymar + yf2*space/2.0)
    ybs = ytotal/(float(sizes.nby)*yf)

    yoff = 0.0

  endif else begin

    xbs  = sizes.bs
    xoff = sizes.xoff
    ybs  = sizes.bs
    yoff = sizes.yoff

  endelse

  xpos0 = float(xipos) * (xbs+space)*xf + xoff + xf2*space/2.0
  xpos1 = float(xipos) * (xbs+space)*xf + xoff + xf2*space/2.0 + xbs*xf

  ypos0 = (1.0-yf2*space/2) - (yipos * (ybs+space)*yf + ybs*yf) - yoff
  ypos1 = (1.0-yf2*space/2) - (yipos * (ybs+space)*yf) - yoff

  pos= [xpos0,ypos0,xpos1,ypos1]

  RETURN

END

