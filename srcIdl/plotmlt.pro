;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
pro plotmlt, maxran, white = white, black = black, 		$
      no00 = no00, no06 = no06, no12 = no12, no18 = no18,	$
      longs = longs, dash = dash

  if n_elements(white) gt 0 then color = 255
  if n_elements(black) gt 0 then color = 0
  if n_elements(color) eq 0 then begin
    if !d.name eq 'PS' then color = 0 else color = 255
  endif

  if n_elements(no00) eq 0 then no00 = 0
  if n_elements(no06) eq 0 then no06 = 0
  if n_elements(no12) eq 0 then no12 = 0
  if n_elements(no18) eq 0 then no18 = 0

  if n_elements(dash) eq 0 then dash = 1 else dash = 2

  if n_elements(longs) eq 0 then begin
    p00 = '00'
    p06 = '06'
    p12 = '12'
    p18 = '18'
  endif else begin
    p00 = '00'
    p06 = '90'
    p12 = '180'
    p18 = '270'
  endelse

  t = findgen(182.0)*2.0*!pi/180.0
  xp = cos(t)
  yp = sin(t)

  plots, maxran*xp, maxran*yp, color = color
  for i=10,maxran, 10 do					$
    oplot, float(i)*xp, float(i)*yp,linestyle=dash, color = color

  oplot, [-maxran,maxran],[0.0,0.0], linestyle =dash, color = color
  oplot, [0.0,0.0], [-maxran,maxran], linestyle = dash, color = color

  xs  = float(!d.x_size)
  ys  = float(!d.y_size)
  pos = float(!p.clip(0:3))
  pos([0,2]) = pos([0,2])/xs
  pos([1,3]) = pos([1,3])/ys

  mid_x = (pos(2) + pos(0))/2.0
  mid_y = (pos(3) + pos(1))/2.0

  ch = 0.8

  y_ch = ch*float(!d.y_ch_size)/ys
  x_ch = ch*float(!d.x_ch_size)/xs

  if no00 eq 0 then 							$
    xyouts, mid_x, pos(1) - y_ch*1.1, p00, alignment=0.5, 		$
      charsize=ch, /norm

  if no06 eq 0 then 							$
    xyouts, pos(2)+x_ch*0.15, mid_y - y_ch/2.0, p06, 			$
      charsize=ch, /norm

  if no12 eq 0 then 							$
    xyouts, mid_x, pos(3) + y_ch*0.15, p12, alignment=0.5, 		$
      charsize=ch, /norm

  if no18 eq 0 then 							$
    xyouts, pos(0)-x_ch*0.15, mid_y - y_ch/2.0, p18, 			$
      charsize=ch, /norm, alignment = 1.0

  return

end
