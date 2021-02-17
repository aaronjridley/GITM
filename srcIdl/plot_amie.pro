
;--------------------------------------------------------------
; Get Inputs from the user
;--------------------------------------------------------------

if (n_elements(initial_guess) eq 0) then begin
   initial_guess = findfile('-t b*')
   iFile = 0
   help, initial_guess,iFile
   while ((strpos(initial_guess(iFile),"data") gt 0 or $
           strpos(initial_guess(iFile),"ps") gt 0 or $
           strpos(initial_guess(iFile),"sum") gt 0) and $
          iFile lt n_elements(initial_guess)-1) do begin
      iFile = iFile + 1
      print, iFile
   endwhile
   print, iFile
   initial_guess = initial_guess(iFile)
   if strlen(initial_guess) eq 0 then initial_guess='b970101'
endif

amie_file = ask('AMIE binary file name',initial_guess)
psfile = ask('ps file',amie_file+'.ps')

ppp = ask('plots per page','12')

read_amie_binary, amie_file, data, lats, mlts, time, fields, 		$
                  imf, ae, dst, hp, cpcp, version

nFields = n_elements(fields)
nLats = n_elements(lats)
nMlts = n_elements(mlts)
nTimes = n_elements(time)

; data is defined a:

; data(number of time, number of fields, number of mlts, number of latitudes)
; so, if you want the potential:

; potential = reform(data(*,0,*,*))

; if you want a specific time, then :

; pot_at_time_4 = reform(potential(4,*,*))

; if you want a time range (say from time 10 to time 100):

; pot_range = potential(10:100,*,*)

;--------------------------------------------------------------
; figure out grid:
;--------------------------------------------------------------

nlats = n_elements(lats)
nmlts = n_elements(mlts)

lat2d = fltarr(nmlts,nlats)
lon2d = fltarr(nmlts,nlats)

for i=0,nlats-1 do lon2d(*,i) = mlts*!pi/12.0 - !pi/2.0
for i=0,nmlts-1 do lat2d(i,*) = lats

x = (90.0-lat2d)*cos(lon2d)
y = (90.0-lat2d)*sin(lon2d)

;--------------------------------------------------------------
; Need to figure out what to plot, so list to fields to the user
;--------------------------------------------------------------

for i=0,nfields-1 do print, tostr(i+1)+'. '+fields(i)

;--------------------------------------------------------------
; Get field to be contoured
;--------------------------------------------------------------

type_1 = fix(ask('field to plot','1'))-1
if (type_1 lt 0) or (type_1 gt nfields-1) then type_1 = 0

;--------------------------------------------------------------
; Get field to be under contour as an image
;--------------------------------------------------------------

type_2 = fix(ask('second field to plot (under first, 0 for none)',tostr(type_1+1)))-1
if (type_2 gt nfields-1) then type_2 = -1

if (type_2 gt -1) then $
   if ((strpos(mklower(fields(type_1)),'east') gt -1) or 		$
       (strpos(mklower(fields(type_1)),'north') gt -1)) and		$
      ((strpos(mklower(fields(type_2)),'east') gt -1) or 		$
       (strpos(mklower(fields(type_2)),'north') gt -1)) then begin
        print,'You have selected 2 vector quantities.'
        vect = ask('whether you would like vectors or contours (v/c)','v')
        if (strpos(mklower(vect),'v') gt -1) then vect = 1 else vect = 0

        if vect eq 1 then begin
          print, 'You can select a quantity to plot under the vectors.'
          print, 'Enter 0 for a magnitude of the vectors.'
    	  type_3 = fix(ask('-1 for nothing','0'))-1
          if (type_3 gt nfields-1) then type_3 = -1
          if (strpos(mklower(fields(type_1)),'east') gt -1) then east = 1 $
            else east = 0
        endif
   endif else begin
     vect = 0
     type_3 = -1
   endelse $
else begin
   vect = 0
   type_3 = -1
endelse

;--------------------------------------------------------------
; need to sort of rearrange things if vectors are selected.
; Want to move color image into data_2, east into data_1,
; and north into data_3.  What a pain.
;--------------------------------------------------------------

if (vect) then begin

  temp_2 = type_2
  type_2 = type_3
  type_3 = temp_2

  if (not east) then begin
    temp_1 = type_1
    type_1 = type_3
    type_3 = temp_1
  endif

endif

;--------------------------------------------------------------
; Get start time and end time, with defaults as the file times
;--------------------------------------------------------------

c_r_to_a, itime_start, min(time)
c_r_to_a, itime_end, max(time)

c_a_to_s, itime_start, stime_start
c_a_to_s, itime_end, stime_end

start_time = ask('start time of plotting',strmid(stime_start,0,15))

if (strlen(start_time) lt 9) then $
  start_time = strmid(stime_start,0,9)+' '+start_time

;--------------------------------------------------------------
; I got sick of typing in the ending date, so if the date is
; the same, assume the user just wants to enter the time
;--------------------------------------------------------------

sdate = strmid(start_time,0,9)
if strpos(stime_end,sdate) gt -1 then 					$
  end_time_default = strmid(stime_end,10,5)				$
else end_time_default = strmid(stime_end,0,15)

end_time   = ask('end time of plotting',end_time_default)

;--------------------------------------------------------------
; If the user entered a short string, assume it is just a time
; and add the date on the front
;--------------------------------------------------------------

if (strlen(end_time) lt 9) then end_time = strmid(start_time,0,9)+' '+end_time

;--------------------------------------------------------------
; Now figure out where in the file these things are, with the
; default to give the user everything
;--------------------------------------------------------------

c_s_to_a, itime_start, start_time
c_a_to_r, itime_start, rtime_start

c_s_to_a, itime_end, end_time
c_a_to_r, itime_end, rtime_end

n_start = where(time ge rtime_start)
if n_start(0) ge 0 then n_start = n_start(0) else n_start = 0

n_end = where(time ge rtime_end)
if n_end(0) ge 0 then n_end = n_end(0) else n_end = n_elements(time)-1

;--------------------------------------------------------------
; Now, allow the user to skip times
;--------------------------------------------------------------

print, 'You have selected '+tostr(n_end-n_start+1)+' times to plot.'
nskip = ask('number of times to skip between each plot','0')+1

;--------------------------------------------------------------
; See if the user wants to plot data locations on the background
;--------------------------------------------------------------

plot_data = 0

if (not vect) then begin

  que=''
  read, 'Would you like Data to be plotted under the convection? [n] ',que

  if (mklower(que) eq 'y') then begin
    datafile = ask('data file (with full path)',amie_file+'_data')
    dt = float(ask('half window size (minutes)','05'))*60.0
    plot_data = 1
  endif

endif

;--------------------------------------------------------------
; See if the user wants to subtract off a background pattern
;--------------------------------------------------------------

back = ask('whether you want to subtract/plot a background pattern (y/n)','n')
if (strpos(mklower(back),'n') ge 0) then back = 0 else back = 1

if (back eq 1) then begin

  back = ask('whether to plot out only the background pattern (y/n)','n')
  if (strpos(mklower(back),'n') ge 0) then back = 1 else back = 2

  start_time = ask('start time of bakground',strmid(stime_start,0,15))

  ;--------------------------------------------------------------
  ; I got sick of typing in the ending date, so if the date is
  ; the same, assume the user just wants to enter the time
  ;--------------------------------------------------------------

  sdate = strmid(start_time,0,9)
  if strpos(stime_end,sdate) gt -1 then 				$
    end_time_default = strmid(stime_end,10,5)				$
  else end_time_default = strmid(stime_end,0,15)

  end_time   = ask('end time of background',end_time_default)

  ;--------------------------------------------------------------
  ; If the user entered a short string, assume it is just a time
  ; and add the date on the front
  ;--------------------------------------------------------------

  if (strlen(end_time) lt 9) then 					$
    end_time = strmid(start_time,0,9)+' '+end_time

  ;--------------------------------------------------------------
  ; Now figure out where in the file these things are, with the
  ; default to give the user everything
  ;--------------------------------------------------------------

  c_s_to_a, itime_start, start_time
  c_a_to_r, itime_start, rtime_start

  c_s_to_a, itime_end, end_time
  c_a_to_r, itime_end, rtime_end

  ut_title = 'Background time : '+strmid(start_time,10,5)+' to '+	$
             strmid(end_time,10,5)+' UT'

  n_back_start = where(time ge rtime_start)
  if n_back_start(0) ge 0 then n_back_start = n_back_start(0) 		$
  else n_back_start = 0

  n_back_end = where(time ge rtime_end)
  if n_back_end(0) ge 0 then n_back_end = n_back_end(0) 		$
  else n_back_end = n_elements(time)-1

  back_1 = fltarr(nmlts,nlats)
  back_2 = fltarr(nmlts,nlats)
  back_3 = fltarr(nmlts,nlats) 

  for n=n_back_start, n_back_end do begin

    back_1 = back_1 + reform(data(n,type_1,*,*))
    if type_2 ge 0 then back_2 = back_2 + reform(data(n,type_2,*,*))
    if type_3 ge 0 then back_3 = back_3 + reform(data(n,type_3,*,*))

  endfor

  back_1 = back_1 / (n_back_end - n_back_start + 1)
  back_2 = back_2 / (n_back_end - n_back_start + 1)
  back_3 = back_3 / (n_back_end - n_back_start + 1)

  if (back eq 2) then begin
    back_1  = -1.0*back_1
    back_2  = -1.0*back_2
    back_3  = -1.0*back_3
    data    = data*0.0
    n_start = 0
    n_end   = 0
  endif

endif else begin

  back_1 = fltarr(nmlts,nlats)*0.0
  back_2 = fltarr(nmlts,nlats)*0.0
  back_3 = fltarr(nmlts,nlats)*0.0

endelse

;--------------------------------------------------------------
; Put the contour data into data_1 array and get field name
;--------------------------------------------------------------

data_1 = reform(data(*,type_1,*,*))
field_1 = strcompress(fields(type_1))

if (vect) then begin
  data_3 = reform(data(*,type_3,*,*))
  field_1 = strmid(field_1,0,strpos(field_1,'('))
  data_2 = sqrt(data_1^2 + data_3^2)
endif

;--------------------------------------------------------------
; Put the image data into data_2, get field name, and try to
; figure out the units on the field
;--------------------------------------------------------------

joule = 0

if (type_2 ge 0 or vect) then begin

  if (vect) then begin
      type_2  = type_1
      field_2 = 'over Magnitude of '+field_1
  endif else begin
      data_2  = reform(data(*,type_2,*,*))
      field_2 = 'over '+fields(type_2)
  endelse

  if (strpos(mklower(fields(type_2)),'potential') gt -1) then 		$
    units = '(kV)' 
  if (strpos(mklower(fields(type_2)),'cond')  gt -1) then 		$
    units = '(mhos)' 
  if (strpos(mklower(fields(type_2)),'electric field') gt -1) then 	$
    units = '(mV/m)' 
  if (strpos(mklower(fields(type_2)),'current') gt -1 and $
      strpos(mklower(fields(type_2)),'field') lt 0) then $
    units = '(A/m!E2!N)' 
  if (strpos(mklower(fields(type_2)),'current') gt -1 and $
      strpos(mklower(fields(type_2)),'field') gt -1) then $
    units = '(!Mm!XA/m!E2!N)' 
  if (strpos(mklower(fields(type_2)),'joule') gt -1) then begin
;    units = 'log(mJ/m!E2!N)' 
    units = '(mJ/m!E2!N)' 
;    data_2 = alog10(data_2+0.001)
;    loc = where(data_2 lt -1.0, count)
;    if count gt 0 then begin
;        data_2(loc) = -1.0
;        if (type_2 eq type_1) then begin
;            data_1 = alog10(data_1+0.001)
;            data_1(loc) = -1.0
;        endif
;    endif
    joule = 1
  endif
endif else begin
    data_2  = data_1*0.0
    field_2 = ' '
    units = ' '
endelse

;--------------------------------------------------------------
; figure out contour levels:
;--------------------------------------------------------------

maxi_array = fltarr(n_end-n_start+1+nskip)
mini_array = fltarr(n_end-n_start+1+nskip)

iE = 0

if (not vect) then begin
    for i = n_start, n_end,nskip do begin
        mini_array(i-n_start:i-n_start+nskip-1) = min(data_1(i,*,*)-back_1)
        maxi_array(i-n_start:i-n_start+nskip-1) = max(data_1(i,*,*)-back_1)
        iE = i-n_start+nskip-1
    endfor
endif else begin
    for i = n_start, n_end,nskip do begin
        mini_array(i-n_start:i-n_start+nskip-1) = min(data_2(i,*,*)-back_2)
        maxi_array(i-n_start:i-n_start+nskip-1) = max(data_2(i,*,*)-back_2)
        iE = i-n_start+nskip-1
    endfor
endelse

mini_array = mini_array(0:iE)
maxi_array = maxi_array(0:iE)

nLevels = 31

mini_1  = min(mini_array)
maxi_1  = max(maxi_array)
if (mini_1 ge 0.0) then begin
  mini_1  = 0.0
  range_1 = maxi_1
endif else begin
    if (not joule) then range_1 = max([abs(mini_1),maxi_1])*2.0 $
      else range_1 = maxi_1 - mini_1
endelse
dc      = 10.0^fix(alog10(range_1/100.0))
factor  = 1.0
while (range_1 gt dc*(nLevels-1)*factor) do factor=factor+0.1
dc = factor*dc
dc = float(ask('contour level for '+field_1,strcompress(string(dc))))
if (mini_1) eq 0.0 then 						$
  levels_1 = findgen(nLevels)*dc					$
else levels_1 = (findgen(nLevels) - (nLevels-1)/2)*dc

;--------------------------------------------------------------
; color stuff is a tiny bit more complicated, since we can have
; nothing to plot, so we have nothing to plot, then we set the
; levels to all be below the data level.  This way we "max out"
; the color scale and have an all white area.
;--------------------------------------------------------------

for i = n_start, n_end, nskip do begin
  mini_array(i-n_start:i-n_start+nskip-1) = min(data_2(i,*,*)-back_2)
  maxi_array(i-n_start:i-n_start+nskip-1) = max(data_2(i,*,*)-back_2)
endfor
mini_2  = min(mini_array)
maxi_2  = max(maxi_array)
if (mini_2 ge 0.0) then begin
  mini_2 = 0.0
  range_2 = maxi_2
endif else begin
    if (not joule) then range_2 = max([abs(mini_2),maxi_2])*2.0 $
      else range_2 = maxi_2 - mini_2
endelse

if (range_2 gt 0.0 and not vect) then begin
  dc      = 10.0^fix(alog10(range_2/100.0))
  factor  = 1.0
  while (range_2 gt dc*(nLevels-1)*factor) do factor=factor+0.05
  dc = factor*dc
  dc = float(ask('contour level for '+strcompress(fields(type_2)),	$
                 strcompress(string(dc))))
  if (mini_2) eq 0.0 then 						$
    levels_2 = findgen(nLevels)*dc					$
  else levels_2 = (findgen(nLevels) - (nLevels-1)/2)*dc
endif else levels_2 = levels_1

;--------------------------------------------------------------
; Determine maximum range
;--------------------------------------------------------------

mr = fix(max(90.0-lats)/10.0)*10.0 + 10.0

mr = 90.0-fix(ask('minimum latitude to plot',tostr(90.0-mr)))

;--------------------------------------------------------------
; Set up device
;--------------------------------------------------------------

setdevice, psfile, 'p', 4, 0.95

;--------------------------------------------------------------
; Read color table. Blue to white to red for -/+ data,
; red to white for + only data.
;--------------------------------------------------------------

if (mini_2 lt 0.0 and not joule) then				$
  makect,'mid' else makect,'wyr'
ncolors = 255
;  readct, ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"	$
;else readct, ncolors, getenv("IDL_EXTRAS")+"white_red.ct"

clevels = (ncolors-1)*findgen(nLevels)/(nLevels-1) + 1

;--------------------------------------------------------------
; Set up plot sizes for the circles
;--------------------------------------------------------------

space = 0.005
pos_space, ppp, space, sizes
total_plots = sizes.nby*sizes.nbx

;--------------------------------------------------------------
; determine how much slop we have for the color table.
;--------------------------------------------------------------

get_position, ppp, space, sizes, 0, pos
x_left = pos(0)
get_position, ppp, space, sizes, sizes.nbx-1, pos
x_right = pos(2)
slop = 1.0 - (x_right - x_left)

;--------------------------------------------------------------
; If there is too little room for the color table, we need
; to make a little more room.
;--------------------------------------------------------------

if (slop gt 0.02) then begin
  if (x_right lt 0.98) then shift_left = 0.0 				$
  else shift_left = 0.02 - (1.0-x_right)
  ct_left  = 0.98
  ct_right = 1.00
  d_shift_left = 0.0
endif else begin
  space = space*2.0
  pos_space, ppp, space, sizes
  get_position, ppp, space, sizes, 0, pos
  x_left = pos(0)
  get_position, ppp, space, sizes, sizes.nbx-1, pos
  x_right = pos(2)
  shift_left = x_left
  d_shift_left = space*0.5
  ct_right = 1.00
  ct_left  = ct_right - shift_left - d_shift_left*(sizes.nbx-1)
endelse

if (back eq 2) then begin
  get_position, ppp, space, sizes, 0, pos
  ct_left = pos(2)+space
  ct_right = ct_left + 0.02
endif

;--------------------------------------------------------------
; Determine the character size in normalized coordinates
;--------------------------------------------------------------

dy  = float(!d.y_ch_size)/float(!d.y_size)

;--------------------------------------------------------------
; If we have a whole bunch of plots, scale the size down.
;--------------------------------------------------------------

if (ppp gt 12) then ch_size = 0.65 else ch_size = 1.0

;--------------------------------------------------------------
; Want to make sure we start with the 0th plot 
;--------------------------------------------------------------

plot_number = -1

;--------------------------------------------------------------
; Begin main loop
;--------------------------------------------------------------

for n=n_start, n_end, nskip do begin

;--------------------------------------------------------------
; Determine current plot number 
;--------------------------------------------------------------

  plot_number = (plot_number + 1) mod ppp

;--------------------------------------------------------------
; Determine current time
;--------------------------------------------------------------

  c_r_to_a, itime, time(n)
  c_a_to_s, itime, stime

;--------------------------------------------------------------
; If we are on a new page, plot the date and the title
;--------------------------------------------------------------

  if plot_number eq 0 then begin

    plotdumb
    plot_string = field_1 + field_2
    xyouts, 0.0, 1.04, plot_string,/norm
    xyouts, 0.0, 1.01, strmid(stime,0,9),/norm
    xyouts, 1.0, 1.04, "AMIE Version "+string(Version,format="(f5.2)"), $
      alignment = 1.0, /norm

    if (back eq 1) then							$
      xyouts, 1.0, 1.01, ut_title, /norm, alignment = 1.0

    if (back eq 2) then							$
      xyouts, mean(pos([0,2])), pos(1)-3.0*dy, ut_title, 		$
        /norm, alignment = 0.5

  endif

;--------------------------------------------------------------
; Figure out where on the page we should be
;--------------------------------------------------------------

  get_position, ppp, space, sizes, plot_number, pos

;--------------------------------------------------------------
; Move over to the left a little bit for the color bar
;--------------------------------------------------------------

  pos([0,2]) = pos([0,2]) - shift_left - 				$
                            d_shift_left*float(plot_number mod sizes.nbx)

;--------------------------------------------------------------
; Color the back ground first
;--------------------------------------------------------------

  loc = where(lat2d(0,*) gt 90.0-mr)

  d = reform(data_2(n,*,loc))-back_2(*,loc)
  l = where(d ge levels_2(nLevels-2),c)
  if (c gt 0) then d(l) = levels_2(nLevels-2)

;  contour, reform(data_2(n,*,loc))-back_2(*,loc),x(*,loc),y(*,loc),	$

  contour, d, x(*,loc),y(*,loc),	$
	/follow, xstyle = 5, ystyle = 5,				$
	xrange = [-mr,mr],yrange=[-mr,mr], levels = levels_2, 		$
	pos = pos, /noerase, /cell_fill, c_color = clevels

;--------------------------------------------------------------
; The plot the contours
;--------------------------------------------------------------

  if (not vect) then begin

    contour, reform(data_1(n,*,loc))-back_1(*,loc),x(*,loc),y(*,loc),	$
	/follow, xstyle = 5, ystyle = 5,				$
	xrange = [-mr,mr],yrange=[-mr,mr], levels = levels_1, 		$
	pos = pos, /noerase, c_linestyle = 3.0*(levels_1 lt 0.0), thick = 4

    if (plot_data) then							$
      plot_data_new, datafile, 0, 0, dt, 90.0-mr, /nodata, ut=time(n)

  endif else begin

    length = (levels_1(1) - levels_1(0))*5.0
    title = string(length, format='(f5.2)')
    normal = 2.0/(levels_1(1) - levels_1(0))

    plot_vectors_polar, reform(data_1(n,*,loc))-back_1(*,loc),		$
	reform(data_3(n,*,*))-back_3(*,loc), 				$
	90.0-lat2d(*,loc), lon2d(*,loc), normal, $
	/arrow, /no0, length = length, pos = [mr,mr], title = title

  endelse

;--------------------------------------------------------------
; Figure out where we are on the page, and whether we need to
; labels or not for the MLT grid
;--------------------------------------------------------------

  no00 = 1
  no06 = 1
  no12 = 1
  no18 = 1

  if (n_end-n lt sizes.nbx) or 						$
     (plot_number+1 gt total_plots-sizes.nbx) then no00 = 0
  if (plot_number mod sizes.nbx eq sizes.nbx-1) or 			$
     (n_end-n eq 0) then no06 = 0
  if (plot_number lt sizes.nbx) then no12 = 0
  if (plot_number mod sizes.nbx eq 0) or (n eq n_start) then no18 = 0

;--------------------------------------------------------------
; Don't draw the 06 if there is a color bar right next to you
;--------------------------------------------------------------

  if (plot_number+1 eq total_plots) and 				$
     (type_2 ge 0) then no06 = 1

  if (back eq 2) and (type_2 ge 0) then no06 = 1

;--------------------------------------------------------------
; Draw the MLT grid
;--------------------------------------------------------------

  plotmlt, mr, no00 = no00, no06 = no06, no12 = no12, no18 = no18

;--------------------------------------------------------------
; Tell us what time it is
;--------------------------------------------------------------

  if (time(n_end)-time(n_start) lt 48.0*3600.0) then begin
      xyouts, pos(0), pos(3)-0.75*dy, strmid(stime,10,5)+' UT',$
        /norm, charsize = ch_size
  endif else begin
      xc = (pos(2)+pos(0))/2
      xr = (pos(2)-pos(0))/2 * 1.01
      yc = (pos(3)+pos(1))/2
      yr = (pos(3)-pos(1))/2 * 1.01
      xp = xc - xr*sin(!pi/4)
      yp = yc + yr*sin(!pi/4)
      xyouts, xp, yp, strmid(stime,0,9), $
          /norm, charsize = ch_size*0.9, align = 0.5, orient = 45

      xp = xc + xr*sin(!pi/4)
      yp = yc + yr*sin(!pi/4)
      xyouts, xp, yp, strmid(stime,10,5)+' UT', $
          /norm, charsize = ch_size*0.9, align = 0.5, orient = -45
  endelse

;--------------------------------------------------------------
; Figure out whether we need to put the max and min in 
; exponential form or not.
;--------------------------------------------------------------

  mini = min(data_1(n,*,*)-back_1)
  maxi = max(data_1(n,*,*)-back_1)

  if (vect) then begin
    v = ((data_1(n,*,*)-back_1)^2.0 + (data_3(n,*,*)-back_3)^2.0)^0.5
    mini = min(v)
    maxi = max(v)
  endif

  order_of_mag = alog10(max([abs(mini),maxi]))

;--------------------------------------------------------------
; If we don't, put the max and min words plus the values.
; If we do, just put the values.
;--------------------------------------------------------------

  if order_of_mag ge 0.0 and order_of_mag lt 3.0 then begin
    mini_string = 'Min:'+tostr(round(mini))
    maxi_string = 'Max:'+tostr(round(maxi))
    up_shift    = 0.5*dy
  endif else begin
    if mini lt 0.0 then mini_string = string(mini,format = '(E9.2)')	$
    else mini_string = string(mini,format = '(E8.2)')
    maxi_string = string(maxi,format = '(E8.2)')
    up_shift    = 0.25*dy
  endelse

;--------------------------------------------------------------
; Output the max and min strings we just created
;--------------------------------------------------------------

  xyouts, pos(0), pos(1)+up_shift, mini_string, /norm, 			$
          charsize = 0.8*ch_size
  xyouts, pos(2), pos(1)+up_shift, maxi_string, /norm, 			$
          alignment=1.0, charsize = 0.8*ch_size

;--------------------------------------------------------------
; Plot the color bar if we need it.
;--------------------------------------------------------------

  if (plot_number+1 eq ppp) or						$
     (n+nskip gt n_end) and (type_2 ge 0) then begin
    ctpos = [ct_left, pos(1), ct_right, pos(3)]
    plotct, ncolors, ctpos, mm(levels_2), units, /right
  endif

endfor

closedevice

end



