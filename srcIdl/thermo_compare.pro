
if (n_elements(filelist) eq 0) then begin
   filelist = findfile('*.bin')
   filelist = filelist(0)
endif

filelist = ask('perturbation filename to plot',filelist(0))

if (n_elements(filelist_base) eq 0) then filelist_base = filelist(0)
filelist_base = ask('baseline filename to plot',filelist_base(0))

filelistp = findfile(filelist)
filelist0 = findfile(filelist_base)

nFiles  = n_elements(filelistp)
nFiles0 = n_elements(filelist0)

if (nFiles ne nFiles0) then begin
   print, "Number of files are not equal"
   stop
endif

for iFile = 0,nFiles-1 do begin

   file      = filelistp(iFile)
   file_base = filelist0(iFile)

   print, "Reading file : ",file
   read_thermosphere_file, file, nvars, nalts, nlats, nlons, $
                           vars, data, rb, cb, bl_cnt, itime1

   print, "Reading file : ",file_base
   read_thermosphere_file, filelist_base, nvars_base, $
                           nalts_base, nlats_base, nlons_base, $
                           vars_base, data_base, $
                           rb_base, cb_base, bl_cnt_base, itime2

   filename = filelist(0)

   alt = reform(data(2,*,*,*)) / 1000.0
   lat = reform(data(1,*,*,*)) / !dtor
   lon = reform(data(0,*,*,*)) / !dtor

   vars = [vars,'TEC']
   nVars++

   if (iFile eq 0) then begin

      for i=0,nvars-1 do print, tostr(i)+'. '+vars(i)
      if (n_elements(sel) eq 0) then sel = 3
      sel = fix(ask('which var to plot',tostr(sel)))

      if (n_elements(plotlog) eq 0) then plotlog = '0'
      plotlog = fix(ask('whether you want log or not (0/1)',tostr(plotlog)))

      if (n_elements(psfile) eq 0) then psfile = filename+'.ps'
      psfile = ask('ps file name',psfile)

      for i=0,nalts-1 do print, tostr(i)+'. '+string(alt(2,2,i))
      if (n_elements(selset) eq 0) then begin
         d = abs(reform(alt(2,2,*)-400.0))
         l = where(d eq min(d))
         selset = l(0)
      endif
      selset = fix(ask('which altitude to plot',tostr(selset)))

      if (n_elements(percentage) eq 0) then percentage = 0
      percentage = ask('percentage variation (0 for no, 1 for yes)',tostr(percentage))

      if (n_elements(smini) eq 0) then smini = '0.0'
      if (n_elements(smaxi) eq 0) then smaxi = '0.0'
      smini = ask('minimum (0.0 for automatic)',smini)
      smaxi = ask('maximum (0.0 for automatic)',smaxi)

      if (n_elements(plotvector) eq 0) then plotvector = 1
      plotvector = fix(ask('whether you want vectors or not (0/1)',tostr(plotvector)))

      if (n_elements(vi_cnt) eq 0) then vi_cnt = 0
      if (plotvector) then begin

         vi_cnt = fix(ask('ions (1) or neutrals (0)',tostr(vi_cnt))) 
         ; vi_cnt is whether to plot vectors of Vi
         ; vn_cnt is whether to plot vectors of Vn
         vn_cnt = 1-vi_cnt
         print,'-1  : automatic selection'
         factors = [1.0, 5.0, 10.0, 20.0, 25.0, $
                    50.0, 75.0, 100.0, 150.0, 200.0]
         nfacs = n_elements(factors)
         for i=0,nfacs-1 do print, tostr(i)+'. '+string(factors(i)*10.0)
         vector_factor = fix(ask('velocity factor','-1'))
      endif else vector_factor = 0

; cursor position variables, which don't matter at this point
      cursor_x = 0.0
      cursor_y = 0.0
      strx = '0.0'
      stry = '0.0'

;cnt1 is a lat/lon plot
      cnt1 = 1

;cnt2 is a lat/alt plot
      cnt2 = 0

;cnt3 is a lon/alt plot
      cnt3 = 0

; yes is whether ghostcells are plotted or not:
      yes = 0
      no  = 1

; yeslog is whether variable should be logged or not:
      if (plotlog) then begin 
         yeslog = 1
         nolog  = 0
      endif else begin
         yeslog = 0
         nolog = 1
      endelse

; yeswrite_cnt is whether we have to output to a ps file or not.
      yeswrite_cnt = 1

; polar is variable to say whether we have polar plots or not
      polestring = mklower(ask('whether you would like polar or flat (p or f)','p'))
      if (strpos(polestring,'p') eq 0) then polar = 1 else polar = 0

; npolar is whether we are doing the northern or southern hemisphere
      if (n_elements(npolar) eq 0) then npolar = 1
      if (polar) then begin
         npolar = mklower(ask('north (1) or south (0) pole',tostr(npolar)))
      endif

; MinLat is for polar plots:
      if (n_elements(MinLat) eq 0) then MinLat = 40.0
      if (polar) then begin
         MinLat = float(ask('minimum latitude to plot',tostr(minlat)))
      endif

; showgridyes says whether to plot the grid or not.
      showgridyes = 0

;plotvectoryes says whether to plot vectors or not
      plotvectoryes = plotvector

; number of points to skip when plotting vectors:
      step = 2

      cursor_cnt = 0

      xrange = [0.0,0.0]

      yrange = [0.0,0.0]

      if (vi_cnt eq 0) then test = "V!Dn!N(east)" $
      else test = "V!Di!N(east)"

      for i=0,nvars-1 do begin
         var=strcompress(vars[i],/remove_all)
         tes = strcompress(test,/remove_all)
         result= STRCMP( var, tes,8 )
         if (result eq 1) then vneast_index = i
      endfor

   sel_save = sel

   endif

   sel=sel_save

   data_sub = data
;   if (min(data_sub(sel,*,*,*)) gt 0.0) then $
;      data_sub(sel,*,*,*) = 100.0*(data_sub(sel,*,*,*) - data_base(sel,*,*,*))/data_base(sel,*,*,*) $
;   else data_sub(sel,*,*,*) = data_sub(sel,*,*,*) - data_base(sel,*,*,*)

   if (sel eq nVars-1) then begin

      dalt = (alt(*,*,1:nAlts-1) - alt(*,*,0:nAlts-2))*1000.0

      e   = reform(data_base(33,*,*,*))
      dalte = dalt * e(*,*,0:nAlts-2)

      tec_base = fltarr(nLons, nLats)
      for i=2,nAlts-3 do tec_base = tec_base + dalte(*,*,i)/ 1.0e16

      e   = reform(data_sub(33,*,*,*))
      dalte = dalt * e(*,*,0:nAlts-2)

      tec_sub = fltarr(nLons, nLats)
      for i=2,nAlts-3 do tec_sub = tec_sub + dalte(*,*,i)/ 1.0e16

      sel = 33
      vars(sel) = 'TEC'

      for i=0,nAlts-1 do begin
         data_sub(sel,*,*,i) = tec_sub
         data_base(sel,*,*,i) = tec_base
      endfor

   endif

   data_sub(sel,*,*,*) = data_sub(sel,*,*,*) - data_base(sel,*,*,*)

   if (percentage eq 1) then $
      data_sub(sel,*,*,*) = 100.0*data_sub(sel,*,*,*)/data_base(sel,*,*,*)

   data_sub(vneast_index,*,*,*) = $
      data_sub(vneast_index,*,*,*) - data_base(vneast_index,*,*,*)
   data_sub(vneast_index+1,*,*,*) = $
      data_sub(vneast_index+1,*,*,*) - data_base(vneast_index+1,*,*,*)

   if (nFiles gt 1) then begin
      p = strpos(psfile,'.ps')
      if (p gt -1) then psfile = strmid(psfile,0,p)
      psfile_final = psfile+'_'+chopr('000'+tostr(iFile),4)+'.ps'
   endif else begin
      psfile_final = psfile
   endelse
   
   thermo_plot,cursor_x,cursor_y,strx,stry,step,nvars,sel,nfiles,	$
               cnt1,cnt2,cnt3,yes,no,yeslog,  	  $
               1-yeslog,nalts,nlats,nlons,yeswrite_cnt,$
               polar,npolar,MinLat,showgridyes,	  $
               plotvectoryes,vi_cnt,vn_cnt,vector_factor,	  $
               cursor_cnt,data_sub,alt,lat,lon,	  $
               xrange,yrange,selset,smini,smaxi,	  $
               filename,vars, psfile_final, 0, 'mid', itime1, $
               /portrait
;, $
;                     plotpole = 1


;       thermo_plot_new,cursor_x,cursor_y,strx,stry,step,nvars,sel,nfiles,$
;                   cnt1,cnt2,cnt3,yes,no,yeslog,  	  $
;                   1-yeslog,nalts,nlats,nlons,yeswrite_cnt,$
;                   polar,npolar,MinLat,showgridyes,	  $
;                   plotvectoryes,vi_cnt,vn_cnt,vector_factor,	  $
;                   cursor_cnt,data_sub,alt,lat,lon,	  $
;                   xrange,yrange,selset,smini,smaxi,	  $
;                   filename,vars, psfile_final, 0, 'mid', itime1 ;, $
;;                       plotsquare = plotsquare



endfor

end
