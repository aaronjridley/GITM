
filelist = findfile("-t *.save")
if (strlen(filelist(0)) eq 0) then filelist = findfile("-t *.dat")

filelist = ask('filelist to plot',filelist(0))

filelist = findfile(filelist)
nFiles = n_elements(filelist)

for iFile = 1, nFiles-1 do begin

   filename = filelist(iFile)
   print, 'Reading file ',filename

   read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $
                           vars, data, rb, cb, bl_cnt, itime1

   if (iFile eq 1) then begin

      filename = filelist(iFile-1)
      print, 'Reading file ',filename

      read_thermosphere_file, filename, nvars_base, $
                              nalts_base, nlats_base, nlons_base, vars_base, $
                              data_base, rb_base, cb_base, bl_cnt_base, itime2

   endif

   alt = reform(data(2,*,*,*)) / 1000.0
   lat = reform(data(1,*,*,*)) / !dtor
   lon = reform(data(0,*,*,*)) / !dtor

   if (iFile eq 1) then begin

      for i=0,nvars-1 do print, tostr(i)+'. '+vars(i)
      sel = fix(ask('which var to plot','9'))

      plotlog = ask('whether you want log or not (y/n)','n')
      if (strpos(plotlog,'y') eq 0) then plotlog = 1 else plotlog = 0

      psfile = filename+'.ps'
      psfile = ask('ps file name',psfile)

      for i=0,nalts-1 do print, tostr(i)+'. '+string(alt(2,2,i))
      selset = fix(ask('which altitude to plot','0'))

      smini = ask('minimum (0.0 for automatic)','0.0')
      smaxi = ask('maximum (0.0 for automatic)','0.0')

      plotvector = ask('whether you want vectors or not (y/n)','y')
      if strpos(plotvector,'y') eq 0 then plotvector=1 else plotvector = 0

      if (plotvector) then begin
         print,'-1  : automatic selection'
         factors = [1.0, 5.0, 10.0, 20.0, 25.0, $
                    50.0, 75.0, 100.0, 150.0, 200.0]
         nfacs = n_elements(factors)
         for i=0,nfacs-1 do print, tostr(i)+'. '+string(factors(i)*10.0)
         vector_factor = fix(ask('velocity factor','-1'))
      endif else vector_factor = -1

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
      npolar = 0
      if (polar) then begin
         npolestring = mklower(ask('north or south pole (n or s)','s'))
         if (strpos(npolestring,'n') eq 0) then npolar = 1 else npolar = 0
      endif

; MinLat is for polar plots:
      MinLat = 40.0

; showgridyes says whether to plot the grid or not.
      showgridyes = 0

;plotvectoryes says whether to plot vectors or not
      plotvectoryes = plotvector

; number of points to skip when plotting vectors:
      step = 2

; vi_cnt is whether to plot vectors of Vi
      vi_cnt = 0

; vn_cnt is whether to plot vectors of Vn
      vn_cnt = 1

      cursor_cnt = 0

      xrange = [0.0,0.0]

      yrange = [0.0,0.0]

      test = "V!Dn!N(east)"
      for i=0,nvars-1 do begin
         var=strcompress(vars[i],/remove_all)
         tes = strcompress(test,/remove_all)
         result= STRCMP( var, tes,8 )
         if (result eq 1) then vneast_index = i
      endfor

   endif

   c_a_to_r, itime1, t1
   c_a_to_r, itime2, t2
   dt = t1-t2
   dlon = dt * 15.0 / 3600.0
   print, dlon
   dlon_in_model = (data(0,1,0,0)-data(0,0,0,0))/!dtor
   sh = fix(dlon/dlon_in_model)

print, sh

   sh = -1

   data_sub = data
   dbs = shift(reform(data_base(sel           ,*,*,*)),sh,0,0)
   dbe = shift(reform(data_base(vneast_index  ,*,*,*)),sh,0,0)
   dbn = shift(reform(data_base(vneast_index+1,*,*,*)),sh,0,0)

   print, dlon, dlon_in_model, sh

   data_sub(sel,*,*,*) = 100.0*(data(sel,*,*,*)-dbs) / dbs

   data_sub(vneast_index,*,*,*) = $
      data(vneast_index,*,*,*) - dbe
   data_sub(vneast_index+1,*,*,*) = data(vneast_index+1,*,*,*) - dbn

   p = strpos(psfile,'.ps')
   if (p gt -1) then psfile = strmid(psfile,0,p)
   psfile_final = psfile+'_'+chopr('000'+tostr(iFile),4)+'.ps'

   smini_final = smini
   smaxi_final = smaxi

   thermo_plot,cursor_x,cursor_y,strx,stry,step,nvars,sel,nfiles,  $
               cnt1,cnt2,cnt3,yes,no,yeslog,  	  $
               1-yeslog,nalts,nlats,nlons,yeswrite_cnt,$
               polar,npolar,MinLat,showgridyes,	  $
               plotvectoryes,vi_cnt,vn_cnt,vector_factor,	  $
               cursor_cnt,data_sub,alt,lat,lon,	  $
               xrange,yrange,selset,smini_final,smaxi_final,	  $
               filename,vars, psfile_final, 0, 'mid', itime1 ;, $
;                     plotpole = 0


   filename_base = filename
   nvars_base = nVars
   nAlts_base = nAlts
   nLons_base = nLons
   nLats_base = nLats
   Vars_base = Vars
   Data_base = Data
   itime2 = itime1

endfor

end
