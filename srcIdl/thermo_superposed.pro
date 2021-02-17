
filelist = findfile("-t *.save")
if (strlen(filelist(0)) eq 0) then filelist = findfile("-t *.dat")

filelist = ask('filelist to plot',filelist(0))

filelist = findfile(filelist)
nFiles = n_elements(filelist)

for iFile = 0, nFiles-1 do begin

   filename = filelist(iFile)
   print, 'Reading file ',filename

   read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $
                           vars, data, rb, cb, bl_cnt, itime1

   alt = reform(data(2,*,*,*)) / 1000.0
   lat = reform(data(1,*,*,*)) / !dtor
   lon = reform(data(0,*,*,*)) / !dtor

   if (iFile eq 0) then begin

      for i=0,nvars-1 do print, tostr(i)+'. '+vars(i)
      sel = fix(ask('which var to plot','3'))

      psfile = filename+'.ps'
      psfile = ask('ps file name',psfile)

      for i=0,nalts-1 do print, tostr(i)+'. '+string(alt(2,2,i))
      selset = fix(ask('which altitude to plot','0'))

      smini = ask('minimum (0.0 for automatic)','0.0')
      smaxi = ask('maximum (0.0 for automatic)','0.0')

   endif

   


endfor
