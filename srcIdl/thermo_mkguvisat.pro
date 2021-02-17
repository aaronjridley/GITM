
filelist = findfile("*.L1B")
nfiles = n_elements(filelist)
for i=0,nfiles-1 do print, tostr(i),". ",filelist(i)
if (n_elements(ifile) eq 0) then ifile = 0
ifile = fix(ask('file to create GUVI sat file from',tostr(ifile)))

ofile = ask('4-letter output code for satellite file','guvi')

file = filelist(iFile)

read_guvi, file, time, data, latitude, longitude

nPts  = n_elements(time)
nSwath = n_elements(latitude(0,*,0))
nSkip  = nSwath/20

c_r_to_a, itime, time(0)
syear  = tostr(itime(0))
smonth = chopr('0'+tostr(itime(1)),2)
sday   = chopr('0'+tostr(itime(2)),2)
shour  = chopr('0'+tostr(itime(3)),2)

ofile = ofile+'.'+syear+smonth+sday+shour+'.in'

openw,1,ofile

printf,1,''
printf,1,'#START'

for i=0,nPts-1 do begin
  c_r_to_a, itime, time(i)
  for j=0,nSwath-1, nSkip do begin
      printf,1,format='(7i5,3f8.2)',itime, 0, $
        longitude(0,j,i), latitude(0,j,i), 0.0
  endfor
endfor

close,1

end
