

if (n_elements(dirin) eq 0) then dirin = '/bigdisk1/Data6/Data/CHAMP'
dirin = ask('directory for CHAMP data',dirin)

if (n_elements(year) eq 0) then year = '2003'
year = ask('year',year)

if (n_elements(month) eq 0) then month = '10'
month = string(fix(ask('month',month)), format='(I02)')

if (n_elements(day) eq 0) then day = '29'
day = string(fix(ask('day',day)), format='(I02)')

julian = string(format='(I03)',jday(fix(year),fix(month),fix(day)))

filein = dirin+'/'+year+'/Density_3deg_'+chopr(year,2)+'_'+julian+'.ascii'

fileout = 'champ'+year+month+day+'.dat'

line = ''

openr,1,filein
openw,2,fileout

printf,2,''
printf,2,'Position taken from file ',filein
printf,2,''
printf,2,'#START'

readf,1,line
readf,1,line

t = ''

while not eof(1) do begin

    readf,1,t
    tmp = float(strsplit(t,/extract))
    hour = floor(tmp(2)/3600.)
    min = floor((tmp(2)/3600.-hour)*60)
    sec = floor(((tmp(2)/3600.-hour)*60-min)*60)
    itime = [year,month,day,hour,min,sec]
    c_a_to_r, itime, rtime
    c_r_to_a, itime, rtime
    lon = float(tmp(5))
    if lon lt 0 then lon = 360. + lon
    printf, 2, itime, 0, lon, tmp(4), tmp(6), format = '(7i5,3f8.2)'


endwhile

close,1,2

end

