
filelist1 = findfile('UA/data.RestartTest.FullRun/log*.dat')
filelist2 = findfile('UA/data.RestartTest.Restarted/log*.dat')

thermo_readlog, filelist1, data1, time1
thermo_readlog, filelist2, data2, time2

; data( 8,*) = dt
; data( 9,*) = min t
; data(10,*) = max t
; data(11,*) = mean t
; data(12,*) = min v
; data(13,*) = max v
; data(14,*) = mean v

diff = data1 - data2

;plot, time1-time1(0), data(

end


