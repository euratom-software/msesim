pro n2,nel_space,tel_space,n2
@com_bloc

; File provided by Manfred von Hellermann (or ADAS) containing the fraction of beam 
; neutrals that are excited (n=2) as function of beam energy and electron density
n2file = 'mvhcxnpoptest.out'

; show what this procedure does:
print, 'Calculating the fraction of the beam density that is in the n=2 excited state:'
print, '=> interpolation of ADAS n=2-fraction table in "'+n2file+'"'
print, '   to the relevant beam energies and electron densities'

; Read the data from the file
n2file = 'mvhcxnpoptest.out'
data=fltarr(7,12)	; file contains a 7x12 matrix of data points
get_lun,unit1
openr,unit1, n2file
readf,unit1,data
close,unit1
free_lun,unit1
mv_e=data[0,1:11]	; the first column contains a range of beam energies [keV]
mv_n=data[1:6,0]	; the first row contains a range of electron densities [m^{-3}]
s=data[1:6,1:11]	; the rest contains the n=2 fraction

; Now interpolate this data to the beam energy and electron density of the code 
n2 = fltarr(xdata[2], ydata[2], zdata[2], n_elements(tcalc),Enumber)
for tidx=0,n_elements(tcalc)-1 do begin
  for gridi=0,xdata(2)-1 do begin
    for gridj=0,ydata(2)-1 do begin  
      for gridk=0,zdata(2)-1 do begin 
        for fr=0,Enumber-1 do begin $
          n2[gridi,gridj,gridk,tidx,fr] = $
          inter_2d(s, mv_n, mv_e, nel_space[gridi,gridj,gridk,tidx], ebeam[fr]/abeam*1e-3)
        endfor
      endfor
    endfor
  endfor
endfor



;stop

end