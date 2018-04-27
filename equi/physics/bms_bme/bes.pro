pro bes,nel_space,tel_space,ratebes
@com_bloc
;restore,f='7107sw1x1bes.xdr',/verbose
;Enumber=3
;restore, f='bms_h1.dat', /verbose
;* **********************************************
ratebes=fltarr(n_elements(nel_space[*,0,0,0]),n_elements(nel_space[0,*,0,0]), $
               n_elements(nel_space[0,0,*,0]),n_elements(nel_space[0,0,0,*]),Enumber)
;***********************************************

for nn=0,n_elements(nel_space[0,0,0,*])-1 do begin
  ebeam = rebin(beamvoltage[nn]/(findgen(Enumber)+1), Enumber, n_elements(nel_space[*,0,0,0]), /sample)
  for zz=0,n_elements(nel_space[0,*,0,0])-1 do begin
    for ff=0,n_elements(nel_space[0,0,*,0])-1 do  begin
      for gg=0,Enumber-1 do begin
        ratebes[*,zz,ff,nn,gg] = $
                   inter_3d(sigmavbes, energy, density, temperature, $
                   ebeam[gg,*]/abeam, nel_space[*,zz,ff,nn]/1e6, tel_space[*,zz,ff,nn])$
                   *nel_space[*,zz,ff,nn]/1e6
      endfor
    endfor
  endfor
endfor


end
