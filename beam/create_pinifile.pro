pro create_pinifile
; creates an ASCII spreadsheet file with the y- and z-coordinates
; of the PINI-holes in the beam coordinate system used by the STARK-code
; (this is the coordinate system with the x-axis on the beam axis and the y-axis
;  in the xy-plane of the machine coordinate system).
; The coordinates are given in [m].
; The data is collected from the 'GRID_MKIII.dat' data-file.

pini=fltarr(2,131)
commentstr=''
get_lun,funit
openr,funit, 'GRID_MKIII.dat'
  readf,funit, commentstr
  readf,funit, a
  readf,funit, yoffset
  readf,funit, zoffset
  readf,funit, d
  readf,funit, e
  readf,funit, f
  readf,funit, commentstr
  readf,funit, commentstr
  readf,funit, commentstr
  readf,funit, pini
free_lun,funit
close,funit

; get the y- and z-coordinates and convert them from mm to m
ypini=(pini[1,*]-yoffset)*1e-3
zpini=(pini[0,*]-zoffset)*1e-3

; mirror the pini-grid to get the lower half
ypini=[[ypini],[ypini]]
zpini=[[zpini],[-zpini]]
pini =fltarr(2,262)
pini = [ypini,zpini]

; write this data to the file 'pini.dat'
get_lun,funit
openw,funit, 'pini.dat'
  printf,funit,'; This ASCII spreadsheet file contains the y- and z-coordinates'
  printf,funit,'; of the PINI-holes in the beam coordinate system used by the STARK-code'
  printf,funit,'; (this is the coordinate system with the x-axis on the beam axis and'
  printf,funit,';  the y-axis in the xy-plane of the machine coordinate system).'
  printf,funit,'; The coordinates are given in [m].'
  printf,funit,'; The data was collected from the "GRID_MKIII.dat" data-file.'
  printf,funit,';'
  printf,funit,';     y [m]        z[m]'
  printf,funit,pini
free_lun,funit
close,funit

end
