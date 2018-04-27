PRO wr2matlab,x,dataname,filename1, append=append

;	write numerical array x to directly executable matlabfile

;	.compile wr2matlab

filename=strcompress(string(filename1)+'.m',/remove_all)

if file_test(filename) then begin
  if file_test(filename, /write) then begin
    if keyword_set(append) then begin
      openw, lun, filename,/append, /get_lun
    endif else begin
      prompt = "'"+filename+"' already exists! Overwrite or cancel? (o/c)"
      answer = ''
      read, prompt=prompt, answer
      if strcmp(answer,'o',/fold_case) then begin
        openw, lun, filename, /get_lun
      endif else begin
        return
      endelse
    endelse
  endif else begin
    print, "'"+filename+"' already exists and is read-only!"
    return
  endelse
endif else begin
  openw, lun, filename, /get_lun
endelse

dim=size(x,/dimensions)

ndim=n_elements(dim)

if ndim eq 1 then begin
if dim eq 0 then begin
dim=dim+1
endif
endif

printf,lun,STRLOWCASE(dataname)+'=['

; 1D
if ndim eq 1 then begin
;for j=0,dim[0]-1 do begin 
form='('+string(dim[0])+'(E,"  ")/)'
printf,lun,x,FORMAT=form
;endfor;j
endif

; 2D arrays
if ndim eq 2 then begin
for j=0,dim[0]-1 do begin 
form='('+string(dim[1])+'(E,"  ")/)'
printf,lun,x[j,*],FORMAT=form
endfor;j
endif

; 3D arrays
if ndim eq 3 then begin
for j=0,dim[0]-1 do begin 
for k=0,dim[1]-1 do begin 
form='('+string(dim[2])+'(E,"  ")/)'
printf,lun,x[j,k,*],FORMAT=form
endfor ;k
endfor;j
endif

; 4D arrays
if ndim eq 4 then begin
for j=0,dim[0]-1 do begin 
for k=0,dim[1]-1 do begin 
for t=0,dim[2]-1 do begin 
form='('+string(dim[3])+'(E,"  ")/)'
printf,lun,x[j,k,t,*],FORMAT=form
endfor ;t
endfor ;k
endfor;j
endif


printf,lun,'];'

if ndim eq 3 then begin
printf,lun,STRLOWCASE(dataname)+'=reshape('+STRLOWCASE(dataname)+','+string(dim[1])+','+string(dim[0])+','+string(dim[2])+');'
printf,lun,STRLOWCASE(dataname)+'=permute('+STRLOWCASE(dataname)+',[2,1,3]);'
endif

if ndim eq 4 then begin
printf,lun,STRLOWCASE(dataname)+'=reshape('+STRLOWCASE(dataname)+','+string(dim[2])+','+string(dim[1])+','+string(dim[0])+','+string(dim[3])+');'
printf,lun,STRLOWCASE(dataname)+'=permute('+STRLOWCASE(dataname)+',[3,2,1,4]);'
endif

close,lun
free_lun, lun

end
