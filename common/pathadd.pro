pro pathadd, pth, quiet=quiet
;
; PATHADD, PATH, QUIET=QUIET
;
;  PATH:  string with the path to add.
;  QUIET: if set no info is printed to the screen
;
; Routine adds 'pth' to the idl-path, but only if it's not already there.
; All subdirectories of 'pth' that contain '.pro' and '.sav' files will be added.
;

; read the current path into an array of strings 
p0=strsplit(!path,':',/extr)

; Find all subdirectories of 'pth' that contain '.pro' and '.sav' files
; and put them into an array of strings.
expan=expand_path('+'+pth,/array)
nexp=n_elements(expan)
; Loop through these directories and add them to the path if they aren't
; already a part of the path
for i=0,nexp-1 do begin
  fnd=where(p0 eq expan[i])
  if fnd[0] ne -1 then begin
    if ~keyword_set(quiet) then print,'path '+expan[i]+' already there, not adding'
    continue
  endif
  ; add to the path
  !path=expan[i]+':'+!path
  if ~keyword_set(quiet) then print,'added ',expan[i]
endfor

end

