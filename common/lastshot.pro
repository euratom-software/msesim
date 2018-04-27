function lastshot
;+
; out = lastshot()
;
; Function returns the lastest shotnumber
;
; :Returns:
;   out : out, required, type=long, scalar
;        the latest MAST shotnumber
;-
  mshot='/funsrv1/home/schedule/Datac/Mast/mshot.dat'
  openr, fileID, mshot, /get_lun
  shot = 0l
  readf, fileID, shot
  close,fileID
  free_lun, fileID
  return, shot
end
