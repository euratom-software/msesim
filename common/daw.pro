pro daw
;
; deletes all open windows
;

winidx = !D.window
while winidx ne -1 do begin
  wdelete
  winidx = !D.window
endwhile

end