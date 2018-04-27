pro error_message, name

  if size(name, /type) eq 7 then begin
    print, name, format='(a)'
    print
  endif

  help, /traceback, output=output

  ord=where(strcmp(output, '%', 1))
  print, output(ord(2):*), format='(a)'
  retall

end
