; Generator which returns the numbers in a file
;
; Useage:
;   First call with the file unit number:
;     status = next_double(fid)
;
;   Then repeatedly call
;     d = next_double()

FUNCTION next_double, fid=fid
  COMMON parse_com, fp, nleft, data
  
  IF KEYWORD_SET(fid) THEN BEGIN
    ; New file - reset
    data = 0
    nleft = 0
    fp = fid
    RETURN, 0
  ENDIF
  
  IF nleft EQ 0 THEN BEGIN
    ; None left, so read in more
    line = ' '
    READF, fp, line
    data = regex_extract(line, '[+-]?([0-9])+($| |(\.[0-9]+([eE][+-]?[0-9]*)?))')
    nleft = N_ELEMENTS(data)
  ENDIF

  num = DOUBLE(0.0)
  
  READS, data[0], num
  
  nleft = nleft - 1
  IF nleft GT 0 THEN data = data[1:*]
  RETURN, num
END
