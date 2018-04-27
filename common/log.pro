;- routine to open new logfile -
pro newlog, fname,fid
  openw, fid, fname, /get_lun
end
;-------------------------------

;- routine to append to the logfile -
pro appendlog, fname, fid
  openu, fid, fname, /get_lun
end
;------------------------------------

;- routine to close logfile -
pro closelog, fname, fid
  close, fid
  free_lun, fid
  file_chmod,fname,/g_write
end
;----------------------------

;- routine to write the screen (standard output), the logfile or both -
pro printlog, data0,data1,data2,data3,data4,data5,data6,data7,data8,data9, format=format, fid=fid, so=so

  if ~keyword_set(format) then format='(10(A,:))'

  nparams = n_params()
  case nparams of
    0   : message, 'PRINTLOG needs at least one argument'
    1   : data = {data0: data0}
    2   : data = {data0: data0,data1: data1}
    3   : data = {data0: data0,data1: data1,data2: data2}
    4   : data = {data0: data0,data1: data1,data2: data2,data3: data3}
    5   : data = {data0: data0,data1: data1,data2: data2,data3: data3,data4: data4}
    6   : data = {data0: data0,data1: data1,data2: data2,data3: data3,data4: data4,data5: data5}
    7   : data = {data0: data0,data1: data1,data2: data2,data3: data3,data4: data4,data5: data5,$
                  data6: data6}
    8   : data = {data0: data0,data1: data1,data2: data2,data3: data3,data4: data4,data5: data5,$
                  data6: data6,data7: data7}
    9   : data = {data0: data0,data1: data1,data2: data2,data3: data3,data4: data4,data5: data5,$
                  data6: data6,data7: data7,data8: data8}
    10  : data = {data0: data0,data1: data1,data2: data2,data3: data3,data4: data4,data5: data5,$
                  data6: data6,data7: data7,data8: data8,data9: data9}
    else: message, 'PRINTLOG can have maximum 10 arguments'
  endcase

  ; print to standard output
  if keyword_set(so) then print, format=format, data

  ; print to file
  if keyword_set(fid) then printf,fid, format=format, data

end
;----------------------------------------------------