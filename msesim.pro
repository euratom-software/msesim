;; Copyright 2018 Euratom/CCFE/TuE

;; Permission is hereby granted, free of charge, to any person obtaining a copy of this
;; software and associated documentation files (the "Software"), to deal in the Software
;; without restriction, including without limitation the rights to use, copy, modify,
;; merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
;; permit persons to whom the Software is furnished to do so, subject to the following
;; conditions:

;; The above copyright notice and this permission notice shall be included in all copies
;; or substantial portions of the Software.

;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
;; INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
;; PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
;; HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
;; CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
;; OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

pro msesim, run, batchfile=batchfile
; This is the main STARK-code program.
;
; The input parameter 'run' is a string containing the name of one of the subdirectories
; of the 'runs'-directory. This program will then:
;  - if the keyword 'batchfile' is set, read the xml-file indicated by 'batchfile' in that directory
;    else it reads the 'batch.xml' file in that directory
;  - execute the desired commands in the batch xml-file by calling the routines:
;        spectrum_calc
;        spectrum_filter
;        spectrum_merge
;        spectrum_plot
;        spectrum_compare2
;        spectrum_compareN
;  - save the data in the output/data (xdr-files) and data/pictures (eps-files) directories
;

;----------------------------------------------------------------------------------
; Distinguish between Windows and UNIX directory separators
;----------------------------------------------------------------------------------
if strcmp( !version.os,'Win32',/fold_case) then sep='\' else sep='/'
;
; Get the source code path and change the working directory to that path
proname = 'msesim'
help, /source, name=proname,output=output
outputstr=''
for i=0,n_elements(output)-1 do outputstr+=output[i]
first = strpos(outputstr,  strupcase(proname) )+strlen(proname)+1
last  = strpos(outputstr,'Compiled Functions:',/reverse_search)
length= last-first
sourcefile = strtrim(strmid(outputstr,first,length),2)
sourcedir  = strmid(sourcefile,0,strpos(sourcefile, sep, /reverse_search))
cd, sourcedir


; set standard 'batch.xml' if no custom one is specified by the batchfile-keyword
if ~keyword_set(batchfile) then batchfile='batch.xml'
batchfile = sourcedir+sep+'runs'+sep+run+sep+batchfile
; read in the command-structure array from the batch xml-file
cmdarr = read_batch(batchfile)
ncmd   = n_elements(cmdarr)

; get some paths:
settingdir = 'runs' + sep + run + sep + 'settings' + sep
print,'settingdir',settingdir
datadir    = 'runs' + sep + run + sep + 'output' + sep + 'data' + sep
pictdir    = 'runs' + sep + run + sep + 'output' + sep + 'pictures' +sep

print, ''
; loop through the commands and execute them
for i=0,ncmd-1 do begin

  case cmdarr[i].name of
    'calc'    : begin
                  ; get the command parameters
                  settingfile = settingdir + cmdarr[i].setting
                  outputfile  = datadir    + cmdarr[i].output
                  print,'settingfile',settingfile
                  print,'outputfile',outputfile
                  ; run the command
                  spectrum_calc, settingfile, outputfile
                end
    'filter'  : begin
                  ; get the command parameters
                  settingfile = settingdir + cmdarr[i].setting
                  inputfile   = datadir    + *cmdarr[i].input	; remember: cmdarr[i].input is a pointer!
                  outputfile  = datadir    + cmdarr[i].output
                  ; run the command
                  spectrum_filter, settingfile, inputfile, outputfile
                end
    'merge'   : begin
                  ; get the command parameters
                  inputfiles  = *cmdarr[i].input		; remember: cmdarr[i].input is a pointer!
                  for k=0,n_elements(inputfiles)-1 do begin
                    inputfiles[k] = datadir +inputfiles[k]
                  endfor
                  outputfile  = datadir + cmdarr[i].output
                  ; run the command
                  spectrum_merge, inputfiles, outputfile
                end
    'plot'    : begin
                  ; get the command parameters
                  settingfile = settingdir + cmdarr[i].setting
                  inputfile   = datadir    + *cmdarr[i].input	; remember: cmdarr[i].input is a pointer!
                  prefix      = cmdarr[i].prefix
                  ; run the command
                  spectrum_plot, settingfile, inputfile, prefix, pictdir
                end
    'compare2': begin
                  ; get the command parameters
                  settingfile = settingdir + cmdarr[i].setting
                  inputfiles  = *cmdarr[i].input		; remember: cmdarr[i].input is a pointer!
                  for k=0,n_elements(inputfiles)-1 do begin
                    inputfiles[k] = datadir +inputfiles[k]
                  endfor
                  labels      = cmdarr[i].labels
                  prefix      = cmdarr[i].prefix
                  ; run the command
                  spectrum_compare2, settingfile, inputfiles, labels, prefix, pictdir
                end
    'comparen': begin
                  ; get the command parameters
                  settingfile = settingdir + cmdarr[i].setting
                  paramfile   = settingdir + cmdarr[i].param
                  inputfiles  = *cmdarr[i].input		; remember: cmdarr[i].input is a pointer!
                  for k=0,n_elements(inputfiles)-1 do begin
                    inputfiles[k] = datadir +inputfiles[k]
                  endfor
                  prefix      = cmdarr[i].prefix
                  ; run the command
                  ;spectrum_compareN, settingfile, paramfile, inputfiles, prefix, pictdir
                end
  endcase
  ; free the pointer of the input field!
  ptr_free, cmdarr[i].input

  ; print a separator
  print, ''
  print, '-----------------------------------------------------------------------------'
  print, ''
endfor

end
