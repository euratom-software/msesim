function read_batch, batchfile
;
; reads the 'batchfile' xml-file in the 'run' directory and returns the command-structure array.
;
; The command-structure array is an array, with as elements a structure.
; This structure has following fields
;   - struct.name    : string containing the name of the command (i.e. 'calc', 'filter', 'merge', 'plot', 'compare2', 'compareN')
;   - struct.setting : string containing the name of the setting file (if any)
;   - struct.input   : pointer to a string or string array containing the name(s) of the input file(s)
;                      (because we don't know how many input files there are this has to be a pointer)
;   - struct.output  : string containing the name of the output file (if any)
;   - struct.prefix  : string containing the prefix for figure-titles (if any)
;   - struct.labels  : string-array (size 2) containing the 2 labels used when 2 spectra are compared (if any)
;   - struct.param   : string containing the parameter file used when N spectra are compared (if any)
;
;
; REMARK: pointers are used. This means one has to be careful to free pointers that are no longer used! 
;         (otherwise you end up spoiling all your memory!)
; v1.0, mdebock 24/07/2007
;


; open document
;--------------

oDocument = obj_new('IDLffXMLDOMDocument', filename=batchfile)

; get the top node. REMARK: this has to be called 'stark_run'!
;-------------------------------------------------------------
nodelist0 = oDocument->GetElementsByTagName('stark_run')
topnode   = nodelist0->Item(0)

; get the 1st level nodes
;------------------------
nodelist = topnode->GetChildNodes()
n_nodes  = nodelist->GetLength()

; and loop through them
;----------------------
cmdarr = create_struct('name','','setting','','input',ptr_new(''),$		; dummy first element of
                       'output','','prefix','','labels',['',''],'param','')	; the command-structure array

for i=0,n_nodes-1 do begin
  ; get the node name
  node = nodelist->Item(i)
  name = node->GetNodeName()

  ; There are only 6 node names accepted:  'calc', 'filter', 'merge', 'plot', 'compare2' and 'compareN'
  name = strlowcase(name)
  skip = 0
  case name of
    'calc'    : begin
                  ; get the node's attributes
                  attributes = node->GetAttributes()
                  ; the get attributes 'settingfile' and 'outputfile' 
                  settingAtt = attributes->getNamedItem('settingfile')
                  outputAtt  = attributes->getNamedItem('outputfile')
                  ; get string values of each these 2 attributes
                  setting = settingAtt->getNodeValue()
                  setting = strtrim(setting,2)
                  output  = outputAtt->getNodeValue()
                  output  = strtrim(output,2)
                  ; the other fields in the command-structure are empty for the 'calc'-command
                  input   = ptr_new('')
                  prefix  = ''
                  labels  = ['','']
                  param   = ''
                end
    'filter'  : begin
                  ; get the node's attributes
                  attributes = node->GetAttributes()
                  ; the get attributes 'settingfile', 'inputfile' and 'outputfile' 
                  settingAtt = attributes->getNamedItem('settingfile')
                  inputAtt   = attributes->getNamedItem('inputfile')
                  outputAtt  = attributes->getNamedItem('outputfile')
                  ; get string values of each these 3 attributes
                  setting = settingAtt->getNodeValue()
                  setting = strtrim(setting,2)
                  output  = outputAtt->getNodeValue()
                  output  = strtrim(output,2)
                  inputstr= inputAtt->getNodeValue()
                  inputstr= strtrim(inputstr,2)
                  input   = ptr_new(inputstr)
                  ; the other fields in the command-structure are empty for the 'filter'-command
                  prefix  = ''
                  labels  = ['','']
                  param   = ''
                end
    'merge'   : begin
                  ; get the node's attributes
                  attributes = node->GetAttributes()
                  ; the get attributes 'inputfiles' and 'outputfile' 
                  inputAtt   = attributes->getNamedItem('inputfiles')
                  outputAtt  = attributes->getNamedItem('outputfile')
                  ; get string values of each these 2 attributes
                  output  = outputAtt->getNodeValue()
                  output  = strtrim(output,2)
                  inputstr= inputAtt->getNodeValue()
                  inputstr= strtrim(strsplit(inputstr,',',/extract),2)
                  input   = ptr_new(inputstr)
                  ; the other fields in the command-structure are empty for the 'filter'-command
                  setting = ''
                  prefix  = ''
                  labels  = ['','']
                  param   = ''
                end
    'plot'    : begin
                  ; get the node's attributes
                  attributes = node->GetAttributes()
                  ; the get attributes 'settingfile','inputfile' and 'prefix' 
                  settingAtt = attributes->getNamedItem('settingfile')
                  inputAtt   = attributes->getNamedItem('inputfile')
                  prefixAtt  = attributes->getNamedItem('prefix')
                  ; get string values of each these 2 attributes
                  setting = settingAtt->getNodeValue()
                  setting = strtrim(setting,2)
                  inputstr= inputAtt->getNodeValue()
                  inputstr= strtrim(inputstr,2)
                  input   = ptr_new(inputstr)
                  prefix  = prefixAtt->getNodeValue()
                  prefix  = strtrim(prefix,2)
                  ; the other fields in the command-structure are empty for the 'filter'-command
                  output  = ''
                  labels  = ['','']
                  param   = ''
                end
    'compare2': begin
                  ; get the node's attributes
                  attributes = node->GetAttributes()
                  ; the get attributes 'settingfile','inputfiles', 'labels' and 'prefix' 
                  settingAtt = attributes->getNamedItem('settingfile')
                  inputAtt   = attributes->getNamedItem('inputfiles')
                  labelsAtt  = attributes->getNamedItem('labels')
                  prefixAtt  = attributes->getNamedItem('prefix')
                  ; get string values of each these 2 attributes
                  setting = settingAtt->getNodeValue()
                  setting = strtrim(setting,2)
                  inputstr= inputAtt->getNodeValue()
                  inputstr= strtrim(strsplit(inputstr,',',/extract),2)
                  input   = ptr_new(inputstr)
                  labels  = labelsAtt->getNodeValue()
                  labels  = strtrim(strsplit(labels,',',/extract),2)
                  prefix  = prefixAtt->getNodeValue()
                  prefix  = strtrim(prefix,2)
                  ; the other fields in the command-structure are empty for the 'filter'-command
                  output  = ''
                  param   = ''
                end
    'comparen': begin
                  ; get the node's attributes
                  attributes = node->GetAttributes()
                  ; the get attributes 'settingfile','paramfile', 'inputfiles' and 'prefix' 
                  settingAtt = attributes->getNamedItem('settingfile')
                  inputAtt   = attributes->getNamedItem('inputfiles')
                  paramAtt  = attributes->getNamedItem('paramfile')
                  prefixAtt  = attributes->getNamedItem('prefix')
                  ; get string values of each these 2 attributes
                  setting = settingAtt->getNodeValue()
                  setting = strtrim(setting,2)
                  inputstr= inputAtt->getNodeValue()
                  inputstr= strtrim(strsplit(inputstr,',',/extract),2)
                  input   = ptr_new(inputstr)
                  param   = paramAtt->getNodeValue()
                  param   = strtrim(param,2)
                  prefix  = prefixAtt->getNodeValue()
                  prefix  = strtrim(prefix,2)
                  ; the other fields in the command-structure are empty for the 'filter'-command
                  output  = ''
                  labels  = ['','']
                end
    else:       skip=1
  endcase
  ; if name is not 'calc', 'filter', etc. then immediately skip to the next node
  if skip then continue

  ; create the command-structure
  cmd = create_struct('name',name,'setting',setting,'input',input,$
                      'output',output,'prefix',prefix,'labels',labels,'param',param)
  ; add it to the command-structure array
  cmdarr = [cmdarr, cmd]
endfor

; remove the first 'dummy' element of the command-structure array
; (DON'T forget to free the pointer first!)
ptr_free, cmdarr[0].input
cmdarr = cmdarr[1:*]

; clean up all xml-objects (i.e. free memory)
obj_destroy, oDocument

; return the setting structure
return, cmdarr

end
