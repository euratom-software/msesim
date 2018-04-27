;+
; WRITEXML, structure, filename [, name=name, /string]
;
; writes a nested structure to an xml-file and
;
; Consider following structure:
; s = {nodename0: 3                                              ,$   ; integer
;      nodename1: { nodename10: ['abc','def']                    ,$   ; string array
;                   nodename11: [1e4, 0.01, 0.003] }             ,$   ; float array
;      nodename2: { nodename20: 30000000                         ,$   ; long
;                   nodename21: { nodename210: 1e100             ,$   ; double
;                                 nodename211: [0, 1, 2, 3] } } }     ; byte array
;
; Then 'writexml, s, 'testwrite.xml', name='test' will produce following xml-file:
; <?xml version="1.0" encoding="UTF-8" standalone="no" ?>
; <test>
;   <nodename0
;       type    = 'integer'
;       value   = '3'
;   />
;   <nodename1>
;     <nodename10
;       type    = 'string'
;       value   = 'abc, def'
;     />
;     <nodename11
;       type    = 'float'
;       value   = '1e4, 0.01, 0.003'
;     />
;   </nodename1>
;   <nodename2>
;     <nodename20
;       type    = 'long'
;       value   = '30000000'
;     />
;     <nodename21>
;       <nodename210
;         type    = 'double'
;         value   = '1e100'
;       />
;       <nodename211
;         type    = 'byte'
;         value   = '0,1,2,3'
;       />
;     </nodename21>
;   </nodename2>
; </test>
;
; the keyword "string" is set, no output file is written.
; Instead the "filename"-argument now is a string variable to which the xml-structure is written.
;
; :Params:
;   structure: in, required, type=structure, scalar
;             IDL structure that needs to be saved into the xml-file
;   filename : in, required, type=string, scalar
;             the name of the xml-file or the xml-string to be read
; :Keywords:
;   name     : in, optional, type=string, scalar
;             name of the structure (top node xml-file/-string)
;   string   : in, optional, type=byte, scalar
;             if set no output file is written. instead the "filename" argument is
;             now a string variable to which the xml-structure is written.
;-

;-------------------------------------
; recursive routine that goes through
; the structure and writes the nodes
;-------------------------------------
pro writenodes, doc, node, struct

ntags    = n_tags(struct)

; if there are no tags, then struct isn't a structure but a variable.
; so then we write a node with the attributes 'type' and 'value'
if ntags eq 0 then begin
  ; get the type (if the type is unknown, then it is considered to be a 'float')
  value = struct
  type  = size(value, /type)
  case type of
    1   : typeStr='byte'
    2   : typeStr='integer'
    3   : typeStr='long'
    4   : typeStr='float'
    5   : typeStr='double'
    6   : typeStr='float'   ; complex float  => make 2 xml-entries!
    7   : typeStr='string'
    9   : typeStr='double'  ; complex double => make 2 xml-entries!
    else: typeStr='float'
  endcase

  ; write the type attribute
  typeNode = doc->CreateAttribute('type')
  typeNode->SetValue, typeStr
  void     = node->SetAttributeNode(typeNode)

  ; write the value attribute
  valueNode = doc->CreateAttribute('value')
  ; this type of xml-file can only handle 1D arrays:
  value = reform(value,n_elements(value))
  ; If value is complex: just save the real part at show a warning
  if (type eq 6) || (type eq 9) then begin
    value=real_part(value)
    print, 'WARNING: Complex numbers detected in the structure! Only the real part will be saved to the XML file!'
  endif
  ; if the data type is 'byte', then convert to integer first, because otherwise the 'string' function
  ; will generate the ASCII character associated with the byte value, rather than the actual number
  if type eq 1 then value=fix(value)
  ; Loop through value and make the valueStr
  if n_elements(value) eq 1 then begin
    valueStr=string(value[0])
  endif else begin
    valueStr=''
    for i=0,n_elements(value)-2 do begin
      valueStr=valueStr+string(value[i])+","
    endfor
    valueStr=valueStr+string(value[n_elements(value)-1])
  endelse
  valueNode->SetValue, valueStr
  void     = node->SetAttributeNode(valueNode)

; if there are tags, then recursively create child nodes
endif else begin
  tagnames = tag_names(struct)
  for i=0,ntags-1 do begin
    name  = strlowcase(tagnames[i])
    childnode = doc->CreateElement(name)
    void      = node->AppendChild(childnode)
    writenodes, doc, childnode, struct.(i)
  endfor
endelse

end


;------------------
; the main routine
;------------------
pro writexml, struct, file, name=name, string=string

if n_params() eq 0 then begin
  print, 'USAGE: writexml, structure, filename (or string variable) [,name=name, /string]'
  return
endif


; open document
oDocument = OBJ_NEW('IDLffXMLDOMDocument')

; create the top node
if ~keyword_set(name) || ~is_string(name) then name='structure'
topnode = oDocument->CreateElement(name)
void    = oDocument->AppendChild(topnode)

; recursively go through the structure and make nodes
writenodes, oDocument, topnode, struct

if ~keyword_set(string) then begin
  ; save the xml-file
  oDocument->Save, filename=file,/pretty_print
endif else begin
  oDocument->Save, string=file,/pretty_print
endelse

; get rid of the object pointers
obj_destroy, oDocument
heap_gc

end
