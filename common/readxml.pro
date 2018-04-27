;+
; a = readxml( filename [, name=name, /string] )
;
; reads an xml-file and returns a nested structure.
;
; e.g. an xml-file with following structure:
; <?xml version="1.0" encoding="UTF-8" standalone="no" ?>
; <structurename>
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
; </structurename>
;
; Then the readxml-function will return following structure:
; {nodename0: 3                                              ,$   ; integer
;  nodename1: { nodename10: ['abc','def']                    ,$   ; string array
;               nodename11: [1e4, 0.01, 0.003] }             ,$   ; float array
;  nodename2: { nodename20: 30000000                         ,$   ; long
;               nodename21: { nodename210: 1e100             ,$   ; double
;                             nodename211: [0, 1, 2, 3] } } }     ; byte array
;
; the keyword "string" is set, the "filename" argument is considered to be an
; xml-structured string that will be translated in an IDL structure
;
; :Params:
;   filename: in, required, type=string, scalar
;            the name of the xml-file or the xml-string to be read
; :Keywords:
;   name  : out, optional, type=string, scalar
;          returns the name of the structure (top node xml-file/-string)
;   string: in, optional, type=byte, scalar
;          if set the "filename" argument is considered to be an xml-structured
;          string that will be translated in an IDL structure
; :Returns:
;   out   : out, required, type=structure, scalar
;          IDL output structure containing the xml-data
;-



;-------------------------------------
; recursive routine that gets all the
; nodes and puts them in a structure
;-------------------------------------
function getnodes, struct, parent

; get the parent name

parentname = parent->GetNodeName()

; get the list with children
children = parent->GetChildNodes()

n_child  = children->GetLength()

; if this particular node has no children: 
;  - get its type and value from the attributes
;  - write a value field to the structure (struct)
if n_child eq 0 then begin
  ; get the node's attributes
  attributes = parent->GetAttributes()

  ; the get attributes 'value' and 'type' 
  valueNode = attributes->getNamedItem('value')
  typeNode  = attributes->getNamedItem('type')
  ; get string values of each these 2 attributes
  valueStr = valueNode->getNodeValue()
  typeStr  = typeNode->getNodeValue()

  ; convert 'valueStr' based upon 'typeStr'. if the type is unknown, then it is considered to be a 'float'
  typeStr  = strlowcase(typeStr) ; avoid errors due to lower and uppercase differences
  case typestr of
     'string'  : data = stringscan(valueStr)
     'float'   : data = float(numscan(valueStr))
     'single'  : data = float(numscan(valueStr))
     'double'  : data = double(numscan(valueStr))
     'integer' : data = fix(numscan(valueStr))
     'int16'   : data = fix(numscan(valueStr))
     'long'    : data = long(numscan(valueStr))
     'int32'   : data = long(numscan(valueStr))
     'int64'   : data = long(numscan(valueStr))
     'byte'    : data = byte(numscan(valueStr))
     'uint8'   : data = byte(numscan(valueStr))
     'complex' : data = complex(numscan(valueStr))
     'scomplex': data = complex(numscan(valueStr))
     'dcomplex': data = numscan(valueStr)
     else      : data = float(numscan(valueStr))
  endcase
  ; make scalars of arrays with size 1
  if n_elements(data) eq 1 then data = data[0]

; if this particular node does have children:
;  - make a substructure
endif else begin
  ; loop recursively through the child nodes (and ignore 'empty' nodes)
  ntag = 0
  for i=0,n_child-1 do begin
    ; get the node name
    child     = children->Item(i)
    childname = child->GetNodeName()
    ; if the node name is '#text', then it isn't a 'real' node
    ; (i.e. text not embedded into the <> </> delimiters, most likely white space)
    if strcmp(childname,'#text') then continue
    ; create a sub structure with empty tags
    if ntag eq 0 then substruct = create_struct(childname,'') else substruct = create_struct(substruct, childname,'')
    ; fill that substructure with the data of the children
    substruct = getnodes(substruct, child)
    ntag = ntag+1
  endfor
  data = substruct

endelse

;  if size(data,/type) eq 8 then  helpstr, data else help, data

  ; save and return the data (either value or structure)
  if size(struct,/type) eq 8 then begin
    struct = replace_tag(struct,parentname,data)
  endif else begin
    struct = data
  endelse  
  return, struct

end


;------------------
; the main function
;------------------
function readxml, file, name=name, string=string
print,file
; open document
if ~keyword_set(string) then begin
  ; get the full filename (xml-routines don't work with things like '~mdebock') and check if it exists
  filesearch = file_search(file[0]+file[1])
  print,'filesearch', filesearch

  file = file[1]
  print,'filename',file

  if file eq '' then begin
    print, format='("ERROR: file <",A,"> not found!")',file
    return, -1
  endif
  ; load the xml-file
  oDocument = obj_new('IDLffXMLDOMDocument', filename=file)
endif else begin
  ; load the xml-string
  oDocument = obj_new('IDLffXMLDOMDocument', string=file)
endelse

topnode = oDocument->GetFirstChild()
print,'topnode', topnode
structName = topnode->GetNodeName()

; start with an empty structure
struct = 0
; recursively get all the nodes under the top node
struct = getnodes(struct, topnode)

; get rid of the object pointers
obj_destroy, oDocument
heap_gc

return, struct

end
