function read_setting, file
;
; reads the xml-setting file and returns a nested structure.
; The 1st level fields are: equi, beam, integration, plot, ...
; Each of these fields is a structure that contains the parameters for
; the equilibrium, the beam, the integration, the plotting, ...


  ; get the full filename (xml-routines don't work with things like '~mdebock') and check if it exists
  file = file_search(file,/fully_qualify_path)
  print,file
  file = file[0]
  if file eq '' then begin
    print, format='("ERROR: file <",A,"> not found!")',file
    return, -1
  endif

  ; load the xml-file
  oDocument = obj_new('IDLffXMLDOMDocument', filename=file)

  ; get the top node. REMARK: this has to be called 'stark_settings'!
	nodelist0 = oDocument->GetElementsByTagName('stark_settings')
	topnode   = nodelist0->Item(0)
	
	; get the 1st level nodes 
	nodelist1 = topnode->GetChildNodes()
	n_nodes1  = nodelist1->GetLength()
	; and loop through them
	k=0
	for i=0,n_nodes1-1 do begin
		; get the node name
		node1 = nodelist1->Item(i)
		name1 = node1->GetNodeName()
		; if the node name is '#text', we ignore the node (i.e. text not 
		; embedded into the <> </> delimiters, most likely white space)
		if strcmp(name1,'#text') then continue
	
		; if the node name differs from '#text', we make up a list of its children
		nodelist2 = node1->GetChildNodes()
		n_nodes2  = nodelist2->GetLength()
		l=0
		for j=0,n_nodes2-1 do begin
			; get the node name
			node2 = nodelist2->Item(j)
			name2 = node2->GetNodeName()
	
			; again ignore the nodes named '#text'
			if strcmp(name2,'#text') then continue
	
	
			; get the node's attributes
			attributes = node2->GetAttributes()
			; the get attributes 'value' and 'type' 
			; ('description' and 'unit' are ignored. they are only present to
			;  indicate what value the user should put into the setting file)
			valueNode = attributes->getNamedItem('value')
			typeNode  = attributes->getNamedItem('type')
			; get string values of each these 2 attributes
			valueStr  = valueNode->getNodeValue()
			typeStr   = typeNode->getNodeValue()
	
			; convert 'valueStr' based upon 'typeStr',
			; we use strsplit for this, which returns a 'comma'-separated list into an array.
			; if there's only one element in the array, the array is converted in a scalar.
			; if the type is unknown, then it is considered to be a 'float'
			typeStr = strlowcase(typeStr)	; avoid errors due to lower and uppercase differences
			case typestr of
				 'string':  begin
											value = strtrim(strsplit(valueStr,',',/extract),2)
											if n_elements(value) eq 1 then value=value[0]
										end
				 'float':   begin
											value = float(strsplit(valueStr,',',/extract))
											if n_elements(value) eq 1 then value=value[0]
										end
				 'double':  begin
											value = double(strsplit(valueStr,',',/extract))
											if n_elements(value) eq 1 then value=value[0]
										end
				 'integer': begin
											value = float(strsplit(valueStr,',',/extract))
											if n_elements(value) eq 1 then value=value[0]
											value = fix(value)
										end
				 'long':    begin
											value = float(strsplit(valueStr,',',/extract))
											if n_elements(value) eq 1 then value=value[0]
											value = long(value)
										end
				 'byte':    begin
											value = float(strsplit(valueStr,',',/extract))
											if n_elements(value) eq 1 then value=value[0]
											value = byte(value)
										end
				 else:      begin
											value = float(strsplit(valueStr,',',/extract))
											if n_elements(value) eq 1 then value=value[0]
										end
			endcase
	
			; create the 1st level structure
			if l eq 0 then begin
				 s1 = create_struct(name2,value)
			endif else begin
				 s1 = create_struct(name2,value,s1)
			endelse
			l++
		endfor
	
		; create the setting structure
		if k eq 0 then begin
			 setting = create_struct(name1,s1)
		endif else begin
			 setting = create_struct(name1,s1,setting)
		endelse
		k++
	endfor
	
	; clean up all xml-objects (i.e. free memory)
	obj_destroy, oDocument
	
	; return the setting structure
	return, setting

end
