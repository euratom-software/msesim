;+
; cp_shaded, data, x, y,   xrange=xrange, yrange=yrange, zrange=zrange, isotropic=isotropic,
;                          xtitle=xtitle, ytitle=ytitle, ztitle=ztitle, title=title,
;                          color=color, ccolor=ccolor, ctable=ctable, invert=invert,
;                          autocontour=autocontour, cvalues=cvalues, clabel=clabel, downhill=downhill,
;                          showscale=showscale, scalewidth=scalewidth, scalelabel=scalelabel, scalepos=scalepos,
;                          position=position, noerase=noerase, charsize=charsize, error=error
;
;   Routine: cp_shaded
;   Version: 1.22
;   Author : M. De Bock
;   Date   : 12.03.10
; based on:
;   Routine: cp_shade
;   Version: 1.20
;   Author : R.Martin
;   Date   : 8.12.03
;
; Compond plot to display data as a contour and/or shaded 2D plot. Optional scale bar.
;
; Overlays a series of plots:
; (1) If requested a scalebar to the right of the plot.
; (2) A shaded surface plot, viewed from above, to produce a 2D colour
;     graphic.
; (3) A contour plot, with label option selected.
; (4) A plot command to add axis and plot titles
;
; Calling sequence.
;
; :Params:
;   data       : in, required, type=numeric, 2D array
;               The data to be plotted
;   x          : in, optional, type=numeric, 1D array
;               The (monotonically increasing) x-vector for the data. Should have the
;               same number of elements as data[*,0]
;   y          : in, optional, type=numeric, 1D array
;               The (monotonically increasing) y-vector for the data. Should have the
;               same number of elements as data[0,*]
; :Keywords:
;   xrange     : in, optional, type=numeric, 1D array
;               x-range of the plot
;   yrange     : in, optional, type=numeric, 1D array
;               y-range of the plot
;   zrange     : in, optional, type=numeric, 1D array
;               z-range of the plot
;   isotropic  : in, optional, type=byte, scalar
;               If set the x-/y-axes will be isotropic
;   xtitle     : in, optional, type=string, scalar
;               xtitle for the plot
;   ytitle     : in, optional, type=string, scalar
;               ytitle for the plot
;   ztitle     : in, optional, type=string, scalar
;               ztitle for the plot (only shown when the /showscale keyword is used)
;   ytitle     : in, optional, type=string, scalar
;               Overall title for the plot
;   color      : in, optional, type=long, scalar
;               Color for the plot axis
;   ccolor     : in, optional, type=long, scalar
;               Color for the contour lines
;   ctable     : in, optional, type=integer, scalar
;               Color table for the shaded plot. If set to -1 then no shaded plot will be shown (just contours are shown)
;   invert     : in, optional, type=byte, scalar
;               If set the color table used for shading is inverted
;   autocontour: in, optional, type=integer, scalar
;               Sets the number of contours to be generated automatically
;   cvalues    : in, optional, type=numeric, 1D array
;               Sets the ontour levels to be displayed (overrules the autocontour keyword)
;               If both autocontour and cvalues are omitted no contours are drawn (just the shaded plot is shown)
;   clabel     : in, optional, type=byte, scalar
;               If set the contour lines are labelled with their value.
;   downhill   : in, optional, type=byte, scalar
;               If set each contour will have short, perpendicular tick marks that point in the "downhill"
;               direction, making the direction of the grade readily apparent.
;   showscale  : in, optional, type=byte, scalar
;               If set the scale bar is shown
;   scalewidth : in, optional, type=float, scalar
;               Sets the width of the scale bar as fraction of width of plot.
;               Its range is 0.005 to 0.9 . This keyword is ignored is the scalepos keyword is defined [default=0.05]
;   scalelabel : in, optional, type=byte, scalar
;               Set the position of the scale bar and the scale label (0: bottom, 1: top, 2: left, 3: right)
;               If the scalepos keyword is defined, this will only set the position of the scale label [default=3]
;   scalepos   : in, optional, type=float, 4-element array
;               Sets the custom, normalised co-ordinates of the scale bar [x0,y0,x1,y1]
;   position   : in, optional, type=float, 4-element array
;               Position of plot on screen in normalised co-ordinates [x0,y0,x1,y1]
;   noerase    : in, optional, type=float, 4-element array
;               If set the currently existing plot is not erased (Just there for compatability)
;   error      : out, optional, type=string, scalar
;               Returns the error message if any. No IDL error is generated is this keyword is set.
;-

;****************************************
; routine to save/reset the color scheme
;****************************************
pro color_save, ctable, reset=reset

  common colors, r0, g0, b0, r1, g1, b1

  common loc_color_save, decomp, info

  if (!d.name ne 'X') then return

  if keyword_set(reset) then begin
    r0=info.red
    g0=info.green
    b0=info.blue
    tvlct, r0, g0, b0
    r1=info.red
    g1=info.green
    b1=info.blue
    !p.color=info.color
    !p.background=info.background
    device, decomp=decomp
    return
  endif

  device, get_decomp=decomp
  if undefined(r0) then begin
    r0=indgen(256)
    g0=indgen(256)
    b0=indgen(256)
    r1=indgen(256)
    g1=indgen(256)
    b1=indgen(256)
  endif
  info={red:r0,      $
        green:g0,    $
        blue:b0,     $
        color:!p.color, $
        background:!p.background}

  if (decomp eq 0) then begin
    i=(!p.background mod 256)
    !p.background=65536L*b0[i]+256*g0[i]+r0[i]
    i=(!p.color mod 256)
    !p.color=65536L*b0[i]+256*g0[i]+r0[i]
    device, decomp=1
  endif

end

;***************************************
; routine to load/reset the color table
;***************************************
pro cp_shaded_loadct, ctable, reset=reset
  if (!d.name eq 'X') then device, decomp=keyword_set(reset)
  if is_integer(ctable) then loadct, abs(ctable(0)) mod 40, /silent
end


;**************
; Main routine
;**************
pro cp_shaded,data, x, y,              $;
              xtitle=xtitle,           $;
              xrange=xrange,           $;
              ytitle=ytitle,           $;
              yrange=yrange,           $;
              title=title,             $;
              color=color,             $;
              ccolor=ccolor,           $;
              ctable=ctable,           $;
              invert=invert,           $;
              isotropic=isotropic,     $;
              cvalues=cvalues,         $;
              clabel=clabel,           $;
              autocontour=autocontour, $;
              zrange=zrange,           $;
              ztitle=ztitle,           $;
              downhill=downhill,       $;
              showscale=showscalearg,  $;
              scalewidth=scalewidtharg,$;
              scalepos=scaleposarg,    $;
              scalelabel=scalelabelarg,$;
              noerase=noerase,         $;
              position=position,       $;
              charsize=charsizearg,    $;
              error=error


  ; Set the character size
  oldcharsize=!p.charsize
  if keyword_set(charsizearg) then charsize=float(charsizearg) else charsize=0
  if charsize ne 0 then !p.charsize=(oldcharsize eq 0) ? charsize : oldcharsize*charsize $
                   else if (oldcharsize eq 0) then !p.charsize=1.

  ; We'll erase the existing plot - if required - by using the PLOT command
  ; with the NOERASE option, while suppressing the data and the axes.
  ; This overcomes the problem of using the erase command on a PS plot.
  ; This will also set the correct margins for character size
  xticks=!x.ticks                       ;Store X-tick value
  yticks=!y.ticks                       ;Store Y-tick value
  plot, [0,1], [0,1], /nodata, xs=4, ys=4, noerase=noerase
  !x.ticks=xticks                       ;Restore X-tick value
  !y.ticks=yticks                       ;Restore Y-tick value


  ; check if whether or not we need a shade plot
  if keyword_set(ctable) && (ctable lt 0) then begin
    showscalearg=0
    ctable=0
    noshade=1
  endif else noshade=0
  color_save, ctable

  ; check wether we need to show the scale bar
  showscalearg=keyword_set(showscalearg) or is_real(scalewidth)

  ; check whether we need to show contours (automatically generated or user-defined)
  customcontours=0B
  if (n_elements(cvalues) gt 0) then begin
    if not_real(cvalues) then begin
      error_message='CP_SHADED: C_LEVEL must be real-array'
      goto, errorcatch
    endif
    zvalues=cvalues
    customcontours=1B
  endif
  if is_integer(autocontour) && (autocontour gt 0) && ~customcontours then begin
    nvalues = autocontour
    autocontour=1B
  endif else autocontour=0B



  ; Check Arguement-1, must be 2D real array
  if not_real(data) then begin
    error_message='CP_SHADED: Arg(1) must be a real-array'
    goto, errorcatch
  endif
  dsize=size(data)
  if (dsize[0] ne 2) then begin
    error_message='CP_SHADED: Arg(1) must be a 2D real-array'
    goto, errorcatch
  endif
  dsize=size(data)
  if (dsize[1] le 1) or (dsize(2) le 1) then begin
    error_message='CP_SHADED: Arg(1) X/Y dimensions must be .gt. 1'
    goto, errorcatch
  endif

  ; Check Arguement-2, if defined must be 1D real array. 
  ; Length must be longer than X-dimension of data array.
  if (n_elements(x) eq 0) then begin
    ix=findgen(dsize[1])
  endif else begin
    if not_real(x) then begin
      error_message='CP_SHADED: Arg(2) must be a real-array'
      goto, errorcatch
    endif
    if (n_elements(x) lt dsize(1)) then begin
      error_message='CP_SHADED: Arg(2) size not compatible with Arg(1)'
      goto, errorcatch
    endif
    if (monotonic(x, /inc, /strict) eq 0) then begin
      error_message='CP_SHADED: Arg(2) must be monotonically increasing'
      goto, errorcatch
    endif
    ix=x[0:(dsize[1]-1)]
  endelse

  ; Check Arguement-3, if defined must be 1D real array.
  ; Length must be longer than Y-dimension of data array.
  if (n_elements(y) eq 0) then begin
    iy=findgen(dsize[2])
  endif else begin
    if not_real(y) then begin
      error_message='CP_SHADED: Arg(3) must be a real-array'
      goto, errorcatch
    endif
    if (n_elements(y) lt dsize(2)) then begin
      error_message='CP_SHADED: Arg(3) size not compatible with Arg(1)'
      goto, errorcatch
    endif
    if (monotonic(y, /inc, /strict) eq 0) then begin
      error_message='CP_SHADED: Arg(3) must be monotonically increasing'
      goto, errorcatch
    endif
    iy=y[0:(dsize[2]-1)]
  endelse

  ;Set XRANGE.
  ;     - IDL doesn't clip SHADE_SURF plots so that data array needs
  ;       to be cut down in size, if the xrange is set.
  ;     - Otherwise the XRANGE is calculated from the x-coordinate
  ;     - If no data is in the given XRANGE a blank plot is shown
  idata=data
  if is_real(xrange) then begin
    xmin=min(xrange, max=xmax)
    if (xmin eq xmax) then begin
      error_message='CP_SHADED: X-Range badly defined (xmin=xmax)'
      goto, errorcatch
    endif

    ixmin=min(ix, max=ixmax)
    if (xmin gt ixmax) or (xmax lt ixmin) then begin
      ymin=0.0
      ymax=1.0
      zmin=0.0
      zmax=1.0
      goto, endplot
    endif

    if (xmin gt ixmin) then begin
      i=value_locate(ix, xmin)
      idata[i, *]=idata[i, *]+(idata[i,*]-idata[i+1,*])*(xmin-ix[i])/(ix[i]-ix[i+1])
      idata=idata[i:*,*]
      ix=ix[i:*]
      ix[0]=xmin
    endif

    if (xmax lt ixmax) then begin
      i=value_locate(ix, xmax)
      idata[i+1, *]=idata[i,*]+(idata[i,*]-idata[i+1,*])*(xmax-ix[i])/(ix[i]-ix[i+1])
      idata=idata[0:(i+1), *]
      ix=ix[0:(i+1)]
      ix[i+1]=xmax
    endif
  endif else begin
    xmin=min(ix, max=xmax)
    xrange=[xmin,xmax]
  endelse

  ;Set YRANGE.
  ;     - IDL doesn't clip SHADE_SURF plots so that data array needs
  ;       to be cut down in size, if the yrange is set.
  ;     - Otherwise the YRANGE is calculated from the y-coordinate
  ;     - If no data is in the given YRANGE a blank plot is shown
  if is_real(yrange) then begin
    ymin=min(yrange, max=ymax)
    if (ymin eq ymax) then begin
      error_message='CP_SHADED: Y-Range badly defined (xmin=xmax)'
      goto, errorcatch
    endif

    iymin=min(iy, max=iymax)
    if (ymin gt iymax) or (ymax lt iymin) then goto, endplot

    if (ymin gt iymin) then begin
      i=value_locate(iy, ymin)
      idata[*, i]=idata[*, i]+(idata[*, i]-idata[*, i+1])*(ymin-iy[i])/(iy[i]-iy[i+1])
      idata=idata[*, i:*]
      iy=iy[i:*]
      iy[0]=ymin
    endif

    if (ymax lt iymax) then begin
      i=value_locate(iy, ymax)
      idata[*, i+1]=idata[*, i]+(idata[*, i]-idata[*, i+1])*(ymax-iy[i])/(iy[i]-iy[i+1])
      idata=idata[*, 0:(i+1)]
      iy=iy[0:(i+1)]
      iy[i+1]=ymax
    endif
  endif else begin
    ymin=min(iy, max=ymax)
    yrange=[ymin,ymax]
  endelse

  ;Check ZRANGE
  ;     - If ZRANGE is not defined a dummy call to PLOT is used
  ;       to determine the extent of the ZRANGE.
  ;     - The NOERASE option is used with the PLOT command, this over comes
  ;       the problem of using the erase command on a PS plot.
  if is_real(zrange) then begin
    zmin=min(zrange, max=zmax)
    if (zmin eq zmax) then begin
      error_message='CP_SHADED: Z-Range badly defined (xmin=xmax)'
      goto, errorcatch
    endif
  endif else begin
    zmin=min(idata, max=zmax, /nan)
    deltaz = zmax-zmin
    if deltaz eq 0 then begin
      print, 'CP_SHADED: WARNING: min(data) eq max(data)! Z-Range set to [min(data)-1,min(data)+1]'
      zmin -= 1.
      zmax += 1.
    endif else begin
      ; extend the z range by 5% at each end
      zmin -= 0.05*deltaz
      zmax += 0.05*deltaz
    endelse
  endelse


  ;Check position of plot on screen. 
  ;       - If ARG(position) isn't defined get the position of the dummy plot
  ;         from the system variables
  ;
  if not_real(position) then begin
    ; get the standard/previously set main plot position in normalised coordinates
    plotpos    = [!x.window, !y.window]
    plotpos    = plotpos[[0,2,1,3]]
  endif else begin
    ; set the custom main plot position in normalised coordinates
    if (n_elements(position) ne 4) then begin
      error_message='CP_SHADED: Position-arg wrong size position=[X0,Y0,X1,Y1]'
      goto, errorcatch
    endif
    if (min(position) lt 0) or (max(position) gt 1) then begin
      error_message='CP_SHADED: Position-arg wrong plot outside window'
      goto, errorcatch
    endif
    if (position[0] eq position[2]) or (position[1] eq position[3]) then begin
      error_message='CP_SHADED: Position-arg badly defined plot region'
      goto, errorcatch
    endif
    plotpos=position
    if (plotpos[0] ge plotpos[2]) then plotpos[[0,2]]=plotpos[[2,0]]
    if (plotpos[1] ge plotpos[3]) then plotpos[[1,3]]=plotpos[[3,1]]
  endelse
  ; initialy we set the scale position to the plot position (will be changed further on)
  scalepos  = plotpos

  ; if the plot is isotropic then recalculate the real position it will take
  if keyword_set(isotropic) then begin
    plot, [xmin, xmax], [ymin, ymax], $
          xstyle=5, ystyle=5, /isotropic, /nodata, /noerase
    isopos=[!x.window, !y.window]
    isopos=isopos[[0,2,1,3]]
    ; centre this with respect to the old plotpos
    xshift = 0.5*(plotpos[0]+plotpos[2])-0.5*(isopos[0]+isopos[2])
    yshift = 0.5*(plotpos[1]+plotpos[3])-0.5*(isopos[1]+isopos[3])
    isopos = isopos
    isopos[0] +=xshift
    isopos[2] +=xshift
    isopos[1] +=yshift
    isopos[3] +=yshift
    plotwidth  = (plotpos[2]- plotpos[0])
    isowidth   = (isopos[2] - isopos[0])
    plotheight = (plotpos[3]- plotpos[1])
    isoheight  = (isopos[3] - isopos[1])
    isoscale   = isoheight/isowidth
  endif

  ; get the standard character width and height
  ch_w  = !d.x_ch_size * !p.charsize  ; in device coordinates
  ch_h  = !d.y_ch_size * !p.charsize
  tmp   = convert_coord(ch_w,ch_h, /device, /to_normal)
  ch_w  = tmp[0]                      ; in normalised coordinates
  ch_h  = tmp[1]

  ;Position of ScaleBar:
  ;     - If the /SHOWSCALE option is choosen a scale will be placed on the 
  ;       right hand side of the plot. The default width is 1/20th the width of
  ;       the plot.
  ;     - The width can be changed using the SCALEWIDTH option
  ;     - SCALELABEL determines whether the colorbar is vertical (1,3) or horizontal (0,2),
  ;       where the label is placed and - if the SCALE_POSITION keyword is not set -
  ;       where the scale is placed  (0:bottom,1: top, 2: left, 3:right [default])
  ;     - The scale bar can be placed anywhere on the plot using the
  ;       SCALE_POSTION option
  if not_real(scalewidtharg) then scalewidth=0.05 else scalewidth=scalewidtharg[0]
  if n_elements(scalelabelarg) eq 0 then scalelabel=3 else scalelabel=floor(abs(scalelabelarg mod 4))
  if keyword_set(showscalearg) then begin
    if not_real(scaleposarg) then begin
      scalewith = 0.005>scalewidth<0.9
      case scalelabel of
      0: begin  ; bottom scale
           ; to make room for the plot label we add 4 characters of space in between the scale and the plot
           space_margin = 4*ch_h
           if keyword_set(isotropic) then begin
             totalheight = space_margin + isoheight*(1.+scalewidth)
             if totalheight gt plotheight then begin                  ; the isotropic plot including scalebar is
               isoheight = (plotheight-space_margin)/(1.+scalewidth)  ; too big for the original plot region
               isowidth  = isoheight/isoscale
             endif
             ; set the plot and scale regions (center the isotropic plot in the original plot region)
             scalepos[3] = scalepos[1] + scalewidth*isoheight
             xshift      = (0.5*(plotpos[0]+plotpos[2]))                      $ ; x-centre of the original plot region
                          -(plotpos[0]+0.5*isowidth)                            ; x-centre of the isotropic plot
             yshift      = (0.5*(plotpos[1]+plotpos[3]))                      $ ; y-centre of the original plot region
                          -(scalepos[3]+space_margin+0.5*isoheight)             ; y-centre of the isotropic plot
             xshift      = xshift>0
             yshift      = yshift>0
             plotpos[0]  = plotpos[0] + xshift
             plotpos[2]  = plotpos[0] + isowidth
             plotpos[1]  = plotpos[1] + (scalewidth)*isoheight + space_margin + yshift
             plotpos[3]  = plotpos[1] + isoheight
           endif else begin
             ; scaleheight = scalewidth*new_plotheight
             ; <=> new_plotpos[1] = old_plotpos[1] + scalewidth*(old_plotpos[3]-new_plotpos[1]) + space_margin
             plotpos[1]  = (plotpos[1] + scalewidth*plotpos[3] + space_margin)/(1.+scalewidth)
             scalepos[3] = scalepos[1] + scalewidth*(plotpos[3]-plotpos[1])
           endelse
         end
      1: begin  ; top scale
           ; to make room for the scale label we add 2 extra characters for the top margin. We also allow
           ; for 3 extra characters space in between the plot and the scale for the main plots title
           top_margin   = 2.0*ch_h
           space_margin = 3.0*ch_h
           if keyword_set(isotropic) then begin
             totalheight = top_margin + space_margin + isoheight*(1.+scalewidth)
             if totalheight gt plotheight then begin                             ; the isotropic plot including scalebar is
               isoheight = (plotheight-space_margin-top_margin)/(1.+scalewidth)  ; too big for the original plot region
               isowidth  = isoheight/isoscale
             endif
             ; set the plot and scale regions (center the isotropic plot in the original plot region)
             scalepos[3] = plotpos[3] - top_margin
             scalepos[1] = scalepos[3] - scalewidth*isoheight
             xshift      = (0.5*(plotpos[0]+plotpos[2]))                      $ ; x-centre of the original plot region
                          -(plotpos[0]+0.5*isowidth)                            ; x-centre of the isotropic plot
             yshift      = (0.5*(plotpos[1]+plotpos[3]))                      $ ; y-centre of the original plot region
                          -(scalepos[1]-space_margin-0.5*isoheight)             ; y-centre of the isotropic plot
             xshift      = xshift>0
             yshift      = yshift<0
             plotpos[0]  = plotpos[0] + xshift
             plotpos[2]  = plotpos[0] + isowidth
             plotpos[3]  = scalepos[1] - space_margin + yshift
             plotpos[1]  = plotpos[3] - isoheight
           endif else begin
             ; scaleheight = scalewidth*new_plotheigth
             ; <=> old_plotpos[3] = new_plotpos[3] + space_margin + scalewidth*(new_plotpos[3]-old_plotpos[1]) + top_margin
             plotpos[3]  = (plotpos[3] + scalewidth*plotpos[1] - space_margin - top_margin)/(1.+scalewidth)
             scalepos[1] = plotpos[3] + space_margin
             scalepos[3] = scalepos[1] + scalewidth*(plotpos[3]-plotpos[1])
           endelse
         end
      2: begin  ; left scale
           ; we add another 4 characters to the left margin to allow for the scale label
           ; to make room for the plot label we add 11 characters of space in between the scale and the plot
           left_margin  =  4.*ch_w
           space_margin = 11.*ch_w
           if keyword_set(isotropic) then begin
             totalwidth = left_margin + space_margin + isowidth*(1.+scalewidth)
             if totalwidth gt plotwidth then begin                               ; the isotropic plot including scalebar is
               isowidth  = (plotwidth-space_margin-left_margin)/(1.+scalewidth)  ; too big for the original plot region
               isoheight = isowidth*isoscale
             endif
             ; set the plot and scale regions (center the isotropic plot in the original plot region)
             scalepos[0] = scalepos[0] + left_margin
             scalepos[2] = scalepos[0] + scalewidth*isowidth
             xshift      = (0.5*(plotpos[0]+plotpos[2]))                      $ ; x-centre of the original plot region
                          -(scalepos[2]+space_margin+0.5*isowidth)              ; x-centre of the isotropic plot
             yshift      = (0.5*(plotpos[1]+plotpos[3]))                      $ ; y-centre of the original plot region
                          -(plotpos[1]+0.5*isoheight)                           ; y-centre of the isotropic plot
             xshift      = xshift>0
             yshift      = yshift>0
             plotpos[0]  = scalepos[2] + space_margin + xshift
             plotpos[2]  = plotpos[0]  + isowidth
             plotpos[1]  = plotpos[1] + yshift
             plotpos[3]  = plotpos[1] + isoheight
           endif else begin
             ; scalewidth = scalewidth*new_plotwidth
             ; <=> new_plotpos[0] = old_plotpos[0] + left_margin + scalewidth*(old_plotpos[2]-new_plotpos[0]) + space_margin
             plotpos[0]  = (plotpos[0] + scalewidth*plotpos[2] + left_margin + space_margin)/(1.+scalewidth)
             scalepos[0] = scalepos[0] + left_margin
             scalepos[2] = scalepos[0] + scalewidth*(plotpos[2]-plotpos[0])
           endelse
         end
      3: begin  ; right scale
           ; to make room for the scale label we add 9 extra characters for the right margin. We also allow 
           ; for 1 character space in between the plot and the scale
           right_margin = 9.0*ch_w
           space_margin = 1.0*ch_w
           if keyword_set(isotropic) then begin
             totalwidth = right_margin + space_margin + isowidth*(1.+scalewidth)
             if totalwidth gt plotwidth then begin                                ; the isotropic plot including scalebar is
               isowidth  = (plotwidth-space_margin-right_margin)/(1.+scalewidth)  ; too big for the original plot region
               isoheight = isowidth*isoscale
             endif
             ; set the plot and scale regions (center the isotropic plot in the original plot region)
             scalepos[2] = plotpos[2] - right_margin
             scalepos[0] = scalepos[2] - scalewidth*isowidth
             xshift      = (0.5*(plotpos[0]+plotpos[2]))                      $ ; x-centre of the original plot region
                          -(scalepos[0]-space_margin-0.5*isowidth)              ; x-centre of the isotropic plot
             yshift      = (0.5*(plotpos[1]+plotpos[3]))                      $ ; y-centre of the original plot region
                          -(plotpos[1]+0.5*isoheight)                           ; y-centre of the isotropic plot
             xshift      = xshift<0
             yshift      = yshift>0
             plotpos[2]  = scalepos[0] - space_margin + xshift
             plotpos[0]  = plotpos[2] - isowidth
             plotpos[1]  = plotpos[1] + yshift
             plotpos[3]  = plotpos[1] + isoheight
           endif else begin
             ; scalewidth = scalewidth*new_plotwidth
             ; <=> old_plotpos[2] = new_plotpos[2] + space_margin + scalewidth*(new_plotpos[2]-old_plotpos[0]) + right_margin
             plotpos[2]  = (plotpos[2] + scalewidth*plotpos[0] - space_margin - right_margin)/(1.+scalewidth)
             scalepos[0] = plotpos[2] + space_margin
             scalepos[2] = scalepos[0] + scalewidth*(plotpos[2]-plotpos[0])
           endelse
         end
      endcase
    endif else begin
      ;is an isotropic plot is required, then use the isotropic plot region
      if keyword_set(isotropic) then plotpos=isopos
      ; use a custom scale position (plot position is not adjusted; it is assumed that
      ; the plot position was adjusted correctly with !p.position or the position keyword)
      scalepos=scaleposarg
      if (n_elements(scalepos) ne 4) then begin
        error_message='CP_SHADED: Scale Position-arg wrong size scale_position=[X0,Y0,X1,Y1]'
        goto, errorcatch
      endif
      if (min(scalepos) lt 0) or (max(scalepos) gt 1) then begin
        error_message='CP_SHADED: Scale Position-arg wrong plot outside window'
        goto, errorcatch
      endif
      if (scalepos[0] eq scalepos[2]) or (scalepos[1] eq scalepos[3]) then begin
        error_message='CP_SHADED: Scale Position-arg bad size'
        goto, errorcatch
      endif

      if (scalepos[0] ge scalepos[2]) then scalepos[[0,2]]=scalepos[[2,0]]
      if (scalepos[1] ge scalepos[3]) then scalepos[[1,3]]=scalepos[[3,1]]
    endelse

    ; make the scale array
    if (scalelabel ge 2) then begin
      bx=[0,1]
      by=[zmin, zmax]
      bdata=[1,1]#by
    endif else begin
      bx=[zmin, zmax]
      by=[0,1]
      bdata=bx#[1,1]
    endelse
    ; plot the scale
    cp_shaded_loadct, ctable
    if keyword_set(invert) then begin
      shade_surf, bdata, bx, by,                               $
                position=scalepos,                             $
                shades=255B-bytscl(bdata, min=zmin, max=zmax), $
                ax=90, az=0, xstyle=5, ystyle=5, zstyle=5,     $
                /noerase, ztick_get=tickvalues
    endif else begin
      shade_surf, bdata, bx, by,                               $
                position=scalepos,                             $
                shades=bytscl(bdata, min=zmin, max=zmax),      $
                ax=90, az=0, xstyle=5, ystyle=5, zstyle=5,     $
                /noerase, ztick_get=tickvalues
    endelse
    cp_shaded_loadct, /reset
    tmpcolor=!p.color
    if exists(color) then !p.color=truecolor(color)
    ; print the scale label
    case scalelabel of 
      0:axis, xaxis=0, xstyle=1, xticklen=0.3, xtitle=ztitle
      1:axis, xaxis=1, xstyle=1, xticklen=0.3, xtitle=ztitle
      2:axis, yaxis=0, ystyle=1, yticklen=0.3, ytitle=ztitle
      3:axis, yaxis=1, ystyle=1, yticklen=0.3, ytitle=ztitle
    endcase
    ; plot a box around the scale
    plots, scalepos[[0,0,2,2,0]], scalepos[[1,3,3,1,1]], /norm
    !p.color=tmpcolor
  endif else begin
    ;is an isotropic plot is required, then use the isotropic plot region
    if keyword_set(isotropic) then plotpos=isopos
  endelse

  ; Plot the axes to the main plot
  if exists(color) then !p.color=truecolor(color)
  plot, [xmin, xmax], [ymin, ymax], /nodata, $
        position=plotpos,                    $
        xtitle=xtitle, xstyle=1,             $
        ytitle=ytitle, ystyle=1,             $
        title=title, /noerase

  ;Add shade-plot
  if ~noshade then begin
    cp_shaded_loadct, ctable

    if keyword_set(invert) then begin
      shade_surf, idata, ix, iy, /noerase,                       $
                  position=plotpos,                              $
                  shades=255B-bytscl(idata, min=zmin, max=zmax), $
                  ax=90, az=0,                                   $
                  xrange=xrange, xstyle=5,                       $
                  yrange=yrange, ystyle=5,                       $
                  zstyle=5, zrange=[zmin, zmax]
    endif else begin
      shade_surf, idata, ix, iy, /noerase,                       $
                  position=plotpos,                              $
                  shades=bytscl(idata, min=zmin, max=zmax),      $
                  ax=90, az=0,                                   $
                  xrange=xrange, xstyle=5,                       $
                  yrange=yrange, ystyle=5,                       $
                  zstyle=5, zrange=[zmin, zmax]
    endelse
    cp_shaded_loadct, /reset

  endif

  ;Add contours, define contour levels
  tmpcolor=!p.color
  if autocontour then begin
    if exists(ccolor) then !p.color=truecolor(ccolor)
    contour, idata, ix, iy, /overplot,  $
             nlevels=nvalues,           $
             downhill=downhill,         $
             c_label=replicate(keyword_set(clabel), nvalues+1),$
             c_charsize=0.75*!p.charsize
  endif
  if customcontours then begin
    if exists(ccolor) then !p.color=truecolor(ccolor)
    zvalues = zvalues[uniq(zvalues, sort(zvalues))]
    ord=where(zvalues gt zmin and zvalues lt zmax, num)
    if (num gt 0) then begin
      contour, idata, ix, iy, /overplot,  $
               level=zvalues[ord],        $
               downhill=downhill,         $
               c_label=replicate(keyword_set(clabel), num+1),$
               c_charsize=0.75*!p.charsize
    endif
  endif
  if is_integer(zerocontour) then begin
    if (zmin lt 0) and (zmax gt 0) then begin
      contour, idata, ix, iy, /overplot,  $
               level=0,                   $
               thick=1>zerocontour[0]<10,$
               c_label=1, col=0
    endif
  endif
  !p.color=tmpcolor

  ;Replot the axes to main-plot. This is done again at the end so that
  ;the graphics system variables contain the correct information.
endplot:

  if exists(color) then !p.color=truecolor(color)
  plot, [xmin, xmax], [ymin, ymax], /nodata, $
        position=plotpos,                    $
        xtitle=xtitle, xstyle=1,             $
        ytitle=ytitle, ystyle=1,             $
        title=title,                         $
        /noerase

  if exists(zmin) then !z.crange=[zmin, zmax] else !z.crange=0
  color_save, /reset

  ; reset the character size
  !p.charsize = oldcharsize

return

  errorcatch:

  color_save, /reset

  if arg_present(error) then begin
    error=error_message
    return
  endif

  error_message, error_message


end