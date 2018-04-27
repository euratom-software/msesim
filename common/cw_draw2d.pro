;*****************************************************************************************************
;*                                                                                                   *
;* A compound widget that allows to represent 2D data - e.g. z(time,radius) - as either:             *
;*     - a scrollable set of 'data as function of x' (e.g. time traces)                              *
;*     - a scrollable set of 'data as function of y' (e.g. radial profiles)                          *
;*     - a contour plot with color bar                                                               *
;*     - a surface plot with color bar                                                               *
;* 1D data can also be plotted, but obviously the data(y), contour                                   *
;* and surface plots are then disabled. (selecting, zooming, panning and saving still works)         *
;*                                                                                                   *
;* Usage:                                                                                            *
;* ------                                                                                            *
;* Creating the widget                                                                               *
;*   cw_id = cw_draw2d(base, xsize=xsize, ysize=ysize, frame=frame, event_pro=event_pro,             *
;*                      background=background, crossh=crossh, ctable=ctable)                         *
;*   with:                                                                                           *
;*     base      : the widget's parent                                                               *
;*     xsize     : xsize of the draw area in device coordinates (default=640)                        *
;*     ysize     : ysize of the draw area in device coordinates (default=480)                        *
;*     frame     : draws a frame around the widget (default=1)                                       *
;*     event_pro : the event routine called by the widget (default='').                              *
;*                 WARNING don't use 'widget_control' to set the event routine!                      *
;*     background: name of the standard background color as accepted by 'truecolor()')               *
;*                 [default='white']                                                                 *
;*     crossh    : switch on crossh [default=0]                                                      *
;*     ctable    : color table to use, -1 means no color shading is used by the contour              *
;*                 and surface plots [default=5]                                                     *
;*     cw_id     : the id of the cw_draw2d widget                                                    *
;*                                                                                                   *
;* Plottting the data                                                                                *
;*   widget_control, cw_id, set_value=plotdata                                                       *
;*                                                                                                   *
;*   with 'plotdata' a structure with following fields                                               *
;*     required fields:                                                                              *
;*       For 2D-data                                                                                 *
;*         .x   : [mx1]-array containing the x-dimension (e.g. time)                                 *
;*         .y   :    [nx1]-array containing the y-dimension (e.g. radius)                            *
;*                or [mxn]-array containing the y-dimension (e.g. time-dependent radius)             *
;*         .z   : [mxn]-array containing the data to be plotted                                      *
;*       For 1D-data                                                                                 *
;*         .x   : [mx1]-array containing the x-dimension (e.g. time)                                 *
;*         .z   : [mx1]-array containing the data to be plotted                                      *
;*     optional fields:                                                                              *
;*         .dz     : [mxn]-array containing the error on the data                                    *
;*         .xrange : sets the range of the x-axis                                                    *
;*         .yrange : sets the range of the y-axis                                                    *
;*         .zrange : sets the range of the z-axis                                                    *
;*         .xtitle : label of the x-axis (uses TeXtoIDL so you can use (some) LaTeX commands)        *
;*         .ytitle : label of the y-axis (uses TeXtoIDL so you can use (some) LaTeX commands)        *
;*         .ztitle : label of the z-axis (uses TeXtoIDL so you can use (some) LaTeX commands)        *
;*         .title  : title above the plot                                                            *
;*         .xsellbl: string that labels the 'data as function of x' selection button                 *
;*                   if omitted: plot.xselect = 'data(x)'                                            *
;*         .ysellbl: string that labels the 'data as function of y' selection button                 *
;*                   if omitted: plot.yselect = 'data(y)'                                            *
;*         .clevels: a string indicating either the number of levels used in the contour plot        *
;*                   - e.g. clevels='#10' - or the values of the contour levels                      *
;*                   - e.g. clevels='0.5,1.3,2.7'. (default='#20')                                   *
;*         .clabels: switch on contour labels                                                        *
;*         .xidx   : initial x-selection [default=0]                                                 *
;*         .yidx   : initial y-selection [default=0]                                                 *
;*         .idxonly: if set the set_value routine will only change the indices                       *
;*         .ptype  : plot type (0:data(x), 1:data(y), 2:contour, 3:surface) [default=0]              *
;*         .xgrid  : x-coordinates for x-gridlines to be plotted                                     *
;*         .ygrid  : y-coordinates for y-gridlines to be plotted                                     *
;*         .zgrid  : z-coordinates for z-gridlines to be plotted                                     *
;*         .linestyle: linestyle for the data(x) and data(y) plots [default=0]                       *
;*         .linethick: thickness of the lines in the data(x) and data(y) plots  [default=1]          *
;*                     (crosshair and gridlines always they have the standard thickness)             *
;*         .psym     : plot symbol for the lines in the data(x) and data(y) plots [default=0]        *
;*         .symsize  : size of the plot symbols [default=1]                                          *
;*         .errstyle : error style for the data(x) and data(y) plots:                                *
;*                      0=shaded area, 1=error bars [default=0]                                      *
;*         .charsize : character size for all plots symbol [default=1]                               *
;*         .bcolor : name of the background color (as accepted by 'truecolor()') [default='white']   *
;*         .fcolor : name of the foreground color (as accepted by 'truecolor()') [default='black']   *
;*         .lcolor : name of the line color (as accepted by 'truecolor()') [default='black']         *
;*         .ecolor : name of the color of the shaded error area if errstyle=0 [default='gray']       *
;*         .overplot : if set the data of this structure is overplotted to the data already there.   *
;*                     If no original data is present this field is ignored. This only has an effect *
;*                     on the data(x) and data(y)-plots, contour and surface plots will only show    *
;*                     the original data. The sliders and mouse selection will also only return the  *
;*                     indices of the original data.                                                 *
;*         .lgdpos : set the position of the automatic legend. 0=top-left, 1=top-right,              *
;*                   2=bottom-right, 3=bottom-left [default=1]                                       *
;*         .extralgd: settings for an extra legend. This only works on the first plot; it is not     *
;                     compatible with the overplot setting. The settings for the extra legend are    *
;                     given in the fields of the extralgd-structure                                  *
;*                    .text      = string array with the legend text                                 *
;*                    .color     = string array with the line colors for the legend                  *
;*                    .textcolor = string array with the text colors for the legend                  *
;*                    .linestyle = byte array with the linestyles for the legend                     *
;*                    .psym      = byte array with the symbols for the legend                        *
;*                    .symsize   = float array with the symbol size for the legend                   *
;*                    .thick     = float array with thickness of the legend                          *
;*                    .box       = byte scalar that draws a box around the legend                    *
;*                    .clear     = byte scalar that clears the area of the legend                    *
;*                    .position  = position of the legend. 0=top-left, 1=top-right,                  *
;*                                 2=bottom-right, 3=bottom-left                                     *
;*                    .horizontal= set a horizontal legend (columns) rather than a vertical (rows)   *
;*                    .charsize  = character size                                                    *
;*                    .charthick = character thickness                                               *
;*                    .spacing   = line spacing of the legend                                        *
;*                                                                                                   *
;*                                                                                                   *
;*                                                                                                   *
;* The widget has following functionallity:                                                          *
;* ----------------------------------------                                                          *
;*   TOP ROW CONTROLS                                                                                *
;*     - plot selector  : select different representations: data(x), data(y), contour or surface     *
;*       This generates an event with the index of the selected x and y (i.e. slide values)          *
;*       and the plot-type as fields:                                                                *
;*           EVENT={ID: 0l, TOP: 0l, HANDLER: 0l, XIDX: 0l, YIDX: 0l, PLOTTYPE:0l}                   *
;*     - croshair switch: plot a crosshair on currently selected data point                          *
;*     - color table    : select the color table to used for contour and surface                     *
;*     - #contours      : the number of contour lines to plot in contour and surface plots           *
;*     - labels switch  : print contour labels in contour and surface plots                          *
;*     - stick button   : makes current data(x) or datay(y) plot 'sticky'                            *
;*                        and gives it a different colour                                            *
;*     - clear button   : clears all sticky plots                                                    *
;*     - save button    : saves the current view to an EPS (vector) or to a PNG (bitmap) file        *
;*                                                                                                   *
;*   SLIDERS                                                                                         *
;*     - for a data(x) plot: the y-slider changes the y-value for which data(x) is plotted           *
;*                           the x-slider moves the crosshair if the crosshair switch is on          *
;*     - for a data(y) plot: the x-slider changes the x-value for which data(y) is plotted           *
;*                           the y-slider moves the crosshair if the crosshair switch is on          *
;*     - for a contour and surface plot: both x- and y-sliders moves the crosshair                   *
;*                                       if the crosshair switch is on                               *
;*     - will return an event with the index of the selected x and y (i.e. slide values)             *
;*       and the plot-type as fields:                                                                *
;*           EVENT={ID: 0l, TOP: 0l, HANDLER: 0l, XIDX: 0l, YIDX: 0l, PLOTTYPE:0l}                   *
;*                                                                                                   *
;*   MOUSE                                                                                           *
;*     - left mouse click: selects closest x and y value in data(x), data(y) and contour plots.      *
;*       This generates an event with the index of the selected x and y (i.e. slide values)          *
;*       and the plot-type as fields:                                                                *
;*           EVENT={ID: 0l, TOP: 0l, HANDLER: 0l, XIDX: 0l, YIDX: 0l, PLOTTYPE:0l}                   *
;*     - left mouse hold and drag: pans a zoomed in data(x), data(y) and contour plots               *
;*                                 rotates the surface plot                                          *
;*     - middle mouse hold and drag: zooms in to the selected area.                                  *
;*       In a contour plot both the main contour plot and the colorbar can be zoomed.                *
;*       In a surface plot only the colorbar can be zoomed                                           *
;*     - right mouse click: reset zoom to [min(x),max(x)], [min(x),max(x)], [min(z),max(z)]          *
;*                                                                                                   *
;*****************************************************************************************************


; some functions used in the widget:
;-----------------------------------
@legend      ; to plot the legend
@cp_shaded   ; the contour plot
@cp_surf     ; the surface plot
@save_dialog ; the save figure dialog


;******************************************************************************************************
; LINE COLORS FUNCTION
;******************************************************************************************************
function cw_draw2d_linecolors
  ; returns a reduced set of nice line colors, the lighter ones separated by 10 from the darker ones.
  ; this makes it easy to print the legend and to plot the 'errorbars' as a shaded area in the lighter color.
  colors=[!p.color                   ,$ ; 10 main line colours (first color => black/white)
          truecolor('crimson')       ,$
          truecolor('blue')          ,$
          truecolor('green')         ,$
          truecolor('darkcyan')      ,$
          truecolor('magenta')       ,$
          truecolor('orangered')     ,$
          truecolor('darkolivegreen'),$
          truecolor('indigo')        ,$
          truecolor('darkgoldenrod') ,$
          truecolor('gray')          ,$ ; 10 corresponding light colours
          truecolor('lightcoral')    ,$
          truecolor('lightblue')     ,$
          truecolor('lightgreen')    ,$
          truecolor('paleturquoise') ,$
          truecolor('lightpink')     ,$
          truecolor('lightsalmon')   ,$
          truecolor('darkseagreen')  ,$
          truecolor('lightsteelblue'),$
          truecolor('gold')           ]
  return, colors
end

;******************************************************************************************************
; EXTRA LEGEND ROUTINE
;******************************************************************************************************
pro cw_draw2d_extralegend, legend_settings, charsize_in=charsize_in
  ; plots an extra legend to the plot

  ; get all the (default if necessary) settings
  if ~tag_exists(legend_settings,'text') then return else text=legend_settings.text
  if ~tag_exists(legend_settings,'color') then color=0 else color=truecolor(legend_settings.color)
  if ~tag_exists(legend_settings,'textcolor') then textcolor=0 else textcolor=truecolor(legend_settings.textcolor)
  if ~tag_exists(legend_settings,'psym')  then psym=0 else psym=legend_settings.psym
  if ~tag_exists(legend_settings,'symsize') then symsize=replicate(1.,n_elements(text)) else symsize=legend_settings.symsize
  if ~tag_exists(legend_settings,'thick') then thick=1. else thick=legend_settings.thick
  if ~tag_exists(legend_settings,'box') then box=1 else box=legend_settings.box
  if ~tag_exists(legend_settings,'clear') then clear=1 else clear=legend_settings.clear
  if ~tag_exists(legend_settings,'horizontal') then horizontal=0 else horizontal=legend_settings.horizontal
  if ~tag_exists(legend_settings,'spacing') then spacing=2. else spacing=2.*legend_settings.spacing
  if ~tag_exists(legend_settings,'charthick') then charthick=1. else charthick=legend_settings.charthick
  if ~tag_exists(legend_settings,'charsize') then charsize=charsize_in $
                                             else charsize=legend_settings.charsize*charsize_in

  top=1 & bottom=0 & left=1 & right=0
  if tag_exists(legend_settings,'position') then begin
    case legend_settings.position of
    0: begin & top=1 & bottom=0 & left=1 & right=0 & end
    1: begin & top=1 & bottom=0 & left=0 & right=1 & end
    2: begin & top=0 & bottom=1 & left=0 & right=1 & end
    3: begin & top=0 & bottom=1 & left=1 & right=0 & end
    endcase
  endif

  if tag_exists(legend_settings,'linestyle') then begin
    legend, text, color=color, textcolor=textcolor,linestyle=legend_settings.linestyle,psym=psym,$
                  symsize=symsize, thick=thick, box=box, clear=clear, horizontal=horizontal,$
                  spacing=spacing, charthick=charthick, charsize=charsize,$
                  top=top,bottom=bottom,left=left,right=right
  endif else begin
    if psym eq 0 then begin
      legend, text, textcolor=textcolor, thick=thick, box=box, clear=clear, horizontal=horizontal,$
                    spacing=spacing, charthick=charthick, charsize=charsize,$
                    top=top,bottom=bottom,left=left,right=right
    endif else begin
      legend, text, textcolor=textcolor,psym=psym, color=color,$
                    symsize=symsize, thick=thick, box=box, clear=clear, horizontal=horizontal,$
                    spacing=spacing, charthick=charthick, charsize=charsize,$
                    top=top,bottom=bottom,left=left,right=right
    endelse
  endelse
end

;******************************************************************************************************
; CHECK PLOT DATA FUNCTION
;******************************************************************************************************
function cw_draw2d_check, data, istag=tag
  ; returns  1 if the data-structure is the correct 2D structure for this widget,
  ; returns -1 if the data-structure is the correct 1D structure for this widget,
  ; returns  0  if the data hasn't the correct structure
  ; if the keyword 'istag' is set to a potential tag-name: returns 1 if the tag is present, 0 if not

  if (size(data,/type) ne 8) then return,0

  if ~keyword_set(tag) then begin
    ; check the necessary tag names:
    required = ['x','z']
    check    = [ 0 , 0 ]
    for i=0,1 do begin
      if tag_exists(data,required[i]) then check[i]=1
    endfor
    if array_equal(check,[1,1]) then begin ; so the data is at least 1D
      check    = 0
      ; loop through the tag names and check the 'y'-tag required for 2D data is present
      if tag_exists(data,'y') then begin ; the data seems to be 2D ... do the sizes fit (and are ther bigger than 2)?
        sz = size(data.z, /dimensions)
        if size(data.y, /n_dim) le 1 then begin  ; y-dimension is constant
          if (n_elements(data.x) eq sz[0]) && (n_elements(data.y) eq sz[1]) then begin
            if (sz[0] gt 1) && (sz[0] gt 1) then return, 1 else begin
              print, 'ERROR: both x- and y-vectors should have more than 1 element!'
            endelse
          endif else begin
              print, 'ERROR: n_elements(x) ne n_elements(z[*,0]) or n_elements(y) ne n_elements(z[0,*])'
          endelse
        endif else begin                         ; y-dimension is dependent on x
          if n_elements(data.x) ne n_elements(data.y[*,0]) then begin
            print, ' ERROR: y is x-dependent, but n_elements(y[*,0]) ne n_elements(x)'
          endif else begin
            if (n_elements(data.x) eq sz[0]) && (n_elements(data.y[0,*]) eq sz[1]) then begin
              if (sz[0] gt 1) && (sz[0] gt 1) then return, 1 else begin
                print, 'ERROR: both x- and y-vectors should have more than 1 element!'
              endelse
            endif else begin
              print, 'ERROR: n_elements(x) ne n_elements(z[*,0]) or n_elements(y[0,*]) ne n_elements(z[0,*])'
            endelse
          endelse
        endelse
      endif else begin                   ; the data seems to be 1D ... do the sizes fit?
        if (n_elements(data.x) eq n_elements(data.z)) then begin
          if (n_elements(data.x) gt 1) then return, -1 else begin
            print, 'ERROR: both x-vector should have more than 1 element!'
          endelse
        endif else begin
            print, 'ERROR: n_elements(x) ne n_elements(z)!'
        endelse
      endelse
    endif
  endif else begin
    ; check whether the tag name is present
    return, tag_exists(data,tag)
  endelse
  ; no valid plot data or tag not found:
  return, 0
end


;******************************************************************************************************
; cw_draw2d_datax (plot data as function of x)
;******************************************************************************************************
pro cw_draw2d_datax, state

  ; get the necessary data
  id = state.id
  if state.ndata eq 0 then return
  data = *state.dataptr[0]

  ; get the settings of the draw widget
  widget_control, id.draw, get_uvalue=info
  ; get the sticky plots
  widget_control, id.stckopt.stick, get_uvalue=sticky

  ; get nice cw_draw2d_linecolors
  lincol = cw_draw2d_linecolors()

  ; get the x-, y-, z- and overall title
  if cw_draw2d_check(data,istag='xtitle') then xtitle=data.xtitle else xtitle='x'
  if cw_draw2d_check(data,istag='ytitle') then ytitle=data.ytitle else ytitle='y'
  if cw_draw2d_check(data,istag='ztitle') then ztitle=data.ztitle else ztitle='z'
  if cw_draw2d_check(data,istag='title')  then  title=data.title  else  title=''
  xtitle = xtitle
  ytitle = ytitle
  ztitle = ztitle
  title  = title
  ; get the plot styles
  if cw_draw2d_check(data,istag='linestyle') then linestyle=data.linestyle  else linestyle=0
  if cw_draw2d_check(data,istag='linethick') then linethick=data.linethick  else linethick=1
  if cw_draw2d_check(data,istag='psym')      then      psym=data.psym       else      psym=0
  if cw_draw2d_check(data,istag='symsize')   then   symsize=data.symsize    else   symsize=1
  if cw_draw2d_check(data,istag='errstyle')  then  errstyle=data.errstyle   else  errstyle=0
  if cw_draw2d_check(data,istag='lcolor')    then  lcolor  =truecolor(data.lcolor) else lcolor=lincol[0]
  if cw_draw2d_check(data,istag='ecolor')    then  ecolor  =truecolor(data.ecolor) else ecolor=lincol[10]
  if cw_draw2d_check(data,istag='charsize')  then  charsize=data.charsize   else  charsize=1.0


  ; plot the axis
  plot, info.xrange, info.zrange, /nodata, xs=1,ys=1,$
        xtitle=xtitle, ytitle=ztitle, title=title, charsize=charsize

  ; plot the sticky plots if any:
  lgdtxt = strarr(1)
  if sticky.yidx[0] ne -1 then begin
    lgdtxt = strarr(n_elements(sticky.yidx)+1)
    for i=0,n_elements(sticky.yidx)-1 do begin
      ; the original data
      ;------------------
      if cw_draw2d_check(data,istag='dz') then begin
        if errstyle eq 1 then begin
        ; plot errorbars
          errplot, data.x, data.z[*,sticky.yidx[i]]+data.dz[*,sticky.yidx[i]],$
                           data.z[*,sticky.yidx[i]]-data.dz[*,sticky.yidx[i]],$
                   linestyle=0,thick=linethick,col=lincol[i+1]
        endif else begin
        ; plot shaded error area
          xtmp = [data.x, reverse(data.x)]
          dtmp = [data.z[*,sticky.yidx[i]]+data.dz[*,sticky.yidx[i]],$
                  reverse(data.z[*,sticky.yidx[i]]-data.dz[*,sticky.yidx[i]])]
          polyfill, xtmp, dtmp, col=lincol[i+11], noclip=0
        endelse
      endif
      ; plot the actual data
      oplot, data.x, data.z[*,sticky.yidx[i]], col=lincol[i+1], linestyle=linestyle,$
             thick=linethick, psym=psym, symsize=symsize
      ; make legend text
      if size(data.y, /n_dim) le 1 then begin  ; y-dimension is constant
        lgdtxt[i+1] = strtrim(string(format='(A," = ",g-10.4)', ytitle,data.y[sticky.yidx[i]]))
      endif else begin                         ; y-dimension is x-dependent
        lgdtxt[i+1] = strtrim(string(format='(A," = ",g-10.4)', ytitle,data.y[info.xidx,sticky.yidx[i]]))
      endelse


      ; the overplot data (if any)
      ;---------------------------
      if state.ndata gt 1 then begin
        for j=1,state.ndata-1 do begin
          ; load overplot data
          odata = *state.dataptr[j]
          ; line thickness and style
          if cw_draw2d_check(odata,istag='linethick') then othick  =odata.linethick  else othick=1
          if cw_draw2d_check(odata,istag='linestyle') then ostyle  =odata.linestyle  else ostyle=2
          if cw_draw2d_check(odata,istag='errstyle')  then oestyle =odata.errstyle   else oestyle=0
          if cw_draw2d_check(odata,istag='psym')      then opsym   =odata.psym       else opsym=0
          if cw_draw2d_check(odata,istag='symsize')   then osymsize=odata.symsize    else osymsize=2
          ; get the yidx in odata closest to yidx in data
          if cw_draw2d_check(odata) eq -1 then yidx=0 else begin
            if size(data.y, /n_dim) le 1 && size(odata.y, /n_dim) le 1 then begin  ; y-dimension is constant for
              tmp=min(abs(odata.y - data.y[sticky.yidx[i]]), yidx)                 ; original and overplot data
            endif
            if size(data.y, /n_dim) le 1 && size(odata.y, /n_dim) gt 1 then begin  ; y-dimension is constant for
              tmp=min(abs(odata.x - data.x[sticky.xidx[i]]), xidx)                 ; original, but x-dependent
              tmp=min(abs(odata.y[xidx,*] - data.y[sticky.yidx[i]]), yidx)         ; for overplot data
            endif
            if size(data.y, /n_dim) gt 1 && size(odata.y, /n_dim) le 1 then begin  ; y-dimension is x-dependent for
              tmp=min(abs(odata.y - data.y[sticky.xidx[i],sticky.yidx[i]]), yidx)  ; original, but constant
            endif                                                                  ; for overplot data
            if size(data.y, /n_dim) gt 1 && size(odata.y, /n_dim) gt 1 then begin  ; y-dimension is x-dependent for
              tmp=min(abs(odata.x - data.x[sticky.xidx[i]]), xidx)                 ; original and overplot data
              tmp=min(abs(odata.y[xidx,*] - data.y[sticky.xidx[i],sticky.yidx[i]]), yidx)
            endif
          endelse

          if cw_draw2d_check(odata,istag='dz') then begin
            if oestyle eq 1 then begin
            ; plot errorbars
              errplot, odata.x, odata.z[*,yidx]+odata.dz[*,yidx],$
                                odata.z[*,yidx]-odata.dz[*,yidx],$
                       linestyle=0,thick=othick,col=lincol[i+1]
            endif else begin
            ; plot shaded error area
              xtmp = [odata.x, reverse(odata.x)]
              dtmp = [odata.z[*,yidx]+odata.dz[*,yidx],$
                      reverse(odata.z[*,yidx]-odata.dz[*,yidx])]
              polyfill, xtmp, dtmp, col=lincol[i+11], noclip=0
            endelse
          endif
          ; plot the actual data
          oplot, odata.x, odata.z[*,yidx], col=lincol[i+1], linestyle=ostyle,$
                 thick=othick, psym=opsym, symsize=osymsize
        endfor
      endif

    endfor
  endif


  ; the original data
  ;------------------
  ; plot the currently selected data
  if cw_draw2d_check(data,istag='dz') then begin
    if errstyle eq 1 then begin
    ; plot errorbars
      errplot, data.x, data.z[*,info.yidx]+data.dz[*,info.yidx],$
                       data.z[*,info.yidx]-data.dz[*,info.yidx],$
               linestyle=0,thick=linethick,col=lcolor
    endif else begin
    ; plot shaded error area
      xtmp = [data.x, reverse(data.x)]
      dtmp = [data.z[*,info.yidx]+data.dz[*,info.yidx],$
              reverse(data.z[*,info.yidx]-data.dz[*,info.yidx])]
      polyfill, xtmp, dtmp, col=ecolor, noclip=0
    endelse
  endif
  ; plot the actual data
  oplot, data.x, data.z[*,info.yidx], col=lcolor, linestyle=linestyle,$
             thick=linethick, psym=psym, symsize=symsize
  ; the overplot data (if any)
  ;---------------------------
  if state.ndata gt 1 then begin
    for j=1,state.ndata-1 do begin
      ; load overplot data
      odata = *state.dataptr[j]
      ; line thickness and style
      if cw_draw2d_check(odata,istag='linethick') then othick  =odata.linethick  else othick=1
      if cw_draw2d_check(odata,istag='linestyle') then ostyle  =odata.linestyle  else ostyle=2
      if cw_draw2d_check(odata,istag='errstyle')  then oestyle =odata.errstyle   else oestyle=0
      if cw_draw2d_check(odata,istag='psym')      then opsym   =odata.psym       else opsym=0
      if cw_draw2d_check(odata,istag='symsize')   then osymsize=odata.symsize    else osymsize=2
      if cw_draw2d_check(odata,istag='lcolor')    then olcolor =truecolor(odata.lcolor) else olcolor=lincol[0]
      if cw_draw2d_check(odata,istag='ecolor')    then oecolor =truecolor(odata.ecolor) else oecolor=lincol[10]
      ; get the yidx in odata closest to yidx in data
      if cw_draw2d_check(odata) eq -1 then yidx=0 else begin
        if size(data.y, /n_dim) le 1 && size(odata.y, /n_dim) le 1 then begin  ; y-dimension is constant for
          tmp=min(abs(odata.y - data.y[info.yidx]), yidx)                      ; original and overplot data
        endif
        if size(data.y, /n_dim) le 1 && size(odata.y, /n_dim) gt 1 then begin  ; y-dimension is constant for
          tmp=min(abs(odata.x - data.x[info.xidx]), xidx)                      ; original, but x-dependent
          tmp=min(abs(odata.y[xidx,*] - data.y[info.yidx]), yidx)              ; for overplot data
        endif
        if size(data.y, /n_dim) gt 1 && size(odata.y, /n_dim) le 1 then begin  ; y-dimension is x-dependent for
          tmp=min(abs(odata.y - data.y[info.xidx,info.yidx]), yidx)            ; original, but constant
        endif                                                                  ; for overplot data
        if size(data.y, /n_dim) gt 1 && size(odata.y, /n_dim) gt 1 then begin  ; y-dimension is x-dependent for
          tmp=min(abs(odata.x - data.x[info.xidx]), xidx)                      ; original and overplot data
          tmp=min(abs(odata.y[xidx,*] - data.y[info.xidx,info.yidx]), yidx)
        endif
      endelse

      if cw_draw2d_check(odata,istag='dz') then begin
        if oestyle eq 1 then begin
        ; plot errorbars
          errplot, odata.x, odata.z[*,yidx]+odata.dz[*,yidx],$
                            odata.z[*,yidx]-odata.dz[*,yidx],$
                    linestyle=0,thick=othick,col=olcolor
        endif else begin
        ; plot shaded error area
          xtmp = [odata.x, reverse(odata.x)]
          dtmp = [odata.z[*,yidx]+odata.dz[*,yidx],$
                  reverse(odata.z[*,yidx]-odata.dz[*,yidx])]
          polyfill, xtmp, dtmp, col=oecolor, noclip=0
        endelse
      endif
      ; plot the actual data
      oplot, odata.x, odata.z[*,yidx], col=olcolor, linestyle=ostyle,$
              thick=othick, psym=opsym, symsize=osymsize
    endfor
  endif

  ; make legend text (only in case of 2D data)
  if cw_draw2d_check(data) eq 1 then begin
    if size(data.y, /n_dim) le 1 then begin  ; y-dimension is constant
      lgdtxt[0] = strtrim(string(format='(A," = ",g-10.4)', ytitle,data.y[info.yidx]))
    endif else begin                         ; y-dimension is x-dependent
      lgdtxt[0] = strtrim(string(format='(A," = ",g-10.4)', ytitle,data.y[info.xidx,info.yidx]))
    endelse
  endif

  ; plot gridlines
  if cw_draw2d_check(data,istag='xgrid') then begin
    for i=0,n_elements(data.xgrid)-1 do oplot, [data.xgrid[i],data.xgrid[i]],info.zrange, col=lincol[0],linestyle=2
  endif
  if cw_draw2d_check(data,istag='zgrid') then begin
    for i=0,n_elements(data.zgrid)-1 do oplot, info.xrange,[data.zgrid[i],data.zgrid[i]], col=lincol[0],linestyle=2
  endif

  ; plot a crosshair if the crosshair switch is on
  widget_control, id.crosshair, get_value=crossh
  if crossh[0] then begin
    oplot, [data.x[info.xidx],data.x[info.xidx]], info.zrange,$
           col=truecolor('blue'), linestyle=2
    oplot, info.xrange, [data.z[info.xidx,info.yidx],data.z[info.xidx,info.yidx]],$
           col=truecolor('blue'), linestyle=2
  endif

  ; plot the legend(only in case of 2D data)
  top=1 & bottom=0 & left=0 & right=1
  if cw_draw2d_check(data,istag='lgdpos') then begin
    case data.lgdpos of
    0: begin & top=1 & bottom=0 & left=1 & right=0 & end
    1: begin & top=1 & bottom=0 & left=0 & right=1 & end
    2: begin & top=0 & bottom=1 & left=0 & right=1 & end
    3: begin & top=0 & bottom=1 & left=1 & right=0 & end
    endcase
  endif
  if cw_draw2d_check(data) eq 1 then $
     legend, lgdtxt, textcolors=[lincol,!p.background],charsize=charsize,$
             top=top,bottom=bottom,left=left,right=right,/box,/clear,spacing=2.0

  ; plot a possible extra legend
  if cw_draw2d_check(data,istag='extralgd') then cw_draw2d_extralegend, data.extralgd, charsize=charsize

end


;******************************************************************************************************
; cw_draw2d_datay (plot data as function of y)
;******************************************************************************************************
pro cw_draw2d_datay, state

  ; get the necessary data
  id = state.id
  if state.ndata eq 0 then return
  data = *state.dataptr[0]

  ; get the settings of the draw widget
  widget_control, id.draw, get_uvalue=info
  ; get the sticky plots
  widget_control, id.stckopt.stick, get_uvalue=sticky

  ; get nice cw_draw2d_linecolors
  lincol = cw_draw2d_linecolors()

  ; get the x-, y-, z- and overall title
  if cw_draw2d_check(data,istag='xtitle') then xtitle=data.xtitle else xtitle='x'
  if cw_draw2d_check(data,istag='ytitle') then ytitle=data.ytitle else ytitle='y'
  if cw_draw2d_check(data,istag='ztitle') then ztitle=data.ztitle else ztitle='z'
  if cw_draw2d_check(data,istag='title')  then  title=data.title  else  title=''
  xtitle = xtitle
  ytitle = ytitle
  ztitle = ztitle
  title  = title
  ; get the plot styles
  if cw_draw2d_check(data,istag='linestyle') then linestyle=data.linestyle  else linestyle=0
  if cw_draw2d_check(data,istag='linethick') then linethick=data.linethick  else linethick=1
  if cw_draw2d_check(data,istag='psym')      then      psym=data.psym       else      psym=0
  if cw_draw2d_check(data,istag='symsize')   then   symsize=data.symsize    else   symsize=1
  if cw_draw2d_check(data,istag='errstyle')  then  errstyle=data.errstyle   else  errstyle=0
  if cw_draw2d_check(data,istag='charsize')  then  charsize=data.charsize   else  charsize=1.0
  if cw_draw2d_check(data,istag='lcolor')    then  lcolor  =truecolor(data.lcolor) else lcolor=lincol[0]
  if cw_draw2d_check(data,istag='ecolor')    then  ecolor  =truecolor(data.ecolor) else ecolor=lincol[10]

  ; plot the axis
  plot, info.yrange, info.zrange, /nodata, xs=1,ys=1,$
        xtitle=ytitle, ytitle=ztitle, title=title, charsize=charsize

  ; plot the sticky plots if any:
  lgdtxt = strarr(1)
  if sticky.xidx[0] ne -1 then begin
    lgdtxt = strarr(n_elements(sticky.xidx)+1)
    for i=0,n_elements(sticky.xidx)-1 do begin
      ; the original data
      ;------------------
      if cw_draw2d_check(data,istag='dz') then begin
        if errstyle eq 1 then begin
        ; plot errorbars
          if size(data.y, /n_dim) le 1 then begin  ; y-dimension is constant
            errplot, data.y, data.z[sticky.xidx[i],*]+data.dz[sticky.xidx[i],*],$
                             data.z[sticky.xidx[i],*]-data.dz[sticky.xidx[i],*],$
                     linestyle=0,thick=linethick,col=lincol[i+1]
          endif else begin                         ; y-dimension is x-dependent
            errplot, data.y[sticky.xidx[i],*], data.z[sticky.xidx[i],*]+data.dz[sticky.xidx[i],*],$
                                               data.z[sticky.xidx[i],*]-data.dz[sticky.xidx[i],*],$
                     linestyle=0,thick=linethick,col=lincol[i+1]
          endelse
        endif else begin
        ; plot shaded error area
          if size(data.y, /n_dim) le 1 then begin  ; y-dimension is constant
            xtmp = [data.y, reverse(data.y)]
          endif else begin                         ; y-dimension is x-dependent
            xtmp = [reform(data.y[sticky.xidx[i],*]), reverse(reform(data.y[sticky.xidx[i],*]))]
          endelse
          dtmp = [reform(data.z[sticky.xidx[i],*]+data.dz[sticky.xidx[i],*]),$
                  reverse(reform(data.z[sticky.xidx[i],*]-data.dz[sticky.xidx[i],*]))]
          polyfill, xtmp, dtmp, col=lincol[i+11], noclip=0
        endelse
      endif
      ; plot the actual data
      if size(data.y, /n_dim) le 1 then begin  ; y-dimension is constant
        oplot, data.y, data.z[sticky.xidx[i],*], col=lincol[i+1], linestyle=linestyle,$
               thick=linethick, psym=psym, symsize=symsize
      endif else begin                         ; y-dimension is x-dependent
        oplot, data.y[sticky.xidx[i],*], data.z[sticky.xidx[i],*], col=lincol[i+1], linestyle=linestyle,$
               thick=linethick, psym=psym, symsize=symsize
      endelse
      ; make legend text
      lgdtxt[i+1] = strtrim(string(format='(A," = ",g-10.4)', xtitle,data.x[sticky.xidx[i]]))

      ; the overplot data (if any)
      ;---------------------------
      if state.ndata gt 1 then begin
        for j=1,state.ndata-1 do begin
          ; load overplot data
          odata = *state.dataptr[j]
          ; line thickness and style
          if cw_draw2d_check(odata,istag='linethick') then othick  =odata.linethick  else othick=1
          if cw_draw2d_check(odata,istag='linestyle') then ostyle  =odata.linestyle  else ostyle=2
          if cw_draw2d_check(odata,istag='errstyle')  then oestyle =odata.errstyle   else oestyle=0
          if cw_draw2d_check(odata,istag='psym')      then opsym   =odata.psym       else opsym=0
          if cw_draw2d_check(odata,istag='symsize')   then osymsize=odata.symsize    else osymsize=2
          ; get the xidx in odata closest to xidx in data
          if cw_draw2d_check(odata) eq -1 then break ; no 2D data!
          tmp=min(abs(odata.x - data.x[sticky.xidx[i]]), xidx)

          if cw_draw2d_check(odata,istag='dz') then begin
            if oestyle eq 1 then begin
            ; plot errorbars
              if size(odata.y, /n_dim) le 1 then begin  ; y-dimension is constant
                errplot, odata.y, odata.z[xidx,*]+odata.dz[xidx,*],$
                                  odata.z[xidx,*]-odata.dz[xidx,*],$
                         linestyle=0,thick=othick,col=lincol[i+1]
              endif else begin                         ; y-dimension is x-dependent
                errplot, odata.y[xidx,*], odata.z[xidx,*]+odata.dz[xidx,*],$
                                          odata.z[xidx,*]-odata.dz[xidx,*],$
                         linestyle=0,thick=othick,col=lincol[i+1]
              endelse
            endif else begin
            ; plot shaded error area
              if size(odata.y, /n_dim) le 1 then begin  ; y-dimension is constant
                xtmp = [odata.y, reverse(odata.y)]
              endif else begin                         ; y-dimension is x-dependent
                xtmp = [reform(odata.y[xidx,*]), reverse(reform(odata.y[xidx,*]))]
              endelse
              dtmp = [reform(odata.z[xidx,*]+odata.dz[xidx,*]),$
                      reverse(reform(odata.z[xidx,*]-odata.dz[xidx,*]))]
              polyfill, xtmp, dtmp, col=lincol[i+11], noclip=0
            endelse
          endif
          ; plot the actual data
          if size(odata.y, /n_dim) le 1 then begin  ; y-dimension is constant
            oplot, odata.y, odata.z[xidx,*], col=lincol[i+1], linestyle=ostyle,$
                   thick=othick, psym=opsym, symsize=osymsize
          endif else begin                         ; y-dimension is x-dependent
            oplot, odata.y[xidx,*], odata.z[xidx,*], col=lincol[i+1], linestyle=ostyle,$
                   thick=othick, psym=opsym, symsize=osymsize
          endelse
        endfor
      endif


    endfor
  endif

  ; plot the currently selected data
  if cw_draw2d_check(data,istag='dz') then begin
    if errstyle eq 1 then begin
    ; plot errorbars
      if size(data.y, /n_dim) le 1 then begin  ; y-dimension is constant
        errplot, data.y, data.z[info.xidx,*]+data.dz[info.xidx,*],$
                         data.z[info.xidx,*]-data.dz[info.xidx,*],$
                 linestyle=0,thick=linethick,col=lcolor
      endif else begin                         ; y-dimension is x-dependent
        errplot, data.y[info.xidx,*], data.z[info.xidx,*]+data.dz[info.xidx,*],$
                                      data.z[info.xidx,*]-data.dz[info.xidx,*],$
                 linestyle=0,thick=linethick,col=lcolor
      endelse
    endif else begin
    ; plot shaded error area
      if size(data.y, /n_dim) le 1 then begin  ; y-dimension is constant
        xtmp = [data.y, reverse(data.y)]
      endif else begin                         ; y-dimension is x-dependent
        xtmp = [reform(data.y[info.xidx,*]), reverse(reform(data.y[info.xidx,*]))]
      endelse
      dtmp = [reform(data.z[info.xidx,*]+data.dz[info.xidx,*]),$
              reverse(reform(data.z[info.xidx,*]-data.dz[info.xidx,*]))]
      polyfill, xtmp, dtmp, col=ecolor, noclip=0
    endelse
  endif
  ; plot the actual data
  if size(data.y, /n_dim) le 1 then begin  ; y-dimension is constant
    oplot, data.y, data.z[info.xidx,*], col=lcolor, linestyle=linestyle,$
           thick=linethick, psym=psym, symsize=symsize
  endif else begin                         ; y-dimension is x-dependent
    oplot, data.y[info.xidx,*], data.z[info.xidx,*], col=lcolor, linestyle=linestyle,$
           thick=linethick, psym=psym, symsize=symsize
  endelse

  ; the overplot data (if any)
  ;---------------------------
  if state.ndata gt 1 then begin
    for j=1,state.ndata-1 do begin
      ; load overplot data
      odata = *state.dataptr[j]
      ; line thickness and style
      if cw_draw2d_check(odata,istag='linethick') then othick  =odata.linethick  else othick=1
      if cw_draw2d_check(odata,istag='linestyle') then ostyle  =odata.linestyle  else ostyle=2
      if cw_draw2d_check(odata,istag='errstyle')  then oestyle =odata.errstyle   else oestyle=0
      if cw_draw2d_check(odata,istag='psym')      then opsym   =odata.psym       else opsym=0
      if cw_draw2d_check(odata,istag='symsize')   then osymsize=odata.symsize    else osymsize=2
      if cw_draw2d_check(odata,istag='lcolor')    then olcolor =truecolor(odata.lcolor) else olcolor=lincol[0]
      if cw_draw2d_check(odata,istag='ecolor')    then oecolor =truecolor(odata.ecolor) else oecolor=lincol[10]
      ; get the yidx in odata closest to yidx in data
      if cw_draw2d_check(odata) eq -1 then break  ; no 2D data!
      tmp=min(abs(odata.x - data.x[info.xidx]), xidx)

      if cw_draw2d_check(odata,istag='dz') then begin
        if oestyle eq 1 then begin
        ; plot errorbars
          if size(odata.y, /n_dim) le 1 then begin  ; y-dimension is constant
            errplot, odata.y, odata.z[xidx,*]+odata.dz[xidx,*],$
                              odata.z[xidx,*]-odata.dz[xidx,*],$
                      linestyle=0,thick=othick,col=olcolor
          endif else begin                         ; y-dimension is x-dependent
            errplot, odata.y[xidx,*], odata.z[xidx,*]+odata.dz[xidx,*],$
                                      odata.z[xidx,*]-odata.dz[xidx,*],$
                      linestyle=0,thick=othick,col=olcolor
          endelse
        endif else begin
        ; plot shaded error area
          if size(odata.y, /n_dim) le 1 then begin  ; y-dimension is constant
            xtmp = [odata.y, reverse(odata.y)]
          endif else begin                         ; y-dimension is x-dependent
            xtmp = [reform(odata.y[xidx,*]), reverse(reform(odata.y[xidx,*]))]
          endelse
          dtmp = [reform(odata.z[xidx,*]+odata.dz[xidx,*]),$
                  reverse(reform(odata.z[xidx,*]-odata.dz[xidx,*]))]
          polyfill, xtmp, dtmp, col=oecolor, noclip=0
        endelse
      endif
      ; plot the actual data
      if size(odata.y, /n_dim) le 1 then begin  ; y-dimension is constant
        oplot, odata.y, odata.z[xidx,*], col=olcolor, linestyle=ostyle,$
                thick=othick, psym=opsym, symsize=osymsize
      endif else begin                         ; y-dimension is x-dependent
        oplot, odata.y[xidx,*], odata.z[xidx,*], col=olcolor, linestyle=ostyle,$
                thick=othick, psym=opsym, symsize=osymsize
      endelse
    endfor
  endif

  ; make legend text
  lgdtxt[0] = strtrim(string(format='(A," = ",g-10.4)', xtitle,data.x[info.xidx]))

  ; plot gridlines
  if cw_draw2d_check(data,istag='ygrid') then begin
    for i=0,n_elements(data.ygrid)-1 do oplot, [data.ygrid[i],data.ygrid[i]],info.zrange,col=lincol[0],linestyle=2
  endif
  if cw_draw2d_check(data,istag='zgrid') then begin
    for i=0,n_elements(data.zgrid)-1 do oplot, info.yrange,[data.zgrid[i],data.zgrid[i]],col=lincol[0],linestyle=2
  endif

  ; plot a crosshair if the crosshair switch is on
  widget_control, id.crosshair, get_value=crossh
  if crossh[0] then begin
    if size(data.y, /n_dim) le 1 then begin  ; y-dimension is constant
      oplot, [data.y[info.yidx],data.y[info.yidx]], info.zrange,$
             col=truecolor('blue'), linestyle=2
    endif else begin                         ; y-dimension is x-dependent
      oplot, [data.y[info.xidx,info.yidx],data.y[info.xidx,info.yidx]], info.zrange,$
             col=truecolor('blue'), linestyle=2
    endelse
    oplot, info.yrange, [data.z[info.xidx,info.yidx],data.z[info.xidx,info.yidx]],$
           col=truecolor('blue'), linestyle=2

  endif

  ; plot the legend
  top=1 & bottom=0 & left=0 & right=1
  if cw_draw2d_check(data,istag='lgdpos') then begin
    case data.lgdpos of
    0: begin & top=1 & bottom=0 & left=1 & right=0 & end
    1: begin & top=1 & bottom=0 & left=0 & right=1 & end
    2: begin & top=0 & bottom=1 & left=0 & right=1 & end
    3: begin & top=0 & bottom=1 & left=1 & right=0 & end
    endcase
  endif
  legend, lgdtxt, textcolors=[lincol,!p.background],charsize=charsize,$
          top=top,bottom=bottom,left=left,right=right,/box,/clear,spacing=2.0

  ; plot a possible extra legend
  if cw_draw2d_check(data,istag='extralgd') then cw_draw2d_extralegend, data.extralgd, charsize=charsize
end


;******************************************************************************************************
; cw_draw2d_contour (plot data as contour)
;******************************************************************************************************
pro cw_draw2d_contour, state

  ; get the necessary data
  id = state.id
  if state.ndata eq 0 then return
  data = *state.dataptr[0]

  ; get the settings of the draw widget
  widget_control, id.draw, get_uvalue=info

  ; get the x-, y-, z- and overall title
  if cw_draw2d_check(data,istag='xtitle') then xtitle=data.xtitle else xtitle='x'
  if cw_draw2d_check(data,istag='ytitle') then ytitle=data.ytitle else ytitle='y'
  if cw_draw2d_check(data,istag='ztitle') then ztitle=data.ztitle else ztitle='z'
  if cw_draw2d_check(data,istag='title')  then  title=data.title  else  title=''
  xtitle = xtitle
  ytitle = ytitle
  ztitle = ztitle
  title  = title
  ; get the plot styles
  if cw_draw2d_check(data,istag='linethick') then linethick=data.linethick  else linethick=1
  if cw_draw2d_check(data,istag='charsize')  then  charsize=data.charsize   else  charsize=1.0

  ; get the contour settings
  widget_control, id.contopt.clevels, get_value=clevels
  clevels=strtrim(clevels[0],2)
  if strcmp(clevels,'#',1) || strcmp(strtrim(clevels,2),'') then begin
    if strcmp(strtrim(clevels,2),'') then nlevels=0 else nlevels=fix(strmid(clevels,1)) <99 >0
    nlevels=fix(strmid(clevels,1)) <99 >0
    if nlevels gt 0 then begin
      cmin = info.zrange[0] + 0.05*(info.zrange[1] - info.zrange[0])
      cmax = info.zrange[1] - 0.05*(info.zrange[1] - info.zrange[0])
      cvalues = cmin + findgen(nlevels)/(nlevels-1.0)*(cmax-cmin)
      widget_control, id.contopt.clabels, get_value=clabels
      clabel=clabels[0]
    endif
    widget_control, id.contopt.clevels, set_value=string(format='("#",I0)', nlevels)
  endif else begin
    cvalues=float(strsplit(clevels,'[, ]', /extract))
    widget_control, id.contopt.clabels, get_value=clabels
    clabel=clabels[0]
  endelse

  ; erase the current plot
  erase

  ; sort the x and y axes
  xsort = sort(data.x)
  xtmp  = data.x[xsort]
  if size(data.y,/n_dim) le 1 then begin  ; y-dimension is constant
    ysort = sort(data.y)
    ytmp  = data.y[ysort]
    ztmp  = data.z[xsort,*]
    ztmp  = ztmp[*,ysort]
  endif else begin                        ; y-dimension is x-dependent
    miny  = min(data.y)
    maxy  = max(data.y)
    n_y   = n_elements(data.y[0,*])
    ytmp  = miny+findgen(n_y)/float(n_y-1.)*(maxy-miny)
    ysort = sort(data.y[0,*])
    y0    = data.y[xsort,*]
    y0    = y0[*,ysort]
    z0    = data.z[xsort,*]
    z0    = z0[*,ysort]
    ztmp  = fltarr(n_elements(xtmp),n_elements(ytmp))
    for i=0,n_elements(xtmp)-1 do begin
      ztmp[i,*] = interpol(z0[i,*],y0[i,*],ytmp)
    endfor
  endelse
  ; make the contour plot
  cp_shaded, ztmp, xtmp, ytmp, $
            xrange=info.xrange, yrange=info.yrange, zrange=info.zrange,$
            xtitle=xtitle, ytitle=ytitle, ztitle=ztitle, title=title, charsize=charsize,$
            /showscale, cvalues=cvalues, clabel=clabel, ctable=info.ctable

  ; plot gridlines
  if cw_draw2d_check(data,istag='xgrid') then begin
    for i=0,n_elements(data.xgrid)-1 do oplot, [data.xgrid[i],data.xgrid[i]],info.yrange,linestyle=2
  endif
  if cw_draw2d_check(data,istag='ygrid') then begin
    for i=0,n_elements(data.ygrid)-1 do oplot, info.xrange,[data.ygrid[i],data.ygrid[i]],linestyle=2
  endif

  ; plot a crosshair if the crosshair switch is on
  widget_control, id.crosshair, get_value=crossh
  if crossh[0] then begin
    oplot, [data.x[info.xidx],data.x[info.xidx]], info.yrange,$
           col=truecolor('blue'), linestyle=2
    if size(data.y, /n_dim) le 1 then begin  ; y-dimension is constant
      oplot, info.xrange, [data.y[info.yidx],data.y[info.yidx]],$
             col=truecolor('blue'), linestyle=2
    endif else begin                         ; y-dimension is x-dependent
      oplot, info.xrange, [data.y[info.xidx,info.yidx],data.y[info.xidx,info.yidx]],$
             col=truecolor('blue'), linestyle=2
    endelse
  endif

end


;******************************************************************************************************
; cw_draw2d_surf (plot data as surface)
;******************************************************************************************************
pro cw_draw2d_surf, state

  ; get the necessary data
  id = state.id
  if state.ndata eq 0 then return
  data = *state.dataptr[0]

  ; get the settings of the draw widget
  widget_control, id.draw, get_uvalue=info

  ; get the x-, y-, z- and overall title
  if cw_draw2d_check(data,istag='xtitle') then xtitle=data.xtitle else xtitle='x'
  if cw_draw2d_check(data,istag='ytitle') then ytitle=data.ytitle else ytitle='y'
  if cw_draw2d_check(data,istag='ztitle') then ztitle=data.ztitle else ztitle='z'
  if cw_draw2d_check(data,istag='title')  then  title=data.title  else  title=''
  xtitle = xtitle
  ytitle = ytitle
  ztitle = ztitle
  title  = title
  ; get the plot styles
  if cw_draw2d_check(data,istag='linethick') then linethick=data.linethick  else linethick=1
  if cw_draw2d_check(data,istag='charsize')  then  charsize=data.charsize   else  charsize=1.0

  ; get the contour settings
  widget_control, id.contopt.clevels, get_value=clevels
  clevels=strtrim(clevels[0],2)
  if strcmp(clevels,'#',1) || strcmp(strtrim(clevels,2),'') then begin
    if strcmp(strtrim(clevels,2),'') then nlevels=0 else nlevels=fix(strmid(clevels,1)) <99 >0
    if nlevels gt 0 then begin
      cmin = info.zrange[0] + 0.05*(info.zrange[1] - info.zrange[0])
      cmax = info.zrange[1] - 0.05*(info.zrange[1] - info.zrange[0])
      cvalues = cmin + findgen(nlevels)/(nlevels-1.0)*(cmax-cmin)
      widget_control, id.contopt.clabels, get_value=clabels
      clabel=clabels[0]
    endif
    widget_control, id.contopt.clevels, set_value=string(format='("#",I0)', nlevels)
  endif else begin
    cvalues=float(strsplit(clevels,'[, ]', /extract))
    widget_control, id.contopt.clabels, get_value=clabels
    clabel=clabels[0]
  endelse

  ; erase the current plot
  erase

  ; sort the x and y axes
  xsort = sort(data.x)
  xtmp  = data.x[xsort]
  if size(data.y,/n_dim) le 1 then begin  ; y-dimension is constant
    ysort = sort(data.y)
    ytmp  = data.y[ysort]
    ztmp  = data.z[xsort,*]
    ztmp  = ztmp[*,ysort]
  endif else begin                        ; y-dimension is x-dependent
    miny  = min(data.y)
    maxy  = max(data.y)
    n_y   = n_elements(data.y[0,*])
    ytmp  = miny+findgen(n_y)/float(n_y-1.)*(maxy-miny)
    ysort = sort(data.y[0,*])
    y0    = data.y[xsort,*]
    y0    = y0[*,ysort]
    z0    = data.z[xsort,*]
    z0    = z0[*,ysort]
    ztmp  = fltarr(n_elements(xtmp),n_elements(ytmp))
    for i=0,n_elements(xtmp)-1 do begin
      ztmp[i,*] = interpol(z0[i,*],y0[i,*],ytmp)
    endfor
  endelse
  ; make the surface plot
  cp_surf, transpose(ztmp), ytmp, xtmp, ax=info.ax, az=info.az,$
            xrange=info.yrange, yrange=info.xrange, zrange=info.zrange,$
            xtitle=ytitle, ytitle=xtitle, ztitle=ztitle, title=title, charsize=charsize,$
            /showscale, cvalues=cvalues, clabel=clabel, ctable=info.ctable

  ; no gridlines in surface plot (gets messy)

  ; plot a crosshair if the crosshair switch is on
  widget_control, id.crosshair, get_value=crossh
  if crossh[0] then begin
    plots, info.yrange,[data.x[info.xidx],data.x[info.xidx]],[info.zrange[1],info.zrange[1]],$
           /data,/T3D,col=truecolor('blue'),linestyle=2, noclip=0
    if size(data.y, /n_dim) le 1 then begin  ; y-dimension is constant
      plots, [data.y[info.yidx],data.y[info.yidx]],info.xrange,[info.zrange[1],info.zrange[1]],$
             /data,/T3D,col=truecolor('blue'),linestyle=2, noclip=0
    endif else begin                         ; y-dimension is x-dependent
      plots, [data.y[info.xidx,info.yidx],data.y[info.xidx,info.yidx]],info.xrange,$
             [info.zrange[1],info.zrange[1]],$
             /data,/T3D,col=truecolor('blue'),linestyle=2, noclip=0
    endelse
    if size(data.y, /n_dim) le 1 then begin  ; y-dimension is constant
      plots, [data.y[info.yidx],data.y[info.yidx]],[data.x[info.xidx],data.x[info.xidx]],$
             [info.zrange[1],data.z[info.xidx,info.yidx]],/data,/T3D,col=truecolor('blue'),$
             linestyle=2, noclip=0
    endif else begin                         ; y-dimension is x-dependent
      plots, [data.y[info.xidx,info.yidx],data.y[info.xidx,info.yidx]],$
             [data.x[info.xidx],data.x[info.xidx]],$
             [info.zrange[1],data.z[info.xidx,info.yidx]],/data,/T3D,col=truecolor('blue'),$
             linestyle=2, noclip=0
    endelse

  endif

end

;******************************************************************************************************
; cw_draw2d_plot ROUTINE (main routine for plotting to the draw widget)
;******************************************************************************************************
pro cw_draw2d_plot, state
  ; get the necessary data
  id = state.id

  ; get the settings of the draw widget
  widget_control, id.draw, get_uvalue=info

  ; select the draw widget from plotting
  widget_control, id.draw, get_value =win
  wset, win
  !x.margin = [14,5]

  ; get the selected plot type
  plotsel = widget_info(id.selector,/droplist_select)
  case plotsel of
    0: begin
         cw_draw2d_datax, state
       end
    1: begin
         cw_draw2d_datay, state
       end
    2: begin
         cw_draw2d_contour, state
       end
    3: begin
         cw_draw2d_surf, state
       end
  endcase
  ; save the scaling and position of the current plot
  info.scale = [!x.s,!y.s]
  pos=[!x.window, !y.window]
  pos=pos[[0,2,1,3]]
  info.pos   = pos
  widget_control, id.draw, set_uvalue=info

end


;******************************************************************************************************
; cw_draw2d_select ROUTINE
;******************************************************************************************************
pro cw_draw2d_select, state

  ; get the necessary data
  id = state.id
  if state.ndata eq 0 then return
  data = *state.dataptr[0]

  ; get the draw info structure
  widget_control, id.draw, get_uvalue=info
  ; get the mouse positions
  mouse = [info.mousepos.x0, info.mousepos.y0]

  ; convert the mouse device coordinates into data coordinates
  !x.s = info.scale[0:1]
  !y.s = info.scale[2:3]
  !y.type = 0
  area = convert_coord(mouse, /device, /to_data)
  area = area[0:1]

  ; get the selected plot type
  plotsel = widget_info(id.selector,/droplist_select)
  case plotsel of
    0: begin
         ; check whether the selection is within the plot area
         if area[0] lt info.xrange[0] then return
         if area[0] gt info.xrange[1] then return
         if area[1] lt info.zrange[0] then return
         if area[1] gt info.zrange[1] then return
         ; get the x-index
         tmp=min(abs(data.x-area[0]),xidx)
         info.xidx=xidx
       end
    1: begin
         ; check whether the selection is within the plot area
         if area[0] lt info.yrange[0] then return
         if area[0] gt info.yrange[1] then return
         if area[1] lt info.zrange[0] then return
         if area[1] gt info.zrange[1] then return
         ; get the x-index
         if size(data.y, /n_dim) le 1 then begin    ; y-dimension is constant
           tmp=min(abs(data.y-area[0]),yidx)
         endif else begin                   ; y-dimension is x-dependent
           tmp=min(abs(data.y[info.xidx,*]-area[0]),yidx)
         endelse
         info.yidx=yidx
       end
    2: begin
         ; check whether the selection is within the plot area
         if area[0] lt info.xrange[0] then return
         if area[0] gt info.xrange[1] then return
         if area[1] lt info.yrange[0] then return
         if area[1] gt info.yrange[1] then return
         ; get the x-index
         tmp=min(abs(data.x-area[0]),xidx)
         if size(data.y, /n_dim) le 1 then begin    ; y-dimension is constant
           tmp=min(abs(data.y-area[1]),yidx)
         endif else begin                   ; y-dimension is x-dependent
           tmp=min(abs(data.y[xidx,*]-area[1]),yidx)
         endelse
         info.xidx=xidx
         info.yidx=yidx
       end
  endcase

  ; save the new range info
  widget_control, id.draw, set_uvalue=info
  ; set the sliders to xidx and yidx
  widget_control, id.xslide, set_value=info.xidx
  widget_control, id.yslide, set_value=info.yidx
end


;******************************************************************************************************
; cw_draw2d_pan ROUTINE
;******************************************************************************************************
pro cw_draw2d_pan, state

  ; get the necessary data
  id = state.id
  if state.ndata eq 0 then return
  data = *state.dataptr[0]

  ; get the draw info structure
  widget_control, id.draw, get_uvalue=info
  ; get the mouse positions
  mouse = [info.mousepos.x0, info.mousepos.y0, info.mousepos.x1, info.mousepos.y1]

  ; convert the mouse device coordinates into data coordinates
  !x.s = info.scale[0:1]
  !y.s = info.scale[2:3]
  !y.type = 0
  area = convert_coord([[mouse[0:1]],[mouse[2:3]]], /device, /to_data)
  xpan = area[0,0]-area[0,1]
  ypan = area[1,0]-area[1,1]

  ; get the selected plot type
  plotsel = widget_info(id.selector,/droplist_select)
  case plotsel of
    0: begin
         ; shift the x- and z-range
         info.xrange+=xpan
         info.zrange+=ypan
       end
    1: begin
         ; shift the x- and z-range
         info.yrange+=xpan
         info.zrange+=ypan
       end
    2: begin
         ; shift the x- and z-range
         info.xrange+=xpan
         info.yrange+=ypan
       end
  endcase

  ; save the new range info
  widget_control, id.draw, set_uvalue=info
  ; and plot
  cw_draw2d_plot, state
end

;******************************************************************************************************
; cw_draw2d_rotate ROUTINE
;******************************************************************************************************
pro cw_draw2d_rotate, state

  ; get the necessary data
  id = state.id
  if state.ndata eq 0 then return
  data = *state.dataptr[0]

  ; get the draw info structure
  widget_control, id.draw, get_uvalue=info
  ; get the mouse positions
  mouse = [info.mousepos.x0, info.mousepos.y0, info.mousepos.x1, info.mousepos.y1]

  ; convert the mouse device coordinates into normalised coordinates
  !x.s = info.scale[0:1]
  !y.s = info.scale[2:3]
  !y.type = 0
  area = convert_coord([[mouse[0:1]], [mouse[2:3]]], /device, /to_normal)
  area = area[0:1,*]
  ;apply rotation
  info.az += 60.*(area[0,1]-area[0,0])
  if info.az gt  90.0 then info.az= 90.0
  if info.az lt   0.0 then info.az=  0.0

  info.ax += 60.*(area[1,0]-area[1,1])
  if info.ax gt 90.0 then info.ax=90.0
  if info.ax lt  0.0 then info.ax= 0.0

  ; save the new range info
  widget_control, id.draw, set_uvalue=info

  ; and replot
  cw_draw2d_plot, state
end

;******************************************************************************************************
; cw_draw2d_zoom ROUTINE
;******************************************************************************************************
pro cw_draw2d_zoom, state, reset=reset

  ; get the necessary data
  id = state.id
  if state.ndata eq 0 then return
  data = *state.dataptr[0]

  ; get the draw info structure
  widget_control, id.draw, get_uvalue=info
  ; get the mouse positions
  mouse = [info.mousepos.x0, info.mousepos.y0, info.mousepos.x1, info.mousepos.y1]

  ; reset zoom to maximum if the reset-keyword is set
  if keyword_set(reset) then begin
    ; reset the x-, y- and z-ranges
    info.xrange = info.fullxrange
    if cw_draw2d_check(data) eq 1 then info.yrange = info.fullyrange
    info.zrange = info.fullzrange
    info.ax = 30.
    info.az = 30.
    widget_control, id.draw, set_uvalue=info
  endif else begin
  ; do the zooming

    ; if the button down position is within a 2 pixel radius
    ; of the button up position (<=> just clicked, not dragged)
    ; then do a standard 2x zoom centered on the clicked point.
    mouse = float(mouse)
    if ((mouse[0]-mouse[2])^2 + (mouse[1]-mouse[3])^2) le 4 then zoom2x=1 else zoom2x=0

    ; convert the mouse device coordinates into normalised coordinates
    !x.s = info.scale[0:1]
    !y.s = info.scale[2:3]
    !y.type = 0
    area = convert_coord([[mouse[0:1]], [mouse[2:3]]], /device, /to_normal)
    area = area[0:1,*]

    ; get the selected plot type
    plotsel = widget_info(id.selector,/droplist_select)
    case plotsel of
    0: begin
         tmp = convert_coord(area, /normal, /to_data)
         if zoom2x then begin
           xrange = 0.5*(info.xrange[1]-info.xrange[0])
           zrange = 0.5*(info.zrange[1]-info.zrange[0])
           xmin   = tmp[0,0]-0.5*xrange
           xmax   = tmp[0,0]+0.5*xrange
           zmin   = tmp[1,0]-0.5*zrange
           zmax   = tmp[1,0]+0.5*zrange
         endif else begin
           xmin = min(tmp[0,*], max=xmax)
           zmin = min(tmp[1,*], max=zmax)
         endelse
         info.xrange = [xmin,xmax]
         info.zrange = [zmin,zmax]
       end
    1: begin
         tmp = convert_coord(area, /normal, /to_data)
         if zoom2x then begin
           yrange = 0.5*(info.yrange[1]-info.yrange[0])
           zrange = 0.5*(info.zrange[1]-info.zrange[0])
           ymin   = tmp[0,0]-0.5*yrange
           ymax   = tmp[0,0]+0.5*yrange
           zmin   = tmp[1,0]-0.5*zrange
           zmax   = tmp[1,0]+0.5*zrange
         endif else begin
           ymin = min(tmp[0,*], max=ymax)
           zmin = min(tmp[1,*], max=zmax)
         endelse
         info.yrange = [ymin,ymax]
         info.zrange = [zmin,zmax]
       end
    2: begin
         if min(area[0,*]) le info.pos[2] then begin  ; zooming the main plot
           tmp = convert_coord(area, /normal, /to_data)
           if zoom2x then begin
             xrange = 0.5*(info.xrange[1]-info.xrange[0])
             yrange = 0.5*(info.yrange[1]-info.yrange[0])
             xmin   = tmp[0,0]-0.5*xrange
             xmax   = tmp[0,0]+0.5*xrange
             ymin   = tmp[1,0]-0.5*yrange
             ymax   = tmp[1,0]+0.5*yrange
           endif else begin
             xmin = min(tmp[0,*], max=xmax)
             ymin = min(tmp[1,*], max=ymax)
           endelse
           info.xrange = [xmin,xmax]
           info.yrange = [ymin,ymax]
         endif else begin                            ; zooming the color bar
           a = (info.zrange[1]-info.zrange[0])/(info.pos[3]-info.pos[1])
           b = info.zrange[0]-a*info.pos[1]
           if zoom2x then begin
             zrange  = 0.5*(info.zrange[1]-info.zrange[0])
             zcentre = a*min(area[1,0])+b
             zmin    = zcentre-0.5*zrange
             zmax    = zcentre+0.5*zrange
           endif else begin
             zmin = a*min(area[1,*])+b
             zmax = a*max(area[1,*])+b
           endelse
           info.zrange = [zmin,zmax]
         endelse
       end
    3: begin
         if min(area[0,*]) le info.pos[2] then return  ; can't zoom 3D surface
         ; zooming the color bar
         a = (info.zrange[1]-info.zrange[0])/(info.pos[3]-info.pos[1])
         b = info.zrange[0]-a*info.pos[1]
         if zoom2x then begin
           zrange  = 0.5*(info.zrange[1]-info.zrange[0])
           zcentre = a*min(area[1,0])+b
           zmin    = zcentre-0.5*zrange
           zmax    = zcentre+0.5*zrange
         endif else begin
           zmin = a*min(area[1,*])+b
           zmax = a*max(area[1,*])+b
         endelse
         info.zrange = [zmin,zmax]
       end
    endcase
  endelse

  ; save the new range info
  widget_control, id.draw, set_uvalue=info

  ; and replot
  cw_draw2d_plot, state

end

;******************************************************************************************************
; SET VALUE ROUTINE
;******************************************************************************************************
pro cw_draw2d_setvalue, base, data

  ; get the state structure from the first child of the top base
  widget_control, widget_info(base,/child), get_uvalue=state
  ; get the widget_hierarchy
  id=state.id
  ; get the draw info structure
  widget_control, id.draw, get_uvalue=info

  ; check is the data is exactly -1 => erase the plot!
  if (size(data,/type) gt 1) && (size(data,/type) lt 6)  && (data[0] eq -1) then begin
    info.xrange = [0.0,1.0]
    info.yrange = [0.0,1.0]
    info.zrange = [0.0,1.0]
    info.ax     = 30.
    info.ax     = 30.
    ; set the labels for 'data(x)' and 'data(y)', select the first plot type ...
    widget_control, id.selector, set_value=['data(x)','data(y)','contour','surface'], set_droplist_select=0
    widget_control, id.selector, sensitive=0
    widget_control, id.contopt.top, sensitive=0
    widget_control, id.stckopt.top, sensitive=0
    ; set the default number of contours
    widget_control, id.contopt.clevels, set_value=string(format='("#",I0)',20)
    ; clear the sticky plots
    sticky={xidx: [-1],$ ; sticky x-indices of the data(y) plots
            yidx: [-1] } ; sticky y-indices of the data(x) plots
    widget_control, id.stckopt.stick, set_uvalue=sticky
    ; set the sliders and the corresponding xidx and yidx in the info structure
    widget_control, id.xslide, set_value=0, set_slider_min=0, set_slider_max=100
    widget_control, id.yslide, set_value=0, set_slider_min=0, set_slider_max=100
    widget_control, id.xslide, sensitive=0
    widget_control, id.yslide, sensitive=0
    info.xidx = 0
    info.yidx = 0
    ; erase the draw widget and set up the color table
    widget_control, id.draw, get_value =win
    wset, win
    erase, truecolor(info.background)
    ; save the draw info
    widget_control, id.draw, set_uvalue=info

    ; free all pointers in the dataptr array
    for i=0,9 do ptr_free, state.dataptr[i]
    state.ndata=0

    ; save the state
    widget_control, widget_info(base,/child), set_uvalue=state
    ; and return
    return
  endif

  ; check if the data has the correct structure
  if ~(cw_draw2d_check(data) eq 0) then begin
    ; find out if this is an overplot or not:
    if cw_draw2d_check(data,istag='overplot') then overplot=abs(fix(data.overplot mod 2)) else overplot=0
    ; if it is: just add a pointer to the data in the dataptr-arrar and increase the ndata field
    if overplot then begin
      if state.ndata eq 10 then begin
        tmp=dialog_message("Cannot overplot: maximum number plots that can be overplot is 10!",$
                           title="max. number of overplots",dialog_parent=id.top)
        return
      endif
      state.dataptr[state.ndata] = ptr_new(data)
      state.ndata+=1

      ; also set the foreground and background colors to that of the 'master' plot data
      data0 = *state.dataptr[0]
      if cw_draw2d_check(data0,istag='fcolor') then !p.color     =truecolor(data0.fcolor) else !p.color     =truecolor('black')
      if cw_draw2d_check(data0,istag='bcolor') then !p.background=truecolor(data0.bcolor) else !p.background=truecolor('white')

    ; if it isn't: remove all existing plot data, set the first element of dataptr to be a pointer
    ;              to the data, set ndata to 1, update ranges, plot type, .... et cetera
    endif else begin

      ; check whether the data is 2D or 1D
      if (cw_draw2d_check(data) eq 1) then data2D=1 else data2d=0

      ; if idxonly is selected skip the range settings, don't clear the sticky settings 
      ; ... just adjust the xidx and yidx and the sliders
      if ~cw_draw2d_check(data,istag='idxonly') || (data.idxonly eq 0) then begin
        ; set and check the x-range
        if cw_draw2d_check(data,istag='xrange') then xrange= data.xrange else begin
          xrange = float([min(data.x), max(data.x)])
          extend = 0.02*[xrange[1]-xrange[0]]
          xrange = [xrange[0]-extend,xrange[1]+extend]
        endelse
        if xrange[0] gt xrange[1] then begin
          print, ' ERROR: xrange[0] > xrange[1] (xrange changed to [min(xrange), max(xrange)])'
          xrange = [min(xrange), max(xrange)]
        endif
        if xrange[0] eq xrange[1] then begin
          print, ' ERROR: xrange[0] = xrange[1] (xrange changed to [xrange[0]-1, xrange[0]+1])'
          xrange = [xrange[0]-1., xrange[0]+1.]
        endif
        if xrange[0] ge max(data.x) then begin
          print, ' ERROR: xrange[0] > max(x) (xrange changed to [min(x), max(x)])'
          xrange = [min(data.x), max(data.x)]
        endif
        if xrange[1] le min(data.x) then begin
          print, ' ERROR: xrange[1] < min(x) (xrange changed to [min(x), max(x)])'
          xrange = [min(data.x), max(data.x)]
        endif
        info.xrange=xrange
        info.fullxrange=xrange

        ; set and check the y-range
        if data2D then begin
          if cw_draw2d_check(data,istag='yrange') then yrange= data.yrange else begin
          yrange = float([min(data.y), max(data.y)])
          extend = 0.02*[yrange[1]-yrange[0]]
          yrange = [yrange[0]-extend,yrange[1]+extend]
        endelse
          if yrange[0] gt yrange[1] then begin
            print, ' ERROR: yrange[0] > yrange[1] (yrange changed to [min(yrange), max(yrange)])'
            yrange = [min(yrange), max(yrange)]
          endif
          if yrange[0] eq yrange[1] then begin
            print, ' ERROR: yrange[0] = yrange[1] (yrange changed to [yrange[0]-1, yrange[0]+1])'
            yrange = [yrange[0]-1., yrange[0]+1.]
          endif
          if yrange[0] ge max(data.y) then begin
            print, ' ERROR: yrange[0] > max(y) (yrange changed to [min(y), max(y)])'
            yrange = [min(data.y), max(data.y)]
          endif
          if yrange[1] le min(data.y) then begin
            print, ' ERROR: yrange[1] < min(y) (yrange changed to [min(y), max(y)])'
            yrange = [min(data.y), max(data.y)]
          endif
        endif else yrange=[0.,1.]
        info.yrange=yrange
        info.fullyrange=yrange

        ; set and check the z-range
        if cw_draw2d_check(data,istag='zrange') then zrange= data.zrange else begin
          fidx = where(finite(data.z),cnt)
          if cnt eq 0 then begin
            ; only NaN data <=> invalid data
            print, 'ERROR: z contains no finite values!'
            return
          endif else begin
            zmin   = min(data.z[fidx], max=zmax)
            zrange = zmax-zmin
            zrange = float([zmin-0.05*zrange, zmax+0.05*zrange])
          endelse
        endelse
        if zrange[0] gt zrange[1] then begin
          print, ' ERROR: zrange[0] > zrange[1] (zrange changed to [min(zrange), max(zrange)])'
          zrange = [min(zrange), max(zrange)]
        endif
        if zrange[0] eq zrange[1] then begin
          print, ' ERROR: zrange[0] = zrange[1] (zrange changed to [zrange[0]-1, zrange[0]+1])'
          zrange = [zrange[0]-1., zrange[0]+1.]
        endif
        if zrange[0] ge max(data.z) then begin
          print, ' ERROR: zrange[0] > max(z) (zrange changed to [min(z), max(z)])'
          fidx = where(finite(data.z),cnt)
          if cnt eq 0 then info.zrange = [0.0,1.0] else begin
            zmin   = min(data.z[fidx], max=zmax)
            zrange = zmax-zmin
            zrange = float([zmin-0.05*zrange, zmax+0.05*zrange])
          endelse
        endif
        if zrange[1] le min(data.z) then begin
          print, ' ERROR: zrange[1] < min(z) (zrange changed to [min(z), max(z)])'
          fidx = where(finite(data.z),cnt)
          if cnt eq 0 then info.zrange = [0.0,1.0] else begin
            zmin   = min(data.z[fidx], max=zmax)
            zrange = zmax-zmin
            zrange = float([zmin-0.05*zrange, zmax+0.05*zrange])
          endelse
        endif
        info.zrange=zrange
        info.fullzrange=zrange

        ; set the foreground and background colors
        if cw_draw2d_check(data,istag='fcolor') then !p.color     =truecolor(data.fcolor) else !p.color     =truecolor('black')
        if cw_draw2d_check(data,istag='bcolor') then !p.background=truecolor(data.bcolor) else !p.background=truecolor('white')

        ; set the labels for 'data(x)' and 'data(y)', select the first plot type ...
        if cw_draw2d_check(data,istag='xsellbl') then  xlbl=data.xsellbl else xlbl='data(x)'
        if cw_draw2d_check(data,istag='ysellbl') then  ylbl=data.ysellbl else ylbl='data(y)'
        if cw_draw2d_check(data,istag='ptype') && data2D then ptype=abs(data.ptype mod 4)  else ptype=0
        widget_control, id.selector, set_value=[xlbl,ylbl,'contour','surface'], set_droplist_select=ptype
        widget_control, id.selector, sensitive=data2D?1:0
        if ptype le 1 then begin
          widget_control, id.contopt.top, sensitive=0
          widget_control, id.stckopt.top, sensitive=data2D?1:0
        endif else begin
          widget_control, id.contopt.top, sensitive=1
          widget_control, id.stckopt.top, sensitive=0
        endelse

        ; set the default number of contours
        if cw_draw2d_check(data,istag='clevels') then clevels=data.clevels else clevels='#20'
        widget_control, id.contopt.clevels, set_value=clevels
        if cw_draw2d_check(data,istag='clabels') then clabels=abs(fix(data.clabels mod 2)) else clabels=0
        widget_control, id.contopt.clabels, set_value=clabels

        ; clear the sticky plots
        sticky={xidx: [-1],$ ; sticky x-indices of the data(y) plots
                yidx: [-1] } ; sticky y-indices of the data(x) plots
        widget_control, id.stckopt.stick, set_uvalue=sticky
      endif

      ; set the sliders and the corresponding xidx and yidx in the info structure
      if cw_draw2d_check(data, istag='xidx') then xidx=data.xidx else xidx=0
      if cw_draw2d_check(data, istag='yidx') && data2D then yidx=data.yidx else yidx=0
      if xidx ge n_elements(data.x) then xidx=0
      if data2D && (yidx ge n_elements(data.y)) then yidx=0
      widget_control, id.xslide, set_value=xidx, set_slider_min=0, set_slider_max=n_elements(data.x)-1
      if data2D then if size(data.y,/n_dim) le 1 then n_ch=n_elements(data.y)-1 else n_ch=n_elements(data.y[0,*])-1
      widget_control, id.yslide, set_value=yidx, set_slider_min=0, set_slider_max=data2D?n_ch:100
      widget_control, id.xslide, sensitive=1
      widget_control, id.yslide, sensitive=data2D?1:0
      info.xidx = xidx
      info.yidx = yidx

      ; save the draw info
      widget_control, id.draw, set_uvalue=info

      ; free all pointers in the dataptr array
      for i=0,9 do ptr_free, state.dataptr[i]
      ; add a pointer to the data
      state.dataptr[0] = ptr_new(data)
      state.ndata=1

    endelse

    ; save the state
    widget_control, widget_info(base,/child), set_uvalue=state

    ; make the plot
    cw_draw2d_plot, state

  endif

end


;******************************************************************************************************
; GET VALUE FUNCTION
;******************************************************************************************************
function cw_draw2d_getvalue, base, data
  ; get widget hierarchy
  widget_control, widget_info(base,/child), get_uvalue=state
  id=state.id
  ; and return the window id of for the drawing area
  widget_control, id.draw, get_value=win
  return, win
end


;******************************************************************************************************
; KILL ROUTINE => runs when the cw_draw2d widget is being destroyed
;******************************************************************************************************
pro cw_draw2d_kill, base
  ; get widget state structure
  widget_control, base, get_uvalue=state
  ; free all pointers in the dataptr array
  for i=0,9 do ptr_free, state.dataptr[i]
end


;******************************************************************************************************
; ZOOM/SELECT/PAN ROUTINE (handles mouse-clicks and drags and converts them in zoom and selections)
;******************************************************************************************************
function cw_draw2d_mouse_event, event, state
  ; returns 0 when no external event needs to be returned (zoom, pan, rotation),
  ; 1 if the event need to be returned (selection)

  ; get the widget hierarchy
  id = state.id

  ; get the draw widget info (e.g. mouse positions)
  widget_control, id.draw, get_uvalue=info
  ; give graphics control to the current draw window
  widget_control, id.draw, get_value=win
  wset, win

  ; flag that determines whether we should report this event to the outside world
  return_event=0

  ; get the event type:
  ; mouse down
  if event.type eq 0 then begin
    ; left button (rotate in surface plot)
    if event.press eq 1 then begin
      ; save the current mouse position
      info.mousepos.x0 = event.x
      info.mousepos.y0 = event.y
      info.mousepos.x1 = event.x
      info.mousepos.y1 = event.y
      ; save the info
      widget_control, id.draw, set_uvalue=info
    endif

    ; middle button (zoom)
    if event.press eq 2 then begin
      ; save the current mouse position
      info.mousepos.x0 = event.x
      info.mousepos.y0 = event.y
      info.mousepos.x1 = event.x
      info.mousepos.y1 = event.y
      ; activate the motion events
      widget_control, id.draw, draw_motion_events=1
      ; and allow overplotting with the inverse color (no replot of the whole plot needed)
      device, get_graphics=oldgraphics, set_graphics=6
      info.graphicsmode=oldgraphics
      ; save the info
      widget_control, id.draw, set_uvalue=info
    endif

  endif

  ; mouse up
  if event.type eq 1 then begin
    ; left button (rotate in surface plot, select/pan in all other plots)
    if event.release eq 1 then begin
     ; get the index of the plot selector: 
      select= widget_info(id.selector,/droplist_select)
      if select eq 3 then begin ; in case of a surface plot
        ; save the current mouse position
        info.mousepos.x1 = event.x
        info.mousepos.y1 = event.y
        ; save the info
        widget_control, id.draw, set_uvalue=info
        ; call the rotate routine
        cw_draw2d_rotate, state
      endif else begin          ; in case of data(x), data(y) and contour plot
        ; save the current mouse position
        info.mousepos.x1 = event.x
        info.mousepos.y1 = event.y
        ; save the info
        widget_control, id.draw, set_uvalue=info
        ; if the 'button up' coordinates are within a 2 pixel area of the
        ; 'button down'coordinates: call the select routine
        if abs( (float(info.mousepos.x0)-float(info.mousepos.x1))^2 $
               +(float(info.mousepos.y0)-float(info.mousepos.y1))^2 ) lt 4 then begin
          cw_draw2d_select, state
          ; set this event to be returned
          return_event=1
        endif else begin
        ; if the 'button up' coordinates differ from the
        ; 'button down'coordinates: call the pan routine
          cw_draw2d_pan, state
        endelse
      endelse
    endif

    ; middle button (zoom)
    if event.release eq 2 then begin
      ; deactivate the motion events
      widget_control, id.draw, draw_motion_events=0
      ; hide the zoom selection box
      plots, /device, [info.mousepos.x0,info.mousepos.x1],$
                      [info.mousepos.y0,info.mousepos.y0],$
                      thick=2, col='ffffff'x
      plots, /device, [info.mousepos.x0,info.mousepos.x1],$
                      [info.mousepos.y1,info.mousepos.y1],$
                      thick=2, col='ffffff'x
      plots, /device, [info.mousepos.x0,info.mousepos.x0],$
                      [info.mousepos.y0,info.mousepos.y1],$
                      thick=2, col='ffffff'x
      plots, /device, [info.mousepos.x1,info.mousepos.x1],$
                      [info.mousepos.y0,info.mousepos.y1],$
                      thick=2, col='ffffff'x
      ; save current mouse position
      info.mousepos.x1 = event.x
      info.mousepos.y1 = event.y
      ; save the info
      widget_control, id.draw, set_uvalue=info
      ; reset the graphics mode
      device, set_graphics=info.graphicsmode
      ; call the zoom routine
      cw_draw2d_zoom, state
    endif

    ; right button (unzoom)
    if event.release eq 4 then begin
      ; call the unzoom routine
      cw_draw2d_zoom, state, /reset
    endif
  endif

  ; mouse motion
  if event.type eq 2 then begin
    ; hide the previous zoom selection box
    plots, /device, [info.mousepos.x0,info.mousepos.x1],$
                    [info.mousepos.y0,info.mousepos.y0],$
                    thick=2, col='ffffff'x
    plots, /device, [info.mousepos.x0,info.mousepos.x1],$
                    [info.mousepos.y1,info.mousepos.y1],$
                    thick=2, col='ffffff'x
    plots, /device, [info.mousepos.x0,info.mousepos.x0],$
                    [info.mousepos.y0,info.mousepos.y1],$
                    thick=2, col='ffffff'x
    plots, /device, [info.mousepos.x1,info.mousepos.x1],$
                    [info.mousepos.y0,info.mousepos.y1],$
                    thick=2, col='ffffff'x
    ; save current mouse position
    info.mousepos.x1 = event.x
    info.mousepos.y1 = event.y
    ; plot the new zoom selection box
    plots, /device, [info.mousepos.x0,info.mousepos.x1],$
                    [info.mousepos.y0,info.mousepos.y0],$
                    thick=2, col='ffffff'x
    plots, /device, [info.mousepos.x0,info.mousepos.x1],$
                    [info.mousepos.y1,info.mousepos.y1],$
                    thick=2, col='ffffff'x
    plots, /device, [info.mousepos.x0,info.mousepos.x0],$
                    [info.mousepos.y0,info.mousepos.y1],$
                    thick=2, col='ffffff'x
    plots, /device, [info.mousepos.x1,info.mousepos.x1],$
                    [info.mousepos.y0,info.mousepos.y1],$
                    thick=2, col='ffffff'x
    ; save the info
    widget_control, id.draw, set_uvalue=info
  endif

  ; return whether an event if requested
  return, return_event

end


;******************************************************************************************************
; cw_draw2d_save_event EVENT ROUTINE
;******************************************************************************************************
pro cw_draw2d_save_event, state

  ; get the widget hierarchy
  id = state.id

  ; check the data
  if state.ndata eq 0 then return
  data = *state.dataptr[0]

  ; Open 'save as ...' dialog
open_save_dialog:
  result = save_dialog(dialog_parent=id.top, filter=['.eps','.png'],file='figure', title='Save figure as ...')
  if strcmp(result[0],'') then return
  ; check if the file already exists
  if file_test(result[0]) then begin
    answer = dialog_message([string(format='("File <",A,"> already exists!")',result[0]),$
                             'Do you want to overwrite it?'], dialog_parent=id.top, /question,$
                            title='File exists')
    if strcmp(answer,'no',/fold_case) then goto, open_save_dialog
  endif

  ; get the settings of the draw widget
  widget_control, id.draw, get_uvalue=info

  ; save as EPS figure
  if strcmp(result[1],'.eps',/fold_case) then begin
    ; set the aspect ratio and the font size
    aspect   = float(info.xsize)/float(info.ysize)
    if cw_draw2d_check(data,istag='charsize') then fontsize = 4.4*data.charsize else fontsize=4.4

    ; open the eps-file for plotting
    sys = init_graphic(aspect=aspect, fontsize=fontsize, eps=result[0], bits_per_pixel=8)

    ; get the selected plot type
    plotsel = widget_info(id.selector,/droplist_select)
    case plotsel of
     0: begin
          cw_draw2d_datax, state
        end
     1: begin
          cw_draw2d_datay, state
        end
     2: begin
          cw_draw2d_contour, state
        end
     3: begin
          cw_draw2d_surf, state
        end
    endcase
    ; close the eps-file
    tmp = reset_graphic(sys, eps=result[0])
  endif

  ; save as PNG figure
  if strcmp(result[1],'.png',/fold_case) then begin
    ; get the current window
    widget_control, id.draw, get_value=win
    wset, win
    ; read in the X-buffer
    img = tvrd(/true)
    ; and write it to the file
    write_png,result[0],img
  endif


end


;******************************************************************************************************
; cw_draw2d_setcolor_event EVENT ROUTINE
;******************************************************************************************************
pro cw_draw2d_setcolor_event, state
  ; get the widget hierarchy
  id = state.id
  ; get the draw info structure
  widget_control, id.draw, get_uvalue=info
  ; get the color table
  ctable = widget_info(id.contopt.color,/droplist_select)-1
  ; set the color table
  info.ctable = ctable
  widget_control, id.draw, set_uvalue=info
end


;******************************************************************************************************
; cw_draw2d_plottype_event EVENT ROUTINE
;******************************************************************************************************
pro cw_draw2d_plottype_event, state
  ; get the widget hierarchy
  id = state.id
  ; get the selected plot type and make the relevant controls active
  ptype = widget_info(id.selector,/droplist_select)
  if ptype le 1 then begin
    widget_control, id.contopt.top, sensitive=0
    widget_control, id.stckopt.top, sensitive=1
  endif else begin
    widget_control, id.contopt.top, sensitive=1
    widget_control, id.stckopt.top, sensitive=0
  endelse
end

;******************************************************************************************************
; cw_draw2d_stick_event EVENT ROUTINE
;******************************************************************************************************
pro cw_draw2d_stick_event, state
  ; get the widget hierarchy
  id = state.id
  ; get the draw info structure
  widget_control, id.draw, get_uvalue=info
  ; get the sticky structure
  widget_control, id.stckopt.stick, get_uvalue=sticky

  ; get the selected plot type and make current plot sticky
  plotsel = widget_info(id.selector,/droplist_select)
  case plotsel of
  0: begin
       if sticky.yidx[0] eq -1 then yidx=info.yidx else begin
         if n_elements(sticky.yidx) eq 9 then begin
           tmp=dialog_message("Cannot stick this plot: maximum number of sticky plots is 9!",$
                              title="max. number of sticky plots",dialog_parent=id.top)
           return
         endif
         tmp = where(sticky.yidx eq info.yidx, exist)
         if exist then return
         yidx=[sticky.yidx, info.yidx]
       endelse
       stickynew = {xidx:sticky.xidx, yidx:yidx}
     end
  1: begin
       if sticky.xidx[0] eq -1 then xidx=info.xidx else  begin
          if n_elements(sticky.xidx) eq 9 then begin
           tmp=dialog_message("Cannot stick this plot: maximum number of sticky plots is 9!",$
                              title="max. number of sticky plots",dialog_parent=id.top)
           return
         endif
         tmp = where(sticky.xidx eq info.xidx, exist)
         if exist then return
         xidx=[sticky.xidx, info.xidx]
       endelse
       stickynew = {xidx:xidx, yidx:sticky.yidx}
     end
  else: return
  endcase
  ; save the new sticky structure
  widget_control, id.stckopt.stick, set_uvalue=stickynew
end


;******************************************************************************************************
; cw_draw2d_unstick_event EVENT ROUTINE
;******************************************************************************************************
pro cw_draw2d_unstick_event, state
  ; get the widget hierarchy
  id = state.id
  ; get the draw info structure
  widget_control, id.draw, get_uvalue=info
  ; get the sticky structure
  widget_control, id.stckopt.stick, get_uvalue=sticky

  ; get the selected plot type and clear the current sticky plots
  plotsel = widget_info(id.selector,/droplist_select)
  case plotsel of
  0: begin
       if sticky.yidx[0] eq -1 then return
       stickynew = {xidx:sticky.xidx, yidx:[-1]}
     end
  1: begin
       if sticky.xidx[0] eq -1 then return
       stickynew = {xidx:[-1], yidx:sticky.yidx}
     end
  else: return
  endcase
  ; save the new sticky structure
  widget_control, id.stckopt.stick, set_uvalue=stickynew
end


;******************************************************************************************************
; GENERAL EVENT FUNCTION
;******************************************************************************************************
function cw_draw2d_event, event

  ; get the state
  widget_control, widget_info(event.handler,/child), get_uvalue=state
  ; widget hierarchy
  id = state.id

  ; get the type of event
  case event.id of
  id.selector       : begin & return_event=1 & replot=1 & cw_draw2d_plottype_event, state & end
  id.crosshair      : begin & return_event=0 & replot=1 & end
  id.contopt.color  : begin & return_event=0 & replot=1 & cw_draw2d_setcolor_event, state & end
  id.contopt.clevels: begin & return_event=0 & replot=1 & end
  id.contopt.clabels: begin & return_event=0 & replot=1 & end
  id.stckopt.stick  : begin & return_event=0 & replot=1 & cw_draw2d_stick_event   , state & end
  id.stckopt.clear  : begin & return_event=0 & replot=1 & cw_draw2d_unstick_event , state & end
  id.save           : begin & return_event=0 & replot=0 & cw_draw2d_save_event    , state & end
  id.draw           : begin & return_event=cw_draw2d_mouse_event(event, state) & replot=2 & end
  id.xslide         : begin & return_event=1 & replot=2 & end
  id.yslide         : begin & return_event=1 & replot=2 & end
  else              : begin & return_event=0 & replot=0 & end
  endcase

  ; get the the position of the x and y sliders
  widget_control, id.xslide, get_value=xidx
  widget_control, id.yslide, get_value=yidx
  ; set the xidx and yidx values to the draw-info structure
  widget_control, id.draw, get_uvalue=info
  info.xidx=xidx
  info.yidx=yidx
  widget_control, id.draw, set_uvalue=info

  ; get the plot type
  ptype = widget_info(id.selector,/droplist_select)

  ; replot if requested
  ; (0: no replot, 1: replot, 2: replot, but not when return_event is 1
  ;                              and a contour or surface plot is selected without crosshair)
  case replot of
  0: ; do nothing
  1: cw_draw2d_plot, state
  2: begin
       if ptype le 1 then begin
         if return_event then cw_draw2d_plot, state
       endif else begin
         widget_control, id.crosshair, get_value=crossh
         if crossh[0] && return_event then cw_draw2d_plot, state
       endelse
     end
  endcase

  ; return the event if requested
  if return_event then begin
    ; creat the new event structure
    ret = {ID:      Event.Handler,$
           TOP:     Event.Top,    $
           HANDLER: 0l,           $
           xidx:    xidx,         $
           yidx:    yidx,         $
           ptype:   ptype         }

    if state.eventpro eq '' then begin
      return, ret               ; call any event function or procedure higher up
    endif else begin
      call_procedure, state.eventpro, ret ; call this specific routine
      return,1
    endelse
  endif else return,1
end


;******************************************************************************************************
; MAIN FUNCTION (draws up the widget)
;******************************************************************************************************
function cw_draw2d, base, xsize=xsize, ysize=ysize, frame=frame, event_pro=epro, background=background,$
                          crossh=crossh, ctable=ctable

  ; set a default size if the xsize and ysize keywords are not present
  if ~keyword_set(xsize) then xsize=640
  if ~keyword_set(ysize) then ysize=480

  ; check whether the frame is present. If not: set frame to 1
  if n_elements(frame) eq 0 then frame=1


  ; set up an internal id-structure for all the controls within the widget
  id={top:0l, draw:0l, selector:0l, crosshair:0l,$
      contopt:{top:0l, color:0l, clabels:0l,     $
               lbllevel:0l, clevels:0l},         $
      stckopt:{top:0l, stick:0l, clear:0l},      $
      save:0l, xslide:0l, yslide:0l              }

  ; top row with plot selector, contour/surface options, stick/clear buttons and 'save to eps' button
  id.top=widget_base(base,/column, frame=frame,$
                    event_func='cw_draw2d_event',$
                    pro_set_value='cw_draw2d_setvalue',$
                    func_get_value='cw_draw2d_getvalue')

  ; the top row: this contains the general controls and its uvalue will alse contain the state
  ; of the widget in it's uvalue (hence we need to know when it's killed) 
  toprow              = widget_base(id.top,/row, kill_notify='cw_draw2d_kill')
  lbl                 = widget_label(toprow, value='plot type:', /align_center)
  id.selector         = widget_droplist(toprow, value=['data(x)','data(y)','contour','surface'],xsize=150, /align_center)
  id.crosshair        = cw_bgroup(toprow, 'crosshair', set_value=keyword_set(crossh),/nonexclusive)
  id.contopt.top      = widget_base(toprow,/row,frame=1)
  loadct, get_names=ctables
  ctables = ['no shade',ctables]
  id.contopt.color    = widget_droplist(id.contopt.top, value=strlowcase(ctables))
  id.contopt.lbllevel = widget_label(id.contopt.top, value='levels:', /align_left)
  id.contopt.clevels  = widget_text(id.contopt.top, value='#20', /editable,xsize=10)
  id.contopt.clabels  = cw_bgroup(id.contopt.top, 'labels', set_value=0,/nonexclusive)
  id.stckopt.top      = widget_base(toprow,/row,frame=1)
  id.stckopt.stick    = widget_button(id.stckopt.top, value='stick', /align_center)
  id.stckopt.clear    = widget_button(id.stckopt.top, value='clear', /align_center)
  spacer              = widget_base(toprow,/col, xsize=15)
  id.save             = widget_button(toprow, value='save...',/align_center, ysize=30,$
                        font='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1')

  ; make the contour/surface options insensitive
  widget_control, id.contopt.top, sensitive=0

  ; main bulk of the widget: plot area and sliders
  drawrow   = widget_base(id.top,/row)
  drawcol   = widget_base(drawrow,/col)
  id.yslide = widget_slider(drawcol,/vertical,minimum=0,maximum=100,ysize=ysize,/drag,scroll=1, /suppress_value)
  drawcol   = widget_base(drawrow,/col)
  id.draw   = widget_draw(drawcol,xsize=xsize,ysize=ysize, /button_events)
  id.xslide = widget_slider(drawcol,minimum=0,maximum=100,xsize=xsize,/drag,scroll=1, /suppress_value)

  ; define data(x) and data(y) sticky plots (none to start with)
  sticky={xidx: [-1],$ ; sticky x-indices of the data(y) plots
          yidx: [-1] } ; sticky y-indices of the data(x) plots
  widget_control, id.stckopt.stick, set_uvalue=sticky

  ; get the default color table for the widget
  if n_elements(ctable) eq 0 then ctable=5l
  widget_control, id.contopt.color, set_droplist_select=ctable+1

  ; get the default background for the widget
  if ~keyword_set(background) || ~is_string(background) then background='white'
  ; set the uvalue of the draw widget <=> graphics mode, mouse positions for zooming and selecting
  draw_info = {graphicsmode: 0l,                     $
               xsize       : xsize,                  $
               ysize       : ysize,                  $
               background  : background,             $
               ctable      : ctable,                 $
               mousepos    :{x0:0, y0:0, x1:0, y1:0},$
               xrange      :[0.0,1.0],               $
               yrange      :[0.0,1.0],               $
               zrange      :[0.0,1.0],               $
               fullxrange  :[0.0,1.0],               $
               fullyrange  :[0.0,1.0],               $
               fullzrange  :[0.0,1.0],               $
               scale       :[1.0,0.0,1.0,0.0],       $
               pos         :[0.0,0.0,1.0,1.0],       $
               ax          : 30.,                    $
               az          : 30.,                    $
               xidx        : 0l,                     $
               yidx        : 0l                      }
  widget_control, id.draw, set_uvalue=draw_info

  ; make an array to hold the pointers to the plot-data structures
  ndata   = 0                   ; number of data structure (max. 10)
  dataptr = ptrarr(10, /nozero) ; 10 unallocated null pointers

  ; set the name of the event_procedure
  if (n_elements(epro) eq 0) then eventpro='' else eventpro=epro

  ; create the state structure
  state = {id      : id     ,$
           ndata   : ndata  ,$
           dataptr : dataptr,$
           eventpro: eventpro}
  ; and set it to the uvalue of the first child of the top base
  ; (i.e. the top row base: we're not using that for anything else anyway)
  widget_control, widget_info(id.top, /child), set_uvalue=state

  ; return the id
  return,id.top
end


