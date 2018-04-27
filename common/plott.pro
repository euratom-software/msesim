;+
; plott,  x, y , z, error=error, xrange=xrange, yrange=yrange, zrange=zrange,
;                   linestyle=linestyle, thick=thick, color=color, psym=psym, symsize=symsize,
;                   errstyle=errstyle, errcolor=errcolor, fcolor=fcolor, bcolor=bcolor,
;                   nlevels=nlevels, clabels=clabels, ctable=ctalble,
;                   xgrid=xgrid, ygrid=ygrid, zgrid=zgrid,
;                   xtitle=xtitle, ytitle=ytitle, ztitle=ztitle, title=title,
;                   xlbl=xlbl, ylbl=ylbl, charsize=charsize, lgdpos=lgdpos, extralgd=extralgd,
;                   crossh=crossh, xidx=xidx, yidx=yidx, ptype=ptype,
;                   xsize=xsize, ysize=ysize,  id=id, overplot=overplot, quiet=quiet,
;                   no_block=no_block, hold=hold
;
; a general plot routine based on the cw_draw2d widget:
;  * allows 1D and 2D plots (data(x) with error bars, data(y) with error bars, contour and surface
;  * allows zooming and panning
;  * allows saving as eps- and png-files
;  * selecting in the plot results in a print of the values on the screen and creating a
;    variable plott_var in the main routine
;
; :Params:
;    x          : in, required, type=any numeric, 1D array
;                [mx1]-array containing the x-dimension (e.g. time)
;   y           : in, required, type=any numeric, 1D or 2D array
;                For 1D data (no "z" input parameter defined)
;                    [mx1]-array containing the data to be plotted
;                For 2D data ("z" input parameter is present)
;                    [nx1]-array containing the y-dimension (e.g. radius)
;                or [mxn]-array containing the y-dimension (e.g. time-dependent radius)
;   z           : in, optional, type=any numeric, 2D array
;                [mxn]-array containing the data to be plotted
; :Keywords:
;   error       : in, optional, type=any numeric, 1D or 2D array
;                Array containing the error on the data. It should have the same size
;                as the data array (i.e. "y" for 1D data, "z" for 2D data)
;   xrange      : in, optional, type=any numeric, 1D array
;                Sets the range of the x-axis
;   yrange      : in, optional, type=any numeric, 1D array
;                Sets the range of the y-axis
;   zrange      : in, optional, type=any numeric, 1D array
;                Sets the range of the z-axis
;   linestyle   : in, optional, type=byte, scalar
;                Sets the linestyle (just like the linestyle-keyword in "plot")
;   thick       : in, optional, type=float, scalar
;                Sets thickness of the line (just like the thick-keyword in "plot")
;                (crosshair and gridlines always they have the standard thickness from !p.thick)
;   color       : in, optional, type=string or integer, scalar or 3-element array
;                Sets the color of the line. WARNING: this is not the same as the color-keyword in "plot"!
;                It will accept either a color name as a string (e.g. 'red') or a 3-element byte array
;                (e.g. [r,g,b]) or an (long) integer with the truecolor index (e.g. '0000ff'x). This keyword
;                will also only set the color of the line! It will not set the color of the axis!
;   psym        : in, optional, type=byte, scalar
;                Sets the plot symbol for the lines (just like the psym-keyword in "plot")
;   symsize     : in, optional, type=float, scalar
;                Sets the size of the plot symbols (just like the symsize-keyword in "plot")
;   errstyle    : in, optional, type=byte, scalar
;                Sets the error style for the data(x) and data(y) plots
;                    0=shaded area, 1=error bars [default=1]
;   errcolor    : in, optional, type=string or integer, scalar or 3-element array
;                Sets the the color of the shaded error area if errstyle=0 [default='gray']
;   fcolor      : in, optional, type=string or integer, scalar or 3-element array
;                Sets the foreground color. This is the color of the axis and the labels [default='black']
;   bcolor      : in, optional, type=string or integer, scalar or 3-element array
;                Sets the background color [default='white']
;   nlevels     : in, optional, type=byte, scalar
;                Number of levels used in the contour plot [default=20, max=100]
;   clabels     : in, optional, type=byte, scalar
;                Switches on contour labels
;   ctable      : in, optional, type=byte, scalar
;                Color table to use for the contour and surface plots [default=5: black-blue-red-yellow-white]
;   xgrid       : in, optional, type=float, 1D array
;                x-coordinates for x-gridlines to be plotted
;   ygrid       : in, optional, type=float, 1D array
;                y-coordinates for y-gridlines to be plotted
;   zgrid       : in, optional, type=float, 1D array
;                z-coordinates for z-gridlines to be plotted
;   xtitle      : in, optional, type=string, scalar
;                Label of the x-axis (uses TeXtoIDL so you can use (some) LaTeX commands. Make sure TeXtoIDL is in your path)
;   ytitle      : in, optional, type=string, scalar
;                Label of the y-axis (uses TeXtoIDL so you can use (some) LaTeX commands. Make sure TeXtoIDL is in your path)
;   ztitle      : in, optional, type=string, scalar
;                Label of the z-axis (uses TeXtoIDL so you can use (some) LaTeX commands. Make sure TeXtoIDL is in your path)
;   xtitle      : in, optional, type=string, scalar
;                Title above the plot (uses TeXtoIDL so you can use (some) LaTeX commands. Make sure TeXtoIDL is in your path)
;   xlbl        : in, optional, type=string, scalar
;                String that labels the 'data as function of x' selection button. if omitted: xlbl = 'data('+xtitle+')'
;   ylbl        : in, optional, type=string, scalar
;                String that labels the 'data as function of y' selection button. if omitted: ylbl = 'data('+ytitle+')'
;   charsize    : in, optional, type=float, scalar
;                Character size for all axes and legends [default=1] 
;   lgdpos      : in, optional, type=byte, scalar
;                Sets the position of the automatically generated legend.
;                0=top-left, 1=top-right, 2=bottom-right, 3=bottom-left [default=1]
;   extralgd    : in, optional, type=structure, scalar
;                Settings for an extra legend. This is NOT compatible with the overplot keywordsetting
;                (Hence create this legend in the first plot before you start overplotting). The settings
;                for the extra legend are given in the fields of the extralgd-structure
;                    .text      = string array with the legend text
;                    .color     = string array with the line colors for the legend
;                    .textcolor = string array with the text colors for the legend
;                    .linestyle = byte array with the linestyles for the legend
;                    .psym      = byte array with the symbols for the legend
;                    .symsize   = float array with the symbol size for the legend
;                    .thick     = float array with thickness of the legend
;                    .box       = byte scalar that draws a box around the legend
;                    .clear     = byte scalar that clears the area of the legend
;                    .position  = position of the legend. 0=top-left, 1=top-right,
;                                 2=bottom-right, 3=bottom-left
;                    .horizontal= set a horizontal legend (columns) rather than a vertical (rows)
;                    .charsize  = character size
;                    .charthick = character thickness
;                    .spacing   = line spacing of the legend
;  crossh       : in, optional, type=byte, scalar
;                Switches on the crosshair [default=0]
;  xidx         : in, optional, type=integer/long, scalar
;                Initial x-selection [default=0]
;  yidx         : in, optional, type=integer/long, scalar
;                Initial y-selection [default=0]
;  ptype        : in, optional, type=byte, scalar
;                Selects the plot type (0:data(x), 1:data(y), 2:contour, 3:surface) [default=0]
;  xsize        : in, optional, type=integer, scalar
;                xsize of the draw area in device coordinates (default=640)
;  xsize        : in, optional, type=integer, scalar
;                ysize of the draw area in device coordinates (default=480)
;  id           : out, optional, type=long, scalar
;                Returns the id of the plott-widget
;  overplot     : in, optional, type=long, scalar
;                If set to a known id for an existing plott-widget, the data will be overplotted onto the data already in
;                that widget. If no original data is present this field is ignored. This only has an effect on the data(x)
;                and data(y)-plots, contour and surface plots will only show the original data. Also the data printed
;                to the screen and saved int he plott_var variable comes from the the original data.
;  quiet        : in, optional, type=byte, scalar
;                If set no data will be printed to the screen and the plott_var variable will not be created
;  block        : in, optional, type=byte, scalar
;                If set, IDL will hold until the plott-widget is closed. If the block-keyword is set, the data will
;                still be printed to the screen, but the plott_var variable will no longer be created!
;  hold         : out, optional, type=byte, scalar
;                If set the widget will not be activated.

@cw_draw2d

; the event
;------------------
pro plott_event, event

  ; get the plot data
  widget_control, event.top, get_uvalue=plotdata

  if plotdata.quiet then return

  ; x and y indices
  xidx= event.xidx
  yidx= event.yidx
  ; create a [x,y,z,dz] variable
  ; 2D plotdata
  if tag_exists(plotdata,'y') then begin
    if tag_exists(plotdata,'dz') then begin
      var=[plotdata.x[xidx],plotdata.y[yidx],plotdata.z[xidx,yidx],plotdata.dz[xidx,yidx]]
      print, format='("  plott_var = [x, y, z, dz] = [",g,", ",g,", ", g,", ",g,"]")', var
    endif else begin
      var=[plotdata.x[xidx],plotdata.y[yidx],plotdata.z[xidx,yidx]]
      print, format='("  plott_var = [x, y, z] = [",g,", ",g,", ", g,"]")', var
    endelse
  ; 1D plotdata
  endif else begin
    if tag_exists(plotdata,'dz') then begin
      var=[plotdata.x[xidx],plotdata.z[xidx],plotdata.dz[xidx]]
      print, format='("  plott_var = [x, y, dy] = [",g,", ",g,", ", g,"]")', var
    endif else begin
      var=[plotdata.x[xidx],plotdata.z[xidx]]
      print, format='("  plott_var = [x, y] = [",g,", ",g,"]")', var
    endelse
  endelse
  ; and save it one level up
  (scope_varfetch('plott_var', /enter, level=-1)) = var
end

; main program
;-------------
pro plott,  x, y , z, error=error, xrange=xrange, yrange=yrange, zrange=zrange,            $
                      linestyle=linestyle, thick=thick, color=color, psym=psym, symsize=symsize,$
                      errstyle=errstyle, errcolor=errcolor, nlevels=nlevels, clabels=clabels,   $
                      xgrid=xgrid, ygrid=ygrid, zgrid=zgrid,                                    $
                      xtitle=xtitle, ytitle=ytitle, ztitle=ztitle, title=title,                 $
                      xlbl=xlbl, ylbl=ylbl, charsize=charsize,                                  $
                      crossh=crossh, xidx=xidx, yidx=yidx, ptype=ptype, ctable=ctalble,         $
                      xsize=xsize, ysize=ysize, fcolor=fcolor, bcolor=bcolor,                   $
                      id=id, overplot=overplot, lgdpos=lgdpos,  extralgd=extralgd,              $
                      block=block, quiet=quiet, hold=hold

  ; set the default settings for the keywords
  default, linestyle, 0
  default, thick, 0
  default, errstyle, 1
  default, psym, 0
  default, symsize, 1
  default, xtitle, 'x'
  default, ytitle, 'y'
  default, ztitle, 'z'
  default, title, ''
  default, xlbl, 'data('+xtitle+')'
  default, ylbl, 'data('+ytitle+')'
  default, xidx, 0
  default, yidx, 0
  default, ptype, 0
  default, ctable, 5
  default, fcolor, 'black'
  default, bcolor, 'white'
  default, xsize, 900
  default, ysize, 675
  default, charsize, 1.0
  default, nlevels, 10
  default, clabels, 0
  default, lgdpos, 1
  default, quiet, 0
  default, block,0
  no_block = ~block
  default, hold, 0

  ; set up the plotdata structure (1D)
  if n_elements(z) eq 0 then begin
    if n_elements(y) eq 0 then begin
      y=x
      x=indgen(n_elements(y))
    endif
    if n_elements(x) ne n_elements(y) then begin
      print, ' ERROR: "x" and "y" do not have the same number of elements!'
      return
    endif
    plotdata={x:x, z:y, xtitle:xtitle, ztitle:ytitle, title:title, xsellbl:xlbl, $
              fcolor:fcolor, bcolor:bcolor, linethick:thick, linestyle:linestyle,$
              psym:psym, symsize:symsize, errstyle:errstyle, charsize:charsize,  $
              xidx:xidx, ptype:ptype, lgdpos:lgdpos, quiet: quiet}
    if n_elements(error) ne 0 then begin
      if n_elements(error) ne n_elements(x) then begin
        print, ' WARNING: "x" and "error" do not have the same number of elements!'
        print, '          the errorbars will not be plotted!'
      endif else begin
        plotdata=create_struct(plotdata,'dz',error)
      endelse
    endif
    if n_elements(xrange)   eq 2 then plotdata=create_struct(plotdata,'xrange',xrange)
    if n_elements(yrange)   eq 2 then plotdata=create_struct(plotdata,'zrange',yrange)
    if n_elements(color)    ne 0 then plotdata=create_struct(plotdata,'lcolor',color)
    if n_elements(errcolor) ne 0 then plotdata=create_struct(plotdata,'ecolor',errcolor)
    if n_elements(xgrid)    ne 0 then plotdata=create_struct(plotdata,'xgrid',xgrid)
    if n_elements(ygrid)    ne 0 then plotdata=create_struct(plotdata,'zgrid',ygrid)

  ; set up the plotdata structure (2D)
  endif else begin
    if n_elements(x) ne n_elements(z[*,0]) then begin
      print, ' ERROR: "x" and "z[*,0]" do not have the same number of elements!'
      return
    endif
    if size(y, /n_dim) le 1. then begin
      if n_elements(y) ne n_elements(z[0,*]) then begin
        print, ' ERROR: "y" and "z[0,*]" do not have the same number of elements!'
        return
      endif
    endif else begin
      if n_elements(x) ne n_elements(y[*,0]) then begin
        print, ' ERROR: "x" and "y[*,0]" do not have the same number of elements!'
        return
      endif
      if n_elements(y[0,*]) ne n_elements(z[0,*]) then begin
        print, ' ERROR: "y[*,0]" and "z[*,0]" do not have the same number of elements!'
        return
      endif
    endelse
    plotdata={x:x, y:y, z:z, xtitle:xtitle, ytitle:ytitle, ztitle:ztitle, title:title,  $
              xsellbl:xlbl, ysellbl:ylbl, fcolor:fcolor, bcolor:bcolor, linethick:thick,$
              linestyle:linestyle, psym:psym, symsize:symsize, errstyle:errstyle,       $
              charsize:charsize, xidx:xidx, yidx:yidx, ptype:ptype, clabels:clabels,    $
              nlevels:nlevels, lgdpos:lgdpos, quiet: quiet}
    if n_elements(error) ne 0 then begin
      if (n_elements(error[*,0]) ne n_elements(z[*,0])) $
         || (n_elements(error[0,*]) ne n_elements(z[0,*]))  then begin
        print, ' WARNING: "z[x,y]" and "error[x,y]" have a different size!'
        print, '          the errorbars will not be plotted!'
      endif else begin
        plotdata=create_struct(plotdata,'dz',error)
      endelse
    endif
    if n_elements(xrange)   eq 2 then plotdata=create_struct(plotdata,'xrange',xrange)
    if n_elements(yrange)   eq 2 then plotdata=create_struct(plotdata,'yrange',yrange)
    if n_elements(zrange)   eq 2 then plotdata=create_struct(plotdata,'zrange',zrange)
    if n_elements(color)    ne 0 then plotdata=create_struct(plotdata,'lcolor',color)
    if n_elements(errcolor) ne 0 then plotdata=create_struct(plotdata,'ecolor',errcolor)
    if n_elements(xgrid)    ne 0 then plotdata=create_struct(plotdata,'xgrid',xgrid)
    if n_elements(ygrid)    ne 0 then plotdata=create_struct(plotdata,'ygrid',ygrid)
    if n_elements(zgrid)    ne 0 then plotdata=create_struct(plotdata,'zgrid',zgrid)
  endelse

  if keyword_set(extralgd) then begin
    plotdata=create_struct(plotdata,'extralgd',extralgd)
  endif

  if keyword_set(overplot) then begin
    ; set the overplot field
    plotdata=create_struct(plotdata,'overplot',1)

    ; and add the data to the plot-tool
    widget_control, overplot, set_value=plotdata
    id = overplot
  endif else begin
    ; create a new window
    top  = widget_base(/column,title='Plot-tool')
    ; create colorselect
    id   = cw_draw2d(top,xsize=xsize,ysize=ysize, frame=0, event_pro='', $
                     crossh=crossh, ctable=ctable, background=bcolor     )
    ; start the widget
    widget_control,top,/realize

    ; and set/plot it
    widget_control, id, set_value=plotdata
    ; and save it to the top level base as well
    widget_control, top, set_uvalue=plotdata
    ; save the top-level base to the uvalue of the cw_draw2d widget
    widget_control, id, set_uvalue=top

  endelse

  ; make it active
  if ~hold then begin
    widget_control, id, get_uvalue=top
    xmanager,'plott',top, no_block=no_block
  endif

end
