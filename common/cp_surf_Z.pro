;---------------------------------------------------------------------------
; Routine: cp_surf_Z
; Version: 1.00
; Author : M.De Bock
; Date   : 19.02.09
;---------------------------------------------------------------------------
;
; almost identical to cp_shade, with the exception that the shade surface
; is seen under the standard 30deg angles, instead of from above
; contour lines are plotted above the surface. This routine uses a Z-buffer.
;
;---------------------------------------------------------------------------
; Routine: cp_shade 
; Version: 1.20
; Author : R.Martin
; Date   : 8.12.03
;---------------------------------------------------------------------------
; cp_shade
;
; Compond plot, display data as a shaded 2D plot. Optional contour lines,
; and scale bar display on right hand of plot.
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
; CP_SHADE, data {,x ,y}
;
; Additional parameters
;   /invert    - reverse colour scheme used for shading
;   /showscale - Display scale bar
;   scalewidth - Real (0.005-0.9) width of scalebar as fraction of
;                width of plot, default 1/20th width of plot
;   /clabel       - Label contours on plot
;   cvalues       - Contour levels to be displayed
;   /autocontour  - 
;   
;
; The following standard IDL graphics parameters are accepted by CP_SHADE
;
;   title, xtitle, ytitle: As for the PLOT command
;   ztitle      : Appears on the scale-bar axis.
;   {xyz}range  : As usual.
;
;   noerase     : Just there for compatability NO LONGER THERE
;   downhill    : As for CONTOUR command 
;   position    : Position of plot on screen in normalised 
;                 co-ordinates [x0,y0,x1,y1]
;
; 19/02/2009, mdebock: minor fixes => keyword charsize introduced
;                                  => keyword scalelabel now used
;                                  => small changes in the positioning,
;                                     such that the scale label fits on
;                                     the plot
;---------------------------------------------------------------------------
;
;

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

pro cp_shade_loadct, ctable, reset=reset

  if (!d.name eq 'X') || (!d.name eq 'Z') then device, decomp=keyword_set(reset)
  if is_integer(ctable) then loadct, abs(ctable(0)) mod 40, /silent

end

pro cp_surf_Z, data, x, y,              $;
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
              ax=axarg, az=azarg,      $;
              error=error

;Check Arguement-1, must be 2D real array

  color_save, ctable

  if not_real(data) then begin
    error_message='CP_SURF: Arg(1) must be a real-array'
    goto, errorcatch
  endif
  
  dsize=size(data)
  if (dsize(0) ne 2) then begin
    error_message='CP_SURF: Arg(1) must be a 2D real-array'
    goto, errorcatch
  endif

  dsize=size(data)
  if (dsize(1) le 1) or (dsize(2) le 1) then begin
    error_message='CP_SURF: Arg(1) X/Y dimensions must be .gt. 1'
    goto, errorcatch
  endif

;Check Arguement-2, if defined must be 1D real array. 
; Length must be longer than X-dimension of data array.

  if (n_elements(x) eq 0) then begin
    ix=findgen(dsize(1))
  endif else begin
    if not_real(x) then begin
      error_message='CP_SURF: Arg(2) must be a real-array'
      goto, errorcatch
    endif
    
    if (n_elements(x) lt dsize(1)) then begin
      error_message='CP_SURF: Arg(2) size not compatible with Arg(1)'
      goto, errorcatch
    endif

    if (monotonic(x, /inc, /strict) eq 0) then begin
      error_message='CP_SURF: Arg(2) must be monotonically increasing'
      goto, errorcatch
    endif

    ix=x(0:(dsize(1)-1))
  endelse

;Check Arguement-2, if defined must be 1D real array.
;  Lenght must be longer than Y-dimension of data array.

  if (n_elements(y) eq 0) then begin
    iy=findgen(dsize(2))
  endif else begin
    if not_real(y) then begin
      error_message='CP_SURF: Arg(3) must be a real-array'
      goto, errorcatch
    endif

    if (n_elements(y) lt dsize(2)) then begin
      error_message='CP_SURF: Arg(3) size not compatible with Arg(1)'
      goto, errorcatch
    endif

    if (monotonic(y, /inc, /strict) eq 0) then begin
      error_message='CP_SURF: Arg(3) must be monotonically increasing'
      goto, errorcatch
    endif

    iy=y(0:(dsize(2)-1))
  endelse

; check the showscale keyword
  showscale=keyword_set(showscalearg) or is_real(scalewidth)  
; set the character size
  oldcharsize=!p.charsize
  if keyword_set(charsizearg) then charsize=float(charsizearg) else charsize=0
  if charsize ne 0 then !p.charsize=(oldcharsize eq 0) ? charsize : oldcharsize*charsize

;Set XRANGE.
;     - IDL doesn't clip SHADE_SURF plots so that data array needs
;       to be cut down in size, if the xrange is set.
;     - Otherwise the XRANGE is calculated from the x-coordinate
;     - If no data is in the given XRANGE a blank plot is shown

  idata=data
  if is_real(xrange) then begin
    xmin=min(xrange, max=xmax)
    if (xmin eq xmax) then begin
      error_message='CP_SURF: X-Range badly defined (xmin=xmax)'
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

      idata(i, *)=idata(i, *)+(idata(i,*)-idata(i+1,*))*(xmin-ix(i))/(ix(i)-ix(i+1))
      idata=idata(i:*,*)
      ix=ix(i:*)
      ix(0)=xmin
    endif

    if (xmax lt ixmax) then begin
      i=value_locate(ix, xmax)
      idata(i+1, *)=idata(i,*)+(idata(i,*)	$
		-idata(i+1,*))*(xmax-ix(i))/(ix(i)-ix(i+1))

      idata=idata(0:(i+1), *)
      ix=ix(0:(i+1))
      ix(i+1)=xmax
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
      error_message='CP_SURF: Y-Range badly defined (xmin=xmax)'
      goto, errorcatch
    endif

    iymin=min(iy, max=iymax)
    if (ymin gt iymax) or (ymax lt iymin) then goto, endplot

    if (ymin gt iymin) then begin
      i=value_locate(iy, ymin)

      idata(*, i)=idata(*, i)+(idata(*, i)-idata(*, i+1))*(ymin-iy(i))/(iy(i)-iy(i+1))
      idata=idata(*, i:*)
      iy=iy(i:*)
      iy(0)=ymin
    endif

    if (ymax lt iymax) then begin
      i=value_locate(iy, ymax)
      idata(*, i+1)=idata(*, i)+(idata(*, i)	$
		-idata(*, i+1))*(ymax-iy(i))/(iy(i)-iy(i+1))

      idata=idata(*, 0:(i+1))
      iy=iy(0:(i+1))
      iy(i+1)=ymax
    endif
  endif else begin
    ymin=min(iy, max=ymax)
    yrange=[ymin,ymax]
  endelse

;Check ZRANGE
;     - If ZRANGE is not defined a dummy call to PLOT is used
;       to determine the extent of the ZRANGE.
;     - The dummy call is also used to calculate the default levels
;       of the contours.
;     - The NOERASE option is used with the PLOT command, this over comes
;       the problem of using the erase command on a PS plot.
;    

  yticks=!y.ticks                       ;Store Y-tick value
  if is_integer(autocontour) then begin
    !y.ticks=(0>autocontour(0))<60

    showcontours=1B
  endif

  if is_real(zrange) then begin
    zmin=min(zrange, max=zmax)
    if (zmin eq zmax) then begin
      error_message='CP_SURF: Z-Range badly defined (xmin=xmax)'
      goto, errorcatch
    endif  

    plot, [0,1], [zmin, zmax], /nodata,   $
          xstyle=5,                       $
          ystyle=5,                       $
          ytick_get=zvalues,              $
          noerase=noerase
  endif else begin
    zmin=min(idata, max=zmax, /nan)
    plot, [0,1], [zmin, zmax], /nodata,   $
          xstyle=5,                       $
          ystyle=4,                       $
          ytick_get=zvalues,              $
          noerase=noerase

    zmin=!y.crange(0)
    zmax=!y.crange(1)
  endelse
  !y.ticks=yticks                       ;Restore Y-tick value

;Check position of plot on screen. 
;       - If ARG(position) isn't defined get the position of the dummy plot
;         from the system variables
;

  if not_real(position) then begin
    plotpos=[!x.window, !y.window]
    plotpos=plotpos([0,2,1,3])
    plotpos(0) = 1.5*plotpos(0)
  endif else begin
    if (n_elements(position) ne 4) then begin
      error_message='CP_SURF: Position-arg wrong size position=[X0,Y0,X1,Y1]'
      goto, errorcatch
    endif

    if (min(position) lt 0) or (max(position) gt 1) then begin
      error_message='CP_SURF: Position-arg wrong plot outside window'
      goto, errorcatch
    endif
    
    if (position(0) eq position(2)) or (position(1) eq position(3)) then begin
      error_message='CP_SURF: Position-arg badly defined plot region'
      goto, errorcatch
    endif      

    plotpos=position
    if (plotpos(0) ge plotpos(2)) then plotpos([0,2])=plotpos([2,0])
    if (plotpos(1) ge plotpos(3)) then plotpos([1,3])=plotpos([3,1])
  endelse

  if keyword_set(isotropic) then begin
    plot, [xmin, xmax], [ymin, ymax], $
          xstyle=5, ystyle=5, /isotropic, /nodata

    plotpos=[!x.window, !y.window]
    plotpos=plotpos([0,2,1,3])
    plotpos(0) = 1.5*plotpos(0)
  endif 

;Position of ScaleBar:
;     - If the /SHOWSCALE option is choosen a scale will be placed on the 
;       right hand side of the plot. The default width is 1/20th the width of
;       the plot.
;     - The width can be changed using the SCALEWIDTH option
;     - The scale bar can be placed anywhere on the plot using the
;       SCALE_POSTION option



; We need the Z-buffer device to deal with 'hidden' surfaces
  ; get the resolution of the current device
  res = convert_coord([1],[1],/norm,/to_device)
  ; the Z-buffer device will always have white (255) as foreground and black (0) as background
  ; so we need to save the desired back and foreground color, so we can replace it at the end
  bcolor = !p.background
  fcolor = !p.color
  ; save the name of the current device and find out whether it is in decomposed mode
  origdevice = !d.name
  device, get_decomp=origdecomp
  ; open Z-buffer device
  set_plot, 'Z'
  device, set_resolution=res[0:1], set_pixel_depth=8, decomposed=0

  if not_real(scalewidtharg) then scalewidth=0.05 else scalewidth=scalewidtharg[0]
  if ~keyword_set(scalelabelarg) then scalelabel=3 else scalelabel=floor(abs(scalelabelarg mod 4))
  if keyword_set(showscale) then begin
    xyouts, 0.0, 0.0, ' ', /normal, width=chwidth

    if not_real(scaleposarg) then begin
      scalepos=plotpos

      if (scalewidth lt 0) or keyword_set(isotropic) then begin
        scalepos(0)=plotpos(2)+0.5*chwidth
        scalepos(2)=plotpos(2)+0.5*chwidth+(plotpos(2)-plotpos(0))*(0.005>abs(scalewidth(0))<0.9)
      endif else begin        
        plotpos(2)=plotpos(2)- 0.5*chwidth - $
                   (plotpos(2)-plotpos(0))*(0.005>scalewidth(0)<0.9)
        scalepos(0)=plotpos(2)+0.5*chwidth
      endelse
    endif else begin
      scalepos=scaleposarg
      if (n_elements(scalepos) ne 4) then begin
        error_message='CP_SURF: Scale Position-arg wrong size scale_position=[X0,Y0,X1,Y1]'
        goto, errorcatch
      endif

      if (min(scalepos) lt 0) or (max(scalepos) gt 1) then begin
        error_message='CP_SURF: Scale Position-arg wrong plot outside window'
        goto, errorcatch
      endif

      if (scalepos(0) eq scalepos(2)) or (scalepos(1) eq scalepos(3)) then begin
        error_message='CP_SURF: Scale Position-arg bad size'
        goto, errorcatch
      endif

      if (scalepos(0) ge scalepos(2)) then scalepos([0,2])=scalepos([2,0])
      if (scalepos(1) ge scalepos(3)) then scalepos([1,3])=scalepos([3,1])

    endelse
    ; make room for the scale label if it is positioned on the y-axis of the scale
    if (scalelabel ge 2) then begin
      plotpos(2) = plotpos(2)- 4.*chwidth
      if scalelabel eq 3 then begin
        scalepos(0) = scalepos(0)- 4.*chwidth
        scalepos(2) = scalepos(2)- 4.*chwidth
      endif
      bx=[0,1]
      by=[zmin, zmax]
      bdata=[1,1]#by
    endif else begin
      bx=[zmin, zmax]
      by=[0,1]
      bdata=bx#[1,1]
    endelse

    cp_shade_loadct, ctable

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
    ntick=n_elements(tickvalues)
    dtick=(tickvalues[ntick-1]-tickvalues[0])/(ntick-1)
    index=where(abs(tickvalues/dtick) lt 5e-7, count)
    if (count gt 0) then tickvalues[index]=0.0
    chwidth=replicate(chwidth, n_elements(tickvalues))
    tickstr=replicate('   ', n_elements(tickvalues))

    cp_shade_loadct, /reset

    tmpcolor=!p.color
    if exists(color) then !p.color=truecolor(color)

    case scalelabel of 
      0:axis, xaxis=0, xstyle=1, xticklen=0.3, xtitle=ztitle
      1:axis, xaxis=1, xstyle=1, xticklen=0.3, xtitle=ztitle
      2:axis, yaxis=0, ystyle=1, yticklen=0.3, ytitle=ztitle, ytickname=tickstr
      3:axis, yaxis=1, ystyle=1, yticklen=0.3, ytitle=ztitle, ytickname=tickstr
    endcase

    if (scalelabel ge 2) then begin
      ytickpos=(tickvalues-zmin)/(zmax-zmin)*(scalepos(3)-scalepos(1))+scalepos(1)

      if (scalelabel eq 2) then xtickpos=scalepos[0]-chwidth*1.5 $
      else xtickpos=scalepos[2]+chwidth*2.5
      xyouts, xtickpos, ytickpos, /normal,      $
              align=0.5, nice_number(float(tickvalues)), orient=90,  $
              charsize=charsize
    endif

    plots, scalepos([0,0,2,2,0]), scalepos([1,3,3,1,1]), /norm
    !p.color=tmpcolor
  endif

;Add shade-surface plot
  ; first plot the axis
  if keyword_set(axarg) then ax=float(axarg) else ax=30.
  if keyword_set(azarg) then az=float(azarg) else az=30.
  plot_3dbox, [xmax,xmax],[ymax,ymax], [zmin, zmax], /noerase, $
              position=plotpos,                                $
              ax=30, az=30, /save,                             $
              xrange=xrange, xstyle=1, xtitle=xtitle,          $
              yrange=yrange, ystyle=1, ytitle=ytitle,          $
              zrange=[zmin, zmax], zstyle=1,ztitle=ztitle,     $
              title=title, charsize=2.0*!p.charsize

  ; then add the contours (on the bottom), define contour levels
  if (n_elements(cvalues) gt 0) then begin
    if not_real(cvalues) then begin
      error_message='CP_SURF: C_LEVEL must be real-array'
      goto, errorcatch
    endif

    zvalues=cvalues
    showcontours=1B
  endif

  tmpcolor=!p.color
  if keyword_set(showcontours) then begin
    if exists(ccolor) then !p.color=truecolor(ccolor)
    ord=where(zvalues gt zmin and zvalues lt zmax, num)
    if (num gt 0) then begin
      contour, idata, ix, iy, /overplot,  $
               level=zvalues, zval=0.0,   $
               downhill=downhill,         $
               c_label=replicate(keyword_set(clabel), num+1),$
               c_charsize=1.5*!p.charsize,/T3D
    endif
  endif

  if is_integer(zerocontour) then begin
    if (zmin lt 0) and (zmax gt 0) then     $
      contour, idata, ix, iy, /overplot,    $
               level=0, zval=0.0,           $
               thick=(1>zerocontour(0))<10, $
               c_label=1, c_charsize=1.5*!p.charsize, /T3D
  endif
  !p.color=tmpcolor

  ; finally plot the surface
  cp_shade_loadct, ctable
  if keyword_set(invert) then begin
    shade_surf, idata, ix, iy, /noerase,                       $
                position=plotpos,                              $
                shades=255B-bytscl(idata, min=zmin, max=zmax), $
                /T3D,                                          $
                xrange=xrange, xstyle=5, xtitle=xtitle,        $
                yrange=yrange, ystyle=5, ytitle=ytitle,        $
                zrange=[zmin, zmax], zstyle=5,ztitle=ztitle,   $
                title=title, charsize=2.0*!p.charsize
  endif else begin
     shade_surf, idata, ix, iy, /noerase,                       $
                 position=plotpos,                              $
                 shades=bytscl(idata, min=zmin, max=zmax),      $
                 /T3D,                                          $
                 xrange=xrange, xstyle=5, xtitle=xtitle,        $
                 yrange=yrange, ystyle=5, ytitle=ytitle,        $
                 zrange=[zmin, zmax], zstyle=5,ztitle=ztitle,   $
                 title=title, charsize=2.0*!p.charsize
  endelse
  cp_shade_loadct, /reset


  ; reset 3D view
  T3D, /reset


endplot:

  ; Read the Z-buffer
  img=tvrd()
  ; find out where the background and foreground pixels are
  bpixel = where(img le !p.background+1, nb)
  fpixel = where(img ge !p.color-1, nf)
  ; set the original foreground and background colors
  img[fpixel]=byte(fcolor)
  img[bpixel]=byte(bcolor)

  ; reload to the original device
  set_plot, origdevice
  ; set it to indexed colors
  device, decomp=0
  ; plot the image
  tv, img
  ; and reset the original decomposed-setting
  device, decomp=origdecomp


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

;---------------------------------------------------------------------------------------
; Modification History
;
; 19.08.08 - Improve colour handling
;             - Add LOADCT option to define colour table for plot
;             - Add COLOR option to define axis colour
;
;
;
;
; 21.03.06 - Add in new ERRORCATCH section
;              - Modify all error trapping to goto error catch section.
;              - Add ERROR option to allow external error trapping
;
; 15.04.05 - Add /nan option to MIN call used to find zmin/zmax
;

