pro pngdump,filename=filename,noreverse=noreverse, raw=raw
; dumps the current window into a png file

if ~keyword_set(filename) then filename='screendump.png'

; read in the current windows image
img=tvrd(true=1)

if ~keyword_set(noreverse) then begin
  ; in normal circumstances we want to make the black background white and the white lines black
  ; (here we inverse all dark colors - coloridx<10 for R,G and B - and all light colors - coloridx>240 for R,G and B)
  im0 =fix(img[0,*,*]) ;R
  im1 =fix(img[1,*,*]) ;G
  im2 =fix(img[2,*,*]) ;B

  im0new = im0
  im1new = im1
  im2new = im2

  idx=where( (im0 lt 10) and (im1 lt 10) and (im2 lt 10))
  if idx[0] ne -1 then begin
    im0new[idx] = -im0[idx]+255
    im1new[idx] = -im1[idx]+255
    im2new[idx] = -im2[idx]+255
  endif

  idx=where( (im0 gt 240) and (im1 gt 240) and (im2 gt 240))
  if idx[0] ne -1 then begin
    im0new[idx] = -(im0[idx]-255)
    im1new[idx] = -(im1[idx]-255)
    im2new[idx] = -(im2[idx]-255)
  endif

  img[0,*,*]=byte(im0new)
  img[1,*,*]=byte(im1new)
  img[2,*,*]=byte(im2new)
endif
tv,img,true=1
write_png,filename,img

end
