;***************************************************
function inter_2d, data, x, y, xnew, ynew
  x_inx = interpol(findgen(n_elements(x)), x, xnew)
  y_inx = interpol(findgen(n_elements(y)), y, ynew)

  inter_2d = interpolate(data, x_inx, y_inx)

  return, inter_2d
end
