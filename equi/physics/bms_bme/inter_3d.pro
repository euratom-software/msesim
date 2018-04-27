;***************************************************
function inter_3d, data, x, y, z, xnew, ynew, znew
  x_inx = interpol(findgen(n_elements(x)), x, xnew)
  y_inx = interpol(findgen(n_elements(y)), y, ynew)
  z_inx = interpol(findgen(n_elements(z)), z, znew)

  inter_3d = interpolate(data, x_inx, y_inx, z_inx,/grid)

  return, inter_3d
end
