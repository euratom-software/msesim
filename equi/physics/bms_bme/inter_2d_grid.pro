;***************************************************
function inter_2d_grid, data, x, y, xnew, ynew
  x_inx = interpol(findgen(n_elements(x)), x, xnew)
  y_inx = interpol(findgen(n_elements(y)), y, ynew)

  inter_2d_grid = interpolate(data, x_inx, y_inx,/grid)

  return, inter_2d_grid
end
