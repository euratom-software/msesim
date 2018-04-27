pro default, var, val
  ; sets default value for a variable if the variable does not yet exist
  if n_elements(var) eq 0 then var=val
end
