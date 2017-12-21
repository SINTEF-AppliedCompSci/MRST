function my_error(msg)
  if ~isempty(msg)
    builtin('error',msg);
  end
end