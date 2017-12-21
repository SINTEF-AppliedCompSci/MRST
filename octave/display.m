function display(msg)
  if(isstruct(msg))
    fname=fieldnames(msg);
    for i=1:numel(fname);
      printf(':%s\n',fname{i})
    end
  else
    %display_org(msg)
    builtin('display',msg)
  end
