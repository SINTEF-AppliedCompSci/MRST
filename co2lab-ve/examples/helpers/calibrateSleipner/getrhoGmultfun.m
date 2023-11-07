function v = getrhoGmultfun(obj, varargin)

   try
       v = getfield(obj, varargin{:});
   catch
       % field did not exist (multipler not yet added to fluid.  Return 1 as default.
       v = 1;
   end
   
end