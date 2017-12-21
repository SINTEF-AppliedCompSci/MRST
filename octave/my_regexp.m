function varargout = my_regexp(str,pat,varargin)
  try
    [varargout{1:nargout}] = builtin("regexp",str,pat,varargin{:});
  catch
    if(nargin>2 && ischar(varargin{1}) && strcmp(varargin{1},'split') )
      assert(pat==pathsep);
      varargout{1}=strsplit(str,pat);
    else
      rethrow(lasterr);
    end
  end
end
