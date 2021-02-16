function res = ls(varargin)
% avoir problem when calling 'ls' on a nonexistant file.  In that case, 'dir'
% returns an empty array, which is what we'd expect from 'ls' too in a Matlab context.
   
   if nargout > 0
      res = dir(varargin{:});
      
      % remove directories '.' and '..'
      remove_ix = arrayfun(@(n) strcmpi(n.name, '.') || strcmpi(n.name, '..'), res);
      res(remove_ix) = [];
      
      res = [arrayfun(@(n) [n.name, ' '], res, 'uniformoutput', false){:}];
      
   else
      dir(varargin{:});
   end
end
