function varargout = my_regexp(str,pat,varargin)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

   try
      [varargout{1:nargout}] = builtin("regexp",str,pat,varargin{:});
   catch
      if (nargin > 2) && ischar(varargin{1}) && strcmp(varargin{1},'split')
         assert(pat==pathsep);
         varargout{1}=strsplit(str,pat);
      else
         rethrow(lasterr);
      end
   end
end
