function fids = listOpenedFiles()
%Get List of Open File Identifiers
%
% SYNOPSIS:
%   fids = listOpenedFiles()
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   fids - Array of open file identifiers.
%
% NOTE:
%   This is an internal function in the ECLIPSE result reader.  Its
%   semantics or existence may change without prior notice.
%
% SEE ALSO:
%   `fopen`.

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

   if exist('openedFiles', 'builtin')
      fids = openedFiles();
   else
      fids = fopen('all');
   end
end
