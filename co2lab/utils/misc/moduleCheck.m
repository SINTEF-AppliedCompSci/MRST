function moduleCheck(varargin)
%Load specified modules unless already active
%
% SYNOPSIS:
%   moduleCheck m1 [m2 ...]
%
% PARAMETERS
%   m1, ... - Module names.  Must be known in mrstPath.  If a module is not
%             alread active, the module will be loaded through mrstModule.
%
% RETURNS:
%   Nothing.
%
% NOTE:
%   Use functional form of moduleCheck if the module name string is
%   contained in a variable.
%
% SEE ALSO:
%   `mrstPath`, `mrstModule`.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

   for module = reshape(varargin, 1, [])
      m = module{1};

      try
         require(m);
      catch
         fprintf('Loading module %s\n', m);
         mrstModule('add', m);
      end
   end
end
