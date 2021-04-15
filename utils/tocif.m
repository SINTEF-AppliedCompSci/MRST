function tocif(bool, varargin)
%Evaluate function TOC if input is true.
%
% SYNOPSIS:
%   tocif(bool)
%   tocif(bool, tstart)
%
% PARAMETERS:
%   bool   - Boolean variable (true/false).
%
%   tstart - Saved start time as defined by functions `tic` or `ticif`.  In
%            this case, function `tocif` measures the elapsed time since
%            `tstart`.  Otherwise, the global `tic`/`toc` timer object is
%            used.
%
% NOTE:
%   Function used for making code cleaner where verbose option is used.
%
% SEE ALSO:
%   `toc`, `ticif`, `dispif`.

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


if bool,
   if nargin > 1,
      toc(varargin{1});
   else
      toc;
   end
end
end
