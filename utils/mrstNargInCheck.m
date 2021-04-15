function mrstNargInCheck(low, high, nargIn)
%Check number of input arguments to function
%
% SYNOPSIS:
%   mrstNargInCheck(low, high, nargIn)
%
% DESCRIPTION:
%   Fails execution (calls `error`) unless actual number of input arguments
%   to calling function is between lower and upper limits inclusive.  This
%   function should usually not be called within a loop as it is
%   implemented in terms of `dbstack`.
%
% PARAMETERS:
%   low    - Minumum number of input arguments needed by calling function.
%            If empty (i.e., if `isempty(low))`, this limit is not checked.
%
%   high   - Maximum number of input arguments allowed by calling function.
%            If empty (i.e., if `isempty(high))`, this limit is not checked.
%
%   nargIn - Actual number of input arguments.  Must be scalar, integral
%            and non-negative.  Should be `nargin` unless there are special
%            circumstances.
%
% SEE ALSO:
%   `nargin`, `nargchk`, `narginchk`, `error`, `dbstack`.

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

   assert (isnumeric(nargIn) && (numel(nargIn) == 1) && ...
           (mod(nargIn, 1) == 0) && (nargIn >= 0), ...
           'Argument count must be scalar, non-negative integer');

   st = dbstack;
   [caller, caller] = fileparts(st(2).file);                    %#ok<ASGLU>

   if ~isempty(low) && (nargIn < low)
      error('ArgCount:Low', ...
           ['Too few input arguments to function ''%''\n', ...
            'Expected at least %d arguments, but got %d'], ...
            caller, low, nargIn);
   end

   if ~isempty(high) && (nargIn > high)
      error('ArgCount:High', ...
           ['Too many input arguments to function ''%''\n', ...
            'Expected at most %d arguments, but got %d'], ...
            caller, high, nargIn);
   end   
end
