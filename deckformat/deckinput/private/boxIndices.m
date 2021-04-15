function [ix, varargout] = boxIndices
%Construct linear (cell) indices corrsponding to current input box
%
% SYNOPSIS:
%   i     = boxIndices
%  [i, m] = boxIndices
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   i - Linear (cell) indices corresponding to the current input box.
%
%   m - Maximum possible cell index in index set 'i'.  Corresponds to
%       PROD(defaultBox) and may exceed MAX(i) (but should never be less
%       than MAX(i)).
%
% SEE ALSO:
%   `private/gridBox`, `private/defaultBox`.

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


   i    = gridBox;
   dims = defaultBox;
   dims = dims(2 : 2 : end);

   [j{1 : 3}] = ndgrid(i(1) : i(2), i(3) : i(4), i(5) : i(6));

   ix = sub2ind(reshape(dims, 1, []), ...
                reshape(j{1}, [], 1), ...
                reshape(j{2}, [], 1), ...
                reshape(j{3}, [], 1));

   if nargout > 1,
      varargout{1} = prod(dims);

      assert (varargout{1} >= ix(end), ...
             ['Maximum possible cell index less than largest ', ...
              'index in current input box?']);
   end
end
