function n = blockNeighbours(conn, b)
%Identify the neighbours of a particular coarse block.
%
% SYNOPSIS:
%   N = blockNeighbours(conn, b)
%
% PARAMETERS:
%   conn - Connection structure as defined by function 'blockConnectivity'.
%
%   b    - Particular block for which to determine the neighbouring coarse
%          blocks.
%
% RETURNS:
%   N - List (represented as an array) of other coarse blocks connected to
%       coarse block 'b'.
%
% SEE ALSO:
%   `blockConnectivity`.

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


   n = conn{b}(conn{b} ~= b);
end
