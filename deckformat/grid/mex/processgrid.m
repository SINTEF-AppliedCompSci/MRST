function G = processgrid(varargin)
%Compute grid topology and geometry from pillar grid description.
%
% SYNOPSIS:
%   G = processgrid(grdecl)
%   G = processgrid(grdecl,ztol)
%
% PARAMETERS:
%   grdecl - Raw pillar grid structure, as defined by function
%            'readGRDECL', with fields COORDS, ZCORN and, possibly, ACTNUM.
%   ztol   - tolerance for unique points
%
% RETURNS:
%   G      - Valid grid definition containing connectivity, cell
%            geometry, face geometry and unique nodes.
%
% EXAMPLE:
%   G = processgrid(readGRDECL('small.grdecl'));
%   plotGrid(G); view(10,45);
%
% SEE ALSO:
%   `processGRDECL`, `readGRDECL`, `splitDisconnectedGrid`.

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

% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% Mex gateway by Jostein R. Natvig, SINTEF ICT.


   G = processgrid_mex(varargin{:});
   G.griddim = 3;
   G = splitDisconnectedGrid(G, 'Verbose', false);
end
