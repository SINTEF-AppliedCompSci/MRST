function N = tpfaExtendedConnections(G, bc, wells)
%Form extended connection (neighbourship) definition by wells and BCs.
%
% SYNOPSIS:
%   N = tpfaExtendedConnections(G, bc, wells)
%
% PARAMETERS:
%   G  - Grid structure.
%
%   bc - Boundary condition structure as defined by function addBC.  If
%        ISEMPTY(bc), then no additional connections will be formed or
%        enumerated for boundary conditions.
%
%   wells -
%        Well data structure as defined by function 'addWell'.  If
%        ISEMPTY(wells), the no additional connections will be formed or
%        enumerated for wells.
%
% RETURNS:
%   N - Connection (neighbourship) definition corresponding to
%       G.faces.neighbors, possibly extended by boundary conditions and/or
%       well perforations.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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


   N   = double(G.faces.neighbors);
   off = G.cells.num;

   % Contributions from boundary conditions
   if ~isempty(bc),
      % Merge cell and bc mobilities: Add entries in N to reflect
      % artificial outer neighbours

      bcf = double(bc.face);
      col = 1 + (N(bcf, 1) ~= 0);

      bc_conn = off + (1 : numel(bcf)).';

      N(sub2ind(size(N), bcf, double(col))) = bc_conn;

      off = off + numel(bcf);
   end

   % Contributions from wells
   if ~isempty(wells),
      wcell = vertcat(wells.cells);
      wconn = off + (1 : numel(wcell)).';

      N = [N; wconn, wcell];  % Pos flux is injection (-> perf in 1st col).
   end
end
