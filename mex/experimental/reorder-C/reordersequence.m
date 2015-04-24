function varargout = reordersequence(G, flux)
%Compute topological sorted strong components of outflow graph (G, flux)
%
% SYNOPSIS:
%   s      = reordersequence(G, flux);
%   [s, b] = reordersequence(G, flux);
%
% Together, the grid G and the face fluxes form a directed graph, where the
% cells a vertices and directed edges exist from vertex i to vertex j if
% there is a positive flux from cell i to cell j.  This graph may be termed
% the outflow graph, and its adjacency matrix is the inflow matrix.
%
% The function REORDERSEQUENCE returns a topological sorted sequence of
% strong components of the inflow graph computed using Tarjans algorithm.
%
% PARAMETERS:
%   G    - Grid structure.  Geometric primitives are not needed by this
%          function.
%   flux - Vector of face fluxes.
%
% RETURNS:
%   s    - vector of cell numbers.
%
%   b    - vector of positions for strong components: the i'th strong
%          component consists of the cells s(b(i):b(i+1)-1),
%          numel(b) = #(strong components) + 1
%
% EXAMPLE:
%   % Make a 3-by-2 grid, reordersequence needs no geometry...
%      G = cartGrid([3, 2]);
%
%   % make fake flux field
%      flux = ones(17, 1);
%      flux(6)  = -1;
%      flux(12) = -1;
%
%   % Plot the grid in 3D-view.
%      [s, b] = reordersequence(G, flux);
%
% SEE ALSO:
%

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


   % ensure fields are int32
   G.cells.faces     = int32(G.cells.faces);
   G.cells.facePos   = int32(G.cells.facePos);
   G.faces.neighbors = int32(G.faces.neighbors);

   % Call MEX edition.
   [varargout{1:nargout}] = reordersequence_mex(G, flux);
end
