function G = emptyGrid(varargin)
%Make empty grid.
%
% SYNOPSIS:
%   G = emptyGrid()
%   G = emptyGrid(numcells)
%   G = emptyGrid(numcells, nodes)
%
% PARAMETERS:
%   numcells - Number of cells.  Used to allocate cells.facePos.  Treated
%              as zero if unspecified.
%
%   nodes    - Node coordinates.  N-by-d array.  Treated as ZEROS([0,3]) if
%              unspecified.
%
% RETURNS:
%   G - Grid structure as described by grid_structure.

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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

   numcells = 0;             if nargin > 0, numcells = varargin{1}; end
   nodes    = zeros(0, 3);   if nargin > 1, nodes    = varargin{2}; end


   % Empty shell:
   cells = struct('num',       double(numcells), ...
                  'facePos',   ones(numcells+1, 1), ...
                  'faces',     zeros(0, 2));
   faces = struct('num',       double(0), ...
                  'nodePos',   1, ...
                  'nodes',     zeros(0, 1), ...
                  'neighbors', zeros(0, 2));
   nodes = struct('num',       double(size(nodes, 1)),...
                  'coords',    nodes);

   G     = struct('cells',     cells, ...
                  'faces',     faces, ...
                  'nodes',     nodes, ...
                  'type',      {mfilename} );
