function G = cartGrid(celldim, varargin)
%Construct 2d or 3d Cartesian grid in physical space.
%
% SYNOPSIS:
%   G = cartGrid(celldim);
%   G = cartGrid(celldim, physdim);
%
% PARAMETERS:
%   celldim  - Vector, length 2 or 3, specifying number of
%              cells in each coordinate direction.
%   physdim  - Vector, length numel(celldim), of physical size in units of
%              meters of the computational domain.
%              OPTIONAL.  Default value == celldim
%              (i.e., each cell has physical dimension 1-by-1-by-1 m).
%
%   cellnodes- OPTIONAL.
%              Default value FALSE.  If TRUE, the corner points of each
%              cell is added as field G.cellNodes.  The field has one row
%              per cell, the sequence of nodes on each is (imin,
%              jmin,kmin), (imax,jmin,kmin), (imin,jmax,kmin), ...
%
% RETURNS:
%   G - Grid structure with a subset of the fields `grid_structure`.
%       Specifically, the geometry fields are missing:
%         - G.cells.volumes
%         - G.cells.centroids
%
%         - G.faces.areas
%         - G.faces.normals
%         - G.faces.centroids
%
%       These fields may be computed using the function `computeGeometry`.
%
%       There is, however, an additional field not described in
%       `grid_structure:
%
%           `cartDims` is a length 2 or 3 vector giving number of cells in
%           each coordinate direction.  In other words 
%
%                      `all(G.cartDims == celldim)`.
%
%       `G.cells.faces(:,2)` contains integers 1-6 corresponding to
%       directions W, E, S, N, T, B respectively.
%
% EXAMPLE:
%   % Make a 10-by-5-by-2 grid on the unit cube.
%      nx = 10; ny = 5; nz = 2;
%      G = cartGrid([nx, ny, nz], [1, 1, 1]);
%
%   % Plot the grid in 3D-view.
%      f = plotGrid(G); view(3);
%
% SEE ALSO:
%   `grid_structure`, `tensorGrid`, `computeGeometry`

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


if ~all(celldim > 0)
   error('CELLDIM must be positive');
end

physdim = celldim;
if nargin > 1 && isnumeric(varargin{1})
   physdim  = varargin{1};
   varargin = varargin(2 : end);
end

dim = numel(celldim);

x = linspace(0, physdim(1), celldim(1)+1);
if dim > 1
    y = linspace(0, physdim(2), celldim(2)+1);
end

switch dim
    case 1
        G = tensorGrid(x, varargin{:});
    case 2
        G = tensorGrid(x, y, varargin{:});
    case 3
        z = linspace(0, physdim(3), celldim(3)+1);
        G = tensorGrid(x, y, z, varargin{:});
    otherwise
        error('Cannot create grid with %d dimensions: Only 1, 2 or 3 is valid.', dim);
end

% Record grid constructur in grid.
G.type    = [G.type, { mfilename }];