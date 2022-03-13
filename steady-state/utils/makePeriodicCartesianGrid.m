function [Gp, bcp] = makePeriodicCartesianGrid(G)
%Undocumented Utility Function

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

% Uses the MRST module 'upscaling'
require upscaling

% Only make the dimensions with cartDims>1 periodic
dims  = find(G.cartDims>1);

bf    = getBoundaryFaces(G);
ndims = numel(dims);
bcl   = cell(1, ndims); % "left" faces
bcr   = cell(1, ndims); % "right" faces
for i = 1:ndims
    d = dims(i);
    bcl{i}.face = bf{d,1};
    bcr{i}.face = bf{d,2};
end

% Set pressure drop to zero for all directions as default
dp = cell(1, ndims);
dp(:) = {0};

% Create periodic grid
[Gp, bcp] = makePeriodicGridMulti3d(G, bcl, bcr, dp, 'dims', dims);
Gp.parent = G;

end

function bfaces = getBoundaryFaces(G)
% Finds the boundary faces of a grid.
% 
% This function finds the boundary faces of the given grid and returns
% the face indecies as a cell structure.

tags = [1 2; 3 4; 5 6];
boundaryFace = any(G.faces.neighbors == 0, 2);
ind = boundaryFace(G.cells.faces(:,1));
faceAndTag = G.cells.faces(ind, :);

bfaces = cell(G.griddim,2);
for i = 1:G.griddim
   bfaces{i,1} = faceAndTag(faceAndTag(:,2) == tags(i,1), 1);
   bfaces{i,2} = faceAndTag(faceAndTag(:,2) == tags(i,2), 1);
end

end
