function Gl = makeLayers(G,nl,flayers)
% This function extrudes the matrix grid and independent fracture grid
% given the number of layers.
%
% SYNOPSIS:
%   Gl = makeLayers(G,nl,flayers)
%
% REQUIRED PARAMETERS:
%
%   G         - Matrix grid structure post fracture processing and
%               gridding.
%
%   nl        - Number of extruded layers.
%
%   flayers   - Indices of extruded layers in which fractures are present.
%
% RETURNS:
%   Gl  - Extruded matrix grid structure with extruded fracture grids
%         stored in Gl.FracGrid.
%
% SEE ALSO:
%   makeLayeredGrid

%{
Copyright 2009-2015: TU Delft and SINTEF ICT, Applied Mathematics.

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


nfl = numel(flayers);

% Layered matrix grid
Gl = computeGeometry(makeLayeredGrid(G,nl));
if isfield(G,'cartDims')
    Gl.cartDims = [G.cartDims nl];
end
FG = G.FracGrid;

% Layered fracture grid
cstart = Gl.cells.num+1; nstart = Gl.nodes.num+1; fstart = Gl.faces.num+1;
for i = 1:numel(fieldnames(FG))
    tempfl = makeLayeredGrid(FG.(['Frac',num2str(i)]),nfl);
    tempfl.nodes.coords(:,3) = tempfl.nodes.coords(:,3) + min(flayers) - 1;
    tempfl = computeGeometry(tempfl);
    tempfl.cells.start = cstart;
    tempfl.faces.start = fstart;
    tempfl.nodes.start = nstart;
    Gl.FracGrid.(['Frac',num2str(i)]) = tempfl;
    cstart = tempfl.cells.start + tempfl.cells.num;
    fstart = tempfl.faces.start + tempfl.faces.num;
    nstart = tempfl.nodes.start + tempfl.nodes.num;
end

return