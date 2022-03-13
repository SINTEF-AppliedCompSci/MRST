function varargout = plotPlume(G, Gt, h, varargin)
%Plot CO2 plume on logically Cartesian grid
%
% SYNOPSIS:
%       plotPlume(G, Gt, h)
%       plotPlume(G, Gt, h, cells)
%       plotPlume(G, Gt, h, 'pn1', pv1, ...)
%       plotPlume(G, Gt, h, 'pn1', pv1, ...)
%   h = plotPlume(...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   Gt      - Corresponding top grid as defined by topSurfaceGrid
%
%   h       - height dataset for the topSurfaceGrid
%
%   cells   - (OPTIONAL) Cells of G to include in plot.
%
%   'pn'/pv - List of other property name/value pairs.  OPTIONAL.
%             This list will be passed directly on to function PATCH
%             meaning all properties supported by PATCH are valid.
%
% RETURNS:
%   h - Handle to resulting PATCH object.  The patch object is added to the
%       current AXES object.
%
% NOTES:
%   Function 'plotPlume' is implemented directly in terms of the low-level
%   function PATCH.  If a separate axes is needed for the graphical output,
%   callers should employ function newplot prior to calling 'plotPlume'.
%
% SEE ALSO:
%   `plotCellData`, `plotGrid`, `newplot`, `patch`, `shading`.

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


% Create interpolant of height
T = TriScatteredInterp(Gt.cells.centroids(:,1), Gt.cells.centroids(:,2), h);
topFaces = Gt.cells.map3DFace;
cells = sum(G.faces.neighbors(topFaces,:), 2);

% Find top and bottom faces for upper layer of cells
fp = mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1);
tmp = G.cells.faces(fp, :);
bottomFaces = tmp(tmp(:,2) == 6, 1);

% Find the nodes of the bottom faces of the layer of cells
nodePos = mcolon(G.faces.nodePos(bottomFaces), G.faces.nodePos(bottomFaces + 1) - 1);
nodeBottom = (G.faces.nodes(nodePos));
c = G.nodes.coords(nodeBottom, :);

% Find the nodes of the top faces of the layer of cells
nodePos = mcolon(G.faces.nodePos(topFaces), G.faces.nodePos(topFaces + 1) - 1);
nodeTop = (G.faces.nodes(nodePos));

% Adjust the nodes of the bottom faces to create the bottom of the
% plume based on the interpolated values
G.nodes.coords(nodeBottom, 3) = G.nodes.coords(nodeTop, 3) + T(c(:,1), c(:,2));

H = h > 1e-3;

plumeCells = sum(G.faces.neighbors(Gt.cells.map3DFace(H),:), 2);
dummy = zeros(G.cells.num, 1);
dummy(plumeCells) = h(H);

ind = false(G.cells.num, 1);
ind(plumeCells) = true;

% Respect user subset selection by going via the and operator
if mod(numel(varargin), 2) == 1
    cc = varargin{1};
    if ~islogical(cc)
        tmp = false(G.cells.num, 1);
        tmp(cc) = true;
        cc = tmp;
    end
    ind = ind & cc;
    varargin = varargin(2:end);
end

% Finally plot and pass on arguments
if any(ind)
    h = plotCellData(G, dummy, ind, varargin{:});
else
    h = nan;
end

if nargout > 0
    varargout{1} = h;
end

end
