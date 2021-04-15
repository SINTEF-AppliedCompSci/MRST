function [cells, segment, param, ptsind, DT] = findWellPathCells(G, wellpath, varargin)
%Convert well path to the intersected cells
%
% SYNOPSIS:
%   cells = findWellPaths(G, wellpath);
%
% DESCRIPTION:
%   By creating a triangulation and mapping perforations to the closest
%   cells, this routine realizes the continuous well path into discrete
%   cells, making it possible to build simulation wells from it.
%
% REQUIRED PARAMETERS:
%  G        - The grid structure we want to realize the wells on.
%
%  wellpath - Well path. See "makeSingleWellpath" for spec.
%
%
% OPTIONAL PARAMETERS:
%
%  interpType - The type of interpolation used to extend the curve between
%               points. Supports the same types as MATLAB builtin interp1.
%               Default: Spline.
%
%  triangulation - The triangulation (typically from delaunayTriangulation)
%                  used to determine proximity in the grid.
%
%  refinement    - Refinement number used to further refine the well curve
%                  before computing which cells it intersect. Default: 100.
% RETURNS:
% 
%  cells      - The list of cells the well intersects. 
%
%  segment    - Segment indicator for each cell, indicating which wellpath
%               segment produced that specific completion. If multiple
%               choices are possible, the segment which comes first in
%               wellpath.points is used.
%
%  ptsind     - Point indicator, indicating which point in the segment was
%               the closest to a given cell. 
%
%  DT         - Triangulation used to produce the results.
%
%
% SEE ALSO:
%   `makeSingleWellpath`, `combineWellPaths`, `getWellFromPath`

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

opt = struct('interpType',      'spline', ...
             'triangulation',   [], ...
             'refinement',      100);
opt = merge_options(opt, varargin{:});

% Grab points
points = wellpath.points;
active = wellpath.active;
% Create a triangulation of the cell centroids. This is a bit
% expensive, but it makes it trivial to find the closest cell for each
% part of the well path.
DT = opt.triangulation;
if isempty(DT)
    DT = delaunayTriangulation(G.cells.centroids);
end

% Intersect is the cells intersected by the well. Note that we use a
% logical map to ensure that segments that are further down the well will
% not be selected to fill a cell over a segment that is higher up.
isIntersect = false(G.cells.num, 1);
% Segment will indicate for each perforated cell which well segment it
% originally came from. ptsind will indicate the same, but in terms of
% the points a cell produced. Param is a entry [0, 1] indicating how
% far along each segment the entries are.
[segment, ptsind, param] = deal(nan(G.cells.num, 1));

% We loop backwards to ensure that the points on the top overwrite
% those at the bottom. Otherwise we would get trouble with segments
% fighting over control over cells, as the well cannot be perforated
% multiple times in the same cell.
for i = numel(points):-1:1
    % Refine the curve to get higher accuracy.
    [pts, v] = refineSpline(points{i}, opt.refinement, opt.interpType);

    % Mask away any points corresponding to inactive subsegments.
    act = active{i};
    segv = min(floor(v), numel(act));
    pts = pts(act(segv), :);

    % Query triangulation to work out where we are in the world.
    cells = DT.nearestNeighbor(pts);

    % Store data
    isIntersect(cells) = true;
    segment(cells) = i;
    ptsind(cells) = round(v);
    param(cells) = v./max(v + sqrt(eps));
end
% Reduce to indices
cells = find(isIntersect);
segment = segment(isIntersect);
ptsind = ptsind(isIntersect);
param = param(isIntersect);

% Global enumeration of points according to parametrization + ordering
% in each segment. This is required for bottom hole pressure
% calculations, as MRST assumed the perforations to be ordered.
sv = param + segment;
[sv, sortInd] = sort(sv);

cells = cells(sortInd);
segment = segment(sortInd);
ptsind = ptsind(sortInd);
end