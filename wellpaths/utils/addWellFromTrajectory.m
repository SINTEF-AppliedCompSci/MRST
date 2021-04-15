function W = addWellFromTrajectory(W, G, rock, traj, varargin)
% This function adds a well based on a piecewise linear trajectory traj.
% All function arguments are the same as for addWell except for in place of
% cellInx, there should be a nx3 matrix of trajectory coordinates.

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

opt = struct('faces',                  [], ...
             'exteriorFaceCorrection', false, ...
             'refDepth',               []);
[opt, other] = merge_options(opt, varargin{:});
assert(G.griddim == 3, 'Function only compatible with 3D grids.');
assert(size(traj,2)==3, 'Trajectory is expected to be nx3 array.');

if ~isfield(G.faces, 'bbox')
    G = addBoundingBoxFields(G);
end

tmp = computeTraversedCells(G, traj, 'faces', opt.faces, ...
            'exteriorFaceCorrection', opt.exteriorFaceCorrection);

if isempty(tmp.cell)
    error('Did not find any traversed cells for trajectory.');
end

if isempty(opt.refDepth) && isfield(G.cells, 'centroids')
    opt.refDepth = G.cells.centroids(tmp.cell(1), 3);
end
% multiply segments by weight (non-unit for segments shared by multiple cells)
seg = bsxfun(@times, tmp.vec, tmp.weight);
W = addWell(W, G, rock, tmp.cell, 'lineSegments', seg, 'refDepth', opt.refDepth, other{:});
end
