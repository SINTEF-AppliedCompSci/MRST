function W = addWellFromTrajectory(W, G, rock, traj, varargin)
% This function adds a well based on a piecewise linear trajectory traj.
% All function arguments are the same as for addWell except for in place of
% cellInx, there should be a nx3 matrix of trajectory coordinates.
opt = struct('faces',                  [], ...
             'exteriorFaceCorrection', false, ...
             'refDepth',               []);
[opt, other] = merge_options(opt, varargin{:});
assert(G.griddim == 3, 'Function only comaptible with 3D grids.');
assert(size(traj,2)==3, 'Trajectory is expected to be nx3 array.');

if ~isfield(G.faces, 'bbox')
    G = addBoundingBoxFields(G);
end

tmp = computeTraversedCells(G, traj, 'faces', opt.faces, ...
            'exteriorFaceCorrection', opt.exteriorFaceCorrection);
        
 
if isempty(opt.refDepth) && isfield(G.cells, 'centroids')
    opt.refDepth = G.cells.centroids(tmp.cell(1), 3);
end
% multiply segments by weight (non-unit for segments shared by multiple cells)
seg = bsxfun(@times, tmp.vec, tmp.weight);
W = addWell(W, G, rock, tmp.cell, 'lineSegments', seg, 'refDepth', opt.refDepth, other{:});
end
