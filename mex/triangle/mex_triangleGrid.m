function G = mex_triangleGrid(points, edges, varargin)
%Construct 2d triangle grid in physical space.
%
% SYNOPSIS:
%   G = mex_triangleGrid(edgepointlist, edgelist)
%   G = mex_triangleGrid(edgepointlist, edgelist, 'pn1', pv1, ...)
%
% PARAMETERS:
%   pointlist - List of node coordinates for edges.
%
%   edgelist  - List of edges defining the boundary of the domain that is
%               to be gridded.
%  'pn'/pv    - List of 'key'/value pairs defining optional parameters.
%               The supported options are:
%
%              maxArea  -- Maximal area of triangle.
%
%              minAngle -- Minimum angle in triangles. Use sensible values
%                          here.  Otherwise, the mex function may not
%                          return.
%              verbose  -- Whether or not to display progress information
%                          Logical.  Default value: Verbose = false.
%
% RETURNS:
%   G - Grid structure mostly as detailed in grid_structure, though lacking
%       the fields
%         - G.cells.volumes
%         - G.cells.centroids
%
%         - G.faces.areas
%         - G.faces.normals
%         - G.faces.centroids
%
%       These fields may be computed using the function computeGeometry.
%
% EXAMPLE:
%   % Make a 10m-by-5m grid.
%      points = [0,0; 5,0; 5,10; 0,10];
%      edges  = [1,2;2,3;3,4;4,1]
%      G = mex_triangleGrid(points, edges, 'maxArea', 0.3);
%
%   % Plot the grid in 3D-view.
%      f = plotGrid(G); axis equal tight; view(2);
%
% SEE ALSO:
%   grid_structure, computeGeometry.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


opt = struct('maxArea',  10, ...
             'minAngle', 32, ...
             'segmentmarkers', [], ...
             'verbose', false);
opt = merge_options(opt, varargin{:});
if isempty(opt.segmentmarkers),
   opt.segmentmarkers = int32(1:size(edges,1));
end
assert(numel(opt.segmentmarkers) == size(edges,1));

assert(size(points, 2)==2);
assert(size(edges, 2)==2);
assert(min(edges(:))>=0 );
assert(max(edges(:))<=size(points,1));
assert(opt.maxArea  > 0 );
assert(opt.minAngle > 0 );
assert(opt.minAngle < 40);

%% domain boundary
in.pointlist   = points';
in.segmentlist = int32(edges'-1);
in.segmentmarkerlist = opt.segmentmarkers;

options = sprintf('Q2pzvAneiq%fa%f', opt.minAngle, opt.maxArea);
if opt.verbose,
   options = [options(2:end), 'V'];
end

%% Construct triangulation
[out, vorout]=mex_triangle(in, options);

G.nodes.num       = size(out.pointlist, 2);
G.nodes.coords    = out.pointlist';
G.faces.num       = size(out.edgelist, 2);
G.faces.numNodes  = int32(repmat(2, [G.faces.num, 1]));
G.faces.neighbors = int32(vorout.edgelist');
G.cells.num       = size(out.trianglelist, 2);
G.cells.numFaces  = int32(repmat(3, [G.cells.num, 1]));
G.faces.nodes       = int32(out.edgelist(:));
G.faces.tag = out.edgemarkerlist;

% Computing cellFaces from triangle list requires some kind of sorting.
%{
e = double(out.edgelist);
M = sparse(e(1,:), e(2,:), 1:size(e, 2));M = M + M';
T = reshape(out.trianglelist([1,2,2,3,3,1], :), 2, [])';
A = int64(T(:,1)); B = int64(T(:,2));
G.cells.faces = int64(full(M(A+(B-1)*G.nodes.num) ));
%}
e = (1:size(out.edgelist, 2))';
n = double(vorout.edgelist)';
n = sortrows([n(:,1), e; n(:,2), e]);
G.cells.faces = n(n(:,1)~=0,2);


% PAtch to convert numFaces and numNodes to facePos and nodePos;
pos = @(n) int32(cumsum([1; double(reshape(n, [], 1))]));
G.faces.nodePos = pos(G.faces.numNodes);
G.cells.facePos = pos(G.cells.numFaces);
G.cells = rmfield(G.cells, 'numFaces');
G.faces = rmfield(G.faces, 'numNodes');

G.griddim = 2;
G.type = { mfilename };
