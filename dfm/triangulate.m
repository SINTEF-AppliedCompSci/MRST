function g = triangulate(vertices,constraints,tags)
% Make a Delaunay Triangulation and return a Mrst grid
%
% SYNOPSIS
% g = triangulate(vertices,edges)
%
% PARAMETERS
%   vertices    - the vertices in the triangulation
%   constraint  - constraint to the triangulation
%  Optinal
%   tags        - one tag for each constraint
%
%
%  RETURNS
%
%   g           - grid structure for mrst as described in grid_structre
%
% Copyright 2011-2012 University of Bergen, 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.


% Create a constrained Delaunay triangulation
% on the point cloud, and create the grid
delTri     = DelaunayTri(vertices, constraints);
g     = triangleGrid(delTri.X, delTri.Triangulation);

% store the nodes of the triangles
g.cells.nodes = reshape(delTri.Triangulation',[],1);
g.cells.nodePos = (1:3:numel(g.cells.nodes)+1)';

% compute the geometry of the grid
g = computeGeometry(g);

% Link the constraints from the Delaunay triangulation  to the tags explicitly,
% in case these have been changed (it shouldn't happen, and will give a warning)
if nargin > 2
    % store the tags in the grid structure.
    g.faces.tags = zeros(g.faces.num,size(tags,2));

    % Recover faceNodes, and sort the nodes column wise (for the use of ismember
    % later)
    faceNodes = sort(reshape(g.faces.nodes,2,[])',2);

    % Use constraints from the Delaunay triangulation, in case these have been
    % changed (it shouldn't happen, and will give a warning)
    % Rowwise sort to prepare for using ismember
    constraints = sort(delTri.Constraints,2);

    % Find the faces linked to the constraint.
    [exists,fracFace] = ismember(constraints,faceNodes,'rows');

    % make sure we have found all constraints
    assert(all(exists))
    % Assign a tag to the newly found fracture face
    g.faces.tags(fracFace,:) = tags;

end


