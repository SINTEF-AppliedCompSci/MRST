function i = testNormals(G)
% Tests if the normals points from neighbor 1 to 2 and returns the index of
% the normals that points in the opposite direction.
%
% SYNOPSIS:
%   i = testNormals(G)
%
% PARAMETERS:
%   G       - Grid data structure.
%   i       - Index of the normals pointing in the wrong direction. Empty
%             if all the normals points correctly.
%
% Copyright 2011-2012 University of Bergen
%
% This file is licensed under the GNU General Public License v3.0.

i1 = G.faces.neighbors(:,1)~=0;

% vector pointing from center of cell 1 to face center
r1 = G.faces.centroids(i1,:) - G.cells.centroids(G.faces.neighbors(i1,1),:);
i2 = G.faces.neighbors(:,2)~=0;
% vector pointing from center of cell 2 to face center
r2 = G.faces.centroids(i2,:)-G.cells.centroids(G.faces.neighbors(i2,2),:);

% normals should point from neighbor 1 to 2
ii1 = find(dot(G.faces.normals(i1,:),r1,2)<0);
ii2 = find(dot(G.faces.normals(i2,:),r2,2)<0);

i=[ii1 ; ii2];





