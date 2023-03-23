function g = makeCoarseDualGrid(vertices,fractures,box,merge,precision)
% Construct the coarse grid from constrains given by
% the fractures
%
% SYNOPSIS:
% g = makeCoarseDualGrid(vertices,fractures,opt)
%
% PARAMETERS
%       vertices  - coordinates of the nodes
%       fractures - the two first columns gives the end nodes
%                   the rest of the columns are tags
%       box       - the bounding box
%       merge     - Boolean telling whether to merge cells or not
%       precision - the relative precision of the grid
%
% OUTPUT
%       g - the coarse dual grid with MRST grid structure
%
% Copyright 2011-2012 University of Bergen, 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.


% add the boundaries to the list of constraints
[vertices, edges,tagmark] = addBondaries(box,vertices,fractures);
clear fractures

% remove duplicated edges
[vertices,edges] = removeDuplicates(vertices,edges);

% create new nodes in the intersection of fractures and
% chop the fractures accordingly
[vertices, edges] = removeFractureIntersections (vertices, edges, box, struct('precision', precision));

% triangulate the domain.
g = triangulate(vertices,edges(:,1:2),edges(:,3:end));

% tags are stored in a matrix tagmarks gives an interpretation of
% the columns of the matrix.
g.faces.tagmark = tagmark;
clear edges vertices

% if true all faces with tags == 0 is removed and the corresponding
% cells are merged.
if merge
    g = mergeCells_loc(g,box);
end

end

function [vertices,edges,tagmark] = addBondaries(box,vertices,fractures)
% we add the boundary as constraints and tag them

corners = [box(1) box(3) ; box(2) box(3) ; box(2) box(4) ; box(1) box(4)];
nVertices = size(vertices,1);
vertices = [vertices; corners];

% Add the boundary. These will have index nbandPoints + 1 : nbandPoints + 4
edges = [fractures(:,1:2) ; nVertices + 1 , nVertices + 2 ; ...
    nVertices + 2 , nVertices + 3 ; ...
    nVertices + 3 , nVertices + 4 ; ...
    nVertices + 4 , nVertices + 1 ];

% Let the boundary be marked according to
%    -1 = SOUTH, -2 = EAST, -3 = NORTH, -4 = WEST
markers = [ zeros(size(fractures,1),1); - (1:4)'];
markers2 = [fractures(:,3) ; zeros(4,1)];
edges = [edges , markers2,markers];

% column one is for the fracture tags and column 2 for the boundary
tagmark.frac = 1;
tagmark.boundary = 2;
end


function g = mergeCells_loc(g,box)
% merge cells for the moment all faces tagged zero are removed
% more advanced strategies may be considered later

% remove all faces with tag zero
remove_faces = find(all(g.faces.tags==0,2) & all(g.faces.neighbors~=0,2));

% we don't want to remove the faces that connected to fractures that ends
% within the domain
endFaces = getEndingFaces(g,box);
remove_faces = setdiff(remove_faces, endFaces);

% remove the face one by one.
n = numel(remove_faces);
count = 0;
while count < n
    [g,remove_faces] = mergeCells(g,remove_faces);
    if isempty(remove_faces)
        return
    end
    count = count + 1;

end

% after cells are merged, faces must be merged to avoid
% cases where two faces has the same neighbors.
g = mergeFaces(g);

end

function endFaces = getEndingFaces(g,box)
% return the face index of faces connected to fractures that ends in the
% interior of the domain

% maps nodes to faces
nodeind = rldecode(1:g.faces.num,diff(g.faces.nodePos),2)';

% find the number of fracture segments meeting at an vertex
eta = accumarray(g.faces.nodes(ismember(nodeind,find(g.faces.tags(:,g.faces.tagmark.frac)>0))),1,[g.nodes.num,1]);

% we are not interested in fractures ending at the boundary
coords = g.nodes.coords;
isb = coords(:,1)==box(1) | coords(:,1) == box(2) |coords(:,2)==box(3) | coords(:,2) == box(4);

% pick the fracture nodes at the interior
endNodes = find(eta <= 2 & ~isb);

% find the corresponding faces
endFaces = nodeind(ismember(g.faces.nodes,endNodes));
end


