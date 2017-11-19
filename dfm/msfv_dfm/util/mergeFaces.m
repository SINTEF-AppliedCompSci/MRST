function g = mergeFaces(g)
% merge faces that has the same pair of neighbors
% for standard grid this do not happen
% but it may occur when cells is merged using mergeCells.m
%
% SYNAPSIS:
% g = mergeFaces(g)
%
% PARAMETERS
%       g - mrst grid structure
%
% OUTPUT
%       g - mrst grid structure where faces with the same pair
%           of neighbors are removed
%
% Copyright 2011-2012 University of Bergen, 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.



% for the moment we are only interested in merging inner faces
isInner = all(g.faces.neighbors~=0,2);

% find one of the faces that shares both neighbors
[~,I] = unique(sort(g.faces.neighbors(isInner,:),2),'rows');
df = find(isInner);
df = setdiff(df,df(I));

% return if no such faces exists
if isempty(df)
    return
end

% find the faces and a link between the faces to be merged
[tf,loc] = ismember(sort(g.faces.neighbors,2),sort(g.faces.neighbors(df,:),2),'rows');

% faces that needs merging
faces = find(tf);

% a map that links the faces
map = loc(tf);

% pick out the associated nodes
nodes = g.faces.nodes(mcolon(g.faces.nodePos(faces),g.faces.nodePos(faces)+1));

% make a node map
nodes_ind = rldecode(map,2);

% store the neighbors
neigh = sort(g.faces.neighbors(tf,:),2);

% store geometry of the faces and remove
% them from the structure.
areas = g.faces.areas(tf,:);
g.faces.areas(tf) = [];
normals = g.faces.normals(tf,:);
g.faces.normals(tf,:) = [];
centroids = g.faces.centroids(tf,:);
g.faces.centroids(tf,:) = [];
tags = g.faces.tags(tf,:);
g.faces.tags(tf,:) = [];
if isfield(g.faces,'fracNodes')
    g.faces.fracNodes(tf,:) = [];
end
% we first remove the faces and then rebuild the new face
g = removeFaces(g,faces);

% build the new face
numfaces = g.faces.num;
numfaces_new = numel(df);
facenodes = [];
neigh_this = zeros(numfaces_new,2);
centroids_new = zeros(numfaces_new,2);
normals_new = zeros(numfaces_new,2);
tags_new = zeros(numfaces_new,size(tags,1));
fracNodes_new = zeros(numfaces_new,2);

for i = 1:numfaces_new;
    % we remove the duplicated nodes ( the one lying between the faces)
    nodes_this = nodes(nodes_ind==i);
    [nodes_this,n] = repcount(nodes_this);
    nodes_this = nodes_this(n==1);

    % update the facenodes
    facenodes = [facenodes;nodes_this];

    % pick out the neighbors
    neigh_this(i,:) = unique(neigh(map==i,:),'rows');

    % get the new geometry, this is ad-hoc for the moment
    % a more correct way is to recalculate the geometry
    % for this new face
    centroids_new(i,:) = mean(centroids(map==i,:));
    normals_new(i,:) = mean(normals(map==i,:));
    tags_new(i,:) = max(tags(map==i,:));

    % store the new face nodes
    fracNodes_new(i,:) = nodes_this';
end

% calculate the new areas
areas_new = accumarray(map,areas);

% in 2d all faces has 2 nodes
numnodes = 2 * ones(numfaces_new,1);

% add the new faces to the structure
g = addFaces(g,facenodes,numnodes,neigh_this);

% remove the unused nodes from the grid
removenodes = setdiff(nodes,facenodes);
nodemap = zeros(g.nodes.num,1);
keepnodes = true(g.nodes.num,1);
keepnodes(removenodes) = false;
nodemap(keepnodes) = (1:g.nodes.num-numel(removenodes))';
g.nodes.coords(removenodes,:) = [];
g.nodes.num = g.nodes.num-numel(removenodes);

% update the face and cell nodes with new node index
g.faces.nodes = nodemap(g.faces.nodes);
[g.cells.nodes,g.cells.nodePos] = removeFromPackedData(g.cells.nodePos,g.cells.nodes,removenodes);
g.cells.nodes = nodemap(g.cells.nodes);

% update the geometry of the structure
g.faces.areas = [g.faces.areas;areas_new];
g.faces.centroids = [g.faces.centroids;centroids_new];
g.faces.tags = [g.faces.tags;tags_new];
g.faces.normals = [g.faces.normals;normals_new];
if isfield(g.faces,'fracNodes')
    g.faces.fracNodes = [g.faces.fracNodes;fracNodes_new];
end

% update the cells with the new faces
cellind = neigh_this(:);
faceind = repmat((1:numfaces_new)',2,1)+numfaces;
cellfaces = g.cells.faces;
facePos = g.cells.facePos;
[cellfaces,facePos] = insertInPackedData(facePos,cellfaces,cellind,faceind);
g.cells.faces = cellfaces;
g.cells.facePos = facePos;
