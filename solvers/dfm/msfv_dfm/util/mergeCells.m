function [g,remove_faces] = mergeCells(g,remove_faces)
% merge two cells separated by the first face in remove_faces
% return the updated remove_faces list as well as a grid
% with one less number of cells.
%
% SYNOPSIS
%     [g,remove_faces] = mergeCells(g,remove_faces)
%
% PARAMETERS
%       g - mrst grid structure
%       remove_faces - a list of faces to be removed
%                       NB Only the first face is removed
%
% OUTPUT
%       g - mrst grid structure with one less cell
%       remove-faces - updated list of faces (- the removed faces)
%
% Copyright 2011-2012 University of Bergen, 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.

% picks the face to be merged.
ind = false(g.faces.num,1);

% we merge two cells at the time starting at the first face in remove_faces
ind(remove_faces(1)) = true;

% the associated cells
cells = g.faces.neighbors(ind,:);
cells = unique(cells)';

% store the face and nodes before removing the cells in
% order to reconstruct the new cell
faceind = rldecode(1:g.cells.num,diff(g.cells.facePos),2)';
nodeind = rldecode(1:g.cells.num,diff(g.cells.nodePos),2)';
faces =  g.cells.faces;
nodes = unique(g.cells.nodes(ismember(nodeind,cells),:));
n = g.faces.neighbors;

% remove the old cells
[g, cellmap, facemap, nodemap] = removeCells_mod(g, cells,ind);

n = n(facemap>0,:);
cellmap(cellmap==0) = g.cells.num+1;
g.faces.neighbors(n(:,1)>0,1) = cellmap(n(n(:,1)>0,1));
g.faces.neighbors(n(:,2)>0,2) = cellmap(n(n(:,2)>0,2));
clear n

% update the node and face indices
remove_faces = facemap(remove_faces);
remove_faces = remove_faces(remove_faces>0);
faces = facemap(faces);
faceind = faceind(faces>0);
faces = faces(faces>0);
faces_this = faces(ismember(faceind,cells));
numFaces = numel(faces_this);

% create new cell
g.cells.faces = g.cells.faces(:,1);
g.cells.faces = [g.cells.faces; faces_this];
g.cells.facePos = [g.cells.facePos; g.cells.facePos(end) + numFaces];
g.cells.num = g.cells.num + 1;

% calculate volume and centroids of the new cell
volume = sum(g.cells.volumes(cells));
centroid = mean(g.cells.centroids(cells,:),1);
nodes = setdiff(nodemap(nodes),0);

% remove old cells
g.cells.volumes(cells,:) = [];
g.cells.centroids(cells,:) = [];

% add the new cell
g.cells.volumes = [g.cells.volumes; volume];
g.cells.centroids = [g.cells.centroids ; centroid];
g.cells.nodes = [g.cells.nodes; nodes];
g.cells.nodePos = [g.cells.nodePos; g.cells.nodePos(end) + numel(nodes)];

end


