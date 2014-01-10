function dual = partitionDualDFM(g)
% Create the dual grid structure as used by MSFVM
%
% SYNOPSIS:
%       dual = partitionDualDFM(g)
%
% PARAMETERS:
%       g - mrst grid structure with additional fields as created
%           by createGridHierarchy.m
%
% OUTPUT:
%       dual - the dual grid structure used by MSFVM
%
% Copyright 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.


% the node cells are represented by hybrid cells of kind 2
% ie. as points in 2d
ncells = find(g.cells.hybrid == 2);

% get the node cells associated with large-scale fracture intersections.
% npick is the intersection between large-scale and large-scale fractures
% while nepick between large-scale and small-scale fractures
% they thus belong among the edge cells
[npick,nepick] = getCoarseNodes(g);

% only hybrid cells of 2nd kind associated with fracture intersection
% between large-scale fractures are classified as node cells
dual.nn = ncells(npick);

% edge cells are all hybrid cells of kind 1
ecells = find(g.cells.hybrid == 1);

% but excluding the fine-scale fractures
epick = g.faces.tags(g.cells.tags(g.cells.hybrid == 1),g.faces.tagmark.map2c)>0;

% the edge cells
dual.ee = [ecells(epick); ncells(nepick)];

% the rest of the cells are inner cells
dual.ii = [find(g.cells.hybrid == 0);ecells(~epick);ncells(~npick&~nepick)];

% make sure we don't lose some cells
assert(numel(dual.nn)+numel(dual.ee)+numel(dual.ii) == g.cells.num)

end


function [n,e] = getCoarseNodes(g)
% identify intersection between large-scale fractures
% and between large and small-scale fractures.


% pick the coarse hybrid cells
isCoarseFace = g.faces.tags(:,g.faces.tagmark.map2c) > 0;

% make sure there is no hybrid face among the faces
assert(all(g.faces.hybrid(isCoarseFace) == 0))

% then we know an edge is between two nodes
edges = [g.faces.nodes(g.faces.nodePos(isCoarseFace)),g.faces.nodes(g.faces.nodePos(isCoarseFace)+1)];

%  edges are duplicated as they are found on both sides
% of the hybrid cells
edges = unique(edges,'rows');

% count number of reference to each node
numberOfNodes = accumarray(edges(:),1);

% check whether the nodes lies on the boundaries.
% 0: interior 1: edge 2: corner
nodeTypes = zeros(numel(numberOfNodes),1);
nodeTypes(numberOfNodes==2) = nodeType(g,find(numberOfNodes==2));

% the intersection of the large-scale fractures should be referred to
% by more than two edges or two edges if we are the boundry
pick_loc_n = numberOfNodes>2 | (numberOfNodes==2 & nodeTypes>0);

% the intersection between the large-scale and small-scale fractures
% should be referred to by two edges.
pick_loc_e = numberOfNodes==2 & nodeTypes==0;

% map the nodes to fine-scale cells
cnodes = g.cells.tags(g.cells.hybrid == 2);
n = ismember(cnodes,find(pick_loc_n));
e = ismember(cnodes,find(pick_loc_e));

end

