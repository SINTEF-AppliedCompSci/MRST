function [fcells,ncells] = getFractureCells(g,tag)
% Return the fracture cells related to a given tag
%
% SYNOPSIS:
% [fcells,ncells] = getFractureCells(g,tag)
%
% PARAMETERS:
%       g - MRST grid structure
%       tag - a list of tags
%
% OUTPUT:
%       fcells - fracture cells
%       ncells - cells in the intersection of the fractures
%
% Copyright 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.



% get the hybrid cells of kind 1 associated with the given tag
hyb1cells = find(g.cells.hybrid == 1);
ftag = g.faces.tags(g.cells.tags(g.cells.hybrid == 1),g.faces.tagmark.frac);
fcells = hyb1cells(ismember(ftag,tag));

% get the hybrid cells of kind 2 associated with the given tag
% i.e the intersection between fractures

isFracFace = ismember(g.faces.tags(:,g.faces.tagmark.frac),tag);

% make sure there is no hybrid face among the fracture face
assert(all(g.faces.hybrid(isFracFace) == 0))

% then we know an edge is between two nodes
edges = [g.faces.nodes(g.faces.nodePos(isFracFace)),g.faces.nodes(g.faces.nodePos(isFracFace)+1)];

% fracture edges are duplicated as they are found on both sides
% of the hybrid fracture cell
edges = unique(edges,'rows');

% count the number of time a node is referred
numberOfNodes = accumarray(edges(:),1);

% find the hybrid cells of second kind and their corresponding nodes
hyb2cells = find(g.cells.hybrid == 2);
cnodes = g.cells.tags(g.cells.hybrid == 2);

% pick out the node cells (cells in the intersection of the fractures)
ncells = hyb2cells(ismember(cnodes,find(numberOfNodes > 0)));


