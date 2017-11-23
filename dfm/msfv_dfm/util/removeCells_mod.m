function [G, cellmap, facemap, nodemap] = removeCells_mod(G, cells,faces)
%Remove cells from grid and renumber cells, faces and nodes.
%
% THIS FILE IS MODIFIED FROM THE ORIGINAL removeCells.m
% The following field are updated from the face field
% - tags, tags2, areas, normals, fracNodes, centroids
% Portions Copyright 2011-2012 University of Bergen.
%
%
% SYNOPSIS:
%   G = removeCells(G, cells);
%   [G, cellmap, facemap, nodemap] = removeCells(G, cells)
%
% PARAMETERS:
%   G          - Valid grid definition
%
%   cells      - list of cell numbers to be removed.
%
%
% RETURNS:
%   G          - Updated grid definition where cells have been removed.
%                In addition, any unreferenced faces or nodes are
%                subsequently removed.
%
% EXAMPLE:
%
%      G = cartGrid([3,5,7]);
%      G = removeCells(G, (1:2:G.cells.num));
%      plotGrid(G);view(-35,20);camlight
%
%
%
% NOTE:
%
%   The process of removing cells is irreversible.
%
% SEE ALSO:
%   `readGRDECL`, `deactivateZeroPoro`.

%{
Copyright 2009, 2010, 2011 SINTEF ICT, Applied Mathematics.

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


   if isempty(cells),
      % Function is no-op if no cells are scheduled for removal.
      %
      cellmap = (1 : G.cells.num) .';
      facemap = (1 : G.faces.num) .';
      nodemap = (1 : G.nodes.num) .';

      return;
   end


  % New numbering of cells
  ind        = false(G.cells.num,1);
  ind(cells) = true;
  cellmap    = mapExcluding(ind);

  % remove and renumber cells in cellFaces
  if isfield(G.cells, 'numFaces'),
     % 'numFaces' exists.  Preserve.
     numFaces                = G.cells.numFaces;
     G.cells.numFaces(cells) = [];
  else
     numFaces = diff(G.cells.facePos);
  end
  G.cells.faces(rldecode(ind, numFaces), :) = [];

  if isfield(G.cells,'nodes')
      numNodes = diff(G.cells.nodePos);
      G.cells.nodes(rldecode(ind, numNodes),:) = [];
      numNodes(cells) = [];
      G.cells.nodePos  = int32(cumsum([1; double(numNodes)]));
  end

  % Alter cell numbering in G.faces.neighbors
  n = G.faces.neighbors;
  G.faces.neighbors(n(:,1)>0,1) = cellmap(n(n(:,1)>0,1));
  G.faces.neighbors(n(:,2)>0,2) = cellmap(n(n(:,2)>0,2));
  clear n

  % Alter cells
  numFaces(cells)         = [];
  G.cells.num             = G.cells.num-numel(cells);
  G.cells.facePos         = int32(cumsum([1; double(numFaces)]));
  if isfield(G.cells, 'indexMap'), G.cells.indexMap(cells) = []; end

  %new numbering of faces.
  ind     = all(G.faces.neighbors(:,1:2)==0,2);
  if ~isempty(faces)
      ind = faces;
  end
  facemap = mapExcluding(ind);

  % remove and renumber faces in faceNodes
  if isfield(G.faces, 'numNodes'),
     % 'numNodes' exists.  Preserve.
     numNodes              = G.faces.numNodes;
     G.faces.numNodes(ind) = [];
  else
     numNodes = diff(G.faces.nodePos);
  end
  G.faces.nodes(rldecode(ind, numNodes)) = [];

  % remove and renumber faces in cellFaces
  G.cells.faces(:,1) = facemap(G.cells.faces(:,1));


  if any(G.cells.faces(:,1)==0),
      error('In removeCells: Too many faces removed!');
  end


  numNodes(ind)            = [];
  G.faces.neighbors(ind,:) = [];
  G.faces.nodePos = int32(cumsum([1; double(numNodes)]));
  if isfield(G.faces, 'tag'),
      G.faces.tag      (ind,:) = [];
  end
  if isfield(G.faces, 'tags'),
      G.faces.tags      (ind,:) = [];
  end
  if isfield(G.faces, 'tags2'),
      G.faces.tags2      (ind,:) = [];
  end
  if isfield(G.faces, 'centroids'),
      G.faces.centroids      (ind,:) = [];
  end
  if isfield(G.faces, 'areas'),
      G.faces.areas      (ind,:) = [];
  end
  if isfield(G.faces, 'fracNodes'),
      G.faces.fracNodes      (ind,:) = [];
  end
  if isfield(G.faces, 'normals'),
      G.faces.normals      (ind,:) = [];
  end


  G.faces.num              = G.faces.num - sum(ind);



  % Construct node map:
  ind = true(G.nodes.num, 1);
  ind(G.faces.nodes) = false;
  nodemap = mapExcluding(ind);

  % Remove nodes
  G.nodes.coords(ind,:) = [];
  G.nodes.num           = G.nodes.num - sum(ind);
  G.faces.nodes           = int32(nodemap(G.faces.nodes));

  if any(G.faces.nodes==0),
      error('In removeCells: Too many nodes removed!');
  end

  if isfield(G, 'type'),
     G.type = [G.type, mfilename];
  end

  if isfield(G.cells,'nodes')
    G.cells.nodes = int32(nodemap(G.cells.nodes));
  end
%  G.cells.actnum(cells)=false;

  function m = mapExcluding(indices)
      n            = numel(indices);
      ind          = ones(n,1);
      ind(indices) = 0;
      m            = cumsum(ind);
      m(indices)   = 0;
