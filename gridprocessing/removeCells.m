function [G, cellmap, facemap, nodemap] = removeCells(G, cells)
%Remove cells from grid and renumber cells, faces and nodes.
%
% SYNOPSIS:
%   H = removeCells(G, cells);
%   [H, cellmap, facemap, nodemap] = removeCells(G, cells)
%
% PARAMETERS:
%   G       - Valid grid definition
%
%   cells   - List of cell numbers (cell IDs) to remove.
%
% RETURNS:
%   H       - Updated grid definition where cells have been removed.
%             Moreover, any unreferenced faces or nodes are subsequently
%             removed.
%
%   cellmap - Cell numbers in `G` for each cell in `H`.  Specifically,
%             `cellmap(i)` is the cell ID of `G` that corresponds to cell
%             `i` in `H`.
%
%   facemap - Face numbers in `G` for each face in `H`.  Specifically,
%             `facemap(i)` is the face ID of `G` that corresponds to face
%             `i` in `H`.
%
%   nodemap - Node numbers in `G` for each node in `H`.  Specifically,
%             `nodemap(i)` is the node ID of `G` that corresponds to node
%             `i` in `H`.
%
% EXAMPLE:
%   G = cartGrid([3, 5, 7]);
%   G = removeCells(G, 1 : 2 : G.cells.num);
%   plotGrid(G), view(-35, 20), camlight
%
% NOTE:
%   The process of removing cells is irreversible.
%
% SEE ALSO:
%   `readGRDECL`, `makeInternalBoundary`

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

% Mark resulting outer boundary with optional tag.
%
%
% If we make tag,
%  m1 = ismember(G.faces.neighbors(:,1), cells);
%  m2 = ismember(G.faces.neighbors(:,2), cells);
%  t  = find((m1 & ~m2) | (~m1 & m2));
%  G.faces.tag(t)=99;

   if isempty(cells)
      % Function is no-op if no cells are scheduled for removal.
      %
      cellmap = (1 : G.cells.num) .';
      facemap = (1 : G.faces.num) .';
      nodemap = (1 : G.nodes.num) .';

      return;
   end

   if islogical(cells)
      assert (numel(cells) == G.cells.num, ...
             ['The length of logical vector differs from ', ...
              'number of grid cells.']);
      cells = find(cells);
   end

  % New numbering of cells
  ind        = false(G.cells.num,1);
  ind(cells) = true;
  cellmap    = mapExcluding(ind);

  % remove and renumber cells in cellFaces
  if isfield(G.cells, 'numFaces')
     warning('MRST:deprecated', ...
            ['The field numFaces will not be supported in ', ...
             'future releases of MRST.'])
     % 'numFaces' exists.  Preserve.
     numFaces                = G.cells.numFaces;
     G.cells.numFaces(cells) = [];
  else
     numFaces = diff(G.cells.facePos);
  end
  G.cells.faces(rldecode(ind, numFaces), :) = [];

  % Alter cell numbering in G.faces.neighbors
  n = G.faces.neighbors;
  G.faces.neighbors(n(:,1)>0,1) = cellmap(n(n(:,1)>0,1));
  G.faces.neighbors(n(:,2)>0,2) = cellmap(n(n(:,2)>0,2));
  clear n

  % Alter cells
  numFaces(cells) = [];
  G.cells.num     = G.cells.num-numel(cells);
  G.cells.facePos = cumsum([1; double(numFaces)]);
  if isfield(G.cells, 'indexMap')
     G.cells.indexMap(cells) = [];
  end

  % new numbering of faces.
  ind     = all(G.faces.neighbors(:,1:2)==0,2);
  facemap = mapExcluding(ind);

  % remove and renumber faces in faceNodes
  if isfield(G.faces, 'nodes')
     if isfield(G.faces, 'numNodes')
        warning('MRST:deprecated', ...
               ['The field ''numNodes'' will be removed in a', ...
                'future release of MRST.'])
        % 'numNodes' exists.  Preserve.
        numNodes              = G.faces.numNodes;
        G.faces.numNodes(ind) = [];
     else
        numNodes = diff(G.faces.nodePos);
     end
     G.faces.nodes(rldecode(ind, numNodes)) = [];
  end

  % remove and renumber faces in cellFaces
  G.cells.faces(:,1) = facemap(G.cells.faces(:,1));

  if any(G.cells.faces(:,1) == 0)
      error('In removeCells: Too many faces removed!');
  end

  assert (isfield(G, 'type'), ...
          'Every grid must record its origin in the field ''type''.');

  % Alter geometry
  if any(strcmp(G.type, 'computeGeometry')) || ...
        any(strcmp(G.type, 'mcomputeGeometry'))
      G.cells.centroids(cells,:) = [];
      G.cells.volumes(cells,:)   = [];
      G.faces.areas(ind)         = [];
      G.faces.centroids(ind,:)   = [];
      G.faces.normals(ind,:)     = [];
  end
  
  G.faces.neighbors(ind,:) = [];

  if isfield(G.faces, 'nodes')
     numNodes(ind)   = [];
     G.faces.nodePos = cumsum([1; double(numNodes)]);
  end

  if isfield(G.faces, 'tag')
     G.faces.tag(ind,:) = [];
  end

  G.faces.num = G.faces.num - sum(ind);

  % Construct node map:
  if isfield(G.faces, 'nodes')
     ind = true(G.nodes.num, 1);
     ind(G.faces.nodes) = false;
     nodemap = mapExcluding(ind);
     
     % Remove nodes
     G.nodes.coords(ind,:) = [];
     G.nodes.num           = G.nodes.num - sum(ind);
     G.faces.nodes         = nodemap(G.faces.nodes);

     if any(G.faces.nodes == 0)
        error('In removeCells: Too many nodes removed!');
     end
  else
     nodemap = [];
  end

  if isfield(G.cells, 'cpnodes')
     pick = cellmap ~= 0;
     G.cells.cpnodes = reshape(nodemap(G.cells.cpnodes(pick, :)), ...
                               sum(pick), []);

     if any(any(G.cells.cpnodes == 0))
        error('CpNodes:ExcessiveNodeRemoval', ...
             ['Too many nodes removed to preserve corner-point ', ...
              'vertex structure']);
     end
  end

  G.type = [G.type, { mfilename }];

  cellmap = find(cellmap > 0);
  facemap = find(facemap > 0);
  nodemap = find(nodemap > 0);
end

%--------------------------------------------------------------------------

function m = mapExcluding(indices)
   n            = numel(indices);
   ind          = ones(n,1);
   ind(indices) = 0;
   m            = cumsum(ind);
   m(indices)   = 0;
end
