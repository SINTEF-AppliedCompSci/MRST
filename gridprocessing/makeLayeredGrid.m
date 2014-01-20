function G = makeLayeredGrid (G, nlayers)
%Extrude 2D grid to layered 3D grid with n layers.
%
% SYNOPSIS:
%   G = makeLayeredGrid(G, n)
%
% PARAMETERS:
%   G      - Valid 2D areal grid.
%
%   n      - Number of layers.  Each layer is 1 meter thick.  (TODO HERE!)
%
%  'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%            supported options are:
%
%              Verbose -- Whether or not to display progress information
%                         Logical.  Default value: Verbose = false.
%
% RETURNS:
%   G      - Valid 3D grid as described in grid_structure.
%
% EXAMPLE:
%   G = makeLayeredGrid(cartGrid([2,2]), 3);
%   plotGrid(G); view(10,45);
%
% SEE ALSO:
%   grid_structure.

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


% Build layered grid from 2D mesh.

assert(G.griddim == 2, ...
   'Areal grid must be 2D.  Field ''griddim'' says otherwise...');
G.type = [G.type, { mfilename }];

% Faces with horizontal normal
%-----------------------------
edges = reshape(G.faces.nodes(:,1), 2, []) .';
dz    = G.nodes.num;
p1    = repmat(edges(:,1), [nlayers,1])+...
        kron((0:nlayers-1)', ones(size(edges,1), 1))*dz;
p2    = repmat(edges(:,2), [nlayers,1])+...
        kron((0:nlayers-1)', ones(size(edges,1), 1))*dz;
hFaces = [p1, p2, p2+dz, p1+dz];
hNumNodes = repmat(4, [size(hFaces, 1), 1]);

% Neighbor relations in horizontal directions
nbrs       = double(G.faces.neighbors);
bdry       = repmat(nbrs==0, [nlayers,1]);
% Find all h-neighbors
hNeighbors = repmat(nbrs, [nlayers, 1]) + ...
             kron((0:nlayers-1)'*G.cells.num, ones(size(nbrs, 1), 2));
% Make sure outer boundary is still there.
hNeighbors(bdry)=0;



% Faces with vertical normal
%----------------------------
cn     = getCellNodes(G);
vFaces = repmat(cn, [nlayers+1, 1]) + ...
         kron((0:nlayers)', ones(size(cn,1), 1))*dz;
vNumNodes = repmat(diff(G.cells.facePos), [nlayers+1, 1]);

% Neighbor relations in vertical direction.
vNeighbors = zeros((nlayers+1)*G.cells.num,2);
cells      = 1:nlayers*G.cells.num;
vNeighbors(cells,       2) = cells;
vNeighbors(cells+G.cells.num,1) = cells;

% Build grid structure
%---------------------------
G.nodes.coords    = [repmat(G.nodes.coords, [nlayers+1,1]), ...
                    kron((0:nlayers)', ones(G.nodes.num,1))];
G.nodes.num       = (nlayers+1)*G.nodes.num;
G.cells.num       = nlayers*G.cells.num;
G.cells.numFaces  = repmat(diff(G.cells.facePos)+2, [nlayers, 1]);
G.faces.neighbors = [vNeighbors; hNeighbors];
G.faces.num       = size(G.faces.neighbors,1);
G.faces.numNodes  = [vNumNodes; hNumNodes];
G.faces.nodes       = [vFaces(:); reshape(hFaces',[], 1)];
pos = @(n) cumsum([1; double(reshape(n, [], 1))]);
G.faces.nodePos = pos(G.faces.numNodes);
G.cells.facePos = pos(G.cells.numFaces);
% Cell topology from G.faces.neighbors
N   = G.faces.neighbors;
tmp = sortrows([N(:,1),(1:G.faces.num)', ones(G.faces.num,1);...
                N(:,2),(1:G.faces.num)',-ones(G.faces.num,1)]);

G.cells.faces       = tmp(tmp(:,1)>0,2);
G.cells.faces       = [G.cells.faces, zeros(size(G.cells.faces))];
G.cells.indexMap    = (1 : G.cells.num) .';

G.numLayers = nlayers;
G.layerSize = G.cells.num / G.numLayers;
if isfield(G, 'cartDims'),
   G = rmfield(G, 'cartDims');
end
G.griddim = 3;
%--------------------------------------------------------------------------

% For 2d only.
function cn = getCellNodes(G)
   % Construct n x 2 table of cell edges with edges oriented the same
   % direction around the cell boundary.
   cellNo    = rldecode(1:G.cells.num, diff(G.cells.facePos), 2).';
   edges     = reshape(G.faces.nodes, 2, [])';
   cellEdges = edges(G.cells.faces(:,1),:);
   ind       = G.faces.neighbors(G.cells.faces(:,1), 1) ~= cellNo;
   cellEdges(ind, :) = cellEdges(ind, [2,1]);

   % Sort edges in each cell:
   for c = 1 : G.cells.num,
      ind = G.cells.facePos(c) : G.cells.facePos(c + 1) - 1;
      cellEdges(ind, :) = sortEdges(cellEdges(ind,:));
   end
   cn = reshape(cellEdges(:,1), 1, [])';

%--------------------------------------------------------------------------

% For 2d only.
function edges = sortEdges(edges)
   % Assume edges vectors are oriented in the same direction around cell.
   % Sort edges such that they are back-to-back.
   % Then cellNodes are edges(:,1).

   for i = 1 : size(edges, 1) - 1,
      for j = i + 1 : size(edges,1),
         if edges(i,2) == edges(j,1), break; end
      end

      % Swap edges i+1 and j
      tmp = edges(i+1,:);
      edges(i+1, :) = edges(j,:);
      edges(j,   :) = tmp;
   end
