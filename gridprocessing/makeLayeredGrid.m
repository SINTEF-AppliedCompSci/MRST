function G = makeLayeredGrid(G, layerSpec)
%Extrude 2D Grid to layered 3D grid with specified layering structure
%
% SYNOPSIS:
%   G = makeLayeredGrid(G, layerSpec)
%
% PARAMETERS:
%   G         - Valid 2D areal grid.
%
%   layerSpec - Layering structure.  Interpreted as the number of layers in
%               an extruded grid (uniform thickness of 1 meter) if positive
%               scalar, otherwise as vector of layer thicknesses if array.
%               Elements should be positive numbers in the latter case.
%
% RETURNS:
%   G - Valid 3D grid as described in grid_structure.
%
% EXAMPLE:
%   % 1) Create three layers of uniform thickness (1 meter).
%   Gu = makeLayeredGrid(cartGrid([2, 2]), 3);
%   figure, plotGrid(Gu), view(10, 45)
%
%   % 2) Create five layers of increasing thickness.
%   Gi = makeLayeredGrid(cartGrid([2, 2]), convertFrom(1:5, ft));
%   figure, plotGrid(Gi), view(10, 45)
%
%   % 3) Add extra layers to a 'topSurfaceGrid' from the co2lab module.
%   %    Probably not too useful in practice.
%   G = tensorGrid(0 : 10, 0 : 10, [0, 1], ...
%                  'depthz', repmat(linspace(0, 1, 10 + 1), [1, 11]));
%   [Gt, G] = topSurfaceGrid(G);
%   Gt_extra = Gt;
%   Gt_extra.nodes.coords = [ Gt_extra.nodes.coords, Gt_extra.nodes.z ];
%   Gt_extra.nodes = rmfield(Gt_extra.nodes, 'z');
%
%   thickness = repmat(convertFrom(1:5, ft), [1, 3]);
%   thickness(6:10) = thickness(10:-1:6);  % Because 'why not?'
%   Gt_extra = makeLayeredGrid(Gt_extra, thickness);
%   figure, plotGrid(Gt_extra), view(10, 45)
%
%   % 4) Create a single layer of non-unit thickness.  This requires extra
%   %    manual steps due to a quirk of the calling interface of function
%   %    'makeLayeredGrid'.
%   G1 = makeLayeredGrid(cartGrid([2, 2]), 1);
%   k  = G1.nodes.coords(:, 3) > 0;
%   G1.nodes.coords(k, 3) = 1234*milli*meter;
%   figure, plotGrid(G1), view(10, 45)
%
% NOTE:
%   The special treatment of a *scalar* `layerSpec` parameter, to preserve
%   backwards compatibility with the original semantics of this function,
%   means that it is not possible to specify a single layer of non-unit
%   thickness.  If you need a single layer of non-unit thickness then you
%   need to manually update the third column of `G.nodes.coords`.
%
% SEE ALSO:
%   `grid_structure`

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


% Build layered grid from 2D mesh.

assert(G.griddim == 2, ...
   'Areal grid must be 2D.  Field ''griddim'' says otherwise...');
G.type = [G.type, { mfilename }];

% Determine layering structure of the final extruded grid.
%-----------------------------
if numel(layerSpec) == 1
   % User passed number of layers.
   nlayers   = fix(layerSpec);
   thickness = ones([nlayers, 1]);
else
   % User passed vector of layer thicknesses.
   nlayers   = numel(layerSpec);
   thickness = reshape(layerSpec, [], 1);
end

% Faces with horizontal normal
%-----------------------------
edges = reshape(double(G.faces.nodes(:,1)), 2, []) .';
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
vFaces = repmat(double(cn), [nlayers + 1, 1]) + ...
         kron((0:nlayers)', ones(size(cn,1), 1))*dz;
vNumNodes = repmat(diff(G.cells.facePos), [nlayers+1, 1]);

% Neighbor relations in vertical direction.
vNeighbors = zeros((nlayers+1)*G.cells.num,2);
cells      = 1:nlayers*G.cells.num;
vNeighbors(cells,       2) = cells;
vNeighbors(cells+G.cells.num,1) = cells;

% Build grid structure
%---------------------------
G.nodes.coords    = extrude_coordinates(G.nodes.coords, thickness);

G.nodes.num       = (nlayers + 1) * G.nodes.num;
G.cells.num       = nlayers * G.cells.num;

G.faces.neighbors = [vNeighbors; hNeighbors];
G.faces.num       = size(G.faces.neighbors,1);
G.faces.nodes     = [vFaces(:); reshape(hFaces .', [], 1)];

numFaces = repmat(diff(G.cells.facePos) + 2, [nlayers, 1]);
numNodes = [vNumNodes; hNumNodes];
pos      = @(n) cumsum([1; double(reshape(n, [], 1))]);

G.faces.nodePos = pos(numNodes);
G.cells.facePos = pos(numFaces);

% Cell topology from G.faces.neighbors
N   = G.faces.neighbors;
tmp = sortrows([N(:,1), (1:G.faces.num).',  ones([G.faces.num, 1]); ...
                N(:,2), (1:G.faces.num).', -ones([G.faces.num, 1])]);

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

function coords = extrude_coordinates(coords, thickness)
   nlayers  = numel(thickness);
   numnodes = size(coords, 1);

   zcoord = reshape(rldecode(cumsum([0 ; thickness]), numnodes), ...
                    [], nlayers + 1);

   if size(coords, 2) == 3
      % Possibly a topSurfaceGrid from the co2lab module.  Grid possibly
      % has non-constant surface depth.  Adjust 'zcoord' accordingly.
      zcoord = bsxfun(@plus, coords(:,3), zcoord);
   end

   coords = ...
      [repmat(coords(:, [1, 2]), [nlayers + 1, 1]), reshape(zcoord, [], 1)];

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
