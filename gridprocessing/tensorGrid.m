function G = tensorGrid(x, varargin)
%Construct Cartesian grid with variable physical cell sizes.
%
% SYNOPSIS:
%   G = tensorGrid(x)
%   G = tensorGrid(x, y)
%   G = tensorGrid(x, y, 'depthz', dz)
%   G = tensorGrid(x, y, z)
%   G = tensorGrid(x, y, z, 'depthz', dz)
%
% PARAMETERS:
%   x,y,z    - Vectors giving cell vertices, in units of meters, of individual
%              coordinate directions.  Specifically, the grid cell at
%              logical location (I,J,K) will have a physical dimension of
%              [x(I+1)-x(I), y(J+1)-y(J), z(K+1)-z(K)] (meters).
%
%   dz       - Depth, in units of meters, at which upper reservoir nodes
%              are encountered.  Assumed to be a
%              NUMEL(x)-by-NUMEL(y) array of nodal depths.
%
%              OPTIONAL.
%              Default value: depthz = ZEROS([numel(x), numel(y)])
%                             (i.e., top of reservoir at zero depth).
%
%   cellnodes- OPTIONAL.
%              Default value FALSE.  If TRUE, the corner points of each
%              cell is added as field G.cellNodes.  The field has one row
%              per cell, the sequence of nodes on each is (imin,
%              jmin,kmin), (imax,jmin,kmin), (imin,jmax,kmin), ...
%
% RETURNS:
%   G - Grid structure with a subset of the fields `grid_structure`.
%       Specifically, the geometry fields are missing:
%         - G.cells.volumes
%         - G.cells.centroids
%
%         - G.faces.areas
%         - G.faces.normals
%         - G.faces.centroids
%
%       These fields may be computed using the function `computeGeometry`.
%
%       There is, however, an additional field not described in
%       `grid_structure:
%
%           `cartDims` is a length 1, 2 or 3 vector giving number of cells
%           in each coordinate direction.  In other words 
%
%                      `all(G.cartDims == celldim)`.
%
%       `G.cells.faces(:,2)` contains integers 1-6 corresponding to
%       directions W, E, S, N, T, B respectively.
%
%
% SEE ALSO:
%   `grid_structure`, `computeGeometry`

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

   % Find first char option
   firstOpt = find(cellfun(@ischar, varargin), 1, 'first');
   dim = nargin();
   if ~isempty(firstOpt)
       % If option N is char, then the dimension is N-1 since the first N-1
       % entries correspond to coordinates. Note that there is an
       % off-by-one offset here since we always assume at least one
       % coordinate is given.
       dim = firstOpt;
   end

   switch dim
       case 1
           G = tensorGrid1D(x, varargin{:});
       case 2
           G = tensorGrid2D(x, varargin{:});
       case 3
           G = tensorGrid3D(x, varargin{:});
       otherwise
           error('Invalid grid dimension "%d"', dim);
   end

   % Record grid constructor in grid.
   G.type    = { mfilename };
   G.griddim = numel(G.cartDims);
end

function G = tensorGrid3D(x, y, z, varargin)
   mrstNargInCheck(3, 5, nargin);

   %-----------------------------------------------------------------------
   % Check input data -----------------------------------------------------
   x = x(:); dx=diff(x);
   if ~all(dx > 0)
      warning('tensorGrid:xData', 'Nonmonotone x-data, truncating..');
      while ~all(dx>0)
         i = find(dx>0); x = x([1; i+1]); dx = diff(x);
      end
   end
   y=y(:); dy = diff(y);
   if ~all(dy > 0)
      warning('tensorGrid:yData', 'Nonmonotone y-data, truncating..');
      while ~all(dy>0)
         i = find(dy>0); y = y([1; i+1]); dy = diff(y);
      end
   end
   z=z(:); dz = diff(z);
   if ~all(dz > 0)
      warning('tensorGrid:yData', 'Nonmonotone z-data, truncating..');
      while ~all(dz>0)
         i = find(dz>0); z = z([1; i+1]); dz = diff(z);
      end
   end

   sx = numel(x)-1; sy = numel(y)-1; sz = numel(z)-1;
   celldim = [sx, sy, sz];

   opt = struct('depthz', zeros(sx+1, sy+1), 'cellnodes', false);
   opt = merge_options(opt, varargin{:});
   if numel(opt.depthz) ~= (sx+1) * (sy+1)
      error(msgid('DepthZ:WrongSize'), ...
         'Input argument ''depthz'' is wrongly sized.')
   end

   numC = sx * sy * sz;             % Number of cells.
   numN = (sx+1) * (sy+1) * (sz+1); % Number of nodes.

   numFX = (sx+1) * sy * sz;        % Number of faces parallel to yz-plane.
   numFY = sx * (sy+1) * sz;        % Number of faces parallel to xz-plane.
   numFZ = sx * sy * (sz+1);        % Number of faces parallel to xy-plane.
   numF  = numFX + numFY + numFZ;

   %--------------------------------------------------------------------------
   % Nodes/Coordinates -------------------------------------------------------

   [xCoord, yCoord, zCoord] = ndgrid(x, y, z);

   zCoord = bsxfun(@plus, zCoord, reshape(opt.depthz, [sx, sy] + 1));

   coords = [xCoord(:), yCoord(:), zCoord(:)];

   %--------------------------------------------------------------------------
   % Generate face-edges ----------------------------------------------------

   % Node index matrix
   N = reshape(1 : numN, [sx+1, sy+1, sz+1]);

   %--------------------------------------------------------------------------
   % x-faces -----------------------------------------------------------------
   %
   NF1 = reshape(N(1:sx+1, 1:sy  , 1:sz  ), 1, []);
   NF2 = reshape(N(1:sx+1, 2:sy+1, 1:sz  ), 1, []);
   NF3 = reshape(N(1:sx+1, 2:sy+1, 2:sz+1), 1, []);
   NF4 = reshape(N(1:sx+1, 1:sy  , 2:sz+1), 1, []);

   faceNodesX = reshape([NF1; NF2; NF3; NF4], [], 1);

   %--------------------------------------------------------------------------
   % y-faces -----------------------------------------------------------------
   %
   NF1 = reshape(N(1:sx  , 1:sy+1, 1:sz  ), 1, []);
   NF2 = reshape(N(1:sx  , 1:sy+1, 2:sz+1), 1, []);
   NF3 = reshape(N(2:sx+1, 1:sy+1, 2:sz+1), 1, []);
   NF4 = reshape(N(2:sx+1, 1:sy+1, 1:sz  ), 1, []);

   faceNodesY = reshape([NF1; NF2; NF3; NF4], [], 1);

   %--------------------------------------------------------------------------
   % z-faces -----------------------------------------------------------------
   %
   NF1 = reshape(N(1:sx  , 1:sy,   1:sz+1), 1, []);
   NF2 = reshape(N(2:sx+1, 1:sy,   1:sz+1), 1, []);
   NF3 = reshape(N(2:sx+1, 2:sy+1, 1:sz+1), 1, []);
   NF4 = reshape(N(1:sx  , 2:sy+1, 1:sz+1), 1, []);

   faceNodesZ = reshape([NF1; NF2; NF3; NF4], [], 1);

   %--------------------------------------------------------------------------
   % Assemble grid_structure faceNodes structure -----------------------------
   %
   faceNodes = [faceNodesX; ...
      faceNodesY; ...
      faceNodesZ];

   clear -regexp ^faceNodes. ^N ^NF.
   %--------------------------------------------------------------------------
   % Generate cell-faces ----------------------------------------------------

   foffset = 0;
   % Face index matrices
   FX = reshape(foffset + (1:numFX), sx+1, sy  , sz  ); foffset = foffset + numFX;
   FY = reshape(foffset + (1:numFY), sx  , sy+1, sz  ); foffset = foffset + numFY;
   FZ = reshape(foffset + (1:numFZ), sx  , sy  , sz+1);

   F1 = reshape(FX(1:sx  , :, :), 1, []); %W == 1
   F2 = reshape(FX(2:sx+1, :, :), 1, []); %E == 2

   F3 = reshape(FY(:, 1:sy  , :), 1, []); %S == 3
   F4 = reshape(FY(:, 2:sy+1, :), 1, []); %N == 4

   F5 = reshape(FZ(:, :, 1:sz  ), 1, []); %T == 5
   F6 = reshape(FZ(:, :, 2:sz+1), 1, []); %B == 6

   cellFaces = [reshape([F1; F2; F3; F4; F5; F6], [], 1), ...
      kron(ones([numC, 1]), [ 1, 2, 3, 4, 5, 6]')];

   clear -regexp ^F.

   %--------------------------------------------------------------------------
   % Generate neighbors -----------------------------------------------------

   % Cell index matrix
   C = zeros([sx+2, sy+2, sz+2]);
   C(2:sx+1, 2:sy+1, 2:sz+1) = reshape(1:numC, [sx, sy, sz]);

   NX1 = reshape(C( 1:sx+1, 2:sy+1, 2:sz+1), [], 1);
   NX2 = reshape(C( 2:sx+2, 2:sy+1, 2:sz+1), [], 1);

   NY1 = reshape(C( 2:sx+1, 1:sy+1, 2:sz+1), [], 1);
   NY2 = reshape(C( 2:sx+1, 2:sy+2, 2:sz+1), [], 1);

   NZ1 = reshape(C( 2:sx+1, 2:sy+1, 1:sz+1), [], 1);
   NZ2 = reshape(C( 2:sx+1, 2:sy+1, 2:sz+2), [], 1);

   neighbors = [ [NX1, NX2]; ...
      [NY1, NY2]; ...
      [NZ1, NZ2] ];

   clear -regexp ^N..

   %-----------------------------------------------------------------------
   % Generate cell nodes -------------------------------------------------
   if opt.cellnodes
      % Index to first node in each cell
      k  = firstnodeindex([sx+1, sy+1, sz+1], 1:sx, 1:sy, 1:sz);
      di = 1;
      dj = sx+1;
      dk = (sx+1)*(sy+1);
      cNodes = [k, k+di, k+dj, k+di+dj, k+dk, k+di+dk, k+dj+dk, k+di+dj+dk];
   end
   %--------------------------------------------------------------------------
   % Assemble structure -----------------------------------------------------

   G.cells = struct('num',      numC,                   ...
                    'facePos',  (1:6:(numC+1)*6)', ...
                    'indexMap', (1 : numC)');

   G.faces = struct('num',       numF,                   ...
                    'nodePos',   (1:4:(numF+1)*4)', ...
                    'neighbors', (neighbors),       ...
                    'tag',       zeros(numF, 1));

   G.nodes = struct('num', numN, 'coords', coords);

   G.cells.faces = cellFaces;
   G.faces.nodes = faceNodes;
   if opt.cellnodes
      G.cellNodes = cNodes;
   end
   G.cartDims = celldim;
end


function G = tensorGrid2D(x, y, varargin)
   mrstNargInCheck(2, 4, nargin);
   opt = struct('cellnodes', false);
   opt = merge_options(opt, varargin{:});

   %-----------------------------------------------------------------------
   % Check input data -----------------------------------------------------
   x = x(:); dx=diff(x);
   if ~all(dx > 0)
      warning('tensorGrid:xData', 'Nonmonotone x-data, truncating..');
      while ~all(dx>0)
         i = find(dx>0); x = x([1; i+1]); dx = diff(x);
      end
   end
   y=y(:); dy = diff(y);
   if ~all(dy > 0)
      warning('tensorGrid:yData', 'Nonmonotone y-data, truncating..');
      while ~all(dy>0)
         i = find(dy>0); y = y([1; i+1]); dy = diff(y);
      end
   end

   sx = numel(x)-1; sy = numel(y)-1;
   celldim = [sx, sy];

   numC = sx * sy;             % Number of cells.
   numN = (sx+1) * (sy+1);     % Number of nodes.s

   numFX = (sx+1) * sy;        % Number of faces parallel to y-axis.
   numFY = sx * (sy+1);        % Number of faces parallel to x-axis.
   numF  = numFX + numFY;

   %-----------------------------------------------------------------------
   % Nodes/Coordinates ----------------------------------------------------

   [xCoord, yCoord] = ndgrid(x, y);
   coords = [xCoord(:), yCoord(:)];

   %-----------------------------------------------------------------------
   % Generate face-edges --------------------------------------------------

   % Node index matrix
   N = reshape(1 : numN, [sx+1, sy+1]);

   %-----------------------------------------------------------------------
   % x-faces --------------------------------------------------------------
   %
   NF1 = reshape(N(1:sx+1, 1:sy  ), 1, []);
   NF2 = reshape(N(1:sx+1, 2:sy+1), 1, []);

   faceNodesX = reshape([NF1; NF2], [], 1);

   %-----------------------------------------------------------------------
   % y-faces --------------------------------------------------------------
   %
   NF1 = reshape(N(1:sx  , 1:sy+1), 1, []);
   NF2 = reshape(N(2:sx+1, 1:sy+1), 1, []);

   faceNodesY = reshape([NF2; NF1], [], 1);
   % Note:
   % Nodes need to be reversed to obtain normals
   % pointing in positive i-direction in computeGeometry.

   %-----------------------------------------------------------------------
   % Assemble grid_structure faceNodes structure --------------------------
   %
   faceNodes = [faceNodesX; ...
                faceNodesY];

   clear -regexp ^faceNodes. ^N ^NF.
   %-----------------------------------------------------------------------
   % Generate cell-faces --------------------------------------------------

   foffset = 0;
   % Face index matrices
   FX = reshape(foffset + (1:numFX), sx+1, sy  ); foffset = foffset + numFX;
   FY = reshape(foffset + (1:numFY), sx  , sy+1);

   F1 = reshape(FX(1:sx  , :), 1, []); %W == 1
   F2 = reshape(FX(2:sx+1, :), 1, []); %E == 2

   F3 = reshape(FY(:, 1:sy  ), 1, []); %S == 3
   F4 = reshape(FY(:, 2:sy+1), 1, []); %N == 4

   cellFaces = [reshape([F1; F3; F2; F4], [], 1), ...
                kron(ones([numC, 1]), [ 1, 3, 2, 4]')];

   clear -regexp ^F.

   %-----------------------------------------------------------------------
   % Generate neighbors ---------------------------------------------------

   % Cell index matrix
   C = zeros([sx+2, sy+2]);
   C(2:sx+1, 2:sy+1) = reshape(1:numC, [sx, sy]);

   NX1 = reshape(C( 1:sx+1, 2:sy+1), [], 1);
   NX2 = reshape(C( 2:sx+2, 2:sy+1), [], 1);

   NY1 = reshape(C( 2:sx+1, 1:sy+1), [], 1);
   NY2 = reshape(C( 2:sx+1, 2:sy+2), [], 1);

   neighbors = [[NX1, NX2]; [NY1, NY2]];

   clear -regexp ^N..

   %-----------------------------------------------------------------------
   % Generate cell nodes --------------------------------------------------
   if opt.cellnodes
      % Index to first node in each cell
      k  = firstnodeindex([sx+1, sy+1], 1:sx, 1:sy);
      di = 1;
      dj = sx+1;
      cNodes = [k, k+di, k+dj, k+di+dj];
   end
   %-----------------------------------------------------------------------
   % Assemble structure ---------------------------------------------------

   G.cells = struct('num',      numC,                   ...
                    'facePos',  (1:4:(numC+1)*4)', ...
                    'indexMap', (1 : numC)');

   G.faces = struct('num',       numF,                   ...
                    'nodePos',   (1:2:(numF+1)*2)', ...
                    'neighbors', neighbors,       ...
                    'tag',       zeros(numF, 1));

   G.nodes = struct('num', numN, 'coords', coords);

   G.cells.faces = cellFaces;
   G.faces.nodes = faceNodes;
   if opt.cellnodes
      G.cellNodes = cNodes;
   end
   G.cartDims = celldim;
end

function G = tensorGrid1D(x, varargin)
   mrstNargInCheck(1, 3, nargin);
   opt = struct('cellnodes', false);
   opt = merge_options(opt, varargin{:});

   %-----------------------------------------------------------------------
   % Check input data -----------------------------------------------------
   x = x(:); dx=diff(x);
   if ~all(dx > 0)
      warning('tensorGrid:xData', 'Nonmonotone x-data, truncating..');
      while ~all(dx>0)
         i = find(dx>0); x = x([1; i+1]); dx = diff(x);
      end
   end
   sx = numel(x)-1;
   celldim = sx;

   numC = sx;         % Number of cells.
   numN = (sx+1);     % Number of nodes.
   numF = numN;       % Number of faces.

   %-----------------------------------------------------------------------
   % Nodes/Coordinates ----------------------------------------------------

   coords = x;
   v = (0:numC)';
   neighbors = [v, circshift(v, -1)];

   G.cells = struct('num',      numC,                   ...
                    'facePos',  (1:2:(numC+1)*2)', ...
                    'indexMap', (1 : numC)');

   G.faces = struct('num',       numF,                   ...
                    'nodePos',   (1:numF+1)', ...
                    'neighbors', neighbors,       ...
                    'tag',       zeros(numF, 1));

   G.nodes = struct('num', numN, 'coords', coords);

   cellFaces = [1; reshape(repmat(2:numC, 2, 1), [], 1); numC+1];
   cellFaceTag = repmat([1; 2], numC, 1);
   % Faces are equal to nodes for this type of grid
   G.cells.faces = [cellFaces, cellFaceTag];
   G.faces.nodes = (1:numF)';
   if opt.cellnodes
      G.cellNodes = (1:numC)';
   end
   G.cartDims = celldim;
end


function k = firstnodeindex(sz, varargin)
   assert(numel(sz) == numel(varargin));
   assert(all(cellfun(@min, varargin) > 0));
   assert(all(cellfun(@max, varargin)<= sz));

   r       = cell(size(varargin));
   [r{:}]  = ndgrid(varargin{:});
   k       = reshape(sub2ind(sz, r{:}), [], 1);
end
