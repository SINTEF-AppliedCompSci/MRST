function G = processGRDECL(grdecl, varargin)
%Compute grid topology and geometry from pillar grid description.
%
% SYNOPSIS:
%   G = processGRDECL(grdecl)
%   G = processGRDECL(grdecl, 'pn1', pv1, ...)
%
% DESCRIPTION:
%
%  This code is designed to compute connectivity of fairly general
%  cornerpoint grids described in Eclipse `SPECGRID/COORDS/ZCORN` format.
%  In short, the algorithm consist of
%
%  a) Compute 8 node coordinates of each grid block `buildCornerPtNodes`.
%
%  b) Find unique points by comparing point coordinates and make matrix `P`
%     of point numbers for each grid block.
%
%  c) Add auxillary top and bottom layer to grid to ease processing of
%     outer boundary in the presence of faults.
%
%  d) Compute connectivity and corresponding face topology in i- j-
%     directions by considering pillar pairs.
%
%     Faulted pillar pairs, those where there is at least one non-matching
%     cell pair (non-neighboring connection), are processed separately by
%     `findFaults`, the remaining pillar pairs with with matching cells are
%     processed by `findFaces`.
%
%     Connectivity computation does not (ever) use the coordinates of
%     inactive cells as they are undefined by the grid format.  The
%     boundary between an active and an inactive region is considered as
%     outer boundary.
%
%     Collapsed faces and cells resulting form pinched layers are removed.
%
%  e) Compute connectivity in the k-direction by `findVerticalFaces`. The
%     connectivity over pinched layers is restored.
%
%  f) Build grid struct `buildGrid`, remove auxillary layers and pinched
%     cells (`removeCells`), change cell and point numbering accordingly.
%
%  g) Check if grid is connected and reasonable.
%
%  For each pair of pillars (1) and (2), compute geometric neighbor cells
%  and geometry of corresponding intersection of cell faces::
%
%
%                (1)                        (2)
%                 |                          |
%        b(1,1)   *                          |
%                 |  *                       |
%        a(1,1)   o-----*--------------------o  a(1,2)
%                 |        *                 |
%                 |           *              |
%        b(2,1)   *              *           |
%                 |   *             *        |
%        a(2,1)   o-------*------------*-----o  a(2,2)
%                 |           *           *  |
%                 |               *          *  b(1,2)
%                 |                   *      |
%                 |                       *  |
%                 |                          *  b(2,2)
%                 |                          |
%        b(3,1)   * * * * * * * * * * * * * *|  b(3,2)
%                 |                          |
%                 |                          |
%        a(3,1)   o--------------------------o  a(3,2)
%                 |                          |
%                 |                          |
%                 |                          |
%        b(4,1)   * * * * * * * * * * * * * *|  b(4,2)
%                 |                          |
%                 |                          |
%
%
%
%
%  Each row in point lists a and b correspond to a line as in the figure.
%  For cornerpoint grids, internal lines are repeated.  In the code below,
%  the odd spaces correspond to actual cells while the even spaces
%  correspond to inaccuracies in the format (or void space in the grid) and
%  are marked as inactive (with cell number 0)
%
%  Once the lines that corespond to active cells are identified, the
%  findConnections function loops through the z-coordinates of a- and
%  b-points to find the geometrical neighbors.  The face geometries are
%  computed last.
%
% PARAMETERS:
%   grdecl - Raw pillar grid structure, as defined by function
%            `readGRDECL`, with fields `COORDS`, `ZCORN` and, possibly,
%            `ACTNUM`.
%
% KEYWORD ARGUMENTS:
%
%   Verbose           - Whether or not to display progress information
%                       Logical.  Default value: `Verbose = mrstVerbose()`.
%
%   Tolerance         - Minimum distinguishing vertical distance for
%                       points along a pillar.  Specifically, two points
%                       (x1,y1,z1) and (x2,y2,z2) are considered separate
%                       only if ABS(z2 - z1) > Tolerance.
%                       Non-negative scalar.
%                       Default value: Tolerance = 0.0 (distinguish all
%                       points along a pillar whose z coordinate differ
%                       even slightly).
%
%   CheckGrid         - Whether or not to perform basic consistency
%                       checks on the resulting grid.
%                       Logical.  Default value: CheckGrid = true.
%
%   SplitDisconnected - Whether or not to split disconnected grid
%                       components into separate grids/reservoirs.
%                       Logical.  Default value: SplitDisconnected=true.
%
%   RepairZCORN       - Make an effort to detect and repair artifacts
%                       that may occur in the corner-point depth
%                       specification.  Specifically, detect and repair
%                       the following, rare, conditions:
%                           - Upper corners of a cell below lower corners
%                             of that same cell
%                           - Lower corners of a cell below that cell's
%                             lower neighbour's upper corners.
%                       Logical.  Default value: RepairZCORN = false.
%
%   PreserveCpNodes   - Whether or not to capture the vertex indicies of
%                       each cell's original corner-point nodes.  If true,
%                       an additional G.cells.num-by-8 array named
%                       'cpnodes' containing vertex indicies will be stored
%                       in the 'cells' substructure of the fully
%                       constructed grid.  Columns 1:4 are the vertices of
%                       the "minimum K" surface while columns 5:8 are the
%                       vertices of the "maximum K" surface.
%
%                       The columns are stored in natural ordering of the
%                       vertices, as demonstrated in the figure below.
%
%                            +-----------> I
%                           /|
%                          / |     1 --------- 2
%                         /  |    /|          /|   Column numbers mapped
%                      J v   |   / |         / |   to vertices in ECLIPSE's
%                            |  3 --------- 4  |   default right-handed
%                          K v  |  |        |  |   coordinate system
%                               |  5 -------|- 6   (origin is top, left,
%                               | /         | /    back vertex with Z-axis
%                               |/          |/     pointing down.)
%                               7 --------- 8
%
% RETURNS:
%   G - Valid grid definition containing connectivity, cell geometry, face
%       geometry and unique nodes.
%
%       If the pillar grid structure contains a pinch-out definition field
%       `PINCH`, then the grid structure will contain a separate field
%       `nnc` in the top-level structure if any explicit non-neighbouring
%       connections are deemed to exist according to the tolerances
%       specified in `PINCH`.  See function `processPINCH` for a more
%       detailed description of the `nnc` field.
%
% EXAMPLE:
%   G = processGRDECL(readGRDECL('small.grdecl'));
%   plotGrid(G); view(10,45);
%
% SEE ALSO:
%   `readGRDECL`, `processPINCH`, `removeCells`, `checkAndRepairZCORN`.

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

% Relies on buildCornerPtNodes, rlencode, rldecode, removeCells, dispif,
% tocif
opt = struct('CheckGrid'        , true       , ...
             'RepairZCORN'      , false      , ...
             'SplitDisconnected', true       , ...
             'Tolerance'        , 0.0        , ...
             'Verbose'          , mrstVerbose, ...
             'PreserveCpNodes'  , false);
opt = merge_options(opt, varargin{:});


if opt.Tolerance > 0
   grdecl.ZCORN = opt.Tolerance * round(grdecl.ZCORN / opt.Tolerance);
end

[X, Y, Z, reverseaxis, actnum, numAuxillaryCells] = ...
   build_coordinates(grdecl, opt.RepairZCORN);

dispif(opt.Verbose,                                   ...
       'Adding %d artificial cells at top/bottom\n\n', ...
       numAuxillaryCells);

%  -----------------------------------------------------------------------
% Replace nan coordinates by inf to avoid special-purpose code below
%
% Nan may occur if top and bottom pillar coordinates coincide.  This
% should only occur in inactive cells (assert this!).  Point comparison is
% used to find unique points and ultimately to detect faults.
%
% nan==nan is always false, while inf==inf is true. Replacing nan with inf
% will result in fewer unique points and fewer faults to process and
% finally less special purpose code to handle nan coordinates.

% TODO:  Assert that NaN junk occurs only in inactive cells

X(isnan(X)) = inf;
Y(isnan(Y)) = inf;

[G, P, B] = create_initial_grid(X, Y, Z, opt.PreserveCpNodes);  clear X Y Z

% -------------------------------------------------------------------------
% Process faces with constant i-index
% Note
% In all caluculations we pass arrays of point numbers P, block numbers
% B and the point coordinates (G.nodes.coords)

tags = [1, 2]; % i.e., [West, East]
G = process_pillar_faces(G, P, B, actnum, tags, 'i', opt);

% ------------------------------------------------------------------------
% Process faces with constant j-index
% Switch i and j indices in order to reuse code for i-faces.  This will
% require some face surgery below.
P = permute(P, [2, 1, 3]);
B = permute(B, [2, 1, 3]);

tags = [3, 4]; % i.e., [South, North]
G = process_pillar_faces(G, P, B, actnum, tags, 'j', opt);


% Due to permutation of i and j indices, the orientation of face nodes is
% reversed. To fix the orientation, reverse the sequence of nodes for each
% j-face:
pos = G.faces.nodePos - 1;
i   = find(any(G.faces.cellTags == 3 | G.faces.cellTags == 4, 2));
j   = mcolon(pos(i)+1, pos(i+1));
k   = mcolon(pos(i+1), pos(i)+1, repmat(-1, [numel(i), 1]));

G.faces.nodes(j)=G.faces.nodes(k);
clear i j k pos B

% To be on the safe side, switch back.
P = permute(P, [2, 1, 3]);


% ------------------------------------------------------------------------
dispif(opt.Verbose, 'Processing regular k-faces\n');

tags = [5,6]; % i.e., [Top, Bottom]
t0   = ticif(opt.Verbose);
G    = findVerticalFaces(G, P, actnum, tags, opt);
tocif(opt.Verbose, t0);
dispif(opt.Verbose, '\n');
clear P

% Deactivate auxillary layers as they are no longer needed
actnum(:,:,1)   = false;
actnum(:,:,end) = false;


% ------------------------------------------------------------------------
dispif(opt.Verbose, 'Building grid structure\n');
G   = buildCellFaces(G, G.faces.cellTags);
clear tags

% pinched cells == cells with no faces
dispif(opt.Verbose, 'removing %d artificial cells at top/bottom\n', ...
       numAuxillaryCells);
dispif(opt.Verbose, 'removing %d inactive and pinched cells\n', ...
       sum(~actnum(:)|diff(G.cells.facePos)==0) - numAuxillaryCells);


% Remove auxillary and inactive cells
G = removeCells(G, find(~actnum(:)|diff(G.cells.facePos)==0));

% When the auxillary layers are removed, we must change the index map
% and Cartesian dimensions
G.cells.indexMap = G.cells.indexMap-prod(G.cartDims(1:2));
G.cartDims(3)    = G.cartDims(3)-2;
G.faces = rmfield(G.faces, 'cellTags');

% ------------------------------------------------------------------------
% A sane grid cell has at least four faces
if opt.CheckGrid
   assert(all(diff(G.cells.facePos)>3));
   assert(all(diff(G.faces.nodePos)>2));
   assert(all(all(~isinf(G.nodes.coords))));
end


% If coordinate axes have been reversed, restore original directions

G.nodes.coords(:, reverseaxis) = -G.nodes.coords(:, reverseaxis);
if mod(sum(reverseaxis), 2)
   G = reverseFaceNodes(G, 1:G.faces.num);
end

if isfield(grdecl, 'PINCH')
   nnc = processPINCH(grdecl, G);

   if ~ isempty(nnc)
      G.nnc = nnc;
   end
end

% Check if grid is connected
if opt.SplitDisconnected
   G = splitDisconnectedGrid(G, 'verbose', opt.Verbose);
end

[G.type]    = deal({ mfilename });
[G.griddim] = deal(3);
end

%==========================================================================
%
%     HELPERS FOLLOW
%
%==========================================================================

function [X, Y, Z, reverseaxis, actnum, numAuxillaryCells] = ...
      build_coordinates(grdecl, repair)

   if repair
      args = {};

      if isfield(grdecl, 'ACTNUM')
         args = [ args, { 'Active', grdecl.ACTNUM } ];
      end

      grdecl.ZCORN = ...
         checkAndRepairZCORN(grdecl.ZCORN, grdecl.cartDims, args{:});
   end

   [X, Y, Z] = buildCornerPtNodes(grdecl);

   % Assume right-handed coordinate system
   reverseaxis = false([1, 3]);

   % Check that ZCORN increases along each pillar.
   % Expand ACTNUM by 2 in each grid dimension
   if isfield(grdecl, 'ACTNUM')
      a = reshape(grdecl.ACTNUM ~= 0, grdecl.cartDims);
   else
      % No ACTNUM.  Assume all cells active.
      a = true(grdecl.cartDims);
   end

   a  = rldecode(a, 2, 1);
   a  = rldecode(a, 2, 2);
   a  = rldecode(a, 2, 3);
   z  = Z; z(~a) = nan;
   dz = diff(z, 1, 3);
   if ~any(dz(:) > 0)
      Z     = -Z;
      z     = Z;
      z(~a) = nan;
      dz    = diff(z, 1, 3);
      reverseaxis(3) = true;
   end

   if any(dz(:) < 0)
      if ~ repair
         warning('processGRDECL:ZCORN', ...
                ['Non-monotonous corner-point depths detected.\n', ...
                 'Possibly wrong results/geometry in grid.\n',     ...
                 'Re-run with option ''RepairZCORN'' enabled to fix.']);
      else
         error(['Option ''RepairZCORN'' failed to correct non-', ...
                'monotonous corner-point depths. Programming error?']);
      end
   end

   % Check if we have a right handed system
   if is_lefthanded(X, Y, Z)
      Y              = -Y;
      reverseaxis(2) = true;
   end

   %-------------------------------------------------------------------
   % Add auxillary top and bottom layer to ensure correct processing of
   % outer boundary at faults.
   minz = min(Z(:));
   maxz = max(Z(:));

   e      = zeros(size(Z(:,:,1)));
   top    = cat(3, minz - 2 + e, minz - 1 + e);
   bottom = cat(3, maxz + 1 + e, maxz + 2 + e);

   Z      = cat(3, top, Z, bottom);
   X      = X(:, :, [1, 1, 1:end, end, end]);
   Y      = Y(:, :, [1, 1, 1:end, end, end]);

   % Mark auxillary layers as active
   if ~isfield(grdecl, 'ACTNUM')
      grdecl.ACTNUM = ones([prod(grdecl.cartDims), 1]);
   end
   actnum   = reshape(double(grdecl.ACTNUM), grdecl.cartDims);
   extLayer = true(size(actnum(:,:,1)));
   actnum   = cat(3, extLayer, actnum, extLayer);

   numAuxillaryCells = 2 * numel(extLayer);
end

%--------------------------------------------------------------------------

function bool = is_lefthanded(X, Y, Z)
   V1 = [X(end,1,1) - X(1,1,1), Y(end,1,1) - Y(1,1,1), 0                  ];
   V2 = [X(1,end,1) - X(1,1,1), Y(1,end,1) - Y(1,1,1), 0                  ];
   V3 = [0,                     0,                     Z(1,1,end)-Z(1,1,1)];

   if (norm(V1) == 0) || (norm(V2) == 0) || (norm(V3) == 0)
       warning('processGRDECL:HandTestFailure', ...
              ['Left-handed test is ill-posed. No way of verifying ', ...
               'right-handed coordinate system']);
   end

   % If triple product is negative we have left-handed coordinates.
   bool = dot(cross(V1, V2), V3) < 0;
end

%--------------------------------------------------------------------------

function [G, P, B] = create_initial_grid(X, Y, Z, preserve_cpnodes)
   % Find unique points in input.  Unique sorts with most significant
   % number in first column. We must use column order Z, Y, X to obtain
   % natural ordering of points.

   G.type = { 'INVALID' };
   G.griddim = 3;

   [G.nodes.coords, to, from] = unique([Z(:), Y(:), X(:)], 'rows');    %#ok
   G.nodes.coords             = fliplr(G.nodes.coords);
   G.nodes.num                = size(G.nodes.coords, 1);
   G.cartDims                 = size(X) / 2;
   G.cells.num                = prod(G.cartDims);
   G.cells.indexMap           = (1 : prod(G.cartDims)) .';

   if preserve_cpnodes
      G.cells.cpnodes = from(create_initial_cpnodes(size(X)));
   end

   % Empty structure used to hold grid
   G.faces.nodes      = [];
   G.faces.nodePos    = [];
   G.faces.neighbors  = zeros(0, 2);
   G.faces.tag        = [];
   G.faces.cellTags   = [];

   P = reshape(from, 2 * G.cartDims);
   B = reshape(1 : prod(G.cartDims), G.cartDims);
end

%--------------------------------------------------------------------------

function G = process_pillar_faces(G, P, B, actnum, tags, name, opt)
   dispif(opt.Verbose, 'Processing regular %s-faces\n', name);

   t0   = ticif(opt.Verbose);
   G    = findFaces(G, P, B, actnum, tags, opt);
   tocif(opt.Verbose, t0);

   dispif(opt.Verbose, '\nProcessing %s-faces on faults\n', name);

   t0 = ticif(opt.Verbose);
   G  = findFaults(G, P, B, actnum, tags, opt);
   tocif(opt.Verbose, t0);

   dispif(opt.Verbose, '\n');
end

%--------------------------------------------------------------------------

function G = reverseFaceNodes(G, f)
   ix1 = mcolon(G.faces.nodePos(f),     G.faces.nodePos(f+1)-1);
   ix2 = mcolon(G.faces.nodePos(f+1)-1, G.faces.nodePos(f)    , -1);
   G.faces.nodes(ix1) = G.faces.nodes(ix2);
end

%--------------------------------------------------------------------------

function cpnodes = create_initial_cpnodes(cp_cartdims)
   % cp_cartdims == 2 .* G.cartDims
   [i, j, k] = ndgrid(1 : 2);
   ijk = [ i(:), j(:), k(:) ];

   ix = reshape(1 : prod(cp_cartdims), cp_cartdims);
   cp = @(nod) reshape(ix(ijk(nod, 1) : 2 : end, ...
                          ijk(nod, 2) : 2 : end, ...
                          ijk(nod, 3) : 2 : end), [], 1);

   cpnodes = zeros([prod(cp_cartdims ./ 2), size(ijk, 1)]);
   for nod = 1 : size(ijk, 1)
      cpnodes(:, nod) = cp(nod);
   end
end

%--------------------------------------------------------------------------

function k = index(i, j, k, sz)
[I,J,K] = ndgrid(i, j, k);
k       = reshape(sub2ind(sz, I, J, K), [], 1);
end

%--------------------------------------------------------------------------

function G = buildCellFaces(G, cellTags)
nf              = size(G.faces.neighbors, 1);
G.faces.num     = nf;
G.nodes.num     = size(G.nodes.coords, 1);
numFaces        = accumarray(G.faces.neighbors(:)+1, 1, [G.cells.num+1, 1]);
G.cells.facePos = cumsum([1; numFaces(2:end)]);

cellTags = fliplr(cellTags);

% Find cellFaces field by sorting
vec         = (double(G.faces.neighbors(:))-1)*nf + repmat((1:nf)', [2,1]);
hf          = find(G.faces.neighbors(:)~=0);
[vec, ind]  = sort(vec(hf));
G.cells.faces = [mod(vec-1, nf) + 1, cellTags(hf(ind))];

% The above code block is equivalent to (but sligthly faster than) the
% following, where sort has been replaced by sortrows.  Speedup is
% attributable to difference between sort and sortrows.

% G.cells.faces = sortrows([G.faces.neighbors(:,1),(1:n)',  ones(n,1),cellTags(:,1);...
%                         G.faces.neighbors(:,2),(1:n)', -ones(n,1),cellTags(:,2)]);
% ind                 = G.cells.faces(:,1)==0;
% G.cells.faces(ind,:)  = [];
% G.cells.faces         = int32(G.cells.faces(:,[2,4]));
end

%--------------------------------------------------------------------------

function G = findFaces(G, P, B, actnum, tags, opt)
%Find REGULAR faces including their nodes, cell-neighbors and direction
%
% See below for details of what is considered a regular face.
%
% SYNOPSIS:
%   G = processGRDECL(G, P, B, actnum, tags, opt)
%
% PARAMETERS:
%   G      - Grid struct that faces will be filled into.
%
%   P      - 3d-array of point numbers for each cell in the grid, .i.e.,
%            the grid cell (i,j,k) has corner points numbers
%                n1 = P(2*i-1,2*j-1,2*k-1)
%                n2 = P(2*i  ,2*j-1,2*k-1)
%                n3 = P(2*i-1,2*j  ,2*k-1)
%                n4 = P(2*i  ,2*j  ,2*k-1)
%                n5 = P(2*i-1,2*j-1,2*k  )
%                n6 = P(2*i  ,2*j-1,2*k  )
%                n6 = P(2*i-1,2*j  ,2*k  )
%                n8 = P(2*i  ,2*j  ,2*k  )
%
%   B      - Block numbers in the grid.  These are included for convenience
%            as they could easily have been computed from (i,j,k).
%
%   actnum - Array of 0/1, 0 indicate inactive cell, 1 active cell.
%
%   tags   - Face tag (numeric) that should be inserted into G.
%
%   opt    - options struct from which opt.verbose is used.
%
% RETURNS:
%   G      - Grid structure with modified G.faces field.  No other field is
%            changed. Postprocessing takes care of changes to cells and
%            nodes.
%
% NODE:
%  Face regularity is only checked by comparing node numbers in B.  A face
%  where node numbers in the two cell neighbors coincide is considered
%  regular.  For simplicity, however, this check is only performed in the
%  i-j directions the following sense.  If the vectors of point numbers
%
%    n2 = P(2*i,    2*j-1, :) and n4 = P(2*i,   2*j  , :)
%
%  of cells B(i,j,:) are identical to the corresponding point numbers in
%  the neighbor cells B(i+1,j,:)
%
%    n1 = P(2*(i+1)-1,  2*j-1, :) and n3 = P(2*(i+1)-1, 2*j  , :),
%
%  All faces are processed as regular. If only one point number differ, the
%  whole stack of faces is processed as faulted.
%
%
% SEE ALSO:
%   `findFaults`

%  Layout of point number k
%     ----> i
%  |
%  |
%  |      o   oo   oo  ...  oo   o
%
%  |      o   oo   oo  ...  oo   o
%         o   oo   oo  ...  oo   o
%             .
%              .
%              .

sz   = size(P);
szB  = size(B);

% Find face corners
di   = 1;
dj   = sz(1);
dk   = sz(1)*sz(2);
k    = index(2:2:sz(1)-1, 1:2:sz(2), 1:2:sz(3)-1, sz);

% Internal faces
f    = [P(k), P(k+dj), P(k+dj+dk), P(k+dk)];
k    = k+di;
g    = [P(k), P(k+dj), P(k+dj+dk), P(k+dk)];

% Neighbor cells
c1   = B(index(1:szB(1)-1, 1:szB(2), 1:szB(3), szB));
c2   = B(index(2:szB(1),   1:szB(2), 1:szB(3), szB));

% Keep only cell pairs (c1,c2) that match completely (no fault)
h    = repmat(all(reshape(all(f==g, 2), sz/2-[1,0,0]), 3), [1,1,sz(3)/2]);
c1   = reshape(c1(h), [], 1);
c2   = reshape(c2(h), [], 1);
f    = f(h,:);
clear g h



% Regular boundary faces
k    = index([1 sz(1)], 1:2:sz(2), 1:2:sz(3)-1, sz);
fB   = [P(k),   P(k+dj), P(k+dj+dk), P(k+dk)];
cB   = B(index([1 szB(1)], 1:szB(2), 1:szB(3), szB));
tagB = repmat (tags(:), [numel(cB)/2, 1]);


% Append boundary faces
f    = [f;    fB]; clear fB
c1   = [c1;   cB(:).*(tagB==tags(2))];
c2   = [c2;   cB(:).*(tagB==tags(1))];



% Filter out inactive and degenerate cell pairs
% Remove inactive cells
i = c1 ~= 0;  c1(i) = c1(i) .* reshape(actnum(c1(i)),[],1);
i = c2 ~= 0;  c2(i) = c2(i) .* reshape(actnum(c2(i)),[],1);  clear i

% Remove faces with no neighbors and pinched faces.
ind = ((c1~=0) | (c2~=0)) & ~(f(:,1)==f(:,4) & f(:,2)==f(:,3)) ;
f    = f(ind,:);
c1   = c1(ind);
c2   = c2(ind);

% Remove zero-area faces
ind = f(:,1)==f(:,4) & f(:,2)==f(:,3);
f   = f(~ind, :);
c1  = c1(~ind);
c2  = c2(~ind);

% Remove repeated nodes that may arise from pinch
f(f(:,1)==f(:,4), 1)=nan;
f(f(:,2)==f(:,3), 2)=nan;

F  = reshape([f, inf(size(f, 1), 1)]', [], 1);%rlencode( reshape([f, inf(size(f, 1), 1)]', [], 1), 1);
nF = diff([0; find(isinf(F(~isnan(F))))])-1;
F  = F(isfinite(F));


% Write result to pre-grid structure
n                  =  numel(c1);
dispif(opt.Verbose, ' Found %d new regular faces\n', n);
G.faces.nodes        = [G.faces.nodes; F];
G.faces.neighbors  = [G.faces.neighbors; [c1(:), c2(:)]];
G.faces.nodePos    = cumsum([1; double(diff(G.faces.nodePos)); nF]);
G.faces.tag        = [G.faces.tag;        zeros(n, 1)];
G.faces.cellTags   = [G.faces.cellTags;   repmat(tags, [numel(c1), 1])];
end

%--------------------------------------------------------------------------

function G = findVerticalFaces(G, P, actnum, tags, opt)
%
%
%
sz         = size(P);
B          = reshape(1:prod(sz/2), sz/2);

% Find faces in k-direction.  k-index running faster than i or j.
[I,J,K]    = ndgrid(1:2:sz(1), 1:2:sz(2), 1:sz(3));
k          = reshape(permute(sub2ind(sz, I, J, K), [3,1,2]), [], 1);
f          = [P(k), P(k+1), P(k+1+sz(1)), P(k+sz(1))];
clear I J K k

% faces contains unique k-faces.
[faces, n] = rlencode(f,1);

% find neighbors in k-direction
c          = zeros(2*numel(actnum)+1,1);
c(2:2:end) = reshape(permute(B.*actnum, [3,1,2]), [], 1);
p          = c([1;1+cumsum(n)]); clear c
new        = [p(1:end-1),p(2:end)];
ind        = any(new(:,1:2)~=0, 2);
c1         = new(ind,1);
c2         = new(ind,2);
f          = faces(ind, :);

n          =  numel(c1);
dispif(opt.Verbose, ' Found %d new regular faces\n', n);
G.faces.nodes        = [G.faces.nodes;        reshape(f', [], 1)];
G.faces.neighbors  = [G.faces.neighbors;  [c1(:), c2(:)]];
G.faces.nodePos    = cumsum([1; double(diff(G.faces.nodePos)); ...
                             repmat(4, [n, 1])]);
G.faces.tag        = [G.faces.tag;        zeros(n, 1)];
G.faces.cellTags   = [G.faces.cellTags;   repmat(tags, [numel(c1), 1])];
end

%--------------------------------------------------------------------------

function G = findFaults(G, P, B, actnum, tags, opt)
%Find FAULTED faces including their nodes, cell-neighbors and direction
%
% See below for details of what is considered a regular face.
%
% SYNOPSIS:
%   G = processGRDECL(G, P, B, actnum, tags, opt)
%
% PARAMETERS:
%   G      - Grid struct that faces will be filled into.
%
%   P      - 3d-array of point numbers for each cell in the grid, .i.e.,
%            the grid cell (i,j,k) has corner points numbers
%                n1 = P(2*i-1,2*j-1,2*k-1)
%                n2 = P(2*i  ,2*j-1,2*k-1)
%                n3 = P(2*i-1,2*j  ,2*k-1)
%                n4 = P(2*i  ,2*j  ,2*k-1)
%                n5 = P(2*i-1,2*j-1,2*k  )
%                n6 = P(2*i  ,2*j-1,2*k  )
%                n6 = P(2*i-1,2*j  ,2*k  )
%                n8 = P(2*i  ,2*j  ,2*k  )
%
%   B      - Block numbers in the grid.  These are included for convenience
%            as they could easily have been compited from (i,j,k).
%
%   actnum - Array of 0/1, 0 indicate inactive cell, 1 active cell.
%
%   tags   - Face tag (numeric) that should be inserted into G.
%
%   opt    - options struct from which opt.verbose is used.
%
% RETURNS:
%   G      - Grid structure with modified G.faces field.  No other field is
%            changed. Postprocessing takes care of changes to cells and
%            nodes.
%
% NODE:
%  Face regularity is only checked by comparing node numbers in B.  A face
%  where node numbers in the two cell neighbors coincide is considered
%  regular.  For simplicity, however, this check is only performed in the
%  i-j directions the following sense.  If the vectors of point numbers
%
%    n2 = P(2*i,    2*j-1, :) and n4 = P(2*i,   2*j  , :)
%
%  of cells B(i,j,:) are identical to the corresponding point numbers in
%  the neighbor cells B(i+1,j,:)
%
%    n1 = P(2*(i+1)-1,  2*j-1, :) and n3 = P(2*(i+1)-1, 2*j  , :),
%
%  All faces are processed as regular. If only one point number differ, the
%  whole stack of faces is processed as faulted.
%
%
% SEE ALSO:
%   `findFaces`.

% Vectorized version of findFaults.
%
% Discover faults and process face geometry for faulted pillar pairs.
%
% G      = early grid structure.
% P      = array of point numbers
% B      = array of cell or block numbers
% actnum = array:is 1 for active cells, 0 otherwise.

% Move z-coordinate to first index
P    = permute(P, [3,1,2]);
B    = permute(B, [3,1,2]);
szP  = size(P); di = szP(1); dj = szP(1)*szP(2); dk = 1; %#ok
szB  = size(B);if numel(szB)==2, szB=[szB, 1];end

% Find index to point numbers belonging to every other pillar (i,j)  (which
% four point set correspond to a pillar?) in i- and j-direction. Point
% numbers are in consequtive order in k-
k    = index(1:szP(1), 2:2:szP(2)-1, 1:2:szP(3), szP);


% Point numbers for faces on side A of pillar pairs
a    = [P(k), P(k+dj)];

% Point numbers for faces on side B of pillar pairs (i
b    = [P(k+di), P(k+di+dj)]; clear k

% Cells associated with each face on side A and ond side B
cA   = B(index(1:szB(1), 1:szB(2)-1, 1:szB(3), szB));
cB   = B(index(1:szB(1), 2:szB(2),   1:szB(3), szB));


% if four point numbers of a match four point numbers of b, i.e., all
% a(2i-1:2i, :)==b(2j-1:2j,:), then cellsA(i) and cellsB(j) match exactly
% along a face.
%
% Below, we construct a logical vector that picks out all pillar pairs
% with ONLY MATCHING cells.
sz2  = [szP(1), szP(2)/2-1, szP(3)/2];
h    = repmat(all(reshape(all(a==b, 2), sz2), 1), [szP(1),1,1]);
dispif(opt.Verbose, ' Found %d faulted stacks\n', sum(~h(:))/szP(1));

% Keep node numbers of faces that DO NOT match exactly
a    = a(~h(:), :);
b    = b(~h(:), :);

% Keep stacks of cells with faces that DO NOT match exactly
cA   = cA(~h(1:2:end));
cB   = cB(~h(1:2:end));

% Make artificial z-increment to separate each stack completely in
% the next function call.
dz   = max(G.nodes.coords(:,3))-min(G.nodes.coords(:,3))+1;
auxz = reshape(1:prod(sz2(2:3)), sz2(2:3))*dz;
dZ   = permute(repmat(auxz, [1,1,szP(1)]), [3,1,2]);
dZ   = repmat(dZ(~h), [1, 2]);


% filter out inactive cells from c1 and c2 and faces belonging to inactive
% cells from a and b. Also, remove entries in dZ.
ind_a= actnum( kron(cA, [1;1]) ) == 1;
ind_b= actnum( kron(cB, [1;1]) ) == 1;
a    = a(ind_a, :);
b    = b(ind_b, :);
cA   = cA(actnum( cA ) == 1);
cB   = cB(actnum( cB ) == 1);
dZa  = repmat(dZ(ind_a), [1,2]);
dZb  = repmat(dZ(ind_b), [1,2]);


% Find z-coordinates + artificial z-increment of each point in a and b.
z    = G.nodes.coords(:,3);
za   = z(a)+dZa;
zb   = z(b)+dZb;
clear z

% Process connectivity across pillar pair, i.e., determine which faces
% b(2j-1:2j,:) on side B overlap with faces a(2i-1:2i, :) each face on side
% A. C returns index pair (2i,2j) for each cell-cell match and a (2i, 2j-1)
% of (2i-1,2j) pair for each cell-void match. These matches translate to
% cell connectivity and pieces of faces that are outer or internal
% boundaries.
%
C    = findConnections(za,zb);

% Construct face-to-cell map NEW; each row corresponds to a face, each of
% the two columns hold cell numbers.  Cell number zero mark outer/inner
% boundary face

ka   = zeros(size(a, 1)-1,1); ka (1:2:end) = cA(:);
kb   = zeros(size(b, 1)-1,1); kb (1:2:end) = cB(:);
new  = [ka(C(:,1)), kb(C(:,2))];
ind  = any(new~=0, 2);  % Keep rows with at least one non-zero
new  = new(any(new~=0, 2), :);

G.faces.neighbors  = [G.faces.neighbors;  new];
G.faces.cellTags   = [G.faces.cellTags;   repmat(tags, [size(new,1), 1])];

% Compute new node coordinates and geometry of the newly found faces
[G.nodes.coords, ncor, cor] = computeFaceGeometry(a, b, C(ind,:), ...
                                                  G.nodes.coords);
G.faces.nodePos  = cumsum([1; double(diff(G.faces.nodePos)); ncor]);
G.faces.nodes    = [G.faces.nodes;      cor];
G.faces.tag      = [G.faces.tag;      ones(numel(ncor), 1)];

dispif(opt.Verbose, ' Found %d new faces\n', numel(ncor));
end

%--------------------------------------------------------------------------

function [points, numNodes, Corners] = computeFaceGeometry(a, b, C, points)
% X-Y directions
%
% Compute intersection of cell faces on side a and b of a pair of pillars.
% The input argument C contains pairs of cell faces (a,b) that has a
% nonzero intersection. The point numbers [a(i, :), a(i+1,:)] define a cell
% face as follows:
%
%                1                2
%                |                |
%   A3 = a(i+1,1)o----------------o A4 = a(i+1,2)
%                |                |
%                |                |
%                |                |
%                |                |
%   A1 = a(i,1)  o----------------o A2 = a(i,2)
%                |                |
%
% The point numbers in b define cell faces in the same manner.
%
% Compute intersection given that these points share pillars, i.e.,
% a(:,1), b(:,1) are collinear (likewise for a(:,2), b(:,2)).  There are 8
% different points that may appear in the face between two cornerpoint
% cells:
%
%   - The 4 points lying on the pillars are the lowest of the top lines
%     and the highest of the bottom lines:
%
%     p1 = argmax (A1, B1)
%     p3 = argmax (A2, B2)
%     p5 = argmin (A3, B3)
%     p7 = argmin (A4, B4)
%
%   - The four points in between are the four possible intersections of the
%     lines A12, B12, A34 and B34.
%
%     p2 = A12 x B12
%     p6 = A34 x B34
%
%     p4 or p8 = A12 x B34
%     p8 or p8 = A34 x B12
%
% The points p1, p2 and p3 form the lower envelope of the face geomtry and
% are independent of the points p5, p6 and p7, forming the upper envelope.
%
% Although all these points will be found in a typical pillar grid, they
% cannot all be present in one face.  The intersection p4 never appear
% together with the points p3 and p5. In the same manner, the intersection
% p8 never appear together with the points p1 and p7.
%
% The numbering of points is will yield faces nodes in a clockwise
% sequence, such that normals computed in computGeometry point from cell a
% to cell b.
%
% The algorithm to obtain the intersection in the fewest numer of steps is
%
%   A) find lower envelope U p1,p2,p3.
%
%   B) find upper envelope L p5,p6,p7.
%
%   C) if upper envelope U is below lower envelope L on pillar 1,
%      find intersection point p8 and remove points p1 and p7 lying on
%      pillar 1.
%
%   D) Likewise, if U is below L on pillar 2, find intersection point p4
%      and remove points p3 and p5 lying on pillar 2.
%
% To obtain face corners in a clockwise sequence, the 8 possible points
% (on pillars or intersections) have designated positions in a row vector of
% length 8:
%
%   * The lower envelope L holds positions 1:3, with 2 being the
%     A12 x B12 intersection if present.
%
%   * The upper envelope U holds positions 5:7, with 6 being the
%     A34 x B34 intersection if present.
%
%   * An intersection between U an L that invalidates points in  positions
%     1 and 7 is placed in position 8. Positions 1 and 7 are set to NaN.
%
%   * Likewise, intersection that invalidates points in positions 3 and
%     5 are placed in position 4.  Position 3 and 5 are set top NaN.
%
% The unassigned positions are NaN.  By extraction all ~isnan positions,
% the point numbers are in the correct sequence.  In code below, J is an
% nx8 table of point numbers with NaN marking unoccupied positions.  The
% number n = size(C, 1) is the number of nonzero face cell face
% intersections.
%
%
% Input:  a(:,1), a(:,2) - point numbers for pillar 1 and 2 for all A-faces
%         b(:,1), b(:,2) - point numbers for pillar 1 and 2 for all B-faces
%         C(k,1), C(k,2) - connection between face a(C(k,1):C(k,1)+1, :)
%                          and face b(C(k,2):C(k,2)+1, :).
%         points         - point coordinates


% For each connection, extract point numbers for a-face and b-face that
% have a nonzero intersection.  The four points are numbered as follows:
%
%           1                2
%           |                |
%   pa(k,3) o----------------o pa(k,4)
%           |                |
%           |                |
%           |                |
%   pa(k,1) o----------------o pa(k,2)
%           |                |
%           |                |
%

pa = [a(C(:,1), :), a(C(:,1)+1, :)];
pb = [b(C(:,2), :), b(C(:,2)+1, :)];

n = size(pa, 1);
if n < 1
   numNodes = [];
   Corners  = [];
   return
end


if all(all (pa == pb))
   % All faces match exactly --> no faults here.
   numNodes = 4*ones(size(pa, 1),1);
   Corners    = reshape(pa', [], 1);


else
    % We only use z-coordinate of corners to determine intersections
    z  = points(:,3);
    az = reshape(z(pa), size(pa));
    bz = reshape(z(pb), size(pb));
    %  clear z;

    % Find possible points along each pillar: For each pair of cells, these
    % are the min of upper cornes and the max of lower corners along each
    % pillar
    i    = [az(:,1:2) < bz(:,1:2), az(:,3:4) > bz(:,3:4)];
    I    = pa;
    I(i) = pb(i);


    % Are there intersections between upper and lower lines?
    PP   = [pa(:,1:2); pa(:,3:4); pa(:,1:2); pa(:,3:4)];
    QQ   = [pb(:,1:2); pb(:,3:4); pb(:,3:4); pb(:,1:2)];

    i    = ( z(PP(:,1))-z(QQ(:,1)) ).*( z(PP(:,2))-z(QQ(:,2)) ) < 0;
    Q    = intersection(PP(i,:), QQ(i,:), points);

    clear PP QQ

    % Remove duplicates and store map:
    [Q, a, b]=unique(Q,'rows');                                      %#ok

    % Store point numbers of intersections in f.  Each column of f
    % represent a type of intersection:
    %
    % f(:,1)  -- A12 x B12  = p2
    % f(:,2)  -- A34 x B34  = p6
    % f(:,3)  -- A12 x B34  = p4 or p8
    % f(:,4)  -- A34 x B12  = p4 or p8
    f    = nan(n,4);
    f(i) = b(1:size( find(i(:)) ))+size(points,1);

    clear i

    % Add intersections to point list,
    points =[points;Q];

    clear Q

    % Point numbers for each face are stored in rows of J.
    % NaN is unassigned or deleted point
    J   = nan(n,8);

    % Points on pillars that may be part of face (p1, p3, p5, p7)
    % Just remove duplicate points first...
    i=I(:,1)==I(:,3); I(i,1)=nan;
    i=I(:,2)==I(:,4); I(i,2)=nan;

    J(:,[1,3,5,7]) = I(:,[1,2,4,3]);

    clear I

    % The simplest cases involve top-top and bottom-bottom line
    % intersection (#):
    %
    %   o---  A
    %   *---  B
    %
    %       1     2            1     2
    %       |     |            |     |
    %       *     |            *-----*
    %       |\    |            |     |
    %  A34  o-#---o            o-----o
    %       |  \  |            |     |
    %       |   \ |            *     |
    %       |    \|            |\    |
    %       |     * B34        | \   |
    %       |     |            |  \  |
    %  A12  o-----o            o---#-o
    %       |     |            |    \|
    %       *-----* B12        |     *
    %       |     |            |     |
    %
    %      A34 x B34           A12 x B12


    % Assign bottom-bottom and top-top intersections. In essence, the upper
    % and lower envelopes of each face are now stored in J.
    %
    J(:, 2)        = f(:,1); % p2
    J(:, 6)        = f(:,2); % p6


    % J(:, 1:3) -- bottom envelope
    % J(:, 5:7) -- top envelope

    % ---------------------------------------------------------------
    % A line in cell A intersects both lines of cell B or vice versa.
    % In other words, top and bottom envelopes intersect on left, right or
    % both sides:
    %
    %
    %   o---  A
    %   *---  B
    %      x  intersection between top and bottom envelope
    %
    %
    %        Case 1         Case2            Case 3          Case 4
    %
    %
    %       1     2        1     2          1     2          1     2
    %       |     |        |     |          |     |          |     |
    %       |     *        *     |          |     |          |     |
    %       |    /|        |\    |          *-----*          *-----*
    %       o---/-o        o-\---o          |     |          |     |
    %       |  /  |        |  \  |          |     |          |     |
    %       o-x---o        o---x-o          |     *          *     |
    %       |/    |        |    \|          |    /|          |\    |
    %       *     |        |     *          o---x-o          o-x---o
    %       |     |        |     |          |  /  |          |  \  |
    %       |     |        |     |          o-/---o          o---\-o
    %       |     |        |     |          |/    |          |    \|
    %       *-----*        *-----*          *     |          |     *
    %       |     |        |     |          |     |          |     |
    %
    %       A34 x B34      A34 x B34        A34 x B12        A34 x B12
    %       A12 x B34      A12 x B34        A12 x B12        A12 x B12
    %
    %       A12 > B34                                        B12 > A34
    %       on pillar 1                                      on pillar 1
    %
    % Combining cases 1 and 3 or cases 2 and 4 yields diamond-shaped faces:
    %
    %       1       2                     1       2
    %       |       |                     |       |
    %       *       |                     |       *
    %       |\      |                     |      /|
    %       | \     |                     |     / |
    %       *  \    |                     |    /  *
    %       |\  \   |                     |   /  /|
    %       o-x--\--o                     o--/--x-o
    %       |  \  \ |                     | /  /  |
    %       o---\--xo                     ox--/---o
    %       |    \  *                     *  /    |
    %       |     \ |                     | /     |
    %       |      \|                     |/      |
    %       |       *                     *       |
    %
    intersect_botA_topB = ~isnan(f(:,3));     %  A12 x B34
    intersect_topA_botB = ~isnan(f(:,4));     %  A34 x B12

    intersect_left = az(:,1) > bz(:,3);       %  A12 > B34 on pillar 1

    % Case 1
    ind = intersect_botA_topB & intersect_left;
    J(ind, [1, 7]) = nan;
    J(ind, 8)      = f(ind,3); % p8

    % Case2
    ind = intersect_botA_topB & ~intersect_left;
    J(ind, [3, 5]) = nan;
    J(ind, 4)      = f(ind,3); % p4

    intersect_left = bz(:,1) > az(:,3);       % B12 > A34 on pillar 1

    % Case 4
    ind = intersect_topA_botB & intersect_left;
    J(ind, [1, 7]) = nan;
    J(ind, 8)      = f(ind,4); % p8

    % Case 3
    ind = intersect_topA_botB & ~intersect_left;
    J(ind, [3, 5]) = nan;
    J(ind, 4)      = f(ind,4); % p4


    % ---------------------------------------------------------------
    % Remove repeated points arising from pinch:
    Corners  = [J, inf(size(J,1),1)]';
    Corners  = rlencode(Corners(~isnan(Corners)), 1);

    numNodes = diff([0;find(isinf(Corners))])-1;
    Corners  = Corners(~isinf(Corners));
end
end

%--------------------------------------------------------------------------

function pts = intersection(La, Lb, PTS)
%
% Find coordinates of intersection between lines.
%
% Each row in La and Lb has two point numbers that are start and
% endpoints of line segments.  Point coordinates are specified in PTS.
% This function computes [x,y,z] of intersection between each line pair in
% La,Lb.
%
z  = PTS(:,3);
za = reshape(z(La),size(La));
zb = reshape(z(Lb),size(Lb));
clear z;

% Find parameter t
t  = (zb(:,1)-za(:,1))./(diff(za, 1, 2)-diff(zb,1,2));
clear zb;

x = PTS(:,1);xa=reshape(x(La), size(La));clear x;
y = PTS(:,2);ya=reshape(y(La), size(Lb));clear y;

% Compute coordinates
pts = zeros(size(La,1),3);
pts(:,1) = diff(xa,1,2).*t + xa(:,1);
pts(:,2) = diff(ya,1,2).*t + ya(:,1);
pts(:,3) = diff(za,1,2).*t + za(:,1);
end

%--------------------------------------------------------------------------

function C = findConnections(za, zb)
%
%           (1)                       (2)
%            |                         |
%  za(i+1,1) o-------------------------o za(i+1,2)
%            |                         |
%  zb(j+1,1) * * * *                   |
%            |       * * * *           |            ^
%            |               * * * *   |            |
%            |                       * * zb(j+1,2)  |z positive
%            |                         |            |
%            |                         |            |
%            |                         |
%    za(i,1) o-------------------------o za(i,2)
%            |                         |
%    zb(j,1) * * * * * * * *           |
%            |               * * * * * * zb(j,2)
%            |                         |
%
% Given:
% vectors of point numbers a and b where a(:,1), b(:,1) refer to pillar (1)
% and a(:,2), b(:,2) refer to pillar (2).
% Each cell a_i is assumed to lie between lines a(i,:) and a(i+1,:) and
% each cell b_j is assumed to lie between the lines b(j,:) and b(j+1,:).
%
% Walk along the stack of cells a_i  and find all connections to cell b_j.

C = zeros(0,2);
j1 = 1;
j2 = 1;
for i = 1 : size(za, 1) - 1

    j = min(j1, j2);  % Largest j where both
                      % zb(j,1) < za(i,1) and
                      % zb(j,2) < za(i,2)

    % While  b(j,:) (bottom of cell b_j) is below a(i+1,:) (top of cell
    % a_i)
     while any(zb(j,:) < za(i+1,:), 2)

        % Precise check to avoid adding pinched layers
       if  doIntersect(za(i,:), za(i+1,:), zb(j,:), zb(j+1,:))
           C = [C; i, j]; %#ok

       end

       % Update candidates for next start of j iteration
       if zb(j,1) < za(i+1,1), j1 = j; end
       if zb(j,2) < za(i+1,2), j2 = j; end

       j=j+1;
     end
end
end

%--------------------------------------------------------------------------

function val = doIntersect(za1, za2, zb1, zb2)
%
%  Does cells given by points (a11,a12,a21,a22) and (b11,b12,b21,b22)
%  given along pillar (1) and pillar (2) have an nonzero intersection?
%  We need only the z-ccordinate to check:
%
%        |                 |
% zb1(1) *                 |
% za1(1) o-*---------------o za1(2)
%        |. .*             |
%        | . . *           |
%        |. . . .*         |
%        | . . . . *       |
% za2(1) o-----------*-----o za2(2)
%        |             *   |
%        |               * |
%        |                 * zb1(2)
%        |                 |
% zb2(1) * * * * * * * * * * zb2(2)
%        |                 |
%       (1)               (2)
%

val = overlap (za1(1), za2(1), zb1(1), zb2(1)) | ... % simple connection
      overlap (za1(2), za2(2), zb1(2), zb2(2)) | ... % simple connection
      (za1(1)-zb1(1)).*(za1(2)-zb1(2)) < 0     | ... % Le fix speciale
      (za2(1)-zb2(1)).*(za2(2)-zb2(2)) < 0;          % Le fix speciale

if all(za1-za2==0), val=false; end
if all(zb1-zb2==0), val=false; end
end

%--------------------------------------------------------------------------

function val = overlap(xa1, xa2, xb1, xb2)
%  Do the intervals overlap?
%
%         xa1(               )xa2
%  -----------------------------------
%   xb1(               )xb2
%
val = max(xa1, xb1) < min(xa2, xb2);
end
