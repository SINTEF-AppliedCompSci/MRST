function [g_top, g] = topSurfaceGrid(g, varargin)
%Make a hybrid grid of the top surface of a 3D grid having a logical
%ijk-index (e.g., a corner-point grid or a logically Cartesian grid)
%
% SYNOPSIS:
%   [Gt,G] = topSurfaceGrid(G)
%
% PARAMETERS:
%   G       - 3D grid as described by grid_structure.
%
%  'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%            supported options are:
%              'Verbose' -- boolean, whether to show log messages or not.
%                           Default value: mrstVerbose
%              'AddFaults' -- boolean indicating if faults should be added
%                             to the resulting 2D model.
%
% RETURNS:
%   Gt - structure representing the top-surface grid. The structure
%      consists of two parts:
%
%      (i) the surface projected onto the xy-plane represented as a
%          standard 2D grid that follows the standard grid conventions as
%          outlined in grid_structure,
%
%      (ii) a set of extra fields that represent the surface elevation, the
%          3D surface normals, the elevation of the bottom surface, the
%          height of the surface, the columns of volumetric cells attached
%          to each cell (quadrilateral patch) in the 2D grid, etc.
%
%      Altogether, the hybrid grid structure contains the following fields:
%
%      - cells   -- A structure specifying properties for each individual
%                   cell in the grid. See CELLS below for details.
%
%      - faces   -- A structure specifying properties for each individual
%                   face in the grid. See FACES below for details.
%
%      - nodes   -- A structure specifying properties for each individual
%                   node (vertex) in the grid. See NODES below for details.
%
%      - columns -- A structure specifying properties of each individual
%                   vertically-averaged column. See COLUMNS below for
%                   details.
%
%      - type    -- A cell array of strings describing the history of grid
%                   constructor and modifier functions through which this
%                   particular grid has been defined.
%
%   G -  3D grid with connected columns and no disconnected cells,
%      corresponding to Gt. Important to use this 3D grid when doing
%      comparisons or if 2D data is mapped back to 3D grid.
%
%
% FIELDS in a 2D top-surface grid:
%
%   COLUMNS - Cell structure Gt.columns that represents the column of
%   cells in the 3D grid that are attached to each cell in the top-surface
%   grid. Containts the fields:
%
%     - cells  -- A G.cells.num-by-1 array of cell indices from the
%                 3D grid, which together with Gt.cells.columnPos
%                 define the columns attached to the top surface.
%
%     - dz     -- A G.cells.num-by-1 array of heights of the cells in the
%                 3D grid defined as the vertical difference between
%                 the centroids of the top and bottom surfaces of each
%                 cell (faces labelled 5 and 6).
%
%     - z      -- A G.cells.num-by-1 array of the z-component of the
%                 centroid of the bottom surface of each cell, defined
%                 relative to the height of the centroid of the top
%                 surface, i.e., a column-wise cumsum of columns.dz.
%
%
%   CELLS - Cell structure Gt.cells, as described by grid_structure
%   but with the extra fields:
%
%     - columnPos -- Indirection map of size [num+1,1] into the
%                  'columns.cells' array. Specifically, the cells in the
%                  column of the 3D grid corresponding to cell 'i' in the
%                  top-surface grid are found in the submatrix
%
%                    Gt.columns.cells(columnPos(i):columnPos(i+1)-1,:)
%
%                  The number of cells in the each column can be computed
%                  by the statement DIFF(columnPos)
%
%     - map3DFace -- A cells.num-by-1 array that maps between a cell in the
%                  top-surface grid and the corresponding top face in the
%                  underlying 3D grid, i.e., cell 'i' corresponds to face
%                  number Gt.cells.map3DFace(i).
%
%     - ij      -- A map of size [num,1] giving logical indices of the
%                  cells in the top-surface grid
%
%     - H       -- A  array of formation heights, one for each cell
%
%     - normals -- normal vector of the surface represented by each cell
%
%     - z       -- A cells.num-by-1 array of the third coordinate of the
%                  cell centroids in R^3. The first two coordinates are
%                  given in Gt.cells.
%
%
%   FACES - Cell structure Gt.faces, as described in grid_structure
%   but with the following extra field created by a subsequent call to
%   computeGeometryVE:
%
%     - z       -- A faces.num-by-1 array of the third coordinate of the
%                  face centroids in R^3. The first two coordinates are
%                  given in Gt.faces.centroids.
%
%
%   NODES - Cell structure Gt.nodes, as described in grid_structure
%   but with the following extra field created by a subsequent call to
%   computeGeometryVE:
%
%     - z       -- A nodes.num-by-1 array of the third coordinate of the
%                  node positions in R^3, i.e., the z-coordinate of the
%                  top surface of the original grid. The first two
%                  coordinates are given in Gt.nodes.coords.
%
%
% COMMENTS:
%   PLEASE NOTE:
%   The current implementation does not honour faults and internal
%   boundaries from 3D grid, i.e., it is possible to get false connections
%   over faults if there is no jump in the ij index over the fault in the
%   3D grid.


%{
#COPYRIGHT#
%}

% $Date: 2012-10-05 15:27:55 +0200 (Fri, 05 Oct 2012) $
% $Revision: 10008 $

opt = struct('Verbose', mrstVerbose, ...
             'AddFaults', false,'add_cellnodes',true);
opt = merge_options(opt, varargin{:});

verbose = opt.Verbose;

%%%% Create 2D grid
% remove inactive cells
A = false(g.cartDims);
A(g.cells.indexMap) = true;
B = cumsum(A,3);
C = (B==0);
% find connected columns, remove all cells beneath possible gap
activeCells = (cumprod(A+C, 3)-C);
% remove all disconnected columns
g_old = g;

if ~all(activeCells(:))
   mlist = mrstModule;
   %mrstModule add libgeometry;
   g = removeCells(g, cart2active(g, find(~activeCells)));
   g = remove_disc_cells(g, verbose);
   g = computeGeometry(g);
   mrstModule('reset',mlist{:});
end
A = false(g.cartDims);
A(g.cells.indexMap) = true;
active = any(A, 3);

g_top = cartGrid(g.cartDims(1:2));
g_top.faces.tag=(1:g_top.faces.num)';
g_top = removeCells(g_top, find(~active)');

cartFaceMap=g_top.faces.tag;

g_top.faces=rmfield(g_top.faces,'tag');

B = cumsum(A,3);

topCells3D=find(B(g.cells.indexMap)==1);

%%%% Find top cells/faces
topCFaces = g.cells.faces(mcolon(g.cells.facePos(topCells3D), ...
                                 g.cells.facePos(topCells3D+1)-1),:);

topFaces3D = topCFaces(topCFaces(:,2) == 5,1);
% check that all faces are boundary faces:
assert(all(any(g.faces.neighbors(topFaces3D,:)==0,2)))

% use ijk index to find cell number in 2D grid
[ijk{1:3}] = ind2sub(g.cartDims, g.cells.indexMap);
ijk = [ijk{:}];

% cells numbers in the old 2D grid (before call to removeCells)
oldCellIx = sub2ind(g.cartDims(1:2), ijk(topCells3D,1), ijk(topCells3D,2));

% cell numbers in current 2D grid
cells2D = cart2active(g_top, oldCellIx);

[trash, ix] = sort(cells2D); %#ok backwards compatability

% renumber according to cells in 2D grid:
topCells3D = topCells3D(ix);
topFaces3D = topFaces3D(ix);

%%%% Averaging of node coordinates in the top-surface grid
% Find corresponence between node number in 3D and 2D grid. One node in 2D
% grid can correspond to two nodes in 3D grid if the top cells are not
% neighbors. Must therefore accumulate and average contributions from these
% nodes to make nodes for top surface grid.
%
%  3D-grid viewed in x-z plane, top surface marked with '
%
%          *'''''''*
%          |       |
%          | cell2 |
%   *''''''*-------*
%   |cell1 |
%   *------*
%
%  Corresponding nodes in top surface grid, the node in the middle has been
%  averaged:
%
%                  *
%          *
%   *
%
%---------------------------------------------------------------------

faceNodes = mcolon(g.faces.nodePos(topFaces3D), ...
                   g.faces.nodePos(topFaces3D+1)-1);
% cn_top = [cell local_node_num global_node_num]
cn_top = cellNodes(g_top);

% find node numbers in 3D grid (take top nodes topfaces)
nodes3D = g.faces.nodes(faceNodes);


% Find correspondace between nodes in 3D grid and 2D grid.
% Local numbering of nodes differ:
%
%   faceNodes:      cellNodes:
%
%    4     3        3     4
%    *-----*        *-----*
%    |     |        |     |
%    |     |        |     |
%    *-----*        *-----*
%    1     2        1     2
%
%  Must therefore use index to match correct nodes:
%ix = repmat([1 2 4 3]', g_top.cells.num,1)+(cn_top(:,1)-1)*4; 
% NB: the above does not work when the tiling order of the top faces in the 
% original grid does not match the indexing of the cells in the top surface
% grid.
% Proposed fix in the two codelines below:

% Determine how to orient the 2D cell relative to the 3D face in order to
% preserve the logical ordering (ensuring that neighbor cells share the
% correct nodes)
rco = relative_corner_order(nodes3D); 
ix = repmat(rco, g_top.cells.num,1) + (cn_top(:,1)-1)*4;

nodes2D = cn_top(ix,3);

g_top.nodes.coords(:,1) = accumarray(nodes2D, ...
                        g.nodes.coords(nodes3D, 1))./accumarray(nodes2D,1);
g_top.nodes.coords(:,2) = accumarray(nodes2D, ...
                        g.nodes.coords(nodes3D, 2))./accumarray(nodes2D,1);
g_top.nodes.z = accumarray(nodes2D, ...
                        g.nodes.coords(nodes3D, 3))./accumarray(nodes2D,1);

% NB: Ideally, we should have 4 z-coordinates per grid cell and only use
% average coordinates if columns are connected. An improvement would be to
% make new nodes if we have a fault without connectivity.

%%% Compute column information

% sort cells in 3D grid by wrt to columns (j, i, and k-index)
% i.e: fix j (x-dir) and sort on i index before we sort on k and numbering
% of columns will be the same as the numbering of the top cells:
%  ---------
%  | 3 | 4 |
%  ---------
%  | 1 | 2 |
%  ---------

[mat, cell_ix] = sortrows(ijk, [2 1 3]);
columns.cells = cell_ix;

% compute dz for each 3D cell
columns.dz = compute_dz(g, cell_ix);

% find index to first cell in each column in columns.cells
% stolen from rlencode: compare differences in position to find run length.
columnPos = 1 + ...
   [0; [find(any(mat(1:end-1,1:2) ~= mat(2:end,1:2),2)); size(mat,1)]];

% Compute column height H:
% NB: inaccurate if dz is inaccurate
colNo = rldecode(1:g_top.cells.num, columnPos(2:end)-columnPos(1:(end-1)),2).';
H = accumarray(colNo, columns.dz);

%%% Add fields to struct

g_top.columns = columns;
g_top.cells.columnPos = columnPos;

% save mappings to 3D grid
g_top.cells.map3DFace = topFaces3D;
g_top.cells.ij    = [ijk(topCells3D,1) ijk(topCells3D,2)];
g_top.parent=g;

%%% Add faults
% Not well tested, but is analogous to what is done in examples
if opt.AddFaults && any(g.faces.tag)
    [g_top c_fault] = add_faults(g, g_top, ijk);
    g_top = computeGeometryVE_2D(g_top);
    
    % Adjust the z values of the corresponding cells
    g_top.cells.z(c_fault) = g.faces.centroids(topFaces3D(c_fault),3);
    
    % Note that no extra nodes have been added to the grid, and g.faces.z
    % has not been adjusted, but this (probably) only matters for plotting
    % the grid. g_top gets an internal boundary across sealing faults.
else
    g_top = computeGeometryVE_2D(g_top);
    
end
if(any(g_top.cells.volumes<0))
        disp('Wrong sign of volumes take abs')  
        g_top.cells.volumes=abs(g_top.cells.volumes);
end
if(any(g_top.faces.areas<0))
        disp('Wrong sign of areas take abs')  
        g_top.faces.areas=abs(g_top.cells.areas);
end

%%% Extract geometry from 3D grid and compute remaining fields
g_top.cells.H = H;
g_top.cells.normals = g.faces.normals(topFaces3D, :);

% g_top = computeGeometryVE(g_top);
if(any(diff(g_top.cells.columnPos)>1))   
    g_top.columns.z = cumulativeHeight(g_top);
else
    g_top.columns.z = g_top.cells.z;
end
g_top.type = [g_top.type, mfilename];

%%% Add function handles to pass replacement functions to MRST default
% solvers. This enables use of solveIncompFlow with top grids.
g_top.grav_pressure = @(g, omega) gravPressureVE_s(g, omega);

g_top.primitives = @primitivesMimeticVE_s;

% Store parent grid
g_top.parent = g;

% Report if 3D grid has been changed, important to use new 3D grid for
% cases where 2D/3D comparisons is done, or where 2D data is mapped back to
% 3D grid.
if nargout == 1 && g_old.cells.num ~= g.cells.num
   warning(msgid([' 3D grid has been changed. ' , ...
      'Use 3D grid returned from topSurfaceGrid.m']));
end
g.type = [g.type, mfilename];

if(opt.add_cellnodes)
   % Speed up plotting
  g_top.cells.sortedCellNodes = getSortedCellNodes(g_top);
  % Temporary addition for the gap between initial co2lab release and mrst
                                % 2013a.
  g_top.cells.cellNodes = g_top.cells.sortedCellNodes;
end
% this map is used for transfering face data on g_top to the face data on
% the underlying Cartesian grid. 
g_top.faces.cartFaceMap=cartFaceMap;
end

% ----------------------------------------------------------------------------
function rco = relative_corner_order(co)
% Determine the relative corner order that should be used when reading
% nodes from 3D faces to 2D cells.  As described in the main routine, the
% default permutation to apply per-cell is [1 2 4 3].  However, this
% presumes that cells are ordered such that row-wise, the 3D face
% corresponding to cell (j+1) immediately follows (rather than precedes)
% the one corresponding to cell (j), and that likewise, the 3D faces
% corresponding to cells in row (i) succeeds (rather than precedes) those
% of row (i+1).  If these assumptions do not hold, we must modify the
% permutation to use accordingly.

rco = [1 2 4 3]'; % should be true for most grids

% Function to determine whether 3D face 2 lies to the 'right' or 'left'
% of 3D face 1. (assuming that each face has 4 corners, and that there
% are at least two faces present)
next_cell_comes_before = @(ixs) (ixs(1) == ixs(6)) && (ixs(4) == ixs(7));
    
% Function to determine whether the second row of 3D faces comes 'before'
% or 'after' the first row.  We assume this happens if the lower-left
% corner of the first cell is later repeated (ignoring the two first
% nodes of the second cell), _or_ if the lower-right node is later
% repeated (ignoring the first node of the second cell)
next_row_comes_before = ...
   @(ixs) ~isempty(find(ixs(7:end)==ixs(1), 1, 'first')) || ...
   ~isempty(find(ixs(6:end)==ixs(2), 1, 'first'));
    
if (next_row_comes_before(co))
   rco = flipud(rco);
end

if (next_cell_comes_before(co))
   rco(1:2) = flip(rco(1:2));
   rco(3:4) = flip(rco(3:4));
end
end


function dz = compute_dz(g, cell_ix)
% Assume that all cells have both top and bottom face.
z_top = g.faces.centroids(g.cells.faces(g.cells.faces(:,2)==6,1),3);
z_bottom = g.faces.centroids(g.cells.faces(g.cells.faces(:,2)==5,1),3);

dz = abs(z_top-z_bottom);
dz = dz(cell_ix);
end

function g = remove_disc_cells(G, verbose)
% Make 3D grid suitable for construction of top surface grid.
% Remove disconnected cells from the 3D model

ix = all(G.faces.neighbors~=0, 2);
I  = [G.faces.neighbors(ix,1);G.faces.neighbors(ix,2)];
J  = [G.faces.neighbors(ix,2);G.faces.neighbors(ix,1)];
N  = double(max(G.faces.neighbors(:)));
A  = sparse(double(I),double(J),1,N,N)+speye(N);
clear ix I J
[a,b,c,d]=dmperm(A);                                           %#ok
clear A

if numel(c) > 2,
   dispif(verbose, '\nGrid has %d disconnected components\n', ...
          numel(c)-  1);

   for i = 1:numel(c) - 1,
      g(i)  = extractSubgrid(G, a(c(i):c(i+1)-1)); %#ok
      sz(i) = g(i).cells.num;                      %#ok
   end

   for i = 1 : numel(c)-1,
      g(i).cartDims = G.cartDims;       %#ok
      g(i).type     = { mfilename };    %#ok
   end

   [i,i] = sort(-sz); %#ok<ASGLU>
   G     = g(i);
end

if numel(c) > 2,
dispif(verbose, ['\n 3D grid has %d disconnected components, ', ...
                   'choosing the largest\n'], numel(c)-  1);
end
g = G(1);

end

function [g_top c_full] = add_faults(g, g_top, ijk)

f_fault = (g.faces.neighbors(:,1) == 0 | ...
   g.faces.neighbors(:,2) == 0) & g.faces.tag~=0;
c_fault = g.faces.neighbors(f_fault,:);
c_fault = c_fault(:);
c_fault = c_fault(c_fault>0);

% Find i/j indices of cells belonging to fault
ij_f = unique(ijk(c_fault,1:2), 'rows');

% Find all cells with those i/j indices
f_c = ismember(ijk(:,1:2), ij_f, 'rows');

% Remove any i/j pairs corresponding to columns where not all fine cells
% belong to the fault.
c_partial = ijk(setdiff(find(f_c), c_fault), 1:2);
c_full    = setdiff(ij_f, c_partial, 'rows');
c_full = find(ismember(g_top.cells.indexMap, ...
   sub2ind(g_top.cartDims, c_full(:,1), c_full(:,2))));

% If a face is between two full faulted cells, remove
face_fault = find(ismember(g_top.faces.neighbors(:,1), c_full) & ...
   ismember(g_top.faces.neighbors(:,2), c_full));

% Make an internal boundary
if any(face_fault)
    g_top = makeInternalBoundary(g_top, face_fault);
end
end
