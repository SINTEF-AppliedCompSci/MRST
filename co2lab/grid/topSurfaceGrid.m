function [Gt, G, transMult] = topSurfaceGrid(G)

   %% Identify columns in layered 3D grid; remove unused cells from G
   % cells in G that will not contribute to the top surface grid 
   % will be removed
   [active_cols, col_cells, G] = identify_column_cells(G);
   
   %% Identify connectivity between columns
   % Absent the presence of faults, this task would be trivial on cartesian
   % grids.  However, a general algorithm is needed in order to handle
   % general layered grids.
   col_adj = determine_column_connectivity(G, active_cols, col_cells);
   
   %% Identify top faces
   tcells = col_cells(1, active_cols); % top cells
   tfaces = identify_lateral_faces(G, tcells); % top faces
   
   %% Construct the 2D grid
   % 'orient' informs about edge orientation, and will be used later when
   % identifying cell neighbors.  However, this has to wait until unwanted
   % discontinuities have been removed (see a couple of lines further down)
   [Gt, orient] = construct_grid_from_top_faces(G, tcells, tfaces);
   
   %% Identify fault edges
   ffaces = identify_fault_faces(G, Gt, tfaces);
   
   %% Stitch-up surface discontinuities
   % update fault faces in case grid faces are merged
   [Gt, ffaces] = stitch_surface_discontinuities(Gt, ffaces, col_adj, orient);
   
   %% Compute and add column information
   [Gt.columns, Gt.cells.columnPos, Gt.cells.H] = ...
       compute_column_info(G, col_cells(:, active_cols), tfaces);

   %% Compute transmissibility multipliers for partially overlapping faults
   transMult = computeFaultTrans(Gt, ffaces);
   
   %% Fill in the remaining fields
   Gt.parent          = G;
   Gt.grav_pressure   = @(G, omega) gravPressureVE_s(G, omega);
   Gt.primitives      = @primitivesMimeticVE_s;
   
   if isfield(G, 'cartDims')
      % also add face tags
      fnum = length(Gt.cells.faces);
      tag_seq = determine_facetag_sequence(G);
      Gt.cells.faces = [Gt.cells.faces, repmat(tag_seq, fnum/4, 1)];
   end
   
   %% Compute geometry
   Gt = compute_geometry(Gt);
end

% ----------------------------------------------------------------------------

function adj_mat = determine_column_connectivity(G, active_cols, col_cells)

   num_cols = length(active_cols);
   num_cells = G.cells.num;
   
   %% making map from cells to column
   cmap = [col_cells(:), ...
           reshape((repmat(1:size(col_cells, 2), size(col_cells, 1), 1)), [], 1)];
   cmap = cmap(cmap(:, 1) ~= 0, :); % removing lines corresponding to padding

   % matrix mapping cells to columns
   cell_col_mat = spones(accumarray(cmap, 1, [num_cells, num_cols], @sum, 0, true));
   
   %% matrix mapping connectivity between cells (i.e. neighbors)
   neighs = G.faces.neighbors;
   neighs = neighs(prod(neighs, 2) ~= 0, :); % keep only internal faces
   
   cell_cell_mat = accumarray(neighs, 1, [num_cells, num_cells], @sum, 0, true);
   cell_cell_mat = spones(cell_cell_mat + cell_cell_mat');
   
   %% Matrix mapping columns to columns
   adj_mat = spones(cell_col_mat' * cell_cell_mat * cell_col_mat);
   adj_mat(logical(speye(num_cols))) = 0; % set diagonal to zero
   
   %% Removing empty columns
   acolix = find(active_cols);
   adj_mat = adj_mat(acolix, acolix);
end

% ----------------------------------------------------------------------------

function transMult = computeFaultTrans(Gt, ffaces)
   
   neigh_cols = Gt.faces.neighbors(ffaces, :);
   int_ind    = prod(neigh_cols, 2) ~= 0;
   neigh_cols = neigh_cols(int_ind, :); % only internal pairs
   ffaces     = ffaces(int_ind);
   
   col_tops    = Gt.cells.z(neigh_cols);
   col_heights = Gt.cells.H(neigh_cols);
   col_bots    = col_heights + col_tops;
   
   % Determine vertical extent of column overlap for  each neighbor pair
   min_heights = min(col_heights, [], 2);
   overlap = max(min(col_bots, [], 2) - max(col_tops, [], 2), 0);
   
   transMult = ones(Gt.faces.num, 1);
   transMult(ffaces) = overlap ./ min_heights;
   
end

% ----------------------------------------------------------------------------

function cfaces = find_connecting_faces(Gt, nodes)

   fnum = Gt.faces.num; % total number of faces
   n_ind = false(Gt.nodes.num, 1);
   n_ind(nodes) = true;
   fnodes = Gt.faces.nodes(mcolon(Gt.faces.nodePos(1:fnum), ...
                                  Gt.faces.nodePos(2:fnum+1)-1));
   fnode_ind = reshape(n_ind(fnodes), 2, []);
   
   cfaces = find(prod(fnode_ind))'; % all faces whose endpoints are both in 'nodes'
   
end

% ----------------------------------------------------------------------------

function Gt = compute_geometry(Gt)

   Gt = computeGeometryVE_2D(Gt);
   if(any(Gt.cells.volumes < 0))
      disp('Wrong sign of volumes. Changing to absolute values.')  
      Gt.cells.volumes = abs(Gt.cells.volumes);
   end
   if(any(Gt.faces.areas<0))
      disp('Wrong sign of areas. Changing to absolute values.')  
      Gt.faces.areas=abs(Gt.cells.areas);
   end
end

% ----------------------------------------------------------------------------

function [cols, col_pos, H] = compute_column_info(G, col_cells, tfaces)

   col_sizes = sum(col_cells~=0)';
   col_pos = cumsum([1; col_sizes]);

   start_ix = cumsum([1;repmat(size(col_cells, 1), size(col_cells, 2), 1)]);
   start_ix = start_ix(1:end-1);
   end_ix = start_ix + col_sizes - 1;
   
   cols.cells = col_cells(mcolon(start_ix, end_ix))';
   
   [top_faces, bot_faces] = identify_lateral_faces(G, cols.cells);

   ref_z = rldecode(G.faces.centroids(tfaces, 3), col_sizes);
   cols.z = G.faces.centroids(bot_faces, 3) - ref_z;   
   cols.dz = G.faces.centroids(bot_faces, 3) - G.faces.centroids(top_faces, 3);
   H = accumarray(rldecode((1:length(col_sizes))', col_sizes'), cols.dz);
   
end

% ----------------------------------------------------------------------------

function [Gt, ffaces] = stitch_surface_discontinuities(Gt, ffaces, col_adj, orient)

   % find the adjacencies currently implied by topology of Gt
   neighs = find_neighbors(Gt.cells.faces, diff(Gt.cells.facePos), orient);
   neighs = neighs(prod(neighs, 2) ~= 0, :); % only consider internal faces
   Gt_cell_adj = accumarray(neighs, 1, size(col_adj), @sum, 0, true);
   Gt_cell_adj = spones(Gt_cell_adj + Gt_cell_adj');
   mcon = col_adj - Gt_cell_adj; % matrix representing still missing connections
   assert(min(mcon(:)) >= 0); % there should only be zero and one values in 'mcon'
   [I, J] = ind2sub(size(mcon), find(triu(mcon))); % neighbors not yet accounted for by Gt
   
   if ~isempty(I)
      [Gt, ffaces] = stitch_neighbors(Gt, ffaces, I, J);
   end
   
   % With topology now in place, we can finally compute cell neighbors
   Gt.faces.neighbors = ...
       find_neighbors(Gt.cells.faces, diff(Gt.cells.facePos), orient);
end

% ----------------------------------------------------------------------------

function [Gt, ffaces] = stitch_neighbors(Gt, ffaces, I, J)

   % Each cell neighbor pair will require merging of two node pairs.  We
   % identify the two nodepairs whose members are closest to each other.
   N1 = padded_nodelist(Gt, I, NaN); % padded list of nodes per cell in I
   N2 = padded_nodelist(Gt, J, NaN); % padded list of nodes per cell in J
   
   comb = combination_pairs(1:size(N1,2), 1:size(N2, 2));
   dist2 = nan(length(I), size(comb, 2));
   r = 1;
   for i = comb
      dist2(:,r) = sum((Gt.nodes.coords(N1(:, i(1)),:) - ...
                        Gt.nodes.coords(N2(:, i(2)),:)).^2, 2);
      r = r+1;
   end
   [~, ix] = sort(dist2,2);
   
   N1 = N1'; N1_rix = cumsum([0; repmat(size(N1, 1), size(N1,2)-1, 1)]);
   N2 = N2'; N2_rix = cumsum([0; repmat(size(N2, 1), size(N2,2)-1, 1)]);
      
   npairs = [N1(N1_rix + comb(1, ix(:, 1))'), ...
             N2(N2_rix + comb(2, ix(:, 1))'); ...  % first pair of closest nodes
             N1(N1_rix + comb(1, ix(:, 2))'), ...
             N2(N2_rix + comb(2, ix(:, 2))')];     % second pair of closest nodes
   
   % remove duplicates.  This will now become the list of node pairs to merge
   npairs = unique(sort(npairs, 2), 'rows');
   
   % Compute and set new coordinates for the nodes to merge
   new_coords = (Gt.nodes.coords(npairs(:, 1),:) + Gt.nodes.coords(npairs(:,2),:))/2;
   new_z = (Gt.nodes.z(npairs(:,1)) + Gt.nodes.z(npairs(:,2)))/2;
   
   Gt.nodes.coords([npairs(:,1); npairs(:,2)], :) = [new_coords; new_coords];
   Gt.nodes.z([npairs(:,1); npairs(:,2)]) = [new_z; new_z];

   % Merging nodes and updating indexing
   [unodes, fnodes_ixs] = uniqueNodes([Gt.nodes.coords, Gt.nodes.z]);
   Gt.nodes.num = size(unodes(:,1));
   Gt.nodes.coords = unodes(:,1:2);
   Gt.nodes.z = unodes(:, 3);

   % Updating Gt faces
   Gt.faces.nodes = fnodes_ixs(Gt.faces.nodes);
   [new_fnodes, nfix] = uniqueNodes(reshape(Gt.faces.nodes, 2, [])');
   Gt.faces.nodes = reshape(transpose(new_fnodes), [], 1);
   Gt.faces.num = length(Gt.faces.nodes)/2;
   Gt.faces.nodePos = (1:2:(numel(Gt.faces.nodes)+1))';
   
   % Updating Gt cells and fault face indexing
   Gt.cells.faces = nfix(Gt.cells.faces);
   ffaces = nfix(ffaces); % update indices to fault faces as well
   
end

% ----------------------------------------------------------------------------
function c = combination_pairs(l1, l2)
   % return list of 2D vecs representing all possible combinations of
   % one element from l1 and one element from l2
   [A, B] = meshgrid(1:numel(l1), 1:numel(l2));
   c = [l1(A(:)); l2(B(:))];
end


% ----------------------------------------------------------------------------

function padded =  padded_nodelist(Gt, cell_ix, fill_val)
   faces = Gt.cells.faces(mcolon(Gt.cells.facePos(cell_ix), ...
                                 Gt.cells.facePos(cell_ix+1)-1));
   nodes = Gt.faces.nodes(mcolon(Gt.faces.nodePos(faces), ...
                                 Gt.faces.nodePos(faces+1)-1));
   % An assumption here is that each face has exactly two nodes, which should
   % always be true in 2D
   num_faces      = Gt.cells.facePos(cell_ix+1) - Gt.cells.facePos(cell_ix);
   ctrl           = rldecode((1:numel(cell_ix))', 2 * num_faces); 
   nodes_cells    = unique([nodes, ctrl], 'rows', 'stable');
   [~, num_nodes] = rlencode(nodes_cells(:,2));
   assert(all(num_nodes == num_faces)); % just to be sure our assumption hold
   
   padded = fill_val * ones(max(num_nodes), numel(cell_ix));

   row_start_ix = cumsum([1; max(num_nodes) * ones(numel(cell_ix)-1, 1)]);

   padded(mcolon(row_start_ix, row_start_ix + num_nodes - 1)) = nodes_cells(:,1);

   padded = padded';
end

% ----------------------------------------------------------------------------

function [ffaces, fnodes] = identify_fault_faces(G, Gt, top_faces)

   fnodes = []; ffaces = [];
   if ~isfield(G.faces, 'tag')
      % no fault information present
      return
   end
             
   % Identify 3D nodes linked to top faces in 3D grid
   top_nodes_3D = false(G.nodes.num, 1);
   top_nodes_3D(G.faces.nodes(mcolon(G.faces.nodePos(top_faces), ...
                                     G.faces.nodePos(top_faces+1)-1))) = true;
   
   % Identify 3D nodes associated with fault, and also the top node subset
   fault_nodes_3D = false(G.nodes.num, 1);
   fault_nodes_3D(G.faces.nodes(mcolon(...
       G.faces.nodePos(logical(G.faces.tag)), ...
       G.faces.nodePos(circshift(logical(G.faces.tag), 1))))) = true;

   top_fault_nodes = top_nodes_3D & fault_nodes_3D;
   
   % Searching for the corresponding nodes in the 2D grid
   clist = [[Gt.nodes.coords, Gt.nodes.z], (1:length(Gt.nodes.coords))'; ...
            G.nodes.coords(top_fault_nodes,:), inf(sum(top_fault_nodes), 1)];
   
   clist = sortrows(clist);
   dup_ix = all(diff(clist(:,1:3)) == 0, 2);
   
   fnodes = clist(dup_ix, 4);
   
   % Identifying the corresponding faces
   ffaces = find_connecting_faces(Gt, fnodes);
end

% ----------------------------------------------------------------------------

function [Gt, orient] = construct_grid_from_top_faces(G, tcells, tfaces)
   
   % get a list with (nonunique) coordinates for involved nodes
   fnodes_start = G.faces.nodePos(tfaces);
   fnodes_end   = G.faces.nodePos(tfaces + 1); % really one-past-end
   fnodes_ixs = G.faces.nodes(mcolon(fnodes_start, fnodes_end-1));
   f_nodenum = fnodes_end - fnodes_start;
      
   coords = G.nodes.coords(fnodes_ixs, :);
   
   % Identify a set of unique nodes, and change fnodes_ixs to refer to these instead
   [unodes, fnodes_ixs] = uniqueNodes(coords); 
   
   %% constructing node field of Gt
   Gt.nodes = struct('num', numel(unodes), 'coords', unodes(:,1:2), 'z', unodes(:,3));
   
   %% Constructing face field (i.e. edges) of Gt
   % We start with a set of edges containing duplicates
   edges = [fnodes_ixs, [fnodes_ixs(2:end);0]];
   last_node_ix = cumsum(f_nodenum); % last node per face
   first_node_ix = [1; last_node_ix(1:end-1)+1]; % first node per face
      
   edges(last_node_ix, 2) = edges(first_node_ix, 1);
   
   % orienting each edge according to node index, but keeping track of original
   % (correct) orientation
   [edges, orient] = sort(edges, 2);
   
   % Removing duplicates from edge set
   [uedges, faces] = uniqueNodes(edges);
   
   % Setting the 'faces' struct.  The 'neighbor' field will have to wait
   % until unwanted discontinuities have been removed from the grid.
   Gt.faces = struct('num', size(uedges, 1), 'nodePos', (1:2:numel(uedges)+1)', ...
                     'nodes', reshape(uedges', numel(uedges), []), ...
                     'z', mean([Gt.nodes.z(uedges(:,1)), Gt.nodes.z(uedges(:,2))],2));
   
   %% Constructing cells of Gt (i.e. corresponding to top faces in the 3D grid)
   Gt.cells = struct('num', numel(tfaces), 'faces', faces, ...
                     'facePos', [1; last_node_ix+1], 'z', G.faces.centroids(tfaces, 3));
   
   %% Other immediate fields
   Gt.griddim = 2;
   Gt.type = [G.type, {'topSurfaceGrid'}]; % @ change to mfilename
   Gt.cells.map3DFace = tfaces;
   Gt.cells.normals = G.faces.normals(tfaces, :);
   
   if isfield(G, 'cartDims')
      Gt.cartDims = G.cartDims(1:end-1);
      [I, J, ~] = ind2sub(G.cartDims, G.cells.indexMap(tcells));
      Gt.cells.ij = [I, J];
      Gt.cells.indexMap = sub2ind(Gt.cartDims, I, J);
   end
   
   % Forwarding the essential information about edge orientation 
   orient = orient(:,1);
   orient(orient==2) = -1;
end


% ----------------------------------------------------------------------------
function neigh = find_neighbors(cell_faces, num_cell_faces, orient)

   num_cells = numel(num_cell_faces);
   cellInx = rldecode((1:num_cells)', num_cell_faces);

   mat = accumarray([cell_faces, cellInx], orient, [], [], [], true);

   [r1, n1] = find(mat == 1);
   [r2, n2] = find(mat == -1);
   
   neigh = zeros(size(mat, 1), 2);
   neigh(r1,1) = n1;
   neigh(r2,2) = n2;
end

% ----------------------------------------------------------------------------

function [tfaces, bfaces] = identify_lateral_faces(G, cells)

  % Etablishing total number of faces for each cell to consider
  start_ix = G.cells.facePos(cells);
  end_ix = G.cells.facePos(cells+1); % really one-past end
  fnum = end_ix - start_ix;
  
  rep = [fnum, max(fnum) - fnum]'; % second row represents padding

  % rows represent faces (padded, if necessary), columns cells
  f_mat = repmat(G.cells.faces(start_ix)', max(fnum), 1);
  
  % indices to non-padded elements
  f_mat_ix = reshape(rldecode(repmat([true;false], numel(cells), 1), ...
                              rep(:)), max(fnum), []);
  
  f_mat(f_mat_ix) = G.cells.faces(mcolon(start_ix, end_ix-1));
  
  % To pick the top face, we eliminate faces whose normals have a larger
  % lateral than vertical component, and pick the 'shallowest' of the
  % remaining faces
  
  % depth value for each face
  z_mat = reshape(G.faces.centroids(f_mat(:), end), max(fnum), []);
  
  nz_mat = reshape(G.faces.normals(f_mat(:), end) ./ ...
                   sqrt(sum(G.faces.normals(f_mat(:), 1:end-1).^2, 2)), ...
                   max(fnum), []);
  
  z_mat(nz_mat < 1/2) = inf; % these sides are more vertical than horizontal,
                             % disregard them by giving them infinite depth
  
  % picking the shallowest face
  [~, row] = min(z_mat);
  tfaces = f_mat(row(:) + (0:numel(cells)-1)' * max(fnum));

  % picking the deepest faces
  z_mat(nz_mat < 1/2) = -inf;
  [~, row] = max(z_mat);
  bfaces = f_mat(row(:) + (0:numel(cells)-1)' * max(fnum));
  
end

% ----------------------------------------------------------------------------

function [active_cols, col_cells, G] = identify_column_cells(G)
% 'active_cols': 'true' for all columns with at least one active cell.
% 'col_cells': Matrix containing the indices of active cells for each column.
%              Each column in this matrix correspond to a vertical stack of cells.
% 'G'        : The input grid G is modified if cells have to be removed.
   
   [lsize, lnum] = layer_size_and_number(G);
   
   col_cells = zeros(lsize, lnum);
   col_cells(G.cells.indexMap) = 1:G.cells.num;
   col_cells = col_cells'; % to make matrix columns correspond to grid columns   
   
   %% First, we peel away inactive cells _above_ the top surface
   
   % following line gives us the number of inactive cells at the top of each pillar
   inactive_upper = sum(cumsum(col_cells) == 0)'; 
   
   % remove the inactive top cells from matrix
   keeps = lnum(ones(lsize, 1)) - inactive_upper; % how many elements to keep per pillar
   c_ixs = rldecode((1:lsize)', keeps); % column indices
   r_ixs = mcolon(inactive_upper+1, lnum * ones(lsize, 1))';
   vals  = col_cells(sub2ind([lnum, lsize], r_ixs, c_ixs));
   r_ixs_new = mcolon(ones(size(inactive_upper)), lnum - inactive_upper)';
   col_cells = accumarray([r_ixs_new, c_ixs], vals, size(col_cells));
      
   %% Then, we peel away any cells _below_ remaining, inactive cells
   
   % for each row, true until first zero is encountered
   keeps = logical(cumprod(col_cells)); 
   discard_cells = reshape(col_cells(~keeps), [], 1);
   discard_cells = discard_cells(discard_cells ~= 0);
   col_cells(~keeps) = 0;
                          
   % Determine which pillars contain at least one cell
   active_cols = logical(sum(col_cells));
   
   %% Remove discarded cells from parent grid
   if ~isempty(discard_cells)
      G = computeGeometry(removeCells(G, discard_cells));
      % call function again, this time with a grid that does not need any
      % further removals.  (This is the easiest way to deal with the new indexing).
      [active_cols, col_cells, G] = identify_column_cells(G);
   end
end

% ----------------------------------------------------------------------------

function [lsize, lnumber] = layer_size_and_number(G)

% We assume the grid is either cartesian, or a layered grid with the fields
% 'numLayers' and 'layerSize'.  What is returned is the _logical_ layer
% size.  It does not take into account whether cells are active or not.
   
   if isfield(G, 'layerSize') 
      lsize   = G.layerSize;
      lnumber = G.numLayers;
   else
      lsize   = prod(G.cartDims(1:2));
      lnumber = G.cartDims(3);
   end
end

% ----------------------------------------------------------------------------

function seq = determine_facetag_sequence(G)
% determine the W/E/S/N sequence following the edges of the top surface of a
% cell.  Assuming all cells have the same layout, we determine this sequence
% looking only at the first grid cell.
   assert(isfield(G, 'cartDims')); % should only be attempted for cartesian grids!
   
   faces = G.cells.faces(G.cells.facePos(1):G.cells.facePos(2)-1, :);
   faces = sortrows(faces, 2); % now should be sorted according to logical direction
   faces = faces(1:5,1); % discard tags, and bottom face (which we do not need)
   
   nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
   assert(numel(nodes) == 20);
   nodes = reshape(nodes, 4, 5); % column 5 corresponds to the nodes of the top face

   ixfun = @(u, v) find(arrayfun(@(x) numel(intersect(nodes(:,x), [u, v]))==2, 1:4));
   
   seq = zeros(4,1);
   seq(1) = ixfun(nodes(1, end), nodes(2, end)); % cardinal dir. for first top edge
   seq(2) = ixfun(nodes(2, end), nodes(3, end)); % cardinal dir. for 2nd top edge
   seq(3) = ixfun(nodes(3, end), nodes(4, end)); % cardinal dir. for 3rd top edge
   seq(4) = ixfun(nodes(4, end), nodes(1, end)); % cardinal dir. for 4th top edge
   
end
