function [Gt, G] = topSurfaceGrid2(G, varargin)

   opt = merge_options(struct('AddFaults', false), varargin{:});
   
   %% Identify columns in layered 3D grid
   [active_cols, col_cells] = identify_column_cells(G);
   
   %% Identify top faces
   tcells = col_cells(1, active_cols); % top cells
   tfaces = identify_lateral_faces(G, tcells); % top faces
   
   %% Construct the 2D grid
   % 'orient' informs about edge orientation, and will be used later when
   % identifying cell neighbors.  However, this has to wait until unwanted
   % discontinuities have been removed (see a couple of lines further down)
   [Gt, orient] = construct_grid_from_top_faces(G, tcells, tfaces);
   
   %% Identify fault nodes
   fault_nodes = identify_fault_nodes(G, Gt, tfaces, opt.AddFaults);
   
   %% Stitch-up non-fault discontinuities
   % also identify cell neighbors.
   Gt = stitch_surface_discontinuities(Gt, fault_nodes, orient);

   %% Compute and add column information
   [Gt.columns, Gt.cells.columnPos, Gt.cells.H] = ...
       compute_column_info(G, col_cells(:, active_cols));
   
   %% Fill in the remaining fields
   Gt.parent          = G;
   Gt.grav_pressure   = @(G, omega) gravPressureVE_s(G, omega);
   Gt.primitives      = @primitivesMimeticVE_s;
end

% ----------------------------------------------------------------------------

function [cols, col_pos, H] = compute_column_info(G, col_cells)

   col_sizes = sum(col_cells~=0)';
   col_pos = cumsum([1; col_sizes]);

   start_ix = cumsum([1;repmat(size(col_cells, 1), size(col_cells, 2), 1)]);
   start_ix = start_ix(1:end-1);
   end_ix = start_ix + col_sizes - 1;
   
   cols.cells = col_cells(mcolon(start_ix, end_ix))';
   
   [top_faces, bot_faces] = identify_lateral_faces(G, cols.cells);

   cols.z = G.faces.centroids(bot_faces, 3);   
   cols.dz = cols.z - G.faces.centroids(top_faces, 3);
   H = accumarray(rldecode((1:length(col_sizes))', col_sizes'), cols.dz);
   
end

% ----------------------------------------------------------------------------

function Gt = stitch_surface_discontinuities(Gt, fault_nodes, orient)

   % identify nodes that are identical when z-coordinate is ignored, but
   % separate otherwise.
   coords = [Gt.nodes.coords, zeros(length(Gt.nodes.coords), 1)];
   coords(fault_nodes, end) = 1:length(fault_nodes); % prevents fault nodes from merging
   [ncoords, nc_ix] = uniqueNodes(coords);
   
   % Updating Gt
   Gt.nodes.num = size(ncoords, 1);
   Gt.nodes.coords = ncoords(:, 1:2);
   Gt.nodes.z = accumarray(nc_ix, Gt.nodes.z) ./ accumarray(nc_ix, 1);
   
   Gt.faces.nodes = nc_ix(Gt.faces.nodes);
   
   % With topology now in place, we can finally compute cell neighbors
   Gt.faces.neighbors = ...
       find_neighbors(Gt.cells.faces, diff(Gt.cells.facePos), orient);
end

% ----------------------------------------------------------------------------

function fnodes = identify_fault_nodes(G, Gt, top_faces, add_faults)

   fnodes = [];
   if ~add_faults || ~isfield(G.faces, 'tag')
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
   edges = [fnodes_ixs, [0;fnodes_ixs(1:end-1)]];
   last_node_ix = cumsum(f_nodenum); % last node per face
   first_node_ix = [1; last_node_ix(1:end-1)+1]; % first node per face
      
   edges(first_node_ix, 2) = edges(last_node_ix, 1);
   
   % orienting each edge according to node index, but keeping track of original
   % (correct) orientation
   [edges, orient] = sort(edges, 2);
   
   % Removing duplicates from edge set
   [uedges, faces] = uniqueNodes(edges);
   
   % Setting the 'faces' struct.  The 'neighbor' field will have to wait
   % until unwanted discontinuities have been removed from the grid.
   Gt.faces = struct('num', size(uedges, 1), 'nodePos', (1:2:numel(uedges)+1), ...
                     'nodes', reshape(uedges', numel(uedges), []), ...
                     'z', mean([Gt.nodes.z(uedges(:,1)), Gt.nodes.z(uedges(:,2))],2));
   
   %% Constructing cells of Gt (i.e. corresponding to top faces in the 3D grid)
   Gt.cells = struct('num', numel(tfaces), 'faces', faces, ...
                     'facePos', [1; last_node_ix+1], 'z', G.faces.centroids(tfaces, 3));
   
   %% Other immediate fields
   Gt.griddim = 2;
   Gt.type = [G.type, {'topSurfaceGrid'}];
   Gt.cells.map3DFace = tfaces;
   
   if isfield(G, 'cartDims')
      Gt.cartDims = G.cartDims(1:end-1);
      [I, J] = ind2sub(G.cartDims, G.cells.indexMap(tcells));
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

function [active_cols, col_cells] = identify_column_cells(G)
% 'active_cols': 'true' for all columns with at least one active cell.
% 'col_cells': Matrix containing the indices of active cells for each column.
%              Each column in this matrix correspond to a vertical stack of cells.
   
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
   col_cells(~keeps) = 0;
                          
   % Determine which pillars contain at least one cell
   active_cols = logical(sum(col_cells));
   
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