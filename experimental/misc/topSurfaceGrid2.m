function [Gt, G] = topSurfaceGrid2(G, varargin)

   opt = merge_options(struct('AddFaults', false), varargin{:});
   
   %% Identify columns in layered 3D grid
   [active_cols, col_cells] = identify_column_cells(G);
   
   %% Identify top faces
   tfaces = identify_top_faces(G, col_cells(1, active_cols));
   
   %% Construct the 2D grid
   Gt = top_faces_to_grid(G, tfaces);
   
   %% Identify fault nodes
   fault_nodes = identify_fault_nodes(G, Gt, opt.AddFaults);
   
   
   %% Stitch-up non-fault discontinuities
   Gt = stitch_surface_discontinuities(Gt, fault_nodes);
   
   %% Compute extra fields
   
end

% ----------------------------------------------------------------------------

function Gt = stitch_surface_discontinuities(Gt, fault_nodes)

   % identify nodes that are identical when z-coordinate is ignored, but
   % separate otherwise  
   [ncoords, ixs] = uniqueNodes(Gt.nodes.coords);
   
end

% ----------------------------------------------------------------------------

function fnodes = identify_fault_nodes(G, Gt, add_faults)

   fnodes = [];
   if ~add_faults
      return
   end
   error('unimplemented');
   
end

% ----------------------------------------------------------------------------

function Gt = top_faces_to_grid(G, tfaces)
   
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
   
   % Setting the 'faces' struct.
   Gt.faces = struct('num', size(uedges, 1), 'nodePos', [1:2:numel(uedges)+1], ...
                     'nodes', reshape(uedges', numel(uedges), []), ...
                     'neighbors', find_neighbors(faces, f_nodenum, orient), ...
                     'z', mean([Gt.nodes.z(uedges(:,1)), Gt.nodes.z(uedges(:,2))],2));
   
   %% Constructing cells of Gt (i.e. corresponding to top faces in the 3D grid)
   Gt.cells = struct('num', numel(tfaces), 'faces', faces, ...
                     'facePos', [1; last_node_ix+1], 'z', G.faces.centroids(tfaces, 3));
   
   %% Other immediate fields
   Gt.griddim = 2;
   Gt.type = [G.type, {'topSurfaceGrid'}];
end


% ----------------------------------------------------------------------------
function neigh = find_neighbors(cell_faces, num_cell_faces, orient)

   num_cells = numel(num_cell_faces);
   cellInx = rldecode([1:num_cells]', num_cell_faces);

   mat = accumarray([cell_faces, cellInx], 1, [], [], [], true);
   assert(all(sum(mat,2) <= 2)); % max 2 neighbors per face
   assert(all(sum(mat,2) > 0)); % no face without at least one neighbor
   
   cs = cumsum(mat, 2);
   [r1, n1] = find(mat & (cs==1));
   [r2, n2] = find(mat & (cs==2));
   
   neigh = zeros(size(mat, 1), 2);
   neigh(r1,1) = n1;
   neigh(r2,2) = n2;

   % determining orientation
   dir_mat = accumarray([cell_faces, cellInx], orient(:,1));
   flip = dir_mat(sub2ind(size(dir_mat), r1, n1)) == 2;
   neigh(flip,:) = fliplr(neigh(flip,:));
end

% % ----------------------------------------------------------------------------
% function neigh = find_neighbors(cell_faces, num_cell_faces, orient)

%    num_cells = numel(num_cell_faces);
%    cellInx = rldecode([1:num_cells]', num_cell_faces);

%    mat = accumarray([cell_faces, cellInx], 1, [], [], [], true);
%    assert(all(sum(mat,2) <= 2)); % max 2 neighbors per face
%    assert(all(sum(mat,2) > 0)); % no face without at least one neighbor
   
%    cs = cumsum(mat, 2);
%    [r1, n1] = find(mat & (cs==1));
%    [r2, n2] = find(mat & (cs==2));
   
%    neigh = zeros(size(mat, 1), 2);
%    neigh(r1,1) = n1;
%    neigh(r2,2) = n2;

%    % determining orientation
%    dir_mat = accumarray([cell_faces, cellInx], orient(:,1));
%    flip = dir_mat(sub2ind(size(dir_mat), r1, n1)) == 2;
%    neigh(flip,:) = fliplr(neigh(flip,:));
% end

% ----------------------------------------------------------------------------

function tfaces = identify_top_faces(G, cells)

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
  
  tfaces = f_mat(row(:) + [0:numel(cells)-1]' * max(fnum));
  
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
   c_ixs = rldecode([1:lsize]', keeps); % column indices
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