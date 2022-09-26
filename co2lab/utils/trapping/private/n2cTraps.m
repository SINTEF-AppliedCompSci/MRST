function [ctraps, ctrap_zvals, ctrap_regions, csommets, ctrap_connectivity, crivers] = ...
	 n2cTraps(Gt, ntrap_regions, ntrap_zvals, ntrap_dstr_neigh, ntrap_connectivity, ntraps, erivers)
% Function converting traps and spill field information from a node-based 
% representation to a cell-based one.
%
% SYNOPSIS
%   [ctraps, ctrap_zvals, ctrap_regions, ctrap_connectivity] = 
%      n2cTraps(Gt, ntrap_regions, ntrap_zvals, ntrap_connectivity)
%
% PARAMETERS:
%   Gt                 - 2D grid structure 
%   ntrap_regions      \ 
%   ntrap_dstr_neigh   | node-based trap information, as computed by the
%   ntrap_zvals        | 'compute-trap-regions' function                
%   ntrap_connectivity / 
% 
% RETURNS:
%   The cell-based versions of the trap-information data matrices provided
%   as input parameters to this function:
%     ctraps         - One value per grid cell; zero for nontrap-cells,
%                      trap number from 1 upwards) for trap cells.
%     ctrap_zvals    - vector with one element per trap, giving the
%                      z-spill value for that trap
%     ctrap_regions  - one value per grid cell.  Gives the trap number of
%                      the trap that the node spills into, or zero if the
%                      cell belongs to the spill region of the 'exterior'
%                      of the domain.
%     csommets       - indices to all cells that represent local maxima
%                      (NB: these are all trap cells, but there may be more
%                      than one sommet per trap)
%     ctrap_connectivity - (Sparse) adjacency matrix with one row/column
%                      per trap. Row 'i' is nonzero only for columns 'j'
%                      where trap 'i' spills directly into trap 'j' when
%                      overflowing.
%     lost_regions  - indices of those edge-based regions whose traps did
%                      not get projected to any cell (e.g. because they
%                      were too small)
%     crivers        - one cell array per trap, containing the 'rivers'
%                      exiting that trap.  A river is presented as a
%                      sequence of consecutive grid cells that lie
%                      geographically along the river.  A river starts in
%                      the trap, and ends either in another trap or at the
%                      boundary of the domain.
%
% SEE ALSO:
%   `computeNodeTraps`

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


% As long as traps are not eliminated, connectivities, z-values and spill
% regions are identical to the node-based version...
  ctrap_connectivity = ntrap_connectivity;
  ctrap_zvals = ntrap_zvals;
  initial_num_regions = numel(ctrap_zvals);
  if isfield(Gt.cells, 'z')
     cell_centroids = Gt.cells.z;
  else
      cell_centroids = compute_cell_z_vals(Gt);
  end

  % identify sommets 
  node_sommets = zeros(size(ntrap_regions));
  node_sommets(intersect(find(ntrap_dstr_neigh == 0), find(ntrap_regions ~= 0))) = 1;
  
  % projecting spill regions and sommets from nodes to cells
  ctrap_regions = nodeFieldToCellField(Gt, ntrap_regions);
  cell_sommets  = n2c_sommets(Gt, node_sommets);

  % isolating traps
  ctraps = zeros(size(ctrap_regions));
  for i = 1:initial_num_regions
    ctrap_region_indices = find(ctrap_regions == i);
    ctrap_cells = intersect(find(cell_centroids < ctrap_zvals(i)), ctrap_region_indices);
    ctraps(ctrap_cells) = i;
  end

  % if any region was eliminated (did not end up to be projected to a single cell), 
  % formally remove it and update the adjacency matrix accordingly.
  lost_regions = identify_lost_regions(initial_num_regions, ctraps);
  
  % regions that are associated to lost traps are re-associated to remaining
  % traps below
  for i = lost_regions
      ixs = (ctrap_regions == i);
      ctrap_regions(ixs) = next_remaining_downstream_trap(i, ...
                                                        ctrap_connectivity, ...
                                                        lost_regions);
  end
  
  % Removing indexes to lost traps, and re-indexing remaining traps
  ctrap_regions = re_index(ctrap_regions, lost_regions);
  ctraps        = re_index(ctraps,        lost_regions);
  ctrap_zvals(lost_regions) = [];  

  % keeping cell sommets that still correspond to traps, and sorting them
  cell_sommets = cell_sommets .* ctrap_regions; % values change from 0/1 to trap index
  %cell_sommets_ixs = unique_cell_sommet_ixs(cell_sommets, cell_centroids);
  cell_sommet_ixs = find(cell_sommets);
  [~ , order] = sort(cell_sommets(cell_sommet_ixs));
  csommets = cell_sommet_ixs(order);
  
  % Updating connectivity matrix and river structure by removing lost regions one by one
  for i = lost_regions 
    % NB: for this to work, it's important that 'lost_regions' is sorted in 
    % descending order
    tnodes = find(ntraps == i);
    
    % Update connectivity matrix
    rmerge = sparse(eye(size(ctrap_connectivity)));  

    lmerge = sparse(eye(size(ctrap_connectivity)));
    lmerge(:,i) = ctrap_connectivity(:,i); %#ok

    lmerge(i,:) = []; %#ok
    rmerge(:,i) = []; %#ok
    
    ctrap_connectivity = spones(lmerge * ctrap_connectivity * rmerge);

    % Update rivers
    upstream_regs = find(ntrap_connectivity(:,i));
    if ~isempty(upstream_regs)
        for j = upstream_regs'
            for k = 1:numel(erivers{j})
                if ntrap_regions(erivers{j}{k}(end)) == i
                    for l = 1:numel(erivers{i})
                        %fprintf('%d %d %d %d', j, k, i, l);
                        %fprintf('(%d, %d) - (%d, %d)\n', size(erivers{j}{k}), size(erivers{i}{l}))
                        %erivers{j}{k} = [erivers{j}{k}; erivers{i}{l}];
                        erivers{j}{k} = join_rivers(erivers{j}{k}, erivers{i}{l}, Gt, tnodes);
                    end
                end
            end
        end
    end
  end
  % removing rivers corresponding to lost regions (information of these
  % rivers has already been taken care of by merging taking place in the
  % previous loop.
  erivers(lost_regions) = [];
  crivers = project_rivers_to_cells(Gt, erivers);
  crivers = truncate_rivers_inside_traps(crivers, ctraps);
end
%===============================================================================

function river = join_rivers(r1, r2, Gt, fillnodes)
    if r1(end) == r2(1)
        % we know these are connected.  Join directly.
        river = [r1(1:end-1); r2];
    else
        % processing endpoints to make a clean search
        eix = find(ismember(r1, fillnodes), 1, 'last');
        six = find(ismember(r2, fillnodes), 1, 'first');
        r1 = r1(1:eix);
        r2 = r2(six:end);
        
        assert(ismember(r1(end), fillnodes));
        assert(ismember(r2(1), fillnodes));
        
        start  = find(fillnodes == r1(end));
        target = find(fillnodes == r2(1));

        M = subregionConnectivityMatrix(Gt, fillnodes, true);
        path = shortestPath(M, start, target);
        path = fillnodes(path);
        river = [r1(1:end-1); path(1:end-1); r2];
    end        
end

% ============================================================================
function cell_sommets = n2c_sommets(Gt, node_sommets)
% For each sommet, associate a corresponding grid cell

% Determining a 'representative' node (the shallowest one) for each cell
    cnodes = activeCellNodes(Gt);
    num_cells = size(cnodes,2);
    [min_z, representative] = min(Gt.nodes.z(cnodes)); %#ok<ASGLU>
    cnodes_rep = cnodes(sub2ind(size(cnodes), representative, 1:num_cells));
    
    % giving an unique value to each sommet, making it possible to recognize
    % it from others
    ns_ix = find(node_sommets);
    num_sommets = numel(ns_ix);
    node_sommets(ns_ix) = 1:num_sommets;
    
    % determine which cells have a given node sommet as its representative
    % (usually more than one!)
    cell_sommets = node_sommets(cnodes_rep);
    
    % choosing unique cell to represent each node
    cwsn = find(cell_sommets);  % indices of all cells with a sommet node
    M = [cwsn, cell_sommets(cwsn), Gt.cells.z(cwsn)];
    M = sortrows(M, [2, 3]); % sort wrt. sommet ix, then wrt. cell depth
    [~, k] = unique(M(:,2), 'first');
    cell_sommets(:) = 0;
    cell_sommets(M(k,1)) = 1;
end

% function csixs = unique_cell_sommet_ixs(cell_sommets, cell_centroids)
% % For each sommet, determine one unique cell to represent it.  This should be
% % the 'highest' of the cells it is projected to, i.e. the one whose centroid
% % has the lowest depth-value.
%     ixs = find(cell_sommets); % not necessarily unique
%     candidate_centroids = cell_centroids(ixs);
%     [~, I] = sort(candidate_centroids); % lowest z-values go first
%     ixs = ixs(I); % ixs now arranged so that cells with lowest x-values
%                   % referenced first
%     [~, keep_ix] = unique(ixs, 'first');
%     csix = ixs(keep_ix);
% end

%===============================================================================
function lost_regions = identify_lost_regions(num_regions, mapped_regions)
  % Determine which regions were 'eliminated' when passing from node-based to 
  % cell-based representation (i.e., not a single cell was attributed to it.

  % 'fliplr' ensures reverse ordering, so that the vector starts with the
  % highest index.
  lost_regions = fliplr(setdiff(1:num_regions, unique(mapped_regions)));
end

%===============================================================================
function mat = re_index(mat, indices)
  for i = indices
      ixs = mat > i;
      mat(ixs) = mat(ixs) - 1;
  end
end	 

%===============================================================================
function ctoids = compute_cell_z_vals(Gt)
  % Compute z-value of each cell as the average of the z-values of its four
  % nodes.
  %
  % We presuppose here that each cell has exactly four nodes.  If this is not 
  % a safe assumption, a more involved implementation is needed below.
  ctoids = mean(Gt.nodes.z(activeCellNodes(Gt)));

end


% %===============================================================================
% function cell_rivers = project_rivers_to_cells(Gt, edge_rivers)
    
%     cell_rivers = cell(size(edge_rivers, 1), 1);
%     cellnodes = activeCellNodes(Gt); % (4 x m)-sized matrix; m = # of active cells
%                                 % in Gt.  Each col. holds the indices of the
%                                 % 4 corner nodes of the corresponding cell.
%     for trap_ix = 1:size(edge_rivers, 1)
%         for r_ix = 1:numel(edge_rivers{trap_ix})
%             nodes_ix = edge_rivers{trap_ix}{r_ix};
%             cells_ix = [];
%             for i = 1:4
%                 % finding index of cells having these nodes as corner 'i'.
%                 cells_ix = [cells_ix, find(ismember(cellnodes(i, :), nodes_ix))]; %#ok

%                 % remove nodes that have been mapped, and iterate on the next corner
%                 % nodes_ix = nodes_ix(find(~ismember(nodes_ix, cellnodes(i, cells_ix))));
%                 if (isempty(nodes_ix)); break; end;
%             end
%             cell_rivers{trap_ix}{r_ix} = cells_ix;
%         end
%     end
% end



%===============================================================================
function cell_rivers = project_rivers_to_cells(Gt, edge_rivers)
    
    nnum = Gt.nodes.num;
    
    % Assume exactly two nodes per edge
    assert(unique(diff(Gt.faces.nodePos)) == 2); 
                                                 
    % one row per edge, giving its end node indices (in ascending order)
    enodes_sorted = sort(reshape(Gt.faces.nodes, 2, []), 1, 'ascend')';
    enodes_lookup = sparse(enodes_sorted(:,1), enodes_sorted(:,2), ...
                           1:Gt.faces.num, nnum, nnum);
    
    % make tables with the diagonals of each cell (needed in case the river
    % goes diagonally across a cell)
    cellnodes = activeCellNodes(Gt)';
    diag1 = sort([cellnodes(:,1), cellnodes(:,4)], 2, 'ascend');
    diag2 = sort([cellnodes(:,2), cellnodes(:,3)], 2, 'ascend');
    diag_lookup = sparse(diag1(:,1), diag1(:,2), (1:Gt.cells.num)', nnum, nnum) + ...
                  sparse(diag2(:,1), diag2(:,2), (1:Gt.cells.num)', nnum, nnum);
    
    cell_rivers = cell(size(edge_rivers, 1), 1);
    if numel(edge_rivers) == 0
       return; % nothing to do
    end
    
    for trap_ix = 1:size(edge_rivers, 1)
        for r_ix = 1:numel(edge_rivers{trap_ix})
            nodes_ix = edge_rivers{trap_ix}{r_ix};
            
            % removing consecutive duplicates 
            nodes_ix = nodes_ix([1; diff(nodes_ix)]~=0);

            cells_ix = [];
            num_edges = numel(nodes_ix) - 1;
            assert(num_edges > 0); % a river should have at least one edge...
            for i = 1:num_edges
                enodes  = sort([nodes_ix(i) nodes_ix(i+1)], 2, 'ascend');
                edge_ix = enodes_lookup(enodes(1), enodes(2));

                if full(edge_ix) == 0
                    % No such edge.  The river is here going diagonally
                    % across a cell.  Determine this cell, and add it as a
                    % river cell
                    c_ix = full(diag_lookup(enodes(1), enodes(2)));
                    assert(c_ix ~= 0); % should be one of the two diagonals
                    cells_ix = [cells_ix, c_ix]; %#ok<AGROW>
                else
                    % Determining cells bordering this edge.   If there is only
                    % one cell, 'project' the edge to it; if there are two,
                    % 'project the edge to the shallowest of them
                    neigh_cells = sort(Gt.faces.neighbors(edge_ix,:), 'ascend');
                    if neigh_cells(1) == 0
                        cells_ix = [cells_ix, neigh_cells(2)]; %#ok<AGROW>
                    else
                        [dummy, min_ix] = min(Gt.cells.z(neigh_cells)); %#ok<ASGLU>
                        cells_ix = [cells_ix, neigh_cells(min_ix)]; %#ok<AGROW>
                    end
                end
            end
            cell_rivers{trap_ix}{r_ix} = cells_ix;
        end
    end
end

%===============================================================================

function crivers =  truncate_rivers_inside_traps(crivers, ctraps)
    
    for t_ix = 1:numel(crivers)
        for r_ix = 1:numel(crivers{t_ix})
            river = crivers{t_ix}{r_ix};
            start_ix = max(1, find(ctraps(river)==0,1) - 1);
            end_ix   = min(numel(river), find(ctraps(river)==0,1, 'last') +1);
            
            crivers{t_ix}{r_ix} = river(start_ix:end_ix);
        end
    end
end

%===============================================================================
function ixs = next_remaining_downstream_trap(trap_ixs, connectivity, lost_regions)
% Return next remaining downstream trap, or '0' if there is none.

    downstream_traps = find(sum(connectivity(trap_ixs, :),1));
    
    % If there are no downstream traps at all, it measn that the flow will
    % exit the domain.  Return '0' to indicate this.
    if isempty(downstream_traps)
        ixs = 0;
        return;
    end
    
    % Are any of the downstream traps remaining in the system after lost
    % regions have been eliminated?
    rem_traps = setdiff(downstream_traps, lost_regions);
    
    if ~isempty(rem_traps)
        % 'rem_traps' contain one or more downstream traps that will remain
        % after those in 'lost_regions' have been removed.  One is not more
        % correct than another, so we arbitrarily pick the first one.
        ixs = rem_traps(1);
    else
        ixs = ...
            next_remaining_downstream_trap(intersect(downstream_traps, lost_regions), ...
                                           connectivity, ... 
                                           lost_regions);
    end
end
