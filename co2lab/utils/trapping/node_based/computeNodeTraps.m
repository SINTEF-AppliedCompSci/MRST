function res = computeNodeTraps(Gt, closed_bnodes, closed_fnodes)
%
% Compute tracks and spill regions for the 2D surface grid 'Gt'.  The
% computation is node-based, not cell based.  Cell-based results can be
% obtained by subsequently applying the function n2cTraps ("node-to-cell
% traps"). 
% 
% SYNOPSIS
%    res = computeNodeTraps(Gt);
% 
% PARAMETERS:
%    Gt            - 2D-grid with z-values assigned to the nodes (Gt.nodes.z
%                    should exist). 
%    closed_bnodes - vector with indices of nodes on the closed part of the
%                    boundary 
%    closed_fnodes - vector with indices of nodes along closed fault lines
% 
% RETURNS: 
%    res - structure with the following fields:
%          traps           - one value per grid node. Zero for non-trap 
%                            nodes, trap number (from 1 upwards) for
%                            trap nodes. 
%          trap_zvals      - vector with one element per trap, giving
%                            the z-spill value for that trap
%          trap_sp_edge_ix - index(ices) of the edge(s) at the spill
%                            point(s) of each trap. There is typically
%                            only one spill point, but might be more in
%                            degenerate cases.  Therefore, the data is
%                            returned as a cell array rather than a
%                            vector or matrix. 
%          dstr_neigh      - one value per grid node.  Gives index of
%                            'most upstream' neighbor node, or '0' if
%                            there is no upstream neighbor (i.e. the
%                            node in question is a sommet or on the
%                            boundary. 
%          trap_regions    - one value per grid node.  Gives the number
%                            of the trap that the node spills into, or
%                            zero if the node spills out of the surface
%                            domain. 
%          connectivity    - (Sparse) adjacency matrix with one
%                            row/column per trap.  Row 'i' is nonzero
%                            only for columns 'j' where trap 'i' leads
%                            directly into trap 'j'.
%          rivers          - One cell array per trap, containing the
%                            'rivers' exiting that trap.  A river is
%                            presented as a sequence of consecutive grid
%                            nodes that lie geographically on the river.
%                            A river starts in the trap, and ends either
%                            in another trap or at the boundary of the
%                            domain. 
% SEE ALSO:
% `n2cTraps`

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

    % Computing spill field
    [res.dstr_neigh, regions, spill_edges] = ...
        nodeSpillField(Gt, closed_bnodes, closed_fnodes);
    num_regions = max(regions); % There is one region per local minimum.  
                                % Will be properly aggregated into traps
                                % below. 

    % Computing adjacency matrix between all (unclustered) regions
    [immediate_adjacencies, all_z_vals, spoint_edges] = ...
        assemble_connectivity_matrix(eye(num_regions), spill_edges);

    % Clustering spill regions for local minima into spill regions for traps 
    % (each trap might contain several local minima)
    [reg_groups, res.trap_zvals, res.trap_sp_edge_ix, res.connectivity, ~] = ...
       track_stream(immediate_adjacencies, spill_edges, all_z_vals, spoint_edges);

    num_traps = size(reg_groups, 2);
    res.trap_regions = zeros(Gt.nodes.num, 1);
    res.traps = res.trap_regions;
    
    % Mapping trap regions to nodes, and determining what parts of the
    % trap spill regions that actually belong to the traps.
    for i = 1:num_traps % only loop over interior regions, so we start a
                        % 1, not at 0.  (trap_regions initialized to 0
                        % anyway, so nodes belonging to the exterior
                        % spill region already have the correct value of
                        % 0). 

	trap_region_indices = find(ismember(regions, find(reg_groups(:,i))));
        res.trap_regions(trap_region_indices) = i;

	trap_nodes = intersect(find(Gt.nodes.z(:) <= res.trap_zvals(i)), ...
                               trap_region_indices); 
	res.traps(trap_nodes) = i;
    end
    
    res.rivers = trace_rivers(Gt, res.traps, res.trap_sp_edge_ix, ...
                              res.dstr_neigh); 

end

function [reg_groupings, grouping_zvals, grouping_sp_ix, cmat, icmat] = ...
        track_stream(cmat, spill_edges, zvals, spoint_edges)
% Based on an initial partition of a surface into spill regions based on
% local surface minima, compute traps, trap spill regions and
% trap-connections by coalescing into clusters. 
% PARAMETERS:
%   cmat         - directed adjacency matrix (square, sparse) between
%                  local minima-based spill regions
%   spill_edges  - complete list of spill edges and associated
%                  information.   The information is presented
%                  as a matrix with one line per spill edge and seven 
%                  columns, giving the following information about the
%                  spill edge: 
%                    Col 1: index of spill edge (into the underlying
%                           grid structure) 
%                    Col 2: region of edge's first node
%                    Col 3: region of edge's second node
%                    Col 4: z-value of edge's first node
%                    Col 5: z-value of edge's second node
%                    Col 6: index of edge's first node
%                    Col 7: index of edge's second node
%                  This matrix can be produced by the 'nodeSpillField'
%                  function.  For the purposes of the present function
%                  ('track_stream'), it is used to recalculate
%                  conenctions when regions merge into clusters.                    
%   zvals        - z-value of spill point(s) for each local-minima
%                  region  
%   spoint_edge  - index(ices) of the edge(s) at the spill point of each 
%                  grouping.  Since there might be a variable number of
%                  spill points for a grouping (typically only one, but
%                  might be more in degenerate cases), the information
%                  is provided as a cell array.
% 
% RETURNS:
% reg_groupings  - sparse matrix expressing which nodes of 'cmat' have
%                  been grouped together.  (This happens when 'cycles'
%                  are detected, so that no of the involved nodes can be
%                  said to be upstream from the others.)  The matrix
%                  has one row per node and one column per grouping. 
% icmat          - "iterated" cmat, i.e. cmat + cmat^2 + cmat^3 + ... until cmat^n = 0.
%                  This matrix thus expresses, for a given grouping
%                  (line), all the groupings that are 'upstream' of it
%                  (columns). 
% cmat           - adjacency (connectivity) matrix between traps, and
%                  between traps and the exterior.
% grouping_zvals - z-value of the spill point of each grouping (one
%                  value per grouping). (Can be thought of as an
%                  'updated' version of the argument zvals.)
% grouping_sp_ix - index(ices) of the spill point edge(s) for each
%                  grouping.  (Can be thought of as an 'updated' version
%                  of the argument 'spoint_edge'. 


    num_nodes = size(cmat, 1); % square, so could be size(cmat,2) as well.

    reg_groupings   = sparse(eye(num_nodes)); % to start off, each region constitutes a group
    icmat           = sparse(zeros(num_nodes));  % will accumulate to cmat | cmat^2 | cmat^3 + ...
    tmp             = sparse(cmat);               % temp. matrix to incrementally compute powers
                                                  % of cmat, and convert nonzero elements to 1.
    grouping_zvals  = zvals;
    grouping_sp_ix  = spoint_edges;
    
    % keep looping until cmat becomes the zero matrix (this should eventually happen)
    while any(any(tmp)) % true if tmp has nonzero elements
	
        icmat = icmat | tmp;

        % first, remove any loops that might have become apparent at this iteration.
        % Loops are revealed by nonzero diagonal elements in 'icmat'.
        % (Loops really represent regions that get clustered together as a result of spilling
        % into each others before spilling out of the total area).
        while (any(diag(icmat)))
            dnz = find(diag(icmat)); % indices of nonzero diagonal elements
            i = dnz(1);
            % nodes participating in the cycle are those both upstream AND upstream
            % of node i.  In other words, the intersection of nonzeros on the i-th row
            % and i-th column of 'tmp'.
            cycle_nodes = find(icmat(i,:) & icmat(:,i)');
            assert(numel(cycle_nodes) > 1); % a cycle should always contain >1 node.
            
            % All the cycle nodes should be merged to one grouping, to be represented
            % by the line/column of the node with the lowest index in the cycle.  Other
            % lines/columns should be removed
            new_node_ix = cycle_nodes(1);
            reg_groupings = merge_cols     (reg_groupings, cycle_nodes);
            cmat          = merge_rows_cols(cmat,  cycle_nodes);  cmat(new_node_ix, new_node_ix) = 0;
            icmat         = merge_rows_cols(icmat, cycle_nodes); icmat(new_node_ix, new_node_ix) = 0;
            tmp           = merge_rows_cols(tmp,   cycle_nodes);   

	        % determine new upstream connection(s) from this cluster 
            [neigh_clusters_ix, new_zval, sp_ixs] = ...
                locate_neighbor_traps(new_node_ix, reg_groupings, spill_edges);
            cmat(new_node_ix, neigh_clusters_ix) = 1;

            grouping_zvals(new_node_ix) = new_zval;
            grouping_zvals(cycle_nodes(2:end)) = [];
            
            grouping_sp_ix{new_node_ix} = sp_ixs;
            grouping_sp_ix(cycle_nodes(2:end)) = [];
        end
        tmp = 1.0 * (tmp * cmat ~= 0);
    end
end

function res = merge_cols(mat, cols_to_merge)
    % Merge a set of columns of booleans as if using 'or', keep a merged column, and remove all the other
    % involved columns.
    mat(:,cols_to_merge(1)) = any(mat(:, cols_to_merge), 2); % Merging columns into one
    mat(:, cols_to_merge(2:end)) = []; % removing all involved columns but the merged one
    res = sparse(mat);
end

function res = merge_rows_cols(mat, cols_rows_to_merge)
    % Merge a set of columns, and corresponing rows, of booleans as if
    % using 'or', keep a merged row and column, and remove all the other
    % involved columns/rows. 

    res = merge_cols(mat, cols_rows_to_merge); % merging columns
    res = merge_cols(res', cols_rows_to_merge); % merging rows
    res = res';
end


function [res_mat, z_heights, spillpoint_edge] = ...
        assemble_connectivity_matrix(reg_groupings, spill_edges)
% 
% SYNOPSIS:
% Compute the adjacencies between traps (i.e. which other trap(s) are
% spilled into when a given trap overflows.  The adjacencies are
% expressed as an adjacency matrix ('res_mat'). Each trap is represented
% as an group of local minima (of the top surface) that are  clustered
% together.  The grouping is expressed by a matrix ('reg_groupings')
% where the regions associated with each local minima are represented by
% the lines, and the clusters that constitutes the traps are represented
% by the columns. In order to determine the connections (which trap
% region spills into another), all spill edges are provided (c.f. the
% 'nodeSpillField' function for more information). 
%  
% PARAMETERS:
% reg_groupings - sparse matrix providing the information on which
%                 regions are grouped together to clusters.  Columns
%                 represent clusters, rows represent regions. Nonzero
%                 elements indicate which regions are included in each
%                 cluster.  Each region should belong to exactly one
%                 cluster.  If no regions are clustered together, this
%                 matrix should thus be the identity.
% spill_edges   - Information all spill edges resulting from the
%                 partitioning of the 2D surface grid.  These are all
%                 edges in the grid whose two nodes belong to two
%                 separate spill regions.  The information is presented
%                 as a matrix with one line per spill edge and seven
%                 column, giving the following information about the
%                 spill edge: 
%                    Col 1: index of spill edge (into the underlying
%                           grid structure) 
%                    Col 2: region of edge's first node
%                    Col 3: region of edge's second node
%                    Col 4: z-value of edge's first node
%                    Col 5: z-value of edge's second node
%                    Col 6: index of edge's first node
%                    Col 7: index of edge's second node
%                 This matrix can be produced by the 'nodeSpillField'
%                 function.                     
% 
% RETURNS:
% res-mat         - adjacency matrix 
% z_heights       - heights of the spill points leading to the
%                   adjacencies between clusters (traps) (one per cluster)
% spillpoint_edge - for each cluster, gives the index(ices) of the edges 
%                   located at the spill point(s).  (Typically only one,
%                   but can be multiple in 'degenerate' cases).  Since
%                   there might be more than one, this structure is
%                   presented as a cell array rather than as a plain
%                   vector. 
% 
% SEE ALSO:
% `nodeSpillField`, `computeNodeTraps`

    num_clusters = size(reg_groupings, 2);
    res_mat = sparse(num_clusters, num_clusters);
    z_heights = zeros(num_clusters, 1); 
    spillpoint_edge = cell(num_clusters,1);

    for cl = 1:num_clusters
	
	[neigh_clusters, highest_spill_z_value, spillpoint_edge{cl}] = ...
            locate_neighbor_traps(cl, reg_groupings, spill_edges);

        res_mat(cl, neigh_clusters) = 1;
        z_heights(cl) = highest_spill_z_value;
    end
end

%===============================================================================
function rivers = trace_rivers(Gt, traps, trap_spoint_edges, upstream_neighbor)
    num_traps = size(trap_spoint_edges,1);
    
    % Assumption: each face (edge) has exactly two nodes, which are stored
    %             consecutively in the Gt.faces.nodes vector.
    edge_nodes = reshape(Gt.faces.nodes, 2, []);
    
    rivers = cell(num_traps, 1);
    
    % computing the rivers exiting each trap
    for trap_ix = 1:num_traps
        % looping over rivers exiting this trap (typically 1, but can be more)
        for r_ix = 1:numel(trap_spoint_edges{trap_ix}) 
            % each river will at least contain the two nodes of the
            % spillpoint edge
            river_nodes = edge_nodes(:, trap_spoint_edges{trap_ix}(r_ix));

            % growing river at 'end' until a sommet is reached
            first_node = river_nodes(1);
            last_node  = river_nodes(2);
            % while ( upstream_neighbor(last_node) ~= 0 ) && ...
            %       ( traps(last_node) == 0)
            while (upstream_neighbor(last_node) ~= 0)
                last_node = upstream_neighbor(last_node);
                river_nodes = [river_nodes; last_node];                    %#ok
            end
            
            % growing river at 'start' until a sommet is reached
            % while ( upstream_neighbor(first_node) ~= 0) && ...
            %       ( traps(first_node) == 0)
            while (upstream_neighbor(first_node) ~= 0)
                first_node = upstream_neighbor(first_node);
                river_nodes = [first_node; river_nodes];                   %#ok
            end
            % ensuring that the river is oriented so that it starts in 'this'
            % trap
            assert((traps(river_nodes(1)) == trap_ix) || ...
                   (traps(river_nodes(end)) == trap_ix));
            if (traps(river_nodes(1)) ~= trap_ix)
                river_nodes = flipud(river_nodes);
            end
            rivers{trap_ix}{r_ix} = river_nodes;
        end
    end
end

function [neigh_traps, highest_spill_z_value, highest_spill_edge_ixs] ...
        = locate_neighbor_traps(trap_ix, all_traps, spill_edges)
% [This function is written for, and called by, the functions 'assemble_connectivity_matrix' 
%  and 'computeNodeTraps']
% For a given trap (as indicated by 'trap_ix'), determine which other traps
% it spills into (i.e. immediate "upstream" neighbor(s)). Each trap is 
% defined by a number of "local-minima" spill regions, as expressed by the 
% 'all_traps' matrix.
% 
% SYNOPSIS:
% [neigh_traps, highest_spill_z_value] = ...
%   locate_neighbor_traps(trap_ix, all_traps, spill_edges)
%
% PARAMETERS:
% trap_ix      - Index of the trap for which we want to identify the immediate upstream traps
% all_traps    - Matrix expressing how individual 'local-minima' spill regions have been grouped 
%                together as 'traps'.  Each line represent a local spill region, each column 
%                represents a trap.  Elements are either '1' (region 'i' is part of trap 'j') or 
%                '0' (region is not part of the trap).  Each region should belong to exactly 1 trap,
%                and each trap should contain one or more regions.
% spill_edges  - Information all spill edges resulting from the partitioning
%                of the 2D surface grid.  These are all edges in the grid
%                whose two nodes belong to two separate spill regions.  The
%                information is presented as a matrix with one line per
%                spill edge and seven columns, giving the following
%                information about the spill edge:
%                   Col 1: index of spill edge (into the underlying grid structure)
%                   Col 2: region of edge's first node
%                   Col 3: region of edge's second node
%                   Col 4: z-value of edge's first node
%                   Col 5: z-value of edge's second node
%                   Col 6: index of edge's first node
%                   Col 7: index of edge's second node
%                This matrix can be produced by the 'nodeSpillField' function.                    
%
% RETURNS:
% neigh_traps            - Index (or vector of indices) to the trap(s) immediate upstream of the 
%                          trap in question.  Can be empty if there is no trap upstream (i.e. 
%                          trap spills out of domain).
% highest_spill_z_value  - z-value of the trap's spill point.
% highest_spill_edge_ixs - Index(ices) of the corresponding spill edge(s)
%                          (typically one, but can theoretically be more)
%
% SEE ALSO:
% `computeNodeTraps`, `assemble_connectivity_matrix`, `nodeSpillField` 

    num_regions = size(all_traps, 1);
    member_regions = find(all_traps(:, trap_ix));
    [all_spills, local_ixs] = trap_spill_edges(spill_edges, member_regions, true);
    
    % highest spill spill is usually one single edge, but can be multiple edges
    highest_spill_z_value = all_spills(1,2);
    highest_spill_ixs = local_ixs(all_spills(:,2) == highest_spill_z_value,1);
    
    % remove cases where multiple spill edges leave from the _same_ node in
    % the region (only keep the steepest) 
    highest_spill_ixs = ...
        purge_multiple_spills_same_node(highest_spill_ixs, spill_edges, member_regions);

    neigh_regions = zeros(1, num_regions);
    neigh_regions(regions_not_in_trap(member_regions, spill_edges(highest_spill_ixs, :))) = 1;
    neigh_traps = find(neigh_regions * all_traps);
    highest_spill_edge_ixs = spill_edges(highest_spill_ixs, 1);

end

function spillp_ixs = purge_multiple_spills_same_node(spillp_ixs, spill_edges, ...
                                                      combined_regions)

    sp_e = [spill_edges(spillp_ixs, :), spillp_ixs]; % all spill point edges
    
    % orienting data of spill point edges so that the _first_ node belongs to
    % the region.
    sp_e_ori = ...
        [sp_e(ismember(sp_e(:,2), combined_regions), [6 5 8]); ...
         sp_e(ismember(sp_e(:,3), combined_regions), [7 4 8])];
    
    % sp_e_ori should now has the columns: 
    %  1) ix of node in _this_ region; 
    %  2) z-val of opposing node (in _other_ region);
    %  3) spillpoint ix

    sp_e_ori = sortrows(sortrows(sp_e_ori, 2), 1); % sort by node ix, then by
                                                   % z-value of opposing node
    spill_node_ixs = unique(sp_e_ori(:, 1));
    
    spillp_ixs = zeros(numel(spill_node_ixs), 1);
    
    for i = 1:numel(spillp_ixs)
        spillp_ixs(i) = sp_e_ori(find(sp_e_ori == spill_node_ixs(i), 1), 3);
    end
end


% function spillp_ixs = purge_multiple_spills_same_node(spillp_ixs, spill_edges)

%     nnh = [spill_edges(spillp_ixs, [6 5]), spillp_ixs];
%     nnh = sortrows(sortrows(nnh, 2), 1); % sort by node ix, then by z-value
%                                          % of opposing node
%     [~, retain_rows] = unique(nnh(:,1), 'stable');
%     spillp_ixs = nnh(retain_rows, 3);
    
% end

function res = regions_not_in_trap(trap_regions, involved_edges)
% trap_regions         - index of nodes belonging to the trap
% involved_spill_edges - index(ices) of the spill edge(s) corresponding to the highest spill value
%	                 of the region.

    % For each edge, column 2 and 3 in involved_edges gives the indices of the regions
    % it connects.
    involved_regions  = unique(reshape(involved_edges(:, 2:3), [], 1));

    % Return indices of the regions involved that are _not_ among the trap regions. (zero tacked
    % on to the group of region indices to exclude, since we only want indices of interior regions)
    res = setdiff(involved_regions, [trap_regions;0]);

end

function [res, loc_ix] = trap_spill_edges(s_edges, combined_regions, sort)
% Return the indices of all spill edges that belongs to the trap
% specified by the combined regions in 'regions'.  In other words, the indices
% of all edges that lead 'in/out of' the trap.
%
% PARAMETERS:
% s_edges          - vector of all spill edges in grid
% combined_regions - vector of indices identifying the local spill regions that belong
%                    to the trap in question	
% sort             - If 'sort' is true, edges are sorted according to z-value.
%
% RETURNS:
%   res   - a matrix with 2 columns, and one line per spill edge of
%           the combined regions.  The first column gives the global
%           index of the edge, the second gives the z-height of that edge,
%           defined as the z-value of the deepest of its two end nodes.
%  loc_ix - a vector of the corresponding local (row) indices into the 
%           s_edges matrix.
% 
% A spill edge is bordering the combined region specified by
% 'regions' if one of its nodes is inside one of the regions and
% the other node is _not_ inside any of the regions.  Hence the use
% of xor in the expression below.
    loc_ix = find(xor(ismember(s_edges(:,2), combined_regions), ...
                      ismember(s_edges(:,3), combined_regions)));

    reg_edges = s_edges(loc_ix, :);

    %     -glob. edge ix.-  -------------- z-value -------------
    res = [reg_edges(:,1), max(reg_edges(:,4), reg_edges(:,5))];

    if sort
        % We must sort res and loc_ix together, to ensure that rows correspond
        tmp = sortrows([res, loc_ix], 2);
        res = tmp(:, 1:2);
        loc_ix = tmp(:, 3);
    end
end
