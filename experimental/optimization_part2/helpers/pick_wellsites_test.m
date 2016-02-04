function [wc, qt] = pick_wellsites_test(Gt, rock2D, co2, ta, num_wells, ...
                                   rhoW, sw, seafloor_temp, seafloor_depth, ...
                                   tgrad, buf_dist, maximize_boundary_distance, ...
                                   domain_buf_dist, pick_highest_pt, ...
                                   ExternalTrapRegion_buf_dist, varargin)

% TEST well buffer dist, maximize bdry dist
% options:
%   a) keep wells away from external domain edge
%   b) keep wells away from internal trap_region bdrys (i.e., trap_region
%   edge that doesn't have another trap_region adjacent to it)
%   c) place wells furthest downslope as possible (while maintaining other
%   constraints) for enhanced residual trapping
%   d) place wells at highest elevation in trap_region (i.e., the top of
%   the structural trap
%   e) in case above constaints fail, do not place well in the particular
%   trap_region, but place it in next trap_region in tree if it exists,
%   otherwise do not place it in tree.
%

% OTHER options to try:
%   a) use an array of wells, only in trap regions. Maintain a trap_region
%   bdry buffer distance, and an outer domain buffer distance
%   b) place wells along edges of adjancet trap regions
%   c) qt should be computed based on shared volume capacities.
%   d) use multiple wells per spill-tree, in the case of large trapping
%   capacities that will be difficult to exploit in a 50 year injection
%   time frame.
%   e) place wells in the tops of each trap, regardless of the spill tree
%
%

% % Trap index of the root trap of each tree
% trees = maximizeTrapping(Gt, 'res', ta); % 'calculateAll', 'removeOverlap' ?
% treeRoots = [trees.root];
% 
% % get the trap indexes of all traps downstream from a tree root
% % (including the trap index of the root itself). Here, we pass in
% % poro * ntg as the porosity, to account for possible ntg.
% [trap_ixs, trap_vols] = downstreamTraps(Gt, rock2D.poro.*rock2D.ntg, ta, treeRoots(i));



    opt.inspectWellPlacement = false;
    opt = merge_options(opt, varargin{:});

    gravity on;

    % Computing local pressure and temperature conditions (in order to
    % compute CO2 densities) 
    P = rhoW .* norm(gravity) .* Gt.cells.z; % hydrostatic pressure
    T = seafloor_temp + (Gt.cells.z - seafloor_depth) .* tgrad ./ 1000;
    
    % Computing local CO2 densities
    rhoCO2 = co2.rho(P, T);

    % Computing trap volumes
    tcells               = find(ta.traps);
    trap_heights         = zeros(Gt.cells.num, 1);
    trap_heights(tcells) = ta.trap_z(ta.traps(tcells)) - Gt.cells.z(tcells);
    trap_heights         = min(trap_heights,Gt.cells.H);
    if ~isfield(rock2D,'ntg')
        rock2D.ntg = ones(Gt.cells.num,1); % account for possible NTG
    end 
    strap_vol = Gt.cells.volumes .* trap_heights .* rock2D.poro .* (1-sw) .* rock2D.ntg;

    assert(all(trap_heights <= Gt.cells.H));
    
    % Computing trap capacities in mass terms and suggesting wellsites
    trapcap = accumarray(ta.traps(tcells), strap_vol(tcells) .* rhoCO2(tcells));
    
    % Suggest wellsites based on maximum cumulative trap capacity reached

    bdist = sqr_distance_from_boundary(Gt); % used below to discourage placement
    bdist = bdist/max(bdist(:));            % of wells close to boundary

    wc = zeros(num_wells, 1);
    qt = wc;

    if opt.inspectWellPlacement
        figure
        plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1)
    end
    for i = 1:num_wells
        cumcap = [0; cumulative_capacity(trapcap, ta.trap_adj)];
        field  = cumcap(ta.trap_regions + 1);

        if maximize_boundary_distance
           % Pick well site as far as possible away from boundary.
        
           % Adding an insignificant quantity so that if everything else is
           % equal, the cell farthest from the boundary is chosen
%            grid_eps = 10 * eps(max(cumcap(:)));
%            field    = field + grid_eps * bdist;
%            [qt(i), wc(i)] = max(field);
           
           % The above approach may fail to capture the farthest cell due
           % to round-off (i.e., when values are the same within machine
           % precision, max() returns the first cell index of the max
           % values). Thus another approach is:
           [~, wc(i)] = max( field + bdist );
           qt(i) = field(wc(i));
           
        else
           % Call to select_wellcell: picks location according to user
           % options and while respecting buffer distances set
           
           % A while loop is used: candidate cells are located as per the
           % best cumulative capacity of a trap region, then a well cell is
           % attempted to be placed. If select_wellcell was successful in
           % placing a well s.t. it satisfies all buffer constaints, the
           % well is placed, otherwise wc(i) is empty and new candidate
           % cells are located in the next best trap region (in terms of
           % cumulative capacity). Once again, select_wellcell tries to
           % place a well within the constraints... Process continues until
           % all the wells are placed, or a well is assigned a rate of zero
           % (which implies no more wells should be placed).
           wc_tmp = [];
           while isempty(wc_tmp)
               
               qt_tmp = max(unique(field)); % field can be array of 0's at some point
               if qt_tmp > 0
                    candidate_cells = find(field == qt_tmp);
                    wc_tmp = select_wellcell(Gt, candidate_cells, buf_dist, domain_buf_dist, ...
                            pick_highest_pt, ExternalTrapRegion_buf_dist, ta, ...
                            opt.inspectWellPlacement);
               
                    field( field == qt_tmp ) = 0;
               else
                    wc_tmp = 0;
               end
               % if wc(i) is empty, none of the candidate cells satisfied the
               % buffer distance constaints. Thus, the field values of
               % these cells were set to 0 so as to remove them from the
               % potential placement cells, and we set the cells with the
               % next highest cumulative capacity to be the new candidate
               % cells.
           end
           qt(i) = qt_tmp;
           wc(i) = wc_tmp;
           
           if qt(i) == 0
                fprintf('No more wells to be placed.\n')
                % keep non-zero values computed so far
                wc = wc(1:i-1); qt = qt(1:i-1);
                % ensure only unique well cell indexes
                [wc, inx] = unique(wc,'stable'); % unsorted
                qt        = qt(inx);
                % return from this function
                return;
           end
           
        end
        
        % Setting available capacity of used traps to zero
        trapcap = set_used_to_zero(trapcap, ta.trap_adj, ta.trap_regions(wc(i)));
        
        % check if there is no more avaible capacity for the wells (which
        % can occur if all traps have been used)
        if all(trapcap == 0) || qt(i)/qt(1) < 0.01
            if all(trapcap == 0)
                fprintf(['You wanted to place %d wells, '...
                'but all traps were used after placing %d wells.\n'], num_wells, i)
            elseif qt(i)/qt(1) < 0.01
                fprintf(['You wanted to place %d wells, '...
                'but the injected volume of well number %d was '...
                'less than 1 percent of well number 1''s volume.\n'], num_wells, i)
            end
            % keep values computed so far
            wc = wc(1:i); qt = qt(1:i);
            % ensure only unique well cell indexes
            [wc, inx] = unique(wc,'stable'); % unsorted
            qt        = qt(inx);
            % return from this function
            return;
        end
    end
    
    % take only the unique well cell indexes
    [wc, inx] = unique(wc,'stable');
    qt = qt(inx);
    num_wells = numel(wc); % updated number of wells
    assert( numel(wc) == numel(qt) )

%     % Add outlines of trap regions
%     if opt.inspectWellPlacement
%     for i = 1:numel(unique(ta.trap_regions( ta.trap_regions > 0 )))
%         plotFaces(Gt, boundaryFaces(Gt, ta.trap_regions == i), 'EdgeColor','r')
%     end
%     end
    

end

% HELPER FUNCTIONS:

% ----------------------------------------------------------------------------
function wc = select_wellcell(Gt, candidates, buffer, domain_buffer, pick_highest_pt, ...
    ExternalTrapRegion_buf_dist, ta, inspectWellPlacement)

   % Identify cells at the boundary of the selected region
   num    = Gt.cells.num;
   n_rels = Gt.faces.neighbors;
   n_rels = double(n_rels(~any(n_rels == 0,2),:)); % remove exterior face relations
   adj    = sparse(n_rels(:,1), n_rels(:,2), 1, num, num, 2 * size(n_rels,1));
   adj    = adj + adj';
   boundary_candidates = candidates(sum(adj(candidates, candidates)) < 4);
   interior_candidates = setdiff(candidates, boundary_candidates);
   
   % compute smallest distance for each interior candidate to boundary
   distances = inf(numel(interior_candidates), 1);
   for bix = boundary_candidates(:)'
      dvec = bsxfun(@minus, ...
                    Gt.cells.centroids(bix, 1:2), ...
                    Gt.cells.centroids(interior_candidates, 1:2));
      dist = sqrt(sum(dvec.^2, 2));
      distances = min(distances, dist);
   end
   
   % Only interior candidates are assessed (not boundary candidates)
   inside_buffer_candidates = interior_candidates(distances > buffer);
   
   % Additional check to ensure inside_buffer_candidates are within:
   % (1) a buffer distance to formation boundary, and
   % (2) a buffer distance to the interior trap_region boundary (i.e., the
   %     boundary between the trap_regions and no trap_regions).
   bdist = sqr_distance_from_boundary(Gt);
   bdist = sqrt(bdist(inside_buffer_candidates));
   inside_buffer_candidates = inside_buffer_candidates(bdist > domain_buffer);
   if ~isempty(inside_buffer_candidates)
        bdist = distance_from_TrapRegionBoundary(Gt, ta);
        bdist = bdist(inside_buffer_candidates);
        inside_buffer_candidates = inside_buffer_candidates(bdist > ExternalTrapRegion_buf_dist);
        %
        if inspectWellPlacement
            plotCellData(Gt, Gt.cells.volumes, boundary_candidates)
            plotCellData(Gt, Gt.cells.volumes, interior_candidates, 'FaceColor','r')
            plotCellData(Gt, Gt.cells.volumes, inside_buffer_candidates, 'FaceColor','y')
        end
   end
   
   % If inside_buff_candidates is not empty, the well is placed using
   % either the highest point if true, otherwise using the furthest
   % downslope. If inside_buffer_candidates are empty, an empty wc is
   % passed out, and the next best leaf-node will be assessed.

   if isempty(inside_buffer_candidates)
      % none of the candidates were far enough from the boundary. An empty
      % wc is passed out, to flag the assessement of the next leaf-node in
      % the spill-tree.
      wc = [];
      
   elseif pick_highest_pt
      % choose the highest point
      candidate_z = Gt.cells.z(candidates);
      [~, min_ix] = min(candidate_z);
      wc = candidates(min_ix);
      
   else
      % the furthest downslope of the inside_buffer_candidates is chosen
      inside_buffer_z = Gt.cells.z(inside_buffer_candidates);
      [~, max_z_ix] = max(inside_buffer_z);
      
      wc = inside_buffer_candidates(max_z_ix);
      %
      if inspectWellPlacement
        plotCellData(Gt, Gt.cells.volumes, [wc wc], 'FaceColor','k')
        plotFaces(Gt, boundaryFaces(Gt, ta.trap_regions == ta.trap_regions(wc)), ...
            'EdgeColor','k','LineWidth',3)
      end
   end
end


% ----------------------------------------------------------------------------

function bdist = distance_from_TrapRegionBoundary(Gt, ta)

   % The candidate cells are all the cells belonging to a trap region:
   candidates = find( ta.trap_regions > 0 );

   % Identify cells at the boundary of the candidate cells (will it work
   % for disconnected trap regions?) @@
   num    = Gt.cells.num;
   n_rels = Gt.faces.neighbors;
   n_rels = double(n_rels(~any(n_rels == 0,2),:)); % remove exterior face relations
   adj    = sparse(n_rels(:,1), n_rels(:,2), 1, num, num, 2 * size(n_rels,1));
   adj    = adj + adj';
   boundary_candidates = candidates(sum(adj(candidates, candidates)) < 4);
   interior_candidates = setdiff(candidates, boundary_candidates);
   
   % compute smallest distance for each interior candidate to boundary
   distances = inf(numel(interior_candidates), 1);
   for bix = boundary_candidates(:)'
      dvec = bsxfun(@minus, ...
                    Gt.cells.centroids(bix, 1:2), ...
                    Gt.cells.centroids(interior_candidates, 1:2));
      dist = sqrt(sum(dvec.^2, 2));
      distances = min(distances, dist);
   end

   % Pass out bdist in an array for each Gt.cell (bdist will be 0 for any
   % cell that is not part of a trap region)
   bdist = zeros(Gt.cells.num,1);
   bdist(interior_candidates) = distances;

end


% ----------------------------------------------------------------------------

function cumcap = cumulative_capacity(trapcap, adj)
    % Compute cumulative capacity of each trap and those it is connected to upstream
    cumcap = trapcap;
    adj_orig = adj;
    while nnz(adj)
        cumcap = cumcap + adj * trapcap;
        adj = spones(adj * adj_orig);
    end
end


% ----------------------------------------------------------------------------

function trapcap = set_used_to_zero(trapcap, adj, trap_ix)
    % Set capacity of trap with index 'trap_ix' and all its ustream traps to zero
    trapcap(trap_ix) = 0;
    ustr_traps = find(adj(trap_ix,:));
    while ~isempty(ustr_traps)
        trapcap(ustr_traps) = 0;
        ustr_traps = find(sum(adj(ustr_traps, :), 1));
    end
end

% ----------------------------------------------------------------------------

function [field, closest_bnd_edge] = sqr_distance_from_boundary(Gt)

    field = inf * ones(Gt.cells.num, 1);
    closest_bnd_edge = zeros(Gt.cells.num, 1);
    
    % Looping over boundary edges
    bnd_edges = find(prod(Gt.faces.neighbors, 2) == 0);
    for b = bnd_edges'
        dx = Gt.cells.centroids(:,1) - Gt.faces.centroids(b,1);
        dy = Gt.cells.centroids(:,2) - Gt.faces.centroids(b,2);
        
        d2 = dx.*dx + dy.*dy;
        replace = (d2 < field);
        field(replace) = d2(replace);
        closest_bnd_edge(replace) = b;
    end
    
end