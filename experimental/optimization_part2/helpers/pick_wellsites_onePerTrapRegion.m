function [wc, qt] = pick_wellsites_onePerTrapRegion(Gt, rock2D, co2, ta, rhoW, ...
                        sw, seafloor_temp, seafloor_depth, tgrad, ...
                        max_num_wells, pick_highest_pt, ...
                        well_buffer_dist, well_buffer_dist_domain, ...
                        well_buffer_dist_catchment, varargin)
                    
                    
% Wells are placed one by one in each trap region, and injection volumes
% (or masses) are determined by the trap region's capacity, not by the
% cumulative capacity of the connecting upstream traps.

% In the trap region, a well can be placed either at the highest elevation
% point, or furthest downslope within the buffer distances. For
% both options, wells are placed only if the buffer distances are
% satisfied, otherwise no well is placed in the trap region.

% Wells are placed one by one starting with the trap region that has the
% highest capacity. No more wells are placed once the max well limit has
% been reached, no more trap regions exist, or the trap region capacity is
% less than 1 % of the highest trap region's capacity (in which case it is
% not considered profitable to place any more wells). @@ this threshold
% could be increased!

    opt.inspectWellPlacement = false;
    opt = merge_options(opt, varargin{:});
    
    
    % Capacity of each trap region
    cells_trapcap = struct_trap_cap(Gt, rock2D, ta, rhoW, sw, seafloor_temp, ...
                                        seafloor_depth, tgrad, co2); % kg
    tcells  = find(ta.traps);
    trapcap = accumarray(ta.traps(tcells), cells_trapcap(tcells));
    % trapcap(1) corresponds to trap_region 1, etc...
    
    % Loop through every trap region, and place a well if buffer distances
    % are met. Otherwise, no well is placed. Iterations go until all trap
    % regions assessed, or until max number of wells placed, or until
    % further well injection masses do not meet threshold.
    wc = []; qt = [];
    for i = 1:numel(trapcap)

        % Trap_region ID with the highest available capacity
        [qt_tmp, trapID] = max(trapcap);

        % Candidate cells from that trap_region
        candidates = find( ta.trap_regions == trapID );
        [bdist, ~] = distance_from_candidates_to_their_bdry(Gt, candidates);
        inside_trapRegion_buffer = candidates(bdist > well_buffer_dist);

        candidates = find( ta.trap_regions > 0 );
        [bdist, ~] = distance_from_candidates_to_their_bdry(Gt, candidates);
        inside_catchment_buffer = candidates(bdist > well_buffer_dist_catchment);    

        candidates = [1:Gt.cells.num]';
        [bdist, ~] = distance_from_candidates_to_their_bdry(Gt, candidates);
        inside_domain_buffer = candidates(bdist > well_buffer_dist_domain);

        % take the candidates that satify all buffer constraints:
        inside_buffer_candidates = intersect(inside_trapRegion_buffer, ...
                                             inside_catchment_buffer);
        inside_buffer_candidates = intersect(inside_buffer_candidates, ...
                                             inside_domain_buffer);

        % place a well if candidates are not empty, and assign it's
        % injection mass
        if ~isempty(inside_buffer_candidates)
            
            % Choose either the highest point of the candidates or the
            % furthest downslope
            candidate_z = Gt.cells.z(inside_buffer_candidates);
            if pick_highest_pt
                [~, min_ix] = min(candidate_z);
                wc_tmp = inside_buffer_candidates(min_ix);
            else
                [~, max_ix] = max(candidate_z);
                wc_tmp = inside_buffer_candidates(max_ix);
            end
            
            % Store results, and set the used trap region's capacity to 0
            wc = [wc; wc_tmp];
            qt = [qt; qt_tmp];
            trapcap(trapID) = 0;
            
            % Determine whether to continue iteration loop
            if numel(wc) > max_num_wells
                fprintf('No more wells placed: max number of wells have been placed.\n')
                % return from this function
                return;
            elseif qt_tmp/qt(1) < 0.01
                fprintf(['No more wells placed: injected mass of the last well was '...
                    'less than 1 percent of the first well.\n'])
                % remove the last wc and its injection mass
                wc = wc(1:end-1); qt = qt(1:end-1);
                % return from this function
                return;
            elseif i == numel(trapcap)
                fprintf(['All trap regions were assessed. ',...
                    'A well was placed in %d out of %d trap regions.\n'], ...
                    numel(wc), numel(trapcap))
            end
            
        else
            % Even though no well was placed in the trap region, we still
            % will set it's capacity to 0 so it will not be considered
            % again for a well placement
            trapcap(trapID) = 0;
        end
        
    end


end

% HELPER FUNCTIONS:
% -------------------------------------------------------------------------

function cap = struct_trap_cap(Gt, rock2D, ta, rhoW, sw, seafloor_temp, seafloor_depth, tgrad, co2)

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
    
    % Trap capacities in mass terms, for all cells (may be 0)
    cap = strap_vol .* rhoCO2; % kg

end


% -------------------------------------------------------------------------

function [bdist, distances_full] = distance_from_candidates_to_their_bdry(Gt, candidates)

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


   % Pass out distances in an array for each Gt.cell (distances_full will
   % be 0 for any cell that is not part of the candidate region)
   distances_full = zeros(Gt.cells.num,1);
   distances_full(interior_candidates) = distances;
   distances_full(boundary_candidates) = 0;
   
   % Pass out distances for each candidate cell (which includes the 0
   % distance along the boundary of the candidates):
   bdist = zeros(numel(candidates),1);
   bdist = distances_full(candidates);
   assert( numel(bdist) == numel(candidates) );

end
