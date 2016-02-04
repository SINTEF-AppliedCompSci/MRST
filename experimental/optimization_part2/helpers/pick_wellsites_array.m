function [wc, qt] = pick_wellsites_array(Gt, rock2D, co2, ta, ...
                                   rhoW, sw, seafloor_temp, seafloor_depth, tgrad, ...
                                   domain_buffer, catchment_buffer, ...
                                   DX, DY, max_num_wells, varargin)

% OTHER options to try:
%   a) use an array of wells, only in trap regions. Maintain a trap_region
%   bdry buffer distance, and an outer domain buffer distance


    % DX and DY is spacing for array of inj wells

    opt.inspectWellPlacement = false;
    opt = merge_options(opt, varargin{:});
    
    %% Determine candidates
    % Candidates are those cells which are located in a trap region and
    % satisfy the buffer distances:
    
    candidates = find( ta.trap_regions > 0 );
    [bdist, ~] = distance_from_candidates_to_their_bdry(Gt, candidates);
    inside_catchment_buffer = candidates(bdist > catchment_buffer);
    
    candidates = [1:Gt.cells.num]';
    [bdist, ~] = distance_from_candidates_to_their_bdry(Gt, candidates);
    inside_domain_buffer = candidates(bdist > domain_buffer);
    
    % take the candidates that satify both buffer constraints:
    inside_buffer_candidates = intersect(inside_catchment_buffer, ...
                                         inside_domain_buffer);
    % @@ debug for situation where inside_buffer_candidates is empty!

    % An array of wells is placed in these candidates by first placing an
    % array of wells all over the formation top surface, and then keeping
    % only the wc that lie in the candidate cells. We ensure a wc is
    % separated from its adjacent wc by at least 2 cells. This can be done
    % by ensuring the well cells make up less than 30% of the total number
    % of candidate cells. Alternatively, well_density_limit could be passed
    % in as a varargin.
    
    added_spacing = 0;
    while 1
        
        % wc covering top surface
        wc = wells_covering_formation(Gt, DX + added_spacing, DY + added_spacing);
        % keep only the wc located in candidate cells
        wc = intersect( wc, inside_buffer_candidates);
        
        % If the distribution of wells in the candidate region is too
        % dense, or more than the well limit have been placed, the spacing
        % between wells is increased for the next iteration. Otherwise, the
        % while loop is exited
        well_density = (numel(wc)/numel(inside_buffer_candidates));
        if well_density > 0.30 || numel(wc) >= max_num_wells
            added_spacing = added_spacing + 1 * kilo * meter;
        else
            break
        end

    end
    
    if opt.inspectWellPlacement
       figure
       plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1)
       plotCellData(Gt, Gt.cells.volumes, wc, 'FaceColor','k')
    end
    
    
    %% Assign initial volume to inject:
    % Loop through each trap region, and count the number of wells placed.
    % If more than 0, assign a well rate based on the trap region capacity
    % and the number of wells. If 0, go to next trap region
    
    cells_trapcap = struct_trap_cap(Gt, rock2D, ta, rhoW, sw, seafloor_temp, ...
                                        seafloor_depth, tgrad, co2); % kg
    
    qt = [];
    tmp = 1:Gt.cells.num;
    for i = 1:numel(unique(ta.trap_regions( ta.trap_regions > 0 )))
        
        trap_region_cap = sum( cells_trapcap( ta.trap_regions == i ) ); % kg
        num_wells = numel(intersect(wc, tmp(ta.trap_regions == i)));
        
        if num_wells > 0
            qt_per_well = trap_region_cap / num_wells;
            qt = [qt; repmat(qt_per_well, [num_wells, 1])]; % kg
        else
            continue 
        end
        
    end
    assert( numel(qt) == numel(wc) )
    
    
end

function wc = wells_covering_formation(Gt, DX, DY)

    all_xcoords = Gt.cells.centroids(:,1);
    all_ycoords = Gt.cells.centroids(:,2);
    wX = [min(all_xcoords) : DX : max(all_xcoords)]';
    wY = [min(all_ycoords) : DY : max(all_ycoords)]';
    Xcoord = repmat(wX, numel(wY),1);
    Ycoord = [];
    for i = 1:numel(wY)
        Ycoord = [Ycoord; repmat(wY(i), numel(wX),1)];
    end
    wc = getCellIndex(Gt, Xcoord, Ycoord); wc = unique(wc);
    
end

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

% -------------------------------------------------------------------------

function cellIndex = getCellIndex(Gt, Xcoord, Ycoord)
% Closest cell index of grid Gt corresponding to physical coordinate (X,Y)

    assert(numel(Xcoord) == numel(Ycoord));
    cellIndex = zeros(numel(Xcoord),1);
    
    for i = 1:numel(Xcoord) % @@ remove for loop
        
        dv        = bsxfun(@minus, Gt.cells.centroids(:,1:2), [Xcoord(i), Ycoord(i)]);
        [v, ind]  = min(sum(dv.^2, 2));
        cellIndex(i) = ind; 

    end
    
end

% -------------------------------------------------------------------------

function [ Xcoord, Ycoord ] = getXYcoords(Gt, cellIndex)

    Xcoord = Gt.cells.centroids(cellIndex,1);
    Ycoord = Gt.cells.centroids(cellIndex,2);

end