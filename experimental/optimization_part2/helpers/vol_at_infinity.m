function [ ult_vol_remaining, ult_vol_leaked ] = vol_at_infinity( Gt, rock2D, state, ta, sw )
% Determine co2 vol that will ultimately remain in formation at time infinity.

% This calculation of co2 vol at time infinity is based on the formations'
% co2 vol at a point in time in which it is safe to assume the flow
% dynamics are gravity-dominated. Thus, the co2 vol in a particular spill
% tree at time infinity is based on that spill tree's capacity, and any
% extra co2 vol is assumed to have been leaked by time infinity.

% % state - a structure of states, including co2 saturation (i.e., states.s(:,2))
% % sw - 
% sw = 0.11; % @@ pass in as input variable
% 
% % Gt
% % rock2D
% % state
% % ta
% % sw
% state = states{end};
% rock2D = other.rock;
% sw = other.residual(1);

    figure;
    plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1);
 
    % 1) get and show current co2 vol: 
    curr_vol = (state.s(:,2) .* Gt.cells.volumes .* Gt.cells.H .* rock2D.poro); % m3 pore vol
    cells = 1:Gt.cells.num;
    plotCellData(Gt, curr_vol, cells(curr_vol>0.1), 'EdgeAlpha',0.1);
    title('current co2 vol [m3]')
    colorbar; axis equal tight off
    
    % show trap regions:
    figure;
    mapPlot(gcf, Gt, 'traps', ta.traps, ...
            'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
            'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
            'maplines', 20);
    colorizeCatchmentRegions(Gt, ta); % this actually colorizes regions around each trap
    axis equal tight off
    title('trap regions')
    
    % get and show structural trap vol:
    % Computing trap volumes
    tcells               = find(ta.traps);
    trap_heights         = zeros(Gt.cells.num, 1);
    trap_heights(tcells) = ta.trap_z(ta.traps(tcells)) - Gt.cells.z(tcells);
    trap_heights         = min(trap_heights, Gt.cells.H);
    strap_vol            = Gt.cells.volumes .* trap_heights .* rock2D.poro .* (1-sw);
    
    figure; plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1);
    plotCellData(Gt, strap_vol, 'EdgeAlpha',0.1)
    colorbar; axis equal tight off
    title('structural trap pore vol [m3]')
    
    assert(all(trap_heights <= Gt.cells.H));
    % Computing trap capacities in vol terms
    trapcap = accumarray(ta.traps(tcells), strap_vol(tcells)); % m3 %@@ result is different from vol returned using downstreamTraps()... look at computeTrapVolume()
    % show:
    figure; plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1);
    num_traps = max(unique(ta.traps));
    for i=1:num_traps
        % trap vol (of a single trap)
        tc = sum(strap_vol(ta.traps == i)); % m3
        plotCellData(Gt, tc*ones(Gt.cells.num,1), ta.traps == i, 'EdgeAlpha',0.1)
        %pause(1)
    end
    title('trap capacities [m3]'); colorbar; axis equal tight off

    
    % 2) current co2 vol in each trap region
    trcells = find(ta.trap_regions);
    curr_vol_per_trap_region = accumarray(ta.trap_regions(trcells), curr_vol(trcells));
    
    figure; plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1);
    num_traps = max(unique(ta.traps));
    for i=1:num_traps
        % curr vol per trap region
        plotCellData(Gt, curr_vol_per_trap_region(i)*ones(Gt.cells.num,1), ...
            ta.traps == i, 'EdgeAlpha',0.1)
        %pause(1)
    end
    title('curr vol per trap region [m3]'); colorbar; axis equal tight off
    
    
    % 3) get the trap index of the root trap of each tree
    trees = maximizeTrapping(Gt, 'res', ta); % 'calculateAll', 'removeOverlap' ?
    % NB: sth buggy in the order of traps for Sto's first spill-pathway, so
    % the trap indexes of traps along a spill-path (i.e., what trees.traps
    % should have contained) are obtained using another function:
    treeRoots = [trees.root];
    
    % 4) here we determine the vol leaked per tree by spilling any vols a
    % trap cannot retain down the tree, until the vol is ultimately spilled
    % out of the domain. During this calculation, we also determine the
    % co2 vol ultimately remaining in the domain
    vol_leaked_per_tree = zeros(numel(treeRoots),1);
    ult_vol_per_trap_region = curr_vol_per_trap_region;
    for i = 1:numel(treeRoots)
        % get the trap indexes of all traps downstream from a tree root
        % (including the trap index of the root itself)
        [trap_ixs, trap_vols] = downstreamTraps(Gt, rock2D.poro, ta, treeRoots(i));
        % NB: trap_vols do not match with trapcap (m3) since trapcap
        % accounted for residual water (sw)! @@ but still there is a
        % discrepancy
        
        for j = 1:numel(trap_ixs)
            % calculate the vol that will spill out of trap_ixs(j)
            spill_vol = max(0, - trapcap(trap_ixs(j)) + ult_vol_per_trap_region(trap_ixs(j)));
            
            if j < numel(trap_ixs) % haven't reach end of spill-tree yet
                % add this spill volume to the next trap down the
                % spill-tree (i.e., trap_ixs(j+1))
                ult_vol_per_trap_region(trap_ixs(j+1)) = ult_vol_per_trap_region(trap_ixs(j+1)) + spill_vol;
            else
                % we've reached the end of the spill-tree, thus the
                % spill_vol goes out of the spill-tree (i.e., out of the
                % domain)
                vol_leaked_per_tree(i) = spill_vol;
            end
            % NB: ult_vol_per_trap_region gets updated each time a
            % spill_vol is added to the next (downsteam) trap. After
            % looping through all tree roots, the final array should
            % represent the amount of co2 per trap region at stable
            % conditions. No more migration should occur under
            % gravity-dominated dynamics.
        end
    end
    ult_vol_remaining = sum(ult_vol_per_trap_region);
    
    % 5) In addition to the vols that will spill-out of each tree, any co2
    % vol that existed outside of the trap regions making up the spill
    % trees is bound to leak from the domain.
    curr_vol_outside_trap_regions = sum(curr_vol) - sum(curr_vol_per_trap_region);
    
    % 6) The final total amount of co2 that will ultimately leak from the
    % domain as stable conditions are reached, or at time infinity:
    ult_vol_leaked = sum(vol_leaked_per_tree) + curr_vol_outside_trap_regions; % m3
    
    % 7) 
    assert( sum(curr_vol) == (ult_vol_remaining + ult_vol_leaked) ) % @@ mis-match in vols!
    
end

