function [ ult_vol_remaining, ult_vol_leaked ] = vol_at_infinity( Gt, rock2D, sG, sw, sr, varargin )
% Determine co2 vol that will ultimately remain in formation at time infinity.
%
% SYNOPSIS:
%   future_vol = vol_at_infinity(Gt, rock2D, states{i}.s(:,2), 0.11)
%   future_vol = vol_at_infinity(Gt, rock2D, states{i}.s(:,2), 0.11, 'plotsOn', true)
%
% DESCRIPTION:
%   This calculation of co2 vol at time infinity is based on the
%   formations' co2 vol at a point in time in which it is safe to assume
%   the flow dynamics are gravity-dominated. Thus, the co2 vol in a
%   particular spill tree at time infinity is based on that spill tree's
%   capacity, and any extra co2 vol is assumed to have been leaked by time
%   infinity.
%
%   NB: this routine was developed to handle ADI variables, thus the syntax
%   used here respects ADI variables, and in some instances a cell array of
%   ADI variables.

% INPUTS    - Gt, top surface grid
%           - rock2D, includes porosity (and possibly ntg) data
%           - sG, current (or i-th) saturation of co2, i.e., states{i}.s(:,2)
%           - sw, residual saturation of water, i.e., 0.11
%           - sr, residual saturation of co2, i.e., 0.21
%           - 'plotsOn', true or false for optional plotting
%           - (optionally) ta, as given by trapAnalysis. If not passed in,
%           will be computed inside this routine.

% RETURNS   - ult_vol_remaining, a total volume in units of m3 co2. This is
%           the amount of co2 remaining in the formation at time infinity
%           - ult_vol_leaked, a total volume in units of m3 co2, of the
%           amount leaked from the current state (i.e., sG)

% Needs more testing with different formations!
% Assess sG versus sG .* model.fluid.bG(p)
% Assess rock2D.poro versus rock2D.poro .* model.fluid.pvMultR(p)

    opt.plotsOn = false;
    opt.ta = [];
    opt = merge_options(opt, varargin{:});

    if isempty(opt.ta)
        fprintf('Obtaining trap structure using trapAnalysis...\n')
        ta = trapAnalysis(Gt, false); % @@ implement option to pass in closed boundary faces.
        fprintf('trapAnalysis done.\n')
    else
        ta = opt.ta;
    end
    
    % To account for possible NTG data:
    poro = rock2D.poro;
    ntg  = ones(Gt.cells.num,1);
    if isfield(rock2D,'ntg')
        ntg = rock2D.ntg;
    end
    
    
    % 1) Get current co2 vol: 
    curr_vol = ( sG .* Gt.cells.volumes .* Gt.cells.H .* poro .* ntg ); % m3 pore vol
    tot_vol = sum(curr_vol);
    cells = 1:Gt.cells.num;
    
    % get structural trap co2 vol:
    % Computing trap co2 volume capacities
    tcells               = find(ta.traps);
    trap_heights         = zeros(Gt.cells.num, 1);
    trap_heights(tcells) = ta.trap_z(ta.traps(tcells)) - Gt.cells.z(tcells);
    trap_heights         = min(trap_heights, Gt.cells.H);
    strap_vol            = Gt.cells.volumes .* trap_heights .* poro .* (1-sw) .* ntg; % without the (1-sw), strap_vol is the structural trap pore volume
    
    assert(all(trap_heights <= Gt.cells.H));
    % Computing trap capacities in vol terms
    %trapcap = accumarray(ta.traps(tcells), strap_vol(tcells)); % m3 %@@ result is different from vol returned using downstreamTraps()... look at computeTrapVolume()
    %@@ accumarray isn't defined for an ADI variable
    num_traps = max(unique(ta.traps));
    %trapcap = zeros(num_traps,1);
    %trapcap = 0*tot_vol;
    trapcap = cell(num_traps,1); % use cell array since trapcap is an ADI var when strap_vol is an ADI
    for i=1:num_traps
        trapcap{i} = sum(strap_vol(ta.traps == i)); % m3
    end
    

    % 2) Current co2 vol in each trap region
    assert( max(unique(ta.traps)) == max(unique(ta.trap_regions)) )
    trcells = find(ta.trap_regions);
    %curr_vol_per_trap_region = accumarray(ta.trap_regions(trcells), curr_vol(trcells));
    %@@ accumarray isn't defined for an ADI variable
    num_traps = max(unique(ta.traps));
    %curr_vol_per_trap_region = zeros(num_traps,1);
    curr_vol_per_trap_region = cell(num_traps,1); % use cell array since trapcap is an ADI var when strap_vol is an ADI
    for i=1:num_traps
        curr_vol_per_trap_region{i} = sum(curr_vol(ta.trap_regions == i)); % m3
    end
    
    
    % 3) Get the trap index of the root trap of each tree
    trees = maximizeTrapping(Gt, 'res', ta); % 'calculateAll', 'removeOverlap' ?
    % NB: sth buggy in the order of traps for Sto's first spill-pathway, so
    % the trap indexes of traps along a spill-path (i.e., what trees.traps
    % should have contained) are obtained using another function:
    treeRoots = [trees.root];
    
    
    % 4) Here we determine the vol that is ultimately leaked per tree by
    % spilling any vols a trap cannot retain down the tree, and out of the
    % domain. During this calculation, we also determine the co2 vol
    % ultimately remaining in the domain
    %vol_leaked_per_tree = zeros(numel(treeRoots),1);
    vol_leaked_per_tree = cell(numel(treeRoots),1);
    ult_vol_per_trap_region = curr_vol_per_trap_region;
    clear i
    %fprintf('There are %d trees. \n\n', numel(treeRoots))
    dispif(mrstVerbose, 'There are %d trees. \n\n', numel(treeRoots));
    for i = 1:numel(treeRoots)
        dispif(mrstVerbose, '\nLooking at tree %d... \n', i);

        % get the trap indexes of all traps downstream from a tree root
        % (including the trap index of the root itself). Here, we pass in
        % poro * ntg as the porosity, to account for possible ntg.
        [trap_ixs, trap_vols] = downstreamTraps(Gt, poro.*ntg, ta, treeRoots(i));
        % NB: trap_vols (m3 'reservoir' pore space) do not match with
        % trapcap (m3 co2) since trapcap is vol of co2; to get the same
        % vaule, we must account for residual water (sw).
        trap_co2_vols = trap_vols' .* (1-sw);
        trapcap_tmp = cellfun( @(x) double(x), {trapcap{trap_ixs}})';
        assert( all( abs( trap_co2_vols - trapcap_tmp )./abs(trap_co2_vols) < 1e-8 ) ) % syntax handles cell array trapcap
        
        dispif(mrstVerbose, 'The trap indexes in this tree are: %s \n\n', ...
            num2str(trap_ixs));
        for j = 1:numel(trap_ixs)
            dispif(mrstVerbose, 'Looking at trap %d... \n', trap_ixs(j));
            
            dispif(mrstVerbose, 'Trap %d has a capacity of %d (m3 co2). \n', ...
                trap_ixs(j), double(trapcap{trap_ixs(j)}) );
            %fprintf('Trap %d initially had %d (m3 co2). \n', trap_ixs(j), trap_co2_vols(j) )
            dispif(mrstVerbose,'Trap %d currently has %d (m3 co2). \n', ...
                trap_ixs(j), double(ult_vol_per_trap_region{trap_ixs(j)}) );
            
            
            % calculate the vol that will spill out of trap_ixs(j)
            spill_vol = max(0, - trapcap{trap_ixs(j)} + ult_vol_per_trap_region{trap_ixs(j)}); % may be an ADI var
            assert( spill_vol >= 0 )
            
            % and remove this vol from trap_ixs(j)
            ult_vol_per_trap_region{trap_ixs(j)} = ult_vol_per_trap_region{trap_ixs(j)} - spill_vol; % may be an ADI var
            
            % spill this "spill vol" into the next trap downstream or out of
            % the tree
            if j < numel(trap_ixs)
                % haven't reach end of spill-tree yet, so add this spill
                % volume to the next trap down the spill-tree (i.e.,
                % trap_ixs(j+1))
                dispif(mrstVerbose, 'Thus %d (m3 co2) spills out of trap %d, and into trap %d.\n', ...
                    double(spill_vol), trap_ixs(j), trap_ixs(j+1));
                ult_vol_per_trap_region{trap_ixs(j+1)} = ult_vol_per_trap_region{trap_ixs(j+1)} + spill_vol; % may be an ADI var
                dispif(mrstVerbose, 'Trap %d now has %d (m3 co2). \n', ...
                    trap_ixs(j+1), double(ult_vol_per_trap_region{trap_ixs(j+1)}) );
            else
                % we've reached the end of the spill-tree, thus the
                % spill_vol goes out of the spill-tree (i.e., out of the
                % domain)
                dispif(mrstVerbose, 'Thus %d (m3 co2) spills out of trap %d, and out of tree %d.\n', ...
                    double(spill_vol), trap_ixs(j), i);
                vol_leaked_per_tree{i} = spill_vol; % may be an ADI var
                
            end
            % NB: ult_vol_per_trap_region gets updated each time a
            % spill_vol is removed from a trap and added to the next
            % (downsteam) trap. After looping through all tree roots, the
            % final array should represent the amount of co2 per trap
            % region at stable conditions. No more migration should occur
            % under gravity-dominated dynamics.
        end
    end
    
    % a for loop is used to take the sum of a cell array of ADI variables:
    % @@ is there a better way to do this?
    ult_vol_remaining = ult_vol_per_trap_region{1};
    for i=2:numel(ult_vol_per_trap_region)
        ult_vol_remaining = ult_vol_remaining + ult_vol_per_trap_region{i};
    end
    %ult_vol_remaining = sum( cellfun( @(x) (x), ult_vol_per_trap_region ) ); % syntax handles a cell array, but loses ADI jacobians
    %ult_vol_remaining = ones(1, Gt.cells.num) * ult_vol_per_trap_region;
    
    
    % 5) In addition to the vols that will spill-out of each tree, any co2
    % vol that existed outside of the trap regions making up the spill
    % trees is bound to leak from the domain. However, we account for an
    % amount of residually trapped co2 that will remain and is therefore
    % not bound to leak.
    total_trap_region_vol = curr_vol_per_trap_region{1};
    for i=2:numel(curr_vol_per_trap_region)
        total_trap_region_vol = total_trap_region_vol + curr_vol_per_trap_region{i};
    end
    curr_vol_outside_trap_regions = sum(curr_vol) - total_trap_region_vol;
    resid_trap_vol_outside_trap_regions = curr_vol_outside_trap_regions * sr;
    
    % 6) The final total amount of co2 that will ultimately leak from the
    % domain as stable conditions are reached, or at time infinity:
    vol_leaked_total_trees = vol_leaked_per_tree{1};
    for i=2:numel(vol_leaked_per_tree)
        vol_leaked_total_trees = vol_leaked_total_trees + vol_leaked_per_tree{i};
    end
    ult_vol_leaked = vol_leaked_total_trees + ...
        curr_vol_outside_trap_regions - resid_trap_vol_outside_trap_regions; % m3, may be an ADI var
    %assert( ult_vol_leaked >= 0 )
    assert( ult_vol_leaked > -1e-5 ) % to account for possible numerical round-off error due to spilling of vols
    ult_vol_remaining = sum(curr_vol) - ult_vol_leaked;
    
    % 7) Check that vols balance out, within some tolerance
    % @@ or they should balance out exactly?
    assert( abs( sum(curr_vol)-(ult_vol_remaining + ult_vol_leaked) ) ./ abs(sum(curr_vol)) < 1e-8 )
    
    
    % 8) Plotting (optional)
    if opt.plotsOn
        
        figure;
        plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1);
        plotCellData(Gt, curr_vol, cells(curr_vol>0.1), 'EdgeAlpha',0.1);
        title('current co2 vol [m3]')
        colorbar; axis equal tight off

        figure;
        mapPlot(gcf, Gt, 'traps', ta.traps, ...
                'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
                'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
                'maplines', 20);
        colorizeCatchmentRegions(Gt, ta); % this actually colorizes regions around each trap
        axis equal tight off
        title('trap regions')
        
        figure; set(gcf,'Position',[1 1 1647 526])
        subplot(1,3,1); plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1);
        num_traps = max(unique(ta.traps));
        for i=1:num_traps
            % trap vol (of a single trap)
            tc = sum(strap_vol(ta.traps == i)); % m3
            plotCellData(Gt, tc*ones(Gt.cells.num,1), ta.trap_regions == i, 'EdgeAlpha',0.1)
            %pause(1)
        end
        title('capacity per trap region [m3 co2]'); colorbar; axis equal tight off

        subplot(1,3,2); plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1);
        num_traps = max(unique(ta.traps));
        for i=1:num_traps
            % curr vol per trap region
            plotCellData(Gt, curr_vol_per_trap_region{i}*ones(Gt.cells.num,1), ...
                ta.trap_regions == i, 'EdgeAlpha',0.1)
            %pause(1)
        end
        title('curr vol per trap region [m3 co2]'); colorbar; axis equal tight off
        
        subplot(1,3,3); plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1);
        num_traps = max(unique(ta.traps));
        for i=1:num_traps
            % curr vol per trap region
            plotCellData(Gt, ult_vol_per_trap_region{i}*ones(Gt.cells.num,1), ...
                ta.trap_regions == i, 'EdgeAlpha',0.1)
            %pause(1)
        end
        title('final vol per trap region [m3 co2]'); colorbar; axis equal tight off
    
    end
end

