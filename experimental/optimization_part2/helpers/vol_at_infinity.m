function [ ult_vol_remaining, ult_vol_leaked ] = vol_at_infinity( Gt, rock2D, sG, fluid, p_curr, p_future, varargin )
% Determine co2 vol (at ref. depth) that will ultimately remain in formation at time infinity.
%
% @@ test impact of using cell-based versus node-based trapAnalysis.
%
%
% SYNOPSIS:
%   future_vol = vol_at_infinity(Gt, rock2D, sG, fluid, p_curr, p_future)
%   future_vol = vol_at_infinity(Gt, rock2D, sG, fluid, p_curr, p_future, ...
%                               'ta', trapAnalysis(Gt,false), 'plotsOn', true)
%
% DESCRIPTION:
%   This calculation of co2 vol at time infinity is based on the
%   formations' co2 vol at a point in time in which it is safe to assume
%   the flow dynamics are gravity-dominated. Thus, the co2 vol in a
%   particular spill tree at time infinity is based on that spill tree's
%   capacity, and any extra co2 vol is assumed to have been leaked by time
%   infinity.
%
%   NB: this routine was developed to handle ADI variables (sG, p_curr,
%   ult_vol_remaining), thus the syntax used here respects ADI variables,
%   and in some instances a cell array of ADI variables.
%
%   The ult_vol_remaining is a volume of co2 at reference conditions. It
%   may be converted into a co2 mass simply by multiplication with the co2
%   density at reference conditions (i.e., fluid.rhoGS)

% INPUTS    - Gt, top surface grid
%           - rock2D, includes porosity (and possibly ntg) data
%           - sG, saturation of co2, i.e., states.s(:,2)
%           - fluid, structure containing:
%               - res_water, res_gas, rhoGS, bG(p), pvMultR(p)
%               - (NB: fluid is assumed to be same fluid structure used in
%                 obtaining sG by simulation)
%
% (optional) - plotsOn, true or false for plotting
%            - ta, as given by trapAnalysis. If not passed in, will be
%              computed. NB: used for traps, trap_regions, spill-depth @@

% RETURNS   - co2 volumes in terms of m3 at ref. depth:
%               1. ult_vol_remaining, amount remaining at time infinity
%               2. ult_vol_leaked, amount leaked at time infinity


    opt.plotsOn = false;
    opt.ta = [];
    opt.time = [];
    opt.trap_method = false;
    opt = merge_options(opt, varargin{:});

    if isempty(opt.ta)
        fprintf('Obtaining trap structure using trapAnalysis...\n')
        fprintf('Method is (1=cell-based, 0=node-based): %1.0f \n', opt.trap_method)
        ta = trapAnalysis(Gt, opt.trap_method); % @@ implement option to pass in closed boundary faces.
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
    
    sw = fluid.res_water;
    sr = fluid.res_gas;
    

    % To handle fluid vols, we convert co2 vol in reservoir to a volume at
    % the reference depth. We distinguish between current and future
    % pressure states since the current vols exist at the current pressure,
    % while the capacity of the domain is in terms of hydrostatic (future)
    % pressure.
    pvMult_curr   = fluid.pvMultR(p_curr); % can be ADI
    pvMult_future = fluid.pvMultR(p_future);
    bG_curr       = fluid.bG(p_curr); % can be ADI
    bG_future     = fluid.bG(p_future);

    cells = 1:Gt.cells.num;
    
    
    % 1.a) Get capacity of structural traps (in m3 co2 at ref. depth)
    
        % computing structual trap heighs (H1) for each cell
        tcells               = find(ta.traps);
        trap_heights         = zeros(Gt.cells.num, 1);
        trap_heights(tcells) = ta.trap_z(ta.traps(tcells)) - Gt.cells.z(tcells);
        trap_heights         = min(trap_heights, Gt.cells.H);
        H2                   = Gt.cells.H - trap_heights;
        assert(all(trap_heights <= Gt.cells.H));
    
        % computing structural trap pore vol (m3)
        strap_pvol           = Gt.cells.volumes .* trap_heights .* poro .* ntg .* pvMult_future;

        % computing free co2 capacity in structural traps (m3 co2, @ ref. depth)
        strap_co2_vol        = strap_pvol .* (1-sw) .* bG_future;               % per cell
        num_traps            = max(unique(ta.traps));
        trapcap              = cell(num_traps,1);
        for i = 1:num_traps
            trapcap{i}       = sum(strap_co2_vol(ta.traps == i));               % per trap
        end
        dispif(mrstVerbose, ...
            'The total structural trapping capacity is %5.3f Mt co2. \n\n', ...
            sum(cell2mat(trapcap))*fluid.rhoGS/1e9);
        
        
    % 1.b) Get maximum residual capacity in formation outside of straps
    
        % computing pore vol outside of straps (m3)
        btrap_pvol         = Gt.cells.volumes .* H2 .* poro .* ntg .* pvMult_future;
        
        % computing capacity for residual co2 (m3 at ref. depth)
        btrap_res_co2_vol  = btrap_pvol .* sr .* bG_future;
        
        dispif(mrstVerbose, ...
            'The total residual trapping capacity is %5.3f Mt co2. \n\n', ...
            sum(btrap_res_co2_vol)*fluid.rhoGS/1e9);
        

    % 2. Get current co2 vol in each trap region (m3 at ref. depth)
    
        % computing distribution of current co2 vol (m3 at ref. depth)
        pvol = Gt.cells.volumes .* Gt.cells.H .* poro .* pvMult_curr .* ntg; % can be ADI
        curr_vol = ( sG .* bG_curr .* pvol );                                % can be ADI
    
        % computing current vol per trap region (m3 at ref. depth)
        assert( max(unique(ta.traps)) == max(unique(ta.trap_regions)) )
        curr_vol_per_trap_region = cell(num_traps,1);
        for i = 1:num_traps
            curr_vol_per_trap_region{i} = sum(curr_vol(ta.trap_regions == i)); % can be ADI
        end
    
    
    % 3. Get the trap index of the root trap of each tree
    
        trees = maximizeTrapping(Gt, 'res', ta); % 'calculateAll', 'removeOverlap' ?
        % NB: sth buggy in the order of traps for Sto's first spill-pathway, so
        % the trap indexes of traps along a spill-path (i.e., what trees.traps
        % should have contained) are obtained using another function:
        treeRoots = [trees.root];
    
    
    % 4. Get the vol that is ultimately leaked per tree:
    
        % We do this by spilling any vols a trap cannot retain down the
        % tree, and out of the domain. During this calculation, we also
        % determine the co2 vol ultimately remaining in the domain
    
        vol_leaked_per_tree = cell(numel(treeRoots),1);
        ult_vol_per_trap_region = curr_vol_per_trap_region;                 % can be ADI
        clear i
        dispif(mrstVerbose, 'There are %d trees. \n\n', numel(treeRoots));
        for i = 1:numel(treeRoots)
            dispif(mrstVerbose, '\nLooking at tree %d... \n', i);

            % Get the trap indexes of all traps downstream from a tree
            % root, including the trap index of the root itself. (NB: While
            % trap_vols is not required, we use its value to compare
            % against the trapcap computed above. We pass in poro * ntg *
            % pvMult_future as the porosity, to account for the impact of
            % pressure and possible ntg data.)
            [trap_ixs, trap_vols] = downstreamTraps(Gt, ...
                poro.*ntg.*pvMult_future.*bG_future, ta, treeRoots(i));
            
            % Compare trap_vols (m3 pore space of traps) with trapcap_vol
            % (m3 reservoir co2). NB: we must account for residual water (sw).
            trap_co2_vols = trap_vols' .* (1-sw);
            trapcap_tmp = cellfun( @(x) double(x), {trapcap{trap_ixs}})';
            assert( all( abs( trap_co2_vols - trapcap_tmp )./abs(trap_co2_vols) < 1e-8 ) ) % syntax handles cell array trapcap
        
            % Spill vols along a tree...
            dispif(mrstVerbose, 'The trap indexes in this tree are: %s \n\n', num2str(trap_ixs));
            for j = 1:numel(trap_ixs)
                dispif(mrstVerbose, 'Looking at trap %d... \n', trap_ixs(j));
                dispif(mrstVerbose, 'Trap %d has a capacity of %d (m3 co2). \n', ...
                    trap_ixs(j), double(trapcap{trap_ixs(j)}) );
                dispif(mrstVerbose,'Trap %d currently has %d (m3 co2). \n', ...
                trap_ixs(j), double(ult_vol_per_trap_region{trap_ixs(j)}) );
            
            
                % Computing the vol that will spill out of trap_ixs(j)
                spill_vol = max(0, - trapcap{trap_ixs(j)} + ult_vol_per_trap_region{trap_ixs(j)}); % may be an ADI var @@ is taking max of ADI ok?
                assert( spill_vol >= 0 )
            
                % and removing this vol from trap_ixs(j)
                ult_vol_per_trap_region{trap_ixs(j)} = ult_vol_per_trap_region{trap_ixs(j)} - spill_vol; % may be an ADI var
            
                % and spilling this "spill vol" into the next trap
                % downstream or out of the tree
                if j < numel(trap_ixs)
                    % haven't reach end of spill-tree yet, so add this
                    % spill volume to the next trap down the spill-tree
                    % (i.e., trap_ixs(j+1))
                    dispif(mrstVerbose, 'Thus %d (m3 co2) spills out of trap %d, and into trap %d.\n', ...
                        double(spill_vol), trap_ixs(j), trap_ixs(j+1));
                    ult_vol_per_trap_region{trap_ixs(j+1)} = ult_vol_per_trap_region{trap_ixs(j+1)} + spill_vol; % may be an ADI var
                    dispif(mrstVerbose, 'Trap %d now has %d (m3 co2). \n', ...
                        trap_ixs(j+1), double(ult_vol_per_trap_region{trap_ixs(j+1)}) );
                else
                    % we've reached the end of the spill-tree, thus the
                    % spill_vol goes out of the spill-tree (i.e., out of
                    % the domain)
                    dispif(mrstVerbose, 'Thus %d (m3 co2) spills out of trap %d, and out of tree %d.\n', ...
                        double(spill_vol), trap_ixs(j), i);
                    vol_leaked_per_tree{i} = spill_vol;                     % may be an ADI var
                
                end
                % NB: ult_vol_per_trap_region gets updated each time a
                % spill_vol is removed from a trap and added to the next
                % (downsteam) trap. After looping through all tree roots,
                % the final array should represent the amount of co2 per
                % trap region at stable conditions. No more migration
                % should occur under gravity-dominated dynamics.
            end
        end
    
        % co2 ultimately remaining in the domain (m3 at ref. depth)
        ult_vol_remaining = ult_vol_per_trap_region{1};
        for i = 2:numel(ult_vol_per_trap_region)
            ult_vol_remaining = ult_vol_remaining + ult_vol_per_trap_region{i}; % can be ADI
        end
    
    
    % 5. Computing mass outside of catchment regions destined to leak:
    
        % In addition to the vols that will spill-out of each tree, any co2
        % vol that existed outside of the trap regions making up the spill
        % trees is bound to leak from the domain. (NB: residually trapped
        % co2 that will remain is accounted for in another step)
        total_trap_region_vol = curr_vol_per_trap_region{1};                % can be ADI
        for i = 2:numel(curr_vol_per_trap_region)
            total_trap_region_vol = total_trap_region_vol + curr_vol_per_trap_region{i};
        end
        curr_vol_outside_trap_regions = sum(curr_vol) - total_trap_region_vol; % can be ADI
        %resid_trap_vol_outside_trap_regions = curr_vol_outside_trap_regions * sr;
        
        
    % 6. Computing vol co2 destined to be residually trapped:
    
        % A migrating plume leaves behind an amount of residually trapped
        % co2. While we do not predict the volume that the free plume
        % migrates through, we can compute a minimum amount of co2 destined
        % to be residually trapped. That is, the curr vol of co2 lying
        % outside of any structural traps is bound to migrate either into a
        % trap or out of the domain, and thus will leave behind residually
        % trapped co2.
        
        % computing co2 heights under caprock
        % (NB: we do not consider any previous residual trapping thus set
        % the historical max saturation equal to the current saturation)
        sGmax       = sG;
        [h, h_max] = upscaledSat2height(sG, sGmax, Gt, 'resSat', [sw sr]);  % may be ADI
        assert( all(h <= Gt.cells.H) )
        assert( all(h_max <= Gt.cells.H) )
        
        % computing trap's spill point depth, zt
        % NB: same as trap_heights
        tr   = ta.trap_regions;
        tr_z = [0; ta.trap_z];
        tr   = 1 + tr;
        tr(tr>numel(tr_z)) = 1;
        zt   = max(tr_z(tr) - Gt.cells.z, 0);
        zt   = min(zt, Gt.cells.H);
        
        % computing vol of free co2 plume (m3 at ref. depth) given future
        % pressure
        pv_area  = Gt.cells.volumes .* poro .* ntg .* pvMult_future;     % effective pore area, per cell
        strucVol = sum( min(zt, h) .* bG_future .* pv_area );               % may be ADI. if h is ADI, ok to take min??
        plumeVol = sum( h .* bG_future .* pv_area ) - strucVol;             % may be ADI. if h is ADI, can be ADI
        
        % computing amount of residually trapped co2 within the free plume
        % that is destined to migrate away
        resid_trap_vol  = plumeVol * sr;                                    % may be ADI
        
    
    % 7. Computing the final total amount of co2 that will ultimately leak:
    
        vol_leaked_total_trees = vol_leaked_per_tree{1};                    % may be an ADI
        for i = 2:numel(vol_leaked_per_tree)
            vol_leaked_total_trees = vol_leaked_total_trees + vol_leaked_per_tree{i};
        end
        ult_vol_leaked = vol_leaked_total_trees + curr_vol_outside_trap_regions ...
            - resid_trap_vol;                                               % m3, may be an ADI var
        %assert( ult_vol_leaked >= 0 )
        assert( ult_vol_leaked > -1e-5 ) % to account for possible numerical round-off error due to spilling of vols
        ult_vol_remaining = sum(curr_vol) - ult_vol_leaked;                 % may be an ADI

        
    % 8. Check that vols balance out, within some tolerance
    
        % these should balance out to machine precision
        assert( abs( sum(curr_vol)-(ult_vol_remaining + ult_vol_leaked) ) ./ abs(sum(curr_vol)) < eps ) %< 1e-8 )   
        
        
        tol = 1e-8;
        % ultimate volume remaining cannot be greater than the maximum
        % retaining capacity of structural + residual trapping
        max_cap = sum(strap_co2_vol) + sum(btrap_res_co2_vol);
        assert( ult_vol_remaining <= max_cap + max_cap*tol )
        
        % ultimate volume remaining in straps cannot be greater than the
        % total structural capacity
        assert( ult_vol_remaining - resid_trap_vol <= sum(strap_co2_vol) + sum(strap_co2_vol)*tol )
        
        % ultimate volume remaining as residually trapped cannot be greater
        % than the total residual capacity
        assert( resid_trap_vol <= sum(btrap_res_co2_vol) + sum(btrap_res_co2_vol)*tol )
        
        
    % 9. Print out summary
    
        if opt.plotsOn
        % force summary to be displayed
        switchBackAfterPlots = false;
        if ~mrstVerbose
            switchBackAfterPlots = true;
            mrstVerbose on;
        end
        
        dispif(mrstVerbose, ...
            '\nThe total structural trapping capacity is %5.3f Mt co2. \n', ...
            sum(cell2mat(trapcap))*fluid.rhoGS/1e9); 
        dispif(mrstVerbose, ...
            'The total residual trapping capacity is %5.3f Mt co2. \n\n', ...
            sum(btrap_res_co2_vol)*fluid.rhoGS/1e9);
    
        dispif(mrstVerbose, ...
            '%5.3f Mt co2 is currently within formation. \n', sum(curr_vol)*fluid.rhoGS/1e9);
        dispif(mrstVerbose, ...
            '   --- %5.3f Mt will remain: \n', ult_vol_remaining*fluid.rhoGS/1e9)
        dispif(mrstVerbose, ...
            '       --- %5.3f Mt structurally (%5.3f perc. of strap capacity) \n', ...
            (ult_vol_remaining - resid_trap_vol)*fluid.rhoGS/1e9, (ult_vol_remaining - resid_trap_vol)/sum(cell2mat(trapcap))*100);
        dispif(mrstVerbose, ...
            '       --- %5.3f Mt residually (%5.3f perc. of resid capacity) \n', ...
            resid_trap_vol*fluid.rhoGS/1e9, resid_trap_vol/sum(btrap_res_co2_vol)*100);
        dispif(mrstVerbose, ...
            '   --- %5.3f Mt will leak \n', ult_vol_leaked*fluid.rhoGS/1e9);
        
        if switchBackAfterPlots
            mrstVerbose off;
        end
        end
    
        
    % 10. Plotting (optional)
    
        % NB: these co2 vols are wrt the ref. depth
        if opt.plotsOn

            if ~ishandle(45)
                figure(45); set(gcf,'Position',[2607 182 560 420])
            else
                % Avoid stealing focus if figure already exists
                set(0, 'CurrentFigure', 45); clf(45)
            end
            plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1);
            plotCellData(Gt, curr_vol, cells(curr_vol>0.1), 'EdgeAlpha',0.1);
            title(['current (t=',num2str(opt.time),') co2 vol [m3]'])
            colorbar; axis equal tight off
            %drawnow;
            
            if ~ishandle(46)
                figure(46); set(gcf,'Position',[4163 917 560 420])
            else
                % Avoid stealing focus if figure already exists
                set(0, 'CurrentFigure', 46); clf(46)
            end
            mapPlot(gcf, Gt, 'traps', ta.traps, ...
                    'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
                    'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
                    'maplines', 20);
            colorizeCatchmentRegions(Gt, ta); % this actually colorizes regions around each trap
            axis equal tight off
            title('trap regions')
            %drawnow;
            
            % pie to show breakdown of leaked and remaining:
            if ~ishandle(48)
                figure(48); set(gcf,'Position',[2684 901 560 420])
            else
                % Avoid stealing focus if figure already exists
                set(0, 'CurrentFigure', 48); clf(48)
            end
            % @@ hack to make pie chart plot "0" values
            ph = pie(max([ult_vol_leaked, (ult_vol_remaining - resid_trap_vol), ...
                resid_trap_vol], eps));
            patches = findobj(ph,'Type','Patch');
            set(patches(1),'FaceColor','r');
            set(patches(2),'FaceColor','y');
            set(patches(3),'FaceColor','g');
            lh = legend('leaked','remaining (in straps)','remaining (residually)');
            set(lh,'Location','SouthEastOutside')
            title(['Future of ',num2str(sum(curr_vol)*fluid.rhoGS/1e9),' Mt co2'])
            %drawnow;
            
            
            if ~ishandle(47)
                figure(47); set(gcf,'Position',[3236 140 1655 529])
            else
                % Avoid stealing focus if figure already exists
                set(0, 'CurrentFigure', 47); clf(47)
            end
            subplot(1,3,1); plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1);
            num_traps = max(unique(ta.traps));
            for i=1:num_traps
                % trap vol (of a single trap)
                tc = sum(strap_co2_vol(ta.traps == i)); % m3
                plotCellData(Gt, tc*ones(Gt.cells.num,1), ta.trap_regions == i, 'EdgeAlpha',0.1)
            end
            title('capacity per trap region [m3 co2]'); colorbar; axis equal tight off

            
            subplot(1,3,2); plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1);
            num_traps = max(unique(ta.traps));
            for i=1:num_traps
                % curr vol per trap region
                plotCellData(Gt, curr_vol_per_trap_region{i}*ones(Gt.cells.num,1), ...
                    ta.trap_regions == i, 'EdgeAlpha',0.1)
            end
            title(['curr (t=',num2str(opt.time),') vol per trap region [m3 co2]']); colorbar; axis equal tight off

            
            subplot(1,3,3); plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1);
            num_traps = max(unique(ta.traps));
            for i=1:num_traps
                % curr vol per trap region
                plotCellData(Gt, ult_vol_per_trap_region{i}*ones(Gt.cells.num,1), ...
                    ta.trap_regions == i, 'EdgeAlpha',0.1)
            end
            title('final vol per trap region [m3 co2]'); colorbar; axis equal tight off
            %drawnow;
            
            % plot amount of trapping capacity that is utilized by time infinity
            if ~ishandle(49)
                figure(49); set(gcf,'Position',[3005 712 926 590])
            else
                % Avoid stealing focus if figure already exists
                set(0, 'CurrentFigure', 49); clf(49)
            end
            % first make mapPlot
            mapPlot(gcf, Gt, 'traps', ta.traps, ...
                    'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
                    'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
                    'maplines', 20);
            colorizeCatchmentRegions(Gt, ta); % this actually colorizes regions around each trap
            hold on
            % then add percentage of trap region utilized
            for i=1:num_traps
                % trap vol (of a single trap)
                tc = sum(strap_co2_vol(ta.traps == i)); % m3
                used_cap(i) = ult_vol_per_trap_region{i} / tc * 100; % percentage 
            end
            % top cell of each trap = ta.top
            % @@ ta.top produced by node-based trap analysis is wrong, so
            % use cell-based ta.top values (as these should be similar)
            ta_cell = trapAnalysis(Gt, true);
            text(Gt.cells.centroids(ta_cell.top,1), Gt.cells.centroids(ta_cell.top,2), cellstr([num2str(used_cap', '%.0f'), repmat(' ',num_traps,1), repmat(char(37),num_traps,1)]), ...
                'HorizontalAlignment','center', 'FontSize',16, 'FontWeight', 'bold')
            axis equal tight off
            %title('Percentage of structural trapping capacity utilized')
            
            %hold on
            %text(Gt.cells.centroids(ta_cell.top,1), Gt.cells.centroids(ta_cell.top,2), Gt.cells.z(ta_cell.top)-5, cellstr(num2str([1:num_traps]')), ...
            %    'HorizontalAlignment','center', 'FontSize',16, 'FontWeight', 'bold')
            %tc = sum(strap_co2_vol(ta.traps == 12)); % m3
            
            
            if ~ishandle(50)
                figure(50); set(gcf,'Position',[3005 712 926 590])
            else
                % Avoid stealing focus if figure already exists
                set(0, 'CurrentFigure', 50); clf(50)
            end
            % Remaining capacity available (in Mt)
            % vol (m3) is at ref depth
            for i=1:num_traps
                % trap vol (of a single trap)
                tc = sum(strap_co2_vol(ta.traps == i)); % m3
                avail_mass_cap(i) = (tc - ult_vol_per_trap_region{i}) * fluid.rhoGS/1e9; % Mt 
            end
            mapPlot(gcf, Gt, 'traps', ta.traps, ...
                    'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
                    'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
                    'maplines', 20);
            colorizeCatchmentRegions(Gt, ta); % this actually colorizes regions around each trap
            hold on
            
            % plot text with 3 decimal places when number is non-zero
            [x, y] = deal(Gt.cells.centroids(ta_cell.top(avail_mass_cap'>0),1), Gt.cells.centroids(ta_cell.top(avail_mass_cap'>0),2));
            data = avail_mass_cap(avail_mass_cap>0)';
            text(x, y, cellstr([num2str(data, '%.3f'), repmat(' Mt',numel(data),1)]), ...
                'HorizontalAlignment','center', 'FontSize',16, 'FontWeight', 'bold')
            
            % plot text with 0 decimal places when number is zero
            % (or simply do not plot these text values)
            %[x, y] = deal(Gt.cells.centroids(ta_cell.top(avail_mass_cap'==0),1), Gt.cells.centroids(ta_cell.top(avail_mass_cap'==0),2));
            %data = avail_mass_cap(avail_mass_cap==0)';
            %text(x, y, cellstr([num2str(data, '%.0f'), repmat(' Mt',numel(data),1)]), ...
            %    'HorizontalAlignment','center', 'FontSize',16, 'FontWeight', 'bold')
            axis equal tight off
        end
end

