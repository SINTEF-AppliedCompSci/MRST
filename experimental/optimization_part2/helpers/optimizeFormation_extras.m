function [Gt, optim, init, history, other] = optimizeFormation_extras(varargin)

% extra parts added to handle other formations (without subscale trapping
% files), pre-defined well locations, initial rates, schedule, etc.

   moduleCheck('ad-core');
   
   gravity on;
   opt = merge_options(opt_defaults(), varargin{:});
   
   %% Physical grid and rock
   if ~strcmpi(opt.modelname, 'Synthetic')
        [Gt, rock2D, ~] = getFormationTopGrid(opt.modelname, opt.coarse_level);
   else
        Gt          = dipped_perturbed_grid('originalGrid',false, 'dz_e',1.05);
        rock2D.poro = NaN * ones(Gt.cells.num,1);
        rock2D.perm = NaN * ones(Gt.cells.num,1);
   end

   if any(isnan(rock2D.perm))
        rock2D.perm = 500*milli*darcy * ones(Gt.cells.num,1);
   end
   if any(isnan(rock2D.poro))
        rock2D.poro = 0.25 * ones(Gt.cells.num,1); 
   end
   % @@ however, we do not want to overwrite all perm, poro, only the NaN
   % entries...
   
   %% Spill-point analysis object
   ta = trapAnalysis(Gt, opt.trap_method);
   
   %% CO2 property object
   co2 = CO2props(); % used in pick well sites
   
   %% Load subscale trapping function, if present
   dh = [];
   if ~isempty(opt.trapfile_name)
      if ~exist(opt.trapfile_name, 'file')
         warning(['Subtrap datafile not found.  Proceeding without subtraps.']);
         if ~strcmpi(opt.surf_topo,'smooth')
            warning('surf_topo option is being changed to ''smooth''.')
            opt.surf_topo = 'smooth';
         end
      else
         dh = computeToptraps(load(opt.trapfile_name), Gt, true);
      end
   elseif ~strcmpi(opt.surf_topo,'smooth')
         warning('surf_topo option is being changed to ''smooth''.')
         opt.surf_topo = 'smooth';
   end

   %% Defining fluid
   % computing fixed temperature field
   % (NB: when ref_temp and ref_depth are the quantities at the seafloor,
   % T_ref is the caprock temperature.)
   T_ref = opt.ref_temp + opt.temp_grad * (Gt.cells.z - opt.ref_depth) / 1000;   
   
   fluid = makeVEFluid(Gt, rock2D, 'sharp interface'                   , ...
                       'fixedT'       , T_ref                          , ...
                       'residual'     , [opt.sw        , opt.sr]       , ...
                       'wat_rho_ref'  , opt.rhoW                       , ...
                       'co2_rho_ref'  , opt.refRhoCO2                  , ... % does this have to relate to mean(caprock_pressure) and T_ref?
                       'wat_rho_pvt'  , [opt.c_water   , opt.refPress] , ...
                       'co2_rho_pvt'  , [opt.p_range   , opt.t_range]  , ...
                       'wat_mu_ref'   , opt.muBrine                    , ...
                       'co2_mu_pvt'   , [opt.p_range   , opt.t_range]  , ...
                       'pvMult_fac'   , opt.pvMult                     , ...
                       'pvMult_p_ref' , opt.refPress                   , ... % mean(caprock_pressure)
                       'dissolution'  , opt.dissolution                , ...
                       'dis_rate'     , opt.dis_rate                   , ...
                       'dis_max'      , opt.dis_max                    , ...
                       'surf_topo'    , opt.surf_topo                  , ...
                       'top_trap'     , dh);
                   
   %% Defining initial state
   initState.pressure = compute_hydrostatic_pressure(Gt, fluid.rhoWS, opt.surface_pressure);
   initState.s = repmat([1 0], Gt.cells.num, 1);
   initState.sGmax = initState.s(:,2);
   
   if opt.dissolution
        initState.rs = zeros(Gt.cells.num, 1);
   end
   
   %% Compute overburden pressure (used when pressure is penalized)
   [P_over, ~] = computeOverburdenPressure(Gt, rock2D, opt.ref_depth, fluid.rhoWS);
   
   %% Compute storage efficiency for closed system (optional)
   % E is computed if system is closed, and is used to ensure the injection
   % masses are not larger than the formation's closed-system capacity,
   % given by:    M_co2 = E_closed * rhoCO2(P,T) * pore_volume
   E = [];
   if strcmpi(opt.btype,'flux') && opt.adjustClosedSystemRates
        E  = closed_system_storage_efficiency( opt.c_water, opt.pvMult, P_over, initState.pressure );
   end
   
    %% Construct schedule
    use_default_schedule = isempty(opt.schedule);
    if use_default_schedule
       
        %%% 1) Place Wells: (qt is in kg)----------------------------------
        if strcmpi(opt.well_placement_type,'use_array')

            [wc, qt] = pick_wellsites_array(Gt, rock2D, co2, ta, opt.rhoW, ...
                        opt.sw, opt.ref_temp, opt.ref_depth, opt.temp_grad, ...
                        opt.well_buffer_dist_domain, opt.well_buffer_dist_catchment, ...
                        opt.DX, opt.DY, opt.max_num_wells, 'inspectWellPlacement',false, ...
                        'E_closed', E);

        elseif strcmpi(opt.well_placement_type,'one_per_trap')

            [wc, qt] = pick_wellsites_onePerTrapRegion(Gt, rock2D, co2, ta, ...
                        opt.rhoW, opt.sw, opt.ref_temp, opt.ref_depth, ...
                        opt.temp_grad, opt.max_num_wells, opt.pick_highest_pt, ...
                        opt.well_buffer_dist, opt.well_buffer_dist_domain, ...
                        opt.well_buffer_dist_catchment, 'inspectWellPlacement',false, ...
                        'E_closed', E);

        elseif strcmpi(opt.well_placement_type,'one_per_path')

            [wc, qt] = pick_wellsites_test(Gt, rock2D, co2, ta, opt.max_num_wells, opt.rhoW, ...
                        opt.sw, opt.ref_temp, opt.ref_depth, opt.temp_grad, ...
                        opt.well_buffer_dist, opt.maximise_boundary_distance, ...
                        opt.well_buffer_dist_domain, opt.pick_highest_pt, ...
                        opt.well_buffer_dist_catchment, ...
                        'inspectWellPlacement',false, 'E_closed', E);
        end
        assert(~isempty(wc), 'No well was placed.')
        
        % Discard wells (optional):
        % @@ should injection masses, qt, be recalculated?
        if ~isempty(opt.well_nums_to_keep)
            wc = wc(opt.well_nums_to_keep);
            qt = qt(opt.well_nums_to_keep);
        end
        
        %qt = qt * 2;
        %qt(1:2) = qt(1:2) * 2.5;% @@
        %qt(3) = qt(3) * 0.5;
        % make a call to a function that adjusts initial rates such that
        % savings (assuming a certain amount of leakage) is positive
        %perc_leakage = 0*ones(numel(wc),1);
        %[M_crit, M_crit_l, M_crit_u] = critical_injection_masses( opt.well_initial_cost, ...
        %    opt.well_operation_cost, opt.leak_penalty*1, ...
        %    opt.co2_tax_credit, 'perc_leakage', perc_leakage );
        %qt(1:2) = M_crit_l;
        %qt(3) = M_crit_l;
        %M_crit(M_crit < 0) = M_crit_u;
        %qt(1:end) = M_crit;
        
        %qt(4) = M_crit_l;
        %qt(1:2) = M_crit_l;
        %qt(3) = 0;
        
        % Set qt to critical mass required for well to be worthwhile (optional):
        if opt.useBreakEvenRate
            clear qt
            qt = ones(numel(wc),1) * opt.well_initial_cost / (opt.co2_tax_credit - opt.well_operation_cost); % tonnes per well
            qt = qt * 1e3; % kg per well
        end
        
      
        %%% 2) Create schedule based on control type: -----------------------
        if strcmpi(opt.well_control_type,'rate')
            opt.schedule = setSchedule_extras( Gt, rock2D, wc, 'rate', ...
                                opt.isteps, opt.itime, opt.msteps, opt.mtime, ...
                                'wqtots', qt/opt.refRhoCO2, 'minval',sqrt(eps));

        elseif strcmpi(opt.well_control_type,'bhp')
            opt.schedule = setSchedule_extras( Gt, rock2D, wc, 'bhp', ...
                                opt.isteps, opt.itime, opt.msteps, opt.mtime, ...
                                'initState', initState); % @@ test out this call

        end
      
        %%% 3) Inspect well placement and schedule (optional): --------------
        if opt.inspectWellPlacement
            cinx_inj = [opt.schedule.control(1).W.cells];

            %close all
            figure(100); set(gcf,'Position',[2929 666 953 615])
            subplot(2,2,[1 3])
            mapPlot(gcf, Gt, 'traps', ta.traps, ...
                'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
                'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
                'maplines', 20, 'wellcells', cinx_inj, 'well_numbering', true);
            colorizeCatchmentRegions(Gt, ta);
            plotGrid(Gt, 'FaceColor','none','EdgeAlpha',0.1)
            axis equal tight off

            subplot(2,2,2); bar(qt);
            xlabel('well number', 'FontSize',16);
            ylabel('mass to inject [kg]', 'FontSize',16); % i.e., mass capacity of empty traps along spill path

            rates = [opt.schedule.control(1).W.val].*opt.refRhoCO2.*(1*year)./10^9; % Mt/yr
            subplot(2,2,4); bar(rates);
            xlabel('well number', 'FontSize',16);
            ylabel({'initial rate [Mt/yr]';['for ',num2str(convertTo(opt.itime,year)),' yrs inj.']}, 'FontSize',16);

            % or in 3D view:
            figure(101);
            plotGrid(Gt.parent, 'facecolor','none', 'edgealpha',0.1); view(3); daspect([1,1,0.02])
            plotCellData(Gt, Gt.cells.z, 'edgealpha',0.5);
            plotWell(Gt.parent, opt.schedule.control(1).W)
            
            % Exit before entering optimization routine
            init = []; optim = []; history = [];
            other.inspectWellPlacement = opt.inspectWellPlacement;
            return
        end

    else
        fprintf('\n Using schedule passed in as varargin. \n')
        if isa(opt.schedule,'char')
            tmp = load(opt.schedule);
            if isfield(tmp,'optim') % extend to handle init.schedule, etc.
                tmp = tmp.optim;
            end
            name = fields(tmp);
            opt.schedule = tmp.(name{1}); % add check to ensure wells are located in grid
            clear tmp
        elseif isa(opt.schedule,'struct')
            % nothing required
        end
        % ensure migration rates are the minimum (as they may have changed
        % slightly when optimal rates were obtained previously).
        for i=1:numel([opt.schedule.control(2).W.val])
           opt.schedule.control(2).W(i).val = sqrt(eps); 
        end
        % ensure injection rates are equal or greater than the minimum (as
        % they may have changed slightly as above)
        for i=1:numel([opt.schedule.control(1).W.val])
           opt.schedule.control(1).W(i).val = max(sqrt(eps),opt.schedule.control(1).W(i).val); 
        end
        %
        % Extend migration period (optional)
        if ~isempty(opt.extended_mig_time)
           opt.schedule = extend_migration_sch( opt.schedule, opt.extended_mig_time ); 
           opt.mtime = opt.extended_mig_time;
        end
        %
        % Set rates during migration to be zero
        if opt.set_mig_rates_to_0
            for i=1:numel([opt.schedule.control(2).W.val])
                opt.schedule.control(2).W(i).val = 0; 
            end
        end
    end
   
    
    %% Add boundary conditions (either constant pressure, or no-flow)
    bfaces = identifyBoundaryFaces(Gt);
    if strcmpi(opt.btype,'pressure')
        bdryVal     = Gt.faces.z(bfaces) * fluid.rhoWS * norm(gravity) + opt.surface_pressure;
    elseif strcmpi(opt.btype,'flux')
        bdryVal     = zeros(numel(bfaces),1);
    end
    bc = addBC( [], bfaces, opt.btype, bdryVal, 'sat', [1 0] );
    for i = 1:numel(opt.schedule.control)
        opt.schedule.control(i).bc = bc;
    end
   
    %% Assign well control limits:
    % either set here internally, or use a pre-defined rate limit
    P_limit = [];
    
    % for injection period
    if strcmpi(opt.well_control_type,'rate')
        min_wvals = sqrt(eps) * ones(numel(opt.schedule.control(1).W), 1);
        if isempty(opt.max_wvals)
            max_wvals = opt.rate_lim_fac * ...
                max([opt.schedule.control(1).W.val]) * ones(numel(opt.schedule.control(1).W), 1);
        else
            assert( numel(opt.max_wvals) == numel(opt.schedule.control(1).W) )
            max_wvals = opt.max_wvals;
        end
    elseif strcmpi(opt.well_control_type,'bhp')
        min_wvals = initState.pressure( [opt.schedule.control(1).W.cells] ); % @@ but will this cause small injection to occur as global formation pressure rises?

        % the max well values are set to be 90% of the overburden pressure
        % at the well cell location
        [P_over, ~] = computeOverburdenPressure(Gt, rock2D, opt.ref_depth, fluid.rhoWS);
        P_limit     = P_over * 0.9;
        max_wvals   = P_limit([opt.schedule.control(1).W.cells]); 
    end
    
    if strcmpi(opt.penalize_type,'pressure')
        [P_over, ~] = computeOverburdenPressure(Gt, rock2D, opt.ref_depth, fluid.rhoWS);
        P_limit     = P_over * opt.p_lim_factor;
    end
    
    % for migration period
    min_mig_rates = sqrt(eps) * ones(numel(opt.schedule.control(1).W), 1);
    
    % checks:
    assert( all(max_wvals >= [opt.schedule.control(1).W.val]') , ...
            'An initial well value exceeds the maximum constraint.')
    assert( all(min_wvals <= [opt.schedule.control(1).W.val]') , ...
            'An initial well value is less than the minimum constraint.')
    if ~opt.set_mig_rates_to_0
        assert( all(min_mig_rates <= [opt.schedule.control(2).W.val]') , ...
            'An initial migration rate is less than the minimum rate constraint.')
    end

   %% Set up model and run optimization
   model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);

   [optim, init, history] = ...
       optimizeRates_extra(initState, model, opt.schedule, ...
                     min_wvals, max_wvals, min_mig_rates, ...
                     'dryrun',              opt.dryrun, ...
                     'last_control_is_migration', true, ...
                     'penalize_type',       opt.penalize_type, ...
                     'leak_penalty',        opt.leak_penalty, ...
                     'pressure_penalty',    opt.pressure_penalty, ...
                     'pressure_limit',      P_limit, ...
                     'account_well_cost',   opt.account_well_cost, ...
                     'well_initial_cost',   opt.well_initial_cost, ...
                     'well_operation_cost', opt.well_operation_cost, ...
                     'co2_tax_credit',      opt.co2_tax_credit, ...
                     'nonlinear_well_cost', opt.nonlinear_well_cost, ...
                     'alpha',               opt.alpha, ...
                     'rho_water',           fluid.rhoWS, ...
                     'surface_pressure',    opt.surface_pressure, ...
                     'lineSearchMaxIt',     opt.lineSearchMaxIt, ...
                     'gradTol',             opt.gradTol, ...
                     'objChangeTol',        opt.objChangeTol, ...
                     'obj_scaling',         []);
                 
   %% Store data
   %other.fluid = fluid; % @@ fluid structure is ~400MB when saved, but required for results
   %other.fluid.rhoGS = fluid.rhoGS;
   other.rock = rock2D;
   other.residual = [opt.sw, opt.sr];
   other.traps = ta;
   other.dh = dh;
   other.initState = initState;
   other.opt = opt;
   other.dissolution = opt.dissolution;
   other.P_over = P_over;
   other.E_closed = E;
   other.adjustClosedSystemRates = opt.adjustClosedSystemRates;
   
end

% ----------------------------------------------------------------------------
function bfaces = identifyBoundaryFaces(Gt)
   
    bfaces = find(any(Gt.faces.neighbors==0,2)); 

end

% ----------------------------------------------------------------------------
function wc = select_wellcell(Gt, candidates, buffer, domain_buffer)

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
   
   inside_buffer_candidates = interior_candidates(distances > buffer);
   
   % Additional check to ensure inside_buffer_candidates are within a
   % buffer distance to formation boundary.
   bdist = sqr_distance_from_boundary(Gt);
   bdist = sqrt(bdist(inside_buffer_candidates));
   inside_buffer_candidates = inside_buffer_candidates(bdist > domain_buffer);
   % However, if inside_buffer_candidates are empty, the highest point in
   % the spill region is selected. This point may not satisfy this domain
   % buffer constaint. Instead, an empty wc should be passed out, and the
   % next best leaf-node is assessed.

   if isempty(inside_buffer_candidates)
      % none of the candidates were far enough from boundary.  We will just
      % choose the highest point
      candidate_z = Gt.cells.z(candidates);
      [~, min_ix] = min(candidate_z);
      wc = candidates(min_ix);
   else
   
      inside_buffer_z = Gt.cells.z(inside_buffer_candidates);
      [~, max_z_ix] = max(inside_buffer_z);
      
      wc = inside_buffer_candidates(max_z_ix);
   end
end

% ----------------------------------------------------------------------------
function [wc, qt] = pick_wellsites(Gt, rock2D, co2, ta, num_wells, ...
                                   rhoW, sw, seafloor_temp, seafloor_depth, ...
                                   tgrad, buf_dist, maximize_boundary_distance, ...
                                   domain_buf_dist)

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
           % Pick well site as far downslope as possible (while keeping a
           % buffer distance from boundary
           
           qt(i) = max(field); 
           candidate_cells = find(field == qt(i));
           wc(i) = select_wellcell(Gt, candidate_cells, buf_dist, domain_buf_dist);
        end
        
        % Setting available capacity of used traps to zero
        trapcap = set_used_to_zero(trapcap, ta.trap_adj, ta.trap_regions(wc(i)));
        
        % check if there is no more avaible capacity for the wells (which
        % can occur if all traps have been used)
        if all(trapcap == 0) || qt(i)/qt(1) < 0.01
            if all(trapcap == 0)
                warning(['You wanted to place %d wells, '...
                'but all traps were used after placing %d wells.'], num_wells, i)
            elseif qt(i)/qt(1) < 0.01
                warning(['You wanted to place %d wells, '...
                'but the injected volume of well number %d was '...
                'less than 1 percent of well number 1''s volume'], num_wells, i)
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

function opt = opt_defaults()

    opt.modelname = 'utsirafm';
    opt.schedule = [];
    opt.extended_mig_time = [];
    opt.set_mig_rates_to_0 = false;
    opt.coarse_level = 3;
    
    % Trapping Structure method:
    opt.trap_method = false; % false for node-based, true for cell-based
    
    % Number of timesteps (injection and migration)
    opt.isteps = 10;
    opt.msteps = 31;

    % Durations of injection an migration phases
    opt.itime  = 50   * year;
    opt.mtime  = 3000 * year;
    
    opt.surface_pressure = 1 * atm;
    
    % Default fluid properties
    [rho_default, mu_default, sr, sw] = getValuesSPE134891();
    opt.refRhoCO2 = rho_default(2);
    opt.rhoW      = rho_default(1);
    opt.muBrine   = mu_default(1);
    opt.muCO2     = 0; % zero value will trigger use of variable viscosity 
    opt.pvMult    = 1e-5/barsa; % pore volume (rock) compressibility
    opt.refPress  = 100 * barsa;
    opt.c_water   = 4.3e-5/barsa; % water compressibility
    opt.p_range   = [0.1, 400] * mega * Pascal; % pressure span of sampled property table
    opt.t_range   = [4 250] + 274; % temperature span of sampled property table 

    % Residual saturations
    opt.sr = sr; % gas
    opt.sw = sw; % brine

    % Temperature regime control parameters
    opt.ref_temp = 7 + 273.15;
    opt.ref_depth = 80; 
    opt.temp_grad = 35.6;

    % Dissolution parameters (default zero to trigger instantaneous
    % dissolution, and suggested nonzero values in comments)
    opt.dissolution = false;
    opt.dis_rate    = 0;    % 0.44 * kilogram / rho / poro / (meter^2) / year = 8.6924e-11;
    opt.dis_max     = 0;    % 53/760 = 0.07

    % Subtrap file
    opt.trapfile_name = '';      % 'utsira_subtrap_function_3'
    opt.surf_topo = 'inf_rough'; % will be set to 'smooth' if trapfile_name is empty
    
    % Directory for saving reslts
    opt.report_basedir = './simulateUtsira_results/';
    
    % Control whether to enter optimization
    opt.dryrun = false;

    % Penalization details:
    opt.penalize_type = 'leakage'; % 'leakage','leakage_at_infinity','pressure'
    opt.leak_penalty  = 10;
    opt.pressure_penalty = [];
    opt.p_lim_factor = 0.9; % factor to apply to P_overburden
    
    % Accounting for well cost:
    opt.account_well_cost = false;
    opt.well_initial_cost = [];
    opt.well_operation_cost = [];
    opt.co2_tax_credit = [];
    opt.nonlinear_well_cost = true;    % if true, opt.alpha must not be empty
    opt.alpha = [];
   
    
    % Boundary type:
    opt.btype = 'pressure'; % 'pressure' or 'flux'
   
    % Well control details:
    opt.well_control_type = 'rate'; % 'rate','bhp'
    opt.rate_lim_fac = 10; % factor for setting upper limit of injection rates
    opt.max_wvals = [];    % or max well vals can be passed in explicitly
    
    % Well placement details:
    opt.well_placement_type = 'use_array'; % 'use_array', 'one_per_trap', 'one_per_path'
    opt.max_num_wells = 40;
    opt.well_nums_to_keep = []; % ignored if empty
    opt.maximise_boundary_distance = false;
    opt.well_buffer_dist = 1 * kilo*meter; % dist from edge of internal catchment
    opt.well_buffer_dist_domain = 5 * kilo*meter;    % dist from edge of domain 
    opt.well_buffer_dist_catchment = 2 * kilo*meter; % dist from edge of external catchment
    opt.pick_highest_pt = false; % otherwise farthest downslope
    opt.DX = 1 * kilo*meter;
    opt.DY = 1 * kilo*meter;
    opt.inspectWellPlacement = false;
    opt.adjustClosedSystemRates = false;
    opt.useBreakEvenRate = false;
    
    % Convergence details:
    opt.lineSearchMaxIt = 10;
    opt.gradTol         = 1e-4;
    opt.objChangeTol    = 1e-4;


end

% ----------------------------------------------------------------------------

function p = compute_hydrostatic_pressure(Gt, rho_water, surface_pressure)
    gravity on;
    p = rho_water * norm(gravity()) * Gt.cells.z + surface_pressure;
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