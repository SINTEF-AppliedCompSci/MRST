function [Gt, optim, init, history, other] = optimizeFormation2(varargin)

% Extra capabilities added to place wells as an array within trap region(s)

   moduleCheck('ad-core');
   
   gravity on;
   opt = merge_options(opt_defaults(), varargin{:});
   
   %% Physical grid and rock
   [Gt, rock2D, ~] = getFormationTopGrid(opt.modelname, opt.coarse_level);
   
   %% Spill-point analysis objec
   ta = trapAnalysis(Gt, false);
   
   %% CO2 property object
   co2 = CO2props();
   
   %% Load subscale trapping function, if present
   dh = [];
   if ~isempty(opt.trapfile_name)
      if ~exist(opt.trapfile_name, 'file')
         warning(['Subtrap datafile not found.  Proceeding without subtraps.']);
      else
         dh = computeToptraps(load(opt.trapfile_name), Gt, true);
      end
   end

   %% Defining fluid
   
   % computing fixed temperature field
   T_ref = opt.ref_temp + opt.temp_grad * (Gt.cells.z - opt.ref_depth) / 1000;   
   
   fluid = makeVEFluid(Gt, rock2D, 'sharp interface'                   , ...
                       'fixedT'       , T_ref                          , ...
                       'wat_rho_pvt'  , [4.3e-5/barsa  , opt.refPress] , ...
                       'residual'     , [opt.sw        , opt.sr]       , ...
                       'wat_rho_ref'  , opt.rhoW                       , ...
                       'co2_rho_ref'  , opt.refRhoCO2                  , ... 
                       'wat_rho_pvt'  , [opt.c_water   , opt.refPress] , ...
                       'co2_rho_pvt'  , [opt.p_range   , opt.t_range]  , ...
                       'co2_mu_pvt'   , [opt.p_range   , opt.t_range]  , ...
                       'wat_mu_ref'   , opt.muBrine                    , ...
                       'pvMult_fac'   , opt.pvMult                     , ...
                       'dissolution'  , false                          , ...
                       'pvMult_p_ref' , opt.refPress                   , ...
                       'surf_topo'    , 'smooth'                    , ...
                       'top_trap'     , dh);
   
   %% Construct schedule
   use_default_schedule = isempty(opt.schedule);
   if use_default_schedule
      [wc, qt] = pick_wellsites(Gt, rock2D, co2, ta, opt.num_wells, opt.rhoW, ...
                                opt.sw, opt.ref_temp, opt.ref_depth, opt.temp_grad, ...
                                opt.well_buffer_dist); 
      qt = qt / 100; %@@
      opt.schedule = setSchedule(Gt, rock2D, wc, qt/opt.refRhoCO2, opt.isteps, ...
                                     opt.itime, opt.msteps , opt.mtime, true, ...
                                     'minval', sqrt(eps));
    
   else
%       tmp = load(opt.schedule);
%       name = fields(tmp);
%       opt.schedule = tmp.(name{1});
        fprintf('\n Using schedule varargin. \n')
        opt.num_wells   = numel(opt.schedule.control(1).W);
   end
   
   % Ensure no repeated well cells
   assert( ~any(diff([opt.schedule.control(1).W.cells])==0) , ...
       'Repeated well cells found, likely due to opt.num_wells > number of catchment areas.')
   
   % Add constant pressure boundary conditions
   % @@ add option to set no-flow boundaries
   bfaces = identifyBoundaryFaces(Gt);
   for i = 1:numel(opt.schedule.control)
      opt.schedule.control(i).bc = addBC([], bfaces, 'pressure', ...
                                         Gt.faces.z(bfaces) * fluid.rhoWS * ...
                                         norm(gravity) + opt.surface_pressure, ...
                                         'sat', [1 0]);
   end
   
   %% Defining initial state
   initState.pressure = computeHydrostaticPressure(Gt, fluid.rhoWS, opt.surface_pressure);
   initState.s = repmat([1 0], Gt.cells.num, 1);
   initState.sGmax = initState.s(:,2);
   
   
   other.fluid = fluid;
   other.rock = rock2D;
   other.residual = [opt.sw, opt.sr];
   other.traps = ta;

   %% Set up model and run optimization
   model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);
   
   min_rates       = sqrt(eps) * ones(opt.num_wells, 1);
   max_rates       = ...
    opt.lim_fac * sum([opt.schedule.control(1).W.val]) * ones(opt.num_wells, 1);

   % Ensure inj well rates are within min and max rates
   vals = [opt.schedule.control(1).W.val]';
   if any(vals<=min_rates)
        warning('An initial well rate is below the minimum rate. Updating initial well rates...')
        vals(vals <= min_rates) = min_rates(1);
   elseif any(vals>=max_rates)
        warning('An initial well rate is above the maximum rate. Updating initial well rates...')
        vals(vals >= max_rates) = max_rates(1);
   end
   
   % Use assert to double-check initial rates within min/max
   assert( all(vals >= min_rates), 'An initial well rate(s) is less than the minimum set rate.')
   assert( all(vals <= max_rates), 'An initial well rate(s) is more than the maximum set rate.')
   vals = num2cell(vals);
   [opt.schedule.control(1).W.val] = vals{:};
   
   [optim, init, history] = ...
       optimizeRates(initState, model, opt.schedule, min_rates, max_rates, ...
                     'last_control_is_migration', true, ...
                     'leak_penalty', opt.leakPenalty, ...
                     'dryrun', opt.dryrun);   %@@
                 
   % pass out initState for later reference.
   other.initState = initState;
   other.min_rates = min_rates;
   other.max_rates = max_rates;
   other.leakPenalty = opt.leakPenalty;
end

% ----------------------------------------------------------------------------
function bfaces = identifyBoundaryFaces(Gt)
   
    bfaces = find(any(Gt.faces.neighbors==0,2)); 

end

% ----------------------------------------------------------------------------
function wc = select_wellcell(Gt, candidates, buffer)

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
                                   tgrad, buf_dist)

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
    strap_vol            = Gt.cells.volumes .* trap_heights .* rock2D.poro .* (1-sw);

    assert(all(trap_heights <= Gt.cells.H));
    
    % Computing trap capacities in mass terms and suggesting wellsites
    trapcap = accumarray(ta.traps(tcells), strap_vol(tcells) .* rhoCO2(tcells));
    
    % Suggest wellsites based on maximum cumulative trap capacity reached

    % bdist = sqrDistanceFromBoundary(Gt); % used below to discourage placement
    % bdist = bdist/max(bdist(:));         % of wells close to boundary

    wc = zeros(num_wells, 1); 
    qt = wc;

    for i = 1:num_wells
        cumcap = [0; cumulative_capacity(trapcap, ta.trap_adj)];
        field  = cumcap(ta.trap_regions + 1);

        % @@ experiment
        qt(i) = max(field); 
        candidate_cells = find(field == qt(i));
        wc(i) = select_wellcell(Gt, candidate_cells, buf_dist);
        
        % % Adding an insignificant quantity so that if everything else is
        % % equal, the cell farthest from the boundary is chosen
        % grid_eps = 10 * eps(max(cumcap(:)));
        % field    = field + grid_eps * bdist;
    
        % [qt(i), wc(i)] = max(field);
        
        % Setting available capacity of used traps to zero
        trapcap = set_used_to_zero(trapcap, ta.trap_adj, ta.trap_regions(wc(i)));
    end

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
    opt.coarse_level = 3;
    opt.num_wells = 10;
    %opt.subtrap_file = 'utsira_subtrap_function_3.mat';

    opt.surface_pressure = 1 * atm;
    
    % Default fluid properties
    [rho_default, mu_default, sr, sw] = getValuesSPE134891();
    opt.refRhoCO2 = rho_default(2);
    opt.rhoW      = rho_default(1);
    opt.muBrine   = mu_default(1);
    opt.muCO2     = 0; % zero value will trigger use of variable viscosity 
    opt.pvMult    = 1e-5/barsa;
    opt.refPress  = 100 * barsa;
    opt.c_water   = 4.3e-5/barsa; % water compressibility
    opt.p_range   = [0.1, 400] * mega * Pascal; % pressure span of sampled property table
    opt.t_range   = [4 250] + 274; % temperature span of sampled property table 

    % residual saturations
    opt.sr = sr; % gas
    opt.sw = sw; % brine
    
    % number of timesteps (injection and migration)
    opt.isteps = 10;
    opt.msteps = 31;

    % durations of injection an migration phases
    opt.itime  = 50   * year;
    opt.mtime  = 3000 * year;
    
    % Temperature regime control parameters
    opt.ref_temp = 7 + 273.15;
    opt.ref_depth = 80; 
    opt.temp_grad = 35.6;

    % dissolution parameters (default zero - suggested nonzero values in comments)
    opt.dis_rate = 0;  % %0.44 * kilogram / rho / poro / (meter^2) / year = 8.6924e-11;
    opt.dis_max = 0;   % 53/760 = 0.07

    % factor for setting upper limit of injection rates
    opt.lim_fac = 10;

    opt.well_buffer_dist = 5 * kilo * meter; % wells should if possible be within
                                             % this distance of spill region boundary
    
    % directory for saving reslts
    opt.report_basedir = './simulateUtsira_results/';
    opt.trapfile_name = ''; % 'utsira_subtrap_function_3'
    
    % extra
    opt.dryrun = false;
    opt.leakPenalty = 1;
end
