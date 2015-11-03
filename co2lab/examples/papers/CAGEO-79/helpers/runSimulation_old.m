function  [initState, Gt, system, schedule, runADIoptions, fluid] = ...
        runSimulation(Gt, rock, schedule, varargin)
% Runs a full VE simulation on the top-surface grid 'Gt', with the provided
% rock and schedule. 'varargin' allows for additional specification of the
% simulation, as follows:
% 'dh'          - one value per cell, describing the sub-scale trapping
%                 capacity for that cell
% 'report_dir'  - where to write the simulation results (in the format
%                 expected by the 'analyseSimulationResult' routine).
% 'cache_dir'   - where to write timestep files (will allow to start at a
%                 timestep other than 1 for a later run).  If not specified,
%                 timestep files will not be written.
% 'do_plot'     - Whether the function should plot key visuals of the
%                 simulation after each timestep (default: true)
% 'extra_hook'  - additional hook to call after each iterations (in addition
%                 to the plot function, and with same call interface).
% 'dissolution' - dissolution model and parameter to use.
%                 * Finite, positive number - rate-based dissolution (meter
%                                             equivalent per year)
%                 * Inf (infinite)          - instant dissolution
%                 * otherwise               - no dissolution
% 'dis_max'     - Maximum allowed dissolution (default: 0.02).  Only
%                 significant if 'dissolution' is not zero.
% 'start_at'    - which timestep to start at (requires previously computed
%                 cache files in 'cache_dir'.
% 'verbose'     - Whether or not to output auxiliary / diagnostic information
%                 during simulation.
% 'dryrun'      - if 'true', all the data structures are set up and returned,
%                 but the simulation itself will not run.  Used if all you
%                 care about is getting the simulation data structures.
% 'save_result' - if true, results will be saved in 'report_dir'

%% TODO
% plot function
% support for multiple wells
% testing
    
    opt = struct('dh'          , []                           , ...
                 'report_dir'  , '.'                          , ...
                 'cache_dir'   , ''                           , ...
                 'do_plot'     , true                         , ...
                 'extra_hook'  , []                           , ...
                 'dissolution' , 0                            , ...
                 'dis_max'     , 0.02                         , ...
                 'start_at'    , 1                            , ...
                 'dryrun'      , false                        , ...
                 'save_result' , true                         , ...
                 'ref_temp'    , 7+273.15                     , ...
                 'ref_depth'   , 80                           , ...
                 'temp_grad'   , 35.6                         , ...
                 'verbose'     , false                        , ...
                 'refRhoCO2'   , 760 * kilogram/meter^3       , ...
                 'refRhoBrine' , 1020 * kilogram/meter^3      , ...
                 'muCO2'       , 6e-2*milli * Pascal * second , ...
                 'muBrine'     , 8e-4 * Pascal * second       , ...
                 'srCO2'       , 0.21                         , ...
                 'srBrine'     , 0.11);
    
    opt = merge_options(opt, varargin{:});
    
    %% Initial check / fundamental parameters
    moduleCheck('ad-fi', 'ad-props', 'mex', 'co2lab');
    gravity on;

    rho = [opt.refRhoBrine, opt.refRhoCO2];
    mu  =  [opt.muBrine, opt.muCO2];
    sr  = opt.srCO2;
    sw  = opt.srBrine;
    
    surface_pressure   = 1 * atm;
    compressible_rock  = true;
    compressible_fluid = true;
    variable_CO2_mu    = (opt.muCO2 == 0); % a zero value signals that
                                           % variable visosity should be employed
                                           
    dissolution_type   = ifelse(opt.dissolution == 0, ...
                                'none', ...
                                ifelse(opt.dissolution == inf,'instant', 'rate'));
    
    dummy = 0; % placeholder variable for unused arguments in function calls

    %% Set up data structures used for simulation
    
    % Establish basic ADI fluid object with linear rel.perms. (to be replaced
    % by real rel.perms below).
    fluid = initSimpleADIFluid('mu', [dummy, mu], 'rho', [dummy, rho]);
    fluid = purgeUnusedFluidFields(fluid);
    fluid = systemSpecificFluidModifications(fluid, ...
                                             Gt, ...
                                             compressible_rock, ...
                                             compressible_fluid, ...
                                             variable_CO2_mu, ...
                                             dissolution_type, ...
                                             opt.dissolution, ...
                                             opt.dis_max, ...
                                             opt.ref_temp, ...
                                             opt.ref_depth, ...
                                             opt.temp_grad);
    
    % Adding _upscaled_ rel.perm. (and cap.press.) functions to fluid object
    fluid = addVERelperm(fluid, Gt, ...
                         'res_oil'   , sw, ...
                         'res_gas'   , sr, ...
                         'top_trap'  , opt.dh, ...
                         'surf_topo' , ifelse(isempty(opt.dh), 'smooth', 'inf_rough'));

    % establish ADI system structure
    components = ifelse(opt.dissolution > 0, {'Oil', 'Gas', 'DisGas'}, {'Oil', 'Gas'});
    system = initADISystemVE(components, Gt, rock, fluid, ...
                             'simComponents' , setupSimCompVe(Gt ,rock), ...
                             'VE'            , true,                     ...
                             'tol'           , 1e-6);
    
    % Define constant pressure boundary conditions
    bfaces = identifyBoundaryFaces(Gt);
    bcVE = addBC([], bfaces, 'pressure', ...           % Not adjusted for water       
                 Gt.faces.z(bfaces)*fluid.rhoOS * ...   % compressibility, but we      
                 norm(gravity) + surface_pressure);    % assume it is accurate enough.

    % % @@ Lower boundary transmissibilities to emulate 100 km of
    % % external aquifer before reaching ocean bottom (in any direction)
    % system.s.T_all(bfaces) = 0 * system.s.T_all(bfaces) / 133;
    
    % Establish initial state object
    p0        = computeHydrostaticPressure(Gt, fluid.rhoOS, surface_pressure); 
    initState = initResSolVE(Gt, p0, [1 0], 'use_ADI', true);
    if opt.dissolution > 0
        initState.sGmax = zeros(Gt.cells.num,1); 
        initState.rs    = zeros(Gt.cells.num, 1);
    end
    initState = purgeUnusedStateFields(initState);
    
    % computing trapping structure (needed to compute quantity of
    % structurally trapped CO2)
    traps = trapAnalysis(Gt, false);
    
    % Definining plot callback function
    basic_plot_fn = ...
        ifelse(opt.do_plot, ...
               @(state, tstep) simple_panel_plot(state, tstep, Gt, ...
                                                 schedule, sr, sw, fluid, ...
                                                 rock, traps, opt.dh), ...
               []);
    plot_fn = setup_hook_function(basic_plot_fn, opt.extra_hook);

    %% Run the ADI solver to compute all states
    runADIoptions = {'bc'                  , bcVE                    , ...
                     'tstep_hook_fn'       , plot_fn                 , ...
                     'force_step'          , false                   , ...
                     'dt_min'              , 60 * day                , ...
                     'targetIts'           , 6                       , ...
                     'outputName'          , 'cache'                 , ...
                     'outputDir'           , opt.cache_dir           , ...
                     'startAt'             , opt.start_at            , ...
                     'Verbose'             , opt.verbose             , ...
                     'report_all'          , false};

    if opt.dryrun  return; end; % no simulation requested.  We are done.
        
    [wellSols, states] = runMrstADI(initState, Gt, system, schedule, runADIoptions{:});
    
    %% Reporting
    if opt.save_result
        saveResults(opt.report_dir, Gt, states, opt.start_at, schedule, ...
                    traps, rock, fluid, sr, sw, opt.dh, 'initState', initState);
    end
end

% ============================================================================

function fn = setup_hook_function(basic_plot_fn, extra_fn)
    if isempty(basic_plot_fn)
        basic_plot_fn = @(state,tstep) true;
    end
    if isempty(extra_fn)
        extra_fn = @(state,tstep) true;
    end
    fn = @(state, tstep) basic_plot_fn(state,tstep) & extra_fn(state,tstep);
end


% ----------------------------------------------------------------------------
function total = compute_total_vol_injected(tstep, schedule)

    total = 0;

    for i = 1:tstep
        W     = [];
        dt    = schedule.step.val(i);
        wconf = schedule.step.control(i);

        if wconf ~= 0
            W = schedule.W{wconf};
            total = total + sum([W.val]) * dt;
        end
    end
end

% ----------------------------------------------------------------------------
function res = simple_panel_plot(state, tstep, Gt, schedule, sr, sw, fluid, rock, traps, dh)

    persistent fig_h;
    if isempty(fig_h)
        fig_h = figure(93);
    end
    figure(fig_h);
    %% Preparing variables

    % Compute current time
    cur_time = cumsum(schedule.step.val(1:tstep));
    
    % Compute total quantity injected so far
    tot_injected = compute_total_vol_injected(tstep, schedule) * fluid.rhoGS; % ref. rho OK here

    % adding height fields
    [state.h, state.h_max] = computePlumeHeight(Gt, state, sw, sr);
    
    % Compute mass distribution
    %masses = massTrappingDistributionVEADI(Gt, state, rock, fluid, sr, sw, traps, dh);
    masses = massTrappingDistributionVEADI(Gt, state, rock, fluid, traps, dh)
    %% plotting figure
    clf;
    
    % Plotting height of plume on 2D top-surface grid
    subplot(1,2,1);
    cvec = state.h;
    cvec(cvec<0.1) = NaN; % @@ Cutoff value reasonable?
    plotFaces(Gt.parent, Gt.cells.map3DFace, cvec);
    axis tight off;
    colorbar;
    
    % Plotting pie chart
    pa=subplot(1,2,2);
    set(pa, 'Position', [.50 0.2 0.5 0.7]);
    assert(numel(masses) == 7); % Anything else not implemented and tested yet
    
    vplot = max([masses, tot_injected - sum(masses)], eps);
    ph = pie(vplot);
    
    lnames = {'Disolved', ...
              'Struct. residual', ...
              'Residual', ...
              'Free residual', ...
              'Struct. movable', ...
              'Subscale trapped', ...
              'Free movable', ...
              'Leaked'};
    for i=1:numel(lnames)
        legn=sprintf('%s\t %4.1f Mt\n',lnames{i},vplot(i)/1e9)
        lnames{i} = legn;
    end
          
    h = legend(lnames, 'Location','SouthOutside', 'orientation','horizontal');
    set(h, 'Location','SouthEastOutside','orientation','vertical');
    pl = get(h, 'Position'); 
    set(h, 'Position', [.80 pl(2)-.2 pl(3:4)]);
    title(['Total injected mass: ', ...
           num2str(round(convertTo(tot_injected, 1e3*mega))),' M tonn']);
    
    %% General figure properties

    % Set background color to white
    set(gcf, 'color', 'white');
    
    % Set size and position to something useful
    if (tstep == 1)
        scrsz = get(0, 'ScreenSize');
        fsize = min(scrsz(3:4), [1024 700]);
        set(gcf, 'Position',  [scrsz(3)-fsize(1), scrsz(4)-fsize(2)-75, fsize]);
    end
    
    res = gcf; % return handle to graphic window
    drawnow;
end

% ----------------------------------------------------------------------------
function fluid = systemSpecificFluidModifications(fluid, ... 
                                                  Gt, ...
                                                  compressible_rock, ...
                                                  compressible_fluids, ...
                                                  variable_CO2_mu, ...
                                                  dissolution_type, ...
                                                  dis_rate, ...
                                                  dis_max, ...
                                                  ref_temp, ...
                                                  ref_depth,...
                                                  temp_grad)
    % @@ We use the following, hard-coded values for now.
    p_ref = 100 * barsa; % we pick this value as the reference pressure for
                         % the densities (water, CO2) given as initial parameters.
    T_ref = ref_temp + temp_grad * (Gt.cells.z-ref_depth) / 1000;
    % T_ref = (7 + 35.6 * mean(Gt.cells.z)/1000) + 273.15; % thermal gradient from article
    %                                                      % mentioned initially
    rock_mul   = 1e-5   / barsa;
    water_mul  = 4.3e-5 / barsa;
        
    %% Rock compressibility related parameters
    if compressible_rock
        fluid.pvMultR = @(p) 1 + rock_mul * (p - p_ref);
    end
    
    %% Fluid compressibility-related parameters
    if compressible_fluids
        % water
        fluid.bO = @(p, varargin) 1 + water_mul * (p - p_ref);
        fluid.BO = @(p, varargin) 1 ./ fluid.bO(p, varargin{:});
        
       % CO2
        fluid.bG = boCO2(T_ref, fluid.rhoGS, 'sharp_phase_boundary', false);
        fluid.BG = @(p) 1./fluid.bG(p);
    end
    
    %% CO2-visosity-related parameters
    if variable_CO2_mu
        co2props  = CO2props('mufile', 'mu_extensive', 'sharp_phase_boundary', false);
        fluid.muG = @(p, varargin) co2props.mu(p, T_ref);
    end
    
    %% Dissolution-related parameters
    switch dissolution_type
      case 'none'
        % do nothing
      otherwise
        fluid.dis_max = dis_max;
        fluid.muO   = @(p, rs, flag, varargin) fluid.muO(p);
        fluid.rsSat = @(p, rs, flag, varargin) (p*0+1) * dis_max; % necessary?

        if strcmp(dissolution_type, 'rate')
            fluid.dis_rate = dis_rate;
        end
    end
end

% ----------------------------------------------------------------------------
function fluid = purgeUnusedFluidFields(fluid)
% To simplify code understanding, we remove those fields that are
% insignificant for the present model.
    purge = {'krO', 'krW', 'krG', 'pcOG', 'pcOW', 'rhoW', 'rhoWS', ...
             'bW', 'BW', 'muW', 'krOW', 'krOO'};
    for i = 1:numel(purge)
        if isfield(fluid, purge{i})
            fluid = rmfield(fluid, purge{i});
        end
    end
end

% ----------------------------------------------------------------------------
function state = purgeUnusedStateFields(state)
% To simplify code understanding, we remove those fields that are
% insignificant for the present model.
    purge = {'h', 'h_max', 'extSat'};
    for i = 1:numel(purge)
        if isfield(state, purge{i})
            state = rmfield(state, purge{i});
        end
    end
end


% ----------------------------------------------------------------------------
function bfaces = identifyBoundaryFaces(Gt)
   
    bfaces = find(any(Gt.faces.neighbors==0,2)); 

end

% ----------------------------------------------------------------------------
function res = ifelse(cond, a, b)
    if cond
        res = a;
    else
        res = b;
    end
end
