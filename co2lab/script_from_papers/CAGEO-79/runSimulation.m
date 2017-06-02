function  [initState, Gt, schedule, fluid] = runSimulation(Gt, rock, schedule, varargin)
% Runs a full VE simulation on the top-surface grid 'Gt', with the provided
% rock and schedule. 'varargin' allows for additional specification of the
% simulation, as follows:
% 'dh'          - one value per cell, describing the sub-scale trapping
%                 capacity for that cell
% 'report_dir'  - where to write the simulation results (in the format
%                 expected by the 'analyseSimulationResult' routine).
% 'do_plot'     - Whether the function should plot key visuals of the
%                 simulation after each timestep (default: true)
% 'dissolution' - dissolution model and parameter to use.
%                 * Finite, positive number - rate-based dissolution (meter
%                                             equivalent per year)
%                 * Inf (infinite)          - instant dissolution
%                 * zero                    - no dissolution
% 'dis_max'     - Maximum allowed dissolution (default: 0.02).  Only
%                 significant if 'dissolution' is not zero.
% 'dryrun'      - if 'true', all the data structures are set up and returned,
%                 but the simulation itself will not run.  Used if all you
%                 care about is getting the simulation data structures.
% 'save_result' - if true, results will be saved in 'report_dir'
    
   opt = struct('dh'          , []                           , ...
                'report_dir'  , '.'                          , ...
                'do_plot'     , true                         , ...
                'extra_hook'  , []                           , ...
                'dissolution' , 0                            , ...
                'dis_max'     , 0.02                         , ...
                'dryrun'      , false                        , ...
                'save_result' , true                         , ...
                'ref_temp'    , 7+273.15                     , ...
                'ref_depth'   , 80                           , ...
                'temp_grad'   , 35.6                         , ...
                'refRhoCO2'   , 760 * kilogram/meter^3       , ...
                'refRhoBrine' , 1020 * kilogram/meter^3      , ...
                'muCO2'       , 6e-2*milli * Pascal * second , ...
                'muBrine'     , 8e-4 * Pascal * second       , ...
                'pvMult'      , 1e-5 / barsa                 , ...
                'ref_press'   , 100 * barsa                  , ...
                'c_water'     , 4.3e-5/barsa                 , ...
                'srCO2'       , 0.21                         , ...
                'srBrine'     , 0.11);
   
   opt = merge_options(opt, varargin{:});
   
   mrstModule add ad-props mex co2lab
   gravity on;
   
   % specified range for sampled property tables
   p_range = [0.1, 400] * mega * Pascal;
   t_range = [4 250] + 274; % in Kelvin
   
   %% make fluid model
   fluid = makeVEFluid(Gt, rock, 'sharp interface', ...
                       'surf_topo'    , ifelse(isempty(opt.dh) , 'smooth', 'inf_rough')    , ...
                       'top_trap'     , opt.dh                                             , ...
                       'pvMult_fac'   , opt.pvMult                                         , ...
                       'pvMult_p_ref' , opt.ref_press                                      , ...
                       'wat_rho_ref'  , opt.refRhoBrine                                    , ...
                       'co2_rho_ref'  , opt.refRhoCO2                                      , ...
                       'wat_rho_pvt'  , [opt.c_water, opt.ref_press]                       , ...
                       'co2_rho_pvt'  , [p_range t_range]                                  , ...
                       'co2_mu_pvt'   , ifelse(opt.muCO2 == 0                              , ...
                                               [p_range t_range]                           , ...
                                               opt.muCO2)                                  , ...
                       'wat_mu_ref'  , opt.muBrine                                         , ...
                       'residual'    , [opt.srBrine, opt.srCO2]                            , ...
                       'dissolution' , opt.dissolution ~= 0                                , ...
                       'dis_rate'    , ifelse(isinf(opt.dissolution) , 0, opt.dissolution) , ...
                       'dis_max'     , opt.dis_max                                         , ...
                       'fixedT'      , opt.ref_temp + ...
                                       opt.temp_grad * (Gt.cells.z-opt.ref_depth) / 1000);
   
   %% defining initial state
   initState = struct('pressure', compute_hydrostatic_pressure(Gt, fluid.rhoWS, 1 * atm), ...
                      's'       , [ones(Gt.cells.num, 1), zeros(Gt.cells.num, 1)], ...
                      'sGmax'   , zeros(Gt.cells.num, 1), ...
                      'rs'      , zeros(Gt.cells.num, 1));
   
   %% check if we are done
   if opt.dryrun return; end; %#ok
   
   %% setup and run complete model

   traps = trapAnalysis(Gt, false);   
   plot_fn = []; 
   if (opt.do_plot)
      plot_fn = @(model, states, reports, solver, schedule, simtime) ...
                deal(model, states, reports, solver, ...
                     simple_panel_plot(states{numel([states{:}])}, ... % last computed state
                                       numel([states{:}])        , ... % time step
                                       Gt                        , ...
                                       schedule                  , ...
                                       opt.srCO2, opt.srBrine    , ...
                                       fluid                     , ...
                                       rock                      , ...
                                       traps                     , ...
                                       opt.dh));
   end
   
   [wellSols, states] = simulateScheduleAD(initState, ...
                                           CO2VEBlackOilTypeModel(Gt, rock, fluid), schedule, ...
                                           'afterStepFn', plot_fn, ...
                                           'NonLinearSolver', NonLinearSolver('useRelaxation', true));%#ok

   %% storing outcome if requested
   if opt.save_result
      saveResults(opt.report_dir, Gt, states, 1, schedule, traps, rock, ...
                  fluid, opt.srCO2, opt.srBrine, opt.dh, 'initState', initState);
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
    [state.h, state.h_max] = compute_plume_height(Gt, state, sw, sr);
    
    % dissolution info, if any
    if isfield(state, 'rs')
       rs = state.rs;
    else
       rs = 0;
    end

    % Compute mass distribution
    masses = massTrappingDistributionVEADI(Gt, state.pressure, state.s(:,2), ...
                                           state.s(:,1), h, h_max, rock, ...
                                           fluid, traps, dh, 'rs', rs);
    
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
        legn=sprintf('%s\t %4.1f Mt\n',lnames{i},vplot(i)/1e9);
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
function total = compute_total_vol_injected(tstep, schedule)

    total = 0;

    for i = 1:tstep
        W     = [];
        dt    = schedule.step.val(i);
        wconf = schedule.step.control(i);

        if wconf ~= 0
            W = schedule.control(wconf).W;
            total = total + sum([W.val]) * dt;
        end
    end
end

 
 
% ----------------------------------------------------------------------------

 function res = ifelse(cond, yes, no)
    if cond
       res = yes;
    else
       res = no;
    end
 end
 
 % ----------------------------------------------------------------------------
 
 function p = compute_hydrostatic_pressure(Gt, rho_water, surface_pressure)

    p = rho_water * norm(gravity()) * Gt.cells.z + surface_pressure;
    
 end
 
 % ----------------------------------------------------------------------------
 
 function [h, h_max] = compute_plume_height(Gt, state, sw, sr)
    
    if isfield(state, 'sGmax')
       smax = state.sGmax; % we operate with dissolution
    else
       smax = state.smax(:,2); % no dissolution.  Current max = historical max
    end
    [h, h_max] = upscaledSat2height(state.s(:,2), smax, Gt, 'resSat', [sw, sr]);
 end
 