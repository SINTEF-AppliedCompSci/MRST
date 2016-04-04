%% Obtain optimized injection rates by maximizing an objective function.

% All inputs are explicitly defined here and passed into
% optimizeFormation_extras.m. Model, schedule, box limits, etc. are
% computed internally in optimizeFormation_extras.m, before calling
% optimizeRates_extras.m. Objective functions are defined in
% optimizeRates_extras.m

% Wells can be either bhp-controlled or rate-controlled (bhp-controls
% require more testing).

% Well placement can be an array covering the catchment regions, placed in
% the best leaf nodes, or placed in each trap at the highest point. Initial
% injection volumes and rates are determined based on trapping capacities.
% If system is closed, these rates are adjusted (either decreased or
% increased) according to the analytical closed-system's storage efficiency.

% Penalization types are: 'leakage' (which is leakage only), 'pressure'
% (which is leakage and pressure), 'leakage_at_infinity' (requires more
% testing).

% Pressure penalty factor, cp, is iteratively determined. That is, we begin
% with a low number, obtain an optimized solution, assess how much the
% pressure limit was surpassed, and decide whether or not we need to obtain
% another optimized solution given a larger cp value. If the decision is to
% run the optimization framework again with a larger cp, the initial rates
% are the previous iterations optimized rates.

% Pressure limit (plim) is computed as a percentage of the formation
% overburden pressure. This percentage is set using 'p_lim_factor', and is
% typically set as 0.9.

% An open-system is specified by 'btype'='pressure', and a closed-system is
% 'btype'='flux'.

%names = {'Stofm'}
names = {'Synthetic'}

for i=1:numel(names)
    
    fprintf('-------------- FORMATION: %s -----------------\n', names{i})
    fmName      = names{i};
    rhoCref     = 760 * kilogram / meter ^3;

    if ~strcmpi(names{1},'Synthetic')
        coarsening = 6; % @@
        [Gt, ~]    = getFormationTopGrid( fmName, coarsening );
        seainfo    = getSeaInfo(fmName, rhoCref);
    else
        coarsening = [];
        Gt      = dipped_perturbed_grid('Lx', 10000, 'Ly', 5000, 'H', 50);
        seainfo = getSeaInfo('NorthSea', rhoCref); 
    end
    gravity on;
    caprock_pressure = (Gt.cells.z * seainfo.water_density * norm(gravity)) ...
                .* (1 + seainfo.press_deviation/100);

    max_rate_fac = 4;
    if any(strcmpi(fmName,{'Sleipnerfm','Huginfmwest','Ulafm','Huginfmeast','Pliocenesand'}))
        max_rate_fac = 4;  % @@ max rates will depend on the starting rates
    %elseif any(strcmpi(fmName,{''}))
    %    max_rate_fac = 4;
    end

    

    %%% Pass everything in explicitly.
    clear Gt optim init history other
    cp = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];    % pressure penalty factor
    sch = [];       % computed internally according to well_placement_type
    max_wvals = []; % will be computed internally using rate_lim_fac
    for r = 1:numel(cp)
        % Assuming penalize_type is 'pressure', get optimized rates for a
        % given cp. If penalize_type is 'leakage', cp is not required and
        % only one iteration of this r-loop will be executed.
        [Gt, optim, init, history, other] = optimizeFormation_extras(...
        'dryrun'                         , true                       , ... % if true, will not proceed to optimization
        'inspectWellPlacement'           , false                        , ... % if true, will not proceed to any simulation
        'adjustClosedSystemRates'        , true                         , ... % if true, initial schedule rates are adjusted for closed-systems, based on analytical closed-system storage efficiency
        'lineSearchMaxIt'                , 5                           , ...
        'gradTol'                        , 1e-4                         , ...
        'objChangeTol'                   , 1e-5, ...%1e-4                         , ...
        'modelname'                      , fmName                       , ...
            'coarse_level'               , coarsening                   , ...
        'schedule'                       , sch, ...%'/home/rebeccaa/Documents/MATLAB/mrst-bitbucket/Synthetic/InjYrs30_MigYrs500_noResTrap_cp2_futLeak_IR0pt5/optim', ...%sch                          , ... % can also pass in path of a previously saved schedule
            'itime'                      , 30 * year                    , ...
            'isteps'                     , 30                           , ...
            'mtime'                      , 3000 * year                  , ...
            'msteps'                     , 300                          , ... 
            'extended_mig_time'          , 3000 * year                  , ... % if sch is passed in, mig period can be extended
        'well_placement_type'            , 'one_per_path'               , ... % options: 'use_array', 'one_per_trap', 'one_per_path'
            'max_num_wells'              , 2                           , ... % used in use_array, one_per_trap
            'maximise_boundary_distance' , false                        , ... % used in one_per_path
            'well_buffer_dist'           , 0 * kilo * meter, ...%1 * kilo * meter             , ... % dist from edge of internal catchment
            'well_buffer_dist_domain'    , 2 * kilo * meter, ...%5 * kilo * meter             , ... % dist from edge of domain    
            'well_buffer_dist_catchment' , 0 * kilo * meter, ...%3 * kilo * meter             , ... % dist from edge of external catchment
            'pick_highest_pt'            , true                         , ... % otherwise farthest downslope, used in one_per_trap and one_per_path
            'DX'                         , 1 * kilo*meter               , ... % used in use_array
            'DY'                         , 1 * kilo*meter               , ... % used in use_array
        'well_control_type'              , 'rate'                       , ...
            'rate_lim_fac'               , max_rate_fac                 , ...
            'max_wvals'                  , max_wvals                    , ...
        'btype'                          , 'pressure'                   , ... % options: 'pressure', 'flux'
        'penalize_type'                  , 'leakage'                   , ... % options: 'pressure', 'leakage', 'leakage_at_infinity'
            'leak_penalty'               , 2                           , ...
            'pressure_penalty'           , cp(r)                        , ... % appropriate penalty is iteratively determined for closed-systems
            'p_lim_factor'               , 0.9                          , ...
        'surface_pressure'              , 1 * atm                       , ...   
        'refRhoCO2'                     , seainfo.rhoCref               , ...
        'rhoW'                          , seainfo.water_density         , ...
        'muBrine'                       , seainfo.water_mu              , ...
        'muCO2'                         , 0                             , ... % zero value will trigger use of variable viscosity 
        'pvMult'                        , 1e-5/barsa                    , ...
        'refPress'                      , mean(caprock_pressure)        , ... % 
        'c_water'                       , 4.3e-5/barsa                  , ... % water compressibility
        'p_range'                       , [0.1, 400] * mega * Pascal    , ... % pressure span of sampled property table
        't_range'                       , [4 250] + 274                 , ...
        'sr'                            , 0, ...%seainfo.res_sat_co2                   , ... % gas
        'sw'                            , 0, ...%seainfo.res_sat_wat                   , ... % brine
        'ref_temp'                      , seainfo.seafloor_temp + 273.15        , ...
        'ref_depth'                     , seainfo.seafloor_depth                , ... 
        'temp_grad'                     , seainfo.temp_gradient                 , ...
        'dissolution'                   , false                                  , ...
            'dis_rate'                  , 0                             , ... % 0 means instantaneous, 0.44 * kilogram / rho / poro / (meter^2) / year = 8.6924e-11;
            'dis_max'                   , 53/seainfo.rhoCref            , ... % 53/760 = 0.07; % 1 kg water holds 0.07 kg of CO2
        'report_basedir'                , './simulateUtsira_results/'   , ... % directory for saving reslts    
        'trapfile_name'                 , []                            , ... % 'utsira_subtrap_function_3.mat'
        'surf_topo'                     , 'smooth' );

        % Save well inspection figure if it was generated,
        % then go to next formation
        if isfield(other,'inspectWellPlacement')
            % Save figure:
            pause(1)
            saveas(figure(100), [figDirName '/' fmName '_wellPlacement'], 'fig')
            break % exit cp(r) loop
        end

        % keep the 'first' initial schedule for post-processing, and
        % keep the 'first' max well values for further cp iterations
        if r == 1
           init0 = init;
           max_wvals = other.opt.rate_lim_fac * ...
                max([init0.schedule.control(1).W.val]) * ones(numel(init0.schedule.control(1).W), 1);
        end
        if strcmpi(other.opt.penalize_type,'pressure')
            % Was pressure limit plus a tolerance surpassed? If yes, use
            % next higher cp value. If within plim + tolerance, results
            % obtained with cp value were acceptable. If system is closed
            % and pressure is under plim, rates could be higher, thus
            % increase initial rates slightly and go to next cp. If system
            % is closed and pressure is under plim, okay.
            plim = other.opt.p_lim_factor * other.P_over;
            [perc_of_plim_reach, perc_of_Pover_reach] = ...
                report_maxPercentage_plim_reached( optim.states, plim, other.P_over );
            if perc_of_Pover_reach/100 - other.opt.p_lim_factor > 0.02
                % use optimized rates as next iteration's initial rates
                sch = optim.schedule;
            elseif strcmpi(other.opt.btype,'flux') && perc_of_Pover_reach/100 < other.opt.p_lim_factor
                % plim was not reached in closed system, thus rates could
                % be higher. Increase next iteration's initial rates by
                % some percentage, while within the max rate limit
                sch = optim.schedule;
                for wi=1:numel([sch.control(1).W.val])
                    if (sch.control(1).W(wi).val * 1.25) <= max_wvals(wi)
                        sch.control(1).W(wi).val = sch.control(1).W(wi).val * 1.25; % 25 percent increase
                    else
                        sch.control(1).W(wi).val = max_wvals(wi);
                    end
                end
            else
                % closed system: if plim < p < (plim + tolerance), thus optimal rates found
                % open system: if p < (plim + tolerance), optimal rates found
                break % exit cp(r) loop
            end
        else
            % penalize leakage or leakage at infinity (without
            % pressure) doesn't require iteratively increasing cp.
            break % exit cp(r) loop
        end
    end
    %
    % Save variables if they exist
    if exist('other','var')
        subVarDirName = [fmName '/' ...
            'InjYrs',num2str(convertTo(other.opt.itime,year)), ...
            '_MigYrs',num2str(convertTo(other.opt.mtime,year)) '_noResTrap_cp2_futLeak_IR0pt5_longSim'];
        mkdir(subVarDirName);
        save([subVarDirName '/' 'Gt'], 'Gt');
        save([subVarDirName '/' 'optim'], 'optim');
        save([subVarDirName '/' 'init0'], 'init0');
        save([subVarDirName '/' 'init'], 'init');
        save([subVarDirName '/' 'history'], 'history');
        save([subVarDirName '/' 'other'], 'other'); % to reduce size, does not contain fluid structure
    end
    %
    % Save optimization iterations/details:
    % (Avoid stealing focus if figures already exist)
    if ishandle(10)
        set(0, 'CurrentFigure', 10);
        saveas(10,[subVarDirName '/' fmName '_optDetails'], 'fig')
        close(10)
    end
    %
    if ishandle(11)
        set(0, 'CurrentFigure', 11);
        title('Final')
        saveas(11,[subVarDirName '/' fmName '_optDetails_2'], 'fig')
        close(11)
    end
    %
    if ishandle(50)
        set(0, 'CurrentFigure', 50); % only exists if penalizing pressure
        title('Initial')
        saveas(50,[subVarDirName '/' fmName '_optDetails_3'], 'fig')
        close(50)
    end
    %
    close all


end