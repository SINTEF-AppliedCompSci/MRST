%% Obtain optimized injection rates by maximizing an objective function.

% NB: export_fig used below. Ensure you have downloaded it from
% <http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig>
assert(exist('export_fig')~=0, 'Ensure export_fig exists and is on path.')
moduleCheck('mrst-gui')
moduleCheck('ad-core')

% -All inputs are explicitly defined here and passed into
% optimizeFormation_extras.m (i.e., probably no default values are used
% inside optimizeFormation_extras.m)

% -Wells can be either bhp-controlled or rate-controlled.
% -Well placement can be an array covering the whole or part of the
% formation, or placed in best leaf nodes using internal function in
% optimizeFormations_extras.m

figDirName = 'WellPlacementFigs_array';
%figDirName = 'OptimizedRates_figs_500yrsMigration';
%figDirName = 'OptimizedRates_figs_2000yrsMigration';
%figDirName = 'OptimizedRates_figs_1000yrsMigration';
%figDirName = 'OptimizedRates_figs_2000yrsMigration';
%figDirName = 'OptimizedRates_figs_3000yrsMigration';
mkdir(figDirName)

varDirName = 'opt_results_array_in_trap_regions';
% will be put in subdir named according to mtime

names = [getBarentsSeaNames() getNorwegianSeaNames() getNorthSeaNames()];
%names = {'Arefm'; 'Bjarmelandfm'; 'Brentgrp'; 'Brynefm'; 'Garnfm'; 'Ilefm'; 'Stofm'; 'Tiljefm'; 'Tubaenfm'};
% % Remove ones already run:
% names = names(~strcmpi(names,'Arefm'));
names = names(~strcmpi(names,'Bjarmelandfm'));
% names = names(~strcmpi(names,'Brentgrp'));
% names = names(~strcmpi(names,'Brynefm'));
% names = names(~strcmpi(names,'Garnfm'));
% names = names(~strcmpi(names,'Ilefm'));
names = names(~strcmpi(names,'Stofm'));
% names = names(~strcmpi(names,'Tiljefm'));
names = names(~strcmpi(names,'Tubaenfm'));
% 
% names = names(~strcmpi(names,'Fensfjordfm'));
% names = names(~strcmpi(names,'Krossfjordfm'));
% names = names(~strcmpi(names,'Sognefjordfm'));


% Remove certain formation names:
names = names(~strcmpi(names,'Nordmelafm'));
names = names(~strcmpi(names,'Rorfm'));
names = names(~strcmpi(names,'Notfm'));
names = names(~strcmpi(names,'Knurrfm'));       % @@ can be constructed??
names = names(~strcmpi(names,'Fruholmenfm'));   % @@
names = names(~strcmpi(names,'Cookfm'));
names = names(~strcmpi(names,'Dunlingp'));
names = names(~strcmpi(names,'Paleocene'));


for i=1:numel(names)
    
fmName      = names{i};
%coarsening  = 5; % @@ prevent break that occurs when resolution too low.
rhoCref     = 760 * kilogram / meter ^3;
btype       = 'pressure'; % 'pressure' or 'flux'
use_bhp_wells = false;
use_default_schedule = true;

%%% Get formation
clear Gt rock2D
[Gt, rock2D, coarsening] = best_coarsening_factors( fmName );
%[Gt, rock2D]    = getFormationTopGrid(fmName, coarsening);
if any(isnan(rock2D.perm))
    rock2D.perm = 500*milli*darcy * ones(Gt.cells.num,1);
end
if any(isnan(rock2D.poro))
    rock2D.poro = 0.25 * ones(Gt.cells.num,1); 
end
clear seainfo
seainfo         = getSeaInfo(fmName, rhoCref);

%%% get overburden pressure, and a pressure limit
% [P_over, P_limit] = computeOverburdenPressure(Gt, rock2D, seainfo.seafloor_depth, ...
%     seainfo.water_density);

%%% Get well sites and inj vols
clear trapCapacities
trapCapacities  = getTrappingInfo(Gt, rock2D, seainfo, 'plotsOn',false, 'fmName',fmName);

                    
% pass in wellinfo.initOverP and BHP is computed internally
%wellinfo.initOverP = 1 * mega * Pascal;


%%% Get fluid/formation properties
% (take the ones used to obtain trapCapacities, ...
% which were also used to compute well rates)
clear seafloor_depth seafloor_temp temp_gradient_water_density co2_density water_mu
clear res_sat_co2 res_sat_wat caprock_pressure caprock_temperature
seafloor_depth  = trapCapacities.info.seafloor_depth;
seafloor_temp   = trapCapacities.info.seafloor_temp;
temp_gradient   = trapCapacities.info.temp_gradient;
water_density   = trapCapacities.info.water_density;
co2_density     = []; % @@ not used as input?
water_mu        = trapCapacities.info.water_mu;
res_sat_co2     = trapCapacities.info.res_sat_co2;
res_sat_wat     = trapCapacities.info.res_sat_wat;
caprock_pressure = trapCapacities.caprock_pressure;
caprock_temperature = trapCapacities.caprock_temperature;


%%% Set-up schedule (controls, steps, wells)
clear itime isteps mtime msteps
itime           = 50 * year;
%injRates        = wellinfo.vols_inj./itime; % m3/s
isteps          = 50;
mtime           = 1000 * year;
msteps          = 100;
single_control  = true; % not used in setSchedule_extras() yet.


if ~use_default_schedule
    % schedule is constructed here:
    
    % first get well info (can be an array, etc...)
    clear wellinfo
%     assert(~isempty(find(trapCapacities.ta.traps>0,1)), ...
%         'Skipping formation because no structural traps found. Consider refining top surface.')
%     %@@ could update getWellInfo() st well rates are set as default values
%     %if no traps detected.
    wellinfo        = getWellInfo(Gt, trapCapacities, ...
                        'limits','none', ...
                        'prod',false, ...
                        'setInjRates',true, 'default_rate', sqrt(eps), ...
                        'buffer', 5000, ...
                        'DX', 4*5000, ...
                        'DY', 4*5000, ... %'num_wells_per_frm_dim', 4, ...
                        'plotsOn', true);
                    
                    
    % then assign well values
    if use_bhp_wells

        %%% Pressure-controlled wells:
        % bhp of wells are set to be the initial pressure plus an additional
        % 'initOverP'. This equates to an initial rate that can be calculated
        % by convertBHPtoRate().
        clear schedule
        surface_pressure    = 1 * atm;
        initState.pressure  = seainfo.water_density * norm(gravity()) ...
            * Gt.cells.z + surface_pressure; %@@ contribution from surf_press is small
        schedule = setSchedule_extras( Gt, rock2D, wellinfo.cinx_inj, 'bhp', ...
                                        isteps, itime, msteps, mtime, ...
                                        'initState', initState, ...
                                        'initOverP', 10 * mega * Pascal, ...
                                        'minval',    sqrt(eps));

        % ensure W.sign is set correctly, as it is left as sign=0 when well is
        % bhp-controlled (confirm this).
        for i=1:numel(schedule.control(1).W)
            schedule.control(1).W(i).sign = 1;
        end

    else       




        % NB: vols_inj is already in m3, thus no need to divide by co2 density
        clear schedule
        if ~isempty(wellinfo.vols_inj)
            schedule = setSchedule_extras(Gt, rock2D, wellinfo.cinx_inj, 'rate', ...
                                        isteps, itime, msteps, mtime, ...
                                        'wqtots', wellinfo.vols_inj, ... %vols_inj, ...
                                        'minval', sqrt(eps));
        else
            fprintf('Injection vol is zero. Skipping formation case...\n')
            continue % will go back to start of main 'for' loop, and start the next iteration
        end
    end
    
else
    schedule = [];  % will trigger the construction of schedule within 
                    % optimizeFormation_extras routine
end



try

%%% Pass everything in explicitly.
clear Gt optim init history other
%[Gt, optim, init, history, other] = ...
optimizeFormation_extras(...
    'modelname'                  , fmName                       , ...
    'schedule'                   , schedule                     , ...
    'coarse_level'               , coarsening                   , ...
    'max_num_wells'              , 40                           , ... % updated instead of num_wells
    'trapfile_name'              , []                           , ... %'utsira_subtrap_function_3.mat', ... %'subtrap_file'               , 'utsira_subtrap_function_3.mat', ...
    'surf_topo'                  , 'inf_rough'                  , ...
    'maximise_boundary_distance' , false                        , ... % used in pick_wellsites()
    'well_buffer_dist'           , 1 * kilo * meter             , ... % dist from edge of internal catchment
    'well_buffer_dist_domain'    , 5 * kilo * meter            , ... % dist from edge of domain    
    'well_buffer_dist_catchment' , 3 * kilo * meter             , ... % dist from edge of external catchment
    'pick_highest_pt'            , false                         , ... % otherwise farthest downslope
    'use_array'                  , true                         , ...
    'DX'                         , 1 * kilo*meter               , ...
    'DY'                         , 1 * kilo*meter               , ...
    'surface_pressure'           , 1 * atm                      , ...   
    'refRhoCO2'                  , rhoCref                      , ...
    'rhoW'                       , water_density                , ...
    'muBrine'                    , water_mu                     , ...
    'muCO2'                      , 0                            , ... % zero value will trigger use of variable viscosity 
    'pvMult'                     , 1e-5/barsa                   , ...
    'refPress'                   , mean(caprock_pressure)       , ... % @@ ? 
    'c_water'                    , 4.3e-5/barsa                 , ... % water compressibility
    'p_range'                    , [0.1, 400] * mega * Pascal   , ... % pressure span of sampled property table
    't_range'                    , [4 250] + 274                , ...
    'sr'                         , res_sat_co2                  , ... % gas
    'sw'                         , res_sat_wat                  , ... % brine
    'ref_temp'                   , seafloor_temp + 273.15       , ...
    'ref_depth'                  , seafloor_depth               , ... 
    'temp_grad'                  , temp_gradient                , ...
    'dis_rate'                   , 0                            , ...   % %0.44 * kilogram / rho / poro / (meter^2) / year = 8.6924e-11;
    'dis_max'                    , 0                            , ...   % 53/760 = 0.07
    'itime'                      , itime                        , ...
    'mtime'                      , mtime                        , ...
    'isteps'                     , isteps                       , ...
    'msteps'                     , msteps                       , ...
    'lim_fac'                    , 2                            , ... % factor for setting upper limit of injection rates
    'report_basedir'             , './simulateUtsira_results/'  , ... % directory for saving reslts
    'leak_penalty'               , 10                           , ...
    'dryrun'                     , false                        , ...
    'inspectWellPlacement'       , true                         , ...
    'penalize_leakage_at_infinity', false                       , ...
    'penalize_pressure'           , false                       , ...
    'pressure_penalty'            , 10000                       , ... % @@ get appropriate penalty with trial-and-error
    'btype'                       , btype );


    pause(1)
    saveas(figure(100), [figDirName '/' fmName '_wellPlacement'], 'fig')
    export_fig(figure(100),[figDirName '/' fmName '_wellPlacement'], '-png','-transparent')
    
    % to prevent the following execution:
    %clear varDirName
    
    % save variables
    subVarDirName = [varDirName '/' fmName '/' ...
        'InjYrs',num2str(convertTo(itime,year)), '_MigYrs',num2str(convertTo(mtime,year))];
    mkdir(subVarDirName);
    save([subVarDirName '/' 'Gt'], 'Gt'); % '-v7.3'); using v7.3 makes size larger!
    save([subVarDirName '/' 'optim'], 'optim');
    save([subVarDirName '/' 'init'], 'init');
    save([subVarDirName '/' 'history'], 'history');
    save([subVarDirName '/' 'other'], 'other'); % does not contain other.fluid
    
    
     %%% save figures:
     
     % figure 10 should contain optimization iterations/details
     saveas(figure(10),[figDirName '/' fmName '_optDetails'], 'fig')
     close(figure(10))
     
     % figure 100 should be plot of grid with wells, and initial
     % volumes to inject
     %drawnow
     saveas(figure(100), [figDirName '/' fmName '_wellPlacement'], 'fig')
     %export_fig(figure(100),[figDirName '/' fmName '_wellPlacement'], '-png','-transparent')
     close(figure(100))
     
     % make bar plot of initial and optimized rates
     compareWellrates_viaWellSols(init.wellSols, optim.wellSols, init.schedule, other.fluid.rhoGS);
     %drawnow
     saveas(gcf, [figDirName '/' fmName '_initVsOptimRates'], 'fig')
     %export_fig(gcf, [figDirName '/' fmName], '-png','-transparent')
     close(gcf)

     
     % make inventory plot
     reports_optim = makeReports_extras(Gt, {other.initState, optim.states{:}}, other.rock, other.fluid, optim.schedule, ...
                      other.residual, other.traps, other.dh, ...
                      optim.wellSols);
     h = figure; plot(1);
     ax = get(h, 'currentaxes');  
     load all timesteps up to last plotted one (a bit of a hack)
     plotTrappingDistribution(ax, reports_optim, 'legend_location', 'northwest');
     fsize = 24;
     set(get(gca, 'xlabel'), 'fontsize', fsize)
     set(get(gca, 'ylabel'), 'fontsize', fsize)
     set(gca,'fontsize', fsize);
     set(gcf, 'position', [1, 1, 850, 850]);
     export_fig(gcf,[figDirName '/' fmName '_inventory'], '-png','-transparent')

    
    close all
    

    
catch
    % should continue the for loop if code executed under try finished or
    % failed
    
    
end



end



            %%% Rate-controlled wells:
        % injected vols could be reduced to avoid neg compressibility values
        % that result from very high injection rates.
        %wellinfo.vols_inj = wellinfo.vols_inj .* 0.2; % 0.1 worked, not 0.2 or higher

    %     % Set initial rates according to initial bhp's
    %     surface_pressure    = 1 * atm;
    %     inSt.pressure  = seainfo.water_density * norm(gravity()) ...
    %         * Gt.cells.z + surface_pressure; %@@ contribution from surf_press is small
    %     s = setSchedule_extras( Gt, rock2D, wellinfo.cinx_inj, 'bhp', ...
    %                                     isteps, itime, msteps, mtime, ...
    %                                     'initState', inSt, ...
    %                                     'initOverP', 10 * mega * Pascal, ...
    %                                     'minval',    sqrt(eps));
    % 
    %     % we want to convert the bhp of the wells into the corresponding flux.
    %     % To do this, we must compute the mobility of co2 in wells, which is
    %     % the relative permeability in the well (assume = 1) divided by the co2
    %     % viscosity (a function of temperature and bhp). Grab the viscosity
    %     % from makeVEFluid:
    %     fluid = makeVEFluid(Gt, rock2D, 'sharp interface'                   , ...
    %                    'fixedT'       ,  caprock_temperature                         , ...
    %                    'wat_rho_pvt'  , [4.3e-5/barsa  , mean(caprock_pressure)] , ...
    %                    'residual'     , [res_sat_wat   , res_sat_co2]       , ...
    %                    'wat_rho_ref'  , water_density                       , ...
    %                    'co2_rho_ref'  , rhoCref                  , ... 
    %                    'wat_rho_pvt'  , [4.3e-5/barsa   , mean(caprock_pressure)] , ...
    %                    'co2_rho_pvt'  , [[0.1, 400] * mega * Pascal   , [4 250] + 274]  , ...
    %                    'co2_mu_pvt'   , [[0.1, 400] * mega * Pascal   , [4 250] + 274]  , ...
    %                    'wat_mu_ref'   , water_mu                    , ...
    %                    'pvMult_fac'   , 1e-5/barsa                     , ...
    %                    'dissolution'  , false                          , ...
    %                    'pvMult_p_ref' , mean(caprock_pressure)                   , ...
    %                    'surf_topo'    , 'smooth'                  , ...
    %                    'top_trap'     , []);
    %     mu_g        = fluid.muG(caprock_pressure);
    %     relperm_g   = fluid.krG(ones(Gt.cells.num,1), caprock_pressure); % will be = 1
    %     mob_g       = relperm_g(wellinfo.cinx_inj) ./ mu_g(wellinfo.cinx_inj);
    %     
    %     wellfluxes  = mob_g .* [s.control(1).W.WI]' .* ...
    %         ([s.control(1).W.val]' - inSt.pressure(wellinfo.cinx_inj));
    %     wellfluxes  = wellfluxes ./ Gt.cells.H(wellinfo.cinx_inj) ./ rock2D.poro(wellinfo.cinx_inj); % @@??
    %     vols_inj    = wellfluxes .* itime;

% 
%         %%% fluxes (qGr) obtained at first time step in initial solution using
%         %%% BHP wells:
%         rates= [2.7264;    3.6041;    2.1753;    3.3926;    2.1596;    2.5669; ...
%                 1.0710;    3.1972;    1.2037;    1.4002;    0.8188;    1.1484; ...
%                 1.4871;    2.0375;    2.3875;    2.3250;    1.5426;    1.7688; ...
%                 2.4895;    0.6788;    1.0161]; % m3/2
%         vols_inj = rates .* itime; % m3
%         clear rates


    

    

   
