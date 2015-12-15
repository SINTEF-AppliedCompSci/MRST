%% explicitly define all input for running optimizeFormation.m

% set-up wells using BHP.
% Modify obj_fun to include penalty against high overpressures.
% Set max/min BHPs instead of rates as optimization limits.

fmName      = 'Stofm'; %Utsirafm';
coarsening  = 4;
rhoCref     = 760 * kilogram / meter ^3;

%%% get formation
[Gt, rock2D]    = getFormationTopGrid(fmName, coarsening);
seainfo         = getSeaInfo(fmName, rhoCref);

%%% get overburden pressure, and a pressure limit
% [P_over, P_limit] = computeOverburdenPressure(Gt, rock2D, seainfo.seafloor_depth, ...
%     seainfo.water_density);

%%% get well sites and inj vols
trapCapacities  = getTrappingInfo(Gt, rock2D, seainfo, 'plotsOn',false);
wellinfo        = getWellInfo(Gt, trapCapacities, ...
                        'limits','none', ...
                        'prod',false, ...
                        'setInjRates',true);
                    
% pass in wellinfo.initOverP and BHP is computed internally
%wellinfo.initOverP = 1 * mega * Pascal;



itime           = 20 * year;
%injRates        = wellinfo.vols_inj./itime; % m3/s
isteps          = 60;
mtime           = 100 * year;
msteps          = 10;
single_control  = true; % not used in setSchedule_extras() yet.


% put well details into schedule (to be passed in)
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

                            
% schedule = setSchedule(Gt, rock2D, wellinfo.cinx_inj, wellinfo.vols_inj, ...
%                                 isteps, itime, msteps, mtime, single_control, ...
%                                 'minval', sqrt(eps));


%%% get fluid/formation properties
% (take the ones used to obtain trapCapacities, ...
% which were also used to compute well rates)
seafloor_depth  = trapCapacities.info.seafloor_depth;
seafloor_temp   = trapCapacities.info.seafloor_temp;
temp_gradient   = trapCapacities.info.temp_gradient;
water_density   = trapCapacities.info.water_density;
co2_density     = []; % @@ not used as input?
water_mu        = trapCapacities.info.water_mu;
res_sat_co2     = trapCapacities.info.res_sat_co2;
res_sat_wat     = trapCapacities.info.res_sat_wat;
caprock_pressure = trapCapacities.caprock_pressure;


%%% bdrys? (will be added to schedule internally)


%%% others?


%%% Pass everything in explicitly.
[Gt, optim, init, history, other] = optimizeFormation_extras(...
    'modelname'                  , fmName                       , ...
    'schedule'                   , schedule                     , ...
    'coarse_level'               , coarsening                   , ...
    'num_wells'                  , 10 , ... %numel(schedule.control(1).W) , ...% used if wells are picked inside
    'trapfile_name'              , 'utsira_subtrap_function_3.mat', ... %'subtrap_file'               , 'utsira_subtrap_function_3.mat', ...
    'surf_topo'                  , 'inf_rough'                  , ...
    'maximise_boundary_distance' , false                         , ... % used in pick_wellsites()
    'well_buffer_dist'           , 5 * kilo * meter             , ... % used in pick_wellsites()
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
    'dryrun'                     , false ); 

                  
 
%%% Post-processing

reports_init  = makeReports_extras(Gt, {other.initState, init.states{:}}, other.rock, other.fluid, init.schedule, ...
                            other.residual, other.traps, other.dh, ...
                            init.wellSols);

reports_optim = makeReports_extras(Gt, {other.initState, optim.states{:}}, other.rock, other.fluid, optim.schedule, ...
                      other.residual, other.traps, other.dh, ...
                      optim.wellSols);
                  

selectedResultsMultiplot(Gt, reports_optim, [2], ...
                         'plot_plume', false, ...
                         'plot_well_numbering', true, ...
                         'plot_distrib', false);
                     
%compareWellrates(init.schedule, optim.schedule, other.fluid.rhoGS);
compareWellrates_viaWellSols(init.wellSols, optim.wellSols, init.schedule, other.fluid.rhoGS);
                  
% initial rates
selectedResultsMultiplot(Gt, reports_init, [2 4 6], ...
                         'background', 'totalCO2', ...
                         'plot_traps', true);
                     
% optimized rates
selectedResultsMultiplot(Gt, reports_optim, [2 4 6], ...
                         'background', 'totalCO2', ...
                         'plot_traps', true);
  
% overpressure
selectedResultsMultiplot(Gt, reports_optim, [2], ...
                         'background', 'overpressure', ...
                         'init_state', other.initState, ...
                         'plot_traps', true, 'plot_distrib', false, ...
                         'plot_well_numbering', true);

                     
% Also:
figure; plotToolbar(Gt, {other.initState, optim.states{:}})
plotWellSols(optim.wellSols)

figure; plotToolbar(Gt, {other.initState, init.states{:}})
plotWellSols(init.wellSols)
    



    

    

    

   
