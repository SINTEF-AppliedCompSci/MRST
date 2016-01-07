%% Obtain optimized injection rates by maximizing an objective function.

% -All inputs are explicitly defined here and passed into
% optimizeFormation_extras.m (i.e., probably no default values are used
% inside optimizeFormation_extras.m)

% -Wells can be either bhp-controlled or rate-controlled.
% -Well placement can be an array covering the whole or part of the
% formation, or placed in best leaf nodes using internal function in
% optimizeFormations_extras.m

%figDirName = 'WellPlacementFigs';
figDirName = 'OptResultsFigs';

names = [getBarentsSeaNames() getNorwegianSeaNames() getNorthSeaNames()];

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
coarsening  = 5;
rhoCref     = 760 * kilogram / meter ^3;
use_bhp_wells = false;
use_default_schedule = true;

%%% Get formation
[Gt, rock2D]    = getFormationTopGrid(fmName, coarsening);
if any(isnan(rock2D.perm))
    rock2D.perm = 500*milli*darcy * ones(Gt.cells.num,1);
end
if any(isnan(rock2D.poro))
    rock2D.poro = 0.25 * ones(Gt.cells.num,1); 
end
seainfo         = getSeaInfo(fmName, rhoCref);

%%% get overburden pressure, and a pressure limit
% [P_over, P_limit] = computeOverburdenPressure(Gt, rock2D, seainfo.seafloor_depth, ...
%     seainfo.water_density);

%%% Get well sites and inj vols
trapCapacities  = getTrappingInfo(Gt, rock2D, seainfo, 'plotsOn',false, 'fmName',fmName);
wellinfo        = getWellInfo(Gt, trapCapacities, ...
                        'limits','none', ...
                        'prod',false, ...
                        'setInjRates',true, ...
                        'buffer', 5000, ...
                        'DX', 5*5000, ...
                        'DY', 5*5000 );
                    
% pass in wellinfo.initOverP and BHP is computed internally
%wellinfo.initOverP = 1 * mega * Pascal;


%%% Get fluid/formation properties
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
caprock_temperature = trapCapacities.caprock_temperature;



%%% Set-up schedule (controls, steps, wells)
itime           = 50 * year;
%injRates        = wellinfo.vols_inj./itime; % m3/s
isteps          = 50;
mtime           = 100 * year;
msteps          = 10;
single_control  = true; % not used in setSchedule_extras() yet.

if ~use_default_schedule
if use_bhp_wells
    
    %%% Pressure-controlled wells:
    % bhp of wells are set to be the initial pressure plus an additional
    % 'initOverP'. This equates to an initial rate that can be calculated
    % by convertBHPtoRate().
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
    
    
    %%% fluxes (qGr) obtained at first time step in initial solution using
    %%% BHP wells:
    rates= [2.7264;    3.6041;    2.1753;    3.3926;    2.1596;    2.5669; ...
            1.0710;    3.1972;    1.2037;    1.4002;    0.8188;    1.1484; ...
            1.4871;    2.0375;    2.3875;    2.3250;    1.5426;    1.7688; ...
            2.4895;    0.6788;    1.0161]; % m3/2
    vols_inj = rates .* itime; % m3
    clear rates
    
    

    % NB: vols_inj is already in m3, thus no need to divide by co2 density
    schedule = setSchedule_extras(Gt, rock2D, wellinfo.cinx_inj, 'rate', ...
                                    isteps, itime, msteps, mtime, ...
                                    'wqtots', vols_inj, ...
                                    'minval', sqrt(eps));
end
end



try

%%% Pass everything in explicitly.
[Gt, optim, init, history, other] = ...
optimizeFormation_extras(...
    'modelname'                  , fmName                       , ...
    'schedule'                   , [], ... %schedule                     , ...
    'coarse_level'               , coarsening                   , ...
    'num_wells'                  , 50 , ... %numel(schedule.control(1).W) , ...% used if wells are picked inside
    'trapfile_name'              , [], ... %'utsira_subtrap_function_3.mat', ... %'subtrap_file'               , 'utsira_subtrap_function_3.mat', ...
    'surf_topo'                  , 'inf_rough'                  , ...
    'maximise_boundary_distance' , true                         , ... % used in pick_wellsites()
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
    
%     % save fig
%     %pause
%     drawnow % @@
%     %export_fig(gcf,[figDirName '/' fmName], '-png','-transparent')
%     saveas(gcf, [figDirName '/' fmName], 'fig')
    compareWellrates_viaWellSols(init.wellSols, optim.wellSols, init.schedule, other.fluid.rhoGS);
    drawnow
    saveas(gcf, [figDirName '/' fmName '_InitVsOptimRates'], 'fig')
    
catch
    % should continue the for loop if code executed under try finished or
    % failed
    
    
end



end
% 
% 
% %%% Post-processing
% 
% reports_init  = makeReports_extras(Gt, {other.initState, init.states{:}}, other.rock, other.fluid, init.schedule, ...
%                             other.residual, other.traps, other.dh, ...
%                             init.wellSols);
% 
reports_optim = makeReports_extras(Gt, {other.initState, optim.states{:}}, other.rock, other.fluid, optim.schedule, ...
                      other.residual, other.traps, other.dh, ...
                      optim.wellSols);
%                   
% 
selectedResultsMultiplot(Gt, reports_optim, [2], ...
                         'plot_plume', false, ...
                         'plot_well_numbering', true, ...
                         'plot_distrib', false);
                     
% %compareWellrates(init.schedule, optim.schedule, other.fluid.rhoGS);
compareWellrates_viaWellSols(init.wellSols, optim.wellSols, init.schedule, other.fluid.rhoGS);
%                   
% % initial rates
% selectedResultsMultiplot(Gt, reports_init, [2 4 6], ...
%                          'background', 'totalCO2', ...
%                          'plot_traps', true);
%                      
% % optimized rates
selectedResultsMultiplot(Gt, reports_optim, [25 61 71], ...
                         'background', 'totalCO2', ...
                         'plot_traps', true);
%   
% % overpressure
selectedResultsMultiplot(Gt, reports_optim, [61], ...
                         'background', 'overpressure', ...
                         'init_state', other.initState, ...
                         'plot_traps', true, 'plot_distrib', false, ...
                         'plot_well_numbering', true);
% 
%                      
% % Also:
optim.states = computeOverpressure(other.initState, optim.states, optim.schedule);
% figure; plotToolbar(Gt, {other.initState, optim.states{:}})
% plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
% plotWellSols(optim.wellSols); set(gcf,'name','optim: bhp-wells')
% 
% 
init.states = computeOverpressure(other.initState, init.states, init.schedule);
% figure; plotToolbar(Gt, {other.initState, init.states{:}})
% plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
% plotWellSols(init.wellSols); set(gcf,'name','init: bhp-wells')
%     

% sat plots just before migration period
figure
subplot(1,2,1)
plotCellData(Gt, init.states{60}.s(:,2), 'EdgeAlpha',0.1)
plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
colorbar
axis equal tight
title('Initial')

subplot(1,2,2)
plotCellData(Gt, optim.states{60}.s(:,2), 'EdgeAlpha',0.1)
plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
colorbar
axis equal tight
title('Optimized')

% press plots just before migration period (likely when the highest
% overpressure occurred)
figure
subplot(1,2,1)
plotCellData(Gt, init.states{60}.pressDev, 'EdgeAlpha',0.1)
plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
colorbar
axis equal tight
title('Initial')

subplot(1,2,2)
plotCellData(Gt, optim.states{60}.pressDev, 'EdgeAlpha',0.1)
plotGrid(Gt, 'EdgeAlpha',0.1, 'FaceColor','none')
colorbar
axis equal tight
title('Optimized')


[~, P_limit] = computeOverburdenPressure(Gt, rock2D, seainfo.seafloor_depth, ...
    seainfo.water_density);
plotFormationPressureChanges(optim.states, other.initState.pressure, ...
    'P_lim', P_limit)


    

    

    

   
