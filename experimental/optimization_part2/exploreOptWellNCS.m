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


% Main Directory name for saving results (each cases results will be put
% into a subdirectory named according to InjYrs and MigYrs)
varDirName = 'opt_results_Array_in_trap_regions_Pressure_plim90_ClosedBdrys';

% Figure Directory name
figDirName = [varDirName '/' 'WellPlacementFigs_Array_in_trap_regions'];
mkdir(figDirName)

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

% Remove ones already run:
% names = names(~strcmpi(names,'Arefm'));
%names = names(~strcmpi(names,'Bjarmelandfm'));
% names = names(~strcmpi(names,'Brentgrp'));
% names = names(~strcmpi(names,'Brynefm'));
% names = names(~strcmpi(names,'Garnfm'));
% names = names(~strcmpi(names,'Ilefm'));
%names = names(~strcmpi(names,'Stofm'));
%names = names(~strcmpi(names,'Tiljefm'));
%names = names(~strcmpi(names,'Tubaenfm'));
% names = names(~strcmpi(names,'Fensfjordfm'));
% names = names(~strcmpi(names,'Krossfjordfm'));
% names = names(~strcmpi(names,'Sognefjordfm'));

% Ensure no repeting names
names = unique(names,'stable');

% Load res containing formation names and their coarsening levels.
% Or get res from testing_coarsening_levels()
%load coarsening_levels_dx3000meter.mat;  n = {res{:,1}}; c_level = {res{:,2}};
load coarsening_levels_70percent_of_full_StrapCap.mat;
n       = {names_and_cellsizes{:,1}};
c_level = {names_and_cellsizes{:,3}};

shared_names = intersect(names, n, 'stable');
assert( numel(shared_names) >= numel(names) )
assert( all(strcmpi(sort(shared_names),sort(names)))==1 )

clear names
names = {'Brynefm'};

for i=1:numel(names)
    
    fprintf('-------------- FORMATION: %s -----------------\n', names{i})
    fmName      = names{i};
    rhoCref     = 760 * kilogram / meter ^3;

    inx             = find(strcmp(fmName,n));
    coarsening      = c_level{inx};
    [Gt, rock2D]    = getFormationTopGrid( fmName, coarsening );
    if any(isnan(rock2D.perm))
        rock2D.perm = 500*milli*darcy * ones(Gt.cells.num,1);
    end
    if any(isnan(rock2D.poro))
        rock2D.poro = 0.25 * ones(Gt.cells.num,1); 
    end
    
    seainfo = getSeaInfo(fmName, rhoCref);
    gravity on;
    caprock_pressure = (Gt.cells.z * seainfo.water_density * norm(gravity)) ...
                .* (1 + seainfo.press_deviation/100);


    try
    

        %%% Pass everything in explicitly.
        clear Gt optim init history other
       [Gt, optim, init, history, other] = optimizeFormation_extras(...
            'dryrun'                         , false                        , ...
            'inspectWellPlacement'           , false                         , ... % if true, will not continue to simulation
            'modelname'                      , fmName                       , ...
                'coarse_level'               , coarsening                   , ...
            'schedule'                       , [], ... %'opt_results_Array_in_trap_regions_Pressure_plim90/Sleipnerfm/InjYrs50_MigYrs2_DissOn_0_rateLim10_pre/new_schedule'                           , ...
                'itime'                      , 50 * year                    , ...
                'isteps'                     , 50                           , ...
                'mtime'                      , 1000 * year                   , ...
                'msteps'                     , 100                           , ... 
            'well_placement_type'            , 'use_array'                  , ... % 'use_array', 'one_per_trap', 'one_per_path'
                'max_num_wells'              , 40                           , ... % used in use_array, one_per_trap
                'maximise_boundary_distance' , false                        , ... % used in one_per_path
                'well_buffer_dist'           , 1 * kilo * meter             , ... % dist from edge of internal catchment
                'well_buffer_dist_domain'    , 5 * kilo * meter             , ... % dist from edge of domain    
                'well_buffer_dist_catchment' , 3 * kilo * meter             , ... % dist from edge of external catchment
                'pick_highest_pt'            , true                       , ... % otherwise farthest downslope, used in one_per_trap and one_per_path
                'DX'                         , 1 * kilo*meter               , ... % used in use_array
                'DY'                         , 1 * kilo*meter               , ... % used in use_array
            'well_control_type'              , 'rate'                       , ...
                'rate_lim_fac'               , 4                            , ...
            'btype'                          , 'flux'                   , ...
            'penalize_type'                  , 'pressure'                    , ... % 'leakage', 'leakage_at_infinity', 'pressure'
                'leak_penalty'               , 10                           , ...
                'pressure_penalty'           , 100                        , ... % @@ get appropriate penalty with trial-and-error
                'p_lim_factor'               , 0.9                       , ...
            'surface_pressure'              , 1 * atm                       , ...   
            'refRhoCO2'                     , seainfo.rhoCref               , ...
            'rhoW'                          , seainfo.water_density         , ...
            'muBrine'                       , seainfo.water_mu              , ...
            'muCO2'                         , 0                             , ... % zero value will trigger use of variable viscosity 
            'pvMult'                        , 1e-5/barsa                    , ...
            'refPress'                      , mean(caprock_pressure)        , ... % @@ ? 
            'c_water'                       , 4.3e-5/barsa                  , ... % water compressibility
            'p_range'                       , [0.1, 400] * mega * Pascal    , ... % pressure span of sampled property table
            't_range'                       , [4 250] + 274                 , ...
            'sr'                            , seainfo.res_sat_co2                   , ... % gas
            'sw'                            , seainfo.res_sat_wat                   , ... % brine
            'ref_temp'                      , seainfo.seafloor_temp + 273.15        , ...
            'ref_depth'                     , seainfo.seafloor_depth                , ... 
            'temp_grad'                     , seainfo.temp_gradient                 , ...
            'dissolution'                   , false                                  , ...
                'dis_rate'                  , 0                             , ... % 0 means instantaneous, 0.44 * kilogram / rho / poro / (meter^2) / year = 8.6924e-11;
                'dis_max'                   , 53/seainfo.rhoCref            , ... % 53/760 = 0.07; % 1 kg water holds 0.07 kg of CO2
            'report_basedir'                , './simulateUtsira_results/'   , ... % directory for saving reslts    
            'trapfile_name'                 , []                            , ... % 'utsira_subtrap_function_3.mat'
            'surf_topo'                     , 'smooth' );

    
        % Save well inspection figure if it was generated
        if isfield(other,'inspectWellPlacement')
            % Save figure:
            pause(1)
            %saveas(figure(100), [figDirName '/' fmName '_wellPlacement'], 'fig')
            export_fig(figure(100),[figDirName '/' fmName '_wellPlacement'], '-png','-transparent')
            continue
        end
        %
        % Save variables
        subVarDirName = [varDirName '/' fmName '/' ...
            'InjYrs',num2str(convertTo(other.opt.itime,year)), ...
            '_MigYrs',num2str(convertTo(other.opt.mtime,year)), ...
            '_DissOn_',num2str(other.dissolution)];
        mkdir(subVarDirName);
        save([subVarDirName '/' 'Gt'], 'Gt'); % '-v7.3'); using v7.3 makes size larger!
        save([subVarDirName '/' 'optim'], 'optim');
        save([subVarDirName '/' 'init'], 'init');
        save([subVarDirName '/' 'history'], 'history');
        save([subVarDirName '/' 'other'], 'other'); % does not contain other.fluid
        %
        % Save optimization iterations/details:
        saveas(figure(10),[subVarDirName '/' fmName '_optDetails'], 'fig')
        close(figure(10))
        %
        % Save other figure
        saveas(figure(11),[subVarDirName '/' fmName '_optDetails_2'], 'fig')
        close(figure(11))
        %
        close all
        
    catch
        % continue the 'for loop' if code under 'try' either finished or
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