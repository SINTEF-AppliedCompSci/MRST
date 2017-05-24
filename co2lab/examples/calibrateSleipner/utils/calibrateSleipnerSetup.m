function [ smodel, schedule, initState, states, wellSols ] = calibrateSleipnerSetup( varargin )
% Setup the IEAGHG Sleipner model, including grid, entry rates, etc.

    [opt.wellXcoord, opt.wellYcoord] = deal(4.38516e5, 6.47121e6);  % Physical coordinate of injection well (Singh et al. 2010)
    opt.ratecase            = 'ORIGINALmodel';  % well injection rate data ('IEAGHGmodel' or 'ORIGINALmodel')
    opt.num_years           = 10;               % use 10 for 1999-2008, and 12 for 1999-2010
    opt.modelname           = 'IEAGHGmodel';    % grid & rock
    opt.refineLevel         = 1;                % Refine using 2, 3, etc. Coarsen using -2, -3, etc.
    opt.modifyParametersCase = 'none';           % ('none' or 'modify')
    
    opt = merge_options(opt, varargin{:});
    
    moduleCheck('co2lab','ad-core', 'deckformat', ...
        'coarsegrid','upscaling','incomp','optimization', 'mrst-gui');
    clear schedule
    gravity on
    
    
    %% Well info, rate data:
    wellXcoord      = opt.wellXcoord;
    wellYcoord      = opt.wellYcoord;
    
    switch opt.ratecase   
        case 'IEAGHGmodel'
            % See Singh et al 2010 for more info about how they determined
            % these rates. Note: the injection rates were reported as
            % surface rates. Both volume and mass were given, thus surface
            % density can be calculated (=1.87 kg/m3). The CO2 density at
            % reservoir conditions was reported as 760 kg/m3.
            inj_year   = [1999; 2000; 2001; 2002; 2003; 2004; 2005; 2006; 2007; 2008; 2009];
            inj_rates  = [2.91e4; 5.92e4; 6.35e4; 8.0e4; 1.09e5; 1.5e5; 2.03e5; 2.69e5; 3.47e5; 4.37e5; 5.4e5] .* meter^3/day;
            % inj_rates is in meter^3/s
            % Convert to rate at reservoir conditions
            inj_rates  = inj_rates.*(1.87/760);       
        case 'ORIGINALmodel' 
            % Note: the CO2 density at reservoir conditions was reported as
            % 695.89 kg/m3, and the surface density 1.87 kg/m3.
            [ inj_year, inj_rates ] = getSleipnerOriginalInjectionRates();
            % inj_rates is in meter^3/s   
        otherwise
            error('The injection rate option was either invalid or not selected.')     
    end
    
    
    %% Specify and compute time step size for injection period.
    % ***Note: inj_time and inj_steps are applied to each inj_rate given***
    num_years   = opt.num_years;
    inj_time    = 1 * year; % DEFAULT. CAN ONLY ADJUST NUMBER OF STEPS.
    inj_steps   = 1;
    dTi         = inj_time / inj_steps; % timestep size in seconds
    
    
    %% Specify fluid properties:
    % Here, we use the benchmark values except for the residual
    % saturations. Seafloor depth info is not given in benchmark, thus a
    % depth of 100m is assumed. Compressibilities were calculated using
    % Span and Wagner's EOS (via coolprops) for reservoir conditions
    % corresponding to initial pressure and temperature (i.e., hydrostatic
    % and computed via thermal gradient), however we set CO2
    % compressibility to be 10 times that of water.
    [rho, mu, sr, sw]   = getValuesSPE134891();
    sr=0; sw=0;
    water_density       = rho(1) * kilogram / meter ^3;
    rhoCref             = rho(2) * kilogram / meter ^3;
    
    seafloor_temp       = 7; % Celsius
    seafloor_depth      = 100; % meters
    temp_gradient       = 35.6; % Celsius / km
    water_compr_val     = 4.37e-5/barsa; % will convert to compr/Pa
    co2_compr_val       = water_compr_val*10; % 1.24e-3/barsa;
    pvMult              = 1e-5/barsa;
    isDissOn            = false;
    dis_max             = (53 * kilogram / meter^3) / rhoCref; % from CO2store
    

    %% Construct grid
    % We only use a portion of the full grid, to speed up our simulations
    [ G_i, ~, rock_i, ~, ~ ] = makeSleipnerVEmodel('usemex',true, 'assign_coords',true);
    
    dx  = (wellXcoord - G_i.cells.centroids(:,1));
    dy  = (wellYcoord - G_i.cells.centroids(:,2));
    ind = dx<1e3 & dx>-1.5e3 & dy<1.3e3 & dy >-4e3;
    
    [G, cellmap] = removeCells(G_i, find(~ind));
    rock.perm    = rock_i.perm(cellmap,:);
    rock.poro    = rock_i.poro(cellmap);
    Gt           = topSurfaceGrid(G);
    rock2D       = averageRock(rock, Gt);
    clear G_i rock_i G rock;

    % Get boundary faces of formation (or grid region)
    %bf = boundaryFaces(Gt);   
    
    
    %% Initial parameter modification:
    switch opt.modifyParametersCase
        case 'modify'
            % Set parameter modifier factors:
            por_mod     = 0.6;
            perm_mod    = 3;
            rhoCO2_mod  = 2/3;
            
            disp('Original rock parameters are being modified ...')       
            %rock.poro   = rock.poro .* por_mod;
            rock2D.poro = rock2D.poro .* por_mod;
            %rock.perm   = rock.perm .* perm_mod;
            rock2D.perm = rock2D.perm .* perm_mod;
            
            disp('Original CO2 density value is being modified ...')
            rhoCref = rhoCref * rhoCO2_mod;
            %rhoCref = 350;
            
        case 'none'
            fprintf('\nRock and density parameters are not modified.\n')    
    end
    
    
    %% Set up initial state:
    initState.pressure  = Gt.cells.z * norm(gravity) * water_density;   % hydrostatic pressure, in Pa=N/m^2
    initState.s         = repmat([1 0], Gt.cells.num, 1);               % sat of water is 1, sat of CO2 is 0
    initState.sGmax     = initState.s(:,2);                             % max sat of CO2 is initially 0
    initState.rs        = 0 * initState.sGmax;                          % initially 0   
    
    
    %% Index of the closest cell to the physical well location: 
    dv              = bsxfun(@minus, Gt.cells.centroids(:,1:2), [wellXcoord, wellYcoord]);
    [v,i]           = min(sum(dv.^2, 2));   
    wellCellIndex   = i;

    
    %% Create schedule:
    clear schedule;
   
    % Well rates:
    % i.e., [inj period 1; inj period 2; etc...; migration period]
    for i = 1:num_years;
        schedule.control(i).W = addWell([], Gt.parent, rock2D, wellCellIndex, ...
            'name', sprintf('W%i', i), 'Type', 'rate', 'Val', inj_rates(i), 'comp_i', [0 1]);
    end
    
    % Boundaries:
    bdryFaces   = find( Gt.faces.neighbors(:,1).*Gt.faces.neighbors(:,2) == 0 );
    bdryType    = 'pressure';
    bdryVal     = Gt.faces.z(bdryFaces) * water_density * norm(gravity);
    bc          = addBC( [], bdryFaces, bdryType, bdryVal, 'sat', [1 0] );   
    for i = 1:numel(schedule.control)
        schedule.control(i).bc = bc;
    end

    % Time step:
    istepvec                = repmat( ones(inj_steps, 1) * dTi , [num_years 1] );
    %mstepvec                = ones(mig_steps, 1) * dTm;    
    schedule.step.val       = [istepvec];
    schedule.step.control   = [];
    for i = 1:numel(schedule.control)
        schedule.step.control = [schedule.step.control; ones(inj_steps, 1) * i];        
    end   
    
    
    %% Create fluid:
    caprock_temperature = 273.15 + seafloor_temp + ...
        (Gt.cells.z - seafloor_depth) / 1e3 * temp_gradient;    % Kelvin
    ref_p               = mean(initState.pressure);             % ref. for linear compressibilities   
    
    fluid = makeVEFluidSens(Gt, rock2D, 'sharp interface', ...
        'fixedT'      , caprock_temperature, ...
        'wat_mu_ref'  , mu(1), ...
        'co2_mu_ref'  , mu(2), ...
        'wat_rho_ref' , water_density, ...
        'co2_rho_ref' , rhoCref, ...
        'wat_rho_pvt' , [water_compr_val, ref_p], ...
        'co2_rho_pvt' , [co2_compr_val  , ref_p], ...
        'pvMult_p_ref', ref_p, ...
        'pvMult_fac'  , pvMult, ...
        'residual'    , [sw, sr] , ...
        'dissolution' , isDissOn, 'dis_max', dis_max);
    
    
    %% Run initial simulation
    % All parameter multipliers are set as 1 for the initial simulation. 
    clear CO2VEBlackOilTypeModelSens
    smodel = CO2VEBlackOilTypeModelSens(Gt, rock2D, fluid);
    initState.dz        = zeros(Gt.cells.num,1);
    initState.rhofac    = 1;
    initState.permfac   = 1;
    initState.porofac   = 1;
    for i=1:numel(schedule.control)
        schedule.control(i).dz      = zeros(Gt.cells.num,1);
        schedule.control(i).rhofac  = 1;
        schedule.control(i).permfac = 1;
        schedule.control(i).porofac = 1;
    end
    %
    [wellSols, states, sim_report] = simulateScheduleAD(initState, smodel, schedule);
    states = addHeightData(states, Gt, fluid);
    
end