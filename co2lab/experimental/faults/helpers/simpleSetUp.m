function [ initState, fluid, rhow, rhoc ] = simpleSetUp( Gt, rock2D )
% Set up initial conditions and model for CO2 injection/migration scenario

    % Fluid data at p = 300 bar
    muw = 0.30860;  rhow = 975.86; sw    = 0.1;
    muc = 0.056641; rhoc = 686.54; srco2 = 0.2;
    kwm = [0.2142 0.85];

    % initial state
    initState.pressure  = Gt.cells.z * norm(gravity) * rhow*kilogram/meter^3;  % hydrostatic pressure, in Pa=N/m^2
    initState.s         = repmat([1 0], Gt.cells.num, 1);                      % sat of water is 1, sat of CO2 is 0
    initState.sGmax     = initState.s(:,2);                                         % max sat of CO2 is initially 0
    initState.rs        = 0 * initState.sGmax;     % initially 0 

    
    ref_p               = mean(initState.pressure);
    info                = getSeaInfo('NorthSea',760);
    caprock_temperature = 273.15 + info.seafloor_temp + ...
            (Gt.cells.z - info.seafloor_depth) / 1e3 * info.temp_gradient; % Kelvin

    fluid = makeVEFluid(Gt, rock2D, 'sharp interface', ...
              'fixedT'      , caprock_temperature, ...
              'wat_mu_ref'  , muw * centi*poise, ...
              'co2_mu_ref'  , muc * centi*poise, ...
              'wat_rho_ref' , rhow * kilogram/meter^3, ...
              'co2_rho_ref' , rhoc * kilogram/meter^3, ...
              'wat_rho_pvt' , [], ...
              'pvMult_p_ref', ref_p, ...
              'pvMult_fac'  , 0, ...
              'residual'    , [sw, srco2] , ...
              'dissolution' , false, 'dis_max', []);
          
    % alternative ?
    %initstate = initState(Gt, W, p0, s0);

end

