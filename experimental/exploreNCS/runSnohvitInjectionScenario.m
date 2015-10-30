function [ wellSols, states, sim_report, opt, var ] = runSnohvitInjectionScenario( Gt, rock2D, varargin )
% Run the Snohvit injection scenario.
%
%   Overview (Hansen et al 2013):
%       - gas production started in August 2007
%       - co2 injection began in April 2008
%       - injected volume in 2008 was lower than expected. Pressure
%       increased quickly during fall 2008
%       - precipitation of salt near wellbore formation is blamed as cause
%       for reduced injectivity and increased pressure
%       - after 3 attempted perforations (to reduce injection pressure)
%       around April 2011, Tubaen storage site was sealed off, and well
%       perforated at a shallower level (i.e., in the Sto formation)
%
%       - between April 2008 and Aug 2012, 1.6 Mt CO2 was injected (not
%       continuously due to interruptions, plant maintence, etc.)
%
%   Estimated injection (Hansen et al 2013):
%       - 30 year lifetime of Snohvit LNG project
%       - ~ 23 Mt CO2 captured for storage (12 e9 m3 at surface conditions)
%       - inj. rate: 2 kt/day
%       - inj. reservoir rate: 2100 to 2750 m3/day
%       - design capacity: 0.7 Mt/yr
%       - inj rates from 2008-2009 are provided in graph in Hansen et al
%       2011
%
    gravity on;
    moduleCheck('ad-core');

    rhoCref = 760 * kilogram / meter ^3; % an (arbitrary) reference density

    %opt.grid_coarsening   = 1;
    %opt.default_formation = 'Stofm'; % the top formation of Hammerfest Basin Aquifer
    % pass in Gt grid and rock2D properties into function. Already
    % coarsened, refined, etc.

    % Pham et al 2011 states Tubaen formation has an average temp of 98 C,
    % at a depth of 2600 meters. The temp gradient should be used to meet
    % this value. Pham et al 2011 also gives a diffusion coefficient range
    % of 4.5e-4 to 4.7e-4 cm2/s, and reports the planned injection rate and
    % period (23 Mt over 30 years = 0.767 Mt per year).
    opt.seafloor_depth    = 330 * meter;     % NPD
    opt.seafloor_temp     =  4;              % in Celsius
    opt.temp_gradient     = 40;              % degrees per kilometer
    % NOTE: see fig 6 in Halland and Riis 2014 for temp gradient data!

    
    % TO CONFIRM for Snohvit
    [rho, mu, sr, sw]   = getValuesSPE134891(); % from Sleipner benchmark
    opt.water_density   = rho(1); % note: chp 6, pg 128 in Atlas says 1.1 g/cm3
    % and formation water is strongly saline.
    opt.co2_density     = rho(2);
    opt.water_mu        = mu(1);
    opt.co2_mu          = mu(2);
    opt.res_sat_wat     = sw;
    opt.res_sat_co2     = sr;
    opt.water_compr_val = 4.3e-5/barsa;
    opt.pvMult          = 1e-5/barsa; % pore volume multiplier
    opt.isDissOn        = false;
    opt.dis_max         = (53 * kilogram / meter^3) / rhoCref; % value from CO2store

    
    % Snohvit injection scenario details:
    
    % Physical coordinate of Snohvit injection (approx, see chp 6 pg 135 in Atlas):
    opt.wellCoords        = [9.225e5, 7.988e6];
    opt.well_radius       = 0.3;
    % 
    opt.inj_rate_MtperYr  = 0.767; % Mt/yr
    opt.inj_time          = 30 * year;
    opt.inj_steps         = 10;
    opt.mig_time          = 1 * year;
    opt.mig_steps         = 1;
    
    % Boundary condition type
    opt.bdryType    = 'pressure';

    opt = merge_options(opt, varargin{:});
    
    % inj rate in m3/s
    opt.inj_rate          = (opt.inj_rate_MtperYr * mega * kilo / rhoCref) / year; % m3/s
    

    %% Set-up (initial state, well location, schedule)
    
    % pressure and temperature
    initState.pressure  = Gt.cells.z * norm(gravity) * opt.water_density;   % hydrostatic pressure, in Pa=N/m^2
    initState.s         = repmat([1 0], Gt.cells.num, 1);                   % sat of water is 1, sat of CO2 is 0
    initState.sGmax     = initState.s(:,2);                                 % max sat of CO2 is initially 0
    initState.rs        = 0 * initState.sGmax;                              % initially 0   
    var.ref_p               = mean(initState.pressure);
    var.caprock_temperature = 273.15 + opt.seafloor_temp + ...
        (Gt.cells.z - opt.seafloor_depth) / 1e3 * opt.temp_gradient; % Kelvin

    
    
    
    % Set up injection:
    numWells = size(opt.wellCoords,1);
    injWells = [];
    offWells = [];
    for i = 1:numWells
        
        % Cell index(es) of injection point(s) in top grid formation
        var.wellCellIndex(i)  = getCellIndex(Gt, opt.wellCoords(i,1), opt.wellCoords(i,2));
        var.wellCoordSim(i,1) = Gt.cells.centroids(var.wellCellIndex(i),1);
        var.wellCoordSim(i,2) = Gt.cells.centroids(var.wellCellIndex(i),2);

        % Put injection controls into schedule
        schedule.control(1).W = addWell(injWells, Gt.parent, rock2D, var.wellCellIndex(i), ...
            'name',     sprintf('W%i', var.wellCellIndex(i)),  ...
            'Type',     'rate', ...
            'Val',      opt.inj_rate, ...
            'comp_i',   [0 1] );
        
        injWells = schedule.control(1).W;
        
        % Put turned-off well controls into schedule
        schedule.control(2).W  = addWell(offWells, Gt.parent, rock2D, var.wellCellIndex(i), ...
            'name',     'W_off',  ...
            'Type',     'rate', ...
            'Val',      0, ...
            'comp_i',   [0 1] );
        
        offWells = schedule.control(2).W;
        
    end

    % Set up boundary conditions, and put into schedule:
    bdryFaces   = find( Gt.faces.neighbors(:,1).*Gt.faces.neighbors(:,2) == 0 );
    bdryVal     = Gt.faces.z(bdryFaces) * opt.water_density * norm(gravity);
    bc          = addBC( [], bdryFaces, opt.bdryType, bdryVal, 'sat', [1 0] );
    for i = 1:numel(schedule.control)
        schedule.control(i).bc = bc;
    end
    
    % Set up time step, and put into schedule:
    dTi         = opt.inj_time / opt.inj_steps;
    dTm         = opt.mig_time / opt.mig_steps;
    istepvec    = ones(opt.inj_steps, 1) * dTi;
    mstepvec    = ones(opt.mig_steps, 1) * dTm;
    schedule.step.val       = [istepvec; mstepvec];
    schedule.step.control   = [ones(opt.inj_steps, 1); ones(opt.mig_steps, 1) * 2];
    
    

    %% Main parts of model and solver
    
    var.fluid = makeVEFluid(Gt, rock2D, 'sharp interface', ...
                                  'fixedT'      , var.caprock_temperature, ...
                                  'wat_mu_ref'  , opt.water_mu, ...
                                  'co2_mu_ref'  , opt.co2_mu, ...
                                  'wat_rho_ref' , opt.water_density, ...
                                  'co2_rho_ref' , opt.co2_density, ...
                                  'wat_rho_pvt' , [opt.water_compr_val, var.ref_p], ...
                                  'pvMult_p_ref', var.ref_p, ...
                                  'pvMult_fac'  , opt.pvMult, ...
                                  'residual'    , [opt.res_sat_wat, opt.res_sat_co2] , ...
                                  'dissolution' , opt.isDissOn, 'dis_max', opt.dis_max);
    var.model = CO2VEBlackOilTypeModel(Gt, rock2D, var.fluid);

    fprintf('\n\n Proceeding to solver... \n\n')
    [wellSols, states, sim_report] = simulateScheduleAD(initState, var.model, schedule);


end


function cellIndex = getCellIndex(Gt, Xcoord, Ycoord)
% Closest cell index of grid Gt corresponding to physical coordinate (X,Y)

    dv        = bsxfun(@minus, Gt.cells.centroids(:,1:2), [Xcoord, Ycoord]);
    [v, ind]  = min(sum(dv.^2, 2));
    cellIndex = ind; 
    % or Gt.cells.indexMap(i);
    
end

