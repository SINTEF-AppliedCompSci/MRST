%% VE simulation with capillary fringe using black-oil solver
% In this example, we will run a few simulations of injection and migration of
% CO_2 into the Gassum formation.  We will look at the impact of boundary
% conditions on pressure buildup, as well as of including a capillary fringe.
% This example employs the new, class-based framework for setting up
% simulations.  This framework was introduced into MRST-co2lab in the 2015a
% release.  Previous to that release, fully-implicit simulations were carried
% out using the approach demonstrated in
% <http://www.sintef.no/projectweb/mrst/modules/co2lab/tutorials/ve-simulation-in-a-standard-black-oil-solver/
% this tutorial>.

% The following line ensures we have the necessary modules loaded
moduleCheck co2lab ad-core

% Ensure gravity is in effect
gravity on

% Setting initial parameters
opt.sf_temp = 7;    % degrees centigrade
opt.sf_depth = 100; % seafloor depth in meters
opt.tgrad = 35;     % degrees per kilometer depth
opt.wellpos = [7.55e5, 6.39e6]; % injection site coordinates
opt.annual_injection = 6 * mega * 1e3; % inject six megatons per year
opt.injection_time = 50 * year;
opt.migration_time = 1000 * year;
opt.injection_steps = 10;
opt.migration_steps = 20;

%% Load simulation grid and rock structure
coarsening_level = 3; % grid downsampling to speed up simulation/visualization
[Gt, rock2D] = getFormationTopGrid('Gassumfm', coarsening_level);


%% Setup fluid objects (with and without capillary fringe)

% Aquifer temperature field (necessary to compute local CO_2 density)
T = 273.15 + opt.sf_temp + (Gt.cells.z - opt.sf_depth) / 1e3 * opt.tgrad;

% Sharp-interface type fluid
fluid_si = makeVEFluid(Gt, rock2D, 'sharp interface', 'fixedT', T);

% Fluid with linear capillary fringe
fluid_cf = makeVEFluid(Gt, rock2D, 'P-scaled table', 'fixedT', T);

%% Setup initial state in the aquifer

% We use water density at reference conditions ('rhoWS') to compute
% (approximate) hydrostatic aquifer pressure.
initState.pressure = Gt.cells.z * norm(gravity) * fluid_si.rhoWS;
initState.s = repmat([1 0], Gt.cells.num, 1); % [water sat, CO2 sat]
initState.sGmax = initState.s(:,2); % historically maximum CO2 saturation


%% Setup boundary conditions, injection well and schedule

% Identifying well cell
dist2 = bsxfun(@minus, Gt.cells.centroids, opt.wellpos);
dist2 = sum(dist.^2,2); % squared distance cell centre to inj. point
[~, wcell_ix] = min(dist2); % index of cell closest to inj. point

% Setup well 
injection_rate = opt.annual_injection / year / fluid.rhoGS; % volumetric
W = addWell([], Gt, rock2D, wcell_ix , ...
            'type'   , 'rate'        , ...
            'val'    , injection_rate, ...
            'comp_i' , [0 1]);

% Specify boundary conditions
bfaces = find(any(Gt.faces.neighbors==0, 2)); % identify boundary faces
open_bc = addBC([], bfaces, ...
                'pressure', Gt.faces.z(bfaces) * fluid_si.rhoWS * norm(gravity), ...
                'sat', [1 0]);

W_off = ;
istep = ;
mstep = ;

% Define schedule with closed boundaries

schedule_cb.control = [struct('W', W, 'bc', []), ...   % well during injection
                       struct('W', W_off, 'bc', [])];  % well during migration
schedule_cb.step    = struct('control', [1 * ones(size(istep)); ...
                                         2 * ones(size(mstep))], ...
                             'val', [istep; mstep]);

% Define schedule with open boundaries

schedule_ob = schedule_cb;
schedule_ob.control(1).bc = open_bc; % open boundaries during injection
scehdule_ob.control(2).bc = open_bc; % open boundaries during migration


%% Setup simulation models
model_si = CO2VEBlackOilTypeModel(Gt, rock2D, fluid_si);
model_cf = CO2VEBlackOilTypeModel(Gt, rock2D, fluid_cf);

%% Run simulation with open boundaries and sharp interface
[wellSols_si, states_si] = simulateScheduleAD(initState, model_si, schedule_ob);

%% Run simulation with open boundaries and capillary fringe
[wellSols_cf, states_cf] = simulateScheduleAD(initState, model_cf, schedule_ob);

%% Run simulation with closed boundaries
[wellSols_cb, states_cb] = simulateScheduleAD(initState, model_si, schedule_cb);