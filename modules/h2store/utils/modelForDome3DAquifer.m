function [description, options, state0, model, schedule,deck, plotOptions] = ModelForDome3DAquifer(deck, varargin)
% This function sets up a conceptual dome-shaped aquifer model for hydrogen storage, 
% focusing on the differences between immiscible and black-oil modeling approaches. 
%
% The reservoir domain measures 5000 × 3000 × 1000 meters with a depth of 600 m, 
% consisting of three geological layers: caprock, bedrock, and storage rock, 
% each characterized by distinct porosities and permeabilities. 
% The mesh is initially cartesian, with interior cell depths adjusted to form 
% an anticlinal structure, comprising a total of 83,232 cells arranged in a 
% [51 × 51 × 32] grid.
%
% The simulation features a 10-cycle operational schedule, where each cycle includes 
% charging, cushioning, idle, and discharging phases. The durations for these phases 
% are set to 10 days for charging, 3 days for idle periods, and 10 days for 
% discharging. Prior to the operational cycles, a 1-year buildup phase is implemented, 
% followed by a 90-day shut-in period to stabilize the reservoir pressure. 
%
% Injection and discharge rates are maintained at 400,000 m³/day, with the maximum 
% injection rate halved during the buildup phase to prevent overpressurization. 
% The rates gradually increase, initially reaching 200,000 m³/day. 
% PVT tables are generated for an initial reservoir pressure of 120 bar and a 
% temperature of 323.15 K.
%
% For further details on this case, refer to:
% Ahmed, E., et al. (2024). Phase behavior and black-oil simulations of hydrogen storage in saline aquifers. 
% Advances in Water Resources, 191, 104772.
%
% SEE ALSO:
%   `MRSTExample`, `example_template`, `exampleSuiteTutorial`.

%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

% Description
description = 'Conceptual model for hydrogen storage with multiple injection/production cycles';

% Optional input arguments
K0 = 273.15 * Kelvin;    % Reference temperature in Kelvin

% Define simulation options
options = struct( ...
    'rateCharge',       4.0e5 * meter^3 / day,    ... % Charge rate
    'rateIdle',         0.0 * kilogram / day,      ... % Idle rate
    'rateCushion',      2.0e5 * meter^3 / day,    ... % Cushion-H2 injection rate
    'rateDischarge',    4.0e5 * meter^3 / day,    ... % Discharge rate
    'bhp',              35.0 * barsa,              ... % Produced bottom hole pressure (BHP)
    'tempCharge',       K0 + 50 * Kelvin,          ... % Charge temperature
    'tempDischarge',    K0 + 50 * Kelvin,          ... % Discharge temperature
    'tempCushion',      K0 + 50 * Kelvin,          ... % Cushion temperature
    'timeCushion',      365 * day,                 ... % Cushion time
    'timeCharge',       10 * day,                  ... % Charging time
    'timeIdle',         3 * day,                   ... % Idle time
    'timeShut',         90 * day,                  ... % Shut-in time
    'timeDischarge',    10 * day,                  ... % Discharging time
    'dtCharge',         1*day,                ... % Time step during charging
    'dtCushion',        1*day * hour,                ... % Time step during cushioning
    'dtIdle',           8.4 * hour,                ... % Time step during idle
    'dtShut',           8.4 * hour,                ... % Time step during shut-in
    'dtDischarge',      1 * day,                ... % Time step during discharging
    'numCycles',        4,                        ... % Number of operational cycles
    'chargeOnly',       0,                         ... % Flag to simulate only the charging phase
    'cushionOnly',      0,                         ... % Flag to simulate only the cushioning phase
    'dischargeOnly',    0,                         ... % Flag to simulate only the discharging phase
    'useGroupCtrl',     false,                     ... % Use group control for injection/production
    'initPres',         37 * barsa,                ... % Initial reservoir pressure
    'initSat',          [1 0],                     ... % Initial reservoir saturation
    'initTemp',         273.15 + 50,               ... % Initial temperature
    'use_bc',           true,                      ... % Use boundary conditions
    'use_cushion',      true,                      ... % Use cushion phase
    'use_bhp',          false                      ... % Use bottom hole pressure control
    );

% Process optional input arguments
[options, fullSetup, ~] = processTestCaseInput(mfilename, options, description, varargin{:});
options = checkOptions(options); % Validate options

if ~fullSetup, return; end
require ad-core ad-props ad-blackoil spe10

%% Get the dome-shaped grid: low depth case with 40 bar at 600 meter depth
% case can be adjusted for high depth and hig pressure (see paper)
dome3DModel();

%% We set up the Rock properties
% porosity
rock = makeRock(G, [1000, 1000, 100]*milli * darcy, 0.25);
rock.poro(caprock)    = 0.1;    % Porosity of caprock
rock.poro(bedrock)    = 0.1;    % Porosity of bedrock
rock.poro(underrock)  = 0.15;   % Porosity of storage rock
rock.poro(toprock)    = 0.15;   % Porosity of top rock

% Layered permeability 
rock.perm(caprock, :)   = 1.0e-04*milli*darcy;   % Permeability of caprock
rock.perm(bedrock, :)   = 1.0e-02*milli*darcy;   % Permeability of bedrock
rock.perm(toprock, :)   = 1.0e-03*milli*darcy;   % Permeability of top rock
rock.perm(underrock, :) = 1.0e-03*milli*darcy; % Permeability of under rock

%% We set saturation regions
deck = convertDeckUnits(deck);
deck.GRID.ACTNUM = ones(G.cells.num, 1);       % Active cells
deck.REGIONS.SATNUM = repmat(1, G.cells.num, 1);  % Enforces size to match the number of cells
deck.REGIONS.SATNUM(bedrock)    = 2;           % Set SATNUM for bedrock
deck.REGIONS.SATNUM(toprock)    = 3;           % Set SATNUM for top rock
deck.REGIONS.SATNUM(caprock)    = 3;           % Set SATNUM for caprock
deck.REGIONS.SATNUM(underrock)  = 2;           % Set SATNUM for under rock

rock.regions.saturation = deck.REGIONS.SATNUM; 
G.cells.indexMap  = rock.regions.saturation;
deck.GRID.PORO    =  rock.poro;
deck.GRID.PERMX   = rock.perm(:,1);
deck.GRID.PERMY   = rock.perm(:,2);
deck.GRID.PERMZ   = rock.perm(:,3);

%% Deck conversion and initializationclose
[~, model, ~] = initEclipseProblemAD(deck,'G',G,'getSchedule',false,'getInitialState', false);
fluid = model.fluid;
deck = model.inputdata;

%% blackoil model setup
gravity reset on;
model = GenericBlackOilModel(G, rock, fluid, 'water', false, 'disgas', true, 'vapoil', true, ...
                             'gravity', [0, 0, -norm(gravity)], 'inputdata', deck);

%% Setup schedule
schedule = setUpSchedule(G, rock, fluid, options);

%% Setup initial state
state0 = setUpInitialState(model, schedule.control(1).W, options);

%% Plotting options
plotOptions = {'View', [0, 0], ...
               'PlotBoxAspectRatio', [1, 1, 0.25], ...
               'Projection', 'orthographic', ...
               'Size', [800, 300]};

% Uncomment the following line to set a title for the RUNSPEC in the deck
% deck.RUNSPEC.TITLE = 'H2_illustration_storage';

% Uncomment to convert the model back to a deck format
% deck_new = model2Deck(model, schedule, 'deck', deck);
end


function W = setUpWells(G, rock, fluid, options)
% setUpWells initializes the well structure for the reservoir model.
%
% This function adds production wells to the grid based on specified 
% criteria and sets up group control if applicable.
%
% Inputs:
%   G       - Grid structure containing cell information
%   rock    - Rock properties structure
%   fluid   - Fluid properties structure
%   options - Options structure containing simulation parameters
%
% Outputs:
%   W       - Well structure containing the initialized wells

    %% Initialize the well structure
    W = [];

    %% Identify the cell indices for production wells based on coordinates
    yd = 2500;
    zd = 600;
    wc = find(abs(G.cells.centroids(:, 1) - yd) < 0.45 & ...
              abs(G.cells.centroids(:, 2) - yd) < 0.2 & ...
              G.cells.centroids(:, 3) > -zd);
    wc = wc(1:2);
    %% Add  well to the well structure
    W = addWell(W, G, rock, wc, ...
                'Name', 'Prod', ...
                'Radius', 12.5 * centi * meter, ...
                'type', 'rate', ...
                'val', options.rateCharge, ...
                'compi', [0, 1]); % Composition indices for the well

    %% Set groups if group control is enabled
    if options.useGroupCtrl
        [W.group] = deal({'Inj', 'Prod'}); % Assign group names to the wells
    end
end


function bc = setUpBc(G, rock, fluid, options)
% setUpBc initializes boundary conditions for the reservoir model.
%
% This function sets the pressure boundary conditions and saturation 
% for the specified grid and fluid properties.
%
% Inputs:
%   G       - Grid structure containing cell and face information
%   rock    - Rock properties structure (not used in this function)
%   fluid   - Fluid properties structure
%   options - Options structure containing simulation parameters
%
% Outputs:
%   bc      - Boundary conditions structure containing pressure and saturation
    rhoOS = fluid.rhoOS;

    % Determine gravity vector
    if G.griddim < 3   
        grav_ = [0 -norm(gravity)]; % 2D gravity vector
    else
        grav_ = gravity; % 3D gravity vector
    end

    %% Identify boundary faces
    f = boundaryFaces(G);
    f1 = f(abs(G.faces.centroids(f, 1)) < eps);

    %% Calculate pressure gradient based on gravity
    ref_position = -100;
    scaling = 15;
    dx1 = bsxfun(@minus, G.faces.centroids(f1, :), ref_position) / scaling;
    dp1 = rhoOS .* (dx1 * reshape(grav_, [], 1));

    %% Set the hydrostatic pressure boundary condition
    pcmax = options.initPres + 0.5*barsa();
    pressure1 = pcmax - dp1;
    bc = addBC([], f1, 'pressure', pressure1, 'sat', options.initSat);

end


% 
% function schedule = setUpSchedule(G0, rock, fluid, options)
%         
%     W = setUpWells(G0, rock, fluid, options);     
%     if options.use_cushion
%        W(1).type     = 'rate';
%        W(1).name     = 'cushion';
%        W(1).val      = options.rateCushion;
%        W(1).T        = options.tempCushion;
%        W(1).sign     = 1;
%     
%        dtCushions       = rampupTimestepsEnds(options.timeCushion, options.dtCushion);
%        rateCushion = options.rateCushion.*dtCushions./max(dtCushions);
%        for i = 1:9    
%            dtCushion      = dtCushions(i);
%            W(1).val      = rateCushion(10);
%            scheduleCushions{i} = simpleSchedule(dtCushion, 'W', W);
%        end
% 
%        dtCushion       = dtCushions(10:end-10);
%        W(1).val      = rateCushion(10);
%        scheduleCushions{10} = simpleSchedule(dtCushion, 'W', W);
% 
%        
%        for i = 1:9    
%            dtCushion      = dtCushions(end-9+i);
%            W(1).val      = rateCushion(end-9+i);
%            scheduleCushions{i+10} = simpleSchedule(dtCushion, 'W', W);
%        end
%     end
%     W = setUpWells(G0, rock, fluid, options);    
%     W(1).type     = 'rate';
%     W(1).name     = 'charge';    
%     W(1).val      = options.rateCharge;
%     W(1).T        = options.tempCharge;
%     W(1).sign     = 1;
%     
%     dtCharges       = rampupTimestepsEnds(options.timeCharge, options.dtCharge);
%     rateCharge = options.rateCharge.*dtCharges./max(dtCharges);
%     for i = 1:9    
%         dtCharge      = dtCharges(i);
%         W(1).val      = rateCharge(10);
%         scheduleCharges{i} = simpleSchedule(dtCharge, 'W', W);
%     end
% 
%     dtCharge       = dtCharges(10:end-10);
%     W(1).val      = rateCharge(10);
%     scheduleCharges{10} = simpleSchedule(dtCharge, 'W', W);
% 
%        
%     for i = 1:9    
%         dtCharge      = dtCharges(end-9+i);
%         W(1).val      = rateCharge(10);
%         scheduleCharges{i+10} = simpleSchedule(dtCharge, 'W', W);
%     end
%     
%     
%     if options.useGroupCtrl
%         groups = [];
%         scheduleCharge.groups = groups;
%     end
% 
%     W(1).type     = 'rate';
%     W(1).val      = options.rateIdle;
%     W(1).name     = 'shut';        
%     W(1).T        = options.tempCushion;
%     W(1).sign     = -1;
%     
%     dtIdle       = rampupTimestepsEnds(options.timeIdle, options.dtIdle);
%     scheduleIdle = simpleSchedule(dtIdle, 'W', W);
%     if options.useGroupCtrl
%         groups = [];
%         scheduleIdle.groups = groups;
%     end
% 
%     dtShut       = rampupTimestepsEnds(options.timeShut, options.dtShut);
%     scheduleShut = simpleSchedule(dtShut, 'W', W);
% 
%     W(1).type     = 'rate';
%     W(1).name     = 'discharge';    
%     W(1).val      = -options.rateDischarge;
%     W(1).sign     = -1;
%     
%     dtDischarge       = rampupTimestepsEnds(options.timeDischarge, options.dtDischarge);
%     scheduleDischarge = simpleSchedule(dtDischarge, 'W', W);
%     if options.useGroupCtrl
%         groups = [];
%         scheduleCharge.groups = groups;
%     end
%     
%     if options.chargeOnly
%         schedule = combineSchedules(scheduleCharges{:}, 'makeConsistent', false);
%     elseif options.dischargeOnly
%         schedule = scheduleDischarge;
%     elseif options.cushionOnly
%         schedule = combineSchedules(scheduleCushions{:}, 'makeConsistent', false);
%     else 
%         if options.use_cushion
%            schedule = combineSchedules(scheduleCharges{:}, scheduleIdle, scheduleDischarge, 'makeConsistent', false);
%            schedule = repmat({schedule}, 1, options.numCycles);
%            schedule = combineSchedules(scheduleCushions{:}, scheduleShut, schedule{:}, 'makeConsistent', false);
%         else
%             schedule = combineSchedules(scheduleCharges{:}, scheduleIdle, scheduleDischarge, 'makeConsistent', false);
%             schedule = repmat({schedule}, 1, options.numCycles);
%             schedule = combineSchedules(schedule{:}, 'makeConsistent', false);  
%         end
%     end
% 
%         
%     if options.use_bc    
%         bc = setUpBc(G0,rock,fluid,options);        
%         for i = 1:numel(schedule.control)        
%             schedule.control(i).bc = bc;
%         end
%     end
% 
% end
function schedule = setUpSchedule(G0, rock, fluid, options)
% setUpSchedule configures the operational schedule for hydrogen storage 
% in the specified reservoir model.
%
% This function sets up various operational phases including charging, 
% cushioning, idling, and discharging, based on the provided options.
%
% Inputs:
%   G0      - Grid structure representing the reservoir model
%   rock    - Rock properties structure
%   fluid   - Fluid properties structure
%   options - Options structure containing simulation parameters
%
% Outputs:
%   schedule - Schedule structure defining the operational phases

    %% Initialize wells based on grid, rock, fluid, and options
    W = setUpWells(G0, rock, fluid, options);     
    
    %% Cushioning phase configuration
    if options.use_cushion
        W(1).type     = 'rate';
        W(1).name     = 'cushion';
        W(1).val      = options.rateCushion;
        W(1).T        = options.tempCushion;
        W(1).sign     = 1;

        %% Generate timesteps for the cushioning phase
        dtCushions = rampupTimestepsEnds(options.timeCushion, options.dtCushion);
        rateCushion = options.rateCushion .* dtCushions ./ max(dtCushions);
        
        %% Create schedules for the initial cushioning phases
        for i = 1:9    
            dtCushion = dtCushions(i);
            W(1).val  = rateCushion(10);
            scheduleCushions{i} = simpleSchedule(dtCushion, 'W', W);
        end

        %% Handle remaining cushioning phases
        dtCushion = dtCushions(10:end-10);
        W(1).val  = rateCushion(10);
        scheduleCushions{10} = simpleSchedule(dtCushion, 'W', W);
        
        %% Continue with final cushioning phases
        for i = 1:9    
            dtCushion = dtCushions(end-9+i);
            W(1).val  = rateCushion(end-9+i);
            scheduleCushions{i+10} = simpleSchedule(dtCushion, 'W', W);
        end
    end

    %% Set up wells again for the charging phase
    W = setUpWells(G0, rock, fluid, options);    
    W(1).type     = 'rate';
    W(1).name     = 'charge';    
    W(1).val      = options.rateCharge;
    W(1).T        = options.tempCharge;
    W(1).sign     = 1;
    
    %% Generate timesteps for the charging phase
    dtCharges = rampupTimestepsEnds(options.timeCharge, options.dtCharge);
    rateCharge = options.rateCharge .* dtCharges ./ max(dtCharges);
    
    %% Create schedules for the initial charging phases
    for i = 1:9    
        dtCharge = dtCharges(i);
        W(1).val  = rateCharge(10);
        scheduleCharges{i} = simpleSchedule(dtCharge, 'W', W);
    end

    %% Handle remaining charging phases
    dtCharge = dtCharges(10:end-10);
    W(1).val  = rateCharge(10);
    scheduleCharges{10} = simpleSchedule(dtCharge, 'W', W);

    %% Continue with final charging phases
    for i = 1:9    
        dtCharge = dtCharges(end-9+i);
        W(1).val  = rateCharge(10);
        scheduleCharges{i+10} = simpleSchedule(dtCharge, 'W', W);
    end
    
    %% Group control setup
    if options.useGroupCtrl
        groups = [];
        scheduleCharge.groups = groups;
    end

    %% Idle phase configuration
    W(1).type     = 'rate';
    W(1).val      = options.rateIdle;
    W(1).name     = 'shut';        
    W(1).T        = options.tempCushion;
    W(1).sign     = -1;
    
    %% Generate schedule for the idle phase
    dtIdle = rampupTimestepsEnds(options.timeIdle, options.dtIdle);
    scheduleIdle = simpleSchedule(dtIdle, 'W', W);
    if options.useGroupCtrl
        groups = [];
        scheduleIdle.groups = groups;
    end

    %% Shut-in phase configuration
    dtShut = rampupTimestepsEnds(options.timeShut, options.dtShut);
    scheduleShut = simpleSchedule(dtShut, 'W', W);

    %% Discharge phase configuration
    W(1).type     = 'rate';
    W(1).name     = 'discharge';    
    W(1).val      = -options.rateDischarge;
    W(1).sign     = -1;
    
    %% Generate schedule for the discharge phase
    dtDischarge = rampupTimestepsEnds(options.timeDischarge, options.dtDischarge);
    scheduleDischarge = simpleSchedule(dtDischarge, 'W', W);
    if options.useGroupCtrl
        groups = [];
        scheduleCharge.groups = groups;
    end
    
    %% Combine schedules based on the specified options
    if options.chargeOnly
        schedule = combineSchedules(scheduleCharges{:}, 'makeConsistent', false);
    elseif options.dischargeOnly
        schedule = scheduleDischarge;
    elseif options.cushionOnly
        schedule = combineSchedules(scheduleCushions{:}, 'makeConsistent', false);
    else 
        %% Combine all schedules for normal operation
        if options.use_cushion
            schedule = combineSchedules(scheduleCharges{:}, scheduleIdle, scheduleDischarge, 'makeConsistent', false);
            schedule = repmat({schedule}, 1, options.numCycles);
            schedule = combineSchedules(scheduleCushions{:}, scheduleShut, schedule{:}, 'makeConsistent', false);
        else
            schedule = combineSchedules(scheduleCharges{:}, scheduleIdle, scheduleDischarge, 'makeConsistent', false);
            schedule = repmat({schedule}, 1, options.numCycles);
            schedule = combineSchedules(schedule{:}, 'makeConsistent', false);  
        end
    end

    %% Setup boundary conditions if specified
    if options.use_bc    
        bc = setUpBc(G0, rock, fluid, options);        
        for i = 1:numel(schedule.control)        
            schedule.control(i).bc = bc;
        end
    end
end

%-------------------------------------------------------------------------%
function state0 = setUpInitialState(model, W, options)
    % setUpInitialState initializes the reservoir state for the simulation.
    %
    % Inputs:
    %   model - The model object containing grid and rock properties.
    %   W     - Well data structure containing well definitions.
    %   options - A structure containing simulation options.
    %
    % Output:
    %   state0 - Initial state structure with pressure, saturation, and well solutions.

    %% Initialize reservoir solution state with specified initial pressure and saturation
    state0 = initResSol(model.G, options.initPres, options.initSat);

    %% Set residual saturations to zero (not applicable for initial state)
    state0.rs = zeros(size(state0.pressure)); % Residual saturation for water
    state0.rv = zeros(size(state0.pressure)); % Residual saturation for gas

    %% Initialize well solutions based on well definitions
    wellSol = initWellSolAD(W, model, state0);
    
    %% Set initial bottom-hole pressure for wells
    [wellSol.bhp] = deal(options.initPres);
    
    %% Assign well solutions to the state structure
    state0.wellSol = wellSol;
end


%-------------------------------------------------------------------------%
function options = checkOptions(options)
    
    assert(~(options.chargeOnly && options.dischargeOnly), ...
        'Cannot simulate only charge and only discharge at the same time');
    
end

function dT = rampupTimestepsEnds(time, dt, n)
% Create timesteps that ramp up geometrically
%
% SYNOPSIS:
%   dT = rampupTimesteps(1*year, 30*day)
%   dT = rampupTimesteps(1*year, 30*day, 5)
%
% DESCRIPTION:
%   This function generates a timestep sequence for a given total time
%   interval that increases geometrically until it reaches some target
%   timestep. The rest of the interval is then divided into a number of
%   target timesteps.
%
% REQUIRED PARAMETERS:
%   time   - The total simulation time so that sum(dt) = time
%
%   dt     - Target timestep after initial ramp-up
%
%   n      - (OPTIONAL) Number of rampup steps. Defaults to 8.
%
% RETURNS:
%   dt     - Array of timesteps.
%
% NOTE:
%   The final timestep may be shorter than dt in order to exactly reach T.
%

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    if nargin < 3
        n = 8;
    end
    if time == 0
        dT = [];
        return
    end
    % Initial geometric series
    dt_init = (dt./2.^[n n:-1:1])';
    cs_time = cumsum(dt_init);
    if any(cs_time > time)
        dt_init = dt_init(cs_time < time);
    end
    
    % Remaining time that must be discretized
    dt_left = time - sum(2.*dt_init);
    % Even steps
    dt_rem = repmat(dt, floor(dt_left/dt), 1);
    % Final ministep if present
    dt_final = time - sum(2.*dt_init) - sum(dt_rem);
    % Less than to account for rounding errors leading to a very small
    % negative time-step.
    if dt_final <= 0
        dt_final = [];
    end
       
    if dt_final >= dt_init(1)
        dt_final = dt_init(1);
    end
    % Combined timesteps
    dT = [dt_init; dt_rem;sort(dt_init,'descend'); dt_final];
end
