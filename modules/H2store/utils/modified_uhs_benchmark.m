function [description, options, state0, model, schedule, plotOptions] = modified_uhs_benchmark(deck,varargin)
% Example modified from:
% Hogeweg, S., Strobel, G., & Hagemann, B. (2022). Benchmark study for the simulation of underground 
% hydrogen storage operations. Comput Geosci, 26, 1367–1378. 
% https://doi.org/10.1007/s10596-022-10163-5
%
% For further details on this case, refer to:
% Ahmed, E., et al. (2024). Phase behavior and black-oil simulations of hydrogen storage in saline aquifers. 
% Advances in Water Resources, 191, 104772.
%
% SEE ALSO:
%   `MRSTExample`, `example_template`, `exampleSuiteTutorial`
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
% Step 1: Test case description and options
% Description
description = ['Modified benchmark for hydrogen storage with ', ...
               'multiple injection/production cycles using the black-oil model'];

%% Optional input arguments: options corresponds to scenario 2 of the paper
K0 = 273.15 * Kelvin;
options = struct( ...
    'rateCharge',       6.14062 * 10^5 * meter^3/day/2,   ... % Charge rate (m³/day)
    'rateIdle',         0.0 * kilogram/day/2,             ... % Idle rate (kg/day)
    'rateCushion',      6.14062 * 10^5 * meter^3/day/2,   ... % Cushion (H₂) injection rate (m³/day)
    'rateDischarge',    6.14062 * 10^5 * meter^3/day/2,   ... % Discharge rate (m³/day)
    'bhp',              90.0 * barsa,                     ... % Production bottom hole pressure (BHP) in bar
    'tempCharge',       K0 + 60 * Kelvin,                 ... % Charging temperature (K)
    'tempDischarge',    K0 + 60 * Kelvin,                 ... % Discharging temperature (K)
    'tempCushion',      K0 + 60 * Kelvin,                 ... % Cushion temperature (K)
    'timeCushion',      120 * day,                        ... % Duration of cushioning phase (days)
    'timeCharge',       30 * day,                         ... % Duration of charging phase (days)
    'timeIdle',         15 * day,                         ... % Duration of idle phase (days)
    'timeShut',         30 * day,                         ... % Duration of shut phase (days)
    'timeDischarge',    30 * day,                         ... % Duration of discharging phase (days)
    'dtCharge',         3.0 * day,                        ... % Time step during charging (days)
    'dtCushion',        3.0 * day,                        ... % Time step during cushioning (H₂) (days)
    'dtIdle',           3.0 * day,                        ... % Time step during idle (days)
    'dtShut',           3.0 * day,                        ... % Time step during shut phase (days)
    'dtDischarge',      3.0 * day,                        ... % Time step during discharging (days)
    'numCycles',        20,                               ... % Total number of cycles (charging and discharging)
    'numCyclesCushions',6,                                ... % Number of cycles for cushion gas (H₂)
    'chargeOnly',       0,                                ... % Flag to simulate only the charging phase (0: No, 1: Yes)
    'cushionOnly',      0,                                ... % Flag to simulate only the cushioning phase (0: No, 1: Yes)
    'dischargeOnly',    0,                                ... % Flag to simulate only the discharging phase (0: No, 1: Yes)
    'useGroupCtrl',     false,                            ... % Flag to use group control (true/false)
    'initPres',         81.6 * barsa,                     ... % Initial pressure (bar)
    'initSat',          [1, 0],                           ... % Initial saturation
    'initTemp',         273.15 + 60,                      ... % Initial temperature (°C)
    'use_bc',           true,                             ... % Flag to use boundary conditions (true/false)
    'use_cushion',      true,                             ... % Flag to use cushion phase (true/false)
    'use_bhp',          true                              ... % Flag to use bottom hole pressure (true/false)
    );

%% Process optinal input arguments
[options, fullSetup, ~] = processTestCaseInput(mfilename, ...
    options, description, varargin{:});
options = checkOptions(options);
if ~fullSetup, return; end
%% Optional input arguments
options = merge_options(options, varargin{:});
if nargout <= 2, return;
end
% just to be sure on the module dependencies
require ad-core ad-props ad-blackoil
%% perm, poro, and grid as in the benchmark
G = processGRDECL(deck.GRID, 'checkgrid', false);
G = computeGeometry(G);
deck = convertDeckUnits(deck);
PERMX = deck.GRID.PERMX;
PERMY = deck.GRID.PERMY;
PERMZ = deck.GRID.PERMZ;
perm =  [PERMX PERMY PERMZ];
poro =  deck.GRID.PORO; clear p
rock = makeRock(G, perm, poro);

%% Set up schedule and initiate model from deck
schedule = setUpSchedule(G, rock, options);
[~, model, ~] = initEclipseProblemAD(deck,'getSchedule',false,'getInitialState', false);
fluid = model.fluid;
deck = model.inputdata;
rock.regions.saturation=deck.REGIONS.SATNUM;
gravity reset on;
%% We reset a blackoil model and "oil" is the water phase
model = GenericBlackOilModel(G, rock, fluid, 'water', false, 'disgas', true, 'inputdata', deck);
%% Set up initial state
state0 = setUpInitialState(model, schedule.control(1).W, options);

%% Plotting
plotOptions = {'View'              , [0,0]         , ...
    'PlotBoxAspectRatio', [1,1,0.25]    , ...
    'Projection'        , 'orthographic', ...
    'Size'              , [800, 300]    };
% deck.RUNSPEC.TITLE ='H2_illustration_storage';
% deck_new = model2Deck(model, schedule, 'deck', deck);
end

function W = setUpWells(G, rock, options)

    %% We also reset the well in the highest point: we reset Well coordinates (wc) and radii (r)
    % ps this is slightly different from the benchmark
    wc = [5891; 9612; 13333; 17054; 20775; 24496; 28217; 31938; 35659; 39380; 43101];
    r  = [0.9479; 0.7354; 0.7338; 2.2337; 2.2337; 2.6778; 2.6722; 5.1490; 5.1491; 0.7834; 0.7834];
    
    %% Add production wells with specified parameters
    W = addWell([], G, rock, wc, 'Name', 'Prod', ...
                'Radius', r, 'type', 'rate', ...
                'val', options.rateCharge, 'compi', [0, 1]);
    
    %% Set groups if group control is used
    if options.useGroupCtrl
        [W.group] = deal({'Inj', 'Prod'});
    end

end

function schedule = setUpSchedule(G0, rock, options)
    %% Set up initial wells based on provided grid, rock, fluid, and options
    W = setUpWells(G0, rock, options);     
    
    %% If cushion stage is being used
    if options.use_cushion
        if options.use_bhp
            % Set well properties for BHP control during cushion phase
            W(1).type     = 'bhp';
            W(1).name     = 'cushion';
            W(1).val      = options.bhp;
            W(1).T        = options.tempCushion;
            W(1).sign     = 1;

            %% Create time steps for cushion phase and set BHP value
            dtCushions = rampupTimestepsEnds(options.timeCushion, options.dtCushion);
            bhpCushion = options.bhp;

            %% Schedule for the initial 9 cushion steps
            for i = 1:9    
                dtCushion = dtCushions(i);
                W(1).val = bhpCushion;  % Fixed BHP value
                scheduleCushions{i} = simpleSchedule(dtCushion, 'W', W);
            end

            %% Intermediate cushion steps
            dtCushion = dtCushions(10:end-10);
            W(1).val = bhpCushion; 
            scheduleCushions{10} = simpleSchedule(dtCushion, 'W', W);

            %% Final cushion steps
            for i = 1:9    
                dtCushion = dtCushions(end-9+i);
                W(1).val = bhpCushion;  % Fixed BHP value
                scheduleCushions{i+10} = simpleSchedule(dtCushion, 'W', W);
            end
        else
            %% Set well properties for rate control during cushion phase
            W(1).type     = 'rate';
            W(1).name     = 'cushion';
            W(1).val      = options.rateCushion;
            W(1).T        = options.tempCushion;
            W(1).sign     = 1;

            %% Create time steps and rates for cushion phase
            dtCushions = rampupTimestepsEnds(options.timeCushion, options.dtCushion);
            rateCushion = options.rateCushion .* dtCushions ./ max(dtCushions);

            %% Schedule for the initial 9 cushion steps with rate control
            for i = 1:9
                dtCushion = dtCushions(i);
                W(1).val = rateCushion(10);
                scheduleCushions{i} = simpleSchedule(dtCushion, 'W', W);
            end

            %% Intermediate cushion steps with constant rate
            dtCushion = dtCushions(10:end-10);
            W(1).val = rateCushion(10); 
            scheduleCushions{10} = simpleSchedule(dtCushion, 'W', W);

            %% Final cushion steps with decreasing rates
            for i = 1:9
                dtCushion = dtCushions(end-9+i);
                W(1).val = rateCushion(end-9+i);  % Vary rate for last steps
                scheduleCushions{i+10} = simpleSchedule(dtCushion, 'W', W);
            end
        end
    end

    %% Set up wells for charge phase
    W = setUpWells(G0, rock, options);
    W(1).type = 'rate';
    W(1).name = 'charge';
    W(1).val  = options.rateCharge;
    W(1).T    = options.tempCharge;
    W(1).sign = 1;

    %% Create time steps and rates for charge phase
    dtCharges = rampupTimestepsEnds(options.timeCharge, options.dtCharge);
    rateCharge = options.rateCharge .* dtCharges ./ max(dtCharges);

    %% Schedule for the charge phase (initial 9 steps)
    for i = 1:9
        dtCharge = dtCharges(i);
        W(1).val = rateCharge(10);
        scheduleCharges{i} = simpleSchedule(dtCharge, 'W', W);
    end

    %% Intermediate charge steps with fixed rate
    dtCharge = dtCharges(10:end-10);
    W(1).val = rateCharge(10);
    scheduleCharges{10} = simpleSchedule(dtCharge, 'W', W);

    %% Final charge steps with constant rate
    for i = 1:9
        dtCharge = dtCharges(end-9+i);
        W(1).val = rateCharge(10);
        scheduleCharges{i+10} = simpleSchedule(dtCharge, 'W', W);
    end

    %% Handle group control if enabled
    if options.useGroupCtrl
        groups = [];
        scheduleCharge.groups = groups;
    end

    %% Set well properties for idle (shut) phase
    W(1).type = 'rate';
    W(1).val  = options.rateIdle;
    W(1).name = 'shut';
    W(1).T    = options.tempCushion;
    W(1).sign = 1;

    %% Create time steps for idle phase
    dtIdle = rampupTimestepsEnds(options.timeIdle, options.dtIdle);
    scheduleIdle = simpleSchedule(dtIdle, 'W', W);
    if options.useGroupCtrl
        groups = [];
        scheduleIdle.groups = groups;
    end

    %% Set up the shut phase
    dtShut = rampupTimestepsEnds(options.timeShut, options.dtShut);
    scheduleShut = simpleSchedule(dtShut, 'W', W);

    %% Set well properties for discharge phase
    W(1).type = 'rate';
    W(1).name = 'discharge';
    W(1).val  = -options.rateDischarge;
    W(1).sign = -1;
    W(1).lims.bhp = 20*barsa;
    W(1).cstatus(2:end) = 0;

    %% Create time steps for discharge phase
    dtDischarge = rampupTimestepsEnds(options.timeDischarge, options.dtDischarge);
    scheduleDischarge = simpleSchedule(dtDischarge, 'W', W);

    %% Handle group control for discharge
    if options.useGroupCtrl
        groups = [];
        scheduleCharge.groups = groups;
    end

    %% Combine schedules based on user options (charge, discharge, or cushion-only)
    if options.chargeOnly    
        schedule = combineSchedules(scheduleCharges{:}, 'makeConsistent', false);
    elseif options.dischargeOnly
        schedule = scheduleDischarge;
    elseif (options.cushionOnly && false)
        schedule = combineSchedules(scheduleCushions{:}, 'makeConsistent', false);
    else
        if options.use_cushion
            schedule = combineSchedules(scheduleCharges{:}, scheduleIdle, scheduleDischarge, 'makeConsistent', false);
            schedule = repmat({schedule}, 1, options.numCycles);

            %% Set different cushion values and combine
            scheduleCushion1 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false);
            for i = 1:length(scheduleCushions)
                scheduleCushions{i}.control.W.val = 94*barsa();
            end
            scheduleCushion2 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false);
            for i = 1:length(scheduleCushions)
                scheduleCushions{i}.control.W.val = 98*barsa();
            end
            scheduleCushion3 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false);
            for i = 1:length(scheduleCushions)
                scheduleCushions{i}.control.W.val = 102*barsa();
            end
            scheduleCushion4 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false);

            %% Final combination of schedules
            if options.cushionOnly
                schedule = combineSchedules(scheduleCushion1, scheduleCushion2, scheduleCushion3, scheduleCushion4, 'makeConsistent', false);
            else
                schedule = combineSchedules(scheduleCushion1, scheduleCushion2, scheduleCushion3, scheduleCushion4, scheduleCushion4, scheduleCushion4, schedule{:}, 'makeConsistent', false);
            end    
        else
            schedule = combineSchedules(scheduleCharges{:}, scheduleIdle, scheduleDischarge, 'makeConsistent', false);
            schedule = repmat({schedule}, 1, options.numCycles);
            schedule = combineSchedules(schedule{:}, 'makeConsistent', false);  
        end
    end

    %% Apply boundary conditions if enabled
    if options.use_bc
        bc = setUpBc(G0, options);
        for i = 1:numel(schedule.control)
            schedule.control(i).bc = bc;
        end
    end
end

function bc = setUpBc(G, options)
    %% Get boundary faces of the grid
    f = boundaryFaces(G);
    
    %% Select lateral faces
    dis = 2.48;
    f1 = (G.faces.normals(f, 3) > dis);
    faces = f(~f1);

    %% Calculate the displacement from a reference point (1140 units)
    Depth = 1140;
    dx = bsxfun(@minus, G.faces.centroids(faces, :), Depth);
    
    %% Calculate pressure drop due to gravity using omega and displacement
    rhoOS = 9.98e+02;
    dp = rhoOS .* (dx * reshape(gravity, [], 1));
    
    %% Define the pressure at the boundary (init pressure plus the gravity-induced pressure drop)
    pressure = options.initPres + dp;

    %% Set the boundary condition with the computed pressure and saturation
    bc = addBC([], faces, 'pressure', pressure, 'sat', options.initSat);    
end



% function schedule = setUpSchedule(G0, rock, options)
%         
%     W = setUpWells(G0, rock, options);     
%     if options.use_cushion
%        if options.use_bhp
%           W(1).type     = 'bhp';
%           W(1).name     = 'cushion';
%           W(1).val      = options.bhp;
%           W(1).T        = options.tempCushion;
%           W(1).sign     = 1;
%     
%           dtCushions       = rampupTimestepsEnds(options.timeCushion, options.dtCushion);
%           bhpCushion = [options.bhp];
%           for i = 1:9    
%               dtCushion      = dtCushions(i);
%               W(1).val      = bhpCushion(1);
%               scheduleCushions{i} = simpleSchedule(dtCushion, 'W', W);
%           end
% 
%           dtCushion       = dtCushions(10:end-10);
%           W(1).val      = bhpCushion(end);
%           scheduleCushions{10} = simpleSchedule(dtCushion, 'W', W);
% 
%        
%           for i = 1:9    
%               dtCushion      = dtCushions(end-9+i);
%               W(1).val      = bhpCushion(end);
%               scheduleCushions{i+10} = simpleSchedule(dtCushion, 'W', W);
%           end
%        else
%            
%            W(1).type     = 'rate';
%            W(1).name     = 'cushion';
%            W(1).val      = options.rateCushion;
%            W(1).T        = options.tempCushion;
%            W(1).sign     = 1;
%     
%            dtCushions       = rampupTimestepsEnds(options.timeCushion, options.dtCushion);
%            rateCushion = options.rateCushion.*dtCushions./max(dtCushions);
%            for i = 1:9               
%                dtCushion      = dtCushions(i);
%                W(1).val      = rateCushion(10);
%                scheduleCushions{i} = simpleSchedule(dtCushion, 'W', W);
%            end
% 
%            dtCushion       = dtCushions(10:end-10);
%            W(1).val      = rateCushion(10);
%            scheduleCushions{10} = simpleSchedule(dtCushion, 'W', W);
% 
%        
%        
%            for i = 1:9               
%                dtCushion      = dtCushions(end-9+i);
%                W(1).val      = rateCushion(end-9+i);
%                scheduleCushions{i+10} = simpleSchedule(dtCushion, 'W', W);
%            end
%        end
%     end
%    
%     W = setUpWells(G0, rock, options);    
%     W(1).type     = 'rate';
%     W(1).name     = 'charge';    
%     W(1).val      = options.rateCharge;
%     W(1).T        = options.tempCharge;
%     W(1).sign     = 1;
% %     W.lims.rate = options.rateCharge;
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
%     W(1).sign     = 1;
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
%     W(1).lims.bhp = 20*barsa;
%     W(1).cstatus(2:end) = 0;
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
%     elseif (options.cushionOnly&&false)
%         schedule = combineSchedules(scheduleCushions{:}, 'makeConsistent', false);
%     else 
%         if options.use_cushion
%            schedule = combineSchedules(scheduleCharges{:}, scheduleIdle, scheduleDischarge, 'makeConsistent', false);
%            schedule = repmat({schedule}, 1, options.numCycles);
%            scheduleCushion1 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false);
%            for i =1:length(scheduleCushions)
%               scheduleCushions{i}.control.W.val =94*barsa();
%            end
%            scheduleCushion2 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false); 
%             for i =1:length(scheduleCushions)
%               scheduleCushions{i}.control.W.val =98*barsa();
%            end
%            scheduleCushion3 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false);
%            for i =1:length(scheduleCushions)
%               scheduleCushions{i}.control.W.val =102*barsa();
%            end
%            scheduleCushion4 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false);  
%            if options.cushionOnly           
%                schedule = combineSchedules(scheduleCushion1,scheduleCushion2,scheduleCushion3,scheduleCushion4, 'makeConsistent', false);
%            else              
%                schedule = combineSchedules(scheduleCushion1,scheduleCushion2,scheduleCushion3,scheduleCushion4, scheduleCushion4, scheduleCushion4, schedule{:}, 'makeConsistent', false);
%            end    
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


%-------------------------------------------------------------------------%
function state0 = setUpInitialState(model, W, options)
    %% Initialize the reservoir solution based on the model grid and initial pressure
    state0 = initResSol(model.G, options.initPres, options.initSat);    
    state0.rs = zeros(size(state0.pressure));  % Set initial solution gas-to-oil ratio to zero
    %% Initialize well solutions
    wellSol = initWellSolAD(W, model, state0);
    wellSol.bhp = options.initPres;  % Set initial bottom hole pressure for all wells
    state0.wellSol = wellSol;  % Store well solutions in the state structure
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
        n = 1;
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