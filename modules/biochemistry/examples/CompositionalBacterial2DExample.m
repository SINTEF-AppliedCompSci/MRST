clear all;
mrstModule add compositional ad-blackoil ad-core ad-props mrst-gui
gravity reset on

%% Read the Eclipse deck file containing the simulation data
% Change input fil by UHS_BENCHMARK_RS_SALT.DATA for SALT EFFECTS
%deck = readEclipseDeck('/home/elyes/Documents/mrst-2023b/spe11-utils/deck_H2/UHS_benchmark/UHSDATA/UHS_BENCHMARK_RS.DATA');
deck = readEclipseDeck('/home/elyes/Documents/Projects/MRST/modules/H2store/data/Illustrative_example/H2STORAGE_RS.DATA');

%% Set up the simulation parameters and model components
[~, options, state0, model, schedule, ~] = H2_illustration_storage_example(deck);



%state0.components= [0.8, 0.0, 0.006, 0.018, 0.176];
% Define compositional fluid model (with CoolProp library support)
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Nitrogen', 'Methane'}, ...
                                      {'H2O', 'H2', 'CO2', 'N2', 'CH4'});


%%=========initialisation proprietes du fluide
%===Compositional fluid model (initialization with the CoolProp library)=

% %Brine-Gas (H2)
% [rhow,rhog]=deal(1000* kilogram/meter^3,8.1688* kilogram/meter^3); %density kilogram/meter^3;
% [viscow,viscog]=deal(1.0*centi*poise,0.0094234*centi*poise);%viscosity
% [cfw,cfg]=deal(0,8.1533e-3/barsa); %compressibility
%  
% [srw,src]=deal(0.0,0.0);
% Phydro0=rhow*norm(gravity).*G.cells.centroids(:,3);
% [Pmaxz,Pref1,Pminz,Pe]=deal(95*barsa,114*barsa,120*barsa,0.1*barsa); %pressions
% 
% %initialisation fluides incompressibles, Brooks-Corey relperm krw=(Sw)^nw
% fluid=initSimpleADIFluid('phases', 'OG', 'mu',[viscow,viscog],...
%                          'rho',[rhow,rhog],'pRef',Pref1,...
%                          'c',[cfw,cfg],'n',[2,2],'smin',[srw,src]);
% 
% % Pression capillaire
% pcOG = @(so) Pe * so.^(-1/2);
% fluid.pcOG = @(sg) pcOG(max((1-sg-srw)./(1-srw), 1e-5)); %@@
% 
% Fluid density and viscosity (kg/m^3 and cP)
[rhow, rhog] = deal(1000 * kilogram / meter^3, 8.1688 * kilogram / meter^3);
[viscow, viscog] = deal(1.0 * centi * poise, 0.0094234 * centi * poise);

% Compressibility (per bar)
[cfw, cfg] = deal(0, 8.1533e-3 / barsa);

% Relative permeability and initial saturations
[srw, src] = deal(0.0, 0.0);
fluid = initSimpleADIFluid('phases', 'OG', 'mu', [viscow, viscog], ...
                           'rho', [rhow, rhog], 'pRef', 114 * barsa, ...
                           'c', [cfw, cfg], 'n', [2, 2], 'smin', [srw, src]);
model.fluid = fluid;
% T = 100*day;
% pv=sum(poreVolume(G,rock))/T;
% rate = 100*pv;
% niter = 30;

% %%==============Conditions aux limites et puits============
% bc = [];
% %puit d'injection
% W = [];
% W = verticalWell(W, G, rock,1,1,nz, 'comp_i', [0, 1],'Radius',0.5,...
%     'name', 'Injector', 'type', 'rate','Val',rate, 'sign', 1);
% W(1).components = [0.0 0.95 0.05 0.0 0.0];
% % for i = 1:numel(W)
% %     W(i).components = info.injection;
% % end
%eos = initDeckEOSModel(deck);
% schedule = convertDeckScheduleToMRST(model, deck);   
% schedule.step.val = schedule.step.val.*day;
for i=1:length(schedule.control)
    schedule.control(i).W.compi=[0, 1];
    schedule.control(i).W.components  = [0.0, 0.95, 0.05, 0.0, 0.0];
    schedule.control(i).W.T = options.initTemp+40;
end
%%==============model compositionnal================
arg = {model.G, model.rock, model.fluid, compFluid,...
    'water', false, 'oil', true, 'gas', true,... % water-oil system
	'bacteriamodel', true,'diffusioneffect',false,'liquidPhase', 'O',...
    'vaporPhase', 'G'}; % water=liquid, gas=vapor
model = BiochemistryModel(arg{:});
nbact0 = 10^6;

state0 = initCompositionalStateBacteria(model, state0.pressure, options.initTemp+40, [0.8, 0.2], [0.8, 0.0, 0.006, 0.018, 0.176], nbact0);
model.outputFluxes = false;
% %===Conditions initiales=====================
% T0=317.5;
% s0= [0.8 0.2]; %initial saturations  Sw=1
% z0 = [0.8,0.0,0.006,0.018,0.176]; %initial composition: H2O,H2,CO2,N2,CH4.
% 
% %===Bacteria model===========================
% if model.bacteriamodel
%     nbact0=10^6;
%     state0 = initCompositionalStateBacteria(model,Phydro0,T0,s0,z0,nbact0);
% else
%     state0 = initCompositionalState(model, Phydro0, T0, s0, z0);
% end
% 
% 
% %===Ajout d'un terme source====================
% src=[];

%===Resolution pression/transport========================
% deltaT = T/niter;
% schedule = simpleSchedule(repmat(deltaT,1,niter),'bc', bc,'src', src,'W',W);
% nls = NonLinearSolver('useRelaxation', true);
%% Set up the linear and nonlinear solvers
lsolve = selectLinearSolverAD(model);                          % Select the linear solver for the model
nls = NonLinearSolver();                                       % Create a nonlinear solver object
nls.LinearSolver = lsolve;                                     % Assign the linear solver to the nonlinear solver

name = 'UHS_2DCASE_COMPOSITIONAL_BACT';
%% Pack the simulation problem with the initial state, model, and schedule
problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);

%% Run the simulation
simulatePackedProblem(problem,'restartStep',71);
%% gGet reservoir and well states
[ws,states] = getPackedSimulatorOutput(problem);
%===Plottings============================================
time=0;
figure;
for i= 1:niter
    x = G.cells.centroids(:,1);
    z = G.cells.centroids(:,3);
    X = reshape(x, [nx,nz]);
    Z = reshape(z, [nx,nz]);
    zH2 = reshape(states{i}.components(:,2), [nx,nz]);
    zH2O = reshape(states{i}.components(:,1), [nx,nz]);
    zCO2 = reshape(states{i}.components(:,3), [nx,nz]);
    Sw = reshape(states{i}.s(:,1), [nx,nz]);
    Sg = reshape(states{i}.s(:,2), [nx,nz]);
    Pres= reshape(states{i}.pressure, [nx,nz]);
    nbacteria=reshape(states{i}.nbact, [nx,nz]);

    subplot(2,2,1);   
    contourf(X,Z,Sw,60,'EdgeColor','auto');
    clim([0 1])
    axis equal
    axis ([0 Lx  depth_res depth_res+H])
    ylabel('Water saturation','FontSize',15)
    set(gca,'Ydir','reverse')
    colormap('jet')
    colorbar 

    subplot(2,2,2);   
    contourf(X,Z,Sg,60,'EdgeColor','auto');
    clim([0 1])
    axis equal
    axis ([0 Lx  depth_res depth_res+H])
    ylabel('Gas saturation','FontSize',15)
    set(gca,'Ydir','reverse')
    colormap('jet')
    colorbar 

    subplot(2,2,3); 
    contourf(X,Z,nbacteria,60,'EdgeColor','auto');
    axis equal
    axis ([0 Lx  depth_res depth_res+H])
    ylabel('nbact','FontSize',15)
    set(gca,'Ydir','reverse')
    colormap('jet')
    colorbar

    subplot(2,2,4); 
    contourf(X,Z,zH2,60,'EdgeColor','auto');
    clim([0 0.8])
    axis equal
    axis ([0 Lx  depth_res depth_res+H])
    ylabel('z_H2','FontSize',15)
    set(gca,'Ydir','reverse')
    colormap('jet')
    colorbar
   
    time = time + deltaT;
    title(sprintf('injection duration = %.2f days',convertTo(time,day)))
    pause(0.001)
end

    function [description, options, state0, model, schedule, plotOptions] = H2_illustration_storage_example_bacteria(deck,varargin)
% This example simulates the injection and behavior of hydrogen (H₂) in a 2D saline aquifer 
% using the black-oil model. The aquifer has a dome-shaped structure defined by the function 
% F(x) = σ + r*sin(π*x), with parameters σ = 25 and r = 5. 
% 
% The domain is 50 m × 50 m, with caprock and bedrock layers providing containment through 
% high entry (capillary) pressure and low permeability. 
% 
% A single well is placed at the top of the trap, in the uppermost cell beneath the caprock. 
% We apply fixed hydrostatic pressure at the lateral boundaries and enforce no-flux conditions 
% at the top and bottom boundaries.
%
% Relative permeability and capillary pressure follow the Brooks-Corey model, with varying 
% residual saturations for different rock types.
%
% SEE ALSO:
%   `MRSTExample`, `example_template`, `exampleSuiteTutorial`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify it under the terms of the GNU 
General Public License as published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

MRST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with MRST. If not, 
see <http://www.gnu.org/licenses/>.
%}
% Step 1: Test Case Description and Options
%---------------------------------------------------------------------%
% Description of the conceptual model for hydrogen storage with multiple
% injection/production (Inj/Prod) cycles.

description = 'Conceptual model for Hydrogen storage with multiple Inj/Prod cycles';

% Time units
K0    = 273.15 * Kelvin;  % Absolute temperature offset
% Optional input arguments
options = struct( ...
    'rateCharge'   , 18 * kilogram/day   , ... % Hydrogen injection rate during charging
    'rateIdle'     , 0.0 * kilogram/day  , ... % No injection during idle periods
    'rateCushion'  , 10 * kilogram/day   , ... % Cushion gas injection rate (H₂)
    'rateDischarge', 18 * kilogram/day   , ... % Hydrogen production rate during discharge
    'bhp'          , 35.0 * barsa        , ... % Bottom hole pressure during production
    'tempCharge'   , K0 + 40 * Kelvin    , ... % Injection temperature during charging
    'tempDischarge', K0 + 40 * Kelvin    , ... % Production temperature during discharge
    'tempCushion'  , K0 + 40 * Kelvin    , ... % Temperature for cushion gas phase
    'timeCushion'  , 90 * day            , ... % Duration of cushion gas phase
    'timeCharge'   , 30 * day            , ... % Duration of hydrogen charging
    'timeIdle'     , 10 * day            , ... % Idle period between charge/discharge cycles
    'timeShut'     , 30 * day            , ... % Shut-in period (no activity)
    'timeDischarge', 30 * day            , ... % Duration of hydrogen production (discharge)
    'dtCharge'     , 8.4 * hour          , ... % Timestep during charging
    'dtCushion'    , 8.4 * hour          , ... % Timestep during cushion phase
    'dtIdle'       , 8.4 * hour          , ... % Timestep during idle phase
    'dtShut'       , 8.4 * hour          , ... % Timestep during shut-in phase
    'dtDischarge'  , 8.4 * hour          , ... % Timestep during discharge
    'numCycles'    , 10                  , ... % Number of injection/production cycles
    'chargeOnly'   , 0                   , ... % Simulate only charging period
    'cushionOnly'  , 0                   , ... % Simulate only cushion gas phase
    'dischargeOnly', 0                   , ... % Simulate only discharge period
    'useGroupCtrl' , false               , ... % Group control for wells (optional)
    'initPres'     , 37 * barsa          , ... % Initial reservoir pressure
    'initTemp'     , K0 + 40             , ... % Initial reservoir temperature
    'initSat'      ,[1 0]                , ... % Initial reservoir saturation
    'use_bc'       , true                , ... % Use boundary conditions
    'use_cushion'  , true                , ... % Include cushion gas in the simulation
    'use_bhp'      , false                 ... % Use bottom hole pressure control
);

% Process optional input arguments
[options, fullSetup, ~] = processTestCaseInput(mfilename, ...
options, description, varargin{:});  % Process test case inputs
options = checkOptions(options);          % Check the validity of options

if ~fullSetup
    return;  % If setup is incomplete, exit early
end
% Merge optional input arguments with the existing options
options = merge_options(options, varargin{:});
% If the number of output arguments is less than or equal to 2, return early
if nargout <= 2
    return;
end

% Define module dependencies for the simulation
require ad-core ad-props ad-blackoil spe10 upr

% Grid setup
% Generate a constraint path along the x-axis and define the grid shape
x = linspace(0.0, 1.0);       % Generate linearly spaced points for x-axis
y = 25 + 5 * sin(pi * x);     % Define y-axis using a sinusoidal function for dome shape
w = {[50 .* x', y']};         % Combine x and y to form 2D points for grid constraints

% Grid size scaling and dimensions
gS = [0.5, 0.5];              % Scaling factors for grid
pdims = [50, 50];             % Grid dimensions in x and y directions

% Generate two versions of a composite PEBI grid
G1 = compositePebiGrid2D(gS, pdims, 'cellConstraints', w, ...
    'interpolateCC', true, 'protLayer', false);  % Without protection layer

% G2 = compositePebiGrid2D(gS, pdims, 'cellConstraints', w, ...
%     'interpolateCC', true, 'protLayer', true);   % With protection layer

% Compute geometry for the grid
G = computeGeometry(G1);

  
% Slight adjustment to the y-axis values for accuracy
y = y + 0.1;  % Adjust y-values to ensure correct positioning above the grid
% Create a polygon representation of the line
linePolygon = [[0, 50]; [50 .* x', y']; [50, 50]];  % Define the polygon that represents the dome structure

% Define the cell vertices from the grid
cellVertices = G.cells.centroids;  % Extract centroids of the grid cells

% Use inpolygon to find which cells are above the line
aboveLineMask = inpolygon(cellVertices(:, 1), cellVertices(:, 2), linePolygon(:, 1), linePolygon(:, 2)); 

% Create rock properties
rock = makeRock(G, [10 * milli * darcy, 10 * milli * darcy], 0.25);  % Define rock permeability and porosity
caprock = aboveLineMask;  % Assign mask for caprock
bedrock = find(G.cells.centroids(:, 2) < 5);  % Identify bedrock cells based on their centroid positions

% Set permeability for caprock and bedrock
rock.perm(caprock, :) = 1.0e-04 * milli * darcy;  % Caprock permeability
rock.perm(bedrock, :) = 1.0e-02 * milli * darcy;  % Bedrock permeability

% Set porosity for caprock and bedrock
rock.poro(caprock) = 0.1;  % Caprock porosity
rock.poro(bedrock) = 0.1;  % Bedrock porosity

% Convert deck units for simulation
deck = convertDeckUnits(deck);
% Initialize activity and saturation numbers in the deck
deck.GRID.ACTNUM = ones(G.cells.num, 1);  % Set active cells
deck.REGIONS.SATNUM = ones(G.cells.num, 1);  % Set saturation numbers to 1

% Update saturation numbers for specific regions
deck.REGIONS.SATNUM(bedrock) = 2;  % Bedrock saturation number
deck.REGIONS.SATNUM(caprock) = 3;  % Caprock saturation number
rock.regions.saturation = deck.REGIONS.SATNUM;  % Assign saturation regions to rock properties
% We update indexmap and rock in G
G.cells.indexMap = rock.regions.saturation;
deck.GRID.PORO=rock.poro;
deck.GRID.PERMX=rock.perm(:,1);
deck.GRID.PERMY=rock.perm(:,2);
deck.GRID = rmfield(deck.GRID,'PERMZ');
% Initialize Eclipse problem for AD (Automatic Differentiation)
[~, model, ~] = initEclipseProblemAD(deck,'getSchedule',false,'G',G,'getInitialState', false);  

% Extract fluid and input data from the model
fluid = model.fluid;  
deck = model.inputdata;  


% Reset gravity for the model
gravity reset on;

% Create the black-oil model with gravity effects and input data
model = GenericBlackOilModel(G, rock, fluid, ...
    'water', false, 'disgas', true, 'gravity', [0 -norm(gravity)], 'inputdata', deck);

% Set up the simulation schedule based on the grid, rock properties, fluid model, and options
schedule = setUpSchedule(G, rock, fluid, options);

% Set up the initial state for the simulation using the first control well from the schedule
state0 = setUpInitialState(model, schedule.control(1).W, options);

% Define plotting options for visualization
plotOptions = {'View'              , [0,0]         , ...
    'PlotBoxAspectRatio', [1,1,0.25]    , ...
    'Projection'        , 'orthographic', ...
    'Size'              , [800, 300]    };
% Optional: Set title for the simulation run in the deck (commented out)
% deck.RUNSPEC.TITLE = 'H2_illustration_storage';

% Optional: Convert the model to a deck format (commented out)
% deck_new = model2Deck(model, schedule, 'deck', deck);

end

function W = setUpWells(G, rock, fluid, options)
    % setUpWells - Sets up the wells for the simulation
    % 
    % Syntax: W = setUpWells(G, rock, fluid, options)
    %
    % Inputs:
    %   G      - Grid structure
    %   rock   - Rock properties structure
    %   fluid  - Fluid properties structure
    %   options - Simulation options structure
    %
    % Outputs:
    %   W      - Well structure containing the defined wells

    % Initialize the well array
    W = [];  
    
    % Identify cells for the production well based on specified coordinates
    wc = find(abs(G.cells.centroids(:, 1) - 25) < 1e-1 & ...
              abs(G.cells.centroids(:, 2)) > 25 & ...
              abs(G.cells.centroids(:, 2)) < 29.5);
    
    % Add a production well at the identified cells with specified properties
    W = addWell(W, G, rock, wc, ...
                'Name', 'Prod', ...                       % Well name
                'Radius', 5 * centi * meter, ...        % Well radius
                'Type', 'rate', ...                      % Well type (rate control)
                'Val', options.rateCharge, ...           % Production rate
                'Compi', [0, 1]);                        % Component indices

    % Set well groups if group control is enabled
    if options.useGroupCtrl
        [W.group] = deal({'Inj', 'Prod'});         % Assign groups for injection and production
    end
end


function bc = setUpBc(G, rock, fluid, options)
    % setUpBc - Sets up boundary conditions for the simulation.
    %
    % Syntax: bc = setUpBc(G, rock, fluid, options)
    %
    % Inputs:
    %   G      - Grid structure
    %   rock   - Rock properties structure
    %   fluid  - Fluid properties structure
    %   options - Simulation options structure
    %
    % Outputs:
    %   bc     - Boundary condition structure

    % Initialize saturation and fluid density
    sat = options.initSat;
    omega = fluid.rhoOS;  % Water saturation density

    % Determine gravity vector based on grid dimensions
    if G.griddim < 3
        grav_ = [0 -norm(gravity)];  % 2D case
    else
        grav_ = gravity;  % 3D case
    end

    % Identify boundary faces
    f = boundaryFaces(G);
    f1 = f(abs(G.faces.centroids(f, 1)) < eps);     % Left boundary
    f2 = f(abs(G.faces.centroids(f, 1) - 50) < eps); % Right boundary

    % Calculate pressure difference for left boundary
    dx1 = bsxfun(@minus, G.faces.centroids(f1, :), 0);
    dp1 = omega .* (dx1 * reshape(grav_, [], 1));  % Pressure change due to gravity
    p0 = options.initPres;                           % Initial pressure
    pcmax = p0 + 3 * barsa();                       % Maximum pressure condition
    pressure1 = pcmax + dp1;                           % Total pressure for left boundary

    % Add boundary condition for left boundary
    bc = addBC([], f1, 'pressure', pressure1, 'sat', sat);

    % Calculate pressure difference for right boundary
    dx2 = bsxfun(@minus, G.faces.centroids(f2, :), 0);
    dp2 = omega .* (dx2 * reshape(grav_, [], 1));  % Pressure change due to gravity
    pressure2 = pcmax + dp2;                       % Total pressure for right boundary

    % Add boundary condition for right boundary
    bc = addBC(bc, f2, 'pressure', pressure2, 'sat', sat);
end



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
    % setUpSchedule - Sets up the well schedule for the reservoir simulation.
    %
    % Syntax: schedule = setUpSchedule(G0, rock, fluid, options)
    %
    % Inputs:
    %   G0      - Grid structure
    %   rock    - Rock properties structure
    %   fluid   - Fluid properties structure
    %   options - Simulation options structure
    %
    % Outputs:
    %   schedule - Schedule structure containing well control information

    % Set up wells
    W = setUpWells(G0, rock, fluid, options); 

    % Define cushion schedule if applicable
    if options.use_cushion
        W(1).type = 'rate';
        W(1).name = 'cushion';
        W(1).val = options.rateCushion;
        W(1).T = options.tempCushion;
        W(1).sign = 1;

        % Ramp up cushion rates
        dtCushions = rampupTimestepsEnds(options.timeCushion, options.dtCushion);
        rateCushion = options.rateCushion .* dtCushions ./ max(dtCushions);
        
        % Create cushion schedules
        for i = 1:9    
            dtCushion = dtCushions(i);
            W(1).val = rateCushion(10);  % Use the tenth value for the first nine steps
            scheduleCushions{i} = simpleSchedule(dtCushion, 'W', W);
        end

        % Handle remaining cushion schedules
        dtCushion = dtCushions(10:end-10);
        W(1).val = rateCushion(10);
        scheduleCushions{10} = simpleSchedule(dtCushion, 'W', W);

        for i = 1:9    
            dtCushion = dtCushions(end-9+i);
            W(1).val = rateCushion(end-9+i);
            scheduleCushions{i+10} = simpleSchedule(dtCushion, 'W', W);
        end
    end

    % Define charge schedule
    W = setUpWells(G0, rock, fluid, options); 
    W(1).type = 'rate';
    W(1).name = 'charge';    
    W(1).val = options.rateCharge;
    W(1).T = options.tempCharge;
    W(1).sign = 1;

    % Ramp up charge rates
    dtCharges = rampupTimestepsEnds(options.timeCharge, options.dtCharge);
    rateCharge = options.rateCharge .* dtCharges ./ max(dtCharges);
    
    % Create charge schedules
    for i = 1:9    
        dtCharge = dtCharges(i);
        W(1).val = rateCharge(10);  % Use the tenth value for the first nine steps
        scheduleCharges{i} = simpleSchedule(dtCharge, 'W', W);
    end

    % Handle remaining charge schedules
    dtCharge = dtCharges(10:end-10);
    W(1).val = rateCharge(10);
    scheduleCharges{10} = simpleSchedule(dtCharge, 'W', W);

    for i = 1:9    
        dtCharge = dtCharges(end-9+i);
        W(1).val = rateCharge(end-9+i);
        scheduleCharges{i+10} = simpleSchedule(dtCharge, 'W', W);
    end

    % Set up idle schedule
    W(1).type = 'rate';
    W(1).val = options.rateIdle;
    W(1).name = 'shut';        
    W(1).T = options.tempCushion;
    W(1).sign = -1;

    dtIdle = rampupTimestepsEnds(options.timeIdle, options.dtIdle);
    scheduleIdle = simpleSchedule(dtIdle, 'W', W);

    % Set up shut schedule
    dtShut = rampupTimestepsEnds(options.timeShut, options.dtShut);
    scheduleShut = simpleSchedule(dtShut, 'W', W);

    % Define discharge schedule
    W(1).type = 'rate';
    W(1).name = 'discharge';    
    W(1).val = -options.rateDischarge;
    W(1).sign = -1;

    dtDischarge = rampupTimestepsEnds(options.timeDischarge, options.dtDischarge);
    scheduleDischarge = simpleSchedule(dtDischarge, 'W', W);

    % Combine schedules based on options
    if options.chargeOnly
        schedule = combineSchedules(scheduleCharges{:}, 'makeConsistent', false);
    elseif options.dischargeOnly
        schedule = scheduleDischarge;
    elseif options.cushionOnly
        schedule = combineSchedules(scheduleCushions{:}, 'makeConsistent', false);
    else 
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

    % Set up boundary conditions if specified
    if options.use_bc    
        bc = setUpBc(G0, rock, fluid, options);        
        for i = 1:numel(schedule.control)        
            schedule.control(i).bc = bc;  % Assign boundary conditions to each control
        end
    end
end

function state0 = setUpInitialState(model, W, options)
    % setUpInitialState - Initializes the state of the reservoir simulation.
    %
    % Syntax: state0 = setUpInitialState(model, W, options)
    %
    % Inputs:
    %   model  - Simulation model structure containing grid and rock properties
    %   W      - Well structure
    %   options - Simulation options structure
    %
    % Outputs:
    %   state0 - Initial state structure containing pressure and saturation

    % Initialize reservoir solution with specified initial pressure and default saturations
    state0 = initResSol(model.G, options.initPres, options.initSat);  % [1, 0] represents initial saturations

    % Initialize residual saturations (set to zero initially)
    state0.rs = 0 .* state0.pressure;

    % Uncomment and modify the following lines to adjust saturations based on bedrock properties
    % bedrock = find(model.G.cells.centroids(:, 2) < 15);
    % state0.s(bedrock, 1) = 1 - (model.rock.poro(bedrock, 1) .* model.G.cells.centroids(bedrock, 2) ./ 5);
    % state0.s(bedrock, 2) = (model.rock.poro(bedrock, 1) .* model.G.cells.centroids(bedrock, 2) ./ 5);

    % Initialize well solutions
    wellSol = initWellSolAD(W, model, state0);
    wellSol.bhp = options.initPres;  % Set initial bottom hole pressure for wells
    state0.wellSol = wellSol;         % Assign well solutions to the state

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