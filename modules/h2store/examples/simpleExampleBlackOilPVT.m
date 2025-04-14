%% MRST Example: Blackoil Simulation with Tabulated PVT Data for H2-Water System
%
% This example demonstrates how to use the MRST blackoil simulator for
% hydrogen storage simulations with tabulated PVT data.
%
% The physical setup is a simple 2.5D domain. We begin by simulating an
% immiscible H2–water system and then enable miscibility by activating the
% DISGAS and/or VAPOIL options. The example illustrates how to build an MRST
% deck using tabulated data and how to set up solubility and relative
% permeability tables.
%
% Workflow:
%   1. Simulate an immiscible case (H2–water).
%   2. Generate and tabulate relative permeability (SGOF) tables.
%   3. Create solubility tables at specified pressure and temperature conditions.
%   4. Generate and tabulate PVT tables using RK or ePC-SAFT EOS.
%   5. Rebuild the MRST deck with all tables.
%   6. Simulate the miscible case by enabling DISGAS and/or VAPOIL.
%
% This example is focused on H2–water systems but can be extended to
% CO2–brine systems if the appropriate EOS is available for solubility tabulation.
%
% SEE ALSO:
%   `generateComponentTable`, `generateSolubilityTable`, `runDeckSimulation`
%
%
%{
Copyright 2009–2025 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

MRST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST. If not, see <http://www.gnu.org/licenses/>.
%}
%% MRST Example: Dome Grid with Tabulated PVT and Petrophysics
%
% This script sets up a simple 2.5D dome-shaped model with three geological 
% layers (underburden, reservoir, caprock). We define grid geometry, 
% petrophysical properties, and simple boundary conditions for use in black-oil 
% simulations of hydrogen storage.

%% Load Required MRST Modules
mrstModule add ad-core ad-blackoil ad-props mrst-gui spe10 deckformat

%% Initialize Gravity
gravity reset on
grav = [0, 0, -norm(gravity)];

%% Define Grid (2D X-Z Domain)
cartDims = [100, 1, 20];                  % Grid dimensions
physDims = cartDims .* [0.5, 1, 1.25];    % Physical dimensions in ft
G = cartGrid(cartDims, physDims);

%% Apply Dome-Shaped Deformation in Z-Direction
x = G.nodes.coords(:,1);
z = G.nodes.coords(:,3);
G.nodes.coords(:,3) = z + 10 * exp(-((x - physDims(1)/2).^2) / (2 * (physDims(1)/8)^2));
G.nodes.coords(:,3) = G.nodes.coords(:,3) - 35;  % Shift dome downward
G = computeGeometry(G);

%% Define Geological Regions
[~, ~, k] = ind2sub(cartDims, G.cells.indexMap);
underburden = k <= 6;
reservoir   = k > 6 & k <= 14;
caprock     = k > 14;

%% Plot Geological Layers
clf;
plotGrid(G, underburden, 'FaceColor', 'g'); hold on
plotGrid(G, reservoir,   'FaceColor', 'r');
plotGrid(G, caprock,     'FaceColor', 'b');
title('2D Dome Grid with 3 Geological Layers');
view(0, 180), axis tight

%% Petrophysical Properties
rock = makeRock(G, [10, 10, 10]*milli*darcy, 0.25);
rock.poro(caprock)      = 0.10;
rock.poro(underburden)  = 0.10;

% Layered permeability
rock.perm(caprock, :)     = 1.0e-04 * milli * darcy;
rock.perm(underburden, :) = 1.0e-02 * milli * darcy;

%% Assign Rock Regions
ACTNUM = ones(G.cells.num, 1);
SATNUM = ones(G.cells.num, 1);
SATNUM(caprock)      = 3;
SATNUM(underburden)  = 2;
rock.regions.saturation = SATNUM;
G.cells.indexMap = rock.regions.saturation;

%% Plot Log-Scale Permeability Field
clf;
plotCellData(G, log10(convertTo(rock.perm(:,1), milli*darcy)));
view(0, 180), axis tight;
title('Log-Scale Permeability Field')

%% Define Well: Injector at to Cell behind caprock
W = addWell([], G, rock, 1150, ...
    'Comp_i', [0, 1], ...
    'Type', 'rate', ...
    'Val', 2 * kilogram/day, ...
    'Radius', 5 * centi * meter, ...
    'Dir', 'z', ...
    'sign', +1, ...
    'Name', 'I');

%% Boundary Conditions: Hydrostatic Pressure on Lateral Faces
rhoOS = 998 * kilogram/meter^3;
boundFaces = boundaryFaces(G);
maxX = max(G.cells.centroids(:,1));
minX = min(G.cells.centroids(:,1));
lateralfaces = boundFaces( ...
    abs(G.faces.centroids(boundFaces,1)) < minX | ...
    abs(G.faces.centroids(boundFaces,1)) > maxX - eps);

z_face = G.faces.centroids(lateralfaces, 3);
z_ref  = max(G.nodes.coords(:,3));
dp_bc = rhoOS * norm(grav) * (z_face - z_ref);
pInit = 500 * psia;
sat = [1.0, 0.0];
pressure_bc = pInit - dp_bc;

bc = addBC([], lateralfaces, 'pressure', pressure_bc, 'sat', sat);

%% Step 1: Generate Relative Permeability and Capillary Pressure Tables
%
% This step builds region-dependent SGOF tables (relative permeability and capillary pressure)
% for a hydrogen-brine system using the `getFluidH2BrineSGOF` utility. 
% The output is formatted like ECLIPSE's 'SGOF' keyword format.
% For details, refer to: `examplePVTGenerationH2Brine`.

nreg = 3;  % Number of geological regions (e.g., underburden, reservoir, caprock)
outputPathPvt = fullfile(mrstOutputDirectory(), 'UHS_PVT_SIMPLE', 'PVT');
disp('Generating gas-oil flow properties for hydrogen-brine system with three regions...');

getFluidH2BrineSGOF( ...
    'n', 100, ...                                      % Number of saturation points
    'plot', false, ...                                  % Plot generated curves
    'outputPath', outputPathPvt, ...                   % Output directory
    'fileName', 'SGOFH2BRINE', ...                     % File name prefix
    'units', 'metric', ...                             % Use metric units
    'nreg', nreg, ...                                  % Number of regions
    'sw_imm',  [0.05, 0.1, 0.1], ...                    % Irreducible water saturation
    'sg_imm',  [0.05, 0.1, 0.1], ...                    % Residual gas saturation
    'c_a1',    3.5, ...                                 % Corey exponent for krw
    'c_a2',    3.5, ...                                 % Corey exponent for krg
    'c_a3',    1.5, ...                                 % Additional Corey exponent (if used)
    'Pe',      [0.4, 10, 10], ...                         % Entry pressures [Pa]
    'P_c_max', [9.5e4, 9.5e4, 9.5e4], ...               % Max capillary pressure [Pa]
    'reCompute', true);                                % Recompute tables even if they exist

% Read the generated SGOF tables
fid  = fopen(fullfile(outputPathPvt, 'SGOFH2BRINE'), 'r');
SGOF = readSGOFTables(fid, 'SGOF', nreg, 4);  % Read SGOF tables for nreg regions, 4 columns

%% Fluid Properties: H2–Brine System: We define an incompressible, immiscible two-phase fluid system consisting of:
%   - Water phase (representing brine)
%   - Gas phase (representing hydrogen)
%
% MRST's `initSimpleADIFluid` is used for this simplified setup. Relative permeability 
% and capillary pressure data are *not* included here, as they are handled separately 
% via region-dependent SGOF tables.
% Phase ordering: 'O' = water (brine), 'G' = gas (hydrogen)
initFluid = initSimpleADIFluid( ...
    'mu',    [1, 0.0087] * centi * poise, ...          % Viscosities: brine, hydrogen [Pa·s]
    'rho',   [998, 0.084] * kilogram / meter^3, ...    % Densities: brine, hydrogen [kg/m³]
    'c',     [4.14e-10, 8.1533e-3 / barsa], ...        % Compressibilities [1/Pa]
    'cR',    6.0e-4 / barsa, ...                       % Rock compressibility [1/Pa]
    'phases', 'OG');                                   % 'O' = brine, 'G' = hydrogen


% Replace synthetic rel. perm. and capillary pressure with SGOF data
immiscFluid = initfluid;

% Interpolate and assign SGOF-based relative permeability
fluid_kr = assignSGOF(immiscFluid, SGOF, struct( ...
    'sat', nreg, ...                        % Region number
    'interp1d', @interpTable));         % Use linear 1D interpolation

% Extract kr curves for the brine (krO) and hydrogen (krG) phases
immiscFluid.krO = fluid_kr.krOG{1};     % 'O' = water/brine phase
immiscFluid.krG = fluid_kr.krG{1};      % 'G' = gas/hydrogen phase


%% Define Reservoir Model
% Incompressible, immiscible H2/brine system using a two-phase oil/gas model
model = GenericBlackOilModel(G, rock, immiscFluid, ...
    'gravity', grav, 'disgas', false, 'vapoil', false, ...
    'water', false, 'oil', true, 'gas', true);
model.OutputStateFunctions{end+1}= 'ComponentPhaseMass';

%% Initialize Formation State
state0.pressure = pInit * ones(G.cells.num, 1);
state0.s = repmat(sat, G.cells.num, 1);

% Plot initial pressure
clf;
plotCellData(G, convertTo(state0.pressure, psia));
view(3), colorbar, axis tight, grid on;
xlabel('x'), zlabel('Depth'), title('Initial Pressure [Psi]');

%% Simulate Immiscible H₂/Brine Injection
% We simulate a first scenario involving immiscible hydrogen injection into a brine-saturated formation.

% Define simple injection schedule
timesteps = [0.01, 0.025, 0.075, 0.1, repmat(0.5, [1, 6]), ...
             ones(1, 6), repmat(2, [1, 13]), repmat(3, [1, 59])] * day;

% Set up schedule with wells and boundary conditions
schedule = simpleSchedule(timesteps, 'W', W, 'bc', bc);

% Simulate and visualize
clf;
fn = getPlotAfterStep(state0, model, schedule, ...
                      'view', [0, 180], 'field', 's:2');

[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
                                                'afterStepFn', fn);

%% Generate PVT Data and Solubility Tables for Miscible Model
% We simulate a miscible hydrogen-brine model, tabulating solubility and generating PVT data.

% Set parameters for temperature, pressure, and salinity
min_temp = 60; max_temp = 61; min_pressure = 30 * atm; max_pressure = 100 * barsa(); nbp = 100; nbt = 10; ms = 0;

% Define output directories
outputPathSol = fullfile(mrstOutputDirectory(), 'UHS_PVT_SIMPLE', 'H2SolubilityTable');
outputPathPvt = fullfile(mrstOutputDirectory(), 'UHS_PVT_SIMPLE', 'PVT');

% Ensure output directories exist
if ~exist(outputPathSol, 'dir'), mkdir(outputPathSol); end
if ~exist(outputPathPvt, 'dir'), mkdir(outputPathPvt); end

% Generate solubility table for H2O and H2
disp('Generating component table for H2O');
tab_H2O = generateComponentProperties('min_temp', min_temp, 'max_temp', max_temp, 'n_temp', nbt, 'min_press', min_pressure, 'max_press', max_pressure, 'n_press', nbp, 'comp_name', 'H2O', 'outputDisplay', false, 'outputPath', outputPathSol);

pause(0.5);  % Ensure smooth execution between commands
disp('Generating component table for H2');
tab_H2 = generateComponentProperties('min_temp', min_temp, 'max_temp', max_temp, 'n_temp', nbt, 'min_press', min_pressure, 'max_press', max_pressure, 'n_press', nbp, 'comp_name', 'H2', 'outputDisplay', false, 'outputPath', outputPathSol);

% Generate solubility table for H2O/H2 in brine
pressures = reshape(tab_H2.("pressure [Pa]"), [], nbt);
temperatures = reshape(tab_H2.("# temperature [°C]"), [], nbt) + 273.15;  % Convert to Kelvin
tab_sol = generateH2WaterSolubilityTable('min_temp', min_temp, 'max_temp', max_temp, 'n_temp', nbt, 'min_press', min_pressure, 'max_press', max_pressure, 'n_press', nbp, 'ms', ms, 'outputDisplay', false, 'outputPath', outputPathSol, 'reCompute', true);

% Generate PVT table
getFluidH2BrineProps(tab_H2O, tab_H2, tab_sol, 'rs', true, 'rv', true, 'PVTGFile', 'PVTGH2BRINE', 'PVTOFile', 'PVTOH2BRINE', 'outputPath', outputPathPvt, 'reCompute', true, 'plot', false);

%% Pass PVT Tables and Setup Eclipse Deck for Miscible Model
% After generating PVT tables, we pass them through the Eclipse deck for simulation. This ensures proper integration of PVT data for the miscible model.
VAPOIL = true;
DISGAS = true;
% We reset our miscible fluid
model.fluid.krOG = fluid_kr.krOG{1};
model.fluid.krG = fluid_kr.krG{1}; clear fluid_kr
deck = model2Deck(model, schedule);

% Add regions and PVT data
deck.GRID.ACTNUM = SATNUM;        % Active cells
deck.REGIONS.SATNUM = SATNUM;     % Region size must match number of cells

% Read and assign PVT data from files
if DISGAS
    fid = fopen(fullfile(outputPathPvt, 'PVTOH2BRINE'), 'r');
    pvto = readMisciblePVTTable(fid, 3, 3, 'PVTO');
    deck.PROPS.PVTO{1} = pvto{2};
    deck.RUNSPEC.DISGAS = true;
    % Add rs to state0
    state0.rs = state0.pressure.*0;
end

if VAPOIL
    fid = fopen(fullfile(outputPathPvt, 'PVTGH2BRINE'), 'r');
    pvtg = readMisciblePVTTable(fid, 3, 3, 'PVTG');
    deck.PROPS.PVTG{1} = pvtg{2};
    deck.RUNSPEC.VAPOIL = true;
    % Add rv to state0
    state0.rv = state0.pressure.*0;
end

% Convert deck units
deck = rmfield(deck, 'PCUNIT');
deck.RUNSPEC = rmfield(deck.RUNSPEC, 'SI');
deck.RUNSPEC.METRIC = 1;
deck = convertDeckUnits(deck);

% Initialize Eclipse problem
[~, modeldeck, ~] = initEclipseProblemAD(deck, 'G', G, 'getSchedule', false, 'getInitialState', false);
% % let define  model with VAPOIL and DISGAS
miscibModel = GenericBlackOilModel(G, model.rock, modeldeck.fluid, 'gravity', grav, 'disgas', DISGAS,...
     'vapoil', VAPOIL, 'water', false, 'inputdata',modeldeck.inputdata);
miscibModel.OutputStateFunctions{end+1}= 'ComponentPhaseMass';
 close all;
fnRs = getPlotAfterStep(state0, miscibModel, schedule, 'view', [0, 180], ...
                       'field', 's:2');

 [wellSolsRs, statesRs, reportRs] = ...
    simulateScheduleAD(state0, miscibModel, schedule,'afterStepFn', fnRs);
disp(['Successful blackoil simulation with tabulated PVT data: disgas = ', ...
      string(DISGAS), ', vapoil = ', string(VAPOIL)]);
disp('Tip: Try setting vapoil = false to observe the impact of disabling vaporization.');
 % <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
