%% Example Simulation From Tenth SPE CSP
%% We start from the original immiscible setup. We show how to set up misccible cases based on 
%% tabutaled PVT data. We generate H2 solubility tables and then PVT tables, refdefine the MRST model to consider ..

clear all;
close all;
mrstModule add ad-core ad-blackoil ad-props mrst-gui spe10
gravity reset on

%% Model Geometry
% Grid is a 100-by-1-by-20 Cartesian box with equally sized cells of
% dimensions 25-by-25-by-2.5 feet.
cartDims = [ 100, 1, 20 ];
physDims = cartDims .* [ 25, 5, 2.5 ]*ft;
G = cartGrid(cartDims, physDims);
G = computeGeometry(G);

%% Petrophysical Properties
% Porosity is constant (=0.2) throughout the formation.  The permeability
% distribution is a correlated geostatistically generated field stored in a
% file supplied by the SPE.

%rock = getSPE10_model_1_rock();
rock = makeRock(G, [500, 500, 1000]*milli * darcy, 0.25);
clf
plotCellData(G, log10(convertTo(rock.perm(:,1), milli*darcy)))
view(3), axis tight, grid on

%% Define Sources and Sinks
% Model is produced from a single producer located at the far end of the
% model (I==100) constrained at a bottom-hole pressure of 95 Psi.  There is
% a single injector at the near end (I==1) providing pressure support.  The
% injector fills the reservoir with 6.97 cubic metres of gas per day.  Both
% wells have an internal diameter of 1 ft.

IJK = gridLogicalIndices(G);
I   = find(IJK{1}(:,1) == 1);
P   = find(IJK{1}(:,1) == G.cartDims(1));  clear IJK

W = addWell([], G, rock, 1450, 'Comp_i', [ 0, 1 ], 'Type', 'rate', ...
            'Val', 1*meter^3/day, 'Radius', 1*ft, 'Dir', 'z', ...
            'sign', +1, 'Name', 'I', 'refDepth', 0*ft);

%% Official Benchmark Relative Permeability Data
% Build a reduced ECLIPSE-style input deck that contains just enough
% information to construct relative permeability curves based on the
% official benchmark data.  In particular we use the fact that the relative
% permeability data is formatted in the same way as ECLIPSE's 'SGOF'
% keyword data.
kr_deck = getSPE10_model_1_relperm();

clf
plot(kr_deck.PROPS.SGOF{1}(:, 1),            ...
     kr_deck.PROPS.SGOF{1}(:, [2, 3]), '*-', ...
     'LineWidth', 2, 'MarkerSize', 5)
legend({'kr_g', 'kr_o'}, 'Location', 'Best')
xlabel('S_g'), title('Relative Permeability, Model I 10th CSP')

%% Fluid Properties
% The fluids in this simulation model are incompressible and immiscible
% with constant viscosities.  This means we can use MRST's special purpose
% fluid constructor |initSimpleADIFluid| to create the fluid object.  We
% will however need to use sampled relative permeability curves so we do
% not enter any relative permeability data in this call.
initfluid = initSimpleADIFluid('mu'    , [1, 0.0087]*centi*poise, ...
                           'rho'   , [998, 0.084]*kilogram/meter^3, ...
                           'c'   , [4.14e-10, 8.1533e-3/barsa], ...
                           'cR'    , 6.0e-4/barsa, ...
                           'phases', 'OG');

%%
% Replace the synthetic relative permeability curves created through
% function |initSimpleADIFluid| with the real benchmark values.
immiscfluid = initfluid;
fluid_kr = assignSGOF(immiscfluid, kr_deck.PROPS.SGOF, struct('sat', 1, ...
                                               'interp1d', @interpTable));
immiscfluid.krG = fluid_kr.krG{1};
immiscfluid.krO = fluid_kr.krOG{1};

%%
% The <matlab:mrstModule('add','spe10') SPE 10 module> contains the special
% purpose function |getSPE10_model_1_fluid| that performes the above fluid
% manipulations so one would generally not do this manually.

%% Form Reservoir Model
% This is an incompressible, immiscible oil/gas system.
model = GenericBlackOilModel(G, rock, immiscfluid, 'gravity', gravity, 'disgas', false,...
    'vapoil', false, 'water', false, 'oil', true, 'gas', true);
%% Initialise Formation
% Formation is initially filled with oil and the initial pressure at the
% top of the model is 100 Psi.
region = getInitializationRegionsBlackOil(model, 0, 'datum_pressure', 500*psia);
state0 = initStateBlackOilAD(model, region);

clf
plotCellData(G, convertTo(state0.pressure, psia))
view(3), colorbar(), axis tight, grid on
xlabel('x'), zlabel('Depth'), title('Initial Pressure Distribution [Psi]')

%%
% Ten years (3650*day), ramp-up time-steps.
timesteps = [ 0.1, 0.2, 0.3, 0.4, repmat(0.5, [1, 6]), ones([1, 6]), ...
              repmat(4, [1, 5]), repmat(6, [1, 8]), repmat(8, [1, 9]), ...
              repmat(8, [1, 5]), repmat(8, [1, 50]), ...
              repmat(8, [1, 38]) ]*day;

% Set up the schedule containing both the wells and the timesteps
schedule = simpleSchedule(timesteps, 'W', W);
%%
fn = getPlotAfterStep(state0, model, schedule, 'view', [0, 0], ...
                      'field', 's:2', 'wells', W);
hold on;
plotWell(G, W, 'Color', 'k', 'FontSize', 10);
[wellSols, states, report] = ...
         simulateScheduleAD(state0, model, schedule, 'afterStepFn', fn);

%% Now we want to tabulate PVT data and experiment MISCIBLE model.
% Here is easy to tabulate the solubility table using our Henry correlation
min_temp = 15;                % Minimum temperature in Celsius
max_temp = 16;                % Maximum temperature in Celsius
min_pressure = 30 * atm;       % Minimum pressure in Pa
max_pressure = 200 * barsa();      % Maximum pressure in Pa
nbp = 100;                     % Number of pressure points
nbt = 1;                     % Number of temperature points
ms = 0;                       % we chose to not consider salinity


% Define the target output directory relative to the current directory
outputPathSol = fullfile(mrstOutputDirectory(), 'Tenth_SPE_CSP', 'H2SolubilityTable');
% Check if the directory exists; if not, create it
if ~exist(outputPathSol, 'dir')
    mkdir(outputPathSol);
end
% Generate H2O Component Table
comp_name = 'H2O';
disp(['Generating component table for: ', comp_name]);
tab_H2O = generateComponentProperties('min_temp',min_temp, 'max_temp',max_temp, 'n_temp', nbt, 'min_press',min_pressure, 'max_press', max_pressure, 'n_press',nbp, 'comp_name', comp_name,'outputDisplay', false,'outputPath',outputPathSol);
% Generate H2 Component Table
pause(0.5);  % Ensure smooth execution between commands
comp_name = 'H2';
disp(['Generating component table for: ', comp_name]);
tab_H2 = generateComponentProperties('min_temp',min_temp, 'max_temp',max_temp, 'n_temp', nbt, 'min_press',min_pressure, 'max_press', max_pressure, 'n_press',nbp, 'comp_name', comp_name,'outputDisplay', false,'outputPath',outputPathSol);
pressures = reshape(tab_H2.("pressure [Pa]"), [], nbt);
temperatures = reshape(tab_H2.("# temperature [°C]"), [], nbt) + 273.15;  % Convert to Kelvin
% Generate solubility tables using ..
tab_sol = HenrySetschenowH2BrineEos(temperatures, 0, pressures);
tab_sol = generateH2WaterSolubilityTable('min_temp',min_temp, 'max_temp',max_temp, 'n_temp', nbt,'min_press',min_pressure, 'max_press',max_pressure, 'n_press', nbp, 'ms', ms,'outputDisplay', false,'outputPath',outputPathSol,'reCompute', true);
%% Generate PVT table
% Define the target output directory relative to the current directory
outputPathPvt = fullfile(mrstOutputDirectory(), 'Tenth_SPE_CSP', 'PVT');
% Check if the directory exists; if not, create it
if ~exist(outputPathPvt, 'dir')
    mkdir(outputPathPvt);
end
getFluidH2BrineProps(tab_H2O, tab_H2, tab_sol,'rs', true, 'rv', true,'PVTGFile', 'PVTGH2BRINE', 'PVTOFile', 'PVTOH2BRINE', 'PVTGFile', 'PVTGH2BRINE','outputPath', outputPathPvt, 'reCompute', true);

% We need a to pass our PVT tables later
miscibfluid = initfluid;
model.fluid.krG = fluid_kr.krG{1};
model.fluid.krOG = fluid_kr.krOG{1};              clear fluid_kr
deck = model2Deck(model,schedule);

fid=fopen(fullfile(outputPathPvt,'PVTOH2BRINE'),'r');
kwo=readMisciblePVTTable(fid, 3, 3, 'PVTO');
deck.PROPS.PVTO{1} = kwo{2};
fid=fopen(fullfile(outputPathPvt,'PVTGH2BRINE'),'r');
kwg=readMisciblePVTTable(fid, 3, 3, 'PVTG');
% deck.PROPS.PVTG{1} = kwg{2};
deck.RUNSPEC.DISGAS = true;
deck.RUNSPEC.VAPOIL = false;
deck =rmfield(deck,'PCUNIT');
deck.RUNSPEC =rmfield(deck.RUNSPEC,'SI');
deck.RUNSPEC.METRIC =1;
deck = convertDeckUnits(deck);
[~, modeldeck, ~] = initEclipseProblemAD(deck,'G',G,'getSchedule',false,'getInitialState', false);
%PVTG  TABLE are computed in metric
% we need to convert them to SI
% u = unitConversionFactors('metric', 'si');
% % first is the pressure
% % --PRES 
% deck.PROPS.PVTO{1}.data(:,1) = deck.PROPS.PVTO{1}.data(:,1)*u.press;
% % --VFV
% deck.PROPS.PVTO{1}.data(:,2) = deck.PROPS.PVTO{1}.data(:,2).*u.volume;
% % --VISC
% deck.PROPS.PVTO{1}.data(:,3) = deck.PROPS.PVTO{1}.data(:,3).*u.viscosity;
% 
% % deck.PROPS.PVTG{1}.data(:,1) = deck.PROPS.PVTG{1}.data(:,1)*1000;
% % --VFV
% deck.PROPS.PVTG{1}.data(:,2) = deck.PROPS.PVTG{1}.data(:,2).*u.volume;
% % --VISC
% deck.PROPS.PVTG{1}.data(:,3) = deck.PROPS.PVTG{1}.data(:,3).*u.viscosity;
% deck.PROPS.PVTG{1}.key = deck.PROPS.PVTG{1}.key.*u.press;
% deck.PROPS.PVTG{1}.data(:,2) = deck.PROPS.PVTG{1}.data(:,2)./10*3;
% deck.RUNSPEC.VAPOIL = true;
% deck.PROPS.SGOF = sgof;
% deck.REGIONS.SATNUM = ones(G.cells.num, 1);  % Set saturation numbers to 1
% % rock.regions.saturation = deck.REGIONS.SATNUM;  % Assign saturation regions to rock properties

fluid = modeldeck.fluid;
deck = modeldeck.inputdata;
% let define first model with dissolved gas
modelRs = GenericBlackOilModel(G, rock, fluid, 'gravity', gravity, 'disgas', true,...
    'vapoil', false, 'water', false, 'inputdata',deck);
% 
fnRs = getPlotAfterStep(state0, modelRs, schedule, 'view', [0, 0], ...
                      'field', 's:1', 'wells', W);
state0.rs = state0.pressure.*0+0.0;
state0.rv = state0.pressure.*0+0.0;

 [wellSolsRs, statesRs, reportRs] = ...
    simulateScheduleAD(state0, modelRs, schedule,'afterStepFn', fnRs);
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
