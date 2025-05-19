%% Geothermal Example: Liquid-Vapor Phase Change (Boiling/Condensation)
% This example demonstrates simulation of phase change (liquid-vapor transition)
% in a geothermal reservoir using the enthalpy formulation.
%
% - 1D vertical column with hydrostatic pressure gradient
% - Initial state: subcooled liquid water
% - Top boundary: pressure/temperature set to induce boiling
% - Shows evolution of vapor region and enthalpy profile
%
% Requirements:
%   - MRST modules: geothermal, ad-core, ad-props, mrst-gui
%   - Run 'startup' in MRST root before executing this script
%
% Steps:
%   1. Add required modules
%   2. Set up 1D grid and rock/fluid properties
%   3. Define initial state (liquid water)
%   4. Apply boundary conditions to induce boiling
%   5. Run simulation and visualize phase change

%% Add necessary MRST modules
mrstModule add geothermal ad-core ad-props mrst-gui

%% Set up 1D vertical grid
nCells = 40;
height = 100; % meters
G = cartGrid([1, 1, nCells], [1, 1, height]);
G = computeGeometry(G);
rock = makeRock(G, 100*milli*darcy, 0.1);
rock = addThermalRockProps(rock);

%% Fluid: water with phase change (enthalpy formulation)
fluid = initSimpleADIFluid('phases', 'W', 'n', 1, 'mu', 1, 'rho', 1);
fluid = addThermalFluidProps(fluid, 'useEOS', true);

%% Model: use enthalpy formulation to allow phase change
model = GeothermalModel(G, rock, fluid, 'thermalFormulation', 'enthalpy');
model.minimumTemperature = 273.15; % K
model.maximumTemperature = 623.15; % K

%% Initial state: subcooled liquid water (no vapor)
p0 = 20*barsa; % bottom pressure
T0 = 300;      % K (below boiling)
state0 = initResSol(G, p0, 1);
state0.T = T0*ones(G.cells.num,1);

%% Boundary conditions: induce boiling at the top
fTop = find(G.faces.centroids(:,3) == max(G.faces.centroids(:,3)));
pTop = 10*barsa; % lower pressure at top
TTop = 450.15;   % K (boiling point at 1 atm)
bc = addBC([], fTop, 'pressure', pTop, 'sat', 1);
bc = addThermalBCProps(bc, 'T', TTop);

%% Time stepping
dt = repmat(10*day, 10, 1);
schedule = simpleSchedule(dt, 'bc', bc);

%% Simulate
[wellSols, states, reports] = simulateScheduleAD(state0, model, schedule);

%% Visualization: plot enthalpy and phase distribution
figure;
for t = 1:numel(states)
    subplot(1,2,1);
    plot(states{t}.T - 273.15, G.cells.centroids(:,3));
    xlabel('Temperature (C)'); ylabel('Depth (m)'); title('Temperature profile');
    set(gca, 'YDir', 'reverse');
    subplot(1,2,2);
    plot(states{t}.h{1}, G.cells.centroids(:,3));
    xlabel('Enthalpy (J/kg)'); ylabel('Depth (m)'); title('Enthalpy profile');
    set(gca, 'YDir', 'reverse');
    drawnow;
end
