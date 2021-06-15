%% A simple, essentially 1D, validation of the hybrid VE solvers
% Here, we solve a simple problem without layering. We set up fine-scale,
% coarse scale, VE and hybrid VE solvers (with non-VE regions) and compare
% the results.
%
% For more information see:
% Møyner, O., & Nilsen, H. M. (2019). Multiresolution coupled vertical
% equilibrium model for fast flexible simulation of CO2 storage.
% Computational Geosciences, 23(1), 1-20.

%% Add modules
gravity reset on;
mrstModule add hybrid-ve co2lab;
mrstModule add ad-core ad-blackoil ad-props;
mrstModule add matlab_bgl coarsegrid;
mrstModule add mrst-gui;
%% Select linear or non-linear version of case
% By default, we use non-linear fluxes. You can specify the variables
% "useNonLinearFlux = false" to simulate the alternative scenario from the
% paper.
if ~exist('useNonLinearFlux', 'var')
    % Define if not found
    useNonLinearFlux = true;
end
if useNonLinearFlux
    % Setup non-linear relperm curves. We use quadratic and cubic
    % Corey-exponents for CO2 and water, respectively
    n_w = 3;
    n_g = 2;
    name = 'simple_displacement_nonlinear';
else
    % Use the linear relative permeability curves. This makes the
    % non-VE coarse-scale solution much better, as there is no upscaling
    % error from the mobility.
    n_w = 1;
    n_g = 1;
    name = 'simple_displacement_linear';
end
%% Set up solvers
% We set up the four different problems that we want to simulate and
% create a cell array containing the packed problem.
%
% We generate the following simple displacement solvers:
% 1) A conventional solver on the fine mesh (reference)
% 2) VE-hybrid with fine cells in parts of the domain
% 3) VE in the entire domain
% 4) A conventional solver on the coarse mesh

%% 1) Fine scale model
% Setup initial fine scale grid and compute geometry. Grid is one cell
% thick in the y direction, therefore essentially 2D. This grid will also
% be used later on when creating the VE grids.

ny = 1;
nx = 100;
nz = 100;

L = 10000;
G = cartGrid([nx, ny, nz], [L, 1, 10]);
G = computeGeometry(G);

% Add rock and fluid properties
rock = makeRock(G, 300*milli*darcy, 0.3);
mu = [0.30860, 0.056641];

fluid = initSimpleADIFluid('rho',     [975.86, 686.54]*kilogram/meter^3, ...
                           'n',       [n_w, n_g], ...
                           'phases',  'WG',...
                           'mu',      mu*centi*poise, ...
                           'c',       [0, 0]/barsa); % Default: Incompressible

% Specify schedule information
time = 100*year;
pv = poreVolume(G, rock);
inj_rate = sum(pv)/time;
nt = 300;

% Flux injection and hydrostatic pressure b.c. for pressure support
bc = [];
bc = fluxside(bc, G, 'xmin', inj_rate, 'sat', [0, 1]);
bc = pside(bc, G, 'xmax', 100*barsa, 'sat', [1, 0]);
isP = strcmpi(bc.type, 'pressure');
bc.value(isP) = bc.value(isP) + fluid.rhoWS.*G.faces.centroids(bc.face(isP), 3)*norm(gravity());
schedule = simpleSchedule(rampupTimesteps(time, time/nt), 'bc', bc);

% Setup model and initial state
model = TwoPhaseWaterGasModel(G, rock, fluid, 1, 1, 'useCNVConvergence', true);
model.extraStateOutput = true;

state0 = initResSol(G, 0, [1, 0]);

% Setup the NonLinearSolver
nls = NonLinearSolver('useRelaxation', true);

% Pack the fine scale simulation problem
p_fine = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls, 'Name', 'fine', 'Description', 'Fine-scale');

%% 2) Hybrid VE-fine
% This model consists of a rectangular model with sections of VE grid
% interspersed with sections of fine-scale grid. We create this by
% calling convertToMultiVEModel(model, fine) where fine contains the
% indices of the cells which should stay as fine cells and not be converted
% to VE cells.

% Find the desired indices of the fine scale cells.
G = model.G;
[ii, jj, kk] = gridLogicalIndices(G);

start = ii < 0.1*G.cartDims(1);
mid = ii > 0.4*G.cartDims(1) & ii < 0.6*G.cartDims(1);
stop = ii > 0.9*G.cartDims(1);
fine = find(start | mid | stop);

% Convert model to hybrid VE and specify fine cells.
model_hyb = convertToMultiVEModel(model, fine);

% Upscale initial state to corresponding coarse hybrid grid
state0_hyb = upscaleState(model_hyb, model, state0);

% Upscale schedule
schedule_hyb = upscaleSchedule(model_hyb, schedule);

% Pack the hybrid problem
p_hyb = packSimulationProblem(state0_hyb, model_hyb, schedule_hyb, name, 'NonLinearSolver', nls, 'Name', 'hybrid', 'Description', 'Hybrid');

% Plot hybrid VE grid 
clf;
plotGrid(model_hyb.G)
view([0 0])

%% 3) VE everywhere
% This contains only VE cells (no fine scale sections). We create this by
% calling convertToMultiVEModel(model) without specifying any fine scale
% cells.

% Convert model, state and schedule to VE. We also save the coarse model
% model_c used to generate the VE model for comparison in the next section.
[model_ves, model_c] = convertToMultiVEModel(model);
state0_ves = upscaleState(model_ves, model, state0);
schedule_ves = upscaleSchedule(model_ves, schedule);

% Pack the VE problem
p_ves = packSimulationProblem(state0_ves, model_ves, schedule_ves, name, 'Name', 've', 'Description', 'VE');

%% 4) The coarse model
% This is the coarse model used to generate the fully VE simulation above. 
% In this case the grid for the coarse model looks identical to the grid 
% for the fine model as there are no internal baffles or fine scale grid 
% sections. However, the results of the simulation will be different as it
% is a regular simulation as opposed to a VE simulation.

% Add cart grid dimensions since we know the original grid
model_c.G.cartDims = [G.cartDims(1), 1, 1];
model_c.G.cells.indexMap = (1:model_c.G.cells.num)';

% Pack the coarse problem
p_c = packSimulationProblem(state0_ves, model_c, schedule_ves, name, 'Name', 'coarse', 'Description', 'Coarse-scale');

%% Create combined problem structure
% Here we put the simulation problems together in a cell array to make it
% easier to solve them consecutively.
problems = {p_fine, p_hyb, p_ves, p_c};

% Get descriptions from problems to use for plot titles later.
np = numel(problems);
descr = cellfun(@(x) x.Description, problems, 'UniformOutput', false);

%% Plot the simulation grids for comparison
close all
for i = 1:3
    g = problems{i}.SimulatorSetup.model.G;
    subplot(3, 1, i)
    plotGrid(g, 'edgealpha', 1, 'facec', 'w');
    
    bc = problems{i}.SimulatorSetup.schedule.control.bc;
    inflow = strcmpi(bc.type, 'flux');
    outflow = ~inflow;
    plotFaces(g, bc.face(inflow), 'edgec', 'r', 'linewidth', 3)
    plotFaces(g, bc.face(outflow), 'edgec', 'b', 'linewidth', 3)
    
    view(0, 0)
    axis equal tight
    daspect([1, 1, 0.005])
    set(gca, 'Fontsize', 16);
    title(descr{i})
end

%% Simulate
% simulatePackedProblem will run four simulations, one for each entry in
% problems. 

% The fine scale simulation may take some time to run. The VE and coarse 
% scale simulations should run much faster.

% Using the packed problem setup means we can save simulation results
% whilst they are running and reload them later without having to re-run 
% the simulations.
[ok, status] = simulatePackedProblem(problems);

%% Get the output from the packed problems
[allws, allstates, allreports, names, T] = getMultiplePackedSimulatorOutputs(problems);

%% Plot all results interactively
% For VE solvers, we get two results: The coarse scale plots, and the
% reconstructed fine-scale.
% Using plotToolbar allows us to interactively inspect the results in each
% figure window.
for i = 1:np
    m = problems{i}.SimulatorSetup.model;
    name = descr{i};
    states = allstates{i};
    g = m.G;
    if isempty(states)
        continue
    end
    figure;
    plotToolbar(g, states, 'field', 's:2');
    title(name);
    view(0, 0)
    if isa(m, 'WaterGasMultiVEModel')
        states = convertMultiVEStates(m, states);
        g = m.G.parent;
        figure;
        plotToolbar(g, states, 'field', 's:2');
        title(name);
        view(0, 0)
        title([name, ' (reconstructed)']);
    end
end
%% Plot the plume height at different times for different simulations
close all
figure('position', [680   682   851   296]);
hold on

plottimes = [0.01, 0.05, 0.3, 0.66, 1];
steps = zeros(size(plottimes));
for i = 1:numel(steps)
    steps(i) = find(T{1} >= plottimes(i)*T{1}(end), 1);
end
plotheights = [1, 2, 5, 6, 8];
nstep = numel(steps);


styles = {'-', '--', '-.', '.'};
colors = lines(nstep);
x = problems{3}.SimulatorSetup.model.G.cells.centroids(:, 1);
h = zeros(np, 1);
for pNo = 1:np
    for stepNo = 1:nstep
        step = steps(stepNo);
        g = problems{pNo}.SimulatorSetup.model.G;
        if isfield(g.cells, 'columns')
            col = g.cells.columns;
        else
            [col, jj, kk] = gridLogicalIndices(g);
        end
        counts = accumarray(col, 1);
        sum_sat = accumarray(col, allstates{pNo}{step}.s(:, 2), [max(col), 1]);
        sum_sat = sum_sat./counts;
        height = 10*sum_sat;
        hh = plot(x, height, styles{pNo}, 'color', colors(stepNo, :));
        if stepNo == 1
            h(pNo) = hh;
        end
        
        if pNo == np
            [h_i, sub] = min(abs(height - plotheights(stepNo)));
        end
    end
    set(gca, 'FontSize', 14)
end
legend(h, descr)
xlabel('x [m]');
ylabel('Plume height [m]');



%% Plot saturation evolution through time
for step = 1:5:308
    figure(1); clf; hold on
    for i = 1:np
        g = problems{i}.SimulatorSetup.model.G;
        if isfield(g.cells, 'columns')
            col = g.cells.columns;
        else
            [col, jj, kk] = gridLogicalIndices(g);
        end
        counts = accumarray(col, 1);
        sum_sat = accumarray(col, allstates{i}{step}.s(:, 2), [max(col), 1]);
        sum_sat = sum_sat./counts;
        
        plot(sum_sat)
    end
    legend(descr)
    ylim([0, 1])
    drawnow
end

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
