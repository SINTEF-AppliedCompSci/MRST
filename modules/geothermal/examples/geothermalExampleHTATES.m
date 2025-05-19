%% High-Temperature Aquifer Thermal Energy Storage (HT-ATES) Example
% This example simulates high-temperature aquifer thermal energy storage (HT-ATES):
% - Stores excess thermal energy as hot water in an underground aquifer
% - Demonstrates both storage and extraction phases
% - Uses a precomputed grid from the upr module
% - Shows setup, simulation, and visualization steps
%
% Requirements:
%   - MRST modules: geothermal, upr, ad-core, ad-props, ad-blackoil, mrst-gui, compositional
%   - Run 'startup' in MRST root before executing this script
%
% Steps:
%   1. Add required modules
%   2. Set up plotting options
%   3. Load precomputed grid and rock data
%   4. Visualize geological model
%   5. Define fluids, wells, and schedule
%   6. Run simulation and visualize results
%
% Note: This test case can also be set up using the function 'htates_geothermal()'.

% --- MODULES AND VERBOSITY ---
mrstModule add geothermal upr ad-core ad-props ad-blackoil mrst-gui compositional % Load required MRST modules
mrstVerbose on % Enable verbose output for debugging
savepng = @(name) disp(name); % Dummy function for saving PNGs (replace with actual save if needed)
saveeps = @(name) disp(name); % Dummy function for saving EPS (replace with actual save if needed)

%% Set plot options
% Set up figure handles and color maps for 2D/3D plotting
fig2D = @() figure('Position', [0,0,800,350]); % 2D figure window
fig3D = @() figure('Position', [0,0,1300, 650]); % 3D figure window
alpha = 0.6; % Transparency for colormaps
cmapT = hot*alpha + (1-alpha); % Colormap for temperature
cmap  = jet*alpha + (1-alpha); % Colormap for other properties
% Axes properties for consistent 3D/2D plots
setAxProps3D = @(ax) set(ax, ...
                              'View'              , [-10,45]          , ...
                              'PlotBoxAspectRatio', [4.96, 2.10, 1.00], ...
                              'Projection'        , 'Perspective'     , ...
                              'Box'               , 'on'              , ...
                              'XLimSpec'          , 'tight'           , ...
                              'YLimSpec'          , 'tight'           , ...
                              'ZLimSpec'          , 'tight'           );
setAxProps2D = @(ax) set(ax , ...
                               'Box'     , 'on', ...
                               'FontSize', 12  );

%% Generate grid
% Load a precomputed grid, rock, and well data from the upr module
% This grid represents the aquifer domain for HT-ATES
% 'layer' is used for later filtering/visualization
% 'W' contains well definitions
% 'rock' contains permeability/porosity/thermal properties
% 'G' is the grid structure
% The .mat file must be present in the geothermal dataset path

data  = load(fullfile(getDatasetPath('geothermal'), 'htates-reservoir.mat'));
G     = data.G;
rock  = data.rock;
W     = data.W;
layer = data.layer;

%% Plot the model
% Visualize the geological model (permeability field)
fig3D(), plotCellData(G, log10(rock.perm), 'edgealpha', 0.2); % Log-permeability
setAxProps3D(gca); camlight; mrstColorbar(gca, log10(rock.perm)); colormap(cmap)
savepng('htates-setup');

%% Create fluid and rock
% Set up fluid and rock properties for the simulation
% - Fluid: single-phase water, p/T-dependent properties using EOS
% - Rock: add thermal properties
fluid = initSimpleADIFluid('phases', 'W', 'n', 1, 'mu', 1, 'rho', 1); % Basic water fluid
fluid = addThermalFluidProps(fluid, 'useEOS', true); % Add thermal EOS
rock  = addThermalRockProps(rock); % Add thermal properties to rock

%% Assign well schedule
% Define well controls for storage, rest, and extraction periods
% - Hot wells inject hot water (storage), then produce (extraction)
% - Cold wells provide pressure support
% - Well properties (type, value, temperature) are set based on location
% - Schedule alternates between storage, rest, extraction, rest
% - Each cycle is one year (4+2+4+2 months)

p0    = 70*barsa;         % Reference reservoir pressure
K0    = 273.15*Kelvin;    % Zero degrees Celsius
rate  = 1000*meter^3/day; % Injection rate for hot wells
bhpSt = p0;               % Storage BHP for cold wells
bhpEx = 85*barsa;         % Extraction BHP for cold wells
Thot  = K0 + 100*Kelvin;  % Hot water temperature
Tcold = K0 + 30*Kelvin;   % Cold water temperature

[W.compi]  = deal(1);                % Single-phase composition
W          = addThermalWellProps(W, G, rock, fluid); % Add thermal props to wells
WSt        = W; % Storage wells
[hNo, cNo] = deal(0); % Counters for hot/cold wells
xmid       = mean(G.cells.centroids); % Reservoir center for grouping
for wNo = 1:numel(W)
    w     = W(wNo);
    cells = w.cells;
    xw    = mean(G.cells.centroids(cells,:));
    w.refDepth = min(xw(:,3));
    if xw(1) < xmid(1)
        % Hot group: set injection rate, temperature
        hNo = hNo + 1;
        w.type = 'rate';
        w.val  = rate;
        w.sign = 1;
        w.name = ['H', num2str(hNo)];
        w.T    = Thot;
    else
        % Cold group: set BHP, temperature
        cNo    = cNo + 1;
        w.type = 'bhp';
        w.val  = bhpSt;
        w.sign = -1;
        w.name = ['C', num2str(cNo)];
        w.T    = Tcold;
    end
    WSt(wNo) = w;
end
[WSt.components] = deal(1); % Single-phase
% Make well structures for extraction and rest periods
[WEx, WRe] = deal(WSt);
isHot = cellfun(@(n) strcmpi(n(1), 'H'), {WSt.name}); % Hot group
[WEx(isHot).val ] = deal(-rate); % Reverse rate for extraction
[WEx(~isHot).val] = deal(bhpEx); % Set extraction BHP for cold
for wNo = 1:numel(WEx)
    WEx(wNo).sign = -WEx(wNo).sign; % Reverse sign for extraction
end
[WRe.status] = deal(false); % Rest: all wells inactive
% Define time periods for each phase
month  = year/12;
[dtSt, dtEx] = deal(rampupTimesteps(4*month, 20*day, 5)); % Storage/extraction
[nSt , nEx]  = deal(numel(dtSt));
dtRe = rampupTimesteps(2*month, 20*day, 2); % Rest
nRe  = numel(dtRe);
control   = struct('W', {WSt, WRe, WEx}); % Well controls
controlNo = [1*ones(nSt,1); 2*ones(nRe,1); 3*ones(nEx,1); 2*ones(nRe,1)];
dt        = [dtSt; dtRe; dtEx; dtRe];
nCycles   = 4; % Number of annual cycles
controlNo = repmat(controlNo, nCycles, 1);
dt        = repmat(dt, nCycles, 1);
schedule  = struct('control', control, ...
                   'step'   , struct('val', dt, 'control', controlNo));
pw = @() plotWell(G, WSt, 'color', 'k', 'height', 45); % Plot wells
               
%% Create model
% Set up the single-phase geothermal model
% - Specify valid pressure/temperature ranges for EOS
% - Use MEX-accelerated AD backend and AMGCL CPR solver for speed
% - Enable gravity

gravity reset on
model = GeothermalModel(G, rock, fluid); % Create model
model.maximumPressure    = 200e6;    % Max pressure for EOS
model.minimumTemperature = K0;       % Min temperature for EOS
model.maximumTemperature = K0 + 275; % Max temperature for EOS
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true); % Fast AD
mrstModule add linearsolvers, lsolver = AMGCL_CPRSolverAD(); % Fast linear solver

%% Initialize the reservoir state
% Set up initial hydrostatic and thermal equilibrium
% - Use a geothermal gradient for initial temperature
% - Apply pressure BC at top faces and simulate to equilibrium
state0   = initResSol(G, p0, 1); % Initial pressure, single-phase

% Geothermal gradient: 30 K/km
% T0 is a function for initial temperature by depth
% Assign initial temperature to all cells

dTdz     = 30*Kelvin/(kilo*meter); % Geothermal gradient
T0       = @(x) K0 + 20*Kelvin + dTdz*x(:,3);
state0.T = T0(G.cells.centroids);
% To find equilibrium, apply fixed pressure BC to top faces and simulate
f       = boundaryFaces(G);
[~, ix] = sort((G.faces.centroids(f,3)));
nf      = 1; % Number of top faces for BC
bc      = addBC([], f(ix(1:nf)), 'pressure', p0, 'sat', 1);
bc      = addThermalBCProps(bc, 'T', T0(G.faces.centroids(f(ix(1:nf)),:)));
bc.x    = ones(numel(bc.face),1);
dtPre             = rampupTimesteps(50*year, 10*year, 3); % Long timesteps
schedulePre       = simpleSchedule(dtPre, 'bc', bc);
[~, statesPre, ~] = simulateScheduleAD(state0, model, schedulePre, ...
                                                  'linearSolver', lsolver);
state0 = statesPre{end}; % Set initial state to equilibrium
state0 = rmfield(state0, 'wellSol'); % Remove well solution field

%% Plot initial state
% Visualize initial temperature/pressure distribution
fig3D(), plotToolbar(G, state0, 'edgealpha', 0.2); pw();
setAxProps3D(gca); view([-200,45]), camlight, view([-10,45]); colormap(cmap);

%% Add plot hook
% Set up a plotting function to visualize temperature after each timestep
% This helps monitor simulation progress interactively
fn = getPlotAfterStep(state0, model, schedule, ...
                      'useTimesteps', true, 'field', 'T' , ...
                      'plotWellSolArgs', {'field', 'T', 'SelectedWells', 1:numel(W)});
pw(); % Plot wells
plotGrid(G, 'faceColor', 'none', 'edgealpha', 0.2), setAxProps3D(gca); colormap(cmapT);
camlight();

%% Simulate HT-ATES
% Run the full simulation using the packed problem interface
% - Uses the schedule with storage/rest/extraction cycles
% - Uses the plot hook for interactive visualization
problem = packSimulationProblem(state0, model, schedule, 'ht-ates', ...
    'ExtraArguments', {'linearSolver', lsolver, 'afterStepFn', fn});
simulatePackedProblem(problem, 'restartStep', 1);

%% Interactive plot of the results
% Visualize the final results interactively (temperature, pressure, etc.)
[wellSols, states, reports] = getPackedSimulatorOutput(problem);
fig3D(), plotToolbar(G, states, 'edgealpha', 0.2); pw();
setAxProps3D(gca); view([-200,45]), camlight, view([-10,45]); colormap(cmap);
dtime = schedule.step.val;
plotWellSols(wellSols, dtime);

%% Animate thermal plume
% Animate the development of the thermal plume in the reservoir
% - Plots cells where temperature deviates >20% from initial
% - Only shows cells below a certain layer (removes top)
fig3D(); alpha = 0.2;
for t = 1:numel(states)
    clf;
    plotGrid(G, 'edgealpha', 0.2, 'faceColor', 'none')
    setAxProps3D(gca); view([-200,45]), camlight, view([-10,45]); colormap(cmapT);
    c = (states{t}.T-K0 > (state0.T-K0)*(1+alpha)  ... % Warm region
       | states{t}.T-K0 < (state0.T-K0)*(1-alpha)) ... % Cold region
       & layer > 5;                                    % Remove top
    plotCellData(G, states{t}.T, c, 'edgecolor', 'none'); pw(); axis off
    caxis([Tcold, Thot]); colorbar('location', 'south');
    pause(0.1)
end

%% Analyse HT-ATES performance
% Post-process simulation results to analyze energy storage and recovery
% - Compute cumulative injected/extracted energy
% - Plot energy recovery factor and well temperature over time
% - Visualize storage, rest, and extraction periods

time   = cumsum(dtime)/year;
period = sum(bsxfun(@gt, time, year*(1:nCycles)),2) + 1;
T  = getWellOutput(wellSols, 'T');      % Well temperature
p  = getWellOutput(wellSols, 'bhp');    % Bottom-hole pressure
qW = getWellOutput(wellSols, 'qWs');    % Water rate
st = getWellOutput(wellSols, 'status'); % Well status
[volume, energy] = deal(zeros(numel(dtime), numel(W)));
for i = 1:numel(W)
    volume(:,i) = st(:,i).*qW(:,i).*dtime;
    energy(:,i) = volume(:,i).*fluid.rhoW(p(:,i), T(:,i)).*fluid.uW(p(:,i), T(:,i));
end
energySt = cumsum(sum(energy(:,isHot),2).*(schedule.step.control == 1));
energyEx = cumsum(sum(energy(:,isHot),2).*(schedule.step.control == 3));
% Plot cumulative energy recovery factor
fig2D(); setAxProps2D(gca), hold on
ix = find(abs(diff(schedule.step.control)) > 0);
cycle = 2;
clr   = 0.2*ones(4,3);
clr([1,3],:) = flipud(lines(2));
clr = brighten(clr, 0.8);
for i = 1:4
    x = time([ix(4*(cycle-1)+i-1), ix(4*(cycle-1)+i)]);
    y = [0, 1];
    patch(x([1,2,2,1]), y([1,1,2,2]), clr(i,:));
end
for i = 1:nnz(ix)
    line([time(ix(i)) time(ix(i))], [0 1], ...
                         'linestyle', '-', 'color', 'k', 'linewidth', 0.1)
end
plot(time, abs(energyEx./energySt), 'k', 'LineWidth', 2),
xlim([0, time(end)]); hold off
xlabel('Time (years)'), ylabel('Energy recovery factor');
legend({'Storage', 'Rest', 'Extraction'}, 'Location', 'southeast')
savepng('htates-erf');
saveeps('htates-erf');

Tlim = [0, Thot-K0 + 10];
% Plot well temperature as a function of time
fig2D(), setAxProps2D(gca), hold on
cycle = 2;
for i = 1:4
    x = time([ix(4*(cycle-1)+i-1), ix(4*(cycle-1)+i)]);
    y = [0, Thot-K0 + 10];
    patch(x([1,2,2,1]), y([1,1,2,2]), clr(i,:));
end
for i = 1:nnz(ix)
    line([time(ix(i)) time(ix(i))], [0 Thot - K0 + 10], ...
                         'linestyle', '-', 'color', 'k', 'linewidth', 0.1)
end
xlim([0, time(end)]); ylim(Tlim);
lclr = lines(2);
red  = lclr(2,:);
blue = lclr(1,:);
for i = 1:numel(W)/2
    plot(time, T(:,i) - K0, 'LineWidth', 2, 'color', red);
    plot(time, T(:,i+4) - K0, 'LineWidth', 2, 'color', blue);
end
hold off
xlabel('Time (years)'), ylabel('Well temperature (C)');
legend({'Storage', 'Rest', 'Extraction'}, 'Location', 'east');
savepng('htates-wtemp');
saveeps('htates-wtemp');

%% Copyright Notice
%
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
