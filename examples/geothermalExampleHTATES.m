%% High-Temperature Aquifer Thermal Energy Storage
% In this example, we simulate high-temperature aquifer thermal energy
% storage (HT-ATES), where excess thermal energy is stored in the form of
% hot water in an undergroud aquifer in periods of energy excess, and
% extracted for use in e.g., heating systems in periods of demand.
% NOTE: this test case also exists as a test case setup function, and can
% be set up with
%   setup = htates_geothermal();
% For pedagogical purposes, this script nevertheless goes thorugh all the
% steps of setting up the test case.
mrstModule add geothermal upr ad-core ad-props ad-blackoil mrst-gui compositional
mrstVerbose on
savepng = @(name) disp(name); % Dummy function
saveeps = @(name) disp(name); % Dummy function

%% Set plot options
% For pretty plotting
fig2D = @() figure('Position', [0,0,800,350]);
fig3D = @() figure('Position', [0,0,1300, 650]);
alpha = 0.6;
cmapT = hot*alpha + (1-alpha);
cmap  = jet*alpha + (1-alpha);
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
% We use a grid from the upr module. For convenience, this is precomputed
% and stored in the module directory
data  = load(fullfile(getDatasetPath('geothermal'), 'htates-reservoir.mat'));
G     = data.G;
rock  = data.rock;
W     = data.W;
layer = data.layer;

%% Plot the model
fig3D(), plotCellData(G, log10(rock.perm), 'edgealpha', 0.2);
setAxProps3D(gca); camlight; mrstColorbar(gca, log10(rock.perm)); colormap(cmap)
savepng('htates-setup');

%% Create fluid and rock
% Viscosity and density are p/T-dependent, and will be set later
fluid = initSimpleADIFluid('phases', 'W', 'n', 1, 'mu', 1, 'rho', 1);
% Assign thermal properties of the fluid, with equation of state from
% Spivey et. al (2004)
fluid = addThermalFluidProps(fluid, 'useEOS', true);
rock  = addThermalRockProps(rock);            % Thermal rock properties

%% Assign well schedule
% The model has two groups of four wells each, which we refer to as hot and
% cold. The hot group will be used to inject hot water for storage during
% the summer months, and extract the hot water for heating in the winter
% months. The cold well is used to provide pressure support, with a low BHP
% during storage, and a hight BHP during extraction. In between storage and
% extraction, there will be a period of rest.
p0    = 70*barsa;         % Reference reservoir pressure
K0    = 273.15*Kelvin;    % Zero degrees Celcius
rate  = 1000*meter^3/day; % Injection rate
bhpSt = p0;               % Storage BHP
bhpEx = 85*barsa;         % Extraction BHP
Thot  = K0 + 100*Kelvin;  % Temperature of stored water
Tcold = K0 + 30*Kelvin;   % Temperature of injected water

[W.compi]  = deal(1);                % Single-phase
W          = addThermalWellProps(W, G, rock, fluid); % Add thermal properties
WSt        = W;
[hNo, cNo] = deal(0);
xmid       = mean(G.cells.centroids); % Reservoir center
for wNo = 1:numel(W)
    w     = W(wNo);
    cells = w.cells;
    xw    = mean(G.cells.centroids(cells,:));
    w.refDepth = min(xw(:,3));
    if xw(1) < xmid(1)
        % Hot group: set injection rate
        hNo = hNo + 1;
        w.type = 'rate';
        w.val  = rate;
        w.sign = 1;
        w.name = ['H', num2str(hNo)];
        w.T    = Thot;
    else
        % Cold group: set bhp
        cNo    = cNo + 1;
        w.type = 'bhp';
        w.val  = bhpSt;
        w.sign = -1;
        w.name = ['C', num2str(cNo)];
        w.T    = Tcold;
    end
    WSt(wNo) = w;
end
[WSt.components] = deal(1);
% Make well structures for the extraction and rest periods
[WEx, WRe] = deal(WSt);
isHot = cellfun(@(n) strcmpi(n(1), 'H'), {WSt.name}); % Hot group
[WEx(isHot).val ] = deal(-rate); % Reverse rate
[WEx(~isHot).val] = deal(bhpEx); % Set extraction bhp
for wNo = 1:numel(WEx)
    WEx(wNo).sign = -WEx(wNo).sign; % Reverse sign
end
[WRe.status] = deal(false); % Rest wells are not active (status = false)
% Storage and extraction periods are four months long, whereas the rest
% periods are two months
month  = year/12;
% Storage period
[dtSt, dtEx] = deal(rampupTimesteps(4*month, 20*day, 5));
[nSt , nEx]  = deal(numel(dtSt));
% Rest period
dtRe = rampupTimesteps(2*month, 20*day, 2);
nRe  = numel(dtRe);
% Controls: storage = 1, rest = 2, extraction = 3
control   = struct('W', {WSt, WRe, WEx});
controlNo = [1*ones(nSt,1);  % Storage
             2*ones(nRe,1);  % Rest
             3*ones(nEx,1);  % Extraction
             2*ones(nRe,1)]; % Rest
dt        = [dtSt; dtRe; dtEx; dtRe]; % Timesteps
% We simulate eight cycles (= eight years) of HTATES
nCycles   = 4;
controlNo = repmat(controlNo, nCycles, 1);
dt        = repmat(dt, nCycles, 1);
schedule  = struct('control', control, ...
                   'step'   , struct('val', dt, 'control', controlNo));
pw = @() plotWell(G, WSt, 'color', 'k', 'height', 45); % For plotting wells
               
%% Create model
% We use a single-phase geothermal model
gravity reset on
model = GeothermalModel(G, rock, fluid); % Make model
% The EOS is valid for pressure/temperature within a given range. We
% provide these to the model so that pressure/temperature are within these
% during the nonlinear solution step
model.maximumPressure    = 200e6;    % Maximum pressure
model.minimumTemperature = K0;       % Minimum temperature 
model.maximumTemperature = K0 + 275; % Maximum temperature 
% To speed up the simulation, we use mex-accelereated AD backend ...
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
% ... and a compiled iterative Krylov solver with CPR/AMG preconditioning
mrstModule add linearsolvers, lsolver = AMGCL_CPRSolverAD();

%% Initialize the reservoir state
% The reserois should be in thermal and hydrostatic equilibrium
state0   = initResSol(G, p0, 1);
dTdz     = 30*Kelvin/(kilo*meter); % Geothermal gradient
T0       = @(x) K0 + 20*Kelvin + dTdz*x(:,3);
state0.T = T0(G.cells.centroids);
% state0.T = K0 + 20*Kelvin;
% To find the equilibrium state, we simply assign a fixed pressure BC to
% the 100 topmost faces, and simulate for 50 years with looong timesteps
f       = boundaryFaces(G);
[~, ix] = sort((G.faces.centroids(f,3)));
nf      = 1;
bc      = addBC([], f(ix(1:nf)), 'pressure', p0, 'sat', 1);
bc      = addThermalBCProps(bc, 'T', T0(G.faces.centroids(f(ix(1:nf)),:)));
bc.x    = ones(numel(bc.face),1);
% Make schedule
dtPre             = rampupTimesteps(50*year, 10*year, 3);
schedulePre       = simpleSchedule(dtPre, 'bc', bc);
[~, statesPre, ~] = simulateScheduleAD(state0, model, schedulePre, ...
                                                  'linearSolver', lsolver);
state0 = statesPre{end}; % Set initial state
state0 = rmfield(state0, 'wellSol');

%% Plot initial state
fig3D(), plotToolbar(G, state0, 'edgealpha', 0.2); pw();
setAxProps3D(gca); view([-200,45]), camlight, view([-10,45]); colormap(cmap);

%% Add plot hook
% We use afterStepFn to plot the progress during the simulation. Tip: Check
% out the temperature evaluation in the well plot
fn = getPlotAfterStep(state0, model, schedule, ...
                      'useTimesteps', true, 'field', 'T' , ...
                      'plotWellSolArgs', {'field', 'T', 'SelectedWells', 1:numel(W)});
pw(); % Plot wells
plotGrid(G, 'faceColor', 'none', 'edgealpha', 0.2), setAxProps3D(gca); colormap(cmapT);
camlight();

%% Simulate HT-ATES
problem = packSimulationProblem(state0, model, schedule, 'ht-ates', ...
    'ExtraArguments', {'linearSolver', lsolver, 'afterStepFn', fn});
simulatePackedProblem(problem, 'restartStep', 1);

%% Interactive plot of the results
[wellSols, states, reports] = getPackedSimulatorOutput(problem);
fig3D(), plotToolbar(G, states, 'edgealpha', 0.2); pw();
setAxProps3D(gca); view([-200,45]), camlight, view([-10,45]); colormap(cmap);
dtime = schedule.step.val;
plotWellSols(wellSols, dtime);

%% Animate thermal plume
% We visualize development of the reservoir temperature by plotting the
% cells where the temerature deviates more that 20% from the initial
% temperature
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
% Finally, we look at the performance of the HT-ATES system.
time   = cumsum(dtime)/year;
period = sum(bsxfun(@gt, time, year*(1:nCycles)),2) + 1;
% Get well output
T  = getWellOutput(wellSols, 'T');      % Temperature
p  = getWellOutput(wellSols, 'bhp');    % Bottom-hole pressure
qW = getWellOutput(wellSols, 'qWs');    % Water injection/production rate
st = getWellOutput(wellSols, 'status'); % Status (on/off)
% Compute injected/extracted energy
[volume, energy] = deal(zeros(numel(dtime), numel(W)));
for i = 1:numel(W)
    volume(:,i) = st(:,i).*qW(:,i).*dtime;
    energy(:,i) = volume(:,i).*fluid.rhoW(p(:,i), T(:,i)).*fluid.uW(p(:,i), T(:,i));
end
energySt = cumsum(sum(energy(:,isHot),2).*(schedule.step.control == 1));
energyEx = cumsum(sum(energy(:,isHot),2).*(schedule.step.control == 3));
% Plot cumulative energy recovery factor (i.e., cumulative extracted over
% stored energy) as a function of time
fig2D(); setAxProps2D(gca), hold on
ix = find(abs(diff(schedule.step.control)) > 0);
% Add in lines indicating storage/rest/extraction
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
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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