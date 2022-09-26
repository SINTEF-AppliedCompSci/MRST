%% Subset of SEP10 Model 2
% In this example, we simulate injection of hot water into Layer 10 of
% SPE10 Model 2, once with thermally insulated boundary conditions, and
% once with fixed-temperature boundary conditions.

%% Add modules
mrstModule add geothermal upr spe10 mrst-gui ad-props ad-core ad-blackoil compositional

%% Set up fine-scale model
% We extract the lower half of layer 13 of SPE10 2
[state0Ref, modelRef, scheduleRef] = setupSPE10_AD('layers', 10);
GRef    = modelRef.G;
rockRef = modelRef.rock;
WRef    = scheduleRef.control(1).W;

%% Make PEBI grid
% We use the upr module to construct a PEBI grid with refinement around the
% wells
rng(2019) % For reproducibility
n = 40;   % Approx number of cells in each direction
% Get well coordinates
l = max(GRef.nodes.coords(:,1:2));
wellLines = mat2cell(GRef.cells.centroids(vertcat(WRef.cells),1:2), ...
                                                    ones(numel(WRef),1), 2)';
% Construct PEBI grid
G = pebiGrid2D(max(l)/n, l, 'cellConstraints', wellLines, ... % Well coords
                            'CCRefinement'   , true     , ... % Refine
                            'CCFactor'       , 0.4      );
G = computeGeometry(G);       % Compute geometry

%% Sample rock properties
% We assign rock properties in the PEBI grid cells by sampling from the
% fine grid using sampleFromBox
poro = sampleFromBox(G, reshape(rockRef.poro, GRef.cartDims));
perm = zeros(G.cells.num,G.griddim);
for i = 1:G.griddim
    perm(:,i) = sampleFromBox(G, reshape(rockRef.perm(:,i), GRef.cartDims));
end
rock = makeRock(G, perm, poro);
rock = addThermalRockProps(rock);

%% Create model
% We use a single-phase geothermal model. Viscosity and density are
% p/T-dependent, and will be set later
fluid = initSimpleADIFluid('phases', 'W', 'n', 1, 'mu', 1, 'rho', 1);
% Assign thermal properties of the fluid, with equation of state from
% Spivey et. al (2004)
fluid = addThermalFluidProps(fluid, 'useEOS', true);
% Add thermal rock properties, with high thermal conductivity and low
% specific heat capacity
Watt  = joule/second;
rock  = addThermalRockProps(rock, 'lambdaR', 20*Watt/(meter*Kelvin)     , ...
                                  'CpR'    , 100*joule/(kilogram*Kelvin));
model = GeothermalModel(G, rock, fluid); % Make model
% The EOS is valid for pressure/temperature within a given range. We
% provide these to the model so that pressure/temperature are within these
% during the nonlinear solution step
K0 = 273.15*Kelvin;
model.maximumPressure    = 200e6;           % Maximum pressure
model.minimumTemperature = K0;              % Minimum temperature 
model.maximumTemperature = K0 + 275*Kelvin; % Maximum temperature
model.extraStateOutput   = true;     % Output density and mobility to state
model.outputFluxes       = true;     % Output fluxes to state

%% Set up schedule and initial state
W        = WRef;
schedule = scheduleRef;
x  = G.cells.centroids;
xwR = GRef.cells.centroids(vertcat(WRef.cells),1:2);
% Slick oneliner to find corresponding cells in the new grid
[~, c] = min(sum(bsxfun(@minus, reshape(xwR, [], 1 , G.griddim), ...
                             reshape(x , 1 , [], G.griddim)).^2,3), [], 2);
for i = 1:numel(W)
    W(i).cells = c(i);
    W(i).components = 1;
    if strcmpi(W(i).type, 'rate')
        W(i).val = 2*W(i).val; % Increase injection rate
    end
end
Tinj = (273.15 + 100)*Kelvin;             % 100 degrees Celsius
W    = addThermalWellProps(W, G, rock, fluid, 'T', Tinj); % Add temperature field T to wells
xw = G.cells.centroids(vertcat(W.cells),1:2);
% [W.components] = deal([0,1]);
schedule.control.W = W;
% Set initial state
state0   = initResSol(G, state0Ref.pressure(1), 1);
Tres     = (273.15 + 30)*Kelvin;
state0.T = repmat(Tres, G.cells.num, 1);

%% Inspect the geological model
figure('Position', [0,0,500,800])
K = convertTo(rock.perm(:,1),milli*darcy);
plotCellData(G,log10(K),'edgeAlpha',.1); axis tight;
set(gca,'FontSize',12)
mrstColorbar(K,'South',true); axis equal tight
hold on; plot(xw(:,1),xw(:,2),'.r','MarkerSize',18); hold off; drawnow

%% Simulate
% No boundary conditions means a completely insulated, closed flow
% compartment
[wellSols, states, reports] = simulateScheduleAD(state0, model, schedule);

%% Interactive plot of the results
figure('Position', [0,0,500,800], 'Name', 'Insulated BCs')
plotToolbar(model.G, states, 'field', 'T'); colormap(hot); axis equal tight
set(gca,'FontSize',12);
plotWellSols(wellSols, schedule.step.val, 'field', 'T');

%% Simulate with fixed-temperature boundary conditions
bc = addBC([], boundaryFaces(G), 'flux', 0);
bc = addThermalBCProps(bc, 'T', Tres);
scheduleBC = schedule;
scheduleBC.control(1).bc = bc;
[wellSolsBC, statesBC, reportsBC] = simulateScheduleAD(state0, model, scheduleBC);

%% Compare with insulated reservoir results
figure('Position', [0,0,500,800], 'Name', 'Fixed-temperature BCs')
dT = cellfun(@(st, stBC) st.T - stBC.T, states, statesBC, 'UniformOutput', false);
t = linspace(0,1,100)'; cmap = vertcat([0,0,1].*(1-t) + t, (1-t) + [1,0,0].*t);
plotToolbar(model.G, dT, 'field', 'T'); colormap(cmap); axis equal tight
set(gca,'FontSize',12)
plotWellSols({wellSols, wellSolsBC}, schedule.step.val, ...
    'dataSetNames', {'Insulated', 'Fixed-temp'}, 'field', 'T', 'SelectedWells', 1:4);

%% Compare advective and conductive heat fluxes
% Compute cell velocities from face fluxes through mimetic mapping
vAdv  = cellfun(@(st) sqrt(sum(faceFlux2cellVelocity(G, st.heatFluxAdv ).^2,2)), states, 'UniformOutput', false);
vCond = cellfun(@(st) sqrt(sum(faceFlux2cellVelocity(G, st.heatFluxCond).^2,2)), states, 'UniformOutput', false);
ratio = cellfun(@(vAdv, vCond) vAdv./vCond, vAdv, vCond, 'UniformOutput', false);
[~, remove] = boundaryFaces(G); remove = vertcat(remove, W.cells);
cells = true(G.cells.num,1); cells(remove) = false;

% Plot
figure('Position', [0,0,500,800], 'Name', 'Convective vs advective heat flux')
plotToolbar(G, ratio, cells, 'log10', true); colormap(jet), axis equal tight

vAdvBC  = cellfun(@(st) sqrt(sum(faceFlux2cellVelocity(G, st.heatFluxAdv ).^2,2)), statesBC, 'UniformOutput', false);
vCondBC = cellfun(@(st) sqrt(sum(faceFlux2cellVelocity(G, st.heatFluxCond).^2,2)), statesBC, 'UniformOutput', false);
ratioBC = cellfun(@(vAdv, vCond) vAdv./vCond, vAdvBC, vCondBC, 'UniformOutput', false);
% Plot
figure('Position', [0,0,500,800], 'Name', 'Convective vs advective heat flux')
plotToolbar(G, ratioBC, cells, 'log10', true); colormap(jet), axis equal tight

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