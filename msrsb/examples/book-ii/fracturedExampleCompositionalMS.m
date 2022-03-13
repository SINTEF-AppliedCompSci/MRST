%% Compositional fractured example
% This example is similar to the example from Section 5.2 of Moyner &
% Tchelepi, SPE J, 23(6), 2018, except that we use a simpler 3-component
% fluid model, somewhat higher fracture permeability, and slightly
% different well positions. See fracturedExampleMS.m for more details on
% the setup.

mrstModule add ad-core ad-blackoil deckformat msrsb
mrstModule add ad-props mrst-gui sequential
mrstModule add coarsegrid
mrstModule add compositional
mrstModule add linearsolvers

%% Notice on Computational Cost
warning('ComputationalCost:High', ...
       ['Please be advised that this example often takes a long time ', ...
        'to run (e.g., more than 1 hour of CPU time)']);
pause(10)

%% Build grid
% Name tags for fluid model and case
fluid_name = 'simple';
caseName = ['FractureMS', '_', fluid_name];

% Load the grid and set up petrophysical model
pth = fullfile(getDatasetPath('MSFractures'), 'setup_fracture.mat');
load(pth);
rdim  = [1000 500];
G.nodes.coords = bsxfun(@times, G.nodes.coords, rdim);
G    = computeGeometry(G);
rock = makeRock(G, perm*milli*darcy, 0.3);
rock.perm(G.cells.tag > 0) = 10*darcy;
pv   = poreVolume(G, rock);

% Define fluid and PVT model
[cf, info] = getBenchmarkMixture(fluid_name);
eos   =  EquationOfStateModel(G, cf);
fluid = initSimpleADIFluid('rho', [1000, 500, 500], ...
                       'mu', [1, 1, 1]*centi*poise, ...
                       'n', [2, 2, 2], ...
                       'c', [1e-5, 0, 0]/barsa);

% Setup the wells
minP    = 50*barsa;
resP    = info.pressure;
totTime = 7*year;
irate   = 0.25*sum(pv)/totTime;
icell   = findEnclosingCell(G, [0.025, 0.05].*rdim);
pcell   = findEnclosingCell(G, [0.975, 0.95].*rdim);
W = addWell([], G, rock, icell, 'comp_i', [0, 1], ...
            'name', 'I', 'Type', 'rate', 'Val', irate);
W = addWell(W, G, rock, pcell, 'comp_i', [0.5, 0.5], ...
            'Name', 'P', 'Val', minP, 'sign', -1);
for i = 1:numel(W)
    W(i).components = info.injection;
end

% Build the model and set accelerated AD backend
model = OverallCompositionCompositionalModel(G, rock, fluid, eos, 'water', false);
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);

% Set the initial state, which must be expanded with information about
% components and temperature
state0   = initResSol(G, resP, [1 0]);
state0.T = repmat(info.temp, G.cells.num, 1);    
state0.components = repmat(info.initial, G.cells.num, 1);

% Build the simulation schedule: we run with uniform time steps of 20 days,
% which give a reasonable compromise between having too high CFL numbers in
% the fractures and too low CFL numbers in the background matrix. In
% addition, we add a standard rampup to stabilize the displacement fronts
% as they move into the reservoir.
dt       = rampupTimesteps(totTime, 20*day);
schedule = simpleSchedule(dt, 'W', W);

% Build packing utility for storing simulation setup
packer = @(model, name, varargin) ...
    packSimulationProblem(state0, model, schedule, caseName, 'name', name, varargin{:});

%% Define fully-implicit solver
nls = NonLinearSolver();
nls.LinearSolver = selectLinearSolverAD(model);
nls.useRelaxation = true;
if isa(nls.LinearSolver, 'AMGCL_CPRSolverAD')
    nls.LinearSolver.strategy = 'amgcl';
end
base = packer(model, 'Fully-implicit', 'NonLinearSolver', nls);

%% Setup baseline sequential simulator
% The baseline simulator uses external linear solvers from the AMGCL header
% library to solve the pressure equation and the transport equation. The
% pressure solver uses CPR, whereas for the transport equation, we must
% reorder the unknowns from unknown-first to cell-first ordering.
ord    = getCellMajorReordering(G.cells.num, size(state0.components, 2));
psolve = AMGCLSolverAD('preconditioner', 'amg', 'tolerance', 1e-5);
tsolve = AMGCLSolverAD('preconditioner', 'relaxation', 'relaxation', 'ilu0', ...
                       'tolerance', 1e-4, 'maxIterations', 25);
tsolve.equationOrdering = ord;
tsolve.variableOrdering = ord;

seqmodel = getSequentialModelFromFI(model);
seqmodel.transportModel.AutoDiffBackend         = model.AutoDiffBackend;
seqmodel.transportNonLinearSolver.useRelaxation = true;
seqmodel.transportNonLinearSolver.LinearSolver  = tsolve;
seqmodel.pressureNonLinearSolver.LinearSolver   = psolve;

seq = packer(seqmodel, 'Sequential');

%% Build coarse grid and multiscale solver
% Based on our observations from the two-phase version of the problem, we
% use a simple 10x10 partition instead of a more unstructured partition
% that adapts to the fractures. For the first multiscale solver, we set the
% convergence criterion to be a factor 50 times looser than the default of
% 1e-3. The linear solver is two orders of magnitude more relaxed than the
% strict AMG settings used above.

pr = sampleFromBox(G, reshape(1:100,[10 10]));
pr = processPartition(G, compressPartition(pr));
CG = generateCoarseGrid(G, compressPartition(pr));
CG = coarsenGeometry(CG);
CG = storeInteractionRegion(CG);
CG = setupMexInteractionMapping(CG);

msmodel = getSequentialModelFromFI(model);
msmodel.transportModel.AutoDiffBackend = model.AutoDiffBackend;
msmodel.pressureModel.incTolPressure = 0.05;
msmodel = addMultiscaleSolverComp(msmodel, CG, 'maxIterations', 50, ...
                                               'useGMRES', true, ...
                                               'tolerance', 1e-3);
msmodel.transportNonLinearSolver.useRelaxation = true;
msmodel.transportNonLinearSolver.LinearSolver  = tsolve;

ms = packer(msmodel, 'Multiscale (Relaxed)');

%% Solve a stricter multiscale model with fine-scale tolerances
strictmsmodel = getSequentialModelFromFI(model);
strictmsmodel.transportNonLinearSolver.useRelaxation = true;
strictmsmodel.transportNonLinearSolver.LinearSolver = tsolve;

strictmsmodel.transportModel.AutoDiffBackend = model.AutoDiffBackend;
strictmsmodel.pressureModel.incTolPressure = 1e-3;

strictmsmodel = addMultiscaleSolverComp(strictmsmodel, CG, ...
                                          'maxIterations', 50, ...
                                          'useGMRES', true, ...
                                          'tolerance', 1e-5);
strictms = packer(strictmsmodel, 'Multiscale (Strict)');

%% Simulate the problems and get the output
problems = {base, seq, ms, strictms};
simulatePackedProblem(problems, 'continueOnError', false);%, 'restartStep', 1);
[ws, states, reports, names, T] = ...
    getMultiplePackedSimulatorOutputs(problems, 'readFromDisk', false, ...
    'readWellSolsFromDisk', true, 'readReportsFromDisk', true);

%% Plot well solutions
plotWellSols(ws, T, 'datasetnames', names, 'linestyles', {'-'}, ...
    'SelectedWells', 2, 'field', 'qGs');

%% Plot the iterations taken in different parts of the solver
ns = numel(reports);
stats = cell(1, ns);
for i = 1:ns
    stats{i} = getPressureTransportIterations(reports{i});
end

total = zeros(ns, 3);
for i = 1:ns
    s = stats{i};
    total(i, 1) = sum(s.pressure);
    total(i, 2) = sum(s.transport);
    total(i, 3) = sum(s.outer);
end
figure;
bar(total)
set(gca, 'XTickLabel', names, 'XTickLabelRotation',20);
legend('Pressure', 'Transport', 'Outer');

%% Prepare for animation of solutions: plot initial state
% We plot hydrocarbon saturation and component concentration
G.cells.sortedCellNodes = getSortedCellNodes(G);

afig = figure('Position',[580 200 860 550]);
flds = {@(state) state.s(:, 2),  @(state) state.components};
nd = numel(flds);
hp = zeros(ns, nd);
for data = 1:nd
    getData = flds{data};
    for i = 1:ns
        subplot(nd,ns, i + ns*(data-1))
        hp(i, data) = plotCellData(G, getData(state0),'EdgeColor','none');
        view(-90,90), axis tight, caxis([0 1]), colormap(flipud(winter(10).^1.5))
        plotGrid(G, vertcat(schedule.control(1).W.cells),'FaceColor','r','EdgeColor','r');
        if data == 1
            caxis([0, 1]);
            title(names{i},'FontSize',10,'FontWeight','normal');
        end
        set(gca, 'XTickLabel', [], 'YTickLabel', [], 'FontSize',12);
    end
end

% Add a colorbar and a ternary diagram outside the respective axes
cb = colorbar('position',[.93 .58 .014 .32]);
tp = axes('position',[.915 .13 .07 .1]);
patch('Vertices', [0 0; 1 0; 0.5 0.5], 'Faces', 1:3, ...
    'FaceVertexCData',[1 0 0; 0 1 0; 0 0 1],'FaceColor','interp');
text(.1,.05,'CH4','Color','w','FontSize',6);
text(.9,.05,'CO2','Color','w','FontSize',6,'HorizontalAlignment','right');
axis off

%% Run the animation
figure(afig);
for stepNo = 1:5:numel(schedule.step.val)
    for i = 1:ns
        s = states{i}{stepNo};
        for data = 1:nd
            getData = flds{data};
            v = getData(s);
            subplot(nd,ns, i + ns*(data-1))
            if size(v, 2) == 1
                fn = 'CData';
            else
                fn = 'FaceVertexCData';
            end
            set(hp(i, data), fn, v);
            if data == 1
                caxis([0, 1]);
            end
        end
    end
    drawnow
end

% Run the following code to get the plots closer to each other
%{
for i=1:ns 
    subplot(2,4,i),    ax=gca;
    ax.OuterPosition = ax.OuterPosition+[0 -.04 0 .04];
    subplot(2,4,i+ns), ax=gca; 
    ax.OuterPosition = ax.OuterPosition+[0 0 0 .04];
end
%}

%% Compute and plot the volume error for all solvers
% Exact mass conservation in transport often implies a volume discrepancy.
% We plot this error to verify that it is within acceptable ranges.
pv = model.operators.pv;
s_err = nan(numel(schedule.step.val), numel(states));
for step = 1:numel(schedule.step.val)
    for i = 1:numel(states)
        state = states{i}{step};
        e = abs(1 - sum(state.s, 2));
        s_err(step, i) = sum(e.*pv)/sum(pv);
    end
end
figure;
plot(T{1}/day, s_err,'LineWidth',2)
set(gca, 'YScale', 'log')
xlabel('Time [days]')
axis tight
legend(names), title('Volume errors');

%% Plot component production
% We plot the component production for each of the four simulations. To be
% able to include two legends, one for the simulation method and one for
% the component, we need to do a trick and plot the first component twice
% for each simulation and add an extra invisible axis.
styles = {'-', '--', '.', ':'};
colors = lines(3);
figure; hold on
start = 20;
h = zeros(3,4);
for c = 1:3
    ind = [1 1:numel(names)];
    for i = 1:numel(ind)
        cprod = -cellfun(@(x) x(2).components(c),  ws{ind(i)});
        h(c,i)=plot(T{ind(i)}(start:end)/day, cprod(start:end)*day, ...
            styles{ind(i)}, 'color', colors(c, :),'LineWidth',2);
    end
end
hold off
axis tight, xlabel('Time [days]'), ylabel('kg/day')
ax1 = gca;
ax2 = axes('position',get(gca,'position'),'visible','off');
legend(ax1, h(1,2:end), names,'Location','East')
legend(ax2, h(:,1), model.EOSModel.getComponentNames,'Location','SouthEast');
set([ax1 ax2],'FontSize',12); title('Component production');

%% Estimate and plot the maximum CFL numbers of the whole simulation
m = GenericOverallCompositionModel(G, rock, fluid, eos, 'water', false);
m = m.validateModel();
cfl_s = zeros(G.cells.num, 1);
cfl_z = zeros(G.cells.num, m.getNumberOfComponents());
for i = 1:numel(schedule.step.val)
    s = states{1}{i};
    dt = schedule.step.val(i);
    cs = estimateSaturationCFL(m, s, dt);
    cz = estimateCompositionCFL(m, s, dt);
    cfl_s = max(cfl_s, cs);
    cfl_z = max(cfl_z, cz);
end

% Plot saturation CFL
figure; 
plotCellData(G, log10(cfl_s),'EdgeAlpha',.1)
ch = colorbar;
axis tight off
set(ch, 'YTick', -2:2, 'YTickLabel', ...
    arrayfun(@(x) sprintf('10^{%d}', x), -2:2, 'Unif', false));
% colormap(parula.^2);
title('Saturation CFL');

% Plot composition CFL
figure; 
plotCellData(G, log10(max(cfl_z, [], 2)),'EdgeAlpha',.1)
ch = colorbar;
axis tight off
set(ch, 'YTick', -3:2, 'YTickLabel', ...
    arrayfun(@(x) sprintf('10^{%d}', x), -3:2, 'Unif', false));
% colormap(parula.^2);
title('Composition CFL');

%%
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
