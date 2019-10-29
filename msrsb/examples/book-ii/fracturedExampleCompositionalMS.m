mrstModule add ad-core ad-blackoil deckformat msrsb
mrstModule add ad-props mrst-gui blackoil-sequential
mrstModule add coarsegrid
mrstModule add compositional
mrstModule add linearsolvers
%% Build grid
% This example is similar to the example from Section 5.2 of Moyner &
% Tchelepi, SPE J, 23(6), 2018. However, we use a simpler three-component
% fluid model, somewhat higher fracture permeability and slightly different
% well positions. See fracturedExampleMS for more details on the setup.

fluid_name = 'simple';
caseName = ['FractureMS', '_', fluid_name];

load setup_fracture

rdim = [1000 500];
G.nodes.coords = G.nodes.coords.*rdim;
G = computeGeometry(G);
rock = makeRock(G, perm*milli*darcy, 0.3);
rock.perm(G.cells.tag > 0) = 10*darcy;
pv = poreVolume(G, rock);
[cf, info] = getCompositionalFluidCase(fluid_name);
eos = EquationOfStateModel(G, cf);
minP = 50*barsa;
resP = info.pressure;
totTime = 7*year;
irate = 0.25*sum(pv)/totTime;

icell = findEnclosingCell(G, [0.025, 0.05].*rdim);
pcell = findEnclosingCell(G, [0.975, 0.95].*rdim);

W = [];
W = addWell(W, G, rock, icell, 'comp_i', [0, 1], 'name', 'Injector', 'Type', 'rate', 'Val', irate);
W = addWell(W, G, rock, pcell, ...
    'comp_i', [0.5, 0.5], 'Name', 'Producer', 'Val', minP, 'sign', -1);

fluid = initSimpleADIFluid('rho', [1000, 500, 500], ...
                       'mu', [1, 1, 1]*centi*poise, ...
                       'n', [2, 2, 2], ...
                       'c', [1e-5, 0, 0]/barsa);
model = OverallCompositionCompositionalModel(G, rock, fluid, eos, 'water', false);
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
s0 = [1, 0];
state0 = initResSol(G, resP, s0);

state0.T = repmat(info.temp, G.cells.num, 1);    
state0.components = repmat(info.initial, G.cells.num, 1);
for i = 1:numel(W)
    W(i).components = info.injection;
end
dt = rampupTimesteps(totTime, 20*day);
schedule = simpleSchedule(dt, 'W', W);
% Build packing utility for storing simulation setup
packer = @(model, name, varargin) packSimulationProblem(state0, model, schedule, caseName, 'name', name, varargin{:});
%% Define fully-implicit solver
nls = NonLinearSolver();
nls.LinearSolver = selectLinearSolverAD(model);
nls.useRelaxation = true;
if isa(nls.LinearSolver, 'AMGCL_CPRSolverAD')
    nls.LinearSolver.strategy = 'amgcl';
end
base = packer(model, 'Fully-implicit', 'NonLinearSolver', nls);
%% Setup transport models
ord = getCellMajorReordering(G.cells.num, size(state0.components, 2));

psolve = AMGCLSolverAD('preconditioner', 'amg', 'tolerance', 1e-5);
tsolve = AMGCLSolverAD('preconditioner', 'relaxation', 'relaxation', 'ilu0', ...
                       'tolerance', 1e-4, 'maxIterations', 25);
tsolve.equationOrdering = ord;
tsolve.variableOrdering = ord;

seqmodel = getSequentialModelFromFI(model);
seqmodel.transportModel.AutoDiffBackend = model.AutoDiffBackend;
seqmodel.transportNonLinearSolver.useRelaxation = true;

seqmodel.transportNonLinearSolver.LinearSolver = tsolve;
seqmodel.pressureNonLinearSolver.LinearSolver = psolve;

seq = packer(seqmodel, 'Sequential');
%% Build coarse grid and multiscale solver
pr     = sampleFromBox(G, reshape(1:100,[10 10]));
pr     = processPartition(G, compressPartition(pr));
CG    = generateCoarseGrid(G, compressPartition(pr));
CG    = coarsenGeometry(CG);
CG    = storeInteractionRegion(CG);
CG    = setupMexInteractionMapping(CG);

msmodel = getSequentialModelFromFI(model);
msmodel.transportModel.AutoDiffBackend = model.AutoDiffBackend;
% We the convergence criterion to be a factor 50 looser than the default of
% 1e-3. The linear solver is two orders of magnitude more relaxed than the
% strict AMG settings used above.
msmodel.pressureModel.incTolPressure = 0.05;
msmodel = addMultiscaleSolverComp(msmodel, CG, 'maxIterations', 50, ...
                                               'useGMRES', true, ...
                                               'tolerance', 1e-3);
msmodel.transportNonLinearSolver.useRelaxation = true;
msmodel.transportNonLinearSolver.LinearSolver = tsolve;

ms = packer(msmodel, 'Multiscale (Relaxed)');
%% Solve a stricter multiscale model with fine-scale tolerances
strictmsmodel = getSequentialModelFromFI(model);
strictmsmodel.transportNonLinearSolver.useRelaxation = true;
strictmsmodel.transportNonLinearSolver.LinearSolver = tsolve;

strictmsmodel.transportModel.AutoDiffBackend = model.AutoDiffBackend;
strictmsmodel.pressureModel.incTolPressure = 1e-3;

strictmsmodel = addMultiscaleSolverComp(strictmsmodel, CG, 'maxIterations', 50, ...
                                                    'useGMRES', true, 'tolerance', 1e-5);
strictms = packer(strictmsmodel, 'Multiscale (Strict)');
%% Simulate the problems and get the output
problems = {base, seq, ms, strictms};
simulatePackedProblem(problems, 'continueOnError', false);%, 'restartStep', 1);
[ws, states, reports, names, T] = getMultiplePackedSimulatorOutputs(problems, 'readFromDisk', false, 'readWellSolsFromDisk', true);
%% Plot well solutions
plotWellSols(ws, T, 'datasetnames', names, 'linestyles', {'-'}, 'SelectedWells', 2, 'field', 'qGs');
%% Plot the iterations taken in different parts of the solver
ns = numel(reports);
stats = cell(1, ns);
for i = 1:ns
    stats{i} = getPressureTransportIterations(reports{i});
end

colors = lines(ns);
total = zeros(ns, 3);
for i = 1:ns
    s = stats{i};
    total(i, 1) = sum(s.pressure);
    total(i, 2) = sum(s.transport);
    total(i, 3) = sum(s.outer);
end
figure;
bar(total)
set(gca, 'XTickLabel', names);
legend('Pressure', 'Transport', 'Outer');
%% Plot initial state
G.cells.sortedCellNodes = getSortedCellNodes(G);

clf;
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
            title(names{i},'FontSize',12,'FontWeight','normal');
        else
            set(gca, 'XTickLabel', []);
            set(gca, 'YTickLabel', []);
        end
    end
end
%% Plot the solution (gas saturation and compositions)
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