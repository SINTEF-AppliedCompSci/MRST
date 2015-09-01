mrstModule add mrst-gui ad-props ad-core ad-blackoil blackoil-sequential new-multiscale spe10 coarsegrid
%% Set up simulation problem
totTime = 2000*day;

layers = 5;
nLayers = numel(layers);

% Offsets to move the wells if needed
ofs = [0, 0];
wloc     = [  1 + ofs(1),   60 - ofs(1),     1 + ofs(1),   60 - ofs(1),   30 ;
              1 + ofs(2),    1 + ofs(2),   220 - ofs(2),  220 - ofs(2),  110 ];

[G, W, rock] = SPE10_setup(layers, wloc);

mp = 0.01;
rock.poro(rock.poro < mp) = mp;

pv = poreVolume(G, rock);

for i = 1:numel(W)
    isinj = W(i).val > 400*barsa;
    if isinj
        W(i).val = 3*sum(pv)/totTime;
        W(i).type = 'rate';
    else
        W(i).val = 4000*psia;
        W(i).type = 'bhp';
    end
    W(i).compi = [1 0];
    W(i).sign = 1 - 2*~isinj;
end

% Set up model with quadratic relperms and constant compressibility in oil
p0 = 300*barsa;
fluid = initSimpleADIFluid('mu', [1, 5, 1]*centi*poise, 'rho', [1000, 700, 1], 'n', [2 2 2]);
c = 0.001/barsa;
fluid_comp.bO = @(p) exp((p - p0)*c);

% Fully implicit model
modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
% Sequential pressure-transport model with same type
model = getSequentialModelFromFI(modelfi);


figure;
plotCellData(G, log10(rock.perm(:, 1)), 'edgec', 'none');
plotWell(G, W)
%% Set up the schedule
nstep = 50;

val = repmat(totTime/nstep, nstep, 1);

clear schedule;
schedule.step.val = val;
schedule.step.control = ones(size(val));
schedule.control.W = W;

state = initResSol(G, p0, [.2, .8]);


state.wellSol = initWellSolAD(W, modelfi, state);

% lim = 2;
% schedule.step.val     = schedule.step.val(1:lim);
% schedule.step.control = schedule.step.control(1:lim);
%%
% split, baseline
timer = tic();
[ws_split, states_split, report_split] = simulateScheduleAD(state, model, schedule, ... 
                                         'outputMinisteps', true);
t_split = toc(timer);

schedule_ms = convertReportToSchedule(report_split, schedule);
%%
[wsfi, statesfi] = simulateScheduleAD(state, modelfi, schedule_ms);
%% Plot well curves for fully implicit problem
plotWellSols({ws_split, wsfi}, cumsum(schedule_ms.step.val), 'datasetnames', {'Sequential', 'FI'})

%% Set up coarse grid
mrstModule add multiscale-devel coarsegrid mrst-experimental
coarsefactors = [10, 20, 1];
cdims = ceil(G.cartDims./coarsefactors);

p = partitionUI(G, cdims);
p = compressPartition(p);

CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegionCart(CG);
%% Set up and solve multiscale problem
if 1
    s = getSmootherFunction('type', 'ilu', 'iterations', 1);
else
    s = getSmootherFunction('type', 'jacobi', 'iterations', 10);
end

tols = 0.1;
its  = 5;

multiscaleModels = cell(numel(tols), 1);

for i = 1:numel(tols)
    it = its(i);
    tol = tols(i);
    
    mssolver = MultiscaleVolumeSolverAD(CG, ...
                                            'maxIterations', it, ...
                                            'useGMRES', true, ...
                                            'getSmoother', s, ...
                                            'tolerance', tol, ...
                                            'prolongationType', 'jacobi-mex');

    msmodel = MultiscalePressureModel(G, rock, model.fluid, model.pressureModel, mssolver);

    msmodel.nonlinearTolerance = 1e-6;
    multiscaleModels{i} = msmodel;
end

nms = numel(multiscaleModels);
[wellsols, states, reports] = deal(cell(nms, 1));
timing = nan(nms, 1);
for i = 1:nms
    model2 = model;
    model2.pressureModel = multiscaleModels{i};
    model2.stepFunctionIsLinear = true;
    model2.outerTolerance = model2.nonlinearTolerance;
    % Reset basis
    model2.pressureModel.multiscaleSolver.basis = [];
    
    timer = tic();
    [wellsols{i}, states, reports{i}] = simulateScheduleAD(state, model2, schedule_ms);
    timing(i) = toc(timer);
end
%%

% close all
plotWellSols({ws_split, wellsols{1}}, cumsum(schedule_ms.step.val), 'datasetnames', {'Sequential', 'Multiscale'})

%%
plotWellSols({ws_split, wellsols{1}, wsfi}, cumsum(schedule_ms.step.val), 'datasetnames', {'Sequential', 'Multiscale', 'FI'})

