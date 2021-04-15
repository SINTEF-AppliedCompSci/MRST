mrstModule add deckformat mrst-gui ad-core ad-blackoil ad-props
mrstModule add sequential ad-unittest multiscale-devel
mrstModule add sequential spe10

%%
totTime = 2000*day;
gravity reset on

layers = 1:10;
nLayers = numel(layers);


coarsefactors = [12 20 5];
ofs = ceil(coarsefactors(1:2)/2);
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

fn = mrstPath('query', 'ad-testdata');
deck = readEclipseDeck(fullfile(fn, 'SINTEF', 'simpleOW', 'simple10x1x10.data'));
deck = convertDeckUnits(deck);

fluid = initDeckADIFluid(deck);

fn = @(p) fluid.BOxmuO(p).*fluid.bO(p);

fluid = rmfield(fluid, {'BOxmuO', 'pcOW', 'sWcon'});
fluid.muO = fn;

fluid.krOW = @(s, varargin) s.^2;
fluid.krW  = @(s, varargin) s.^2;


p0 = 300*barsa;
fluid.bO = @(p) 0.5 + (p-p0)./(1000*barsa);

modelfi = TwoPhaseOilWaterModel(G, rock, fluid);


if ispc
    datadir = 'D:\jobb\data\';
else
    datadir = '/data/';
end

cdims = ceil(G.cartDims./coarsefactors);

p = partitionUI(G, cdims);
p = compressPartition(p);

figure;
plotCellData(G, log10(rock.perm(:, 1)), 'edgec', 'none');
outlineCoarseGrid(G, p)
plotWell(G, W)
%%

nstep = 50;

val = repmat(totTime/nstep, nstep, 1);

clear schedule;
schedule.step.val = val;
schedule.step.control = ones(size(val));
schedule.control.W = W;

state = initResSol(G, p0, [.2, .8]);


state.wellSol = initWellSolAD(W, modelfi, state);

% lim = 1;
% schedule.step.val     = schedule.step.val(1:lim);
% schedule.step.control = schedule.step.control(1:lim);
%%
mrstModule add agmg
mrstVerbose on
clear solver


outersolver = NonLinearSolver('enforceResidualDecrease', false, ...
                         'useRelaxation', true, ...
                         'verbose', true);
outersolver.maxIterations = 10;

outersolver.timeStepSelector = IterationCountTimeStepSelector();
outersolver.timeStepSelector.maxTimestep = 90*day;
outersolver.timeStepSelector.firstRampupStep = 1*day;
outersolver.timeStepSelector.targetIterationCount = 5;
  


% mrstModule add agmg
% solverp.LinearSolver = AGMGSolverAD();

mrstVerbose on
model = getSequentialModelFromFI(modelfi);

model.transportNonLinearSolver.useRelaxation = true;
model.pressureNonLinearSolver.useRelaxation = true;

model.transportNonLinearSolver.LinearSolver = GMRES_ILUSolverAD();
% model.pressureNonLinearSolver.LinearSolver = BackslashSolverAD();
model.pressureNonLinearSolver.LinearSolver = AGMGSolverAD();
model.pressureNonLinearSolver.LinearSolver.tolerance = 1e-5;
model.pressureNonLinearSolver.LinearSolver.maxIterations = 100;

model.transportModel.extraStateOutput = false;
model.pressureModel.extraStateOutput = false;

model.pressureModel.nonlinearTolerance = 5e-6;
model.transportModel.nonlinearTolerance = 5e-6;
% model.transportModel.conserveOil = false;
% model.transportModel.conserveWater = true;
%%
% split, baseline
timer = tic();
[ws_split, states_split, report_split] = simulateScheduleAD(state, model, schedule, 'NonLinearSolver', outersolver, ... 
                                         'outputMinisteps', true);
t_split = toc(timer);

schedule_ms = convertReportToSchedule(report_split, schedule);

% save(fullfile(datadir, 'spe10_base.mat'), 'ws_split', 'report_split', 'states_split', '-v7.3')


% schedule_ms.step.val = schedule.step.val(1);
% schedule_ms.step.control = schedule.step.control(1);

%%
mrstModule add multiscale-devel coarsegrid mrst-experimental



CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegionCart(CG);
%%
rehash
clear mssolver msmodel MultiscalePressureModel

if 1
    s = 'ilu';
    smooth = 1;
else
    s = 'jacobi';
    smooth = 10;
end

tols = 0.1;
its  = 100;

multiscaleModels = cell(numel(tols), 1);

for i = 1:numel(tols)
    it = its(i);
    tol = tols(i);
    
    mssolver = MultiscaleVolumeSolverAD(CG, 'useGalerkinRestriction', false,...
                                            'maxIterations', it, ...
                                            'useGMRES', true, ...
                                            'smoother', s, ...
                                            'tolerance', tol, ...
                                            'prolongationType', 'jacobi-mex', ...
                                            'smoothCycleIterations', smooth);

    msmodel = MultiscalePressureModel(G, rock, model.fluid, model.pressureModel, mssolver, 'extraStateOutput', false);

    msmodel.nonlinearTolerance = 1e-6;
    multiscaleModels{i} = msmodel;
end

nms = numel(multiscaleModels);
[wellsols, states, reports] = deal(cell(nms, 1));
timing = nan(nms, 1);
for i = 1:nms
    
    solver = NonLinearSolver('enforceResidualDecrease', false, ...
                             'useRelaxation', true, ...
                             'verbose', true);

    model2 = model;
    model2.pressureModel = multiscaleModels{i};
    % Reset basis
    model2.pressureModel.multiscaleSolver.basis = [];
    
    timer = tic();
    [wellsols{i}, states, reports{i}] = simulateScheduleAD(state, model2, schedule_ms, 'NonLinearSolver', solver);
    timing(i) = toc(timer);
end
%%

% save(fullfile(datadir, 'spe10_out_ms2.mat'), 'wellsols', 'ws_split', 'reports', 'report_split', '-v7.3')
% 
%%
close all
plotWellSols({ws_split, wellsols{1}}, cumsum(schedule_ms.step.val), 'datasetnames', {'Sequential', 'Multiscale'})


t1 = breakDownSolver(report_split);
t2 = breakDownSolver(reports{1});
figure;
subplot(1, 2, 1); 
bar([t1; t2]/3600)
set(gca, 'XTickLabel', {'Sequential', 'Multiscale'})
legend('Pressure', 'Transport', 'Other')
title('Timing for full SPE10 (two-phase, black-oil)')
ylabel('Time taken (hours, wall time)')
grid on

subplot(1, 2, 2); 
d = [sum(t1(1:2))./sum(t2(1:2)), sum(t1)/sum(t2)];
bar(d)
set(gca, 'XTickLabel', {'Solver', 'Simulation'})
title('Speedup')
%%
clc
disp(['Pressure:', formatTimeRange(t1(1)), '| ', formatTimeRange(t2(1))])
disp(['Transport:', formatTimeRange(t1(2)), '| ', formatTimeRange(t2(2))])
disp(['Other:', formatTimeRange(t1(3)), '| ', formatTimeRange(t2(3))])
disp(['Total:', formatTimeRange(sum(t1)), '| ', formatTimeRange(sum(t2))])

%%
clc
disp(['Speedup pressure: ', num2str(sum(t1(1))/sum(t2(1)))])
disp(['Speedup solver: ', num2str(sum(t1(1:2))/sum(t2(1:2)))])
disp(['Speedup simulation: ', num2str(sum(t1(1:3))/sum(t2(1:3)))])

%%
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
