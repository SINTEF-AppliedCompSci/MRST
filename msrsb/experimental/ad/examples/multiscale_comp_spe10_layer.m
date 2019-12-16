mrstModule add mrst-gui ad-props ad-core ad-blackoil blackoil-sequential msrsb spe10 coarsegrid
%% Set up simulation problem
totTime = 2000*day;

layers = 5;
nLayers = numel(layers);

% Offsets to move the wells if needed
ofs = [0, 0];
wloc     = [  1 + ofs(1),   60 - ofs(1),     1 + ofs(1),   60 - ofs(1),   30 ;
              1 + ofs(2),    1 + ofs(2),   220 - ofs(2),  220 - ofs(2),  110 ];

[G, W, rock] = getSPE10setup(layers, wloc);

mp = 0.01;
rock.poro(rock.poro < mp) = mp;

pv = poreVolume(G, rock);

for i = 1:numel(W)
    isinj = W(i).val > 400*barsa;
    if isinj
        W(i).val = 0.25*sum(pv)/totTime;
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
fluid = initSimpleADIFluid('mu', [1, 5, 1]*centi*poise, 'rho', [1000, 700, 1], 'n', [2 2 2], 'c', [1e-6, 1e-6, 1e-6]/barsa);

% Fully implicit model
modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
% Sequential pressure-transport model with same type
model = getSequentialModelFromFI(modelfi);
model.transportModel.useCNVConvergence = modelfi.useCNVConvergence;


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


% state.wellSol = initWellSolAD(W, modelfi, state);

% lim = 2;
% schedule.step.val     = schedule.step.val(1:lim);
% schedule.step.control = schedule.step.control(1:lim);
%%
% split, baseline
timer = tic();
[ws_split, states_split, report_split] = simulateScheduleAD(state, model, schedule);
t_split = toc(timer);

%%
[wsfi, statesfi] = simulateScheduleAD(state, modelfi, schedule);
%% Plot well curves for fully implicit problem
plotWellSols({ws_split, wsfi}, cumsum(schedule.step.val), 'datasetnames', {'Sequential', 'FI'})

%% Set up coarse grid
mrstModule add multiscale-devel coarsegrid mrst-experimental
coarsefactors = [10, 10, 1];
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

tols = 0.01;
its  = 100;

multiscaleModels = cell(numel(tols), 1);

for i = 1:numel(tols)
    it = its(i);
    tol = tols(i);
    
    mssolver = MultiscaleVolumeSolverAD(CG, ...
                                            'maxIterations', it, ...
                                            'useGMRES', true, ...
                                            'getSmoother', s, ...
                                            'tolerance', tol, ...
                                            'useMex', false, ...
                                            'prolongationType', 'msrsb');

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
    [wellsols{i}, states, reports{i}] = simulateScheduleAD(state, model2, schedule);
    timing(i) = toc(timer);
end
%%

% close all
plotWellSols({ws_split, wellsols{1}}, cumsum(schedule.step.val), 'datasetnames', {'Sequential', 'Multiscale'})

%%
plotWellSols({ws_split, wellsols{1}, wsfi}, cumsum(schedule.step.val), 'datasetnames', {'Sequential', 'Multiscale', 'FI'})

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.
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
