useWater = true;
caseName = 'verysimple';
p = 150*barsa;
bhp = 100*barsa;
pvi = 0.2;
nstep = 30;
time = 1*year;
n = 10;

[G, rock, state0, schedule, fluid, eos] = ...
    getTest1D_compositional(n, 'water', useWater, ...
                               'fluid', caseName, ...
                               'nkr', 2, ...
                               'nstep', nstep, ...
                               'pvi', pvi, ...
                               'pressure', p, ...
                               'bhp', bhp, ...
                               'perm', 0.5*darcy, ...
                               'time', time);
G = computeCellDimensions2(G);
G.cells.ghost = false(G.cells.num,1);

args = {G, rock, fluid, eos.fluid, 'water', useWater, 'useIncTolComposition', true};

pmodel = PressureNaturalVariablesModel(args{:});
tmodel = TransportNaturalVariablesModel(args{:});
modelSeq = SequentialPressureTransportModel(pmodel, tmodel);

ix = 1:30;
subschedule = schedule;
subschedule.step.val     = subschedule.step.val(ix);
subschedule.step.control = subschedule.step.control(ix);

%%

[wsFV, stFV, rep] = simulateScheduleAD(state0, modelSeq, subschedule);

%%

close all
plotToolbar(G, stFV, 'plot1d', true)

%%

degree = 0;
[jt, ot, mt] = deal(Inf);
mt = 0.0;
% pmodel = PressureNaturalVariablesModelSemiDG(args{:});
pmodel = PressureNaturalVariablesModel(args{:});
tmodel = TransportNaturalVariablesModelDG(args{:});
modelDG = SequentialPressureTransportModel(pmodel, tmodel);
disc = DGDiscretization(tmodel                                        , ...
                                'degree'               , degree       , ...
                                'basis'                , 'legendre'   , ...
                                'useUnstructCubature'  , true         , ...
                                'jumpTolerance'        , jt           , ...
                                'outTolerance'         , ot           , ...
                                'outLimiter'           , 'orderReduce', ...
                                'meanTolerance'        , mt           , ...
                                'limitAfterConvergence', false        , ...
                                'plotLimiterProgress'  , false        );
modelDG.transportModel = TransportNaturalVariablesModelDG(G, rock, fluid, eos.fluid, ...
                                   'disc'    , disc    , ...
                                   'dsMaxAbs', 0.1     , ...
                                   'nonlinearTolerance', 1e-3, ...
                                   'useIncTolComposition', true);
% modelDG.pressureModel.disc = disc;

state0 = assignDofFromState(disc, state0);
[wsDG, stDG, rep] = simulateScheduleAD(state0, modelDG, subschedule);

%%

close all
figure; plotToolbar(G, stFV, 'plot1d', true)
figure; plotToolbar(G, stDG, 'plot1d', true)

%%
plotWellSols({wsFV, wsDG})

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
