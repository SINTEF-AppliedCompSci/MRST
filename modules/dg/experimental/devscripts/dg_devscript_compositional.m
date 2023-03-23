mrstModule add msrsb compositional dg gasinjection ad-blackoil ...
    blackoil-sequential ad-core ad-props coarsegrid mrst-gui vista ...
    coarsegrid msrsb vemmech vem
mrstVerbose on

%%

% caseName = 'simple_1d_wat';
% caseName = 'simple_3ph';
caseName = 'immiscible_denis';
[state0, modelFI, schedule, CG] = setupCompositionalPaperCase(caseName);
modelFI = NaturalVariablesCompositionalModel(modelFI.G, modelFI.rock, modelFI.fluid, modelFI.EOSModel.fluid, 'water', modelFI.water);
gravity reset off;

%%

G     = modelFI.G;
rock  = modelFI.rock;
fluid = modelFI.fluid;
compFluid = modelFI.EOSModel.fluid;
G = computeCellDimensions2(G);
G.cells.ghost = false(G.cells.num,1);
args = {G, rock, fluid, modelFI.EOSModel.fluid, 'water', modelFI.water, ...
                                                'nonlineartolerance', 1e-4};
useNat = true;
if useNat
    model = NaturalVariablesCompositionalModel(args{:});
    pmodel = PressureNaturalVariablesModel(args{:});
    tmodel = TransportNaturalVariablesModel(args{:});
else
    model = OverallCompositionCompositionalModel(args{:});
    pmodel = PressureOverallCompositionModel(args{:});
    tmodel = TransportOverallCompositionModel(args{:});
end
modelSeq = SequentialPressureTransportModel(pmodel, tmodel);

%%

ix = 1:45;%numel(schedule.step.val);
subschedule = schedule;
subschedule.step.val     = subschedule.step.val(ix);
subschedule.step.control = subschedule.step.control(ix);

%%

[wsFV, stFV, rep] = simulateScheduleAD(state0, modelSeq, subschedule);

%%

close all
plotToolbar(modelSeq.transportModel.G, stFV);
plotWellSols(wsFV);

%%

[jt, ot, mt] = deal(Inf);
mt = 0.0;
% pmodel = PressureNaturalVariablesModelSemiDG(args{:});
pmodel = PressureNaturalVariablesModel(args{:});
pmodel.incTolPressure = 1e-3;
tmodel = TransportNaturalVariablesModelDG(args{:});
degree = 1;
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
modelDG.transportModel = TransportNaturalVariablesModelDG(G, rock, fluid, compFluid, ...
                                   'disc'    , disc    , ...
                                   'dsMaxAbs', 0.1     , ...
                                   'dzMaxAbs', 0.1     , ...
                                   'water', modelFI.water, ...
                                   'nonlinearTolerance', 1e-4, ...
                                   'useIncTolComposition', false, ...
                                   'reduceLinearSystem', false);
% modelDG.pressureModel.disc = disc;

state0 = assignDofFromState(disc, state0);
state0.twoPhase = false(G.cells.num,1);
[wsDG, stDG, rep] = simulateScheduleAD(state0, modelDG, subschedule);

%%

figure
plotToolbar(modelDG.transportModel.G, stDG);
plotWellSols({wsFV, wsDG});

%%

stD = cellfun(@(s1, s2) compareStates(s1, s2), stFV, stDG, 'unif', false)

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.
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
