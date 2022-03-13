mrstModule add dg vem vemmech ad-props ad-core ad-blackoil blackoil-sequential gasinjection mrst-gui reorder matlab_bgl

%%

gravity reset off

n = 10;
l = 1000*meter;
G = computeGeometry(cartGrid([n,1], [1,0.1]*l));
G = computeVEMGeometry(G);
G = computeCellDimensions(G);

rock = makeRock(G, 100*milli*darcy, 1);
fluid = initSimpleADIFluid('phases', 'WO'                   , ...
                           'rho'   , [1000, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise     , ...
                           'n'     , [1, 1] );
fluid.krW = @(s) fluid.krW(s).*(s >= 0 & s <= 1) + (s > 1);
fluid.krO = @(s) min(max(fluid.krO(s),0),1);

modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(modelfi);
modelDG = modelFV;

%%

time = 2*year;
rate = 2*sum(poreVolume(G, rock))/time;
W = [];
W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

dt    = 20*day;
dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W);

sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);

%%

degree = [2, 3, 4, 5];
jt = 0.2;
mt = 0.0;
ot = 0.0;

nls = NonLinearSolver('maxIterations', 25, 'useLinesearch', false);
[wsDG, statesDG, disc] = deal(cell(numel(degree),1));
for dNo = 1:numel(degree)
    disc{dNo} = DGDiscretization(modelDG.transportModel                    , ...
                                    'degree'             , degree(dNo), ...
                                    'basis'              , 'legendre' , ...
                                    'useUnstructCubature', false      , ...
                                    'jumpTolerance'      , jt         , ...
                                    'jumpLimiter'        , 'kill'      , ...
                                    'outTolerance'       , ot         , ...
                                    'outLimiter'         , 'orderReduce', ...
                                    'meanTolerance'      , mt         , ...
                                    'plotLimiterProgress', false       , ...
                                    'limitAfterConvergence', true, ...
                                    'limitAfterNewtonStep', true);
    modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid  , ...
                                       'disc'    , disc{dNo}          , ...
                                       'dsMaxAbs', 0.05                , ...
                                       'nonlinearTolerance', 1e-3);
    modelDG.pressureModel = PressureOilWaterModelSemiDG(G, rock, fluid, ...
                                       'disc', disc{dNo}              , ...
                                       'extraStateOutput', true);
%     modelDG.transportNonLinearSolver = nls;
    modelDG = SequentialPressureTransportModelDG(modelDG.pressureModel, modelDG.transportModel);
                                    
    state0 = assignDofFromState(modelDG.transportModel.disc, state0);
    [wsDG{dNo}, statesDG{dNo}, rep] ...
        = simulateScheduleAD(state0, modelDG, schedule);
    
end

%%

[wsFV, statesFV, rep] = simulateScheduleAD(state0, modelFV, schedule);

%%

dsnDG = cellfun(@(d) ['dG(' num2str(d), ')'] , num2cell(degree), 'unif', false);

close all

for dNo = 1:numel(degree)
    figure('name', dsnDG{dNo});
    plotToolbar(G, statesDG{dNo}, 'plot1d', true);
    colormap(jet);
end

figure('name', 'FV');
plotToolbar(G, statesFV, 'plot1d', true);
colormap(jet)

dsn = horzcat('FV', dsnDG);

plotWellSols({wsFV, wsDG{:}}, schedule.step.val, 'datasetNames', dsn)
% plotWellSols({wsFV, wsDG{:}, wsDGReorder}, schedule.step.val)

%%

close all

ll = cellfun(@(d) ['dG(', num2str(d), ')'], num2cell(degree), 'unif', false);

cNo = round(G.cartDims(1)/2);
s = zeros(numel(schedule.step.val),numel(degree));
for dNo = 1:numel(degree)
    s(:,dNo) = cellfun(@(s) s.s(cNo,1), statesDG{dNo});
end

plot(s, 'linew', 2);
legend(ll, 'location', 'southeast');
    
%%
close all
figure('position', [-1000, 0, 800, 600])
% h = zeros(numel(degree),1);
saturation = cell(numel(degree), 1);
hold on
for dNo = 1:numel(degree)
    [h(dNo), saturation{dNo}, coords, keep, n] = ...
        plotSaturationDG(disc{dNo}, statesDG{dNo}{1}, 'plot1d', true, 'linew', 2);
end
ll = cellfun(@(d) ['dG(', num2str(d), ')'], num2cell(degree), 'unif', false);
legend(ll, 'location', 'southwest');
hold off
axis tight; box on
ylim([-0.1,1.1])

for sNo = 1:numel(statesDG{1})
    for dNo = 1:numel(degree)
        s = saturation{dNo}(statesDG{dNo}{sNo});
        set(h(dNo), 'YData', s);
    end
    pause(0.5);
end

%%

close all

sNo = 10;
hold on
for dNo = 1:numel(degree)
    [h(dNo), saturation{dNo}, coords, keep, n] = ...
        plotSaturationDG(disc{dNo}, statesDG{dNo}{sNo}, 'plot1d', true, 'linew', 2);
end
hold off
legend(ll)

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
