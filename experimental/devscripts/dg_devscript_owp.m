mrstModule add dg vem vemmech ad-props ad-core ad-blackoil ...
    blackoil-sequential gasinjection mrst-gui reorder matlab_bgl ...
    ad-eor trust-region
mrstVerbose on

%%

gravity reset off

n = 50;
l = 1000*meter;
G = computeGeometry(cartGrid([n,1], [1, 0.01]*l));
G = computeVEMGeometry(G);
G = computeCellDimensions(G);

rock  = makeRock(G, 100*milli*darcy, 1);
fluid = initSimpleADIFluid('phases', 'WO'                   , ...
                           'rho'   , [1000, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise     , ...
                           'n'     , [2, 2]                 );
fluid = addSimplePolymerProperties(fluid);

model   = OilWaterPolymerModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(model);
modelDG = modelFV;

%%

time = 0.5*year;
rate = 0.5*sum(poreVolume(G, rock))/time;
W = [];
W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);
[W.c] = deal(0);

dt       = 10*day;
dtvec    = rampupTimesteps(time, dt, 0);
schedule_w = simpleSchedule(dtvec, 'W', W);
schedule_p = schedule_w;
schedule_p.step.control = schedule_p.step.control + 1;
schedule = schedule_w;
schedule.step.val     = [schedule_w.step.val; ...
                         schedule_p.step.val; ...
                         schedule_w.step.val];
schedule.step.control = [schedule_w.step.control; ...
                         schedule_p.step.control; ...
                         schedule_w.step.control];
schedule.control(2) = schedule.control(1);
[schedule.control(2).W.c] = deal(fluid.cmax);
% [schedule.control(1).W.c] = deal(fluid.cmax);

sW          = 0.0;
state0      = initResSol(G, 100*barsa, [sW,1-sW]);
state0.c    = zeros([model.G.cells.num, 1]);
state0.cmax = zeros([model.G.cells.num, 1]);

%%

ix = ':';
% ix = 1:10;
subschedule = schedule;
subschedule.step.val = subschedule.step.val(ix);
subschedule.step.control= subschedule.step.control(ix);

%%

[wsFV, stFV, repFV] = simulateScheduleAD(state0, modelFV, subschedule);

%%

degree = [0];

[jt, ot, mt] = deal(Inf);
% 
jt = 0.2;
% mt = 0.0;
ot = 0.01;

[wsDG, stDG] = deal(cell(numel(degree),1));
for dNo = 1:numel(degree)
    modelDG.transportModel = TransportOilWaterPolymerModelDG(G, rock, fluid, ...
                                             'dsMaxAbs'     , 0.2        , ...
                                             'degree'       , degree(dNo), ...
                                             'meanTolerance', 1e-3       , ...
                                             'outTolerance' , 1e-3       );

    state0 = assignDofFromState(modelDG.transportModel.disc, state0);
    [wsDG{dNo}, stDG{dNo}, rep] = simulateScheduleAD(state0, modelDG, subschedule);
    
end

%%

plotWellSols({wsFV, wsDG{:}});

%%

close all
figure
d = cellfun(@(s1, s2) compareStates(s1, s2), stFV, stDG{1}, 'unif', false);
plotToolbar(G, d, 'plot1d', true)
figure
plotToolbar(G, stFV, 'plot1d', true)
figure
plotToolbar(G, stDG{1}, 'plot1d', true)

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
