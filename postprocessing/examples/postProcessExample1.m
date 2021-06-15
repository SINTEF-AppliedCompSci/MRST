mrstModule add ad-core ad-blackoil postprocessing ad-props
mrstModule add mrst-gui
mrstModule add diagnostics

%% Set up five-spot
dims = [25, 25, 1];
time = 5*year;
G = cartGrid(dims, [1000, 1000, 10]);
G = computeGeometry(G);
rock = makeRock(G, 0.1*darcy, 0.3);

fluid = initSimpleADIFluid('phases', 'WO', 'n', [2, 2], 'mu', [1, 2]*centi*poise);

model = TwoPhaseOilWaterModel(G, rock, fluid);

irate = sum(model.operators.pv)/(4*time);
W = [];
W = verticalWell(W, G, rock, 1, 1, [], 'type', 'rate', 'val', irate, 'comp_i', [1, 0]);
W = verticalWell(W, G, rock, dims(1), 1, [], 'type', 'rate', 'val', irate, 'comp_i', [1, 0]);
W = verticalWell(W, G, rock, 1, dims(2), [], 'type', 'rate', 'val', irate, 'comp_i', [1, 0]);
W = verticalWell(W, G, rock, dims(1), dims(2), [], 'type', 'rate', 'val', irate, 'comp_i', [1, 0]);

W = verticalWell(W, G, rock, ceil(dims(1)/2), ceil(dims(2)/2), [], 'type', 'bhp', 'val', 100*barsa, 'comp_i', [1, 0]);

% Put different tracer in all wells
for i = 1:5
    W(i).tracer = zeros(1, 4);
    if i < 5
        W(i).tracer(i) = 1;
    end
end

state0 = initResSol(G, 150*barsa, [0, 1]);
schedule = simpleSchedule(rampupTimesteps(time, time/50), 'W', W);
[ws, states, report] = simulateScheduleAD(state0, model, schedule);
%% Create a model which simply checks convergence
verifyModel = CheckConvergenceModel(model);
[ws2, states2, report2] = simulateScheduleAD(state0, verifyModel, schedule, 'initialGuess', states);
%% Create and solve tracers
tmodel = TracerModel(G, rock, fluid, 'tracerNames', {'T1', 'T2', 'T3', 'T4'});
dmodel = DiagnosticsPostProcessModel(G, rock, fluid);
postModel = PostProcessModel(verifyModel, {tmodel, dmodel});
[ws3, states3, report3] = simulateScheduleAD(state0, postModel, schedule, 'initialGuess', states);
%% Plot tracer
figure; plotToolbar(G, states3, 'field', 'tracer:1')
%% Plot tracer blend using diagnostics tools
tracer = states3{end}.tracer;
[t, tpart] = max(tracer, [], 2);
figure;
plotTracerBlend(G, tpart, t);

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
