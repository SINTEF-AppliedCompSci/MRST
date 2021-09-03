%% Example demonstrating compositional solvers with K-values
% We solve SPE1 using two compositional models with the black-oil
% properties represented as pressure and composition-dependent
% equilibrium constants.

%% First, setup the base case
mrstModule add ad-core ad-blackoil compositional ad-props mrst-gui example-suite
[G, rock, fluid, deck, state] = setupSPE1();

model = ThreePhaseBlackOilModel(G, rock, fluid, 'disgas', true);
model.extraStateOutput = true;
schedule = convertDeckScheduleToMRST(model, deck);

%% Switch well limits
% The compositional solvers use a different definition of oil-rate
% controls. For this reason, we adjust the control of the producer to use
% BHP.
schedule.control.W(2).type = 'bhp';
schedule.control.W(2).val = 200*barsa;
schedule.control.W(2).lims = [];

%% Convert state and wells
state0 = convertBlackOilStateToCompositional(model, state);
for i = 1:numel(schedule.control(1).W)
    schedule.control(1).W(i).components = schedule.control(1).W(i).compi(2:end);
end
%% Convert model to compositional model
% This is based on interpolation and is not completely robust.
modelNat = convertBlackOilModelToCompositionalModel(model, 'interpolation', 'pchip');
modelNat.FacilityModel = FacilityModel(modelNat);
modelNat.FacilityModel.toleranceWellRate = 1e-3;
modelNat.dsMaxAbs = 0.2;
modelNat.dzMaxAbs = inf;
modelNat.incTolPressure = inf;

%% Run natural variables
[wsNat, statesNat, reportNat] = simulateScheduleAD(state0, modelNat, schedule);
%% Run overall composition
modelOverall = OverallCompositionCompositionalModel(G, rock, fluid, modelNat.EOSModel);
modelOverall.FacilityModel = FacilityModel(modelOverall);
modelOverall.FacilityModel.toleranceWellRate = modelNat.FacilityModel.toleranceWellRate;

modelOverall.dsMaxAbs = modelNat.dsMaxAbs;
modelOverall.dzMaxAbs = modelNat.dzMaxAbs;
modelOverall.incTolPressure = modelNat.incTolPressure;


[wsOverall, statesOverall, reportOverall] = simulateScheduleAD(state0, modelOverall, schedule);
%% Run the reference
model.FacilityModel = FacilityModel(model);
model.FacilityModel.toleranceWellRate = modelNat.FacilityModel.toleranceWellRate;

[wsRef, statesRef, reportRef] = simulateScheduleAD(state, model, schedule);
%% Store component rates for the reference solution in the same format
for i = 1:numel(wsRef)
    for j = 1:numel(wsRef{i})
        wsRef{i}(j).PseudoGas = wsRef{i}(j).qGs.*fluid.rhoGS;
        wsRef{i}(j).PseudoOil = wsRef{i}(j).qOs.*fluid.rhoOS;
    end
end
%% Add interactive plotting
% Note: Gas and oil rates in the producer are defined differently between
% the solvers, and will not match
names = {'MRST-NaturalVariables', 'MRST-MolarVariables', 'MRST-BlackOil'};
figure;
plotToolbar(G, statesNat);
title('Compositional Natural Variables')
axis tight
view(30, 30);
figure;
plotToolbar(G, statesOverall);
title('Compositional Molar Variables')
axis tight
view(30, 30);

figure;
plotToolbar(G, statesRef);
title('Black-oil');
axis tight
view(30, 30);
plotWellSols({wsNat, wsOverall, wsRef}, 'datasetnames', names)
%% Plot nonlinear iterations
figure;
stairs(cumsum([reportNat.Iterations, reportOverall.Iterations, reportRef.Iterations]));
legend(names);
ylabel('Cumulative nonlinear iterations')

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
