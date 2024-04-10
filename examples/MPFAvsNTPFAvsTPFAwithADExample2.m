%% Combining consistent discretizations with AD-OO
% We follow example 6.1.2 in the MRST book (see
% examples/1ph/showInconsistentTPFA in the book module).
% We create a skewed grid with pressure BCs.
clear all
close all

mrstModule add ad-core mpfa ad-blackoil compositional ad-props mrst-gui nfvm

dims = [21, 10];
G = cartGrid(dims, [2, 1]);
makeSkew = @(c) c(:, 1) + .4 * (1 - (c(:, 1) - 1).^2) .* (1 - c(:, 2));
G.nodes.coords(:, 1) = 2 * makeSkew(G.nodes.coords);
G.nodes.coords(:, 1) = G.nodes.coords(:, 1) * 1000;
G.nodes.coords(:, 2) = G.nodes.coords(:, 2) * 1000;
G = twister(G, 0.1);
G = computeGeometry(G);

% Homogeneous reservoir properties
rock = makeRock(G, 100*milli*darcy, .2);
pv = sum(poreVolume(G, rock));

% Pressure BCs
p0 = 10;
p1 = 20;
bc = [];
bc = pside(bc, G, 'Xmin', p0, 'sat', [0, 1]);
bc = pside(bc, G, 'Xmax', p1, 'sat', [0, 1]);

%% We can simulate with either immiscible or compositional fluid physics
% The example uses the general simulator framework and as such we can
% easily simulate the same problem with different underlying physics.

gravity reset off;

fluid = initSimpleADIFluid('cR', 1e-8/barsa, 'rho', [1, 1000, 100]);
if ~exist('useComp', 'var')
    useComp = false;
end

if useComp
    % Compositional, two-component
    [f, info] = getCompositionalFluidCase('verysimple');
    eos = EquationOfStateModel(G, f);
    model = GenericOverallCompositionModel(G, rock, fluid, eos, 'water', false);
    for i = 1:numel(W)
        W(i).components = info.injection;
    end
    z0 = info.initial;
    state0 = initCompositionalState(G, info.pressure, info.temp, [1, 0], z0, eos);
    W(1).val = 100 * W(1).val;
else
    % Immiscible two-phase
    model = GenericBlackOilModel(G, rock, fluid, 'water', true, 'oil', true, 'gas', false);
    state0 = initResSol(G, 1*barsa, [0, 1]);
end
% Schedule
dt = [1; 9; repmat(15, 26, 1)] * day;
schedule = simpleSchedule(dt, 'W', [], 'bc', bc);

%% Simulate the implicit TPFA base case
disp('TPFA implicit')
[wsTPFA, statesTPFA] = simulateScheduleAD(state0, model, schedule);
plotFinalPressure(G, statesTPFA, 'TPFA')

%% Simulate implicit AvgMPFA
disp('AvgMPFA implicit')
mrstModule add nfvm
ratio = [];
model_avgmpfa = setAvgMPFADiscretization(model, 'myRatio', ratio);
[wsAvgMPFA, statesAvgMPFA] = simulateScheduleAD(state0, model_avgmpfa, schedule);
plotFinalPressure(G, statesAvgMPFA, 'AvgMPFA')

%% Simulate implicit NTPFA
disp('NTPFA implicit')
mrstModule add nfvm
ratio = [];
model_ntpfa = setNTPFADiscretization(model, 'myRatio', ratio);

% Illustrating options
model_ntpfa.nonlinearTolerance = 1e-14;
model_ntpfa.toleranceCNV = 1e-14;
model_ntpfa.toleranceMB = 1e-14;
model_ntpfa.verbose = true;

[wsNTPFA, statesNTPFA] = simulateScheduleAD(state0, model_ntpfa, schedule);
plotFinalPressure(G, statesNTPFA, 'NTPFA')

%% Simulate implicit MPFA
% The simulator reuses the multipoint transmissibility calculations from
% the MPFA module. We instantiate a special phase potential difference that
% is computed using MPFA instead of the regular two-point difference for
% each face.
disp('MPFA implicit')
mrstModule add mpfa
model_mpfa = setMPFADiscretization(model);
[wsMPFA, statesMPFA] = simulateScheduleAD(state0, model_mpfa, schedule);
plotFinalPressure(G, statesMPFA, 'MPFA')

%% Simulate and plot with explicit time discretization
model_exp = setTimeDiscretization(model, 'Explicit', 'initialStep', dt(1));
model_avgmpfa_exp = setTimeDiscretization(model_avgmpfa, 'Explicit');
model_ntpfa_exp = setTimeDiscretization(model_ntpfa, 'Explicit');
model_mpfa_exp = setTimeDiscretization(model_mpfa, 'Explicit');

disp('TPFA explicit')
[wsExplicit, statesExplicit] = simulateScheduleAD(state0, model_exp, schedule);
plotFinalPressure(G, statesExplicit, 'TPFA explicit')

disp('AvgMPFA explicit')
[wsAvgMPFAExplicit, statesAvgMPFAExplicit] = simulateScheduleAD(state0, model_avgmpfa_exp, schedule);
plotFinalPressure(G, statesAvgMPFAExplicit, 'AvgMPFA explicit')

disp('NTPFA explicit')
[wsNTPFAExplicit, statesNTPFAExplicit] = simulateScheduleAD(state0, model_ntpfa_exp, schedule);
plotFinalPressure(G, statesNTPFAExplicit, 'NTPFA explicit')

disp('MPFA explicit')
[wsMPFAExplicit, statesMPFAExplicit] = simulateScheduleAD(state0, model_mpfa_exp, schedule);
plotFinalPressure(G, statesMPFAExplicit, 'MPFA explicit');

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
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
