% Load modules
mrstModule add ad-core ad-blackoil compositional ad-props mrst-gui mpsaw mpfa

clear all
close all

%% Setup geometry

dims = [41, 20];
% dims = [2, 2];
G = cartGrid(dims, [2, 1]);
makeSkew = @(c) c(:,1) + .4*(1-(c(:,1)-1).^2).*(1-c(:,2));
G.nodes.coords(:,1) = 2*makeSkew(G.nodes.coords);
G.nodes.coords(:, 1) = G.nodes.coords(:, 1)*1000;
G.nodes.coords(:, 2) = G.nodes.coords(:, 2)*1000;

G = computeGeometry(G);

%% Homogeneous reservoir properties
alpha = 1; % biot's coefficient

rock = makeRock(G, 100*milli*darcy, .2);
rock.alpha = alpha*ones(G.cells.num, 1);
pv = sum(poreVolume(G, rock));

%% setup fluid and wells 

[f, info] = getCompositionalFluidCase('verysimple');
eos = EquationOfStateModel(G, f);

pRef = 50*barsa;
gravity reset off;

fluid = initSimpleADIFluid('cR', 1e-8/barsa, 'rho', [1, 1000, 100]);

% Symmetric well pattern
[ii, jj] = gridLogicalIndices(G);
% Injector + two producers
W = [];
W = addWell(W, G, rock, find(ii == ceil(G.cartDims(1)/2) & jj == G.cartDims(2)), 'comp_i', [1, 0], 'type', 'rate', 'val', pv/year);
W = addWell(W, G, rock, find(ii == G.cartDims(1) & jj == 1), 'comp_i', [1, 0], 'type', 'bhp', 'val', pRef);
W = addWell(W, G, rock, find(ii == 1 & jj == 1), 'comp_i', [1, 0], 'type', 'bhp', 'val', pRef);


for i = 1:numel(W)
    W(i).components = info.injection;
end
z0 = info.initial;
state0 = initCompositionalState(G, info.pressure, info.temp, [1, 0], z0, eos);
W(1).val = 100*W(1).val;

%% setup mechanics mech structure (with field prop and loadstruct)

mu = 1 * giga * Pascal;
lambda = 1e3*mu;

lambda = lambda*ones(G.cells.num, 1);
mu = mu*ones(G.cells.num, 1);
mechprop = struct('lambda', lambda, 'mu', mu);

[tbls, mappings] = setupStandardTables(G);

% We set zero displacement at all external faces

extfaces = find(any(G.faces.neighbors == 0, 2));
nextf    = numel(extfaces);
extfaces = rldecode(extfaces, 2*ones(nextf, 1));
linform  = repmat([[1, 0]; [0, 1]], nextf, 1);
bcvals   = zeros(numel(extfaces), 1);

bc = struct('linform'    , linform , ...
            'extfaces'   , extfaces, ...
            'linformvals', bcvals);

bc = setupFaceBC(bc, G, tbls);

nodefacecoltbl = tbls.nodefacecoltbl;
extforce = zeros(nodefacecoltbl.num, 1);

cellcoltbl = tbls.cellcoltbl;
force = zeros(cellcoltbl.num, 1);

loadstruct.bc = bc;
loadstruct.extforce = extforce;
loadstruct.force = force;

% setup mech structure 
mech.prop = mechprop;
mech.loadstruct = loadstruct;

modelMpfa = BiotCompositionalModel(G, rock, fluid, eos, mech, 'water', false);
modelMpfa.OutputStateFunctions = {'Dilatation', 'Stress'};
modelTpfa = BiotTpfaCompositionalModel(G, rock, fluid, eos, mech, 'water', false);
modelTpfa.OutputStateFunctions = {'Dilatation', 'Stress'};

mechmodel = MechModel(G, mech);
statemech = mechmodel.solveMechanics();

state0.u = statemech.u;
state0.lambdamech = statemech.lambdamech;
state0.biotpressure = state0.pressure;

dt = [1; 9; repmat(15, 22, 1)]*day;
schedule = simpleSchedule(dt, 'W', W);

[wsmpfa, statesMpfa] = simulateScheduleAD(state0, modelMpfa, schedule);
[wstpfa, statesTpfa] = simulateScheduleAD(state0, modelTpfa, schedule);

%% Reformat stress
for i = 1 : numel(statesMpfa)
    stress = modelMpfa.getProp(statesMpfa{i}, 'Stress');
    stress = formatField(stress, G.griddim, 'stress');
    statesMpfa{i}.stress = stress;
    stress = modelTpfa.getProp(statesTpfa{i}, 'Stress');
    stress = formatField(stress, G.griddim, 'stress');
    statesTpfa{i}.stress = stress;
end

%% plot

figure;
plotToolbar(G, statesMpfa);
title('MPFA')

figure
plotToolbar(G, statesTpfa);
title('TPFA')

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% The MPSA-W module is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% The MPSA-W module is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with the MPSA-W module.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>

