clear all
close all

%% Setup geometry

dims = [30, 30];
% dims = [2, 2];
G = cartGrid(dims, [2, 1]);
makeSkew = @(c) c(:,1) + .4*(1-(c(:,1)-1).^2).*(1-c(:,2));
G.nodes.coords(:,1) = 2*makeSkew(G.nodes.coords);
G.nodes.coords(:, 1) = G.nodes.coords(:, 1);
G.nodes.coords(:, 2) = G.nodes.coords(:, 2);

G = computeGeometry(G);

%% Homogeneous reservoir properties
alpha = 1; % biot's coefficient

rock = makeRock(G, 1, 1);
rock.alpha = alpha*ones(G.cells.num, 1);
pv = sum(poreVolume(G, rock));

%% setup fluid and wells 

pRef = 0*barsa;
gravity reset off;
fluid = initSimpleADIFluid('cR', 1, 'pRef', pRef, 'mu', [1, 1, 1]);

% Symmetric well pattern
[ii, jj] = gridLogicalIndices(G);

% Injector + two producers
W = [];
W = addWell(W, G, rock, find(ii == ceil(G.cartDims(1)/2) & jj == G.cartDims(2)), 'radius', 1e-3, 'comp_i', [1, 0], 'type', 'rate', 'val', pv);
W = addWell(W, G, rock, find(ii == G.cartDims(1) & jj == 1), 'radius', 1e-3, 'comp_i', [1, 0], 'type', 'bhp', 'val', pRef);
W = addWell(W, G, rock, find(ii == 1 & jj == 1), 'radius', 1e-3, 'comp_i', [1, 0], 'type', 'bhp', 'val', pRef);

%% setup mechanics mech structure (with field prop and loadstruct)

lambda    = 0;
mu        = 1;

lambda = lambda*ones(G.cells.num, 1);
mu = mu*ones(G.cells.num, 1);
mechprop = struct('lambda', lambda, 'mu', mu);

useVirtual = true;
[tbls, mappings] = setupMpxaStandardTables(G, 'useVirtual', useVirtual);

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

nodefacevectbl = tbls.nodefacevectbl;
extforce = zeros(nodefacevectbl.num, 1);

cellvectbl = tbls.cellvectbl;
force = zeros(cellvectbl.num, 1);

loadstruct.bc = bc;
loadstruct.extforce = extforce;
loadstruct.force = force;

% setup mech structure 
mech.prop = mechprop;
mech.loadstruct = loadstruct;

modelMpfa = BiotBlackOilModel(G, rock, fluid, mech, 'water', true, 'oil', true, 'gas', false);
modelMpfa.OutputStateFunctions = {'Dilatation', 'Stress'};
modelTpfa = BiotTpfaBlackOilModel(G, rock, fluid, mech, 'water', true, 'oil', true, 'gas', false);
modelTpfa.OutputStateFunctions = {'Dilatation', 'Stress'};

mechmodel = MechModel(G, mech);
statemech = mechmodel.solveMechanics();

state0 = initResSol(G, pRef, [0, 1]);
state0.u = statemech.u;
state0.lambdamech = statemech.lambdamech;
state0.biotpressure = state0.pressure;

dt = 1e-1;
nsteps = 10;
schedule.step.val = dt*ones(nsteps, 1);
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control = struct('W', W);

[ws, statesMpfa] = simulateScheduleAD(state0, modelMpfa, schedule);
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

