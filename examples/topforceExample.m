%% Example where we impose a controlled force at the top
% The intensity of the force is gradually increased (see setup of schedule
% below).
clear all
close all

%% Load required modules
mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui mpsaw mpfa

%% Setup grid
physdim = [1, 1] * meter;
nx = 40; ny = 40;
resolution = [nx, ny];
G = cartGrid(resolution, physdim);
G = computeGeometry(G);

% flow parameters
perm = 1;
muW  = 1;
poro = 1;

% Young's modulus
E = 1;
% Poisson's ratio
nu = 0;

% First Lamé parameter
lambda = E*nu/((1 + nu)*(1 - 2*nu));
% Second Lamé parameter (also called Shear modulus and denoted G)
mu = E/(2*(1 + nu));
Gm = mu;
% Bulk modulus K
K = E/(2*(1 - 2*nu)); 

% Biot's coefficient
alpha = 1; 

% force at the top (maximal value)
top_force = 1;

rock.poro  = poro * ones(G.cells.num, 1);
rock.perm  = (perm/muW) * ones(G.cells.num, 1);
rock.alpha = alpha * ones(G.cells.num, 1);

% fluid compressibility
cW = 0; % incompressible

% Reference pressure
pref = 0*barsa;

%% setup mechanics mech structure (with field prop and loadstruct)
lambda = lambda*ones(G.cells.num, 1);
mu = mu*ones(G.cells.num, 1);
mechprop = struct('lambda', lambda, 'mu', mu);

[tbls, mappings] = setupStandardTables(G);

% We recover the top, bottom and lateral faces using function pside
dummy = 0;
bc = pside([], G, 'Ymax', dummy); 
topfaces = bc.face;
bc = pside([], G, 'Xmin', dummy); 
leftfaces = bc.face;
bc = pside([], G, 'Xmax', dummy); 
rightfaces = bc.face;
bc = pside([], G, 'Ymin', dummy); 
bottomfaces = bc.face;

lateralfaces = leftfaces;

%% Setup the Dirichlet boundary condition for the mechanics

% on the bottom, we have rolling condition in x-direction
bottomlinform = repmat([0, 1], numel(bottomfaces), 1);
% on the lateral walls, we have rolling condition in y-direction
laterallinform = repmat([1, 0], numel(lateralfaces), 1);

linform  = [bottomlinform; laterallinform];
extfaces = [bottomfaces; lateralfaces];
bcvals   = zeros(numel(extfaces), 1);

bc = struct('linform'    , linform , ...
            'extfaces'   , extfaces, ...
            'linformvals', bcvals);

bc = setupFaceBC(bc, G, tbls);

%% setup vertical force at the top

topfacetbl.faces = topfaces;
topfacetbl = IndexArray(topfacetbl);
nodefacecoltbl = tbls.nodefacecoltbl;
topnodefacecoltbl = crossIndexArray(topfacetbl, nodefacecoltbl, {'faces'});

cellnodefacetbl = tbls.cellnodefacetbl;
cellnodefacecoltbl = tbls.cellnodefacecoltbl;
facetNormals = computeFacetNormals(G, cellnodefacetbl);

map = TensorMap();
map.fromTbl = cellnodefacecoltbl;
map.toTbl = topnodefacecoltbl;
map.mergefds = {'faces', 'nodes', 'coldim'};
map = map.setup();

topnormals = map.eval(facetNormals);
extforce = -top_force*topnormals;

map = TensorMap();
map.fromTbl = topnodefacecoltbl;
map.toTbl = nodefacecoltbl;
map.mergefds = {'nodes', 'faces', 'coldim'};
map = map.setup();

extforce = map.eval(extforce);

cellcoltbl = tbls.cellcoltbl;
force = zeros(cellcoltbl.num, 1);

loadstruct.bc = bc;
loadstruct.extforce = zeros(nodefacecoltbl.num, 1); % get the dimension right
loadstruct.force = force;

% Setup the mechanical structure (input to to setup model)
mech.prop = mechprop;
mech.loadstruct = loadstruct;

%% Setup flow parameters (with field c and bcstruct)

fluid.c = cW;
fluid.src = []; % no source

%% Setup boundary conditions for flow
% We impose a constant pressure on the faces on the right handside.

bcfaces = rightfaces;
bcvals = pref*ones(numel(bcfaces));

nodefacetbl = tbls.nodefacetbl;

bcfacetbl.faces = bcfaces;
bcfacetbl = IndexArray(bcfacetbl);
bcnodefacetbl = crossIndexArray(nodefacetbl, bcfacetbl, {'faces'});

map = TensorMap();
map.fromTbl = bcfacetbl;
map.toTbl = bcnodefacetbl;
map.mergefds = {'faces'};
map = map.setup();

bcvals = map.eval(bcvals);

bcdirichlet = struct('bcnodefacetbl', bcnodefacetbl, ...
                     'bcvals', bcvals);
bcneumann = [];

bcstruct = struct('bcdirichlet', bcdirichlet, ...
                  'bcneumann'  , bcneumann);

fluid.bcstruct = bcstruct;

%% Setup Biot model

model = BiotModel(G, rock, fluid, mech);
model = model.validateModel();
model.OutputStateFunctions = {'Dilatation', 'Stress'};

%% Setup schedule
% We gradually increase the force exterted at the top
tsteps1 = 10; 
duration1 = 1;
tsteps2 = 10;
duration2 = 4;
val1 = (duration1/tsteps1)*ones(tsteps1, 1);
val2 = linspace(0, 1, tsteps2 + 1)'.^4;
val2 = duration2*diff(val2);

schedule.step.val = [val1; val2];
schedule.step.control = [(1 : tsteps1)'; (tsteps1 + 1)*ones(tsteps2, 1)];
for i = 1 : tsteps1
    coef = i/tsteps1;
    control(i) = struct('W', [], 'extforce', coef*extforce);
end
control(tsteps1 + 1) = struct('W', [], 'extforce', extforce);
schedule.control = control;

%% Setup initial state
clear initState;
% fluid
initState.pressure = zeros(G.cells.num, 1);
nlf = size(bcstruct.bcdirichlet.bcvals, 1);
initState.lambdafluid = zeros(nlf, 1);
% mech
cellcoltbl = tbls.cellcoltbl;
initState.u = zeros(cellcoltbl.num, 1);
nlm = size(loadstruct.bc.linformvals, 1);
initState.lambdamech = zeros(nlm, 1);
initState.extforce = 0*extforce;

[wsol, states] = simulateScheduleAD(initState, model, schedule);

%% Pressure for some selected times
figure
clf
xind = (1 : nx)';
xc = G.cells.centroids(xind, 1);
hold on

inds = 1 : floor(tsteps1/5) : tsteps1; 
legends = {};
tt = cumsum(schedule.step.val);
for i = 1 : numel(inds)
    ind = inds(i);
    p = states{ind}.pressure;
    p = p(xind);
    plot(xc, p);
    legends{end + 1} = sprintf('time %g s', tt(ind));
end
legend(legends{:});
title('Pressure profile at some selected times');

%% Reformat stress
for i = 1 : numel(states)
    stress = model.getProp(states{i}, 'Stress');
    stress = formatField(stress, G.griddim, 'stress');
    states{i}.stress = stress;
end 

%% plotting
figure 
clf
plotToolbar(G, states);

figure
clf
ind = 1;
pmid = cellfun(@(state) state.pressure(ind), states);
tt = cumsum(schedule.step.val);
plot(tt, pmid, '*');
title('Pressure evolution at left edge')

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

