%% Mandel-like problem

clear all
% close all

%% Load required modules

mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui mpsaw mpfa

%% Setup grid

physdim = [1, 1] * meter;
nx = 40;, ny = 20;
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

% consolidation coefficient
cv = perm/muW*(K + 3*Gm); % see reference verruijt2013theory

alpha = 1; % biot's coefficient

% for at top
% top_force = 100 * mega * Pascal;
top_force = 1;

rock.poro  = poro * ones(G.cells.num, 1);
rock.perm  = (perm/muW) * ones(G.cells.num, 1);
rock.alpha = alpha * ones(G.cells.num, 1);

% incompressible fluid
cW = 0; 

% reference pressure on the side
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

% setup vertical force at the top

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

% setup mech structure 
mech.prop = mechprop;
mech.loadstruct = loadstruct;

%% Setup flow parameters (with field c and bcstruct)

fluid.c = cW;
fluid.src = [];

% Setup boundary conditions for flow

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

%% Setup schedule
tsteps = 100;
duration = 1;
schedule.step.val = duration/tsteps * ones(tsteps, 1);
schedule.step.control = (1 : tsteps)';
for i = 1 : tsteps
    coef = i/tsteps;
    control(i) = struct('W', [], 'extforce', coef*extforce);
end
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

solver = NonLinearSolver('maxIterations', 100);
[wsol, states] = simulateScheduleAD(initState, model, schedule, 'nonlinearsolver', solver);

%% Pressure for some selected times
figure
clf
xind = (1 : nx)';
xc = G.cells.centroids(xind, 1);
hold on

inds = 1 : floor(tsteps/5) : tsteps; 
legends = {};
tt = cumsum(schedule.step.val);
for i = 1 : numel(inds);
    ind = inds(i);
    p = states{ind}.pressure;
    p = p(xind);
    plot(xc, p);
    legends{end + 1} = sprintf('time %g s', tt(ind));
end
legend(legends{:});
title('Pressure profile at some selected times');

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

