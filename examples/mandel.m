
clear all
close all

%% Load required modules

mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui

%% Setup grid

nx = 10; ny = 10;
G = cartGrid([nx, ny], [1, 1]);
G = computeGeometry(G);
nc = G.cells.num;

%% setup mechanics mech structure (with field prop and loadstruct)

lambda = 0*ones(nc, 1);
mu     = ones(nc, 1);

mechprop = struct('lambda', lambda, ...
                  'mu', mu);

[tbls, mappings] = setupStandardTables(G);

% We recover the top, bottom and lateral faces using function pside
dummy = 0;
bc = pside([], G, 'Ymax', dummy); 
topfaces = bc.face;
bc = pside([], G, 'Xmin', dummy); 
bc = pside(bc, G, 'Xmax', dummy); 
lateralfaces = bc.face;
bc = pside([], G, 'Ymin', dummy); 
bottomfaces = bc.face;

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
extforce = -topnormals;

map = TensorMap();
map.fromTbl = topnodefacecoltbl;
map.toTbl = nodefacecoltbl;
map.mergefds = {'nodes', 'faces', 'coldim'};
map = map.setup();

extforce = map.eval(extforce);

cellcoltbl = tbls.cellcoltbl;
force = zeros(cellcoltbl.num, 1);

loadstruct.bc = bc;
loadstruct.extforce = extforce;
loadstruct.force = force;

% setup mech structure 
mech.prop = mechprop;
mech.loadstruct = loadstruct;

%% Setup rock parameters (for flow)

% rock.perm = 100*milli*darcy*ones(G.cells.num, 1);
rock.perm = 1*ones(G.cells.num, 1);
% rock.poro = 0.25*ones(G.cells.num, 1);
rock.poro = 1*ones(G.cells.num, 1);
rock.alpha = 1e-1*ones(G.cells.num, 1);

%% Setup flow parameters (with field c and bcstruct)

fluid.c = 0; % incompressible

% setup boundary condition for flow

pref = 0*barsa;

bcfaces = lateralfaces;
bcvals = pref*ones(numel(bcfaces));

useVirtual = false;
[tbls, mappings] = setupStandardTables(G, 'useVirtual', useVirtual);
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

model =  BiotModel(G, rock, fluid, mech);

%% Setup schedule
tottime = 1;
N = 30;
alpha = 20;

dt = 1/N;
t = [0 : dt : 1];
t = tottime/(exp(alpha) - 1)*(exp(alpha*t)- 1);

schedule.step.val = diff(t);
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control = struct('W', []);

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

solver = NonLinearSolver('maxIterations', 100);
[wsol, states] = simulateScheduleAD(initState, model, schedule, 'nonlinearsolver', solver);

%% plot results
plotToolbar(G, states);
colorbar
