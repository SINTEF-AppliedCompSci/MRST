
clear all
close all


%% Load required modules

mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui

%% Setup grid

dimcase = 2;
switch dimcase
  case 2
    nx = 40; ny = 40;
    G = cartGrid([nx, ny]);
  case 3
    nx = 4; ny = 3; nz = 5;
    G = cartGrid([nx, ny, nz]);
end

G = twister(G, 0.1); 
G = computeGeometry(G);
nc = G.cells.num;

%% setup mechanics mech structure (with field prop and loadstruct)


lambda = ones(nc, 1);
mu     = ones(nc, 1);

mechprop = struct('lambda', lambda, ...
                  'mu', mu);

[tbls, mappings] = setupStandardTables(G);
runcase = '2d-compaction';
loadstruct = setupBCpercase(runcase, G, tbls, mappings);
bc = setupFaceBC(loadstruct.bc, G, tbls);
loadstruct.bc = bc;

mech.prop = mechprop;
mech.loadstruct = loadstruct;

%% Setup rock parameters (for flow)

rock.perm = ones(G.cells.num, 1);
% rock.perm = ones(G.cells.num, 1);
rock.poro = 0.3*ones(G.cells.num, 1);
% this is the Biot parameter which will be used to multiply the given dilatation
rock.alpha = ones(G.cells.num, 1);

%% Setup flow parameters (with field c and bcstruct)

fluid.c = 1;

%% setup bc for pressure

pbottom = 0*barsa;

dim = G.griddim;
zface = G.faces.centroids(:, dim);
maxz = max(zface);
minz = min(zface);
bottomfaces = find(G.faces.centroids(:, dim) < (minz + 1e-4));
topfaces = find(G.faces.centroids(:, dim) > (maxz - 1e-4));

nbf = numel(bottomfaces);
bcfaces = bottomfaces; 
bcvals = pbottom*ones(nbf, 1);

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

schedule.step.val = 1*ones(40, 1);
schedule.step.control = ones(numel(schedule.step.val), 1);
extforce = loadstruct.extforce;
schedule.control = struct('W', [], 'extforce', extforce);

%% Setup initial state
clear initState;
% fluid
initState.pressure = pbottom*ones(G.cells.num, 1);
nlf = size(bcstruct.bcdirichlet.bcvals, 1);
initState.lambdafluid = zeros(nlf, 1);
% mech
cellcoltbl = tbls.cellcoltbl;
initState.u = zeros(cellcoltbl.num, 1);
nlm = size(loadstruct.bc.linformvals, 1);
initState.lambdamech = zeros(nlm, 1);
initState.extforce = extforce;

solver = NonLinearSolver('maxIterations', 100);
[wsol, states] = simulateScheduleAD(initState, model, schedule, 'nonlinearsolver', solver);

%% plot results
plotToolbar(G, states);