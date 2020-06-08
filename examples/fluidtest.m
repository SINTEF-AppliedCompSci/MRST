
clear all
close all

%% load modules

mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui mpfa mpsaw

%% Setup grid

dimcase = 2;
switch dimcase
  case 2
    nx = 10; ny = 10;
    G = cartGrid([nx, ny]);
  case 3
    nx = 4; ny = 3; nz = 5;
    G = cartGrid([nx, ny, nz]);
end
G = twister(G, 0.1); 
G = computeGeometry(G);
nc = G.cells.num;

%% Setup rock parameters (for flow)

rock.perm = darcy*ones(G.cells.num, 1);
% rock.perm = ones(G.cells.num, 1);
rock.poro = 0.3*ones(G.cells.num, 1);
% this is the Biot parameter which will be used to multiply the given dilatation
rock.alpha = ones(G.cells.num, 1);

%% Setup fluid parameters from SPE1

fluid.c = 1e-8;

%% Gravity
% The gravity in this option affects only the fluid behavior
gravity off;

ptop    = 2*barsa;
pbottom = 1*barsa;

%% setup bc for pressure
dim = G.griddim;
zface = G.faces.centroids(:, dim);
maxz = max(zface);
minz = min(zface);
bottomfaces = find(G.faces.centroids(:, dim) > (maxz - 1e-4));
topfaces = find(G.faces.centroids(:, dim) < (minz + 1e-4));

nbf = numel(bottomfaces);
ntf = numel(topfaces);
bcfaces = [bottomfaces; topfaces];
bcvals = ptop*ones(nbf + ntf, 1);
bcvals(1 : nbf) = pbottom;

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


%% construct model

fluidmodel = BiotFluidModel(G, rock, fluid, bcstruct);

%% Setup schedule
schedule.step.val = 1e-1*day*ones(10, 1);
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control = struct('W', []);

%% Setup initial state
clear initState;
initState.pressure = pbottom*ones(G.cells.num, 1);
initState.lambda = zeros(bcnodefacetbl.num, 1);

solver = NonLinearSolver('maxIterations', 100);
[wsol, states] = simulateScheduleAD(initState, fluidmodel, schedule, 'nonlinearsolver', solver);

%% plot results
plotToolbar(G, states);
