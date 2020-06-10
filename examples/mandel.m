
clear all
close all

%% Load required modules

mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui mpsaw mpfa

%% Setup grid

physdim = [40, 20] * meter;
resolution = [40, 20];
G = cartGrid(resolution, physdim);
G = computeGeometry(G);


% flow parameters
perm = 100 * milli * darcy;
muW = 0.89 * milli * Pascal / second;
poro = 0.25;

% elastic parameters
young = 1 * giga * Pascal;
poisson = 0.3;

E = young; nu = poisson; % shortcuts
lambda = E*nu/((1 + nu)*(1 - 2*nu));
mu = E/(2*(1 + nu));

alpha = 1; % biot's coefficient
% alpha = 0;

% for at top
top_force = 100 * mega * Pascal;

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
loadstruct.extforce = extforce;
loadstruct.force = force;

% setup mech structure 
mech.prop = mechprop;
mech.loadstruct = loadstruct;

%% Setup flow parameters (with field c and bcstruct)

fluid.c = cW; 

% setup boundary condition for flow


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
tsteps = 100;
% tsteps = 1;
duration = 10 * second;
schedule.step.val = duration/tsteps * ones(tsteps, 1);
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

u = states{end}.u;
u = reshape(u, 2, [])';

figure
plotCellData(G, u(:, 1));
title('displacenent in x-direction')
figure
plotCellData(G, u(:, 2));
title('displacenent in y-direction')

return

%% plot results from first row
figure
hold on
ind = (1 : resolution(1))';
xc = G.cells.centroids(ind, 1);
for i = 1 : numel(states)
    pc = states{i}.pressure(ind);
    plot(xc, pc);
end

%% plot value at middle 
figure
ind = floor(resolution(1)/2);
pmid = cellfun(@(state) state.pressure(ind), states);
