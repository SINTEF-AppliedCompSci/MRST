clear all

mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui mpsaw mpfa


%% Define and process geometry
% Construct a Cartesian grid 
runcases = {'2d-refinement', ...
            '2d-linear'    , ...
            '2d-compaction', ...
            '3d-linear'    , ...
            '3d-compaction' };

runcase = '2d-compaction';

switch runcase
  case '2d-refinement'
    ny = 4;
    dx = 1e-3;
    dy = [dx; ones(ny, 1)];
    y = [0; cumsum(dy)];
    y = 1/max(y)*y;
    dx = [dx; ones(ny, 1); dx];
    x = [0; cumsum(dx)];
    x = 1/max(x)*x;
    G = tensorGrid(x, y);    
  case {'2d-linear', '2d-compaction'}
    N = 40;
    nx = N; ny = N;
    G = cartGrid([nx, ny], [1, 1]);
    % N = 3;
    % Nx = N*ones(1, 2);
    % G = createBisectedTriangleGrid(Nx,1);
    % G = twister(G, 0.1);
  case {'3d-linear', '3d-compaction'}
    nx = 5;
    ny = nx;
    nz = nx;
    G = cartGrid([nx, ny, nz], [1, 1, 1]);
  otherwise
    error('runcase not recognized');
end

% G = twister(G, 0.1);
% compute Grid geometry
G = computeGeometry(G);
dim = G.griddim;

% set material properties
nc = G.cells.num;

lambda = ones(nc, 1);
mu     = ones(nc, 1);

prop = struct('lambda', lambda, ...
              'mu', mu);


[tbls, mappings] = setupStandardTables(G, 'useVirtual', false);
loadstruct = setupBCpercase(runcase, G, tbls, mappings);

mech.prop = prop;
mech.loadstruct = loadstruct;

%% Setup model
model = MechModel(G, mech);

state = model.solveMechanics();

u = model.getProp(state, 'displacement');
u = formatField(u, dim, 'u');

unf = model.getProp(state, 'FaceNodeDisplacement');

cellnodefacecoltbl = tbls.cellnodefacecoltbl;
nodefacecoltbl = tbls.nodefacecoltbl;
cellcoltbl = tbls.cellcoltbl;

map = TensorMap();
map.fromTbl = nodefacecoltbl;
map.toTbl = cellnodefacecoltbl;
map.mergefds = {'nodes', 'faces', 'coldim'};
map = map.setup();
ucnf = map.eval(unf);

map = TensorMap();
map.fromTbl = cellnodefacecoltbl;
map.toTbl = cellcoltbl;
map.mergefds = {'cells', 'coldim'};
map = map.setup();
unf = map.eval(ucnf);

unf = formatField(unf, dim, 'u');

stress = model.getProp(state, 'Stress');
stress = formatField(stress, dim, 'stress');

figure(1)
clf
plotCellData(G, u(:, 1))
title('x-displacement')

figure(2)
clf
plotCellData(G, u(:, 2))
title('y-displacement')

figure(3)
clf
plotCellData(G, stress(:, 1))
title('xx-stress')

figure(4)
clf
plotCellData(G, stress(:, 2))
title('yy-stress')

figure(5)
clf
plotCellData(G, stress(:, 3))
title('xy-stress')

figure
clf
plotCellData(G, unf(:, 1))
title('face node x-displacement')

figure
clf
plotCellData(G, unf(:, 2))
title('face node y-displacement')
