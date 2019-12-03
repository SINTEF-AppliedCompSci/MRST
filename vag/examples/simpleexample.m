%% Load modules

mrstModule add mpfa vem vag vemmech

%% Setup a Cartesian grid

nx = 10;
ny = 10;
nz = 1;
cartDims = [nx, ny, nz];
physDims = [nx*100, ny*50, 20];

G = cartGrid(cartDims, physDims);
G = computeGeometry(G);

nc = G.cells.num;
nn = G.nodes.num;

perm = 0.008527*100;
rock = makeRock(G, perm*ones(nc, 1), 0.1*ones(nc, 1));

%% Compute the VAG transmissibilities

vagstruct = computeVagTrans(G, rock);

%% Setup system matrix

[A, op] = setupSystem(vagstruct, G);

rhsfun = op.rhsfun;
computeCellPressure = op.computeCellPressure;

%% setup source and boundary conditions

% We set up a source. Let us inject at a given rate in the first cell
cellinflux    = zeros(nc, 1);
cellinflux(1) = 100; % some rate
nodeinflux = zeros(nn, 1);
rhs = rhsfun(nodeinflux, cellinflux);


% We define "boundary nodes" where we impose a given pressure. In this case,
% we consider only the last node

bcnodes = nn;
nbcn = numel(bcnodes);

intnodes = true(nn, 1);
intnodes(bcnodes) = false;
intnodes = find(intnodes);

A11 = A(intnodes, intnodes);
A12 = A(intnodes, bcnodes);

% We impose a given pressure at the nodes
pval = 100;
pbc = pval*ones(nbcn, 1);

rhs = rhs(intnodes);
rhs = rhs - A12*pbc;


%% Solve system

x = A11\rhs;

%% Recover nodal  pressures

pn = NaN(nn, 1);
pn(intnodes) = x;
pn(bcnodes)  = pbc;

%% Recover cell  pressures

pc = computeCellPressure(pn, cellinflux);

%% Plot solution

close all

figure
plotNodeData(G, pn)
title('nodal pressure')
axis equal
colorbar

figure
plotCellData(G, pc)
title('cell pressure')
axis equal
colorbar