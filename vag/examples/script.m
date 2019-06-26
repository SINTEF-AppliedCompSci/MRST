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

[Atrans, tbls] = computeVagTrans(G, rock);

cellnode2tbl = tbls.cellnode2tbl;

%% Setup system matrix

[A, op] = setupSystem(Atrans, cellnode2tbl, G);

rhsfun = op.rhsfun;
computeCellPressure = op.computeCellPressure;

%% setup source and boundary conditions

% We set up a source. Let us inject at a given rate in the first cell
cellinflux    = zeros(nc, 1);
cellinflux(1) = 100; % some rate
nodeinflux = zeros(nn, 1);
rhs = rhsfun(nodeinflux, cellinflux);

clear nodetbl;
nodetbl.nodes = (1 : nn)';
nodetbl.num = nn;
clear bcnodetbl;
bcnodetbl.nodes = nn;
bcnodetbl.num = 1;

intnodes = true(nn, 1);
intnodes(bcnodetbl.nodes) = false;
intnodes = find(intnodes);
clear intnodetbl
intnodetbl.nodes = intnodes;
intnodetbl.num = numel(intnodes);

A11 = A(intnodetbl.nodes, intnodetbl.nodes);
A12 = A(intnodetbl.nodes, bcnodetbl.nodes);

pval = 100;
pbc = pval*ones(bcnodetbl.num, 1);

rhs = rhs(intnodetbl.nodes);
rhs = rhs - A12*pbc;


%% Solve system
x = A11\rhs;

%% Recover node  pressures
pn = NaN(nn, 1);
pn(intnodetbl.nodes) = x;
pn(bcnodetbl.nodes)  = pbc;

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