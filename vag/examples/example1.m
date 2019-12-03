%% Use VAG transmissibilities and hybridization to compute solution

%% Load modules

mrstModule add vag vem

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

%% Setup system matrixrhsfun = op.rhsfun;

[A, op] = setupSystem(vagstruct, G);

rhsfun = op.rhsfun;
computeCellPressure = op.computeCellPressure;

%% Setup boundary condition: Given influx at first node and given pressure at
%% last node

% We inject in the first node, no injection elsewhere

cellinflux    = zeros(nc, 1);
nodeinflux    = zeros(nn, 1);
nodeinflux(1) = 1e8*meter^3/day;
rhs = rhsfun(nodeinflux, cellinflux);

% We impose a given pressure p in the last node

p = 100*barsa;
intnodes = 1 : (nn - 1);
A11 = A(intnodes, intnodes);
A12 = A(intnodes, nn);

rhs = rhs(intnodes) - A12*p;

%% Solve system

x = A11\rhs;

%% Recover node  pressures

pn = NaN(nn, 1);
pn(intnodes) = x;
pn(nn)       = p;


%% Recover cell  pressures

pc = computeCellPressure(pn, cellinflux);


%% plotting

figure(1)
clf
plotCellData(G, pc);
colorbar

figure(2)
clf
plotNodeData(G, pn);
colorbar