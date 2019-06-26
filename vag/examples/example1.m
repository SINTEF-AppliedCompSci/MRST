%% Use VAG transmissibilities and hybridization to compute solution

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

A       = vagstruct.A;
cellnode2tbl = vagstruct.cellnode2tbl;
cellnodetbl  = vagstruct.cellnodetbl;

cellnodetbl = addLocInd(cellnodetbl, 'cnind');

[~, cellnode2tbl] = setupTableMapping(cellnode2tbl, cellnodetbl, {'cells', ...
                    {'nodes1', 'nodes'}});
cellnode2tbl = replacefield(cellnode2tbl, {'cnind', 'cnind1'});
[~, cellnode2tbl] = setupTableMapping(cellnode2tbl, cellnodetbl, {'cells', ...
                    {'nodes2', 'nodes'}});
cellnode2tbl = replacefield(cellnode2tbl, {'cnind', 'cnind2'});

tbl = cellnode2tbl; %alias
A = sparse(tbl.cnind1, tbl.cnind2, A, cellnodetbl.num, cellnodetbl.num);

nnc = cellnodetbl.num;
R1 = sparse(cellnodetbl.cnind, cellnodetbl.cells, ones(nnc, 1), nnc, nc);
R2 = sparse(cellnodetbl.cnind, cellnodetbl.nodes, ones(nnc, 1), nnc, nn);

% We defined the matrix opcells with coefficients either
% equal to zero and one We have opcell(i, j) = 1 only if the
% node i (in node-localCell indexing) and the node j (in
% cell-localNode indexing) has a cell in common.
L1 = R1';
L2 = sparse(cellnodetbl.cells, cellnodetbl.cnind, ones(nnc, 1), nc, nnc);
opcells = L2'*L1;

% We defined the matrix opnodes. Same as opcells but for nodes.
L1 = sparse(cellnodetbl.nodes, cellnodetbl.cnind, ones(nnc, 1), nn, nnc);
L2 = R2';
opnodes = L2'*L1;

% We take intersection of the two: It means that S(i, j) = 1
% if node i and j (in their respective indexing) corresponds
% to the same "physical" cell and node.
S = opcells.*opnodes;

% We assemble the helper matrices, as described in slides
M11 = R1'*A*R1;
M12 = -R1'*A*S'*R2;
M21 = M12';
M22 = R2'*S*A*S'*R2;

M = [[M11, M12]; [M21, M22]];
n = nc + nn;

% We use Lagrangian multipliers to enforce the boundary
% conditions for pressure
%ind = nc + 1; % we impose pressure at first node
ind = n  ; % we set pressure at last node
nind = size(ind, 1);
dirop = sparse((1 : nind)', ind, ones(nind, 1), nind, n);

M = [[M, -dirop']; [dirop, zeros(nind)]];

% We can set up the source.
% Let us choose the first cell 
f = zeros(n, 1);
f(1) = 1000; % some rate 
%pval = zeros(nind, 1); % we impose a non zero pressure;
pval=100;


rhs = [f; pval];

x = M\rhs;
u = x(1 : nc);
v = x((nc + 1)' : n);
bhp= x(n+1:n+1);
figure(1)
clf
plotCellData(G, u)
figure(2)
clf
plotNodeData(G, v)
