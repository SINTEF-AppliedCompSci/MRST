%% Assembly of MPSA-weak
%%
%% Reference paper:
%% Finite volume methods for elasticity with weak symmetry
%% Keilegavlen, Eirik and Nordbotten, Jan Martin
%% International Journal for Numerical Methods in Engineering
%% 2017

clear all

tic

% load modules
mrstModule add mimetic mpsaw incomp vemmech

eta = 1/3;

%% Define and process geometry
% Construct a Cartesian grid 
runcases = {'2d-refinement', ...
            '2d-linear'    , ...
            '2d-compaction', ...
            '3d-linear'    , ...
            '3d-compaction' };

runcase = '3d-compaction';

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
    nx = 5; ny = 5;
    G = cartGrid([nx, ny], [1, 1]);
  case {'3d-linear', '3d-compaction'}
    nx = 5; ny = 5; nz = 5;
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

[tbls, mappings] = setupStandardTables(G);
loadstruct = setupBCpercase(runcase, G, tbls, mappings);
assembly = assembleMPSA(G, prop, loadstruct, eta, tbls, mappings);

B   = assembly.B  ;
rhs = assembly.rhs;

sol = B\rhs;

% displacement values at cell centers.
cellcoltbl = tbls.cellcoltbl;
n = cellcoltbl.num;

u = sol(1 : n);

% Force where the Dirichlet BC are imposed (the are given by the lagrange
% multipliers)
lagmult = sol(n + 1: end);

% Compute displacement at nodes

force = assembly.extforce;
matrices = assembly.matrices;
nodaldisp_op = assembly.nodaldisp_op;

invA11 = matrices.invA11;
D      = matrices.D;
A12    = matrices.A12;

% displacement values at facenode
unf = invA11*(force - A12*u + D*lagmult);
un = nodaldisp_op*unf;

%% plotting
% 
close all

unvec = reshape(un, dim, [])';
figure 
coef = 1e0;
plotGridDeformed(G, coef*unvec);

%%
dim = G.griddim;
uvec = reshape(u, dim, [])';

figure
plotCellData(G, uvec(:, 1));
title('displacement - x direction')
colorbar
figure
plotCellData(G, uvec(:, 2));
title('displacement - y direction')
colorbar