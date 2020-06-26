%% Run Biot example
%%
%% Reference paper:
%% Finite volume methods for elasticity with weak symmetry
%% Keilegavlen, Eirik and Nordbotten, Jan Martin
%% International Journal for Numerical Methods in Engineering
%% 2017

clear all
close all

tic

% load modules
mrstModule add mimetic mpsaw incomp vemmech mpfa

eta = 0;

bcetazero = false;

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
    N = 20;
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

useVirtual = false;
[tbls, mappings] = setupStandardTables(G, 'useVirtual', useVirtual);


% Setup mechanical driving forces (volumetric forces and boundary condition)
loadstruct = setupBCpercase(runcase, G, tbls, mappings);

% Setup fluid driving forces (source terms and boundary condition)

extfaces = loadstruct.bc.extfaces; % get "external faces" from setupBCpercase
extfaces = unique(extfaces);

bcfacetbl.faces = extfaces;
bcfacetbl = IndexArray(bcfacetbl);
nodefacetbl = tbls.nodefacetbl;
bcnodefacetbl = crossIndexArray(bcfacetbl, nodefacetbl, {'faces'});
bcvals = zeros(bcnodefacetbl.num, 1);
bcneumann = [];
bcdirichlet = struct('bcnodefacetbl', bcnodefacetbl, ...
                     'bcvals', bcvals);
bcstruct = struct('bcdirichlet', bcdirichlet, ...
                  'bcneumann'  , bcneumann);

nc = G.cells.num;
src = zeros(nc, 1);
src(end) = 1;

fluidforces = struct('bcstruct', bcstruct, ...
                     'src', src);


% Setup fluid parameter properties
cellcolrowtbl = tbls.cellcolrowtbl;
colrowtbl = tbls.colrowtbl;

map = TensorMap();
map.fromTbl = colrowtbl;
map.toTbl = cellcolrowtbl;
map.mergefds = {'coldim', 'rowdim'};
map = map.setup();

K = [1; 0; 0; 1];
K = map.eval(K);
fluidprops.K = K;

eta = 1/3;
assembly = assembleMPFA(G, K, bcstruct, src, eta, tbls, mappings);

B = assembly.B;
rhs = assembly.rhs;

sol = B\rhs;

nc = G.cells.num;
p = sol(1 : nc);

plotCellData(G, p);