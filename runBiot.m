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


% setup mechanical driving forces (volumetric forces and boundary condition)
loadstruct = setupBCpercase(runcase, G, tbls, mappings);



% setup fluid driving forces (source terms and boundary condition)

extfaces = loadstruct.bc.extfaces; % get "external faces" from setupBCpercase

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
src = 0;

fluidforces = struct('bcstruct', bcstruct, ...
                     'sr', src);

% setup Biot's system driving forces
drivingforces = struct('mechanics', loadstruct, ...
                       'fluid', fluidforces);


% Setup mechanical parameter properties
% set mechanical material properties
nc = G.cells.num;
lambda = ones(nc, 1);
mu     = ones(nc, 1);

mechprops = struct('lambda', lambda, ...
                   'mu'    , mu);

% Setup fluid parameter properties
cellcolrowtbl = tbls.cellcolrowtbl;
colrowtbl = tbls.colrowtbl;

map = TensorMap();
map.fromTbl = colrowtbl;
map.toTbl = cellcolrowtbl;
map.mergefds = {'coldim', 'rowdim'};
map = map.setup();

K = [1; 0; 0; 1];
fluidprops = map.eval(K);

% setup Biot parameter properties

props = struct('mechprops' , mechprops, ...
               'fluidprops', fluidprops);

return

B   = assembly.B  ;
rhs = assembly.rhs;

sol = B\rhs;

% we compute the divergence
divop = assembly.divop;
divu = divop(sol);

% displacement values at cell centers.
cellcoltbl = tbls.cellcoltbl;
n = cellcoltbl.num;

u = sol(1 : n);

plotdeformedgrid = true;
if plotdeformedgrid
    lagmult = sol(n + 1 : end);
    unode = assembly.computeNodeDisp(u, lagmult);
    dim = G.griddim;
    unvec = reshape(unode, dim, [])';
end

dim = G.griddim;
uvec = reshape(u, dim, [])';

doplotsol = true;
if doplotsol
    figure
    plotCellData(G, uvec(:, 1));
    titlestr = sprintf('displacement - %s direction, eta=%g', 'x', eta);
    title(titlestr)
    colorbar
    figure
    plotCellData(G, uvec(:, 2));
    titlestr = sprintf('displacement - %s direction, eta=%g', 'y', eta);
    title(titlestr)
    colorbar
    figure
    plotCellData(G, divu);
    titlestr = sprintf('divergence, eta=%g', eta);
    title(titlestr)
    colorbar
end

if plotdeformedgrid
    figure 
    coef = 1e0;
    plotGridDeformed(G, coef*unvec);
end

doplotcontpoints = false;
if doplotcontpoints
    % plot continuity points
    [~, nodefacecents] = computeNodeFaceCentroids(G, eta, tbls, 'bcetazero', ...
                                                  bcetazero);
    figure
    hold on
    plotGrid(G)
    nodefacecents = reshape(nodefacecents, dim, [])';
    if dim == 2
        plot(nodefacecents(:, 1), nodefacecents(:, 2), '*');
    else
        plot3(nodefacecents(:, 1), nodefacecents(:, 2), nodefacecents(:, 3), '*');
    end

end

printcondest = false;
if printcondest

    matrices = assembly.matrices;

    A11 = matrices.A11;
    A12 = matrices.A12;
    A21 = matrices.A21;
    A22 = matrices.A22;
    D   = matrices.D  ;

    Z1 = zeros(size(A22, 1), size(D, 2));
    Z = zeros(size(D', 1), size(D, 2));

    A = [[A11, A12, -D];
         [A21, A22,  Z1];
         [D' , Z1',  Z]];

    fprintf('condest(A): %g\n', condest(A));
    fprintf('condest(A11): %g\n', condest(A11));
    fprintf('condest(A22): %g\n', condest(A22));
    fprintf('condest(B): %g\n', condest(B));
    
end

dplotinzdir = true;
if dplotinzdir
    if strcmp(runcase, '2d-linear')
        % plot displacement
        x = G.cells.centroids(:, 1);
        y = G.cells.centroids(:, 2);
        figure
        plot(x, uvec(:, 1), '*');
        titlestr = sprintf('displacement - %s direction, eta=%g', 'x', eta);
        title(titlestr)
        figure
        plot(y, uvec(:, 2), '*');
        titlestr = sprintf('displacement - %s direction, eta=%g', 'y', eta);
        title(titlestr)
    end
end

