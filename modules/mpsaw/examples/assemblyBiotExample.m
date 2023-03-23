%% Run Biot example
%
% Simple examples

clear all
close all

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
    N = 15;
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

% get some "external faces" from setupBCpercase
bcnodefacetbl = loadstruct.bc.bcnodefacetbl; 
extfaces = bcnodefacetbl.get('faces');
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
src = [];

fluidforces = struct('bcstruct', bcstruct, ...
                     'src', src);

% Setup Biot's system driving forces
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
K = map.eval(K);
fluidprops.K = K;

% setup Biot parameter properties

alpha = 1;
rho = 0;

alpha = alpha*ones(nc, 1);
rho = rho*ones(nc, 1);
coupprops = struct('alpha', alpha, ...
                   'rho', rho);

props = struct('mechprops' , mechprops , ...
               'fluidprops', fluidprops, ...
               'coupprops' , coupprops);

casetype = 'blockassembly';
switch casetype
  case 'blockassembly'
    assembly = blockAssembleBiot(G, props, drivingforces, eta, tbls, mappings, 'blocksize', 20, 'verbose', true);
  case 'standard'
    assembly = assembleBiot(G, props, drivingforces, eta, tbls, mappings);
end

B = assembly.B;
rhs = assembly.rhs;
sol = B\rhs;

ncc  = tbls.cellcoltbl.num;
nc   = tbls.celltbl.num;
n1 = 1; n2 = ncc;
ucc  = sol(n1 : n2);
n1 = n2 + 1; n2 = n1 + nc - 1;
pc = sol(n1 : n2);

doplot = true;
if doplot
    dim = G.griddim;
    uvec = reshape(ucc, dim, [])';
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
    plotCellData(G, pc);
    titlestr = sprintf('pressure');
    title(titlestr)
    colorbar
    
end

plotdeformedgrid = true;
if plotdeformedgrid
    unc = computeNodeDisp(ucc, tbls);
    unvec = reshape(unc, dim, [])';
    figure 
    coef = 1e0;
    plotGridDeformed(G, coef*unvec);
end

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% The MPSA-W module is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% The MPSA-W module is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with the MPSA-W module.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>

