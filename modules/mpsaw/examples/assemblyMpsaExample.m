%% Assembly of MPSA-weak
%
% Reference paper:
% Finite volume methods for elasticity with weak symmetry
% Keilegavlen, Eirik and Nordbotten, Jan Martin
% International Journal for Numerical Methods in Engineering 2017

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
    N = 10;
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

G = twister(G, 0.1);
G = computeGeometry(G);
dim = G.griddim;

% Set material properties
nc = G.cells.num;
lambda = ones(nc, 1);
mu     = ones(nc, 1);

prop = struct('lambda', lambda, ...
              'mu', mu);

useVirtual = false;
[tbls, mappings] = setupStandardTables(G, 'useVirtual', useVirtual);
loadstruct = setupBCpercase(runcase, G, tbls, mappings);

casetype = 'standard';
switch casetype
  case 'blockassembly'
    assembly = blockAssembleMPSA(G, prop, loadstruct, eta, tbls, mappings, 'blocksize', 100, 'verbose', true);
  case 'standard'
    assembly = assembleMPSA(G, prop, loadstruct, eta, tbls, mappings, 'bcetazero', bcetazero, 'extraoutput', true, 'addAdOperators', true); 
  case 'experimental'
    assembly = assembleMPSA2(G, prop, loadstruct, eta, tbls, mappings, 'bcetazero', bcetazero, 'extraoutput', true);     
end

B   = assembly.B  ;
rhs = assembly.rhs;

sol = B\rhs;

% Recover displacement values at cell centers.
cellcoltbl = tbls.cellcoltbl;
n = cellcoltbl.num;
u = sol(1 : n);
% values of Lagrangian corresponding to Dirichlet contraint
lm = sol((n + 1) : end);

isstresscomputed = false;
if strcmp(casetype, 'standard')
    % We compute node-face displacement values
    op = assembly.adoperators;
    unf = op.facenodedispop(u, lm);
    % We compute stress
    stress = op.stressop(unf, u);
    assert(dim == 2, 'only 2d treated here at the moment');
    stress = formatField(stress, dim, 'stress');
    isstresscomputed = true;
end

plotdeformedgrid = true;
if plotdeformedgrid
    lagmult = sol(n + 1 : end);
    unode = computeNodeDisp(u, tbls);
    unvec = formatField(unode, dim, 'displacement');
    coef = 1e0;
    plotGridDeformed(G, coef*unvec);
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
    if isstresscomputed
        figure
        plotCellData(G, stress(:, 1));
        titlestr = sprintf('stress xx', eta);
        title(titlestr)
        colorbar
        figure
        plotCellData(G, stress(:, 2));
        titlestr = sprintf('stress yy', eta);
        title(titlestr)
        colorbar
        figure
        plotCellData(G, stress(:, 3)); 
        titlestr = sprintf('stress xy', eta);
        title(titlestr)
        colorbar
    end
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

%% plot displacement values as function of z. 
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

