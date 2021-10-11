%% VAG example
%  Example of VAG computation, with given rate injection in some cells and given
%  pressure at other cells. It uses the function incompVAG to solve the system

%% Load modules

mrstModule add mimetic incomp mpfa vem vag vemmech

%% Setup a Cartesian grid

nx = 10;
ny = 10;
% nx = 1;
% ny = 1;
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

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
