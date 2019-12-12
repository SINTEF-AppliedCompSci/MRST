%% Linear pressure test
%
% We impose constant pressure on the bottom and top faces with a given pressure
% drop. We consider a constant permeability over the whole domain.
%
% The exact solution for pressure is a linear function of the z-variable.
%
% The VAG scheme is exact for linear pressure field. In this example, we check
% that this property holds. We use a twisted grid.


%% Load modules

mrstModule add mpfa vem vag vemmech

%% Setup a Cartesian grid

nx = 3;
ny = 3;
nz = 5;
cartDims = [nx, ny, nz];
physDims = [1, 1, 1];

G = cartGrid(cartDims, physDims);
% We twist the grid
G = twister(G, 0.1);
G = computeGeometry(G);

nc = G.cells.num;
nn = G.nodes.num;

perm = 1;
rock = makeRock(G, perm*ones(nc, 1), 0.1*ones(nc, 1));

%% Compute the VAG transmissibilities

vagstruct = computeVagTrans(G, rock);

%% Setup system matrix

[A, op] = setupSystem(vagstruct, G);

rhsfun = op.rhsfun;
computeCellPressure = op.computeCellPressure;

%% setup source and boundary conditions
% In this case, constant pressure at top and bottom faces with given pressure
% drop

% No volumetric source in this case
cellinflux = zeros(nc, 1);
nodeinflux = zeros(nn, 1);
rhs = rhsfun(nodeinflux, cellinflux);

% Get nodes on the top and bottom faces to impose boundary conditions
clear facetbl
nf = G.faces.num;
facetbl.faces = (1 : nf)';
facetbl.num = nf;

clear facenodetbl
facenodetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos)');
facenodetbl.nodes = G.faces.nodes;
facenodetbl.num = numel(facenodetbl.faces);

facez = G.faces.centroids(:, 3);
zmax = max(facez);
zmin = min(facez);
zmean = sum(G.cells.centroids(:, 3).*G.cells.volumes)/sum(G.cells.volumes);

topfacetbl.faces = find((facez - zmean) > (1 - 1e-10)*(zmax - zmean));
topfacetbl.num = numel(topfacetbl.faces);
[~, topfacenodetbl] = setupTableMapping(topfacetbl, facenodetbl, {'faces'});
topnodetbl = projTable(topfacenodetbl, {'nodes'});

bottomfacetbl.faces = find((facez - zmean) < (1 - 1e-10)*(zmin - zmean));
bottomfacetbl.num = numel(bottomfacetbl.faces);
[~, bottomfacenodetbl] = setupTableMapping(bottomfacetbl, facenodetbl, {'faces'});
bottomnodetbl = projTable(bottomfacenodetbl, {'nodes'});

ntop = topnodetbl.num;
nbottom = bottomnodetbl.num;
clear bcnodetbl
bcnodetbl.nodes = [topnodetbl.nodes; bottomnodetbl.nodes];
bcnodetbl.num   = numel(bcnodetbl.nodes);

pval = 100;
toppbc = pval*ones(ntop, 1); 
bottompbc = zeros(nbottom, 1);
pbc = [toppbc; bottompbc];

clear nodetbl;
nodetbl.nodes = (1 : nn)';
nodetbl.num = nn;

intnodes = true(nn, 1);
intnodes(bcnodetbl.nodes) = false;
intnodes = find(intnodes);
clear intnodetbl
intnodetbl.nodes = intnodes;
intnodetbl.num = numel(intnodes);

A11 = A(intnodetbl.nodes, intnodetbl.nodes);
A12 = A(intnodetbl.nodes, bcnodetbl.nodes);


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

%% Plot pressure versus z-coordinate
close all
figure
hold on
zn = G.nodes.coords(:, 3);
plot(zn, pn, '*');
zc = G.cells.centroids(:, 3);
plot(zc, pc, '*');
legend({'node pressures', 'cell pressures'});
title('pressure versus z-coordinate');
xlabel('z-coordinate');
ylabel('pressure');

%% Plot solution
figure
plotNodeData(G, pn)
title('nodal pressure')
view(-60, 20);
axis equal
colorbar

figure
plotCellData(G, pc)
title('cell pressure')
view(-60, 20);
axis equal
colorbar

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.
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
