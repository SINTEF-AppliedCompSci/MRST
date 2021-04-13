%% Linear pressure test
%
%
% The MPFA method is exact for linear pressure field. We impose boundary
% condition such that the exact solution is linear. The grid is a twisted
% grid made from a Cartesian grid.

mrstModule add ad-core ad-props incomp mrst-gui mpfa postprocessing mpsaw

clear all
close all

dimcase = 2;
switch dimcase
  case 2
    nx = 10; ny = 10;
    G = cartGrid([nx, ny]);
  case 3
    nx = 4; ny = 3; nz = 5;
    G = cartGrid([nx, ny, nz]);
end
G = twister(G, 0.1); 
G = computeGeometry(G);
nc = G.cells.num;


%% setup bc for pressure
dim = G.griddim;
zface = G.faces.centroids(:, dim);
maxz = max(zface);
minz = min(zface);
bottomfaces = find(G.faces.centroids(:, dim) > (maxz - 1e-4));
topfaces = find(G.faces.centroids(:, dim) < (minz + 1e-4));

nbf = numel(bottomfaces);
ntf = numel(topfaces);
bcfaces = [bottomfaces; topfaces];
bcvals = zeros(nbf + ntf, 1);
bcvals(1 : nbf) = 1*barsa;

useVirtual = false;
[tbls, mappings] = setupStandardTables(G, 'useVirtual', useVirtual);
nodefacetbl = tbls.nodefacetbl;

bcfacetbl.faces = bcfaces;
bcfacetbl = IndexArray(bcfacetbl);
bcnodefacetbl = crossIndexArray(nodefacetbl, bcfacetbl, {'faces'});

map = TensorMap();
map.fromTbl = bcfacetbl;
map.toTbl = bcnodefacetbl;
map.mergefds = {'faces'};
map = map.setup();

bcvals = map.eval(bcvals);

bcdirichlet = struct('bcnodefacetbl', bcnodefacetbl, ...
                     'bcvals', bcvals);
bcneumann = [];

bcstruct = struct('bcdirichlet', bcdirichlet, ...
                  'bcneumann'  , bcneumann);

src = []; % no source in this case

% setup identity K in cellcolrowtbl
cellcolrowtbl = tbls.cellcolrowtbl;
colrowtbl = tbls.colrowtbl;

map = TensorMap();
map.fromTbl = colrowtbl;
map.toTbl = cellcolrowtbl;
map.mergefds = {'coldim', 'rowdim'};
map = map.setup();

K = [1; 0; 0; 1];
K = map.eval(K);

eta = 1/3;

doblock = false;
if doblock
    assembly = blockAssembleMPFA(G, K, bcstruct, src, eta, tbls, mappings, 'blocksize', 20, 'verbose', true);
else
    assembly = assembleMPFA(G, K, bcstruct, src, eta, tbls, mappings);
end

B = assembly.B;
rhs = assembly.rhs;

sol = B\rhs;

p = sol(1 : nc);

%% plotting

figure
plotCellData(G, p);

z = G.cells.centroids(:, dim);
vec = [z, p];
vec = sortrows(vec);
figure
plot(z, p);
xlabel('z');
ylabel('pressure');
title('linear test');

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

