clear all
close all

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
    N = 40;
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

% set material properties
nc = G.cells.num;

lambda = ones(nc, 1);
mu     = ones(nc, 1);

prop = struct('lambda', lambda, ...
              'mu', mu);

useVirtual = true;

[tbls, mappings] = setupMpsaStandardTables(G, 'useVirtual', useVirtual);
loadstruct = setupBCpercase(runcase, G, tbls, mappings, [], 'useVirtual', useVirtual);

mech.prop = prop;
mech.loadstruct = loadstruct;

%% Setup model
model = MechModel(G, mech);

state = model.solveMechanics();

u = model.getProp(state, 'displacement');
u = formatField(u, dim, 'u');

unf = model.getProp(state, 'FaceNodeDisplacement');

cellnodefacevectbl = tbls.cellnodefacevectbl;
nodefacevectbl     = tbls.nodefacevectbl;
cellvectbl         = tbls.cellvectbl;
vectbl             = tbls.vectbl;

map = TensorMap();
map.fromTbl  = nodefacevectbl;
map.toTbl    = cellnodefacevectbl;
map.mergefds = {'nodes', 'faces', 'vec'};

if useVirtual

    nodefacetbl     = tbls.nodefacetbl;
    cellnodefacetbl = tbls.cellnodefacetbl;
    celltbl         = tbls.celltbl;

    map.pivottbl = cellnodefacevectbl;

    [vec, i] = ind2sub([vectbl.num, cellnodefacetbl.num], (1 : cellnodefacevectbl.num)');
    map.dispind1 = sub2ind([vectbl.num, nodefacetbl.num], vec, mappings.nodeface_from_cellnodeface(i));
    map.dispind2 = (1 : cellnodefacevectbl.num)';
    map.issetup = true;
else
    map = map.setup();
end

ucnf = map.eval(unf);

map = TensorMap();
map.fromTbl  = cellnodefacevectbl;
map.toTbl    = cellvectbl;
map.mergefds = {'cells', 'vec'};
if useVirtual

    map.pivottbl = cellnodefacevectbl;

    [vec, i] = ind2sub([vectbl.num, cellnodefacetbl.num], (1 : cellnodefacevectbl.num)');
    map.dispind1 = (1 : cellnodefacevectbl.num)';
    map.dispind2 = sub2ind([vectbl.num, celltbl.num], vec, mappings.cell_from_cellnodeface(i));
    map.issetup = true;
    
else
    map = map.setup();
end

unf = map.eval(ucnf);

unf = formatField(unf, dim, 'u');

stress = model.getProp(state, 'Stress');
stress = formatField(stress, dim, 'stress');

cdiv = model.getProp(state, 'ConsistentDiv');

figure
plotCellData(G, u(:, 1))
title('x-displacement')

figure
plotCellData(G, u(:, 2))
title('y-displacement')

figure
plotCellData(G, stress(:, 1))
title('xx-stress')

figure
plotCellData(G, stress(:, 2))
title('yy-stress')

figure
plotCellData(G, stress(:, 3))
title('xy-stress')

figure
plotCellData(G, cdiv)
title('consistent div')

figure
plotCellData(G, unf(:, 1))
title('face node x-displacement')

figure
plotCellData(G, unf(:, 2))
title('face node y-displacement')

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

