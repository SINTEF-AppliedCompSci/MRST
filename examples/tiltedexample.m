clear all
% close all

mrstModule add mpsaw vemmech mpfa

N = 10;
G = cartGrid([N, N], [1 , 1]);
angle = 10/180*pi;
[G, bcfaces] = rotateGrid(G, angle);
G = computeGeometry(G);

dim = G.griddim;
nc = G.cells.num;

lambda = ones(nc, 1);
mu     = ones(nc, 1);

prop = struct('lambda', lambda, ...
              'mu', mu);

useVirtual = false;
[tbls, mappings] = setupStandardTables(G, 'useVirtual', useVirtual);

extfaces{1} = bcfaces.ymin;
n = numel(extfaces{1});
linforms{1} = repmat([sin(angle), -cos(angle)], n, 1);
linformvals{1} = zeros(n, 1);

extfaces{2} = bcfaces.xmin;
n = numel(extfaces{2});
linforms{2} = repmat([cos(angle), sin(angle)], n, 1);
linformvals{2} = zeros(n, 1);

nodefacecoltbl  = tbls.nodefacecoltbl;
coltbl  = tbls.coltbl;
facetbl.faces = (1 : G.faces.num)';
facetbl = IndexArray(facetbl);
facecoltbl  = crossIndexArray(facetbl, coltbl, {});

extfacetbl.faces = bcfaces.ymax;
extfacetbl = IndexArray(extfacetbl);
extfacecoltbl = crossIndexArray(extfacetbl, coltbl, {});

normals = G.faces.normals;
normals = reshape(normals', [], 1);

map = TensorMap();
map.fromTbl = facecoltbl;
map.toTbl   = extfacecoltbl;
map.mergefds = {'faces', 'coldim'};
map = map.setup();
extFacetNormals = map.eval(normals);

map = TensorMap();
map.fromTbl = extfacecoltbl;
map.toTbl = nodefacecoltbl;
map.mergefds = {'faces', 'coldim'};

map = map.setup();

extforce = map.eval(-extFacetNormals);

cellcoltbl = tbls.cellcoltbl;
force = zeros(cellcoltbl.num, 1);

bc.extfaces    = vertcat(extfaces{:});
bc.linform     = vertcat(linforms{:});
bc.linformvals = vertcat(linformvals{:});

loadstruct.bc = bc;
loadstruct.extforce = extforce;
loadstruct.force = force;

eta = 1/3;
bcetazero = false;

assembly = assembleMPSA(G, prop, loadstruct, eta, tbls, mappings, 'bcetazero', ...
                        bcetazero);

B   = assembly.B  ;
rhs = assembly.rhs;

sol = B\rhs;

% displacement values at cell centers.
n = cellcoltbl.num;

ucell = sol(1 : n);

op = setupCell2NodeDispOperator(G, tbls);
unode = op*ucell;

dim = G.griddim;
ucell = reshape(ucell, dim, [])';
unode = reshape(unode, dim, [])';

figure
plotGrid(G, 'facecolor', 'none', 'edgecolor', 'blue');
plotGridDeformed(G, unode, 'facecolor', 'none', 'edgecolor', 'red');
axis equal

figure
plotCellData(G, ucell(:, 1));
title('disp x direction');
colorbar
axis equal

figure
plotCellData(G, ucell(:, 2));
title('disp y direction');
colorbar
axis equal