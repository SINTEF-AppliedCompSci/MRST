mrstModule add mimetic mpfa incomp

clear all

isverbose = true;
eta       = 1/3;
blocksize = 10;

%% Define and process geometry
% Construct a Cartesian grid of size 10-by-10-by-4 cells, where each cell
% has dimension 1-by-1-by-1. Because our flow solvers are applicable for
% general unstructured grids, the Cartesian grid is here represented using
% an unstructured formate in which cells, faces, nodes, etc. are given
% explicitly.
nx = 10; ny = 10;
nx = 1; ny = 1;
G = cartGrid([nx, ny]);
% G = twister(G, 0.1);
G = computeGeometry(G);
nc = G.cells.num;

%% Set rock and fluid data
% The only parameters in the single-phase pressure equation are the
% permeability $K$, which here is homogeneous, isotropic and equal 100 mD.
% The fluid has density 1000 kg/m^3 and viscosity 1 cP.
% We make a non diagonal rock tensor
rock = makeRock(G, 1e-3*darcy, 1);


nc  = G.cells.num;
nf  = G.faces.num;
dim = G.griddim;

coltbl.coldim = (1 : dim)';
coltbl.num = dim;
rowtbl = coltbl;
rowtbl = replacefield(rowtbl, {'coldim', 'rowdim'});

celltbl.cells = (1 : nc)';
celltbl.num = nc;

cellcoltbl = crossTable(celltbl, coltbl, {}); % ordering is cell - col

cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos)); 
cellfacetbl.faces = G.cells.faces(:, 1);
cellfacetbl.num   = numel(cellfacetbl.cells);

facenodetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos)); 
facenodetbl.nodes = G.faces.nodes;
facenodetbl.num   = numel(facenodetbl.faces);
% We setup the face-node table and it is ordered along ascending node numbers so
% that we will have a block structure for the nodal scalar product.
facenodetbl = sortTable(facenodetbl, {'nodes', 'faces'});

cellnodefacetbl = crossTable(cellfacetbl, facenodetbl, {'faces'});

% We setup the cell-face-node table, cellnodefacetbl. Each entry determine a
% unique facet in a corner
% We order cellnodeface in cell-node-face order. This is node to optimize
% for-end loop below.
cellnodefacetbl = sortTable(cellnodefacetbl, {'cells', 'nodes', 'faces'});

% We setup the cell-node table, cellnodetbl. Each entry determine a unique
% corner
cellnodetbl = projTable(cellnodefacetbl, {'nodes', 'cells'});
cellnodecoltbl = crossTable(cellnodetbl, coltbl, {});
 
% shortcuts:
fno = cellnodefacetbl.faces;
cno = cellnodefacetbl.cells;
nno = cellnodefacetbl.nodes;

cellnodefacecoltbl = crossTable(cellnodefacetbl, coltbl, {});
cellnodefacecoltbl = sortTable(cellnodefacecoltbl, {'cells', 'nodes', ...
                    'faces', 'coldim'});


cellFacetVec = G.faces.centroids(fno, :) - G.cells.centroids(cno, :) + ...
    eta*(G.nodes.coords(nno, :) - G.faces.centroids(fno, :));

cellFacetVec = reshape(cellFacetVec', [], 1);


cellnodefacemat = convertTableToArray(cellnodefacetbl, {'cells', 'nodes', ...
                    'faces'});
cellnodefacefcolmat = [cellnodefacemat, repmat([1 : coltbl.num]', cellnodetbl.num)];
cellnodefacefcoltbl = convertArrayToTable(cellnodefacefcolmat, {'cells', 'nodes', ...
                    'faces', 'fakecoldim'});
cellnodefacefcolcoltbl = crossTable(cellnodefacefcoltbl, coltbl, {});

cellnodefcolcoltbl = projTable(cellnodefacefcolcoltbl, {'cells', 'nodes', ...
                    'fakecoldim', 'coldim'});

% we map cellFacetVec from cellnodefacecoltbl to cellnodefcolcoltbl
% (These mappings probably in practice correspond to identity. When
% optimizing the code we should remove them).
fds = {'cells', 'nodes', 'faces', 'coldim'};
cellFacetVec = tblmap(cellFacetVec, cellnodefacecoltbl, cellnodefacefcolcoltbl, fds);
fds = {'cells', 'nodes', 'fakecoldim', 'coldim'};
cellFacetVec = tblmap(cellFacetVec, cellnodefacefcolcoltbl, cellnodefcolcoltbl, fds);


fds = {{'cells' , {'cells1', 'cells2'}}, ...
       {'nodes' , {'nodes1', 'nodes2'}}, ...
       {'coldim', {'fakecoldim', 'coldim2'}}};

[I, cell2node2fcolcoltbl] = setupIdentity(cellnodecoltbl, fds);

prod = TensorProd();
prod.tbl1 = cellnodefcolcoltbl;
prod.tbl2 = cell2node2fcolcoltbl;
prod.replacefds1 = {{'cells'     , 'cells1'}     , ...
                    {'coldim'    , 'coldim1'}    , ...
                    {'nodes'     , 'nodes1'}};
prod.mergefds  = {'cells1', 'nodes1'};
prod.reducefds = {'fakecoldim'};
prod = prod.setup();

res = prod.evalProd(cellFacetVec, I);
prodtbl = prod.prodtbl;

cellnodecoltbl = addLocInd(cellnodecoltbl, 'locind');

fds = {{'cells' , 'cells1'}, ...
       {'nodes' , 'nodes1'}, ...
       {'coldim', 'coldim1'}};
ind1 = tblmap(cellnodecoltbl.locind, cellnodecoltbl, prodtbl, fds);
fds = {{'cells' , 'cells2'}, ...
       {'nodes' , 'nodes2'}, ...
       {'coldim', 'coldim2'}};
ind2 = tblmap(cellnodecoltbl.locind, cellnodecoltbl, prodtbl, fds);

n = cellnodecoltbl.num;
A = sparse(ind1, ind2, res, n, n);

opt.invertBlocks = 'mex';
bi = blockInverter(opt);

sz = repmat(coltbl.num, cellnodetbl.num, 1);
invA = bi(A, sz);

[locind1, locind2, g] = find(invA);

gtbl.locind1 = locind1;
gtbl.locind2 = locind2;
gtbl.num = numel(gtbl.locind1);

[gtbl, indstruct] = crossTable(gtbl, cellnodecoltbl, {{'locind1', 'locind'}});
g = tbldispatch1(g, indstruct);
gtbl = replacefield(gtbl, {'coldim', 'fakecoldim'});
[gtbl, indstruct] = crossTable(gtbl, cellnodecoltbl, {{'locind2', 'locind'}});
g = tbldispatch1(g, indstruct);

% we map g from  gtbl to cellnodefacecoltbl
% (These mappings probably in practice are equivalent to identity. When
% optimizing the code we should remove them).
fds = {'cells', 'nodes', 'fakecoldim', 'coldim'};
g = tblmap(g, gtbl, cellnodefacefcolcoltbl, fds);
fds = {'cells', 'nodes', 'faces', 'coldim'};
g = tblmap(g, cellnodefacefcolcoltbl, cellnodefacecoltbl, fds);

return

%% some test

clear seltbl
seltbl.cells1 = 1;
seltbl.nodes1 = 1;
seltbl.cells2 = 1;
seltbl.nodes2 = 1;
seltbl.num = 1;

fds    = {'cells1', 'cells2', 'nodes1', 'nodes2'}
seltbl = crossTable(prodtbl, seltbl, fds);

clear sseltbl
sseltbl.cells = 1;
sseltbl.nodes = 1;
sseltbl.num = 1;
sselcoltbl = crossTable(sseltbl, coltbl, {});
fds = {'cells', 'nodes', 'coldim'};
sselcoltbl = crossTable(sselcoltbl, cellnodecoltbl, fds);
sselcolmat = convertTableToArray(sselcoltbl, {fds{:}, 'locind'})

fds = {{'cells' , 'cells1'}, ...
       {'nodes' , 'nodes1'}, ...
       {'coldim', 'coldim1'}};
selcoltbl = crossTable(cellnodecoltbl, seltbl, fds);


fds = {'cells1' , 'cells2', ...
       'nodes1' , 'nodes2', ...
       'coldim1', 'coldim2'};

selmat = convertTableToArray(seltbl, fds);


fds = {{'cells' , 'cells1'}, ...
       {'nodes' , 'nodes1'}, ...
       {'coldim', 'coldim1'}};
tbl1 = crossTable(cellnodecoltbl, seltbl, fds);
tbl1 = replacefield(tbl1, fds);
fds = {{'cells' , 'cells2'}, ...
       {'nodes' , 'nodes2'}, ...
       {'coldim', 'coldim2'}};
tbl2 = crossTable(cellnodecoltbl, seltbl, fds);
tbl2 = replacefield(tbl2, fds);

fds = {'cells1' , 'cells2', ...
       'nodes1' , 'nodes2', ...
       'coldim1', 'coldim2', 'locind'};

tblmat1 = convertTableToArray(tbl1, fds)
tblmat2 = convertTableToArray(tbl2, fds)


%%

fds = {{'cells1', 'cells'}, ...
       {'nodes1', 'nodes'}, ...
       {'faces1', 'faces'}};

[tbl, indstruct] = crossTable(prodtbl, cellnodefacefcoltbl, fds)

return

numnodes = double(diff(G.faces.nodePos));
numnodes = numnodes(fno);
facetNormals = G.faces.normals(fno, :);
facetNormals = bsxfun(@ldivide, numnodes, facetNormals);

sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
facetNormals = sgn.*facetNormals; % Outward normals with respect to cell
                                  % in cellnodeface.

facetNormals = reshape(facetNormals', [], 1);

                    
if opt.verbose
    fprintf('assemble facet K*normals ...\n');
end
% Assemble facePermNormals which corresponds to $Kn$ where n are the *outward*
% normals at the facets.
[perm, r, c] = permTensor(rock, G.griddim);
permmat = perm;
perm = reshape(permmat', [], 1);
% setup cellcolrow table for the vector perm
[~, colrowtbl] = setupTableMapping(coltbl, rowtbl, []);
[~, cellcolrowtbl] = setupTableMapping(colrowtbl, celltbl, []);
cellcolrowtbl = sortTable(cellcolrowtbl, {'cells', 'coldim', 'rowdim'});
cellcolrowtbl = addLocInd(cellcolrowtbl, 'ccrind');

% dispatch perm on cellnodeface

[~, cellnodefacecolrowtbl] = setupTableMapping(cellcolrowtbl, cellnodefacetbl, ...
                                                             {'cells'});
op = setupTableMapping(cellcolrowtbl, cellnodefacecolrowtbl, {'ccrind'}, ...
                                     'fastunstable', true);
perm = op*perm;
% Multiply perm with facetNormals
map1 = setupTableMapping(cellnodefacecoltbl, cellnodefacecolrowtbl, ...
                                       {'cnfind', 'coldim'}, 'fastunstable', ...
                                       true);
Kn = perm.*(map1*facetNormals);

map2 = setupTableMapping(cellnodefacecolrowtbl, cellnodefacecoltbl, ...
                                       {'cnfind', {'rowdim', 'coldim'}}, ...
                                       'fastunstable', true);
Kn = map2*Kn;

% store Kn in matrix form in facePermNormals.
op = setupTableMapping(cellnodefacetbl, cellnodefacecoltbl, {'cnfind'}, ...
                                     'fastunstable', true);
ind1 = cellnodefacetbl.cnfind;
ind1 = op*ind1;
op = setupTableMapping(coltbl, cellnodefacecoltbl, {'coldim'}, 'fastunstable', ...
                               true);
ind2 = (1 : coltbl.num)';
ind2 = op*ind2;    
facePermNormals = sparse(ind1, ind2, Kn, cellnodefacetbl.num, coltbl.num);

% Some shortcuts
cno = cellnodefacetbl.cells;
fno = cellnodefacetbl.faces;
nno = cellnodefacetbl.nodes;



% set up areas and volumes
areas = G.faces.areas(fno);
vols  = G.cells.volumes(cno);

map = setupTableMapping(cellnodetbl, cellnodefacetbl, {'cells', 'nodes'}); 
nfaces = diag(map'*map);
nfaces = map*nfaces;

cnf_i = 1; % start indice for the cellnodefacetbl index
mat_i = 1; % start indice for the mattbl index
    