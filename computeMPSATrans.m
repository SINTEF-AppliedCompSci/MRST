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

nodefacetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos)); 
nodefacetbl.nodes = G.faces.nodes;
nodefacetbl.num   = numel(nodefacetbl.faces);
% We setup the face-node table and it is ordered along ascending node numbers so
% that we will have a block structure for the nodal scalar product.
nodefacetbl = sortTable(nodefacetbl, {'nodes', 'faces'});
nodefacecoltbl = crossTable(nodefacetbl, coltbl, {});

cellnodefacetbl = crossTable(cellfacetbl, nodefacetbl, {'faces'});

% We setup the cell-face-node table, cellnodefacetbl. Each entry determine a
% unique facet in a corner
% We order cellnodeface in cell-node-face order. This is node to optimize
% for-end loop below.
cellnodefacetbl = sortTable(cellnodefacetbl, {'cells', 'nodes', 'faces'});
cellnodefacetbl = addLocInd(cellnodefacetbl, 'cnfind');

% We setup the cell-node table, cellnodetbl. Each entry determine a unique
% corner
cellnodetbl = projTable(cellnodefacetbl, {'nodes', 'cells'});

cellnodecoltbl = crossTable(cellnodetbl, coltbl, {});
cellnodecoltbl = sortTable(cellnodecoltbl, {'cells', 'nodes', 'coldim'});
cellnodecolrowtbl = crossTable(cellnodecoltbl, rowtbl, {});
cellnodecoltbl = addLocInd(cellnodecoltbl, 'cncind');

% shortcuts:
fno = cellnodefacetbl.faces;
cno = cellnodefacetbl.cells;
nno = cellnodefacetbl.nodes;

cellnodefacecoltbl = crossTable(cellnodefacetbl, coltbl, {});

fds = {'cells', 'nodes', 'coldim'};
cellnodefacecoltbl = crossTable(cellnodefacecoltbl, cellnodecoltbl, fds);

% this sorting may be unecessary. We do it to be sure
fds = {'cells', 'nodes', 'faces', 'coldim', 'cnfind', 'cncind'};
cellnodefacecoltbl = sortTable(cellnodefacecoltbl, fds);

cellFacetVec = G.faces.centroids(fno, :) - G.cells.centroids(cno, :) + ...
    eta*(G.nodes.coords(nno, :) - G.faces.centroids(fno, :));

cellFacetVec = reshape(cellFacetVec', [], 1);


fds = {'cells', 'nodes', 'faces'};
cnf = cellnodefacetbl;
ind1 = tblmap(cnf.cnfind, cnf, cellnodefacecoltbl, fds);

fds = {'cells', 'nodes', 'coldim'};
cnc = cellnodecoltbl;
ind2 = tblmap(cnc.cncind, cnc, cellnodefacecoltbl, fds);


n = cellnodecoltbl.num; 
assert(n == cellnodefacetbl.num, ['This implementation of mpsaw cannot handle ' ...
                    'this grid']);

A = sparse(ind1, ind2, cellFacetVec, n, n);

opt.invertBlocks = 'mex';
bi = blockInverter(opt);

sz = repmat(coltbl.num, cellnodetbl.num, 1);
invA = bi(A, sz);

[cncind, cnfind, g] = find(invA);

gtbl.cncind = cncind;
gtbl.cnfind = cnfind;
gtbl.num = numel(gtbl.cncind);

fds = {'cncind', 'cnfind'};
g = tblmap(g, gtbl, cellnodefacecoltbl, fds);

% We clean-up the tables

cellnodefacecoltbl = rmfield(cellnodefacecoltbl, {'cnfind', 'cncind'});
cellnodecoltbl     = rmfield(cellnodecoltbl, {'cncind'});

% The cellnodefacecol part of the grad operator from cellnodefacecoltbl to
% cellnodecolrowtbl is obtained for any u in cellnodefacecoltbl by using v =
% grad_nodeface.evalProd(g, u) where grad_nodeface is defined below
%
prod = TensorProd();
prod.tbl1 = cellnodefacecoltbl;
prod.tbl2 = nodefacecoltbl;
prod.replacefds1 = {'coldim', 'rowdim'};
prod.reducefds = {'faces'};
prod.mergefds = {'nodes'};
prod.prodtbl = cellnodecolrowtbl;
prod = prod.setup();

grad_nodeface = SparseTensor();
grad_nodeface = grad_nodeface.setFromTensorProd(g, prod);

% The cellcol part of the grad operator from cellcoltbl to cellnodecolrowtbl is
% obtained for any u in cellcoltbl by using v = grad_cell.evalProd(greduced, u)
% where greduced and grad_cell are defined below 
%
fds = {'cells', 'nodes', 'coldim'};
greduced = - tblmap(g, cellnodefacecoltbl, cellnodecoltbl, fds); 
prod = TensorProd();
prod.tbl1 = cellnodecoltbl;
prod.tbl2 = cellcoltbl;
prod.replacefds1 = {'coldim', 'rowdim'};
prod.mergefds = {'cells'};
prod = prod.setup();

grad_cell = SparseTensor();
grad_cell = grad_cell.setFromTensorProd(greduced, prod);

%% setup the facet normals
fno = cellnodefacetbl.faces;
cno = cellnodefacetbl.cells;
numnodes = double(diff(G.faces.nodePos));
numnodes = numnodes(fno);
facetNormals = G.faces.normals(fno, :);
facetNormals = bsxfun(@ldivide, numnodes, facetNormals);

sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
facetNormals = sgn.*facetNormals; % Outward normals with respect to cell
                                  % in cellnodefacetbl.
facetNormals = reshape(facetNormals', [], 1);

% The nodefacecol part of the divergence operator from cellnodecolrowtbl to
% nodefacecoltbl is obtained for any u in cellnodecolrowtbl by evaluating the
% expression div_nodeface.evalProd(d, u) where d and div_nodeface are defined
% below
%
d = facetNormals; 
prod = TensorProd();
prod.tbl1 = cellnodefacecoltbl;
prod.tbl2 = cellnodecolrowtbl;
prod.replacefds1 = {'coldim', 'rowdim'};
prod.reducefds = {'rowdim', 'cells'};
prod.mergefds = {'nodes'};
prod.prodtbl = nodefacecoltbl;
prod = prod.setup();

div_nodeface = SparseTensor();
div_nodeface = div_nodeface.setFromTensorProd(d, prod);

% the cellcol part of the divergence operator from cellnodecolrowtbl to
% cellcoltbl is obtained for any u in cellnodecolrowtbl by evaluating the
% expression div_cell.evalProd(dreduced, u) where dreduced and div_cell
% are defined below
%

warning('div_cell to be fixed...')
fds = {'cells', 'nodes', 'coldim'};
% note the minus sign below (see formula in paper)
dreduced = - tblmap(facetNormals, cellnodefacecoltbl, cellnodecoltbl, fds);

prod = TensorProd();
prod.tbl1 = cellnodecoltbl;
prod.tbl2 = cellnodecolrowtbl;
prod.replacefds1 = {'coldim', 'rowdim'};
prod.reducefds   = {'rowdim', 'nodes'};
prod.mergefds    = {'cells'};
prod = prod.setup();

div_cell = SparseTensor();
div_cell = div_cell.setFromTensorProd(dreduced, prod);


divgrad_nodeface = multSparseTensor(div_nodeface, grad_nodeface);

divgradmat = divgrad_nodeface.getMatrix();

return


prod = TensorProd();
prod.tbl1 = nodefacecoltbl;
prod.tbl2 = cellnodefacecoltbl;
prod.mergefds = {'faces', 'nodes', 'coldim'};
prod.prodtbl = cellnodefacecoltbl;
prod = prod.setup();

tens = SparseTensor();
val = ones(cellnodefacecoltbl.num, 1);
tens = tens.setFromTensorProd(val, prod, 'argindex', 2);


divgrad_nodeface = multSparseTensor(divgrad_nodeface, tens);

mat = divgrad_nodeface.getMatrix();