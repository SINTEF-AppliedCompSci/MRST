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
grad_nodeface = TensorProd();
grad_nodeface.tbl1 = cellnodefacecoltbl;
grad_nodeface.tbl2 = cellnodefacecoltbl;
grad_nodeface.replacefds1 = {'coldim', 'rowdim'};
grad_nodeface.reducefds = {'faces'};
grad_nodeface.mergefds = {'cells', 'nodes'};
grad_nodeface = grad_nodeface.setup();

% The cellcol part of the grad operator from cellcoltbl to cellnodecolrowtbl is
% obtained for any u in cellcoltbl by using v = grad_cell.evalProd(greduced, u)
% where greduced and grad_cell are defined below 
%
fds = {'cells', 'nodes', 'coldim'};
greduced = - tblmap(g, cellnodefacecoltbl, cellnodecoltbl, fds); 
grad_cell = TensorProd();
grad_cell.tbl1 = cellnodecoltbl;
grad_cell.tbl2 = cellcoltbl;
grad_cell.replacefds1 = {'coldim', 'rowdim'};
grad_cell.mergefds = {'cells'};
grad_cell = grad_cell.setup();


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
% expression div_facenode.evalProd(d, u) where d and div_facenode are defined
% below
%
d = facetNormals; 
div_facenode = TensorProd();
div_facenode.tbl1 = cellnodefacecoltbl;
div_facenode.tbl2 = cellnodecolrowtbl;
div_facenode.replacefds1 = {'coldim', 'rowdim'};
div_facenode.reducefds = {'rowdim', 'cells'};
div_facenode.mergefds = {'nodes'};
div_facenode.mergefds = {'nodes'};
div_facenode = div_facenode.setup();

% the cellcol part of the divergence operator from cellnodecolrowtbl to
% cellcoltbl is obtained for any u in cellnodecolrowtbl by evaluating the
% expression div_facenode.evalProd(dreduced, u) where dreduced and div_facenode
% are defined below
%
fds = {'cells', 'nodes', 'coldim'};
% note the minus sign below (see formula in paper)
dereduced = - tblmap(facetNormals, cellnodefacecoltbl, cellnodecoltbl, fds);

div_cell = TensorProd();
div_cell.tbl1 = cellnodecoltbl;
div_cell.tbl2 = cellnodecolrowtbl;
div_cell.replacefds1 = {'coldim', 'rowdim'};
div_cell.reducefds   = {'rowdim', 'nodes'};
div_cell.mergefds    = {'cells'};
div_cell = div_cell.setup();



