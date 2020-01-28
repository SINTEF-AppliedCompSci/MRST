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
nn  = G.nodes.num;
dim = G.griddim;

coltbl.coldim = (1 : dim)';
coltbl.num = dim;
rowtbl = coltbl;
rowtbl = replacefield(rowtbl, {'coldim', 'rowdim'});

celltbl.cells = (1 : nc)';
celltbl.num = nc;

nodetbl.nodes = (1 : nn)';
nodetbl.num = nn;

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

fds = {'cells', 'nodes', 'coldim'};
% note the minus sign below (see formula in paper)
dreduced = - tblmap(facetNormals, cellnodefacecoltbl, cellnodecoltbl, fds);

prod = TensorProd();
prod.tbl1 = cellnodecoltbl;
prod.tbl2 = cellnodecolrowtbl;
prod.replacefds1 = {'coldim', 'rowdim'};
prod.replacefds2 = {'coldim', 'rowdim', 'interchange'};
prod.reducefds   = {'rowdim', 'nodes'};
prod.mergefds    = {'cells'};
prod.prodtbl = cellcoltbl;
prod = prod.setup();

div_cell = SparseTensor();
div_cell = div_cell.setFromTensorProd(dreduced, prod);

divgrad_nodeface = multSparseTensor(div_nodeface, grad_nodeface);

colrowtbl = crossTable(coltbl, rowtbl, {});
nodecolrowtbl = crossTable(nodetbl, colrowtbl, {});

col2row2tbl.coldim2 = colrowtbl.coldim;
col2row2tbl.rowdim2 = colrowtbl.rowdim;
col2row2tbl.coldim1 = colrowtbl.rowdim;
col2row2tbl.rowdim1 = colrowtbl.coldim;
col2row2tbl.num = colrowtbl.num;

prod = TensorProd();
prod.tbl1 = col2row2tbl;
prod.tbl2 = nodecolrowtbl;
prod.replacefds1 = {{'coldim1', 'coldim'}, ...
                    {'rowdim1', 'rowdim'}};
prod.replacefds2 = {{'coldim', 'coldim2'}, ...
                    {'rowdim', 'rowdim2'}};
prod.reducefds = {'coldim2', 'rowdim2'};
prod.prodtbl = nodecolrowtbl;
prod = prod.setup();

sym_op = SparseTensor();
sym_op = sym_op.setFromTensorProd(ones(col2row2tbl.num, 1), prod);

% get nodal average of cellnode tensor
prod = TensorProd();
prod.tbl1 = celltbl;
prod.tbl2 = cellnodecolrowtbl;
prod.reducefds = {'cells'};
prod.prodtbl = nodecolrowtbl;

prod = prod.setup();

nodeaverage_op = SparseTensor();
nodeaverage_op = nodeaverage_op.setFromTensorProd(ones(celltbl.num), prod);

symnodeaverage_op = multSparseTensor(sym_op, nodeaverage_op);

% Defines stiffness operator

dim = G.griddim;
vdim = dim*(dim + 1)/2;
voigttbl.voigt = (1 : vdim)';
switch dim
  case 2
    voigttbl.coldim = [1; 2; 2];
    voigttbl.rowdim = [1; 2; 1];
    voigt = [1; 3; 3; 2];
  case 3
    voigttbl.coldim = [1; 2; 3; 3; 3; 2];
    voigttbl.rowdim = [1; 2; 3; 2; 1; 1];
    voigt = [1; 6; 5; 6; 2; 4; 5; 4; 3];
end
voigttbl.num = vdim;

colrowvoigttbl.coldim = colrowtbl.coldim;
colrowvoigttbl.rowdim = colrowtbl.rowdim;
colrowvoigttbl.voigt  = voigt;
colrowvoigttbl.num    = colrowtbl.num;

stensdim = vdim*(vdim + 1)/2;

s = 0;
voigt1 = [];
voigt2 = [];
stensind = [];


for i = 1 : vdim
    for j = i : vdim
        voigt1(end + 1) = i;
        voigt2(end + 1) = j;
        s = s + 1;
        stensind(end + 1) = s;
    end
end


stenstbl.voigt1 = voigt1';
stenstbl.voigt2 = voigt2';
stenstbl.stensind = stensind';
stenstbl.num = stensdim;
fds = {'stensind', 'voigt1', 'voigt2'};
stenstbl = sortTable(stenstbl, fds);

voigt2tbl = crossTable(voigtbl, voigtbl, 'crossextend', {'voigt1', ...
                    'voigt2'});


for i = 1 : voigt2tbl.num
    v1 = voigt2tbl.voigt1;
    v2 = voigt2tbl.voigt2;
    if v1 < v2;
        stensind(i) = 
    
% stensmat = convertTableToArray(stenstbl, {'stensind', 'voigt1', 'voigt2'});


% stensfullmat = [stensmat; ...
%             [stensmat(:, 1), stensmat(:, 3), stensmat(:, 2)]];

% stensfulltbl = convertArrayToTable(stensfullmat, {'stensind', 'voigt1', ...
%                     'voigt2'});

% stensfulltbl = crossTable(stensfulltbl, voigttbl, {{'voigt1', 'voigt'}});
% stensfulltbl = replacefield(stensfulltbl, {{'coldim', 'coldim1'}, {'rowdim', 'rowdim1'}});
% stensfulltbl = crossTable(stensfulltbl, voigttbl, {{'voigt2', 'voigt'}});
% stensfulltbl = replacefield(stensfulltbl, {{'coldim', 'coldim2'}, {'rowdim', 'rowdim2'}});
% stensfulltbl = rmfield(stensfulltbl, {'voigt1', 'voigt2'});

% fds = {'stensind', 'rowdim1', 'coldim1', 'rowdim2', 'coldim2'};
% stensfulltbl = sortTable(stensfulltbl, fds);
% stensfullmat = convertTableToArray(stensfulltbl, fds);

crossfds = {{'coldim', {'coldim1', 'coldim2'}}, ...
            {'rowdim', {'rowdim1', 'rowdim2'}}, ...
            {'voigt', {'voigt1', 'voigt2'}}};

col2row2voigt2tbl = crossTable(colrowvoigttbl, colrowvoigttbl, {}, 'crossextend', ...
                               crossfds);

fds = {'rowdim1', 'coldim1', 'voigt1', ...
       'rowdim2', 'coldim2', 'voigt2', ...
      };
col2row2voigt2mat = convertTableToArray(col2row2voigt2tbl, fds);

Cvoigt = (1 : stenstbl.num)';

% prod = TensorProd();
% prod.tbl1 = stenstbl;
% prod.tbl2 = voigttbl;
% prod.replacefds1 = {'stensind', ''};
% prod.replacefds2 = {{'voigt', 'voigt2'}, {'coldim', ''}, {'rowdim', ''}};
% prod.reducefds = {'voigt2'};
% prod = prod.setup();
% 
% Cvoigtmat = SparseMatrix();
% Cvoigtmat = Cvoigtmat.setFromTensorProd(Cvoigt, prod);
% Cvoigtmat = Cvoigtmat.matrix;

% C = tblmap(Cvoigt, stenstbl, stensfulltbl, {'stensind'});

C = tblmap(Cvoigt, stenstbl, col2row2voigt2tbl, {'voigt1', 'voigt2'});

return

prod = TensorProd();
prod.tbl1 = col2row2tbl;
prod.tbl2 = colrowtbl;
prod.replacefds1 = {{'coldim1', 'coldim'}, {'rowdim1', 'rowdim'}};
prod.replacefds2 = {{'coldim', 'coldim2'}, {'rowdim', 'rowdim2'}};
prod.reducefds = {'coldim2', 'rowdim2'};
prod = prod.setup();

Cmat = SparseMatrix();
Cmat = Cmat.setFromTensorProd(C, prod);
Cmat = Cmat.matrix;


