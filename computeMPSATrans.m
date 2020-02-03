%% Assembly of MPSA-weak
%%
%% Reference paper:
%% Finite volume methods for elasticity with weak symmetry
%% Keilegavlen, Eirik and Nordbotten, Jan Martin
%% International Journal for Numerical Methods in Engineering
%% 2017

clear all

% load modules
mrstModule add mimetic mpfa incomp

eta = 1/3;

%% Define and process geometry
% Construct a Cartesian grid 
dimcase = 2;
switch dimcase
  case 2
    nx = 10; ny = 10;
    G = cartGrid([nx, ny]);
  case 3
    nx = 5; ny = 5; nz = 5;
    G = cartGrid([nx, ny, nz]);
end
% G = twister(G, 0.1);
% compute Grid geometry
G = computeGeometry(G);

%% Basic table setup
%

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

% We setup boundary face table
intfaces = find(all(G.faces.neighbors, 2));
intfacetbl.faces = intfaces;
intfacetbl.num = numel(intfaces);

% We setup the cell-face-node table, cellnodefacetbl. Each entry determine a
% unique facet in a corner
% We order cellnodeface in cell-node-face order. This is node to optimize
% for-end loop below.
cellnodefacetbl = crossTable(cellfacetbl, nodefacetbl, {'faces'});
cellnodefacetbl = sortTable(cellnodefacetbl, {'cells', 'nodes', 'faces'});
cellnodefacetbl = addLocInd(cellnodefacetbl, 'cnfind');

% We setup the cell-node table, cellnodetbl. Each entry determine a unique
% corner
cellnodetbl = projTable(cellnodefacetbl, {'nodes', 'cells'});

cellnodecoltbl = crossTable(cellnodetbl, coltbl, {});
cellnodecoltbl = sortTable(cellnodecoltbl, {'cells', 'nodes', 'coldim'});
cellnodecolrowtbl = crossTable(cellnodecoltbl, rowtbl, {});
cellnodecoltbl = addLocInd(cellnodecoltbl, 'cncind');

cellnodefacecoltbl = crossTable(cellnodefacetbl, coltbl, {});

fds = {'cells', 'nodes', 'coldim'};
cellnodefacecoltbl = crossTable(cellnodefacecoltbl, cellnodecoltbl, fds);

% this sorting may be unecessary. We do it to be sure
fds = {'cells', 'nodes', 'faces', 'coldim', 'cnfind', 'cncind'};
cellnodefacecoltbl = sortTable(cellnodefacecoltbl, fds);


%% Construction of tensor g (as defined in paper eq 4.1.2)
% shortcuts:
%

fno = cellnodefacetbl.faces;
cno = cellnodefacetbl.cells;
nno = cellnodefacetbl.nodes;

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

%% Construction of the gradient operator
%
% We clean-up the tables

cellnodefacecoltbl = rmfield(cellnodefacecoltbl, {'cnfind', 'cncind'});
cellnodecoltbl     = rmfield(cellnodecoltbl, {'cncind'});

% Construction of gradnodeface_op : nodefacecoltbl -> cellnodecolrowtbl
%
% The nodefacecol part of the grad operator from nodefacecoltbl to
% cellnodecolrowtbl is obtained for any u in nodefacecoltbl by using v =
% prod.evalProd(g, u) where prod is defined below
%
prod = TensorProd();
prod.tbl1 = cellnodefacecoltbl;
prod.tbl2 = nodefacecoltbl;
prod.replacefds2 = {'coldim', 'rowdim'};
prod.reducefds = {'faces'};
prod.mergefds = {'nodes'};
prod.prodtbl = cellnodecolrowtbl;
prod = prod.setup();

gradnodeface_T = SparseTensor();
gradnodeface_T = gradnodeface_T.setFromTensorProd(g, prod);

% Construction of gradcell_T : cellcoltbl -> cellnodecolrowtbl
%
% The cellcol part of the grad operator from cellcoltbl to cellnodecolrowtbl is
% obtained for any u in cellcoltbl by using v = prod.evalProd(greduced, u)
% where greduced and prod are defined below 
%
fds = {'cells', 'nodes', 'coldim'};
greduced = - tblmap(g, cellnodefacecoltbl, cellnodecoltbl, fds); 
prod = TensorProd();
prod.tbl1 = cellnodecoltbl;
prod.tbl2 = cellcoltbl;
prod.replacefds2 = {'coldim', 'rowdim'};
prod.mergefds = {'cells'};
prod.prodtbl = cellnodecolrowtbl;
prod = prod.setup();

gradcell_T = SparseTensor();
gradcell_T = gradcell_T.setFromTensorProd(greduced, prod);

%% Construction of the divergence operator
%
% setup the facet normals
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

% divnodeface_T : cellnodecolrowtbl -> nodefacecoltbl
%
% The nodefacecol part of the divergence operator from cellnodecolrowtbl to
% nodefacecoltbl is obtained for any u in cellnodecolrowtbl by evaluating the
% expression divnodeface_T.evalProd(d, u) where d and divnodeface_T are defined
% below
%
d = facetNormals; 
prod = TensorProd();
prod.tbl1 = cellnodefacecoltbl;
prod.tbl2 = cellnodecolrowtbl;
prod.replacefds1 = {'coldim', 'rowdim'};
prod.replacefds2 = {'coldim', 'rowdim', 'interchange'};
prod.reducefds = {'rowdim', 'cells'};
prod.mergefds = {'nodes'};
prod.prodtbl = nodefacecoltbl;
prod = prod.setup();

divnodeface_T = SparseTensor();
divnodeface_T = divnodeface_T.setFromTensorProd(d, prod);

% divcell_T : cellnodecoltbl -> cellcoltbl
%
% the cellcol part of the divergence operator from cellnodecolrowtbl to
% cellcoltbl is obtained for any u in cellnodecolrowtbl by evaluating the
% expression divcell_T.evalProd(dreduced, u) where dreduced and divcell_T
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

divcell_T = SparseTensor();
divcell_T = divcell_T.setFromTensorProd(dreduced, prod);


%% Construction of transpose operator for matrices at nodes (that are
%% elements of nodecolrowtbl)
%
%  trans_T: nodecolrowtbl -> nodecolrowtbl
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

trans_T = SparseTensor();
trans_T = trans_T.setFromTensorProd(ones(col2row2tbl.num, 1), prod);

%% Construction of nodal average for cellnode tensor
%
% transnodeaverage_T : cellnodecolrowtbl -> nodecolrowtbl
%
% (later this operator is dispatched to cells)
prod = TensorProd();
prod.tbl1 = celltbl;
prod.tbl2 = cellnodecolrowtbl;
prod.reducefds = {'cells'};
prod.prodtbl = nodecolrowtbl;

prod = prod.setup();

nodeaverage_T = SparseTensor();
nodeaverage_T = nodeaverage_T.setFromTensorProd(ones(celltbl.num), prod);

transnodeaverage_T = trans_T*nodeaverage_T;

% we need to dispatch this tensor to cellnodecolrowtbl
% now we have
% transnodeaverage_T : cellnodecolrowtbl -> cellnodecolrowtbl

prod = TensorProd();
prod.tbl1 = celltbl;
prod.tbl2 = nodecolrowtbl;
prod.prodtbl = cellnodecolrowtbl;
prod = prod.setup();

celldispatch_T = SparseTensor();
celldispatch_T = celldispatch_T.setFromTensorProd(ones(celltbl.num), prod);

transnodeaverage_T = celldispatch_T*transnodeaverage_T;

%% Construction of the stiffness operator
%
% C_T : cellnodecolrowtbl -> cellnodecolrowtbl
%
dim = G.griddim;
vdim = dim*(dim + 1)/2;
asymdim = dim*dim - vdim;
voigttbl.voigt = (1 : vdim)';
voigttbl.num = vdim;

switch dim
  case 2
    voigttbl.coldim = [1; 2; 2];
    voigttbl.rowdim = [1; 2; 1];
    
    colrowvoigttbl.coldim = colrowtbl.coldim;
    colrowvoigttbl.rowdim = colrowtbl.rowdim;
    colrowvoigttbl.voigt = [1; 3; 3; 2];
    colrowvoigttbl.num    = colrowtbl.num;

    voigt = [1; 0.5; 0.5; 1];
    
    colrowasymtbl.asym = [1; 1];
    colrowasymtbl.coldim = [1; 2];    
    colrowasymtbl.rowdim = [2; 1];
    colrowasymtbl.num = numel(colrowasymtbl.asym);
    
    asym = [1; -1];
    
  case 3
    error('not implemented yet');
    voigttbl.coldim = [1; 2; 3; 3; 3; 2];
    voigttbl.rowdim = [1; 2; 3; 2; 1; 1];
    voigt = [1; 6; 5; 6; 2; 4; 5; 4; 3];
end

clear voigttbl
voigttbl.voigt = (1 : vdim)';
voigttbl.num = vdim;

voigt2tbl = crossTable(voigttbl, voigttbl, {}, 'crossextend', {{'voigt', ...
                    {'voigt1', 'voigt2'}}});

col2row2voigt2tbl = crossTable(colrowvoigttbl, voigt2tbl, {{'voigt', ...
                    'voigt1'}});
fds = {{'voigt', 'voigt1'}, {'coldim', 'coldim1'}, {'rowdim', 'rowdim1'}};
col2row2voigt2tbl = replacefield(col2row2voigt2tbl, fds);

col2row2voigt2tbl = crossTable(col2row2voigt2tbl, colrowvoigttbl, {{'voigt2', 'voigt'}} );
fds = {{'coldim', 'coldim2'}, {'rowdim', 'rowdim2'}};
col2row2voigt2tbl = replacefield(col2row2voigt2tbl, fds);


col2row2tbl = rmfield(col2row2voigt2tbl, {'voigt1', 'voigt2'});
% note that col2row2voigt2tbl and col2row2tbl correspond to the same indexing

C = [9; 3; 2; 3; 10; 1; 2; 1; 11];

C = tblmap(C, voigt2tbl, col2row2voigt2tbl, {'voigt1', 'voigt2'});

fds = {{'voigt', 'voigt1'}, {'coldim', 'coldim1'}, {'rowdim', 'rowdim1'}};
voigt1 = tblmap(voigt, colrowvoigttbl, col2row2voigt2tbl, fds);
fds = {{'voigt', 'voigt2'}, {'coldim', 'coldim2'}, {'rowdim', 'rowdim2'}};
voigt2 = tblmap(voigt, colrowvoigttbl, col2row2voigt2tbl, fds);

C = voigt1.*C.*voigt2;

% note necessary in fact because col2row2voigt2tbl and col2row2tbl correspond
% to the same indexing
C = tblmap(C, col2row2voigt2tbl, col2row2tbl, {'coldim1', 'rowdim1', 'coldim2', ...
                    'rowdim2'});

colrowtbl = addLocInd(colrowtbl, 'colrowdim');

ind1 = tblmap(colrowtbl.colrowdim, colrowtbl, col2row2tbl, {{'coldim', 'coldim1'}, ...
                    {'rowdim', 'rowdim1'}});
ind2 = tblmap(colrowtbl.colrowdim, colrowtbl, col2row2tbl, {{'coldim', 'coldim2'}, ...
                    {'rowdim', 'rowdim2'}});

n = colrowtbl.num;
Cmat = sparse(ind1, ind2, C, n, n);
Cmat = full(Cmat);

asymtbl.asym = (1 : asymdim)';
asymtbl.num = asymdim;

asym2tbl = crossTable(asymtbl, asymtbl, {}, 'crossextend', {{'asym', ...
                    {'asym1', 'asym2'}}});

col2row2asym2tbl = crossTable(colrowasymtbl, asym2tbl, {{'asym', ...
                    'asym1'}});
fds = {{'asym', 'asym1'}, {'coldim', 'coldim1'}, {'rowdim', 'rowdim1'}};
col2row2asym2tbl = replacefield(col2row2asym2tbl, fds);

col2row2asym2tbl = crossTable(col2row2asym2tbl, colrowasymtbl, {{'asym2', 'asym'}} );
fds = {{'coldim', 'coldim2'}, {'rowdim', 'rowdim2'}};
col2row2asym2tbl = replacefield(col2row2asym2tbl, fds);

Casym = [1];
Casym = tblmap(Casym, asym2tbl, col2row2asym2tbl, {'asym1', 'asym2'});

fds = {{'asym', 'asym1'}, {'coldim', 'coldim1'}, {'rowdim', 'rowdim1'}};
asym1 = tblmap(asym, colrowasymtbl, col2row2asym2tbl, fds);
fds = {{'asym', 'asym2'}, {'coldim', 'coldim2'}, {'rowdim', 'rowdim2'}};
asym2 = tblmap(asym, colrowasymtbl, col2row2asym2tbl, fds);

Casym = asym1.*Casym.*asym2;

fds = {'coldim1', 'rowdim1', 'coldim2', 'rowdim2'};
Casym = tblmap(Casym, col2row2asym2tbl, col2row2tbl, fds);
Casymmat = sparse(ind1, ind2, Casym, n, n);
Casymmat = full(Casymmat);

C = Cmat + Casymmat;

return

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

mat = sparse(voigt1, voigt2, stensind);
mat = diag(diag(mat)) + abs(mat - mat');
[i, j, ind] = find(mat);

stensfulltbl.voigt1 = i;
stensfulltbl.voigt2 = j;
stensfulltbl.stensind = ind;
stensfulltbl.num = numel(i);

crossfds = {{'coldim', {'coldim1', 'coldim2'}}, ...
            {'rowdim', {'rowdim1', 'rowdim2'}}, ...
           };
col2row2tbl = crossTable(colrowtbl, colrowtbl, {}, 'crossextend', crossfds);

fds = {{'coldim1', 'coldim'}, {'rowdim1', 'rowdim'}};
col2row2voigt2tbl = crossTable(col2row2tbl, colrowvoigttbl, fds);
col2row2voigt2tbl = replacefield(col2row2voigt2tbl, {'voigt', 'voigt1'});
fds = {{'coldim2', 'coldim'}, {'rowdim2', 'rowdim'}};
col2row2voigt2tbl = crossTable(col2row2voigt2tbl, colrowvoigttbl, fds);
col2row2voigt2tbl = replacefield(col2row2voigt2tbl, {'voigt', 'voigt2'});

col2row2voigt2tbl = crossTable(col2row2voigt2tbl, stensfulltbl, {'voigt1', ...
                    'voigt2'});

% We choose some Stiffness matrix

switch dimcase
  case 2
    temptbl.voigt1 = [1; 2; 3];
    temptbl.voigt2 = [1; 2; 3];
    temptbl.num = 3; 
    C = tblmap([1; 1; 0.5], temptbl, stenstbl, {'voigt1', 'voigt2'});
  case 3
end

C = tblmap(C, stenstbl, col2row2voigt2tbl, {'stensind'});

prod = TensorProd();
prod.tbl1 = col2row2voigt2tbl;
prod.tbl2 = colrowtbl;
prod.replacefds1 = {{'coldim1', 'coldim'}, {'rowdim1', 'rowdim'}, ...
                    {'voigt1', ''}, {'voigt2', ''}};
prod.replacefds2 = {{'coldim', 'coldim2'}, {'rowdim', 'rowdim2'}};
prod.reducefds = {'coldim2', 'rowdim2'};
prod.prodtbl = colrowtbl;
prod = prod.setup();

Cmat = SparseMatrix();
Cmat = Cmat.setFromTensorProd(C, prod);
Cmat = Cmat.matrix;

col2row2tbl = rmfield(col2row2voigt2tbl, {'voigt1', 'voigt2', 'stensind'});
[cellnodecol2row2tbl, indstruct] = crossTable(cellnodetbl, col2row2tbl, {});
C = tbldispatch2(C, indstruct);

prod = TensorProd();
prod.tbl1 = cellnodecol2row2tbl;
prod.tbl2 = cellnodecolrowtbl;
prod.replacefds1 = {{'coldim1', 'coldim'}, {'rowdim1', 'rowdim'}};
prod.replacefds2 = {{'coldim', 'coldim2'}, {'rowdim', 'rowdim2'}};
prod.mergefds = {'cells', 'nodes'};
prod.reducefds = {'coldim2', 'rowdim2'};
prod.prodtbl = cellnodecolrowtbl;
prod = prod.setup();

C_T = SparseTensor();
C_T = C_T.setFromTensorProd(C, prod);

% Cgradnodeface_T = C_T*gradnodeface_T;
Cgradnodeface_T = gradnodeface_T;
transaverCgradnodeface_T = transnodeaverage_T*Cgradnodeface_T;

combCgradnodeface_T = Cgradnodeface_T + transaverCgradnodeface_T;

% Cgradcell_T = C_T*gradcell_T;
Cgradcell_T = gradcell_T;
transaverCgradcell_T = transnodeaverage_T*Cgradcell_T;

combCgradcell_T = Cgradcell_T + transaverCgradcell_T;

A11 = divnodeface_T*combCgradnodeface_T;
A12 = divnodeface_T*combCgradcell_T; 
A21 = divcell_T*combCgradnodeface_T;
A22 = divcell_T*combCgradcell_T; 


%% Setup for Dirichlet condition
%
% Setup projection tensor P_t : nodefacecoltbl -> nodeintfacecoltbl
%
nodeintfacecoltbl = crossTable(intfacetbl, nodefacecoltbl, {'faces'});
fds = {'nodes', 'faces', 'coldim'};
nodeintfacecoltbl = sortTable(nodeintfacecoltbl, fds);

prod = TensorProd();
prod.tbl1 = intfacetbl;
prod.tbl2 = nodeintfacecoltbl;
prod.mergefds = {'faces'};
prod.prodtbl = nodefacecoltbl;
prod = prod.setup();

P_T = SparseTensor();
P_T = P_T.setFromTensorProd(ones(intfacetbl.num, 1), prod);

P_Tt = P_T.transpose();

A11 = P_Tt*A11*P_T;
A12 = P_Tt*A12;
A21 = A21*P_T;

% convert Aij to matrices;

A11 = A11.getMatrix();
A12 = A12.getMatrix();
A21 = A21.getMatrix();
A22 = A22.getMatrix();

% get the block structure
% We count the number of degrees of freedom that are connected to the same
% node.
[nodes, sz] = rlencode(nodeintfacecoltbl.nodes, 1);

invA11 = bi(A11, sz);

A = A22 - A21*invA11*A12;

%% solve a problem

% We setup a source-term
switch dimcase
  case 2
    indcell = floor(nx/2) + nx*floor((ny - 1)/2);
    force = [0; 1]; % force in upward direction
  case 3
    indcell = floor(nx/2 + ny/2*nx + nz/2*nx*ny);
    force = [0; 0; 1]; % force in upward direction    
end

sourcetbl.cells = indcell;
sourcetbl.num = numel(indcell);

sourcetbl = crossTable(sourcetbl, coltbl, {});

force = tblmap(force, coltbl, sourcetbl, {'coldim'});
force = tblmap(force, sourcetbl, cellcoltbl, {'cells', 'coldim'});

u = A\force;

u = reshape(u, dimcase, [])';

return

close all
figure
plotGrid(G)
plotGrid(G, indcell, 'facecolor', 'red');
title('source cell')

figure
plotCellData(G, u(:, 1));
title('displacement - x direction')
figure
plotCellData(G, u(:, 2));
title('displacement - y direction')
