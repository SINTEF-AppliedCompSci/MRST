%% Assembly of MPSA-weak
%%
%% Reference paper:
%% Finite volume methods for elasticity with weak symmetry
%% Keilegavlen, Eirik and Nordbotten, Jan Martin
%% International Journal for Numerical Methods in Engineering
%% 2017

clear all

tic

% load modules
mrstModule add mimetic mpfa incomp

eta = 1/3;

%% Define and process geometry
% Construct a Cartesian grid 
runcases = {'2d-refinement', ...
            '2d-linear'    , ...
            '2d-compaction', ...
            '3d-linear'};

runcase = '2d-linear';

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
    nx = 30; ny = 30;
    G = cartGrid([nx, ny], [1, 1]);
  case '3d-linear'
    nx = 5; ny = 5; nz = 5;
    G = cartGrid([nx, ny, nz]);
end

% G = twister(G, 0.1);
% compute Grid geometry
G = computeGeometry(G);

%% Basic table setup
%
tic

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

% We setup the cell-face-node table, cellnodefacetbl. Each entry determine a
% unique facet in a corner
% We order cellnodeface in cell-node-face order. This is node to optimize
% for-end loop below.
cellnodefacetbl = crossTable(cellfacetbl, nodefacetbl, {'faces'});
cellnodefacetbl = sortTable(cellnodefacetbl, {'cells', 'nodes', 'faces'});

% We setup the cell-node table, cellnodetbl. Each entry determine a unique
% corner
cellnodetbl = projTable(cellnodefacetbl, {'nodes', 'cells'});
cellnodetbl = sortTable(cellnodetbl, {'cells', 'nodes'});

fds = {'cells'};
cell_from_cellnode = getDispatchInd(celltbl, cellnodetbl, fds);
fds = {'nodes'};
node_from_cellnode = getDispatchInd(nodetbl, cellnodetbl, fds);
fds = {'cells', 'faces'};
cellface_from_cellnodeface = getDispatchInd(cellfacetbl, cellnodefacetbl, fds);
fds = {'cells', 'nodes'};
cellnode_from_cellnodeface = getDispatchInd(cellnodetbl, cellnodefacetbl, fds);
fds = {'faces', 'nodes'};
nodeface_from_cellnodeface = getDispatchInd(nodefacetbl, cellnodefacetbl, fds);

cellnodecoltbl    = crossTable(cellnodetbl, coltbl, {});
cellnodecolrowtbl = crossTable(cellnodecoltbl, rowtbl, {});

cellnodefacecoltbl = crossTable(cellnodefacetbl, coltbl, {});
cellnodefacecolrowtbl = crossTable(cellnodefacecoltbl, rowtbl, {});

% fds = {'cells', 'nodes', 'coldim'};
% cellnodefacecoltbl = crossTable(cellnodefacecoltbl, cellnodecoltbl, fds);
% this sorting may be unecessary. We do it to be sure
% fds = {'cells', 'nodes', 'faces', 'coldim', 'cnfind', 'cncind'};
% cellnodefacecoltbl = sortTable(cellnodefacecoltbl, fds);


% some shortcuts
c_num     = celltbl.num;
n_num     = nodetbl.num;
cnf_num   = cellnodefacetbl.num;
cnfc_num  = cellnodefacecoltbl.num;
cn_num    = cellnodetbl.num;
cncr_num  = cellnodecolrowtbl.num;
nf_num    = nodefacetbl.num;
nfc_num   = nodefacecoltbl.num;
cnfcr_num = cellnodefacecolrowtbl.num;
d_num     = coltbl.num;


%% Construction of tensor g (as defined in paper eq 4.1.2)
% shortcuts:
%

fno = cellnodefacetbl.faces;
cno = cellnodefacetbl.cells;
nno = cellnodefacetbl.nodes;

cellFacetVec = G.faces.centroids(fno, :) - G.cells.centroids(cno, :) + ...
    eta*(G.nodes.coords(nno, :) - G.faces.centroids(fno, :));

cellFacetVec = reshape(cellFacetVec', [], 1);

[c, i] = ind2sub([d_num, cnf_num], (1 : cnfc_num)');
ind1 = i;
ind2 = sub2ind([d_num, cn_num], c, cellnode_from_cellnodeface(i));

n = cellnodecoltbl.num; 
assert(n == cellnodefacetbl.num, ['This implementation of mpsaw cannot handle ' ...
                    'this grid']);

A = sparse(ind1, ind2, cellFacetVec, n, n);

opt.invertBlocks = 'mex';
bi = blockInverter(opt);

sz = repmat(coltbl.num, cellnodetbl.num, 1);
invA = bi(A, sz);

[cncind, cnfind, g] = find(invA);
[c, i] = ind2sub([d_num, cn_num], cncind);
ind = sub2ind([d_num, cnf_num], c, cnfind);

g(ind) = g;

%% Construction of the gradient operator
%


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
prod.reducefds   = {'faces'};
prod.mergefds    = {'nodes'};
prod.tbl3 = cellnodecolrowtbl;

[r, c, i] = ind2sub([d_num, d_num, cnf_num], (1 : cnfcr_num)');

prod.dispind1 = sub2ind([d_num, cnf_num], c, i);
prod.dispind2 = sub2ind([d_num, cnf_num], r, nodeface_from_cellnodeface(i));
prod.dispind3 = sub2ind([d_num, d_num, cn_num], r, c, cellnode_from_cellnodeface(i));
prod.issetup = true;

gradnodeface_T = SparseTensor('matlabsparse', true);
gradnodeface_T = gradnodeface_T.setFromTensorProd(g, prod);

% Construction of gradcell_T : cellcoltbl -> cellnodecolrowtbl
%
% The cellcol part of the grad operator from cellcoltbl to cellnodecolrowtbl is
% obtained for any u in cellcoltbl by using v = prod.evalProd(greduced, u)
% where greduced and prod are defined below 
%
prod = TensorProd();
prod.tbl1 = cellnodefacecoltbl;
prod.tbl2 = cellnodefacecoltbl;
prod.tbl3 = cellnodecoltbl;
prod.pivottbl = cellnodefacecoltbl;
prod.reducefds = {'faces'};

prod.dispind1 = (1 : cnfc_num)';
prod.dispind2 = (1 : cnfc_num)';
[c, i] = ind2sub([d_num, cnf_num], (1 : cnfc_num)');
prod.dispind3 = sub2ind([d_num, cn_num], c, cellnode_from_cellnodeface(i));

prod.issetup = true;

greduced = - prod.evalProd(ones(cnfc_num, 1), g);

prod = TensorProd();
prod.tbl1 = cellnodecoltbl;
prod.tbl2 = cellcoltbl;
prod.tbl3 = cellnodecolrowtbl;
prod.replacefds2 = {'coldim', 'rowdim'};
prod.mergefds = {'cells'};

prod.pivottbl = cellnodecolrowtbl;
[r, c, i] = ind2sub([d_num, d_num, cn_num], (1 : cncr_num)');
prod.dispind1 = sub2ind([d_num, cn_num], c, i);
prod.dispind2 = sub2ind([d_num, c_num], r, cell_from_cellnode(i));
prod.dispind3 = (1 : cncr_num);

prod.issetup = true;

gradcell_T = SparseTensor('matlabsparse', true);
gradcell_T = gradcell_T.setFromTensorProd(greduced, prod);

% some test gradnodeface_T and gradcell_T
dotest = false;
if dotest
    fno = nodefacetbl.faces;
    nno = nodefacetbl.nodes;
    facetcent = G.faces.centroids(fno, :) + eta*(G.nodes.coords(nno, :) - ...
                                                 G.faces.centroids(fno, :));
    facetcent = reshape(facetcent', [], 1);

    g1 = gradnodeface_T.getMatrix()*facetcent;

    cellcent = G.cells.centroids(celltbl.cells, :);
    cellcent = reshape(cellcent', [], 1);

    g2 = gradcell_T.getMatrix()*cellcent;

    g = g1 + g2;
    % g should correspond to identity in cellnodecolrowtbl
end


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
prod.tbl3 = nodefacecoltbl;

prod.pivottbl = cellnodefacecolrowtbl;
[r, c, i] = ind2sub([d_num, d_num, cnf_num], (1 : cnfcr_num)');
prod.dispind1 = sub2ind([d_num, cnf_num], r, i);
prod.dispind2 = sub2ind([d_num, d_num, cn_num], c, r, cellnode_from_cellnodeface(i));
prod.dispind3 = sub2ind([d_num, nf_num], c, nodeface_from_cellnodeface(i));

prod.issetup = true;

divnodeface_T = SparseTensor('matlabsparse', true);
divnodeface_T = divnodeface_T.setFromTensorProd(d, prod);


% some test for dinnodeface_T
dotest = false;
if dotest
    % create uniform gradient tensor (take unity)
    assert(dim == 2);
    colrowtbl = crossTable(coltbl, rowtbl, {});
    g = [1; 0; 0; 1];
    g = tblmap(g, colrowtbl, cellnodecolrowtbl, {'coldim', 'rowdim'});
    d = divnodeface_T.getMatrix()*g;
end

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
prod.tbl3 = cellcoltbl;
prod.replacefds1 = {'coldim', 'rowdim'};
prod.replacefds2 = {'coldim', 'rowdim', 'interchange'};
prod.reducefds   = {'rowdim', 'nodes'};
prod.mergefds    = {'cells'};

prod.pivottbl = cellnodecolrowtbl;
[r, c, i] = ind2sub([d_num, d_num, cn_num], (1 : cncr_num)');
prod.dispind1 = sub2ind([d_num, cn_num], r, i);
prod.dispind2 = sub2ind([d_num, d_num, cn_num], c, r, i);
prod.dispind3 = sub2ind([d_num, c_num], c, cell_from_cellnode(i));

prod.issetup = true;

divcell_T = SparseTensor('matlabsparse', true);
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
prod.tbl3 = nodecolrowtbl;
prod.replacefds1 = {{'coldim1', 'coldim'}, ...
                    {'rowdim1', 'rowdim'}};
prod.replacefds2 = {{'coldim', 'coldim2'}, ...
                    {'rowdim', 'rowdim2'}};
prod.reducefds = {'coldim2', 'rowdim2'};

nodecol2row2tbl = crossTable(nodetbl, col2row2tbl, {});
nc2r2_num = nodecol2row2tbl.num; % shortcut

% note the definition of col2row2tbl above
prod.pivottbl = nodecol2row2tbl;
[r, c, i] = ind2sub([d_num, d_num, n_num], (1 : nc2r2_num)');
c2 = c;
r2 = r;
c1 = r;
r1 = c;
prod.dispind1 = sub2ind([d_num, d_num], r, c);
prod.dispind2 = sub2ind([d_num, d_num, n_num], r2, c2, i);
prod.dispind3 = sub2ind([d_num, d_num, n_num], r1, c1, i);

prod.issetup = true;

trans_T = SparseTensor('matlabsparse', true);
trans_T = trans_T.setFromTensorProd(ones(col2row2tbl.num, 1), prod);

%% Construction of nodal average for cellnode tensor
%
% transnodeaverage_T : cellnodecolrowtbl -> nodecolrowtbl
%
% (later this operator is dispatched to cells)
%

% Compute number of cell per node
[~, indstruct] = crossTable(cellnodetbl, nodetbl, {'nodes'});
nnodepercell = tblmap1to2(ones(cellnodetbl.num, 1), indstruct);
coef   = tblmap2to1(1./nnodepercell, indstruct);

% we eliminitate the places (at the boundaries) where the local reconstruction
% is ill-posed: nodes with one cell in 2d (corners of a Cartesian grid) and
% nodes with less the two nodes in 3d (edges of a Cartesian grid);

switch dim
  case 2
    maxnnodepercell = 1;
  case 3
    maxnnodepercell = 2;
end
    
fixnodetbl.nodes = find(nnodepercell <= maxnnodepercell);
fixnodetbl.num = numel(fixnodetbl.nodes);

coef(coef >= 1/maxnnodepercell) = 0;

prod = TensorProd();
prod.tbl1 = cellnodetbl;
prod.tbl2 = cellnodecolrowtbl;
prod.tbl3 = nodecolrowtbl;
prod.reducefds = {'cells'};
prod.mergefds = {'nodes'};

prod.pivottbl = cellnodecolrowtbl;
[r, c, i] = ind2sub([d_num, d_num, cn_num], (1 : cncr_num)');
prod.dispind1 = i;
prod.dispind2 = (1 : cncr_num)';
prod.dispind3 = sub2ind([d_num, d_num, n_num], r, c, node_from_cellnode(i));

prod.issetup = true;

nodeaverage_T = SparseTensor('matlabsparse', true);
nodeaverage_T = nodeaverage_T.setFromTensorProd(coef, prod);

transnodeaverage_T = trans_T*nodeaverage_T;

% we need to dispatch this tensor to cellnodecolrowtbl
% now we have
% transnodeaverage_T : cellnodecolrowtbl -> cellnodecolrowtbl

prod = TensorProd();
prod.tbl1 = celltbl;
prod.tbl2 = nodecolrowtbl;
prod.tbl3 = cellnodecolrowtbl;

prod.pivottbl = cellnodecolrowtbl;
[r, c, i] = ind2sub([d_num, d_num, cn_num], (1 : cncr_num)');
prod.dispind1 = cell_from_cellnode(i);
prod.dispind2 = sub2ind([d_num, d_num, n_num], r, c, node_from_cellnode(i));
prod.dispind3 = (1 : cncr_num)';

prod.issetup = true;

celldispatch_T = SparseTensor('matlabsparse', true);
celldispatch_T = celldispatch_T.setFromTensorProd(ones(celltbl.num), prod);

transnodeaverage_T = celldispatch_T*transnodeaverage_T;

%% we need to multiply by 2 at the place where we discarded the symmetry requirement

fixcellnodecolrowtbl = crossTable(fixnodetbl, cellnodecolrowtbl, {'nodes'});

ind = tblmap(ones(fixnodetbl.num, 1), fixnodetbl, cellnodecolrowtbl, ...
             {'nodes'});

c = ones(cellnodecolrowtbl.num, 1);
c(logical(ind)) = 2;

prod = TensorProd();
prod.tbl1 = cellnodecolrowtbl;
prod.tbl2 = cellnodecolrowtbl;
prod.mergefds = {'cells', 'nodes', 'coldim', 'rowdim'};
prod.tbl3 = cellnodecolrowtbl;

prod.pivottbl = cellnodecolrowtbl;
cncr_num = cellnodecolrowtbl.num; %shortcut
prod.dispind1 = (1 : cncr_num)';
prod.dispind2 = (1 : cncr_num)';
prod.dispind3 = (1 : cncr_num)';

prod.issetup = true;

bcfix_T = SparseTensor('matlabsparse', true);
bcfix_T = bcfix_T.setFromTensorProd(c, prod);

% some test for transnodeaverage_T
dotest = false;
if dotest
    assert(dim == 2);
    colrowtbl = crossTable(coltbl, rowtbl, {});
    g = [1; 2; 3; 1];
    g = tblmap(g, colrowtbl, cellnodecolrowtbl, {'coldim', 'rowdim'});
    g = transnodeaverage_T.getMatrix()*g;
end

%% Construction of the stiffness operator
%
% C_T : cellnodecolrowtbl -> cellnodecolrowtbl
%

switch dim
    
  case 2
    mu = 1;
    lambda = 1;
    Cvoigt = mu*[[2 0 0]; ...
                 [0 2 0]; ...
                 [0 0 1]];
    Z = [0; 0];
    Cvoigt = [[lambda*ones(2, 2), Z]; ...
              [Z', 0] ...
             ] + Cvoigt;
    Casym = mu*1;

    n1 = size(Cvoigt, 1);
    n2 = size(Casym, 1);

    Z = zeros(n1, n2);

    C = [[Cvoigt, Z];...
         [Z', Casym];
        ];

    % convert from colrow to voigt (including asymetric part)
    % follows indexing of colrowtbl
    M = [[1  0 0 0]; ...
         [0  0 0 1]; ...
         [0  1 1 0]; ...
         [0 -1 1 0] ...
        ];

    % We should add an extra multiplication for the coef 2 on diagonal for Voigt.
    V = diag([1; 1; 0.5; 0.5]);

  case 3
    
end

    
C = M'*V*C*M;


fds = {{'rowdim', {'rowdim1', 'rowdim2'}}, ...
       {'coldim', {'coldim1', 'coldim2'}}};
col2row2tbl = crossTable(colrowtbl, colrowtbl, {}, 'crossextend', fds);

C = reshape(C', [], 1);

dotest = false;
if dotest
    prod = TensorProd();
    prod.tbl1 = col2row2tbl;
    prod.tbl2 = colrowtbl;
    prod.replacefds1 = {{'coldim1', 'coldim'}, {'rowdim1', 'rowdim'}};
    prod.replacefds2 = {{'coldim', 'coldim2'}, {'rowdim', 'rowdim2'}};
    prod.reducefds = {'coldim2', 'rowdim2'};
    prod.tbl3 = colrowtbl;
    prod = prod.setup();

    C_T = SparseTensor('matlabsparse', true);
    C_T = C_T.setFromTensorProd(C, prod);
    
    printTensor(C_T); 
end 

[cellnodecol2row2tbl, indstruct] = crossTable(cellnodetbl, col2row2tbl, {});
C = tbldispatch2(C, indstruct);

prod = TensorProd();
prod.tbl1 = cellnodecol2row2tbl;
prod.tbl2 = cellnodecolrowtbl;
prod.replacefds1 = {{'coldim1', 'coldim'}, {'rowdim1', 'rowdim'}};
prod.replacefds2 = {{'coldim', 'coldim2'}, {'rowdim', 'rowdim2'}};
prod.mergefds = {'cells', 'nodes'};
prod.reducefds = {'coldim2', 'rowdim2'};
prod.tbl3 = cellnodecolrowtbl;

prod.pivottbl = cellnodecol2row2tbl;

d = d_num; %shortcut
cnc2r2_num = cellnodecol2row2tbl.num; %shortcut
[r2, c2, r1, c1, i] = ind2sub([d, d, d, d, cn_num], (1 : cnc2r2_num)');
prod.dispind1 = (1 : cnc2r2_num)';
prod.dispind2 = sub2ind([d, d, cn_num], r1, c1, i);
prod.dispind3 = sub2ind([d, d, cn_num], r2, c2, i);

prod.issetup = true;

C_T = SparseTensor('matlabsparse', true);
C_T = C_T.setFromTensorProd(C, prod);

Cgradnodeface_T = bcfix_T*C_T*gradnodeface_T;
% Cgradnodeface_T = bcfix_T*gradnodeface_T;
transaverCgradnodeface_T = transnodeaverage_T*Cgradnodeface_T;

combCgradnodeface_T = Cgradnodeface_T + transaverCgradnodeface_T;

Cgradcell_T = bcfix_T*C_T*gradcell_T;
% Cgradcell_T = gradcell_T;
transaverCgradcell_T = transnodeaverage_T*Cgradcell_T;

combCgradcell_T = Cgradcell_T + transaverCgradcell_T;

A11 = divnodeface_T*combCgradnodeface_T;
A12 = divnodeface_T*combCgradcell_T; 
A21 = divcell_T*combCgradnodeface_T;
A22 = divcell_T*combCgradcell_T; 

dotest = false;
if dotest
    A11mat = A11.getMatrix()
    [nodes, sz] = rlencode(nodefacecoltbl.nodes, 1);
    invA11 = bi(A11mat, sz);
end

A11 = A11.getMatrix();
A12 = A12.getMatrix();
A21 = A21.getMatrix();
A22 = A22.getMatrix();

tbls = struct('nodefacetbl'       , nodefacetbl       , ...
              'nodefacecoltbl'    , nodefacecoltbl    , ...
              'cellnodefacetbl'   , cellnodefacetbl   , ...
              'cellnodefacecoltbl', cellnodefacecoltbl, ...
              'coltbl'            , coltbl);

mappings = struct('nodeface_from_cellnodeface', nodeface_from_cellnodeface);

[D, force] = setupBC(runcase, G, tbls, mappings, 'facetNormals', facetNormals);

% get the block structure
% We count the number of degrees of freedom that are connected to the same
% node.
[nodes, sz] = rlencode(nodefacecoltbl.nodes, 1);
invA11 = bi(A11, sz);

B11 = A22 - A21*invA11*A12;
B12 = A21*invA11*D;
B21 = -D'*invA11*A12;
B22 = D'*invA11*D;

force1 = -A21*invA11*force;
force2 = -D'*invA11*force;

B = [[B11, B12]; ...
     [B21, B22]];

force = [force1; force2];

u = B\force;
n = cellcoltbl.num;
u = u(1 : n); 
u = reshape(u, dim, [])';

%% plotting
% 
toc

close all
figure
plotCellData(G, u(:, 1));
title('displacement - x direction')
colorbar
figure
plotCellData(G, u(:, 2));
title('displacement - y direction')
colorbar