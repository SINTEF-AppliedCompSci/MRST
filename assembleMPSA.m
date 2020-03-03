function [assembly, tbls] = assembleMPSA(G, prop, loadstruct, eta, tbls, mappings, ...
                                         runcase)

%% Assembly of MPSA-weak
%%
%% Reference paper:
%% Finite volume methods for elasticity with weak symmetry
%% Keilegavlen, Eirik and Nordbotten, Jan Martin
%% International Journal for Numerical Methods in Engineering
%% 2017

    coltbl                = tbls.coltbl;
    celltbl               = tbls.celltbl;
    nodetbl               = tbls.nodetbl;
    cellfacetbl           = tbls.cellfacetbl;
    cellnodetbl           = tbls.cellnodetbl;
    nodefacetbl           = tbls.nodefacetbl;
    cellcoltbl            = tbls.cellcoltbl;
    nodecoltbl            = tbls.nodecoltbl;
    nodefacecoltbl        = tbls.nodefacecoltbl;
    cellnodefacetbl       = tbls.cellnodefacetbl;
    cellnodecoltbl        = tbls.cellnodecoltbl;
    cellnodecolrowtbl     = tbls.cellnodecolrowtbl;
    cellnodefacecoltbl    = tbls.cellnodefacecoltbl;
    cellnodefacecolrowtbl = tbls.cellnodefacecolrowtbl;
    colrowtbl             = tbls.colrowtbl;
    nodecolrowtbl         = tbls.nodecolrowtbl;
    col2row2tbl           = tbls.col2row2tbl;
    cellcol2row2tbl       = tbls.cellcol2row2tbl;

    cell_from_cellnode         = mappings.cell_from_cellnode;
    node_from_cellnode         = mappings.node_from_cellnode;
    cellface_from_cellnodeface = mappings.cellface_from_cellnodeface;
    cellnode_from_cellnodeface = mappings.cellnode_from_cellnodeface;
    nodeface_from_cellnodeface = mappings.nodeface_from_cellnodeface;
    
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
    
    dim = coltbl.num;

    %% Construction of tensor g (as defined in paper eq 4.1.2)
    % shortcuts:
    %

    fno = cellnodefacetbl.get('faces');
    cno = cellnodefacetbl.get('cells');
    nno = cellnodefacetbl.get('nodes');

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
    % prod.eval(g, u) where prod is defined below
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
    % obtained for any u in cellcoltbl by using v = prod.eval(greduced, u)
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

    greduced = - prod.eval(ones(cnfc_num, 1), g);

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
    fno = cellnodefacetbl.get('faces');
    cno = cellnodefacetbl.get('cells');
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
    % expression divnodeface_T.eval(d, u) where d and divnodeface_T are defined
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
        g = [1; 0; 0; 1];
        g = tblmap(g, colrowtbl, cellnodecolrowtbl, {'coldim', 'rowdim'});
        d = divnodeface_T.getMatrix()*g;
    end

    % divcell_T : cellnodecoltbl -> cellcoltbl
    %
    % the cellcol part of the divergence operator from cellnodecolrowtbl to
    % cellcoltbl is obtained for any u in cellnodecolrowtbl by evaluating the
    % expression divcell_T.eval(dreduced, u) where dreduced and divcell_T
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

    symcol2row2tbl.coldim2 = colrowtbl.get('coldim');
    symcol2row2tbl.rowdim2 = colrowtbl.get('rowdim');
    symcol2row2tbl.coldim1 = colrowtbl.get('rowdim');
    symcol2row2tbl.rowdim1 = colrowtbl.get('coldim');
    symcol2row2tbl = IndexTable(symcol2row2tbl);

    prod = TensorProd();
    prod.tbl1 = symcol2row2tbl;
    prod.tbl2 = nodecolrowtbl;
    prod.tbl3 = nodecolrowtbl;
    prod.replacefds1 = {{'coldim1', 'coldim'}, ...
                        {'rowdim1', 'rowdim'}};
    prod.replacefds2 = {{'coldim', 'coldim2'}, ...
                        {'rowdim', 'rowdim2'}};
    prod.reducefds = {'coldim2', 'rowdim2'};

    symnodecol2row2tbl = crossTable(nodetbl, symcol2row2tbl, {});
    nc2r2_num = symnodecol2row2tbl.num; % shortcut

    % note the definition of symcol2row2tbl above
    prod.pivottbl = symnodecol2row2tbl;
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
    trans_T = trans_T.setFromTensorProd(ones(symcol2row2tbl.num, 1), prod);

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
    fixnodetbl = IndexTable(fixnodetbl);

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
        g = [1; 2; 3; 1];
        g = tblmap(g, colrowtbl, cellnodecolrowtbl, {'coldim', 'rowdim'});
        g = transnodeaverage_T.getMatrix()*g;
    end

    %% Construction of the stiffness operator
    %
    % C_T : cellnodecolrowtbl -> cellnodecolrowtbl
    %

    C = setupStiffnessTensor(prop, tbls);

    [cellnodecol2row2tbl, indstruct] = crossTable(cellnodetbl, cellcol2row2tbl, {'cells'});
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

    %% Assembly
    
    Cgradnodeface_T = bcfix_T*C_T*gradnodeface_T;
    transaverCgradnodeface_T = transnodeaverage_T*Cgradnodeface_T;

    combCgradnodeface_T = Cgradnodeface_T + transaverCgradnodeface_T;

    Cgradcell_T = bcfix_T*C_T*gradcell_T;
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
    
    % Uses the block structure for the local reduction
    % We count the number of degrees of freedom that are connected to the same
    % node.
    [nodes, sz] = rlencode(nodefacecoltbl.get('nodes'), 1);
    invA11 = bi(A11, sz);


    % We enforce the boundary conditions as Lagrange multipliers

    D = setupBC(loadstruct, G, tbls);
    
    % the solution is given by the system
    %
    % A = [[A11, A12, -D];
    %      [A21, A22,  0];
    %      [D' , 0  ,  0]];
    %
    % u = [u(nodeface);
    %      u(cell);
    %      lagmult];
    %
    % f = [force(extnodeface);
    %      force(cells);
    %      0];
    %
    % A*u = f
    %
    % By construction of the method, the matrix A11 is block-diagonal. Hence,
    % we invert it directly and reduce to a cell-centered scheme.
    
    matrices = struct('A11', A11, ...
                      'A12', A12, ...
                      'A21', A21, ...
                      'A22', A22, ...
                      'D'  , D  , ...
                      'invA11', invA11);
    
    B11 = A22 - A21*invA11*A12;
    B12 = A21*invA11*D;
    B21 = -D'*invA11*A12;
    B22 = D'*invA11*D;

    force = loadstruct.force;
    
    rhs1 = -A21*invA11*force;
    rhs2 = -D'*invA11*force;

    B = [[B11, B12]; ...
         [B21, B22]];

    rhs = [rhs1; rhs2];

    % setup mapping from nodeface to node

    map = TensorMap();
    map.fromTbl = nodefacecoltbl;
    map.toTbl   = nodecoltbl;
    map.mergefds = {'nodes', 'coldim'};
    map = map.setup();

    coef = map.eval(ones(nodefacecoltbl.num));
    coef = 1./coef;

    prod = TensorProd();
    prod.tbl1 = nodecoltbl;
    prod.tbl2 = nodefacecoltbl;
    prod.tbl3 = nodecoltbl;
    prod.mergefds = {'nodes', 'coldim'};
    prod = prod.setup();


    nodaldisp_T = SparseTensor('matlabsparse', true);
    nodaldisp_T = nodaldisp_T.setFromTensorProd(coef, prod);

    nodaldisp_op = nodaldisp_T.getMatrix();

    assembly = struct('B'           , B       , ...
                      'rhs'         , rhs     , ...
                      'force'       , force     , ...
                      'matrices'    , matrices, ...
                      'nodaldisp_op', nodaldisp_op);
    
    tbls = struct('cellcoltbl', cellcoltbl);
end

