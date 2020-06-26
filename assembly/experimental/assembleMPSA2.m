function assembly = assembleMPSA2(G, prop, loadstruct, eta, tbls, mappings, varargin)
%% Assembly of MPSA-weak
%%
%% Reference paper:
%% Finite volume methods for elasticity with weak symmetry
%% Keilegavlen, Eirik and Nordbotten, Jan Martin
%% International Journal for Numerical Methods in Engineering
%% 2017

    % the solution is given by the system
    %
    % A = [[A11, A12, -D];
    %      [A21, A22,  0];
    %      [D' , 0  ,  0]];
    %
    % u = [u (displacement at nodefacecoltbl);
    %      u (displacement at cellcoltbl);
    %      lagmult (forces in the linear directions at the boundary)];
    %
    % f = [extforce (force at nodefacecoltbl);
    %      force    (volumetric force at cellcoltbl);
    %      bcvals   (for the linear form at the boundary)];
    %
    % A*u = f
    %
    % Note: extforce is sparse and should only give contribution at facets
    % that are at the boundary
    %
    % By construction of the method, the matrix A11 is block-diagonal. Hence,
    % we invert it directly and reduce to a cell-centered scheme.
    
    
    dooptimize = true;
    
    opt = struct('bcetazero'       , true , ...
                 'assemblyMatrices', false, ...
                 'extraoutput'     , false);
    opt = merge_options(opt, varargin{:});
    
    % Recover IndexArrays
    coltbl                = tbls.coltbl;
    celltbl               = tbls.celltbl;
    nodetbl               = tbls.nodetbl;
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
    col12row12tbl           = tbls.col2row2tbl;
    cellcol12row12tbl       = tbls.cellcol2row2tbl;
    cellnodecol12row12tbl   = tbls.cellnodecol2row2tbl;
    
    cell_in_cellnode         = mappings.cell_from_cellnode;
    node_in_cellnode         = mappings.node_from_cellnode;
    cellnode_in_cellnodeface = mappings.cellnode_from_cellnodeface;
    nodeface_in_cellnodeface = mappings.nodeface_from_cellnodeface;
    
    node_in_cellnodeface = node_in_cellnode(cellnode_in_cellnodeface);
    cell_in_cellnodeface = cell_in_cellnode(cellnode_in_cellnodeface);
    
    % Some shortcuts
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

    %% Construction of tensor g (as defined in paper eq 4.1.2), belongs to cellnodefacecoltbl
    g = computeConsistentGradient(G, eta, tbls, mappings, 'bcetazero', opt.bcetazero);
    
    %% Setup of the facet normals, belongs to cellnodefacecoltbl
    facetNormals = computeFacetNormals(G, cellnodefacetbl);

    %% Setup of stiffness tensor, belongs to cellcol12row12tbl
    C = setupStiffnessTensor(prop, tbls);

    %% Compute number of cell per node
    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = nodetbl;
    map.mergefds = {'nodes'};
    map = map.setup();
    
    ncellpernode = map.eval(ones(cellnodetbl.num, 1));
    
    % setup cell12nodeface12tbl
    [cell12nodeface12tbl, indstruct] = crossIndexArray(cellnodefacetbl, cellnodefacetbl, {'nodes'}, 'crossextend', {{'cells', ...
                        {'cells1', 'cells2'}}, {'faces', {'faces1', 'faces2'}}});
    if dooptimize
        cell1nodeface1_in_cell12nodeface12 = indstruct{1}.inds;
        cell2nodeface2_in_cell12nodeface12 = indstruct{2}.inds;
    end
    
    % setup cell12nodeface2tbl
    [cell12nodeface2tbl, indstruct] = projIndexArray(cell12nodeface12tbl, {'cells1', 'cells2', 'nodes', 'faces2'});
    if dooptimize
        cell12nodeface2_in_cell12nodeface12 = indstruct.inds;
        
        ind1 = cell2nodeface2_in_cell12nodeface12;
        ind2 = cell12nodeface2_in_cell12nodeface12;
        ind = zeros(cell12nodeface2tbl.num, 1);
        ind(ind2) = (1 : cell12nodeface12tbl.num)';
        cell2nodeface2_in_cell12nodeface2 = ind1(ind);
        
    end
    
    % setup nodeface12tbl
    [nodeface12tbl, indstruct] = projIndexArray(cell12nodeface12tbl, {'nodes', 'faces1', 'faces2'});
    if dooptimize
        
        nodeface12_in_cell12nodeface12 = indstruct.inds;
        
        % setup mapping nodeface1 in nodeface12
        map = TensorMap();
        map.fromTbl = nodefacetbl;
        map.toTbl = nodeface12tbl;
        map.replaceFromTblfds = {{'faces', 'faces1'}};
        map.mergefds = {'nodes', 'faces1'};
        nodeface1_in_nodeface12 = map.getDispatchInd();
        
        % setup mapping nodeface2 in nodeface12
        map = TensorMap();
        map.fromTbl = nodefacetbl;
        map.toTbl = nodeface12tbl;
        map.replaceFromTblfds = {{'faces', 'faces1'}};
        map.mergefds = {'nodes', 'faces1'};
        nodeface1_in_nodeface12 = map.getDispatchInd();
        
    end
    
    % setup cell1nodeface2tbl
    [cell1nodeface2tbl, indstruct] = projIndexArray(cell12nodeface12tbl, {'cells1', 'nodes', 'faces2'});
    if dooptimize
        
        cell1nodeface2_in_cell12nodeface12 = indstruct.inds;
        
        % setup mapping cell1 in cell1nodeface2
        map = TensorMap();
        map.fromTbl = celltbl;
        map.toTbl = cell1nodeface2tbl;
        map.replaceFromTblfds = {{'cells', 'cells1'}};
        map.mergefds = {'cells1'};
        cell1_in_cell1nodeface2 = map.getDispatchInd();
        
        % setup mapping nodeface2 in cell1nodeface2
        map = TensorMap();
        map.fromTbl = nodefacetbl;
        map.toTbl = cell1nodeface2tbl;
        map.replaceFromTblfds = {{'faces', 'faces2'}};
        map.mergefds = {'nodes', 'faces2'};
        nodeface2_in_cell1nodeface2 = map.getDispatchInd();
        
    end
    
    % setup cell2nodeface1tbl
    [cell2nodeface1tbl, indstruct] = projIndexArray(cell12nodeface12tbl, {'cells2', 'nodes', 'faces1'});
    if dooptimize
        
        cell2nodeface1_in_cell12nodeface12 = indstruct.inds;

        % setup mapping cell2 in cell2nodeface1
        map = TensorMap();
        map.fromTbl = celltbl;
        map.toTbl = cell2nodeface1tbl;
        map.replaceFromTblfds = {{'cells', 'cells2'}};
        map.mergefds = {'cells2'};
        cell2_in_cell2nodeface1 = map.getDispatchInd();
        
        % setup mapping nodeface1 in cell2nodeface1
        map = TensorMap();
        map.fromTbl = nodefacetbl;
        map.toTbl = cell2nodeface1tbl;
        map.replaceFromTblfds = {{'faces', 'faces1'}};
        map.mergefds = {'nodes', 'faces1'};
        nodeface1_in_cell2nodeface1 = map.getDispatchInd();
        
    end
    
    % setup cell12tbl
    [cell12tbl, indstruct] = projIndexArray(cell12nodeface12tbl, {'cells1', 'cells2'});
    if dooptimize
        
        cell12_in_cell12nodeface12 = indstruct.inds;

        % setup mapping cell1 in cell12
        map = TensorMap();
        map.fromTbl = celltbl;
        map.toTbl = cell12tbl;
        map.replaceFromTblfds = {{'cells', 'cells1'}};
        map.mergefds = {'cells1'};
        cell1_in_cell12 = map.getDispatchInd();
        
        % setup mapping cell2 in cell12
        map = TensorMap();
        map.fromTbl = celltbl;
        map.toTbl = cell12tbl;
        map.replaceFromTblfds = {{'cells', 'cells2'}};
        map.mergefds = {'cells2'};
        cell2_in_cell12 = map.getDispatchInd();
        
    end
    
    % Extend tables in the spatial dimensions.
    colrow12tbl = crossIndexArray(colrowtbl, colrowtbl, {'coldim'}, 'crossextend', {{'rowdim', {'rowdim1', 'rowdim2'}}});
    colrow12tbl = sortIndexArray(colrow12tbl, {'coldim', 'rowdim1', 'rowdim2'});
    
    cellnodefacecolrow12tbl     = crossIndexArray(cellnodefacetbl    , colrow12tbl, {}, 'optpureproduct', true);
    nodefacecolrow12tbl         = crossIndexArray(nodefacetbl        , colrow12tbl, {}, 'optpureproduct', true);
    cell12nodeface2colrow12tbl  = crossIndexArray(cell12nodeface2tbl , colrow12tbl, {}, 'optpureproduct', true);
    cell12nodeface12colrow12tbl = crossIndexArray(cell12nodeface12tbl, colrow12tbl, {}, 'optpureproduct', true);
    
    col12tbl = crossIndexArray(coltbl, coltbl, {}, 'crossextend', {{'coldim', {'coldim1', 'coldim2'}}});
    col12tbl = sortIndexArray(col12tbl, {'coldim1', 'coldim2'});
    
    cell12nodeface12col12tbl = crossIndexArray(cell12nodeface12tbl, col12tbl, {}, 'optpureproduct', true);
    nodeface12col12tbl       = crossIndexArray(nodeface12tbl      , col12tbl, {}, 'optpureproduct', true);
    cell1nodeface2col12tbl   = crossIndexArray(cell1nodeface2tbl  , col12tbl, {}, 'optpureproduct', true);
    cell2nodeface1col12tbl   = crossIndexArray(cell2nodeface1tbl  , col12tbl, {}, 'optpureproduct', true);
    cell12col12tbl           = crossIndexArray(cell12tbl          , col12tbl, {}, 'optpureproduct', true);
    
    cellnodefacecol12row12tbl = crossIndexArray(cellnodefacetbl, col12row12tbl, {}, 'optpureproduct', true);
    
    % shorcuts
    d = d_num;
    c12nf12_num       = cell12nodeface12tbl.num;
    c12nf12cr12_num   = cell12nodeface12colrow12tbl.num;
    c12nf2_num        = cell12nodeface2tbl.num;
    c12nf2cr12_num    = cell12nodeface2colrow12tbl.num;
    cnfc12r12_num     = cellnodefacecol12row12tbl.num;
    cnfcr12_num       = cellnodefacecolrow12tbl.num;
    
    %% We include the fix at the boundary when the symmetry condition cannot be imposed.
    
    switch dim
      case 2
        minncellpernode = 1;
      case 3
        minncellpernode = 2;
    end
    
    bcfix = ones(nodetbl.num, 1);
    bcfix(ncellpernode <= minncellpernode) = 2;
    
    prod = TensorProd();
    prod.tbl1 = nodetbl;
    prod.tbl2 = cellnodefacecoltbl;
    prod.tbl3 = cellnodefacecoltbl;
    prod.mergefds = {'nodes'};
    
    if dooptimize
        prod.pivottbl = cellnodefacecoltbl;
        [c, i] = ind2sub([d_num, cnf_num], (1 : cnfc_num)');
        prod.dispind1 = node_in_cellnodeface(i);
        prod.dispind2 = (1 : cnfc_num)';
        prod.dispind3 = (1 : cnfc_num)';
        prod.issetup = true;
        assert(prod.checkSetup(), 'problem');
    else
        prod = prod.setup();
    end
    
    gbcfix1 = prod.eval(bcfix, g);

    bcfix = ones(nodetbl.num, 1);
    bcfix(ncellpernode <= minncellpernode) = 0;
    
    gbcfix2 = prod.eval(bcfix, g);

    
    prod = TensorProd();
    prod.tbl1 = cellcol12row12tbl;
    prod.tbl2 = cellnodefacecoltbl;
    prod.tbl3 = cellnodefacecolrow12tbl;
    prod.replacefds1 = {{'coldim1', 'coldim'}};
    prod.replacefds2 = {{'coldim', 'coldim2'}};
    prod.mergefds = {'cells'};
    prod.reducefds = {'coldim2'};
    if dooptimize
        prod.pivottbl = cellnodefacecol12row12tbl;
        [r2, c2, r1, c1, i] = ind2sub([d, d, d, d, cnf_num], (1 : cnfc12r12_num)');
        prod.dispind1 = sub2ind([d, d, d, d, c_num], r2, c2, r1, c1, cell_in_cellnodeface(i));
        prod.dispind2 = sub2ind([d, cnf_num], c2, i);
        prod.dispind3 = sub2ind([d, d, d, cnf_num], r2, r1, c1, i);
        prod.issetup = true;
        assert(prod.checkSetup({{'coldim1', 'coldim'}}), 'problem');
    else
        prod = prod.setup();
    end
    
    Cg = prod.eval(C, gbcfix1);
    
    CAverg = prod.eval(C, gbcfix2);
    
    % We take average over the neighoring cells
    prod = TensorProd();
    prod.tbl1 = nodetbl;
    prod.tbl2 = cellnodefacecolrow12tbl;
    prod.tbl3 = cellnodefacecolrow12tbl;
    prod.mergefds = {'nodes'};
    if dooptimize
        prod.pivottbl = cellnodefacecolrow12tbl;
        [r2, r1, c1, i] = ind2sub([d, d, d, cnf_num], (1 : cnfcr12_num)');
        prod.dispind1 = node_in_cellnodeface(i);
        prod.dispind2 = (1 : cnfcr12_num)';
        prod.dispind3 = (1 : cnfcr12_num)';
        prod.issetup = true;
        assert(prod.checkSetup(), 'problem');
    else
        prod = prod.setup();
    end

    CAverg = prod.eval(1./ncellpernode, CAverg);

    %% We take the transpose of CAverg
    
    map = TensorMap();
    map.fromTbl = cellnodefacecolrow12tbl;
    map.toTbl = cellnodefacecolrow12tbl;
    map.replaceFromTblfds = {{'coldim', 'rowdim1', 'interchange'}};
    map.mergefds = {'cells', 'nodes', 'faces', 'coldim', 'rowdim1', 'rowdim2'};
    if dooptimize
        map.pivottbl = cellnodefacecolrow12tbl;
        [r2, r1, c1, i] = ind2sub([d, d, d, cnf_num], (1 : cnfcr12_num)');
        map.dispind1 = sub2ind([d, d, d, cnf_num], r2, c1, r1, i);
        map.dispind2 = (1 : cnfcr12_num)';
        map.issetup = true;
        assert(map.checkSetup(), 'problem');
    else
        map = map.setup();
    end
    
    CAverg = map.eval(CAverg);
    
    % We map Cg and CAverg in cell12nodeface2colrow12tbl
    
    map = TensorMap();
    map.fromTbl = cellnodefacecolrow12tbl;
    map.toTbl = cell12nodeface2colrow12tbl;
    map.replaceFromTblfds = {{'cells', 'cells2'}, {'faces', 'faces2'}};
    map.mergefds = {'cells2', 'nodes', 'faces2', 'coldim', 'rowdim1', 'rowdim2'};
    if dooptimize
        map.pivottbl = cell12nodeface2colrow12tbl;
        [r2, r1, c1, i] = ind2sub([d, d, d, c12nf2_num], (1 : c12nf2cr12_num)');
        map.dispind1 = sub2ind([d, d, d, cnf_num], r2, r1, c1, cell2nodeface2_in_cell12nodeface2(i));
        map.dispind2 = (1 : c12nf2cr12_num)';
        map.issetup = true;
        assert(map.checkSetup(), 'problem');
    else
        map = map.setup();
    end
    
    CAverg = map.eval(CAverg);
    
    cells1 = cell12nodeface2colrow12tbl.get('cells1');
    cells2 = cell12nodeface2colrow12tbl.get('cells2');
    injcoef = zeros(cell12nodeface2colrow12tbl.num, 1);
    injcoef(cells1 == cells2) = 1;
    
    prod = TensorProd();
    prod.tbl1 = cell12nodeface2colrow12tbl;
    prod.tbl2 = cellnodefacecolrow12tbl;
    prod.tbl3 = cell12nodeface2colrow12tbl;
    prod.replacefds2 = {{'cells', 'cells2'}, {'faces', 'faces2'}};
    prod.mergefds = {'cells2', 'nodes', 'faces2', 'coldim', 'rowdim1', 'rowdim2'};
    if dooptimize
        prod.pivottbl = cell12nodeface2colrow12tbl;
        [r2, r1, c1, i] = ind2sub([d, d, d, c12nf2_num], (1 : c12nf2cr12_num)');
        prod.dispind1 = (1 : c12nf2cr12_num)';
        prod.dispind2 = sub2ind([d, d, d, cnf_num], r2, r1, c1, cell2nodeface2_in_cell12nodeface2(i));;
        prod.dispind3 = (1 : c12nf2cr12_num)';
        prod.issetup = true;
        assert(prod.checkSetup(), 'problem');
    else
        prod = prod.setup();
    end
    
    
    Cg = prod.eval(injcoef, Cg);

    Cg = 0.5*(Cg + CAverg);

    %% We multiply Cg with facetNormals
    prod = TensorProd();
    prod.tbl1 = cellnodefacecoltbl;
    prod.tbl2 = cell12nodeface2colrow12tbl;
    prod.tbl3 = cell12nodeface12col12tbl;
    prod.replacefds1 = {{'cells', 'cells1'}, {'faces', 'faces1'}};
    prod.replacefds2 = {{'rowdim1', 'coldim1'}, {'rowdim2', 'coldim2'}};
    prod.mergefds = {'cells1', 'nodes'};
    prod.reducefds = {'coldim'};
    if dooptimize
        prod.pivottbl = cell12nodeface12colrow12tbl;
        [r2, r1, c1, i] = ind2sub([d, d, d, c12nf12_num], (1 : c12nf12cr12_num)');
        prod.dispind1 = sub2ind([d, cnf_num], c1, cell1nodeface1_in_cell12nodeface12(i));
        prod.dispind2 = sub2ind([d, d, d, c12nf2_num], r2, r1, c1, cell12nodeface2_in_cell12nodeface12(i));
        prod.dispind3 = sub2ind([d, d, c12nf12_num], r2, r1, i);
        prod.issetup = true;
        assert(prod.checkSetup({{'rowdim1', 'coldim1'}, {'rowdim2', 'coldim2'}}), 'problem');
    else
        prod = prod.setup();
    end
    nCg = prod.eval(facetNormals, Cg);

    dooptimize = false;
    
    %% We setup A11
    
    map = TensorMap();
    map.fromTbl = cell12nodeface12col12tbl;
    map.toTbl = nodeface12col12tbl;
    map.mergefds = {'nodes', 'faces1', 'faces2', 'coldim1', 'coldim2'};
    if dooptimize
        % prod.pivottbl = cell12nodeface12col12tbl;
    else
        map = map.setup();
    end
    
    A11 = map.eval(nCg);
    
    prod = TensorProd();
    prod.tbl1 = nodeface12col12tbl;
    prod.tbl2 = nodefacecoltbl;
    prod.tbl3 = nodefacecoltbl;
    prod.replacefds1 = {{'faces1', 'faces'}, {'coldim1', 'coldim'}};
    prod.replacefds2 = {{'faces', 'faces2'}, {'coldim', 'coldim2'}};
    prod.mergefds = {'nodes'};
    prod.reducefds = {'faces2', 'coldim2'};
    if dooptimize
        % prod.pivottbl = nodeface12col12tbl;
    else
        prod = prod.setup();
    end
    

    
    A11_T = SparseTensor();
    A11_T = A11_T.setFromTensorProd(A11, prod);
    A11 = A11_T.getMatrix();
    
    %% We setup A12
    
    map = TensorMap();
    map.fromTbl = cell12nodeface12col12tbl;
    map.toTbl = cell2nodeface1col12tbl;
    map.mergefds = {'cells2', 'nodes', 'faces1', 'coldim1', 'coldim2'};
    if dooptimize
        % prod.pivottbl = cell12nodeface12col12tbl;
    else
        map = map.setup();
    end
    

    
    % note the minus sign
    A12 = map.eval( - nCg);
    
    prod = TensorProd();
    prod.tbl1 = cell2nodeface1col12tbl;
    prod.tbl2 = cellcoltbl;
    prod.tbl3 = nodefacecoltbl;
    prod.replacefds1 = {{'coldim1', 'coldim'}, {'cells2', 'cells'}, {'faces1', 'faces'}};
    prod.replacefds2 = {{'coldim', 'coldim2'}};
    prod.reducefds = {'cells', 'coldim2'};
    if dooptimize
        % prod.pivottbl = cell2nodeface1col12tbl;
    else
        prod = prod.setup();
    end
    

    
    A12_T = SparseTensor();
    A12_T = A12_T.setFromTensorProd(A12, prod);
    A12 = A12_T.getMatrix();
    
    %% We setup A21
    
    map = TensorMap();
    map.fromTbl = cell12nodeface12col12tbl;
    map.toTbl = cell1nodeface2col12tbl;
    map.mergefds = {'cells1', 'nodes', 'faces2', 'coldim1', 'coldim2'};
    if dooptimize
        % prod.pivottbl = cell12nodeface12col12tbl;
    else
        map = map.setup();
    end
    

        
    % note the minus sign
    A21 = map.eval( - nCg);
    
    prod = TensorProd();
    prod.tbl1 = cell1nodeface2col12tbl;
    prod.tbl2 = nodefacecoltbl;
    prod.tbl3 = cellcoltbl;
    prod.replacefds1 = {{'coldim1', 'coldim'}, {'cells1', 'cells'}, {'faces2', 'faces'}};
    prod.replacefds2 = {{'coldim', 'coldim2'}};
    prod.reducefds = {'nodes', 'faces', 'coldim2'};
    if dooptimize
        % prod.pivottbl = cell1nodeface2col12tbl;
    else
        prod = prod.setup();
    end
    

    
    A21_T = SparseTensor();
    A21_T = A21_T.setFromTensorProd(A21, prod);
    A21 = A21_T.getMatrix();
    
    %% We setup A22
    
    map = TensorMap();
    map.fromTbl = cell12nodeface12col12tbl;
    map.toTbl = cell12col12tbl;
    map.mergefds = {'cells1', 'cells2', 'coldim1', 'coldim2'};
    if dooptimize
        % prod.pivottbl = cell12nodeface12col12tbl;
    else
        map = map.setup();
    end
    

    
    A22 = map.eval(nCg);
    
    prod = TensorProd();
    prod.tbl1 = cell12col12tbl;
    prod.tbl2 = cellcoltbl;
    prod.tbl3 = cellcoltbl;
    prod.replacefds1 = {{'coldim1', 'coldim'}, {'cells1', 'cells'}};
    prod.replacefds2 = {{'coldim', 'coldim2'}, {'cells', 'cells2'}};
    prod.reducefds = {'coldim2', 'cells2'};
    if dooptimize
        % prod.pivottbl = cell12col12tbl;
    else
        prod = prod.setup();
    end
    

    
    A22_T = SparseTensor();
    A22_T = A22_T.setFromTensorProd(A22, prod);
    A22 = A22_T.getMatrix();
    
    % Uses the block structure for the local reduction. We count the number of degrees of freedom that are connected to the same node.
    [nodes, sz] = rlencode(nodefacecoltbl.get('nodes'), 1);
    opt.invertBlocks = 'mex';
    bi = blockInverter(opt);
    invA11 = bi(A11, sz);

    % We enforce the boundary conditions as Lagrange multipliers

    bc = loadstruct.bc;
    if ~isfield(bc, 'bcnodefacetbl')
        bc = setupFaceBC(bc, G, tbls);
    end
    [D, bcvals] = setupNodeFaceBc(bc, G, tbls);
    
    extforce = loadstruct.extforce;
    force = loadstruct.force;

    fullrhs{1} = extforce;
    fullrhs{2} = force;
    fullrhs{3} = bcvals;
    
    matrices = struct('A11', A11, ...
                      'A12', A12, ...
                      'A21', A21, ...
                      'A22', A22, ...
                      'D'  , D  , ...
                      'invA11', invA11);
    matrices.fullrhs = fullrhs;
    
    % We reduced the system (shur complement) using invA11
    % We obtain system of the form
    %
    % B*u = rhs
    %
    % where
    %
    % B = [[B11, B12];
    %      [B21, B22]];
    %
    % u = [u (displacement at cellcoltbl);
    %      lagmult];
    %
    % rhs = [-A21*invA11*extforce;  +  [force;
    %        -D'*invA11*extforce  ]     bcvals]
    
    B11 = A22 - A21*invA11*A12;
    B12 = A21*invA11*D;
    B21 = -D'*invA11*A12;
    B22 = D'*invA11*D;

    B = [[B11, B12]; ...
         [B21, B22]];
    
    rhs{1} = -A21*invA11*extforce + force; 
    rhs{2} = -D'*invA11*extforce + bcvals;
    
    rhs = vertcat(rhs{:});

    % Assembly of operator to compute u_{nodefacecoltbl} from solution of the system
    % (which consists of the concatenation of u_{cellcol} and lagmult) and
    % extforce which is a force in nodefacecoltbl
    %
    % We have  u_{nodefacecoltbl} = R1*sol + R2*extforce
    %
    % where R1 = invA11*[-A12, D] and R2 = invA11
    
    R1 = invA11*[-A12, D];
    R2 = invA11;
    
    % The divergence operator (integrated over the volume)
    % is given by 
    %
    %  div[c] = sum (m[f,s] u_[f,n,i] n[c,f,i])
    %
    % where u:solution, n:normal, m:area
    % indices : c:cell, f:face, n:node.
    
    % The facetNormals are already weighted with respect to area
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacecoltbl;
    prod.tbl2 = nodefacecoltbl;
    prod.tbl3 = celltbl;
    prod.reducefds = {'faces', 'nodes', 'coldim'};
    % prod = prod.setup();
    
    prod.pivottbl = cellnodefacecoltbl;
    prod.dispind1 = (1 : cnfc_num)';
    [c, i] = ind2sub([d_num, cnf_num], (1 : cnfc_num)');
    prod.dispind2 = sub2ind([d_num, nf_num], c, nodeface_in_cellnodeface(i));
    prod.dispind3 = cell_in_cellnode(cellnode_in_cellnodeface(i));
    prod.issetup = true;
    assert(prod.checkSetup(), 'problem');
    
    div_T = SparseTensor;
    div_T = div_T.setFromTensorProd(facetNormals, prod);
    div = div_T.getMatrix();
    
    
    assembly = struct('B'       , B       , ...
                      'rhs'     , rhs     , ...
                      'g'       , g       , ...
                      'extforce', extforce, ...
                      'R1'      , R1      , ...
                      'R2'      , R2);
    
    if opt.assemblyMatrices
        assembly.matrices = matrices;
    end
    
    if opt.extraoutput
        assembly.divop = @(sol) mpsaDivOperator(sol, extforce, R1, R2, div);
    end
    
end




