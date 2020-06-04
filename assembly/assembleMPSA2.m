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
    col2row2tbl           = tbls.col2row2tbl;
    cellcol2row2tbl       = tbls.cellcol2row2tbl;
    cellnodecol2row2tbl   = tbls.cellnodecol2row2tbl;
    
    cell_from_cellnode         = mappings.cell_from_cellnode;
    node_from_cellnode         = mappings.node_from_cellnode;
    cellnode_from_cellnodeface = mappings.cellnode_from_cellnodeface;
    nodeface_from_cellnodeface = mappings.nodeface_from_cellnodeface;
    
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

    %% Setup of stiffness tensor, belongs to cellcol2row2tbl
    C = setupStiffnessTensor(prop, tbls);

    %% Compute number of cell per node
    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = nodetbl;
    map.mergefds = {'nodes'};
    map = map.setup();
    
    ncellpernode = map.eval(ones(cellnodetbl.num, 1));
    
    colrow2tbl = crossIndexArray(colrowtbl, colrowtbl, {'coldim'}, 'crossextend', {{'rowdim', {'rowdim1', 'rowdim2'}}});
    cellnodefacecolrow2tbl = crossIndexArray(cellnodefacetbl, colrow2tbl, {});
    nodefacecolrow2tbl = crossIndexArray(nodefacetbl, colrow2tbl, {});
    
    % nodefacecolrowtbl = crossIndexArray(nodecolrowtbl, nodefacetbl, {'nodes'});
    % nodefacecolrow2tbl = crossIndexArray(nodefacecolrowtbl, rowtbl, {}, 'crossextend', {{'rowdim', {'rowdim1', 'rowdim2'}}});
    
    % table for dispatching Caverg
    % cellnodefacecolrow2tbl2 = crossIndexArray(nodefacecolrow2tbl, cellnodetbl, {'nodes'});
    
    % Large sparsity than cellnodefacetbl ! 
    cell2nodefacetbl = crossIndexArray(cellnodetbl, cellnodefacetbl, {'nodes'}, 'crossextend', {{'cells', {'cells1', 'cells2'}}});
    cell2nodefacecolrow2tbl = crossIndexArray(cell2nodefacetbl, colrow2tbl, {});
    
    % cellnodefacetbl2 = crossIndexArray(nodefacetbl, cellnodetbl, {'nodes'});
    % colrow2tbl = crossIndexArray(colrowtbl, colrowtbl, {'coldim'}, 'crossextend', {{'rowdim', {'rowdim1', 'rowdim2'}}});
    % cellnodefacecolrow2tbl2 = crossIndexArray(cellnodefacetbl2, colrow2tbl, {});    
    
    % cellnodeface2tbl = crossIndexArray(cellnodefacetbl2, cellnodefacetbl2, {'cells', 'nodes'}, 'crossextend', { {'faces', ...
                        % {'faces1', 'faces2'}}});
    % col2tbl = crossIndexArray(coltbl, coltbl, {}, 'crossextend', {{'coldim', {'coldim1', 'coldim2'}}});
    % cellnodeface2col2tbl = crossIndexArray(cellnodeface2tbl, col2tbl, {});
    
    % nodeface2col2tbl    = projIndexArray(cellnodeface2col2tbl, {'nodes', 'faces1', 'faces2', 'coldim1', 'coldim2'});
    % cellnodefacecol2tbl = projIndexArray(cellnodeface2col2tbl, {'cells', 'nodes', 'faces1', 'coldim1', 'coldim2'});
    % cellnodefacecol2tbl = replacefield(cellnodefacecol2tbl, {'faces1', 'faces'});
    % cellcol2tbl         = projIndexArray(cellnodeface2col2tbl, {'cells', 'coldim1', 'coldim2'});
    
    prod = TensorProd();
    prod.tbl1 = cellcol2row2tbl;
    prod.tbl2 = cellnodefacecoltbl;
    prod.tbl3 = cellnodefacecolrow2tbl;
    prod.replacefds1 = {{'coldim1', 'coldim'}};
    prod.replacefds2 = {{'coldim', 'coldim2'}};
    prod.mergefds = {'cells'};
    prod.reducefds = {'coldim2'};
    prod = prod.setup();
    
    Cg = prod.eval(C, g);
    
    % We take average over Cg and and obtain CAverg
    prod = TensorProd();
    prod.tbl1 = nodetbl;
    prod.tbl2 = cellnodefacecolrow2tbl;
    prod.tbl3 = cellnodefacecolrow2tbl;
    prod.mergefds = {'nodes'};
    prod = prod.setup();
    
    CAverg = prod.eval(1./ncellpernode, Cg);

    %% We take the transpose of CAverg
    
    map = TensorMap();
    map.fromTbl = cellnodefacecolrow2tbl;
    map.toTbl = cellnodefacecolrow2tbl;
    map.replaceFromTblfds = {{'coldim', 'rowdim1', 'interchange'}};
    map.mergefds = {'cells', 'nodes', 'faces', 'coldim', 'rowdim1', 'rowdim2'};
    map = map.setup();
    
    CAverg = map.eval(CAverg);
    
    %% We include the fix at the boundary when the symmetry condition cannot be imposed.
    
    dobcfix = true;
    if dobcfix
        switch dim
          case 2
            minncellpernode = 1;
          case 3
            minncellpernode = 2;
        end
        
        bcfix = ones(nodetbl.num, 1);
        bcfix(ncellpernode <= minncellpernode) = 0;
        
        prod = TensorProd();
        prod.tbl1 = nodetbl;
        prod.tbl2 = cellnodefacecolrow2tbl;
        prod.tbl3 = cellnodefacecolrow2tbl;
        prod.mergefds = {'nodes'};
        prod = prod.setup();
        
        CAverg = prod.eval(bcfix, CAverg);
        
        bcfix = ones(nodetbl.num, 1);
        bcfix(ncellpernode <= minncellpernode) = 2;
        
        prod = TensorProd();
        prod.tbl1 = nodetbl;
        prod.tbl2 = cellnodefacecolrow2tbl;
        prod.tbl3 = cellnodefacecolrow2tbl;
        prod.mergefds = {'nodes'};
        prod = prod.setup();
        
        Cg = prod.eval(bcfix, Cg);

    end
    
    % We map Cg and CAverg in cell2nodefacecolrow2tbl
    
    map = TensorMap();
    map.fromTbl = cellnodefacecolrow2tbl;
    map.toTbl = cell2nodefacecolrow2tbl;
    map.replaceFromTblfds = {{'cells', 'cells2'}};
    map.mergefds = {'cells2', 'nodes', 'faces', 'coldim', 'rowdim1', 'rowdim2'};
    map = map.setup();
    
    CAverg = map.eval(CAverg);
    
    cells1 = cell2nodefacecolrow2tbl.get('cells1');
    cells2 = cell2nodefacecolrow2tbl.get('cells2');
    injcoef = zeros(cell2nodefacecolrow2tbl.num, 1);
    injcoef(cells1 == cells2) = 1;
    
    prod = TensorProd();
    prod.tbl1 = cell2nodefacecolrow2tbl;
    prod.tbl2 = cellnodefacecolrow2tbl;
    prod.tbl3 = cell2nodefacecolrow2tbl;
    prod.replacefds2 = {{'cells', 'cells1'}};
    prod.mergefds = {'cells1', 'nodes', 'faces', 'coldim', 'rowdim1', 'rowdim2'};
    prod = prod.setup();
    
    Cg = prod.eval(injcoef, Cg);

    Cg = 0.5*(Cg + CAverg);
    % Cg = CAverg;
    
    debug = true;

    if debug
        
        % assemble mappin Cg : nodefacecol -> cellnodecolrow
        % prod = TensorProd();
        % prod.tbl1 = cellnodefacecol2rowtbl;
        % prod.tbl2 = nodefacecoltbl;
        % prod.tbl3 = cellnodecolrowtbl;
        % prod.replacefds1 = {{'coldim1', 'coldim'}};
        % prod.replacefds2 = {{'coldim', 'coldim2'}};
        % prod.reducefds = {'coldim2'};
        % prod.mergefds = {'faces', 'nodes'};
        % prod = prod.setup();
        
        prod = TensorProd();
        prod.tbl1 = cell2nodefacecolrow2tbl;
        prod.tbl2 = nodefacecoltbl;
        prod.tbl3 = cellnodecolrowtbl;
        prod.replacefds1 = {{'cells1', 'cells'}, {'rowdim1', 'rowdim'}};
        prod.replacefds2 = {{'coldim', 'rowdim2'}};
        prod.reducefds = {'faces', 'rowdim2'};
        prod.mergefds = {'nodes'};
        prod = prod.setup();

        C1_T = SparseTensor();
        C1_T = C1_T.setFromTensorProd(Cg, prod);
        C1 = C1_T.getMatrix();

        prod = TensorProd();
        prod.tbl1 = cell2nodefacecolrow2tbl;
        prod.tbl2 = cellcoltbl;
        prod.tbl3 = cellnodecolrowtbl;
        prod.replacefds1 = {{'cells1', 'cells'}, {'rowdim1', 'rowdim'}};
        prod.replacefds2 = {{'cells', 'cells2'}, {'coldim', 'rowdim2'}};
        prod.reducefds = {'cells2', 'rowdim2'};
        prod = prod.setup();

        C2_T = SparseTensor();
        C2_T = C2_T.setFromTensorProd(Cg, prod);
        C2 = C2_T.getMatrix();
        
    end
    
    return
    
    %% We multiply Cg with facetNormals
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacecoltbl;
    prod.tbl2 = cellnodefacecolrow2tbl2;
    prod.tbl3 = cellnodeface2col2tbl;
    prod.replacefds1 = {{'faces', 'faces1'}};
    prod.replacefds2 = {{'rowdim1', 'coldim1'}, {'rowdim2', 'coldim2'}, {'faces', 'faces2'}};
    prod.mergefds = {'cells', 'nodes'};
    prod.reducefds = {'coldim'};
    prod = prod.setup();
    
    nCg = prod.eval(facetNormals, Cg);
    
    %% We setup A11
    
    map = TensorMap();
    map.fromTbl = cellnodeface2col2tbl;
    map.toTbl = nodeface2col2tbl;
    map.mergefds = {'nodes', 'faces1', 'faces2', 'coldim1', 'coldim2'};
    map = map.setup();
    
    A11 = map.eval(nCg);
    
    prod = TensorProd();
    prod.tbl1 = nodeface2col2tbl;
    prod.tbl2 = nodefacecoltbl;
    prod.tbl3 = nodefacecoltbl;
    prod.replacefds1 = {{'faces1', 'faces'}, {'coldim1', 'coldim'}};
    prod.replacefds2 = {{'faces', 'faces2'}, {'coldim', 'coldim2'}};
    prod.mergefds = {'nodes'};
    prod.reducefds = {'faces2', 'coldim2'};
    prod = prod.setup();
    
    A11_T = SparseTensor();
    A11_T = A11_T.setFromTensorProd(A11, prod);
    A11 = A11_T.getMatrix();
    
    %% We setup A12
    
    map = TensorMap();
    map.fromTbl = cellnodeface2col2tbl;
    map.toTbl = cellnodefacecol2tbl;
    map.replaceFromTblfds = {{'faces1', 'faces'}};
    map.mergefds = {'cells', 'nodes', 'faces', 'coldim1', 'coldim2'};
    map = map.setup();
    
    % note the minus sign
    A12 = map.eval( - nCg);
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacecol2tbl;
    prod.tbl2 = cellcoltbl;
    prod.tbl3 = nodefacecoltbl;
    prod.replacefds1 = {{'coldim1', 'coldim'}};
    prod.replacefds2 = {{'coldim', 'coldim2'}};
    prod.reducefds = {'cells', 'coldim2'};
    prod = prod.setup();
    
    A12_T = SparseTensor();
    A12_T = A12_T.setFromTensorProd(A12, prod);
    A12 = A12_T.getMatrix();
    
    %% We setup A21
    
    map = TensorMap();
    map.fromTbl = cellnodeface2col2tbl;
    map.toTbl = cellnodefacecol2tbl;
    map.replaceFromTblfds = {{'faces2', 'faces'}};
    map.mergefds = {'cells', 'nodes', 'faces', 'coldim1', 'coldim2'};
    map = map.setup();
    
    % note the minus sign
    A21 = map.eval( - nCg);
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacecol2tbl;
    prod.tbl2 = nodefacecoltbl;
    prod.tbl3 = cellcoltbl;
    prod.replacefds1 = {{'coldim2', 'coldim'}};
    prod.replacefds2 = {{'coldim', 'coldim1'}};
    prod.reducefds = {'nodes', 'faces', 'coldim1'};
    prod = prod.setup();
    
    A21_T = SparseTensor();
    A21_T = A21_T.setFromTensorProd(A21, prod);
    A21 = A21_T.getMatrix();
    
    %% We setup A22
    
    map = TensorMap();
    map.fromTbl = cellnodeface2col2tbl;
    map.toTbl = cellcol2tbl;
    map.mergefds = {'cells', 'coldim1', 'coldim2'};
    map = map.setup();
    
    A22 = map.eval(nCg);
    
    prod = TensorProd();
    prod.tbl1 = cellcol2tbl;
    prod.tbl2 = cellcoltbl;
    prod.tbl3 = cellcoltbl;
    prod.replacefds1 = {{'coldim1', 'coldim'}};
    prod.replacefds2 = {{'coldim', 'coldim2'}};
    prod.mergefds = {'cells'};
    prod.reducefds = {'coldim2'};
    prod = prod.setup();
    
    A22_T = SparseTensor();
    A22_T = A22_T.setFromTensorProd(A22, prod);
    A22 = A22_T.getMatrix();
    
    % Uses the block structure for the local reduction
    % We count the number of degrees of freedom that are connected to the same
    % node.
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
    prod.dispind2 = sub2ind([d_num, nf_num], c, nodeface_from_cellnodeface(i));
    prod.dispind3 = cell_from_cellnode(cellnode_from_cellnodeface(i));
    prod.issetup = true;
    
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




