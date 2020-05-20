function assembly = assembleMPFA(G, K, bcstruct, src, eta, tbls, mappings, varargin)

    opt = struct('bcetazero', true);
    opt = merge_options(opt, varargin{:});

    cellnodefacecoltbl = tbls.cellnodefacecoltbl;
    cellcolrowtbl      = tbls.cellcolrowtbl;
    
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
    cellnodefacecolrowtbl = tbls.cellnodefacecolrowtbl;
    colrowtbl             = tbls.colrowtbl;
    nodecolrowtbl         = tbls.nodecolrowtbl;
    col2row2tbl           = tbls.col2row2tbl;
    cellcol2row2tbl       = tbls.cellcol2row2tbl;
    cellnodecol2row2tbl   = tbls.cellnodecol2row2tbl;

    %  g belongs to cellnodefacecoltbl;
    g = computeConsistentGradient(G, eta, tbls, mappings);

    % facetNormals belongs to cellnodefacecoltbl;
    normals = computeFacetNormals(G, cellnodefacetbl);

    % K belongs to cellcolrowtbl

    prod = TensorProd();
    prod.tbl1 = cellcolrowtbl;
    prod.tbl2 = cellnodefacecoltbl;
    prod.tbl3 = cellnodefacecoltbl;
    prod.replacefds2 = {{'coldim', 'rowdim'}};
    prod.mergefds = {'cells'};
    prod.reducefds = {'rowdim'};
    prod = prod.setup();
    
    Kg = prod.eval(K, g);
    
    cellnodeface2tbl = crossIndexArray(cellnodefacetbl, cellnodefacetbl, {'cells', 'nodes'}, 'crossextend', {{'faces', {'faces1', 'faces2'}}});
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacecoltbl;
    prod.tbl2 = cellnodefacecoltbl;
    prod.tbl3 = cellnodeface2tbl;
    prod.replacefds1 = {{'faces', 'faces1'}};
    prod.replacefds2 = {{'faces', 'faces2'}};
    prod.mergefds = {'cells', 'nodes'};
    prod.reducefds = {'coldim'};
    prod = prod.setup();
    
    nKg = prod.eval(normals, Kg);
    
    nodeface2tbl = crossIndexArray(nodefacetbl, nodefacetbl, {'nodes'}, 'crossextend', {{'faces', {'faces1', 'faces2'}}});
    
    %% Setup A11 matrix (facenode dof -> facenode dof)
    
    map = TensorMap();
    map.fromTbl = cellnodeface2tbl;
    map.toTbl = nodeface2tbl;
    map.mergefds = {'nodes', 'faces1', 'faces2'};
    map = map.setup();
    
    A11 = map.eval(nKg);
    
    prod = TensorProd();
    prod.tbl1 = nodeface2tbl;
    prod.tbl2 = nodefacetbl;
    prod.tbl3 = nodefacetbl;
    prod.replacefds1 = {{'faces1', 'faces'}};
    prod.replacefds2 = {{'faces', 'faces2'}};
    prod.mergefds = {'nodes'};
    prod.reducefds = {'faces2'};
    prod = prod.setup();
    
    A11_T = SparseTensor();
    A11_T = A11_T.setFromTensorProd(A11, prod);
    A11 = A11_T.getMatrix();
    
    [~, sz] = rlencode(nodefacetbl.get('nodes'), 1);
    opt.invertBlocks = 'mex';
    bi = blockInverter(opt);
    invA11 = bi(A11, sz);

    
    %% Setup A12 matrix (cell dof -> facenode dof)    
   
    map = TensorMap();
    map.fromTbl = cellnodeface2tbl;
    map.toTbl = cellnodefacetbl;
    map.replaceFromTblfds = {{'faces1', 'faces'}};
    map.mergefds = {'cells', 'nodes', 'faces'};
    map = map.setup();
    
    % beware minus sign here
    A12 = - map.eval(nKg);
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacetbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = nodefacetbl;
    prod.reducefds = {'cells'};
    prod = prod.setup();
    
    A12_T = SparseTensor();
    A12_T = A12_T.setFromTensorProd(A12, prod);
    A12 = A12_T.getMatrix();
    
    assembly = struct('A11', A11, ...
                      'A12', A12, ...
                      'invA11', invA11);
    
    
    %% Setup A21 matrix (facenode dof -> cell dof)    
    
   
    map = TensorMap();
    map.fromTbl = cellnodeface2tbl;
    map.toTbl = cellnodefacetbl;
    map.replaceFromTblfds = {{'faces2', 'faces'}};
    map.mergefds = {'cells', 'nodes', 'faces'};
    map = map.setup();
    
    % beware minus sign here
    A21 = -map.eval(nKg);
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacetbl;
    prod.tbl2 = nodefacetbl;
    prod.tbl3 = celltbl;
    prod.reducefds = {'faces', 'nodes'};
    prod = prod.setup();
    
    A21_T = SparseTensor();
    A21_T = A21_T.setFromTensorProd(A21, prod);
    A21 = A21_T.getMatrix();
    
    
    %% Setup A22 matrix (cell dof -> cell dof)    
   
    map = TensorMap();
    map.fromTbl = cellnodeface2tbl;
    map.toTbl = celltbl;
    map.mergefds = {'cells'};
    map = map.setup();
    
    A22 = map.eval(nKg);
    
    prod = TensorProd();
    prod.tbl1 = celltbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = celltbl;
    prod.mergefds = {'cells'};
    prod = prod.setup();
    
    A22_T = SparseTensor();
    A22_T = A22_T.setFromTensorProd(A22, prod);
    A22 = A22_T.getMatrix();

    % We enforce the Dirichlet boundary conditions as Lagrange multipliers
    bcdirichlet = bcstruct.bcdirichlet;
    [D, bcvals] = setupMpfaNodeFaceBc(bcdirichlet, tbls);
    
    % the solution is given by the system
    %
    % A = [[A11, A12, -D];
    %      [A21, A22,  0];
    %      [D' , 0  ,  0]];
    %
    % u = [u  (displacement at nodefacecoltbl);
    %      u  (displacement at cellcoltbl);
    %      lagmult];
    %
    % f = [extforce  (force at nodefacecoltbl);
    %      force  (volumetric force at cellcoltbl);
    %      bcvals (for the linear form at the boundary)];
    %
    % A*u = f
    %
    % Note: extforce is sparse and should only give contribution at facets
    % that are at the boundary
    %
    % By construction of the method, the matrix A11 is block-diagonal. Hence,
    % we invert it directly and reduce to a cell-centered scheme.
    
    matrices = struct('A11', A11, ...
                      'A12', A12, ...
                      'A21', A21, ...
                      'A22', A22, ...
                      'D'  , D  , ...
                      'invA11', invA11);
    
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
    % u = [u (displacement at celltbl);
    %      lagmult];
    %
    % rhs = [-A21*invA11*extflux;  +  [src;
    %        -D'*invA11*extflux  ]     bcvals]
    
    B11 = A22 - A21*invA11*A12;
    B12 = A21*invA11*D;
    B21 = -D'*invA11*A12;
    B22 = D'*invA11*D;

    if ~isempty(bcstruct.bcneumann)
        error('not yet implemented');
    else
        nf_num = nodefacetbl.num;
        extflux = zeros(nf_num, 1);
    end

    B = [[B11, B12]; ...
         [B21, B22]];
    
    rhs{1} = -A21*invA11*extflux + src; 
    rhs{2} = -D'*invA11*extflux + bcvals;
    
    rhs = vertcat(rhs{:});
    
    assembly = struct( 'B'  , B  , ...
                       'rhs', rhs, ...
                       'matrices', matrices);
end

