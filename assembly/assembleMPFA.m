function assembly = assembleMPFA(G, K, bcstruct, src, eta, tbls, mappings, varargin)
%
%   
% Before nodal reduction, the solution is given by the system
%
% A = [[A11, A12, -D];
%      [A21, A22,  0];
%      [D' , 0  ,  0]];
%
% u = [pnf     (pressure at nodefacetbl);
%      pc      (pressure at celltbl);
%      lagmult (flux values Dirichlet boundary)];
%
% f = [flux   (flux at nodefacetbl);
%      src    (source at celltbl);
%      bcvals (pressure values at Dirichlet boundary)];
%
% A*u = f
%
% Note: extforce is sparse and should only give contribution at facets
% that are at the boundary
%
% By construction of the method, the matrix A11 is block-diagonal. Hence,
% we invert it directly and reduce to a cell-centered scheme.
% 
%
% First, we assemble the matrices A11, A12, A21, A22

    
    opt = struct('bcetazero'           , true , ...
                 'addAssemblyMatrices' , false, ...
                 'addAdOperators'      , false, ...
                 'onlyAssemblyMatrices', false, ...
                 'dooptimize'          , true);
    
    opt = merge_options(opt, varargin{:});
    
    dooptimize = opt.dooptimize;
    
    if tbls.useVirtual
        assert(dooptimize, 'We cannot use virtual tables in a non-optimized run of assembleMPFA');
    end
    
    cellcolrowtbl         = tbls.cellcolrowtbl;
    cellnodecolrowtbl     = tbls.cellnodecolrowtbl;
    cellnodeface2coltbl   = tbls.cellnodeface2coltbl;
    cellnodeface2tbl      = tbls.cellnodeface2tbl;
    cellnodefacecolrowtbl = tbls.cellnodefacecolrowtbl;
    cellnodefacecoltbl    = tbls.cellnodefacecoltbl;
    cellnodefacetbl       = tbls.cellnodefacetbl;
    celltbl               = tbls.celltbl;
    coltbl                = tbls.coltbl;
    nodeface2tbl          = tbls.nodeface2tbl;
    nodefacecoltbl        = tbls.nodefacecoltbl;
    nodefacetbl           = tbls.nodefacetbl;
    
    if dooptimize
        % fetch the index mappings to set explictly the tensor products or tensor mappings
        cell_from_cellnodeface     = mappings.cell_from_cellnodeface;
        nodeface_from_cellnodeface = mappings.nodeface_from_cellnodeface;
        cellnodeface_1_from_cellnodeface2 = mappings.cellnodeface_1_from_cellnodeface2;
        cellnodeface_2_from_cellnodeface2 = mappings.cellnodeface_2_from_cellnodeface2;
        nodeface_1_from_nodeface2 = mappings.nodeface_1_from_nodeface2;
        nodeface_2_from_nodeface2 = mappings.nodeface_2_from_nodeface2;
    end
    
    % Some shortcuts
    c_num     = celltbl.num;
    cnf_num   = cellnodefacetbl.num;
    nf_num    = nodefacetbl.num;
    cnfcr_num = cellnodefacecolrowtbl.num;
    d_num     = coltbl.num;
    
    %  g belongs to cellnodefacecoltbl;
    g = computeConsistentGradient(G, eta, tbls, mappings, 'bcetazero', opt.bcetazero);

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
    
    if dooptimize
        prod.pivottbl = cellnodefacecolrowtbl;
        [r, c, i] = ind2sub([d_num, d_num, cnf_num], (1 : cnfcr_num)');
        prod.dispind1 = sub2ind([d_num, d_num, c_num], r, c, cell_from_cellnodeface(i));
        prod.dispind2 = sub2ind([d_num, cnf_num], r, i);
        prod.dispind3 = sub2ind([d_num, cnf_num], c, i);
        prod.issetup = true;
    else
        prod = prod.setup();
    end
    
    Kg = prod.eval(K, g);
    
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacecoltbl;
    prod.tbl2 = cellnodefacecoltbl;
    prod.tbl3 = cellnodeface2tbl;
    prod.replacefds1 = {{'faces', 'faces1'}};
    prod.replacefds2 = {{'faces', 'faces2'}};
    prod.mergefds = {'cells', 'nodes'};
    prod.reducefds = {'coldim'};
    
    if dooptimize
        cnf2_num = cellnodeface2tbl.num;
        cnf2c_num = cellnodeface2coltbl.num;
        prod.pivottbl = cellnodeface2coltbl;
        [c, i] = ind2sub([d_num, cnf2_num], (1 : cnf2c_num)');
        prod.dispind1 = sub2ind([d_num, cnf_num], c, cellnodeface_1_from_cellnodeface2(i));
        prod.dispind2 = sub2ind([d_num, cnf_num], c, cellnodeface_2_from_cellnodeface2(i));
        prod.dispind3 = i;
        prod.issetup = true;
    else
        prod = prod.setup();
    end
    
    nKg = prod.eval(normals, Kg);
    
    
    %% Setup A11 matrix (facenode dof -> facenode dof)
    
    map = TensorMap();
    map.fromTbl = cellnodeface2tbl;
    map.toTbl = nodeface2tbl;
    map.mergefds = {'nodes', 'faces1', 'faces2'};
    map = map.setup(); % not optimized (use generic setup function)
    
    A11 = map.eval(nKg);
    
    prod = TensorProd();
    prod.tbl1 = nodeface2tbl;
    prod.tbl2 = nodefacetbl;
    prod.tbl3 = nodefacetbl;
    prod.replacefds1 = {{'faces1', 'faces'}};
    prod.replacefds2 = {{'faces', 'faces2'}};
    prod.mergefds = {'nodes'};
    prod.reducefds = {'faces2'};
    
    if dooptimize
        prod.pivottbl = nodeface2tbl;
        i = (1 : nodeface2tbl.num)';
        prod.dispind1 = i;
        prod.dispind2 = nodeface_2_from_nodeface2(i);
        prod.dispind3 = nodeface_1_from_nodeface2(i);
        prod.issetup = true;
    else
        prod = prod.setup();
    end
    
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
    
    if dooptimize
        map.pivottbl = cellnodeface2tbl;
        cnf2_num = cellnodeface2tbl.num;
        i = (1 : cnf2_num)';
        map.dispind1 = i;
        map.dispind2 = cellnodeface_1_from_cellnodeface2(i);
        map.issetup = true;
    else
        map = map.setup();
    end
    
    % beware minus sign here
    A12 = - map.eval(nKg);
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacetbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = nodefacetbl;
    prod.reducefds = {'cells'};
    
    if dooptimize
        prod.pivottbl = cellnodefacetbl;
        i = (1 : cellnodefacetbl.num)';
        prod.dispind1 = i;
        prod.dispind2 = cell_from_cellnodeface(i);
        prod.dispind3 = nodeface_from_cellnodeface(i);
        prod.issetup = true;
    else
        prod = prod.setup();
    end
    
    A12_T = SparseTensor();
    A12_T = A12_T.setFromTensorProd(A12, prod);
    A12 = A12_T.getMatrix();

    
    %% Setup A21 matrix (facenode dof -> cell dof)    
    
   
    map = TensorMap();
    map.fromTbl = cellnodeface2tbl;
    map.toTbl = cellnodefacetbl;
    map.replaceFromTblfds = {{'faces2', 'faces'}};
    map.mergefds = {'cells', 'nodes', 'faces'};
    
    if dooptimize
        map.pivottbl = cellnodeface2tbl;
        i = (1 : cellnodeface2tbl.num)';
        map.dispind1 = i;
        map.dispind2 = cellnodeface_2_from_cellnodeface2(i);
        map.issetup = true;
    else
        map = map.setup();
    end
    
    % beware minus sign here
    A21 = -map.eval(nKg);
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacetbl;
    prod.tbl2 = nodefacetbl;
    prod.tbl3 = celltbl;
    prod.reducefds = {'faces', 'nodes'};
    prod = prod.setup();
    
    if dooptimize
        prod.pivottbl = cellnodefacetbl;
        i = (1 : cellnodefacetbl.num)';
        prod.dispind1 = i;
        prod.dispind2 = nodeface_from_cellnodeface(i);
        prod.dispind3 = cell_from_cellnodeface(i);
    else
        prod = prod.setup();
    end
    
    A21_T = SparseTensor();
    A21_T = A21_T.setFromTensorProd(A21, prod);
    A21 = A21_T.getMatrix();
    
    
    %% Setup A22 matrix (cell dof -> cell dof)    
   
    map = TensorMap();
    map.fromTbl = cellnodeface2tbl;
    map.toTbl = celltbl;
    map.mergefds = {'cells'};

    if dooptimize
        map.pivottbl = cellnodeface2tbl;
        i = (1 : cellnodeface2tbl.num);
        map.dispind1 = i;
        map.dispind2 = cell_from_cellnodeface(cellnodeface_1_from_cellnodeface2(i));
        map.issetup = true;
    else
        map = map.setup();
    end
    
    A22 = map.eval(nKg);
    
    prod = TensorProd();
    prod.tbl1 = celltbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = celltbl;
    prod.mergefds = {'cells'};
    
    if dooptimize
        prod.pivottbl = celltbl;
        i = (1 : celltbl.num);
        prod.dispind1 = i;
        prod.dispind2 = i;
        prod.dispind3 = i;
        prod.issetup = true;
    else
        prod = prod.setup();
    end
    
    A22_T = SparseTensor();
    A22_T = A22_T.setFromTensorProd(A22, prod);
    A22 = A22_T.getMatrix();
    

    % We enforce the Dirichlet boundary conditions as Lagrange multipliers
    if ~isempty(bcstruct.bcneumann)
        error('not yet implemented');
    else
        nf_num = nodefacetbl.num;
        extflux = zeros(nf_num, 1);
    end
    

    if isempty(src)
        src = zeros(celltbl.num, 1);
    end
    
    bcdirichlet = bcstruct.bcdirichlet;
    if ~isempty(bcdirichlet)
        [D, bcvals] = setupMpfaNodeFaceBc(bcdirichlet, tbls);
        % scale D
        fac    = max(max(abs(A11)));
        D      = fac*D;
        bcvals = fac*bcvals;
    else
        D      = ones(size(A11, 1), 0);
        bcvals = [];
    end
    
    fullrhs{1} = extflux;
    fullrhs{2} = src;
    fullrhs{3} = bcvals;
    
    if (opt.onlyAssemblyMatrices | opt.addAssemblyMatrices | opt.addAdOperators)
        
        matrices = struct('A11', A11, ...
                          'A12', A12, ...
                          'A21', A21, ...
                          'A22', A22, ...
                          'D'  , D  , ...
                          'invA11', invA11);
        
        matrices.fullrhs = fullrhs;
        assembly.matrices = matrices;
        assembly.nKg = nKg;
        
        if opt.onlyAssemblyMatrices
            return
        end
        
    end
    
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
    % u = [p (pressure at celltbl);
    %      lagmult];
    %
    % rhs = [-A21*invA11*extflux;  +  [src;
    %        -D'*invA11*extflux  ]     bcvals]
    
    B11 = A22 - A21*invA11*A12;
    B12 = A21*invA11*D;
    B21 = -D'*invA11*A12;
    B22 = D'*invA11*D;


    B = [[B11, B12]; ...
         [B21, B22]];
    
    adrhs{1} = -A21*invA11*extflux + src; 
    adrhs{2} = -D'*invA11*extflux + bcvals;
    
    rhs = vertcat(adrhs{:});
    
    assembly.B = B;
    assembly.rhs = rhs;
    
    if opt.addAdOperators
        
        adB = cell(2, 2);
        adB{1, 1} = B11;
        adB{2, 1} = B21;
        adB{1, 2} = B12;
        adB{2, 2} = B22;
        
        adoperators.B     = adB;
        adoperators.rhs   = adrhs;        
        
        % Setup fluid flux operator
        mpfaKgrad = setupMpfaFlux(G, assembly, tbls);
        fluxop = @(p) fluxFunc(p, mpfaKgrad);
        
        adoperators.fluxop = fluxop;
        
        assembly.adoperators = adoperators;
        
    end
    
end

function flux = fluxFunc(p, mpfaKgrad)
   flux = mpfaKgrad*p;
end
