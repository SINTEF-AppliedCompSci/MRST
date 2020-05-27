function assembly = blockAssembleMPFA(G, K, bcstruct, src, eta, globtbls, globmappings, varargin)

    opt = struct('verbose'  , mrstVerbose, ...
                 'blocksize', [], ...
                 'bcetazero', false, ...
                 'useVirtual', false);
    opt = merge_options(opt, varargin{:});
    
    useVirtual = opt.useVirtual;
    blocksize = opt.blocksize;
    
    
    nn = G.nodes.num;
    nblocks = floor(nn/blocksize);
    blocksizes = repmat(blocksize, nblocks, 1);
    if nn > nblocks*blocksize
        blocksizes = [blocksizes; nn - nblocks*blocksize];
    end
    nblocks = numel(blocksizes);
    blockinds = cumsum([1; blocksizes]);

    globcoltbl                = globtbls.coltbl;
    globcolrowtbl             = globtbls.colrowtbl;
    globcol2row2tbl           = globtbls.col2row2tbl;
    globnodetbl               = globtbls.nodetbl;
    globcellcoltbl            = globtbls.cellcoltbl;
    globnodecoltbl            = globtbls.nodecoltbl;
    globcellnodetbl           = globtbls.cellnodetbl;
    globcellnodecoltbl        = globtbls.cellnodecoltbl;
    globnodefacecoltbl        = globtbls.nodefacecoltbl;
    globcellcol2row2tbl       = globtbls.cellcol2row2tbl;
    globcellcolrowtbl         = globtbls.cellcolrowtbl;
    globcoltbl                = globtbls.coltbl;
    globcelltbl               = globtbls.celltbl;
    globnodetbl               = globtbls.nodetbl;
    globcellfacetbl           = globtbls.cellfacetbl;
    globcellnodetbl           = globtbls.cellnodetbl;
    globnodefacetbl           = globtbls.nodefacetbl;
    globcellcoltbl            = globtbls.cellcoltbl;
    globnodecoltbl            = globtbls.nodecoltbl;
    globnodefacecoltbl        = globtbls.nodefacecoltbl;
    globcellnodefacetbl       = globtbls.cellnodefacetbl;
    globcellnodecoltbl        = globtbls.cellnodecoltbl;
    globcellnodecolrowtbl     = globtbls.cellnodecolrowtbl;
    globcellnodefacecolrowtbl = globtbls.cellnodefacecolrowtbl;
    globcolrowtbl             = globtbls.colrowtbl;
    globnodecolrowtbl         = globtbls.nodecolrowtbl;
    globcol2row2tbl           = globtbls.col2row2tbl;
    globcellcol2row2tbl       = globtbls.cellcol2row2tbl;
    globcellnodecol2row2tbl   = globtbls.cellnodecol2row2tbl;
    
    if ~isempty(bcstruct.bcneumann)
        error('not yet implemented');
    else
        nf_num = globnodefacetbl.num;
        globextflux = zeros(nf_num, 1);
    end

    if isempty(src)
        src = zeros(globcelltbl.num, 1);
    end
    
    globsrc = src;
    
    globbcdirichlet = bcstruct.bcdirichlet;
    globbcnodefacetbl = globbcdirichlet.bcnodefacetbl;
    globbcvals        = globbcdirichlet.bcvals;
    
    gnc = globcelltbl.num;
    gnbc = globbcnodefacetbl.num;
    
    B11 = sparse(gnc, gnc);
    B12 = sparse(gnc, gnbc);
    B21 = sparse(gnbc, gnc);
    B22 = sparse(gnbc, gnbc);

    rhsc  = zeros(gnc, 1);
    rhsbc = zeros(gnbc, 1);
    
    for iblock = 1 : nblocks

        %% Construction of tensor g (as defined in paper eq 4.1.2)
        nodes = [blockinds(iblock) : (blockinds(iblock + 1) - 1)]';

        clear nodetbl;
        nodetbl.nodes = nodes;
        nodetbl = IndexArray(nodetbl);

        if opt.verbose
            fprintf('Assembling block %d/%d (%d nodes)\n', iblock, nblocks, nodetbl.num);
        end
        
        [tbls, mappings] = setupStandardBlockTables(G, nodetbl, globtbls, 'useVirtual', useVirtual);

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
        nodecolrowtbl         = tbls.nodecolrowtbl;
        cellcol2row2tbl       = tbls.cellcol2row2tbl;
        cellnodecol2row2tbl   = tbls.cellnodecol2row2tbl;
        cellcolrowtbl         = tbls.cellcolrowtbl;
        
        %  g belongs to cellnodefacecoltbl;
        g = computeConsistentGradient(G, eta, tbls, mappings);

        % facetNormals belongs to cellnodefacecoltbl;
        normals = computeFacetNormals(G, cellnodefacetbl);

        % K belongs to cellcolrowtbl

        % Set up product of K (in cellcolrowtbl) with g (in cellnodefacecoltbl)
        prod = TensorProd();
        prod.tbl1 = cellcolrowtbl;
        prod.tbl2 = cellnodefacecoltbl;
        prod.tbl3 = cellnodefacecoltbl;
        prod.replacefds2 = {{'coldim', 'rowdim'}};
        prod.mergefds = {'cells'};
        prod.reducefds = {'rowdim'};
        prod = prod.setup();
        
        Kg = prod.eval(K, g);
        
        % Set up space for local mapping around the nodes.
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
        bcnodefacetbl = crossIndexArray(globbcnodefacetbl, nodefacetbl, {'nodes', 'faces'});
        
        bcterm_exists = true;
        if bcnodefacetbl.num == 0
            bcterm_exists = false;
        end
        
        if bcterm_exists
            map = TensorMap();
            map.fromTbl = globbcnodefacetbl;
            map.toTbl = bcnodefacetbl;
            map.mergefds = {'nodes', 'faces'};
            bcind = map.getDispatchInd();
        
            clear bcdirichlet;
            bcdirichlet.bcnodefacetbl = bcnodefacetbl;
            bcdirichlet.bcvals = []; % not used locally
            
            [D, ~] = setupMpfaNodeFaceBc(bcdirichlet, tbls);
            
        end
        
        map = TensorMap();
        map.fromTbl = globnodefacetbl;
        map.toTbl = nodefacetbl;
        map.mergefds = {'nodes', 'faces'};
        map = map.setup();
        
        extflux = map.eval(globextflux);

        locB11 = A22 - A21*invA11*A12;
        locrhsc = -A21*invA11*extflux; 
        
        if bcterm_exists
            locB12 = A21*invA11*D;
            locB21 = -D'*invA11*A12;
            locB22 = D'*invA11*D;
            locrhsbc = -D'*invA11*extflux;
        end

        % locB11 : celltbl       -> celltbl
        % locB22 : bcnodefacetbl -> bcnodefacetbl
        % locB12 : bcnodefacetbl -> celltbl
        % locB21 : celltbl       -> bcnodefacetbl 
        % Above, we recall that these index arrays are all local.

        nc = celltbl.num;
        map = TensorMap();
        map.fromTbl = globcelltbl;
        map.toTbl = celltbl;
        map.mergefds = {'cells'};
        cellind = map.getDispatchInd();
        
        B11 = B11 + sparse(repmat(cellind, 1, nc), repmat(cellind', nc, 1), locB11, gnc, gnc);
        rhsc(cellind) = rhsc(cellind) + locrhsc;
        
        if bcterm_exists
            nbc = bcnodefacetbl.num;
            B22 = B22 + sparse(repmat(bcind, 1, nbc),   repmat(bcind', nbc, 1),   locB22, gnbc, gnbc);
            B12 = B12 + sparse(repmat(cellind, 1, nbc), repmat(bcind', nc, 1),    locB12, gnc,  gnbc);
            B21 = B21 + sparse(repmat(bcind, 1, nc),    repmat(cellind', nbc, 1), locB21, gnbc, gnc);
            rhsbc(bcind) = rhsbc(bcind) + locrhsbc;
        end
        
    end
    
    B = [[B11, B12]; ...
         [B21, B22]];

    rhsc = rhsc + globsrc;
    rhsbc = rhsbc + globbcvals;
    
    rhs = [rhsc; ...
           rhsbc];
    
    assembly = struct( 'B'  , B  , ...
                       'rhs', rhs);    
end
