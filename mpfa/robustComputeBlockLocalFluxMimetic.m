function [B, tbls] = robustComputeBlockLocalFluxMimetic(G, rock, opt)

    % Some short aliases 
    nc = G.cells.num;
    nf = G.faces.num;
    dim = G.griddim;
   
    cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos)); 
    cellfacetbl.faces = G.cells.faces(:, 1);
    cellfacetbl.num   = numel(cellfacetbl.cells);

    faces = rldecode((1 : nf)', diff(G.faces.nodePos)); 
    nodes = G.faces.nodes;
    % We setup the face-node table and it is ordered along ascending node numbers so
    % that we will have a block structure for the nodal scalar product.
    collatefacenode = [nodes, faces];
    collatefacenode = sortrows(collatefacenode, 1);
    facenodetbl.nodes = collatefacenode(:, 1);
    facenodetbl.faces = collatefacenode(:, 2);
    facenodetbl.num = numel(facenodetbl.faces);
    
    [op, cellnodefacetbl] = setupTableMapping(cellfacetbl, facenodetbl, {'faces'});

    % We setup the cell-face-node table, cellnodefacetbl. Each entry determine a
    % unique facet in a corner
    % We order cellfacenode in cell-node-face order. This is node to optimize
    % for-end loop below.
    orderingmat = [cellnodefacetbl.cells, cellnodefacetbl.nodes, ...
                   cellnodefacetbl.faces];
    orderingmat = sortrows(orderingmat);
    cellnodefacetbl.cells = orderingmat(:, 1);
    cellnodefacetbl.nodes = orderingmat(:, 2);
    cellnodefacetbl.faces = orderingmat(:, 3);
    
    % Somer shortcuts
    cno = cellnodefacetbl.cells;
    fno = cellnodefacetbl.faces;
    nno = cellnodefacetbl.nodes;
    
    % We setup the cell-node table, cellnodetbl. Each entry determine
    % a unique corner
    cfn = cellnodefacetbl;
    cellnode = [cfn.cells, cfn.nodes];
    cellnode = unique(cellnode, 'rows');
    cellnodetbl.cells = cellnode(:, 1);
    cellnodetbl.nodes = cellnode(:, 2);
    cellnodetbl.num   = numel(cellnodetbl.nodes);
    orderingmat = [cellnodetbl.cells, cellnodetbl.nodes];
    orderingmat = sortrows(orderingmat);
    cellnodetbl.cells = orderingmat(:, 1);
    cellnodetbl.nodes = orderingmat(:, 2);    
    
    % Nodal scalar product is stored in vector nodeM
    % mattbl is the table which specifies how nodeM is stored: a matrix for
    % each "corner" (cell-node pair).
    duplicate = {'faces', {'faces1', 'faces2'}};
    [~, mattbl] = setupTableMapping(cellnodefacetbl, cellnodefacetbl, {'cells', 'nodes'}, ...
                                         'duplicate', {duplicate});
    % We order mattbl in cell-node-face1-face2 order
    % This is done to optimize for-end loop below
    orderingmat = [mattbl.cells, mattbl.nodes, ...
                   mattbl.faces1, mattbl.faces2];
    orderingmat = sortrows(orderingmat);
    mattbl.cells  = orderingmat(:, 1);
    mattbl.nodes  = orderingmat(:, 2);
    mattbl.faces1 = orderingmat(:, 3);    
    mattbl.faces2 = orderingmat(:, 4);    
    
    
    numnodes = double(diff(G.faces.nodePos));
    numnodes = numnodes(fno);
    facetNormals = G.faces.normals(fno, :);
    facetNormals = bsxfun(@ldivide, numnodes, facetNormals);
    
    sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
    facetNormals = sgn.*facetNormals; % outward normals.
    
    % Assemble facePermNormals which corresponds to $Kn$ where n are the *outward*
    % normals at the facets.
    [perm, r, c] = permTensor(rock, G.griddim);
    Kn = bsxfun(@times, perm(cno, :), facetNormals(:, c));
    Kn = rldecode(Kn, dim*ones(numel(cno, 1)));
    coef = sparse(r, (1 : dim*dim)', ones(dim*dim, 1), dim, dim*dim);
    coef = repmat(coef, numel(cno), 1);
    Kn = coef.*Kn;
    Kn = sum(Kn, 2);
    Kn = reshape(Kn, dim, [])';
    
    facePermNormals = Kn;
    
    % Default option (opt.eta = 0): Use original face centroids and cell centroids,
    % NOT actual subface centroids. This corresponds to an MPFA method
    % (O-method) R = G.faces.centroids(fno,:) - G.cells.centroids(cno,:);
    cellFacetVec = G.faces.centroids(fno,:) - G.cells.centroids(cno,:) + ...
        opt.eta*(G.nodes.coords(nno,:) - G.faces.centroids(fno,:));
    
    % set up areas and volumes
    areas = G.faces.areas(fno);
    vols  = G.cells.volumes(cno);
    
    map = setupTableMapping(cellnodetbl, cellnodefacetbl, {'cells', 'nodes'}); 
    nfaces1 = full(diag(map'*map));
    nfaces2 = map*nfaces1;
    
    blocksize = opt.blocksize;
    ncn = cellnodetbl.num;
    nblocks = floor(ncn/blocksize);
    blocksizes = [repmat(blocksize, nblocks, 1); ...
                  ncn - nblocks*blocksize];
    nblocks = numel(blocksizes);
    blockinds = cumsum([1; blocksizes]);
    nodeMs = cell(nblocks, 1);

    cnf_i = 1; % start indice for the cellnodefacetbl index
    
    if opt.verbose
        fprintf('Number of blocks: %d\n', nblocks);
    end
    
    for i = 1 : nblocks
    
        blocksize = blocksizes(i);
        ind = [blockinds(i) : (blockinds(i + 1) - 1)];
        nblockfaces = nfaces1(ind);
        nodeMs{i} = zeros(sum(nblockfaces.^2), 1);
        mat_i = 1;

        if opt.verbose
            t0 = tic;
        end
    
        for j = 1 : blocksize
    
            nface = nfaces2(cnf_i);
            cnfind = cnf_i : (cnf_i + (nface - 1));
    
            N     = facePermNormals(cnfind, :); 
            R     = cellFacetVec(cnfind, :);
            a     = areas(cnfind);
            v     = vols(cnfind);
            faces = cellnodefacetbl.faces(cnfind);
    
            cellno = cellnodefacetbl.cells(cnf_i);

            K = reshape(perm(cellno, :), [dim, dim]);
    
            % Assemble local nodal scalar product ( function node_ip2 below handle case when
            % N is invertible)
            locM = node_ip(a, v, ...
                           full(N), ...
                           full(R), ...
                           K);
            locM = reshape(locM', [], 1);
    
            assert(numel(locM) == nface*nface, 'mismatch');
            matind = mat_i : (mat_i + (nface*nface - 1));
            nodeMs{i}(matind) = nodeMs{i}(matind) + locM;
    
            cnf_i = cnf_i + nface;
            mat_i = mat_i + nface*nface;        
    
        end

        if opt.verbose
            t0 = toc(t0);
            fprintf('Assembly of block %d (block matrix size : %d) done in %g\n', i, sum(nblockfaces), t0);
        end
    
    end

    % Condensate on nodes (sum up cell contributions for give node).
    op = setupTableMapping(facenodetbl, facenodetbl, {'nodes'});
    [colind, rowind] = find(op);
    redmattbl.nodes  = facenodetbl.nodes(colind);
    redmattbl.faces1 = facenodetbl.faces(colind);
    redmattbl.faces2 = facenodetbl.faces(rowind);
    redmattbl.num    = numel(redmattbl.nodes);
    nodeM = zeros(redmattbl.num, 1);
    
    mat_i = 1;
    for i = 1 : nblocks
        if opt.verbose
            t0 = tic;
        end
        blocksize = blocksizes(i);
        ind = [blockinds(i) : (blockinds(i + 1) - 1)];
        nblockfaces = nfaces1(ind);
        nblockfacessq = sum(nblockfaces.*nblockfaces);
        mat_i2 = mat_i + nblockfacessq - 1;
        ind = mat_i : mat_i2;
        blockmattbl.cells  = mattbl.cells(ind);
        blockmattbl.nodes  = mattbl.nodes(ind);
        blockmattbl.faces1 = mattbl.faces1(ind);
        blockmattbl.faces2 = mattbl.faces2(ind);
        sgn1 = 2*(blockmattbl.cells == G.faces.neighbors(blockmattbl.faces1, 1)) - 1;
        sgn2 = 2*(blockmattbl.cells == G.faces.neighbors(blockmattbl.faces2, 1)) - 1;
        nodeMs{i} = nodeMs{i}.*sgn1.*sgn2;
        op = setupTableMapping(blockmattbl, redmattbl, {'nodes', 'faces1', 'faces2'});
        nodeM = nodeM + op*nodeMs{i};
        mat_i = mat_i2 + 1;
        if opt.verbose
            t0 = toc(t0);
            fprintf('reduction of block %d done in %g sec\n', i, t0);
        end
    end

    % Setup matrix
    % First set up facet indices in the redmattbl table
    mat1tbl.nodes  = facenodetbl.nodes;
    mat1tbl.faces1 = facenodetbl.faces;
    mat1tbl.ind    = (1 : facenodetbl.num)';
    op = setupTableMapping(redmattbl, mat1tbl, {'nodes', 'faces1'});
    [colind, rowind] = find(op);
    facetind1 = zeros(redmattbl.num, 1);
    facetind1(rowind) = mat1tbl.ind(colind);
    redmattbl.facetind1 = facetind1;
    
    mat2tbl.nodes  = facenodetbl.nodes;
    mat2tbl.faces2 = facenodetbl.faces;
    mat2tbl.ind    = (1 : facenodetbl.num)';
    op = setupTableMapping(redmattbl, mat2tbl, {'nodes', 'faces2'});
    [colind, rowind] = find(op);
    facetind2 = zeros(redmattbl.num, 1);
    facetind2(rowind) = mat2tbl.ind(colind);
    redmattbl.facetind2 = facetind2;
    
    % Assembly of B
    B = sparse(redmattbl.facetind1, ...
               redmattbl.facetind2, ...
               nodeM, ...
               facenodetbl.num, ...
               facenodetbl.num);
    
    tbls = struct('cellnodefacetbl', cellnodefacetbl, ...
                  'cellfacetbl'    , cellfacetbl    , ...
                  'cellnodetbl'    , cellnodetbl    , ...
                  'facenodetbl'    , facenodetbl);
    
end

function M = node_ip(a, v, N, R, K)
% a : areas of the facets
% v : volume of the corner (for the moment we use volume of cell)
% N : permeability*(facets' normals), corresponds to $\tilde N_c$ in paper of
%     Lipnikov et al (2009)
% R : vector of cell's to facets' centroids, corresponds to $R_c$ in paper of
%     Lipnikov et al (2009)
    
    [U, D, V] = svd(N);
    fnum = size(N, 1);
    dim = size(N, 2);
    
    Dp = D(1 : dim, 1 : dim);
    d = diag(Dp);
    assert(prod(d)>0, 'cannot assemble mpfa, need extra fix'); 
    invd = 1./d;
    H = [diag(invd), zeros(dim, fnum - dim)];
    M = R*V*H*U';
    
    % Add stabilization term S
    S = blkdiag(zeros(dim), eye(fnum - dim));
    U = diag(a)*U;
    t = 6 * sum(diag(K)) / size(K, 2);
    M = M + (t/v)*U*S*U';

end

function M = node_ip2(a, v, N, R, K)
    M = R*inv(N);
end

