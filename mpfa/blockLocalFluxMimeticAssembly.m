function [B, tbls] = blockLocalFluxMimeticAssembly(G, rock, nodes, opt)

    nodetbl.nodes = nodes;
    nodetbl.num = numel(nodes);
    
    nc  = G.cells.num;
    nf  = G.faces.num;
    dim = G.griddim;

    % Setup cellnodetbl for *whole* grid (this could be moved out of this function)
    cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos)); 
    cellfacetbl.faces = G.cells.faces(:, 1);
    cellfacetbl.num   = numel(cellfacetbl.cells);
    facenodetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos)); 
    facenodetbl.nodes = G.faces.nodes;    
    facenodetbl.num = numel(facenodetbl.faces);    
    [~, cellnodefacetbl] = setupTableMapping(cellfacetbl, facenodetbl, ...
                                                          {'faces'});
    cellnodetbl = projTable(cellnodefacetbl, {'cells', 'nodes'});
    
    % Restrict cellnodetbl to the nodes blocks (nodes that are sent as argument)
    [~, cellnodetbl] = setupTableMapping(cellnodetbl, nodetbl, {'nodes'});
    % Restrict cellnodefacetbl to the node blocks. Each entry of cellnodefacetbl
    % determine a unique facet in a corner
    [~, cellnodefacetbl] = setupTableMapping(cellnodetbl, cellnodefacetbl, ...
                                                          {'nodes', 'cells'});
    % We order cellnodeface in cell-node-face order. This is done to optimize the
    % for-end loop below.
    cellnodeface = convertTableToArray(cellnodefacetbl, {'cells', 'nodes', 'faces'});
    cellnodeface = sortrows(cellnodeface);
    cellnodefacetbl = convertArrayToTable(cellnodeface, {'cells', 'nodes', 'faces'});
    cellnodefacetbl = addLocInd(cellnodefacetbl, 'cnfind');
    
    % We set up the facenode table. It is reordered so that we obtain a diagonal
    % matrix in the local indices.
    facenodetbl = projTable(cellnodefacetbl, {'nodes', 'faces'});
    facenode = convertTableToArray(facenodetbl, {'nodes', 'faces'});
    facenode = sortrows(facenode);
    facenodetbl = convertArrayToTable(facenode, {'nodes', 'faces'});
    facenodetbl = addLocInd(facenodetbl, 'fnind');
    
    % col and row tables for matrices in the spatial dimension (dim).
    coltbl.coldim = (1 : dim)';
    coltbl.num = dim;
    rowtbl = coltbl;
    rowtbl = replacefield(rowtbl, 'coldim', 'rowdim');

    % We set up cell table (only cells that belong to the block).
    celltbl = projTable(cellnodetbl, {'cells'});
    
    % Nodal scalar product is stored in vector nodeM mattbl is the table which
    % specifies how nodeM is stored: a matrix for each "corner" (cell-node
    % pair).
    duplicate = {'faces', {'faces1', 'faces2'}};
    [~, mattbl] = setupTableMapping(cellnodefacetbl, cellnodefacetbl, {'cells', 'nodes'}, ...
                                         'duplicate', {duplicate});
    % We order mattbl in cell-node-face1-face2 order. This is done to optimize
    % for-end loop below
    orderingmat = convertTableToArray(mattbl, {'cells', 'nodes', 'faces1', 'faces2'});
    orderingmat = sortrows(orderingmat);
    mattbl = convertArrayToTable(orderingmat, {'cells', 'nodes', 'faces1', 'faces2'});
    
    nodeM = zeros(mattbl.num, 1);

    if opt.verbose
        fprintf('assemble facet normals ...\n');
    end
    
    % We compute the product of the permeability tensor K and the normal at each
    % node-face. The normal is weighted by area divided number of nodes the face
    % is dealt with.
    fno = cellnodefacetbl.faces;
    cno = cellnodefacetbl.cells;
    numnodes = double(diff(G.faces.nodePos));
    numnodes = numnodes(fno);
    facetNormals = G.faces.normals(fno, :);
    facetNormals = bsxfun(@ldivide, numnodes, facetNormals);
    
    sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
    facetNormals = sgn.*facetNormals; % Outward normals with respect to cell
                                      % in cellnodeface.
    
    [~, cellnodefacecoltbl] = setupTableMapping(cellnodefacetbl, coltbl, []);
    a = convertTableToArray(cellnodefacecoltbl, {'cells', 'nodes', 'faces', ...
                        'coldim', 'cnfind'});
    a = sortrows(a);
    cellnodefacecoltbl = convertArrayToTable(a, {'cells', 'nodes', 'faces', ...
                        'coldim', 'cnfind'});
    facetNormals = reshape(facetNormals', [], 1);
    
    if opt.verbose
        fprintf('assemble facet K*normals ...\n');
    end
    % Assemble facePermNormals which corresponds to $Kn$ where n are the *outward*
    % (weighted) normals at the facets
    [perm, r, c] = permTensor(rock, G.griddim);
    permmat = perm(celltbl.cells, :);
    perm = reshape(permmat', [], 1);
    % Setup cellcolrow table for the vector perm
    [~, colrowtbl] = setupTableMapping(coltbl, rowtbl, []);
    [~, cellcolrowtbl] = setupTableMapping(colrowtbl, celltbl, []);
    a = convertTableToArray(cellcolrowtbl, {'cells', 'coldim', 'rowdim'});
    a = sortrows(a);
    cellcolrowtbl = convertArrayToTable(a, {'cells', 'coldim', 'rowdim'});
    cellcolrowtbl = addLocInd(cellcolrowtbl, 'ccrind');
    
    % Dispatch permeability tensor on cellnodeface
    [~, cellnodefacecolrowtbl] = setupTableMapping(cellcolrowtbl, cellnodefacetbl, ...
                                                                 {'cells'});
    op = reduceDispatchMapping(cellcolrowtbl, cellnodefacecolrowtbl, 'ccrind');
    perm = op*perm;
    % Multiply permeability tensor with facetNormals
    map1 = setupTableMapping(cellnodefacecoltbl, cellnodefacecolrowtbl, ...
                                           {'cnfind', 'coldim'}, 'fastunstable', ...
                                           true);
    Kn = perm.*(map1*facetNormals);
    
    map2 = setupTableMapping(cellnodefacecolrowtbl, cellnodefacecoltbl, ...
                                           {'cnfind', {'rowdim', 'coldim'}}, ...
                                           'fastunstable', true);
    Kn = map2*Kn;
    
    % store Kn in matrix form in facePermNormals.
    op = reduceDispatchMapping(cellnodefacetbl, cellnodefacecoltbl, 'cnfind');
    ind1 = cellnodefacetbl.cnfind;
    ind1 = op*ind1;
    op = reduceDispatchMapping(coltbl, cellnodefacecoltbl, 'coldim');
    ind2 = (1 : coltbl.num)';
    ind2 = op*ind2;    
    facePermNormals = sparse(ind1, ind2, Kn, cellnodefacetbl.num, coltbl.num);
    
    % Some shortcuts
    cno = cellnodefacetbl.cells;
    fno = cellnodefacetbl.faces;
    nno = cellnodefacetbl.nodes;
    % Default option (opt.eta = 0): Use original face centroids and cell centroids,
    % NOT actual subface centroids. This corresponds to an MPFA method
    % (O-method) R = G.faces.centroids(fno,:) - G.cells.centroids(cno,:);
    cellFacetVec = G.faces.centroids(fno,:) - G.cells.centroids(cno,:) + ...
        opt.eta*(G.nodes.coords(nno,:) - G.faces.centroids(fno,:));
    
    % Set up areas and volumes
    areas = G.faces.areas(fno);
    vols  = G.cells.volumes(cno);
    
    map = setupTableMapping(cellnodetbl, cellnodefacetbl, {'cells', 'nodes'}); 
    nfaces = diag(map'*map);
    nfaces = map*nfaces;
    
    cnf_i = 1; % start indice for the cellnodefacetbl index
    mat_i = 1; % start indice for the mattbl index
    
    for i = 1 : cellnodetbl.num
        
        % if opt.verbose
            % t0 = tic();
        % end
        
        nface = nfaces(cnf_i);
        cnfind = cnf_i : (cnf_i + (nface - 1));
        
        N     = facePermNormals(cnfind, :); 
        R     = cellFacetVec(cnfind, :);
        a     = areas(cnfind); % areas of the faces the facets belong to
        v     = vols(cnf_i); % volume of the current cell
        faces = cellnodefacetbl.faces(cnfind);
        
        cellno = cellnodefacetbl.cells(cnf_i);
        K = reshape(permmat(cellno, :), [dim, dim]);
        
        % Assemble local nodal scalar product ( function node_ip2 below handle case when
        % N is invertible)
        locM = node_ip(a, v, ...
                       full(N), ...
                       full(R), ...
                       K);
        locM = reshape(locM', [], 1);
        
        matind = mat_i : (mat_i + (nface*nface - 1));
        nodeM(matind) = nodeM(matind) + locM;
        
        % increment start indices
        cnf_i = cnf_i + nface;
        mat_i = mat_i + nface*nface;        
        
        % if opt.verbose 
            % t0 = toc(t0);
            % fprintf('assembly cellnode %d / %d took %g sec\n', i, cellnodetbl.num, ...
                    % t0);
        % end
    end
    
    if opt.verbose
        fprintf('include face sign ...\n');
    end
    sgn1 = 2*(mattbl.cells == G.faces.neighbors(mattbl.faces1, 1)) - 1;
    sgn2 = 2*(mattbl.cells == G.faces.neighbors(mattbl.faces2, 1)) - 1;
    nodeM = nodeM.*sgn1.*sgn2;   
    
    if opt.verbose
        fprintf('Condensate on nodes ...\n');
    end    
    % Condensate on nodes (sum up cell contributions for given node).
    [~, redmattbl] = setupTableMapping(facenodetbl, facenodetbl, {'nodes'}, ...
                                               'duplicate', {{'faces', {'faces1', ...
                        'faces2'}}, {'fnind', {'fnind1', 'fnind2'}}});
    op = setupTableMapping(mattbl, redmattbl, {'nodes', 'faces1', 'faces2'});
    nodeM = op*nodeM;
    
    if opt.verbose
        fprintf('Set up matrix ...\n');
    end    
    
    % Assembly of B
    B = sparse(redmattbl.fnind1, ...
               redmattbl.fnind2, ...
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
    
    k = 6 * sum(diag(K)) / size(K, 2);
    scalfact = k*sum(a)/numel(a);
    N = (1/scalfact)*N;
    [U, D, V] = svd(N);
    D = scalfact*D;
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
    t = 6 * v^(5/3) * sum(diag(K)) / size(K, 2);
    regM = (1/t)*U*S*U';
    % fprintf('norm main: %g, norm reg: %g\n', norm(M), norm(regM));
    M = M + regM;
    % fprintf('condition number: %g\n', condest(M));
end

function M = node_ip2(a, v, N, R, K)
    M = R*inv(N);
end
