function [B, tbls] = robustComputeLocalFluxMimetic(G, rock, opt)

    % Some short aliases 
    nc  = G.cells.num;
    nf  = G.faces.num;
    dim = G.griddim;

    coltbl.coldim = (1 : dim)';
    coltbl.num = dim;
    rowtbl = coltbl;
    rowtbl = replacefield(rowtbl, 'coldim', 'rowdim');

    celltbl.cells = (1 : nc)';
    celltbl.num = nc;

    cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos)); 
    cellfacetbl.faces = G.cells.faces(:, 1);
    cellfacetbl.num   = numel(cellfacetbl.cells);

    faces = rldecode((1 : nf)', diff(G.faces.nodePos)); 
    nodes = G.faces.nodes;
    % We setup the face-node table and it is ordered along ascending node numbers so
    % that we will have a block structure for the nodal scalar product.
    collatefacenode = [nodes, faces];
    collatefacenode = sortrows(collatefacenode, 1);
    facenodetbl = convertArrayToTable(collatefacenode, {'nodes', 'faces'});
    
    [~, cellnodefacetbl] = setupTableMapping(cellfacetbl, facenodetbl, {'faces'});

    % We setup the cell-face-node table, cellnodefacetbl. Each entry determine a
    % unique facet in a corner
    % We order cellnodeface in cell-node-face order. This is node to optimize
    % for-end loop below.
    orderingmat = convertTableToArray(cellnodefacetbl, {'cells', 'nodes', 'faces'});
    orderingmat = sortrows(orderingmat);
    cellnodefacetbl = convertArrayToTable(orderingmat, {'cells', 'nodes', 'faces'});
    
    % We setup the cell-node table, cellnodetbl. Each entry determine
    % a unique corner
    cfn = cellnodefacetbl;
    cellnode = convertTableToArray(cellnodefacetbl, {'nodes', 'cells'});
    cellnode = cellnode(:, [1, 2]);
    cellnode = unique(cellnode, 'rows');
    cellnodetbl = convertArrayToTable(cellnode, {'nodes', 'cells'});
    
    % Nodal scalar product is stored in vector nodeM
    % mattbl is the table which specifies how nodeM is stored: a matrix for
    % each "corner" (cell-node pair).
    duplicate = {'faces', {'faces1', 'faces2'}};
    [~, mattbl] = setupTableMapping(cellnodefacetbl, cellnodefacetbl, {'cells', 'nodes'}, ...
                                         'duplicate', {duplicate});
    % We order mattbl in cell-node-face1-face2 order
    % This is done to optimize for-end loop below
    orderingmat = convertTableToArray(mattbl, {'cells', 'nodes', 'faces1', 'faces2'});
    orderingmat = sortrows(orderingmat);
    mattbl = convertArrayToTable(orderingmat, {'cells', 'nodes', 'faces1', 'faces2'});
    
    nodeM = zeros(mattbl.num, 1);
    
    facetbl.faces = (1 : G.faces.num)';
    facetbl.num   = G.faces.num;
    op = setupTableMapping(facetbl, facenodetbl, {'faces'});
    numnodes = diag(op'*op); % Number of node per face
    op = setupTableMapping(facetbl, cellnodefacetbl, {'faces'});
    numnodes = op*numnodes;
    fno = cellnodefacetbl.faces;
    cno = cellnodefacetbl.cells;
    facetNormals = G.faces.normals(fno, :);
    facetNormals = bsxfun(@ldivide, numnodes, facetNormals);
    sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
    facetNormals = sgn.*facetNormals; % Outward normals with respect to cell
                                      % in cellnodeface.
    [~, cellnodefacecoltbl] = setupTableMapping(cellnodefacetbl, coltbl, []);
    facetNormals = reshape(facetNormals', [], 1);
    
    % Assemble facePermNormals which corresponds to $Kn$ where n are the *outward*
    % normals at the facets.
    [perm, r, c] = permTensor(rock, G.griddim);
    permmat = perm;
    perm = reshape(permmat', [], 1);
    % setup cellcolrow table for the vector perm
    [~, colrowtbl] = setupTableMapping(coltbl, rowtbl, []);
    [~, cellcolrowtbl] = setupTableMapping(colrowtbl, celltbl, []);
    a = convertTableToArray(cellcolrowtbl, {'cells', 'coldim', 'rowdim'});
    a = sortrows(a);
    cellcolrowtbl = convertArrayToTable(a, {'cells', 'coldim', 'rowdim'});
    
    % dispatch perm on cellnodeface
    [~, cellnodefacecolrowtbl] = setupTableMapping(cellcolrowtbl, cellnodefacetbl, ...
                                                                 {'cells'});
    op = setupTableMapping(cellcolrowtbl, cellnodefacecolrowtbl, {'cells', ...
                        'coldim', 'rowdim'});
    perm = op*perm;
    % Multiply perm with facetNormals
    map1 = setupTableMapping(cellnodefacecoltbl, cellnodefacecolrowtbl, {'cells', ...
                        'faces', 'nodes', 'coldim'});
    Kn = perm.*(map1*facetNormals);
    
    map2 = setupTableMapping(cellnodefacecolrowtbl, cellnodefacecoltbl, {'cells', ...
                        'faces', 'nodes', {'rowdim', 'coldim'}});
    Kn = map2*Kn;
    
    % store Kn in matrix form in facePermNormals.
    op = setupTableMapping(cellnodefacetbl, cellnodefacecoltbl, {'cells', ...
                        'faces', 'nodes'});
    ind1 = (1 : cellnodefacetbl.num)';
    ind1 = op*ind1;
    op = setupTableMapping(coltbl, cellnodefacecoltbl, {'coldim'});
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
    
    % set up areas and volumes
    areas = G.faces.areas(fno);
    vols  = G.cells.volumes(cno);
    
    map = setupTableMapping(cellnodetbl, cellnodefacetbl, {'cells', 'nodes'}); 
    nfaces = diag(map'*map);
    nfaces = map*nfaces;
    
    cnf_i = 1; % start indice for the cellnodefacetbl index
    mat_i = 1; % start indice for the mattbl index
    
    for i = 1 : cellnodetbl.num
        
        if opt.verbose
            t0 = tic();
        end
        
        nface = nfaces(cnf_i);
        cnfind = cnf_i : (cnf_i + (nface - 1));
        
        N     = facePermNormals(cnfind, :); 
        R     = cellFacetVec(cnfind, :);
        a     = areas(cnfind); % areas of the faces the facets belong to
        v     = vols(cnf_i); % volume of the current cell
        faces = cellnodefacetbl.faces(cnfind);
        
        cell = cellnodefacetbl.cells(cnf_i);
        node = cellnodefacetbl.nodes(cnf_i);

        K = reshape(permmat(cell, :), [dim, dim]);
        
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

    sgn1 = 2*(mattbl.cells == G.faces.neighbors(mattbl.faces1, 1)) - 1;
    sgn2 = 2*(mattbl.cells == G.faces.neighbors(mattbl.faces2, 1)) - 1;
    nodeM = nodeM.*sgn1.*sgn2;   
    
    % Condensate on nodes (sum up cell contributions for given node).
    [~, redmattbl] = setupTableMapping(facenodetbl, facenodetbl, {'nodes'}, ...
                                               'duplicate', {{'faces', {'faces1', ...
                        'faces2'}}});
    op = setupTableMapping(mattbl, redmattbl, {'nodes', 'faces1', 'faces2'});
    nodeM = op*nodeM;
    
    % Setup matrix
    % First set up facet indices in the redmattbl table
    op = setupTableMapping(facenodetbl, redmattbl, {'nodes', {'faces', ...
                        'faces1'}});
    facesind1 = (1 : facenodetbl.num)';
    facesind1 = op*facesind1;

    op = setupTableMapping(facenodetbl, redmattbl, {'nodes', {'faces', ...
                        'faces2'}});
    facesind2 = (1 : facenodetbl.num)';
    facesind2 = op*facesind2;
    
    % Assembly of B
    B = sparse(facesind1, ...
               facesind2, ...
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
