function [B, tbls] = robustComputeLocalFluxMimetic(G, rock, opt)

    % Some short aliases 
    nc  = G.cells.num;
    nf  = G.faces.num;
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
    % We order cellnodeface in cell-node-face order. This is node to optimize
    % for-end loop below.
    orderingmat = [cellnodefacetbl.cells, cellnodefacetbl.nodes, ...
                   cellnodefacetbl.faces];
    orderingmat = sortrows(orderingmat);
    cellnodefacetbl.cells = orderingmat(:, 1);
    cellnodefacetbl.nodes = orderingmat(:, 2);
    cellnodefacetbl.faces = orderingmat(:, 3);
    
    % Some shortcuts
    cno = cellnodefacetbl.cells;
    fno = cellnodefacetbl.faces;
    nno = cellnodefacetbl.nodes;
    
    % We setup the cell-node table, cellnodetbl. Each entry determine
    % a unique corner
    cfn = cellnodefacetbl;
    cellnode = [cfn.nodes, cfn.cells];
    cellnode = unique(cellnode, 'rows');
    cellnodetbl.nodes = cellnode(:, 1);
    cellnodetbl.cells = cellnode(:, 2);
    cellnodetbl.num   = numel(cellnodetbl.nodes);
    % following reordering should not matter ...
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
    
    nodeM = zeros(mattbl.num, 1);
    
    facetbl.faces = (1 : G.faces.num)';
    facetbl.num   = G.faces.num;
    op = setupTableMapping(facetbl, facenodetbl, {'faces'});
    numnodes = diag(op'*op); % Number of node per face
    op = setupTableMapping(facetbl, cellnodefacetbl, {'faces'});
    numnodes = op*numnodes;
    facetNormals = G.faces.normals(fno, :);
    facetNormals = bsxfun(@ldivide, numnodes, facetNormals);
    sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
    facetNormals = sgn.*facetNormals; % Outward normals with respect to cell
                                      % in cellnodeface.

    % Assemble facePermNormals which corresponds to $Kn$ where n are the *outward*
    % normals at the facets.
    cnfnum = cellnodefacetbl.num; %alias
    cellnodefacerowtbl.cells = rldecode(cno, dim*ones(cnfnum, 1));
    cellnodefacerowtbl.faces = rldecode(fno, dim*ones(cnfnum, 1));
    cellnodefacerowtbl.nodes = rldecode(nno, dim*ones(cnfnum, 1));
    cellnodefacerowtbl.rowdim = repmat((1 : dim)', cnfnum, 1);
    cellnodefacerowtbl.num = numel(cellnodefacerowtbl.cells);
    facetNormals = reshape(facetNormals', [], 1);
    
    % Assemble facePermNormals which corresponds to $Kn$ where n are the *outward*
    % normals at the facets.
    [perm, r, c] = permTensor(rock, G.griddim);
    permmat = perm;
    perm = reshape(permmat', [], 1);
    % setup cellcolrow table for the vector perm
    ind = (1 : dim)';
    col = repmat(ind, dim, 1);
    row = rldecode(ind, dim*ones(dim, 1));
    cellcolrowtbl.cells = rldecode((1 : nc)', dim^2*ones(nc, 1));
    cellcolrowtbl.coldim = repmat(col, nc, 1);
    cellcolrowtbl.rowdim = repmat(row, nc, 1);
    cellcolrowtbl.num = numel(cellcolrowtbl.cells);
    % dispatch perm on cellnodeface
    [~, cellnodefacecolrowtbl] = setupTableMapping(cellcolrowtbl, cellnodefacetbl, ...
                                                                 {'cells'});
    op = setupTableMapping(cellcolrowtbl, cellnodefacecolrowtbl, {'cells', ...
                        'coldim', 'rowdim'});
    perm = op*perm;
    % Multiply perm with facetNormals
    map1 = setupTableMapping(cellnodefacerowtbl, cellnodefacecolrowtbl, {'cells', ...
                        'faces', 'nodes', 'rowdim'});
    Kn = perm.*(map1*facetNormals);
    
    cellnodefacecoltbl = replacefield(cellnodefacerowtbl, 'rowdim', 'coldim');
    map2 = setupTableMapping(cellnodefacecolrowtbl, cellnodefacecoltbl, {'cells', ...
                        'faces', 'nodes', 'coldim'});
    Kn = map2*Kn;
    
    % store Kn in matrix form in facePermNormals.
    op = setupTableMapping(cellnodefacecoltbl, cellnodefacetbl, {'cells', ...
                        'faces', 'nodes'});
    [ind1, ind2] = find(op); 
    ind2 = cellnodefacecoltbl.coldim(ind2);
    facePermNormals = sparse(ind1, ind2, Kn, cnfnum, dim);
    
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
        locM = node_ip2(a, v, ...
                        full(N), ...
                        full(R), ...
                        K);
        locM = reshape(locM', [], 1);
        
        matind = mat_i : (mat_i + (nface*nface - 1));
        nodeM(matind) = nodeM(matind) + locM;
        
        % increment start indices
        cnf_i = cnf_i + nface;
        mat_i = mat_i + nface*nface;        
        
        if opt.verbose
            t0 = toc(t0);
            fprintf('assembly cellnode %d / %d took %g sec\n', i, cellnodetbl.num, ...
                    t0);
        end
    end

    sgn1 = 2*(mattbl.cells == G.faces.neighbors(mattbl.faces1, 1)) - 1;
    sgn2 = 2*(mattbl.cells == G.faces.neighbors(mattbl.faces2, 1)) - 1;
    nodeM = nodeM.*sgn1.*sgn2;   
    
    % Condensate on nodes (sum up cell contributions for give node).
    op = setupTableMapping(facenodetbl, facenodetbl, {'nodes'});
    [colind, rowind] = find(op);
    redmattbl.nodes  = facenodetbl.nodes(colind);
    redmattbl.faces1 = facenodetbl.faces(colind);
    redmattbl.faces2 = facenodetbl.faces(rowind);
    redmattbl.num    = numel(redmattbl.nodes);
    
    op = setupTableMapping(mattbl, redmattbl, {'nodes', 'faces1', 'faces2'});
    nodeM = op*nodeM;
    
    % Setup matrix
    % First set up facet indices in the redmattbl table
    clear tmptbl;
    tmptbl.nodes  = redmattbl.nodes;
    tmptbl.faces  = redmattbl.faces1;
    tmptbl.faces2 = redmattbl.faces2;
    tmptbl.num = numel(tmptbl.nodes);
    op = setupTableMapping(tmptbl, facenodetbl, {'nodes', 'faces'});
    [colind, rowind] = find(op);
    facesind1(rowind) = colind;

    clear tmptbl;
    tmptbl.nodes  = redmattbl.nodes;
    tmptbl.faces1 = redmattbl.faces1;
    tmptbl.faces  = redmattbl.faces2;
    tmptbl.num  = numel(tmptbl.nodes);
    op = setupTableMapping(tmptbl, facenodetbl, {'nodes', 'faces'});
    [colind, rowind] = find(op);
    facesind2(rowind) = colind;
    
    
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
