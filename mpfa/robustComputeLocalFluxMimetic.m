function [B, tbls] = robustComputeLocalFluxMimetic(G, rock, opt)

    % Some short aliases 
    nc  = G.cells.num;
    nf  = G.faces.num;
    dim = G.griddim;

    coltbl.coldim = (1 : dim)';
    coltbl.num = dim;
    rowtbl = coltbl;
    rowtbl = replacefield(rowtbl, {'coldim', 'rowdim'});

    celltbl.cells = (1 : nc)';
    celltbl.num = nc;

    cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos)); 
    cellfacetbl.faces = G.cells.faces(:, 1);
    cellfacetbl.num   = numel(cellfacetbl.cells);

    facenodetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos)); 
    facenodetbl.nodes = G.faces.nodes;
    facenodetbl.num   = numel(facenodetbl.faces);
    % We setup the face-node table and it is ordered along ascending node numbers so
    % that we will have a block structure for the nodal scalar product.
    facenodetbl = sortTable(facenodetbl, {'nodes', 'faces'});
    
    cellnodefacetbl = crossTable(cellfacetbl, facenodetbl, {'faces'});

    % We setup the cell-face-node table, cellnodefacetbl. Each entry determine a
    % unique facet in a corner
    % We order cellnodeface in cell-node-face order. This is node to optimize
    % for-end loop below.
    cellnodefacetbl = sortTable(cellnodefacetbl, {'cells', 'nodes', 'faces'});
    cellnodefacetbl = addLocInd(cellnodefacetbl, 'cnfind');    
    
    % We setup the cell-node table, cellnodetbl. Each entry determine a unique
    % corner
    cellnodetbl = projTable(cellnodefacetbl, {'nodes', 'cells'});
    
    % Nodal scalar product is stored in vector nodeM
    % mattbl is the table which specifies how nodeM is stored: a matrix for
    % each "corner" (cell-node pair).
    crossextend = {'faces', {'faces1', 'faces2'}};
    mattbl = crossTable(cellnodefacetbl, cellnodefacetbl, {'cells', 'nodes'}, ...
                        'crossextend', {crossextend});
    % We order mattbl in cell-node-face1-face2 order
    % This is done to optimize for-end loop below
    mattbl = sortTable(mattbl, {'cells', 'nodes', 'faces1', 'faces2'});
    
    nodeM = zeros(mattbl.num, 1);

    if opt.verbose
        fprintf('assemble facet normals ...\n');
    end
    
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
    cellnodefacecoltbl = sortTable(cellnodefacecoltbl, {'cells', 'nodes', ...
                        'faces', 'coldim', 'cnfind'});
    facetNormals = reshape(facetNormals', [], 1);
    
    if opt.verbose
        fprintf('assemble facet K*normals ...\n');
    end
    % Assemble facePermNormals which corresponds to $Kn$ where n are the *outward*
    % normals at the facets.
    [perm, r, c] = permTensor(rock, G.griddim);
    permmat = perm;
    perm = reshape(permmat', [], 1);
    % setup cellcolrow table for the vector perm
    colrowtbl = crossTable(coltbl, rowtbl, []);
    cellcolrowtbl = crossTable(colrowtbl, celltbl, {});
    cellcolrowtbl = sortTable(cellcolrowtbl, {'cells', 'coldim', 'rowdim'});
    cellcolrowtbl = addLocInd(cellcolrowtbl, 'ccrind');
    
    % Multiply perm with facetNormals
    prod = TensorProd();
    prod.tbl1 = cellcolrowtbl;
    prod.tbl2 = cellnodefacecoltbl;
    prod.replacefds1 = {{'coldim', 'temp'}, {'rowdim', 'coldim'}, {'temp', 'rowdim'}};
    prod.replacefds2 = {'coldim', 'rowdim'};
    prod.mergefds = {'cells'};
    prod.reducefds = {'rowdim'};
    prod.prodtbl = cellnodefacecoltbl;
    prod = prod.setup();
    
    Kn = prod.evalProd(perm, facetNormals);
   
    % store Kn in matrix form in facePermNormals.
    ind1 = (1 : cellnodefacetbl.num)';
    ind1 = tblmap(ind1, cellnodefacetbl, cellnodefacecoltbl, {'cells', 'nodes', 'faces'});
    ind2 = (1 : coltbl.num)';
    ind2 = tblmap(ind2, coltbl, cellnodefacecoltbl, {'coldim'});
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
    
    % number of faces per cell-nodes.
    nfaces = tblmap(ones(cellnodefacetbl.num, 1), cellnodefacetbl, cellnodetbl, {'cells', 'nodes'}); 
    % we setup nfaces indexed along cellnodefacetbl
    nfaces = tblmap(nfaces, cellnodetbl, cellnodefacetbl, {'cells', 'nodes'}); 
    
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
        
        % Assemble local nodal scalar product 
        switch opt.ip_compmethod
          case 'general'
            locM = node_ip(a, v, full(N), full(R), K);
          case 'nicecorner'
            locM = node_ip2(full(N), full(R));
          case 'directinverse'
            volE = G.cells.volumes(cellno);
            locM = node_ip3(full(N), K, G.griddim, volE);
          otherwise
            error('ip_compmethod not recognized');
        end
        
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
    redmattbl = projTable(mattbl, {'nodes', 'faces1', 'faces2'});
    nodeM = tblmap(nodeM, mattbl, redmattbl, {'nodes', 'faces1', 'faces2'});
    
    if opt.verbose
        fprintf('Set up matrix ...\n');
    end    
    % Setup matrix
    % First set up facet indices in the redmattbl table
    fnind = (1 : facenodetbl.num)';
    facesind1 = tblmap(fnind, facenodetbl, redmattbl, {'nodes', {'faces', ...
                        'faces1'}});
    facesind2 = tblmap(fnind, facenodetbl, redmattbl, {'nodes', {'faces', ...
                        'faces2'}});
    
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

function M = node_ip2(N, R)
    M = R*inv(N);
end

function M = node_ip3(N, K, d, volE)
% volE : volume of cell
% d    : Spatial dimension (2 or 3)
% K    : permeability tensor    
    invN = inv(N);
    switch d
      case 2
        mE = 3;
      case 3
        mE = 4;
      otherwise
        error('wrong spatial dimension (should be equal to 2 or 3)');
    end
    
    M = (volE/mE)*(invN'*K*invN);
end

