function cubature  = makeCubature(disc, type)

    G = disc.G;
    degree = disc.degree + 1;
    issurf = strcmp(type, 'surface');
    dim = G.griddim - issurf;
    basis = dgBasis(dim, degree, 'legendre');
    nDof = basis.nDof;
    psi  = basis.psi;
    
    x = getQuadraturePoints(degree+1, dim);
    
    P = reshape(cell2mat(cellfun(@(p) p(x), psi, 'unif', false)), nDof, nDof)';
    
    switch type
        case 'volume'
            [xq, w, ~, ii, jj, cellNo, faceNo] = makeCellIntegrator(G, (1:G.cells.num)', degree, type);
        case 'surface'
            [xq, w, ~, ii, jj, cellNo, faceNo] = makeFaceIntegrator2(G, (1:G.cells.num)', degree);
    end
    
    W = sparse(ii, jj, w);
    
    xq = disc.transformCoords(xq, cellNo);
    
    I = cellfun(@(p) W*p(xq), psi, 'unif', false);
    rhs = zeros(nDof, G.cells.num);
    for dofNo = 1:nDof
        rhs(dofNo, :) = I{dofNo};
    end
    
    w = reshape(P\rhs, [], 1);
    
    nq0 = nDof;
    if issurf
        [x, nq, cellNo, faceNo] = mapToFace(x,disc);
    else
        x = repmat(x, G.cells.num, 1);
        cellNo = reshape(repmat((1:G.cells.num), nq0, 1), [], 1);
        nq = repmat(nq0, G.cells.num, 1);
    end
   x = disc.transformCoords(x, cellNo, true);
    
%     x = repmat(x, G.cells.num, 1) + rldecode(G.cells.centroids, nq, 1);

    cubature = struct('x', x, 'w', w, 'nq', nq, 'nq0', nq0, 'cellNo', cellNo, 'faceNo', faceNo);
    
    end

    
function [x, nq, cellNo, faceNo] = mapToFace(x, disc)

    G = disc.G;

    faces = G.cells.faces(:,1);
    [bf, bc] = boundaryFaces(G);
    ncbf = sum((1:G.cells.num)' == bc',2);
    ncf = diff(G.cells.facePos);
    faces = faces(~ismember(faces, bf));
    
    nq = size(x,1);
    
    if G.griddim == 3
        edges = G.faces.edges(mcolon(G.faces.edgePos(faces), G.faces.edgePos(faces+1)-1));
        nfe   = diff(G.faces.edgePos);
        nfe   = nfe(faces);
        edges = G.faces.edges(cumsum(nfe));
        nodes = G.edges.nodes(mcolon(G.edges.nodePos(edges), G.edges.nodePos(edges+1)-1));
    else
        nfe = diff(G.cells.facePos);
        nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
    end
    
%     edges = G.faces.edges(cumsum(nfe));
%     nodes = G.edges.nodes(mcolon(G.edges.nodePos(edges), G.edges.nodePos(edges+1)-1));
    
    
    vc = cell(G.griddim-1,1);
    v1    = G.nodes.coords(nodes(2:2:end),:) - G.nodes.coords(nodes(1:2:end-1),:);
    v1    = v1./sqrt(sum(v1.^2,2));
    vc{1} = v1;
    if G.griddim == 3
        n     = G.faces.normals./G.faces.areas(faces);
        v2    = cross(v1, n, 2);
        vc{2} = v2;
    end
    
    v = zeros(G.griddim-1, G.griddim*numel(faces));
    for vNo = 1:G.griddim - 1
        v(vNo, :) = reshape(vc{vNo}', [], 1)';
    end
    
    xx = x*v;
    x = zeros(nq*numel(faces), G.griddim);
    for dNo = 1:G.griddim
        
        x(:, dNo) = reshape(xx(:, dNo:G.griddim:end), [], 1);
        
    end
    
    faceNo = reshape(repmat(faces', nq, 1), [], 1);
    cellNo = rldecode((1:G.cells.num)', (ncf - ncbf)*nq, 1);
    nq = sum((1:G.cells.num)' == cellNo', 2);
    
    xf = disc.transformCoords(G.faces.centroids(faceNo,:), cellNo);
    x = x + xf;

end