function l2Err = l2Error3D(G, sol, u, k)
%--------------------------------------------------------------------------
%   Calculates the square of the L^2-norms over each cell K of the
%   difference between the solution to the Laplace equation and the
%   approximated solution using a kth order VEM.
%
%   SYNOPSIS:
%       l2Err = l2Error3D(G, sol, u, k)
%
%   DESCRIPTION:
%       Evaluates the integrals
%   
%           \int_{K} |u-U|^2 \dx
% 
%       for each cell K, where u is the analytical solution to the Laplace
%       equation
%
%       -\Delta u = f,                                                (1)  
%
%       and U is the approximated solution to the same equation, using a
%       kth order VEM.
%
%   REQUIRED PARAMETERS:
%       G       - 2D MRST grid, with sorted edges, G = sortEdges(G),
%                 computed VEM geometry, G = computeVEMGeometry(G), and
%                 cell projectors.
%       sol     - Solution struct obtained from kth order VEM.
%       u       - Analytical solution to (1).
%       k       - Method order. Can be lower than order of VEM used to
%                 obtain sol, but will then use only nodeValues.
%
%   RETURNS:
%       l2Err   - Square of L^2(K)-norm of difference between analytical
%                 solution and approximated solution for (1), for each cell
%                 K.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}

nK = G.cells.num;
nF = G.faces.num;
[m, ~, ~] = retrieveMonomials(3,k);

if k == 2;
    uFaceMoments = polygonInt3D(G,1:nF,u,7)./G.faces.areas;
    uCellMoments = polyhedronInt(G,1:nK, u, 7)./G.cells.volumes;
end

[Xq, w, V, vol] = tetrahedronQuadRule(7);
nq = size(Xq,1);
Vdiff = V(1:end-1,:) - V(2:end,:);

Kc = G.cells.centroids;
hK = G.cells.diameters;

l2Err = zeros(nK,1);

for K = 1:nK
    
    nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);

    if k == 2
        edgeNum = G.cells.edgePos(K):G.cells.edgePos(K+1)-1;
        edges = G.cells.edges(edgeNum);
        faceNum = G.cells.facePos(K):G.cells.facePos(K+1)-1;
        faces = G.cells.faces(faceNum);
    end
    
    tri = delaunay(X);
    nTri = size(tri,1);
    
    X = X(reshape(tri',[],1),:);
    ii = mcolon(1:4:4*nTri, 3:4:4*nTri);

    Xdiff = X(ii,:) - X(ii+1,:);
    
%     Xdiff = mat2cell(X(ii,:) - X(ii+1,:), 3*ones(nTri,1), 3);
    
    ii = repmat((1:3)',3*nTri,1);
    add = repmat((0:1:nTri-1)*3,9,1);
    ii = ii + add(:);
    jj = repmat((1:3), 3, 1);
    jj = repmat(jj(:), nTri,1);
    jj = jj(:) + add(:);
    Vd = sparse(ii,jj,repmat(Vdiff(:),nTri,1), 3*nTri, 3*nTri);
    Xdiff = Xdiff';
    Xd = sparse(jj,ii,Xdiff(:), 3*nTri, 3*nTri);
    
%     Xb = mat2cell(X(1:4:end,:),ones(nTri,1), 3);

    A = Vd\Xd;

    Xb = reshape(X(1:4:end,:)',1,[]);

    b = Xb - repmat(V(1,:),1,nTri)*A;

%     A = cellfun(@(X) Vdiff\X, Xdiff, 'UniformOutput', false);
%     b = cellfun(@(X,Y) X - V(1,:)*Y, Xb, A, 'UniformOutput', false);
    
%     D = cellfun(@(X) abs(det(X)), A);
    [n,~] = size(A);
    
    main = ones(n,1);
    sub1 = main;
    sub1(3:3:end) = 0;
    sub2 = zeros(n,1);
    sub2(1:3:end) = 1;
    Id = spdiags([sub2, sub1, main,sub1(end:-1:1),sub2(end:-1:1)], -2:1:2, n,n);

    Ad = full(A(Id == 1));
    Ad = reshape(Ad,3,3,[]);
    D = squeeze(Ad(1,1,:).*Ad(2,2,:).*Ad(3,3,:) + ...
                Ad(2,1,:).*Ad(3,2,:).*Ad(1,3,:) + ...
                Ad(3,1,:).*Ad(1,2,:).*Ad(2,3,:) - ...
                Ad(3,1,:).*Ad(2,2,:).*Ad(1,3,:) - ...
                Ad(2,1,:).*Ad(1,2,:).*Ad(3,3,:) - ...
                Ad(1,1,:).*Ad(3,2,:).*Ad(2,3,:));
    D = abs(D);
    
    XhatTmp = bsxfun(@plus, repmat(Xq, 1, nTri)*A, b);
    Xhat = zeros(nq*nTri,3);
    Xhat(:,1) = reshape(XhatTmp(:,1:3:3*nTri),[],1);
    Xhat(:,2) = reshape(XhatTmp(:,2:3:3*nTri),[],1);
    Xhat(:,3) = reshape(XhatTmp(:,3:3:3*nTri),[],1);
    
%     Xhat = cell2mat(cellfun(@(X,y) bsxfun(@plus, Xq*X,y), A, b, ...
%                                    'UniformOutput', false));
    Xmon = bsxfun(@minus, Xhat, Kc(K,:))/hK(K);
    
    PNstarNum = G.cells.PNstarPos(K):G.cells.PNstarPos(K+1)-1;
    PNstar = G.cells.PNstarT(PNstarNum,:)';
    
    if k == 1
        UChi = sol.nodeValues(nodes);
    elseif k == 2
        UChi = [sol.nodeValues(nodes) ; sol.edgeValues(edges); ...
                sol.faceMoments(faces); sol.cellMoments(K)   ];
    end
                               
    vals= (m(Xmon)*PNstar*UChi - u(Xhat)).^2;
    detAw = rldecode(D,nq*ones(nTri,1),1).*repmat(vol*w',nTri,1); 
    vals = bsxfun(@times,vals,detAw);
    
    l2Err(K) = sum(vals,1);
        
end

end