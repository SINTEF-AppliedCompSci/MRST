function l2Err = l2Error3D(G, sol, u, k)
%--------------------------------------------------------------------------
%   Calculates the square of the L^2-norms over each cell K of the
%   difference between the solution to the Laplace equation and the
%   approximated solution using a kth order VEM.
%
%   SYNOPSIS:
%       l2Err = l2Error(G, sol, u, k)
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
%       G       - 2D MRST grid, with sorted edges, G = sortEdges(G) and
%                 computed VEM geometry, G = computeVEMGeometry(G).
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
[m, ~, ~] = retrieveMonomials(3,k);

if k == 2;
    uFaceMoments = polygonInt3D(G,1:G.faces.num,u,7)./G.faces.areas;
    uCellMoments = polyhedronInt(G,1:nK, u, 7)./G.cells.volumes;
end

[Xq, w, V, vol] = polyhedronQuadRule(7);
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
    Xdiff = mat2cell(X(ii,:) - X(ii+1,:), 3*ones(nTri,1), 3);
    Xb = mat2cell(X(1:4:end,:),ones(nTri,1), 3);

    A = cellfun(@(X) Vdiff\X, Xdiff, 'UniformOutput', false);
    b = cellfun(@(X,Y) X - V(1,:)*Y, Xb, A, 'UniformOutput', false);
    
    D = cellfun(@(X) abs(det(X)), A);
    
    Xhat = cell2mat(cellfun(@(X,y) bsxfun(@plus, Xq*X,y), A, b, ...
                                   'UniformOutput', false));
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