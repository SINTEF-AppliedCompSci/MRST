function l2Err = l2Error2D(G, sol, u, k)
%   Calculates the square of the L^2-norms over each cell K of the
%   difference between the solution to the Poisson equation and the
%   approximated solution using a kth order VEM.
%
%   SYNOPSIS:
%       l2Err = l2Error(G, sol, u, k)
%
%   DESCRIPTION:
%       Approximates the integrals
%   
%           \int_{K} |u-U|^2 \dx
%           \approx \int_{K}|u-\Pi^\nabla U|^2_{0,K} \dx
% 
%       for each cell K, where u is the analytical solution to the Poisson
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
%       l2Err   - Square of L^2(K) norm of difference between analytical
%                 solution and approximated solution for (1), for each cell
%                 K.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See COPYRIGHT.txt
   for details.
%}

assert(k == 1 || k == 2, 'Method order k must be 1 or 2.')

nK = G.cells.num;
[m, ~, ~] = retrieveMonomials(2,k);

if k == 2;
    uCellMoments = polygonInt(G,1:nK, u, 7)./G.cells.volumes;
end

[Xq, w, ~, vol] = triangleQuadRule(7);
nq = size(Xq,1);

Kc = G.cells.centroids;
hK = G.cells.diameters;

l2Err = zeros(nK,1);

for K = 1:nK

    nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
    nodes   = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    if k == 2
        faceNum = G.cells.facePos(K):G.cells.facePos(K+1)-1;
        faces = G.cells.faces(faceNum);
    end

    tri = delaunay(X);
    nTri = size(tri,1);

    bA = X(tri(:,1),:);
    A = X(tri(:,2:end),:) - repmat(bA,2,1);
    A = A(mcolon(1:nTri,2*nTri,nTri),:);
    A = mat2cell(A,2*ones(nTri,1),2);
    detA = cellfun(@(X) abs(det(X)), A);
    
    Xhat = cell2mat(cellfun(@(X) Xq*X, A, 'uniformOutput', false));
    Xhat = Xhat + rldecode(bA,nq*ones(nTri,1),1);
    Xmon = (Xhat - repmat(Kc(K,:),nq*nTri,1))/hK(K);

    PNstar = G.cells.PNstarT(G.cells.PNstarPos(K):G.cells.PNstarPos(K+1)-1,:)';
    
    if k == 1
        UChi = sol.nodeValues(nodes);
    elseif k == 2
        UChi = [sol.nodeValues(nodes); sol.edgeValues(faces); sol.cellMoments(K)];
    end
    
    mVals = (m(Xmon)*PNstar*UChi - u(Xhat)).^2;
    detAw = rldecode(detA,nq*ones(nTri,1),1).*repmat(vol*w',nTri,1); 
    mVals = bsxfun(@times,mVals,detAw);
    
    l2Err(K) = sum(mVals,1);
        
end

end