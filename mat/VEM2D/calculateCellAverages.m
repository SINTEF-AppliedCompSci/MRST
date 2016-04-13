function sol = calculateCellAverages(G, sol)
%--------------------------------------------------------------------------
%   Calculates cell averages of solution sol obtained from a 1st order VEM.
%
%   SYNOPSIS:
%       sol = calculateCellAverages(G, sol)
%
%   DESCRIPTION:
%       Evaluates the integrals
%   
%           |K|^{-1}\int_{K} u \dx = |K|^{-1}\int_{K} \Pi^{\nabla} u \dx 
% 
%       for each cell K of the solution u. Each cell is trangulated, and a
%       map from reference triangle with vertices (0,0), (1,0) and (0,1)
%       is constructed. A first-order quadrature rule is used, and 
%       \Pi^{\nabla} u is evaluated at each quadrature point through the
%       monomial basis \mathcal{M}_1(K). See [1] for details.
%
%   REQUIRED PARAMETERS:
%       G       - 2D MRST grid, with sorted edges, G = sortEdges(G),
%                 computed VEM geometry, G = computeVEMGeometry(G), and
%                 projectors \Pi^\nabla.
%       sol     - Solution obtained from 1st order VEM.
%
%   RETURNS:
%       sol     - Updates solution structure with cell averages
%                 (cellMoments).
%
%   REFERENCES:
%       [1]     - Thesis title.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}

if isempty(sol.cellMoments)

    nK = G.cells.num;
    [m, ~, ~] = retrieveMonomials(2,1);

    [Xq, w, ~, vol] = triangleQuadRule(1);
    nq = size(Xq,1);

    Kc = G.cells.centroids;
    hK = G.cells.diameters;
    aK = G.cells.volumes;

    cellAverages = zeros(nK, 1);

    for K = 1:nK
                                %   Node data for cell K.
        nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
        nodes   = G.cells.nodes(nodeNum);
        X = G.nodes.coords(nodes,:);

        tri = delaunay(X);
        nTri = size(tri,1);
                                %   Construct map from refrence triangle to
                                %   triangles in triangulation.
        bA = X(tri(:,1),:);
        A = X(tri(:,2:end),:) - repmat(bA,2,1);
        A = A(mcolon(1:nTri,2*nTri,nTri),:);
        A = mat2cell(A,2*ones(nTri,1),2);
        detA = cellfun(@(X) abs(det(X)), A);

                                %   Map from reference triangle to triangels.
        Xhat = cell2mat(cellfun(@(X) Xq*X, A, 'uniformOutput', false));
        Xhat = Xhat + rldecode(bA,nq*ones(nTri,1),1);
                                %   Scale coordinates for use in 2D monomials.
        Xmon = (Xhat - repmat(Kc(K,:),nq*nTri,1))/hK(K);

        PNstar = G.PNstarT(G.PNstarPos(K):G.PNstarPos(K+1)-1,:)';
        uChi = sol.nodeValues(nodes);
        mVals = m(Xmon)*PNstar*uChi; 
                                %   Multilply by wheights and determinants.
        detAw = rldecode(detA,nq*ones(nTri,1),1).*repmat(vol*w',nTri,1); 
        mVals = bsxfun(@times,mVals,detAw);

        cellAverages(K) = sum(mVals,1)/aK(K);

    end
    
    sol.('cellMoments') = cellAverages;
    
else
    
    warning('Cell averages already calculated');
    
end
    
end
    