function sol = calculateFaceAverages(G, sol)
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
%       map from a reference tetrahedron is constructed, see
%       tetrahedronQuadRule for details. A first-order quadrature rule is
%       used, and \Pi^{\nabla} u is evaluated at each quadrature point
%       through the monomial basis \mathcal{M}_1(K). See [1] for details.
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
%       [1]     - Ø. S. Klemetsdal: 'The virtual element method as a common
%                 framework for finite element and finite difference
%                 methods - Numerical and theoretical analysis'. MA thesis.
%                 Norwegian University of Science and Technology.
%-----------------------------------------------------------------ØSK-2016-

%{
    Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See COPYRIGHT.txt
    for details.
%}

if isempty(sol.faceMoments)

    nF = G.faces.num;
    [m, ~, ~] = retrieveMonomials(2,1);

    [Xq, w, ~, vol] = triangleQuadRule(1);
    nq = size(Xq,1);

    Fc = G.faces.centroids;
    hF = G.faces.diameters;
    aF = G.faces.areas;

    faceAverages = zeros(nF, 1);

    for F = 1:nF
                                %   Node data for cell K.
        nodeNum = G.faces.nodePos(F):G.faces.nodePos(F+1)-1;
        nodes   = G.faces.nodes(nodeNum);
        X       = G.nodes.coords(nodes,:);
        
      
        T = G.faces.localCoords(G.faces.TPos(F):G.faces.TPos(F+1)-1,:);
        X = bsxfun(@minus, X, Fc(F,:))*T;

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
        Xmon = Xhat/hF(F);

        PNstar = G.faces.PNstarT(G.faces.PNstarPos(F):G.faces.PNstarPos(F+1)-1,:)';
        uChi = sol.nodeValues(nodes);
        mVals = m(Xmon)*PNstar*uChi; 
                                %   Multilply by wheights and determinants.
        detAw = rldecode(detA,nq*ones(nTri,1),1).*repmat(vol*w',nTri,1); 
        mVals = bsxfun(@times,mVals,detAw);

        faceAverages(F) = sum(mVals,1)/aF(F);

    end
    
    sol.('faceMoments') = faceAverages;
    
else
    
    warning('Face averages already calculated');
    
end
    
end