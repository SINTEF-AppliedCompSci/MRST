function I = polygonInt(G, cells, f, k)
%   Integrates the function f over each cell in cells of grid G, using a
%   quadrature rule of precision k.
%
%   SYNOPSIS:
%       I = polygonInt(G, cells, f, k)
%
%   DESCRIPTION:
%       Approximates the integrals
%           \int_K f \dx
%       over specified cells K of G of using a quadrature rule of
%       precission k. Each cell is trangulated, and a map F from reference
%       triangle with vertices (0,0), (1,0) and (0,1) is constructed. Using
%       that
%
%           \int_K f \dx = |\det(F)|\int_T f(F(y)) \dy,
%
%       the integral can be approximated by the quadrature rule.
%
%   REQUIRED PARAMETERS:
%       G       - MRST grid.
%       cells   - Cells over which to integrate f.
%       f       - Integrand.
%       k       - Precission of quadrature rule.
%
%   RETURNS:
%       I       - Approximated solution to the integral.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}

[Xq, w, ~, vol] = triangleQuadRule(k);

nq = size(Xq,1);

nK = numel(cells);

I = zeros(nK,size(f([0,0]),2));

for i = 1:nK
    
    nodeNum = G.cells.nodePos(cells(i)):G.cells.nodePos(cells(i)+1)-1;
    nodes = G.cells.nodes(nodeNum);
    
    X = G.nodes.coords(nodes,:);
    tri = delaunay(X);
    nTri = size(tri,1);

    bA = X(tri(:,1),:);
    A = X(tri(:,2:end),:) - repmat(bA,2,1);
    A = A(mcolon(1:nTri,2*nTri,nTri),:);
    A = mat2cell(A,2*ones(nTri,1),2);
    D = cellfun(@(X) abs(det(X)), A);
    
    Xhat = cell2mat(cellfun(@(X) Xq*X, A, 'uniformOutput', false));
    Xhat = Xhat + rldecode(bA,nq*ones(nTri,1),1);
    I(i,:) = vol*repmat(w,1,nTri).*(rldecode(D,nq*ones(nTri,1),1))'*f(Xhat);
    
end

end