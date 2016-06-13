function  I = polygonInt3D(G, faces, f, k)
%   Integrates the function f over each face in faces of grid 3D G, using a
%   quadrature rule of precision k.
%
%   SYNOPSIS:
%       I = polygonInt(G, cells, f, k)
%
%   DESCRIPTION:
%       Approximates the integrals
%
%           \int_F f \dx
%
%       over specified faces F of G of using a quadrature rule of
%       precission k. Each face is mapped from 3D to 2D,  triangluated, and
%       a map F from reference tringle Tr with vertices (0,0), (1,0) and
%       (0,1) is constructed. Using that
%
%           \int_T f \dx = |\det(Fr)|\int_Tr f(Fr(y)) \dy,
%
%       the integral can be approximated by the quadrature rule.
%
%   REQUIRED PARAMETERS:
%       G       - 3D MRST grid.
%       faces   - faces over which to integrate f.
%       k       - Precision of quadrature rule.
%
%   RETURNS:
%       I       - Approximated solution to the integral.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See COPYRIGHT.txt
   for details.
%}

[Xq, w, ~, vol] = triangleQuadRule(k);
nq = size(Xq,1);
nF = numel(faces);

I = zeros(nF,size(f([0,0,0]),2));

for i = 1:nF
    
    nodeNum = G.faces.nodePos(faces(i)):G.faces.nodePos(faces(i)+1)-1;
    nodes = G.faces.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    
    %   Construct map from polygon to local face coordinates.
    
    vec1 = X(2,:) - X(1,:);
    vec1 = vec1/norm(vec1,2);
    vec2 = cross(G.faces.normals(faces(i),:),vec1);
    vec2 = vec2/norm(vec2,2);
    bT = X(1,:);
    T = [vec1;vec2];
    
    X = (X - repmat(bT,size(X,1),1))*T';
    
    %   Triangulate face, construct map from triangles to reference
    %   triangle.
    
    tri = delaunay(X);
    nTri = size(tri,1);

    bA = X(tri(:,1),:);
    A = X(tri(:,2:end),:) - repmat(bA,2,1);
    A = A(mcolon(1:nTri,2*nTri,nTri),:);
    
    Ad = reshape(A',2,2,[]);
    D2 = abs(squeeze(Ad(1,1,:).*Ad(2,2,:) - Ad(1,2,:).*Ad(2,1,:)));
    
    ii = repmat((1:2*nTri)',2,1);
    jj = 1:2*nTri;
    jj = reshape(jj,2,[])';
    jj = repmat(jj(:)',2,1);
    jj = jj(:);
    
    A2 = sparse(ii, jj, A(:), 2*nTri, 2*nTri);
    
    Xhat2Tmp = repmat(Xq,1,nTri)*A2 + repmat(reshape(bA',1,[]),nq,1);
    Xhat2 = zeros(nq*nTri,2);
    Xhat2(:,1) = reshape(Xhat2Tmp(:,1:2:end),[],1);
    Xhat2(:,2) = reshape(Xhat2Tmp(:,2:2:end),[],1);
    
    %   Evaluate integrals.
    
    XqF = Xhat2*T + repmat(bT,nTri*nq,1);
    I(i,:) = vol*repmat(w,1,nTri).*(rldecode(D2,nq*ones(nTri,1),1))'*f(XqF);
    
end

end