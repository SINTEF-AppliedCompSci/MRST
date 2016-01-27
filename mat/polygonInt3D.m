    function  I = polygonInt3D(G, faces, f)

                            %   Gauss-Lobatto quadrature point and
                            %   wheights for refenrence triangle.
Xq = [0.0, 0.0; 0.0, 1.0; 0.5, 0.0; 0.5, 0.5; 0.0, 0.5; 0.5, 0.25];
w = [1/36, 1/36, 1/18, 1/18, 1/9, 2/9];
                            %   For each triangle t, evaluate integral.

nq = size(Xq,1);

nF = numel(faces);

I = zeros(nF,size(f([0,0,0]),2));

for i = 1:nF
    
    nodeNum = G.faces.nodePos(faces(i)):G.faces.nodePos(faces(i)+1)-1;
    nodes = G.faces.nodes(nodeNum);
    nNF = size(nodes,1);
    X = G.nodes.coords(nodes,:);
    
    vec1 = X(2,:) - X(1,:);
    vec1 = vec1/norm(vec1,2);
    vec2 = cross(G.faces.normals(faces(i),:),vec1);
    vec2 = vec2/norm(vec2,2);
    bT = X(1,:);
    T = [vec1;vec2];
    
                            %   Map from Polygon to face
    X = (X - repmat(bT,size(X,1),1))*T';
    
                            %   Triangulate face
    tri = delaunay(X);
    nTri = size(tri,1);
    
                            %   Construct map from polygon to local face
                            %   coordinates.
    bA = X(tri(:,1),:);
    A = X(tri(:,2:end),:) - repmat(bA,2,1);
    A = A(mcolon(1:nTri,2*nTri,nTri),:);
    A = mat2cell(A,2*ones(nTri,1),2);
    D = cellfun(@(X) abs(det(X)), A);
    
    Xhat = cell2mat(cellfun(@(X) Xq*X, A, 'uniformOutput', false));
    Xhat = Xhat + rldecode(bA,nq*ones(nTri,1),1);
    XqF = Xhat*T + repmat(bT,nTri*nq,1);
    I(i,:) = repmat(w,1,nTri).*(rldecode(D,nq*ones(nTri,1),1))'*f(XqF);
    
end

end