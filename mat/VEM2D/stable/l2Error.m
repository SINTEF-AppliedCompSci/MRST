function l2Err = l2Error(G, sol, u, k)

nK = G.cells.num;
[m, ~, ~] = retrieveMonomials(k);

if k == 2;
    uCellMoments = polygonInt_v2(G,1:nK, u, 7)./G.cells.volumes;
end

[Xq, w, ~, vol] = triangleQuadRule(7);
nq = size(Xq,1);

Kc = G.cells.centroids;
hK = G.cells.diameters;
aK = G.cells.volumes;

l2Err = zeros(nK,1);

for K = 1:nK
                            %   Node data for cell K.
    nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
    nodes   = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    if k == 2
        faceNum = G.cells.facePos(K):G.cells.facePos(K+1)-1;
        faces = G.cells.faces(faceNum);
        Ec = G.faces.centroids(faces,:);
    end

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
    
    if k == 1
        UChi = sol.nodeValues(nodes);
    elseif k == 2
        UChi = [sol.nodeValues(nodes); sol.edgeValues(faces); sol.cellMoments(K)];
    end
    
    mVals = (m(Xmon)*PNstar*UChi - u(Xhat)).^2;
                            %   Multilply by wheights and determinants.
    detAw = rldecode(detA,nq*ones(nTri,1),1).*repmat(vol*w',nTri,1); 
    mVals = bsxfun(@times,mVals,detAw);
    
    l2Err(K) = sum(mVals,1);
        
end