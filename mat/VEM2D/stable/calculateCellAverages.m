function cellAverages = calculateCellAverages(G, nodeValues, PNstarT)

nK = G.cells.num;
[m, ~, ~] = retrieveMonomials(1);

[Xq, w, ~, vol] = triangleQuadRule(2);
nq = size(Xq,1);

Kc = G.cells.centroids;
hK = G.cells.diameters;
aK = G.cells.volumes;

PNstarPos = [1, cumsum( diff(G.cells.nodePos')) + 1];

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

    PNstar = PNstarT(PNstarPos(K):PNstarPos(K+1)-1,:)';
    uChi = nodeValues(nodes);
    mVals = m(Xmon)*PNstar*uChi; 
                            %   Multilply by wheights and determinants.
    detAw = rldecode(detA,nq*ones(nTri,1),1).*repmat(vol*w',nTri,1); 
    mVals = bsxfun(@times,mVals,detAw);
    
    cellAverages(K) = sum(mVals,1)/aK(K);
        
end
    