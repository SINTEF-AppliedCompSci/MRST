function sol = cellAverages(G, sol, PNstarT)

nK = G.cells.num;
nk = (k+1)*(k+2)/2;

nNK = diff(G.cells.nodePos);
cellAverages = zeros(nK, 1);

[Xq, w, ~, vol] = triangleQuadRule(1);
     
nq = size(Xq,1);

Kc = G.cells.centroids;
hK = G.cells.diameters;

PNstarPos = [1, cumsum( diff(G.cells.nodePos')) + 1];

[m, ~, ~] = retrieveMonomials(1);

for K = 1:nK
                            %   Node data for cell K.
    nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
    nodes   = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);

    tri = delaunay(XF);
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

    PNstar = PNstar(PNstarPos(K):PNst
    
    mVals = m(Xmon)*PNstar;
    
                            %   Multilply by wheights and determinants.
    detAw = repmat(rldecode(detA,nq*ones(nTri,1),1).*repmat(vol*w',nTri,1),nK,1); 
    mVals = bsxfun(@times,mVals,detAw);
    mVals = mat2cell(mVals,nq*nTri*ones(nK,1),nk);
        
    uChi = sol.nodeValues(nodes);
                            %   Evaluate integrals.
    int = cell2mat(cellfun(@(X) X'*mVals, uChi, 'UniformOutput', false));
    
                            %   Construct local to global map.
    if k == 1
        dofVec = faceNodes';
    elseif k == 2
        dofVec = [faceNodes', faceEdges' + G.nodes.num, i + G.nodes.num + G.edges.num];
    end
    
    intNum = mcolon(intPos(cells),intPos(cells+1)-1);
    
    I(intNum, dofVec) = I(intNum, dofVec) + int;
    
    %   Speed: first rows of B can be obatined from F.
    


        
end
    