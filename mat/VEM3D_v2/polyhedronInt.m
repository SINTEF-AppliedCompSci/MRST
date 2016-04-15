function I = polyhedronInt(G,cells,f, k)

%   http://people.sc.fsu.edu/~jburkardt/m_src/tetrahedron_arbq_rule/tetrahedron_arbq_rule.html

[Xq, w, V, vol] = polyhedronQuadRule(k);
    
Vdiff = V(1:end-1,:) - V(2:end,:);
   
nK = numel(cells);
nq = size(Xq,1);
      
I = zeros(nK,size(f([0,0,0]),2));

for i = 1:nK
    
    nodeNum = G.cells.nodePos(cells(i)):G.cells.nodePos(cells(i)+1)-1;
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
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
    
    I(i,:) = vol*repmat(w,1,nTri).*rldecode(D,nq*ones(nTri,1),1)'*f(Xhat);

end

end