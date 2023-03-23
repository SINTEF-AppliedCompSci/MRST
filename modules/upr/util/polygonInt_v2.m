function I = polygonInt_v2(G,cells,f,k)
% This code is written by Oeystein Klemetsdal:
% Øystein Strengehagen Klemetsdal,
% The virtual element method as a common framework for finite element
% and finite difference methods - Numerical and theoretical analysis,
% June 2016,
% Masters thesis,
% Department of Mathematics,
% Norwegian University of Science and Technology.
% git@bitbucket.org:strene/the-virtual-element-method.git
[Xq, w, V, vol] = triangleQuadRule(k);

nq = size(Xq,1);

nK = numel(cells);

I = zeros(nK,size(f([0,0],1),2));

for i = 1:nK
    
    nodeNum = G.cells.nodePos(cells(i)):G.cells.nodePos(cells(i)+1)-1;
    nodes = G.cells.nodes(nodeNum);
    
    X = G.nodes.coords(nodes,:);
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
    I(i,:) = vol*repmat(w,1,nTri).*(rldecode(D,nq*ones(nTri,1),1))'*f(Xhat,i);
    
end

end