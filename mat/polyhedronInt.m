function I = polyhedronInt(G, f)

nK = G.cells.num;
% w        = [1/216, 1/216, ...
%             1/54, 1/54, 1/54, 1/54, 1/54, ...
%             1/36, ...
%             2/27, 2/27, 2/27, 2/27, ...
%             1/6, ...
%             1/9, ...
%             8/27]/6;
% Xq       = [0,0,0;   0,0,1; ...
%             0,0,0.5; 0,0.5,0; 0,0.5,0.5; 0.5,0,0; 0.5,0,0.5; ...
%             0,1,0;
%             0,0.5,0.25; 0.5,0,0.25; 0.5,0.25,0; 0.5,0.25,0.25; ...
%             1,0,0;
%             0.5,0.5,0;
%             0.5,0.25,0.125];

w = [1/216, 1/216, 1/216, 1/216, ...
     1/108, 1/108, 1/108, 1/108, ...
     1/54, 1/54, ...
     1/27, 1/27];
Xq = [0.0,0.0,0.0; 0.0,0.0,1.0; 0.5,0.0,0.0; 0.5,0.0,0.5; ...
      0.0,0.5,0.0; 0.0,0.5,0.5; 0.5,0.25,0.0; 0.5,0.25,0.25; ...
      0.0,0.0,0.5; 0.5,0.0,0.25; ...
      0.0,0.5,0.25; 0.5,0.25,0.125];

nq = size(Xq,1);
      
I = zeros(nK,size(f([0,0,0]),2));

nK = 1;

X = [0,0,0; 1,0,0; 0,1,0; 0,0,1];
% 1,0,0; 0,1,0; 0,0,1];
% X = [0,0,1; 0,0,0; 0,1,0; 1,0,0];

for i = 1:nK
    
    nodeNum = G.cells.nodePos(i):G.cells.nodePos(i+1)-1;
    nodes = G.cells.nodes(nodeNum);
%     X = G.nodes.coords(nodes,:);
    tri = delaunay(X);
    nTri = size(tri,1);
    tri = [1,2,3,4];
%     IT = 0;
%     for k= 1:nTri
%     
%         b = X(tri(k,1),:)
%         A = X(tri(k,2:end),:)-repmat(b,3,1)
%         D = abs(det(A));
%         XhatT = Xq*A + repmat(b,nq,1);
%         IT = IT + D*w*f(XhatT);
%     end
    
    b = X(tri(:,1),:);
    
    A = X(tri(:,2:end),:) - repmat(b,3,1);
    A = A(mcolon(1:nTri,3*nTri,nTri),:);
    A = mat2cell(A,3*ones(nTri,1),3);
    
    D = cellfun(@(X) abs(det(X)), A);
    
    
    
    Xhat = cell2mat(cellfun(@(X) Xq*X, A, 'uniformOutput', false));
    Xhat = Xhat + rldecode(b,nq*ones(nTri,1),1);

    D*w*f(Xhat)
%     I(i,:) = repmat(w,1,nTri)*(bsxfun(@times, rldecode(D,nq*ones(nTri,1),1),f(Xhat)));
    I(i,:) = repmat(w,1,nTri).*rldecode(D,nq*ones(nTri,1),1)'*f(Xhat);
    
    figure();
    view(3);
    tetramesh(tri,X)
    for j = 1:nTri
        hold on;

        Xh = Xhat((j-1)*12+1:j*12,:);
        Xh = Xhat;
        plot3(Xh(:,1), Xh(:,2), Xh(:,3),'or')
        pause;
        
        
    end
end

end