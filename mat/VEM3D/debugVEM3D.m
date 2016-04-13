%%  VEM3D_loc

% Mdb = zeros(nk,nk);
% Mdb(1,:) = cellIntegrals/vol;
% %    1 2 3 4 5   6  7  8   9  10
% %   {1,x,y,z,x^2,xy,xz,y^2,yz,z^2}
% mm = zeros(nk,nk);
% mm(2,2) = vol; mm(3,3) = vol; mm(4,4) = vol;
% mm(2,[5,6,7])  = cellIntegrals([2,3,4]).*[2,1,1];
% mm(3,[6,8,9])  = cellIntegrals([2,3,4]).*[1,2,1];
% mm(4,[7,9,10]) = cellIntegrals([2,3,4]).*[1,1,2];
% mm(5,[5,6,7])  = cellIntegrals([5,6,7]).*[4,2,2]; 
% mm(6,6)        = cellIntegrals(5) + cellIntegrals(8);
% mm(6,[7,8,9])  = cellIntegrals([9,6,7]).*[1,2,1];
% mm(7,7)        = cellIntegrals(5) + cellIntegrals(10);
% mm(7,[9,10])   = cellIntegrals([6,7]).*[1,2];
% mm(8,[8,9])    = cellIntegrals([8,9]).*[4,2];
% mm(9,9)        = cellIntegrals(8) + cellIntegrals(10);
% mm(9,10)       = cellIntegrals(9)*2;
% mm(10,10)      = cellIntegrals(10)*4;
% 
% 
% mm = mm/hK^2;
% 
% Mdb(2:nk,2:nk) = triu(mm(2:nk,2:nk)) + tril(mm(2:nk,2:nk)',-1);
% 
% norm(Mdb-M,'fro')
% (abs(Mdb-M)>10e-6)
% 
% 
% Fc = G.faces.centroids(faces,:);
% edgeNum = mcolon(G.faces.edgePos(faces),G.faces.edgePos(faces+1)-1);
% edges = G.faces.edges(edgeNum);
% edgeNormals = G.faces.edgeNormals(edgeNum,:);
% clf;
% plotGrid(G,K)
% hold on
% plot3(X(:,1), X(:,2), X(:,3),'*')
% plot3(Ec(:,1), Ec(:,2), Ec(:,3), 'o')
% Ec = G.edges.centroids(edges,:);
% % plot3(Ec(:,1) + edgeNormals(:,1), ...
% %       Ec(:,2) + edgeNormals(:,2), ...
% %       Ec(:,3) + edgeNormals(:,3), '+')
% 
% plot3(Fc(:,1),Fc(:,2), Fc(:,3), 'o')
% plot3(Fc(:,1) + faceNormals(:,1), Fc(:,2) ...
%               + faceNormals(:,2), Fc(:,3) ...
%               + faceNormals(:,3), '+')
% 
% hold off
% axis equal

% % Q = orth(eye(NK)-PN);
% q1 = [ 1 -1 -1  1  1 -1 -1  1];
% q2 = [ 1 -1  1 -1 -1  1 -1  1];
% q3 = [ 1  1 -1 -1 -1 -1  1  1];
% q4 = [-1  1  1 -1  1 -1 -1  1];
% 
% hx = abs(max(X(:,1)) - min(X(:,1)))/2;
% hy = abs(max(X(:,2)) - min(X(:,2)))/2;
% hz = abs(max(X(:,3)) - min(X(:,3)))/2;

% delta1 = sqrt(9/(8*hx*hy*hz));
% delta2 = sqrt(27/(8*hx*hy*hz));
% Q = [delta1*q1', delta1*q2', delta1*q3', delta2*q4'];
% P = Q'*Q;
% S = diag(alpha, 0);
     

% g = @(X) X(:,1).^2 + X(:,3).*X(:,2)*1000/3 + 10;
% g = @(X) ones(size(X,1),1);
% 
% gF = polygonInt3D(G,faces,g);
% gK = polyhedronInt(G,K,g);
% gv = [g([X;Ec]); gF./faceAreas; gK/vol];
% 
% er1 = max(abs((gv - PN*gv)./gv))
% er2 = max(abs(g(X) - m3D(X)*PNstar*gv)./g(X))


% D = [monomialNodeVals; monomialEdgeVals;  ...
%      monomialFaceInt ; monomialCellInt/vol]; 
%  
%             m3D =      @(X) [ones(size(X,1),1) , ...
%                 (X(:,1)-Kc(1))/hK              , ...   %   (1,0,0)
%                (X(:,2)-Kc(2))/hK               , ...   %   (0,1,0)
%                (X(:,3)-Kc(3))/hK               , ...   %   (0,0,1)
%                (X(:,1)-Kc(1)).^2/hK^2          , ...   %   (2,0,0)
%                (X(:,1)-Kc(1)).*(X(:,2)-Kc(2))/hK^2 , ...   %   (1,1,0)
%                (X(:,1)-Kc(1)).*(X(:,3)-Kc(3))/hK^2, ...   %   (1,0,1)
%                (X(:,2)-Kc(2)).^2/hK^2, ...   %   (0,2,0) 
%                (X(:,2)-Kc(2)).*(X(:,3)-Kc(3))/hK^2, ...   %   (0,1,1)
%                (X(:,3)-Kc(3)).^2/hK^2];      %   (0,0,2)icenter of K.
%            
% Xm = (X-repmat(Kc,size(X,1),1))./hK;
% Ecm = (Ec-repmat(Kc,size(Ec,1),1))./hK;
% mFI= bsxfun(@rdivide, polygonInt3D(G,faces,m3D),faceAreas);
% mCI = polyhedronInt(G,K,m3D)/vol;
% D-[m3D([X;Ec]);mFI;mCI]


% m3D =      @(X) [ones(size(X,1),1) , ...
%                 X(:,1)              , ...   %   (1,0,0)
%                X(:,2)               , ...   %   (0,1,0)
%                X(:,3)               , ...   %   (0,0,1)
%                X(:,1).^2          , ...   %   (2,0,0)
%                X(:,1).*X(:,2) , ...   %   (1,1,0)
%                X(:,1).*X(:,3), ...   %   (1,0,1)
%                X(:,2).^2, ...   %   (0,2,0) 
%                X(:,2).*X(:,3), ...   %   (0,1,1)
%                X(:,3).^2];      %   (0,0,2)icenter of K.
% 
% 
% Mdb = zeros(nk,nk);
% Mdb(1,:) = cellIntegrals/vol;
% %    1 2 3 4 5   6  7  8   9  10
% %   {1,x,y,z,x^2,xy,xz,y^2,yz,z^2}
% mm = zeros(nk,nk);
% mm(2,2) = vol/(hK^2); mm(3,3) = vol/(hK^2); mm(4,4) = vol/(hK^2);
% mm(2,5) = cellIntegrals(2)*2/(hK^2);
% mm(2,[6,7]) = cellIntegrals([3,4])/hK^2;
% mm(3,[6,8,9]) = cellIntegrals([2,3,4]).*[1/hK^2,2/hK^2, 1/hK^2];
% mm(4,[7,9,10]) = cellIntegrals([2,3,4]).*[1/hK^2,2/hK^2, 1/hK^2];
% mm(5,[5,6,7]) = cellIntegrals([5,6,7]).*[4/hK^2,2/hK^2, 2/hK^2]; 
% mm(6,6) = (cellIntegrals(5) + cellIntegrals(8))/hK^2;
% mm(6,[7,8,9]) = cellIntegrals([9,6,7])/hK^2;
% mm(7,7) = (cellIntegrals(5) + cellIntegrals(10))/hK^2;
% mm(7,[9,10]) = cellIntegrals([6,7]).*[1/hK^2,2/hK^2];
% mm(8,[8,9]) = cellIntegrals([8,9]).*[4/hK^2,2/hK^2];
% mm(9,9) = (cellIntegrals(8) + cellIntegrals(10))/hK^2;
% mm(9,10) = cellIntegrals(9)*2/hK^2;
% mm(10,10) = cellIntegrals(10)*4/hK^2;
% 
% Mdb(2:nk,2:nk) = triu(mm(2:nk,2:nk)) + tril(mm(2:nk,2:nk)',-1);
% 
% norm(Mdb-M,'fro')
% Mdb-M

%%  DELETED STUFF

% faceIntNum = G.cells.faceIntPos(K):G.cells.faceIntPos(K+1)-1;
% monomialFaceInt     = bsxfun(@rdivide, ...
%                      G.cells.monomialFaceIntegrals(faceIntNum,:), ...
%                      faceAreas);
% mIint = G.cells.monomialFaceIntegrals(faceIntNum,:)


% monomialCellInt = G.cells.monomialCellIntegrals(K,:);