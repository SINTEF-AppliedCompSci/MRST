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

%%  FACEPROJECTORS

%%  DELETED STUFF

                 
% grad_m3D = @(X) ...
% [ones(size(X,1),1)   , zeros(size(X,1),1), zeros(size(X,1),1) ;...    %   (1,0,0)
% zeros(size(X,1),1)   , ones(size(X,1),1) , zeros(size(X,1),1) ;...    %   (0,1,0)
% zeros(size(X,1),1)   , zeros(size(X,1),1), ones(size(X,1),1)  ;...    %   (0,0,1)
% X(:,1)*2             , zeros(size(X,1),1), zeros(size(X,1),1) ;...    %   (2,0,0)
% X(:,2)               , X(:,1)            , zeros(size(X,1),1) ;...    %   (1,1,0)
% X(:,3)               , zeros(size(X,1),1), X(:,1)             ;...    %   (1,0,1)
% zeros(size(X,1),1)   , X(:,2)*2          , zeros(size(X,1),1) ;...    %   (0,2,0)
% zeros(size(X,1),1)   , X(:,3)            , X(:,2)             ;...    %  (0,1,1)
% zeros(size(X,1),1)   , zeros(size(X,1),1), X(:,3)*2];    %   (0,0,2)



% ii = repmat([1,2,3,1,2,3],1,nF) + rldecode(0:3:3*nF-1,6,2);
% jj = rldecode(1:2*nF,3,2);
% T = sparse(ii,jj,reshape([vec1;vec2],[],1));

% X = bsxfun(@rdivide, X - repmat(rldecode(Fc,numFaceNodes,1),3,1), ...
%                          repmat(rldecode(hF,numFaceNodes,1),3,1));
% 
%                             %   Apply coordinate transform to coordinates
%                             %   and edgeNormals.
% X = cell2mat(cellfun(@(X,Y,b) (X)*Y, ...
%              mat2cell(X,repmat(diff(G.faces.nodePos),3,1),3), ...
%              mat2cell(repmat(T,3,1),3*ones(3*nF,1),2), ...
%              mat2cell(repmat(b,3,1),ones(3*nF,1),3), ...
%                          'UniformOutput', false));

%PN     = cellfun(@(D, PNstar) D*PNstar, D, PNstar, 'UniformOutput', false); 
%PN = cell2mat(PN);

%   w = [ ...
%     0.22338158967801146570; ...
%     0.22338158967801146570; ...
%     0.22338158967801146570; ...
%     0.10995174365532186764; ...
%     0.10995174365532186764; ...
%     0.10995174365532186764 ];
% 
%   Xq = [ ...
%     0.10810301816807022736,  0.44594849091596488632; ...
%     0.44594849091596488632,  0.10810301816807022736; ...
%     0.44594849091596488632,  0.44594849091596488632; ...
%     0.81684757298045851308,  0.091576213509770743460; ...
%     0.091576213509770743460,  0.81684757298045851308; ...
%     0.091576213509770743460,  0.091576213509770743460 ];

%     faceNormal = bsxfun(@times, repmat(faceNormal,nq*nTri*nK*9,1),...
%                  repmat(rldecode(faceSign,nq*nTri*ones(nK,1),1),9,1));
%     grad_mVals = sum(grad_m3D(XK).*faceNormal,2)./...
%                  repmat(rldecode(hK(cells),nq*nTri*ones(nK,1),1),9,1);
%                         
%     Dw = repmat(rldecode(D,nq*ones(nTri,1),1).*repmat(w,nTri,1),nK,1);
%     grad_mVals = reshape(grad_mVals, nq*nTri*nK,9);
%     grad_mVals = bsxfun(@times,grad_mVals,Dw);
% 
%     grad_mVals = mat2cell(grad_mVals,nq*nTri*ones(nK,1),9);
%     int = cell2mat(cellfun(@(X) X'*mVals, grad_mVals, 'UniformOutput', false));
%  

% Xq = [0.0, 0.0; 0.0, 1.0; 0.5, 0.0; 0.5, 0.5; 0.0, 0.5; 0.5, 0.25];
% w  = [1/36    ; 1/36    ; 1/18    ; 1/18    ; 1/9     ; 2/9      ];


%%  DEBUG

% for i = 1:nF
%     tt = T(3*(i-1)+1:3*i,:)'
%     xx = X(4*(i-1)+1:4*i,:);
%     plot3(xx(:,1), xx(:,2), xx(:,3), '*')
%     hold on
%     plot3(xx(1,1) + tt(1,1), ...
%           xx(1,2) + tt(1,2), ...
%           xx(1,3) + tt(1,3), 'o')
%     pause;
%     plot3(xx(1,1) + tt(2,1), ...
%           xx(1,2) + tt(2,2), ...
%           xx(1,3) + tt(2,3), 'o')
%     plot3(b(i,1), b(i,2), b(i,3), 'sq');
%     pause;
%     hold off
% end

% Xplot = cell2mat(cellfun(@(X,Y,b) (X-repmat(b,size(X,1),1))*Y, ...
%              mat2cell(X,repmat(diff(G.faces.nodePos),3,1),3), ...
%              mat2cell(repmat(T,3,1),3*ones(3*nF,1),2), ...
%              mat2cell(repmat(b,3,1),ones(3*nF,1),3), ...
%                          'UniformOutput', false));
% 
%



% sqrt(sum((Xplot(1:nN,:) - Xplot(nN+1:2*nN,:)).^2,2)) - hE
% 
% for i = 1:nF
%     figure();
%     plot(Xplot(3*(i-1)+1:3*i,1), Xplot(3*(i-1)+1:3*i,2),'*')
%     figure();
%     plot3(Xu(3*(i-1)+1:3*i,1),Xu(3*(i-1)+1:3*i,2), Xu(3*(i-1)+1:3*i,3),'o')
%     pause
%     close all;
% end

% for i = 1:nF
%     xx = Xplot(3*(i-1)+1:3*i,:);
%     plot(xx(:,1), xx(:,2), '*');
%     hold on
%     plot(xx(:,1) + edgeNormals(3*(i-1)+1:3*i,1), xx(:,2) + edgeNormals(3*(i-1)+1:3*i,2), 'o')
%     hold off
% end

% for i = 1:nF
% 
%     tt = T(3*(i-1)+1:3*i,:);
%     
%     nodeNum  =G.faces.nodePos(i):G.faces.nodePos(i+1)-1;
%     nodes = G.faces.nodes(nodeNum);
%     xx = G.nodes.coords(nodes,:);
%     b = xx(1,:);
%     xxu = xx;
%     xxI = (xx-repmat(b,size(xx,1),1))*tt;
%     xx = (xx - repmat(Fc(i,:),size(xx,1),1))/hF(i);
%     xx = xx*tt;
%     
%     edgeNum = G.faces.edgePos(i):G.faces.edgePos(i+1)-1;
%     edges  = G.faces.edges(edgeNum);
%     ee = G.edges.centroids(edges,:);
%     ee = (ee - repmat(Fc(i,:),size(xx,1),1))/hF(i);
%     ee = ee*T(3*(i-1)+1:3*i,:);
% 
%     en = G.faces.edgeNormals(edgeNum,:);
%     en = en*tt;
%     
%     ff = (Fc(i,:) - b)*tt;
%     
%     mm = @(X) [ ones(size(X,1),1), (X(:,1)-ff(1))/hF(i), ...
%                                     (X(:,2)-ff(2))/hF(i), ...
%                                     ((X(:,1)-ff(1))/hF(i)).^2, ...
%                                     (X(:,1)-ff(1))/hF(i).*(X(:,2)-ff(2))/hF(i), ...
%                                     ((X(:,2)-ff(2))/hF(i)).^2]; 
%     %%  Check for D. OK
%     dd = D{i};
%     dd-[m([xx;ee]); polygonInt(xxI, mm)/aF(i)];
%     
%     %%  Check for B.
%     
%     bb = BT{i}';
%     bb;
%     
% end


                            %   Map from Polygon to face
%     Xmon = (X-repmat(Fc(i,:),size(X,1),1))/hF(i);
    
%      XuF = Xu(nodeNum,:);
%     Xmon = (XuF-repmat(Fc(i,:),size(XuF,1),1))/hF(i);
% %     X = (X - repmat(bT,size(X,1),1))*TF;
% 
%     Xmon = (Xmon)*TF;
% % 
%     Ec = G.edges.centroids(faceEdges,:);
%     Ec = (Ec-repmat(bT,size(Ec,1),1))*TF;
% 
%     g = @(X) X(:,1).^2 + X(:,2) + X(:,2).^2;
%     gI = polygonInt(XF,g)./aF(i);
% 
%     gv = [g([XF;Ec]); gI];
%     g(XF) - m(Xmon)*PNFstar*gv
% % 
%     TF

    
% %     xp = (X*TF' + repmat(bT,size(X,1),1)-repmat(Kc(cells,:),size(X,1),1))./hK(i);
%     xp = (XuF-repmat(Kc(cells(1),:),size(XuF,1),1))/hK(cells(1));
%     for j = 1:nTri
% %         xf = XF(nq*(i-1)+1:nq*i,:);
%         xk = XK(nq*(j-1)+1:nq*j,:);
%         plot3(xk(:,1), xk(:,2),xk(:,3),'*')
%         hold on
%         plot3(xp(3*(j-1)+1:3*j,1),xp(3*(j-1)+1:3*j,2),xp(3*(j-1)+1:3*j,3),'o')
%         hold off
%     end

% 
%     xp = (XF-repmat(FcT,size(XF,1),1))/hF(i);
%     for j = 1:nTri
%         xf = XFmon(nq*(j-1)+1:nq*j,:);
% %         xk = XK(nq*(j-1)+1:nq*j,:);
%         plot(xf(:,1), xf(:,2),'*')
%         hold on
%         plot(xp(3*(j-1)+1:3*j,1),xp(3*(j-1)+1:3*j,2),'o')
%         hold off
%         pause
%     end

%     for j = 1:nK
%         c = cells(j);
%         nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1;
%         nn = G.cells.nodes(nodeNum);
%         xx = G.nodes.coords(nn,:);
%         xx = (xx-repmat(Kc(c,:), size(xx,1),1))/hK(c);
%         plot3(xx(:,1), xx(:,2), xx(:,3), 'o')
%         hold on;
%         ii = nq*(j-1)+1:nq*j;
%         plot3(XKmon(ii,1), XKmon(ii,2), XKmon(ii,3), '*')
%         hold off
%         pause
%     end

        
%         xx = (XF - repmat(Fc(i,:),size(XF,1),1))/hF(i);
%         plot(xx(:,1), xx(:,2), '*')
%         hold on
%         plot(XFmon(:,1), XFmon(:,2), 'o');
%         hold off
%         pause