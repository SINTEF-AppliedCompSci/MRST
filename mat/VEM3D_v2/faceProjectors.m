function I = faceProjectors(G, k)

% m      = @(X) ...
%          [ones(size(X,1),1) , ...   %   (0,0)
%           X(:,1)            , ...   %   (1,0)
%           X(:,2)            , ...   %   (0,1)
%           X(:,1).^2         , ...   %   (2,0)
%           X(:,1).*X(:,2)    , ...   %   (1,1)
%           X(:,2).^2         , ...   %   (0,2) 
%                             ];
% 
% grad_m = @(X) ...
%          [ones(size(X,1),1)    , zeros(size(X,1),1); ...
%           zeros(size(X,1),1)   , ones(size(X,1),1) ; ...
%           X(:,1)*2             , zeros(size(X,1),1); ...
%           X(:,2)               , X(:,1)            ; ...
%           zeros(size(X,1),1)   , X(:,2)*2          ];
% 
% int_m  = @(X) ...
%          [X(:,1)                , ...
%           X(:,1).^2/2           , ...
%           X(:,1).*X(:,2)        , ...
%           X(:,1).^3/3           , ...
%           X(:,1).^2.*X(:,2)/2   , ...
%           X(:,1).*X(:,2).^2     ];

[m ,grad_m, int_m] = retrieve2DMonomials(k);
      
m3D    = @(X) ...
         [X(:,1)    , ...
          X(:,2)    , ...
          X(:,3)    ];

nk = (k+1)*(k+2)/2;
      
                            %   Face data.
nF          = G.faces.num;
Fc          = G.faces.centroids;
hF          = G.faces.diameters;
aF          = G.faces.areas;
faceNormals = G.faces.normals;

                            %   Edge data for each face.
edgeNum     = mcolon(G.faces.edgePos(1:end-1),G.faces.edgePos(2:end)-1);
edges       = G.faces.edges(edgeNum);
Ec          = G.edges.centroids(edges,:);
hE          = G.edges.lengths(edges);
edgeNormals = G.faces.edgeNormals(edgeNum,:);

                            %   Node data for each edge of each face.
                            %   Sort nodes to be in counter-clockwise
                            %   order. nodes(1:nN) = all first nodes of
                            %   each edge, nodes(nN+1:2*nN) = all last
                            %   nodes of each edge.
nodeNum = mcolon(G.edges.nodePos(edges),G.edges.nodePos(edges+1)-1);
nodes   = G.edges.nodes(nodeNum);
nodes   = reshape(nodes,2,[])';
nN      = size(nodes,1);
nodes(G.faces.edgeSign(edgeNum) == -1,:) ...
        = nodes(G.faces.edgeSign(edgeNum) == -1,2:-1:1);
nodes   = reshape(nodes,[],1);

X = [G.nodes.coords(nodes,:); Ec];


                            %   Build local coordinate systems.
                            %   x -> Tx + b
vec1 = (X(G.faces.nodePos(1:end-1)+1,:) - X(G.faces.nodePos(1:end-1),:));
vec1 = bsxfun(@rdivide, vec1, sqrt(sum(vec1.^2,2)));
vec2 = cross(faceNormals,vec1,2);
vec2 = bsxfun(@rdivide, vec2, sqrt(sum(vec2.^2,2)));
vec1 = vec1'; vec2 = vec2';
T    = [vec1(:), vec2(:)];
b    = X(G.faces.nodePos(1:end-1),:);

                            %   Map from polygon to local face coords
X = cell2mat(cellfun(@(X,Y,b) (X-repmat(b,size(X,1),1))*Y, ...
             mat2cell(X,repmat(diff(G.faces.nodePos),3,1),3), ...
             mat2cell(repmat(T,3,1),3*ones(3*nF,1),2), ...
             mat2cell(repmat(b,3,1),ones(3*nF,1),3), ...
                         'UniformOutput', false));

Fc = cell2mat(cellfun(@(X,Y,b) (X-b)*Y, ...
             mat2cell(Fc,ones(nF,1),3), ...
             mat2cell(T,3*ones(nF,1),2), ...
             mat2cell(b,ones(nF,1),3), ...
                         'UniformOutput', false));
                     
edgeNormals = ...
    cell2mat(cellfun(@(X,Y) X*Y, ...
    mat2cell(edgeNormals,diff(G.faces.edgePos),3), ...
    mat2cell(T,3*ones(nF,1),2), ...
    'UniformOutput', false));

                            %   Scale monomial coordinates.
numFaceNodes = diff(G.faces.nodePos);
Xmon = bsxfun(@rdivide, X - repmat(rldecode(Fc,numFaceNodes,1),3,1), ...
                            repmat(rldecode(hF,numFaceNodes,1),3,1));                     

                            %   Scale edgeNormals by length.
edgeNormals = bsxfun(@times, edgeNormals, hE);

%%  CALCULATE INTEGRALS FOR D MATRICES                                   %%

if k == 1
    
    D = m(Xmon(1:nN,:));
    
    intB = .5*sum(grad_m(Xmon(2*nN+1:end,:)).*repmat(edgeNormals,2,1),2);
    intB = reshape(intB, nN, 2);

              
    tmp  = intB(G.faces.edgePos(2:end) - 1,:);
    intB2 = zeros(size(intB));
    intB2(2:nN,:) = intB(1:nN-1,:);
    intB2(G.faces.edgePos(1:end-1),:) = tmp;

    intB = bsxfun( ...
                  @rdivide   , ...
                  intB + intB2       , ...
                  rldecode(hF,diff(G.faces.edgePos),1) );    
       
    NF = diff(G.faces.edgePos);
              
    BT = [ones(sum(NF), 1)./rldecode(NF,NF,1), intB];
    

elseif k == 2
    intD = bsxfun(@times, ...
                  (int_m(Xmon(1:nN,:)) + int_m(Xmon(nN+1:2*nN,:)))/6    ... 
                                       + int_m(Xmon(2*nN+1:end,:))*2/3, ...
                   edgeNormals(:,1));

    intD = cell2mat(cellfun(@(X) sum(X,1)                   , ...
                    mat2cell(intD,diff(G.faces.edgePos),6)  , ....
                    'UniformOutput', false));

    intD = bsxfun(@times, intD, hF./aF);

%%  CALCULATE INTEGRALS FOR B MARTICES                                   %%

intB = sum(grad_m(Xmon).*repmat(edgeNormals,5*3 ,1),2);
intB = reshape(intB,3*nN,5);
tmp  = intB(G.faces.nodePos(2:end) - 1 + nN,:);
intB(nN+2:2*nN,:) = intB(nN+1:2*nN-1,:);
intB(G.faces.nodePos(1:end-1)+nN,:) = tmp;

intB = bsxfun( ...
       @rdivide                                                       , ...
       [(intB(1:nN,:) + intB(nN+1:2*nN,:))/6; intB(2*nN+1:end,:)*2/3] , ...
       repmat(rldecode(hF,diff(G.faces.edgePos),1),2,1));

%%  BUILD FACE PROJECTION OPERATORS Pi^\Nabla_{F,*}                      %%
   
diffVec = cumsum(diff(G.faces.nodePos));
ii = [mcolon(G.faces.nodePos(1:end-1)                        ...
             + [0;diffVec(1:end-1) + (1:nF-1)']            , ...
             G.faces.nodePos(2:end)                          ...
             + [0;diffVec(1:end-1) + (1:nF-1)'] -1)        , ...
      mcolon(G.faces.nodePos(1:end-1) + diffVec + (0:nF-1)', ...
             G.faces.nodePos(2:end)   + diffVec + (0:nF-1)' - 1)];
NF = 2*diff(G.faces.nodePos) + 1;
N = sum(2*diff(G.faces.nodePos) + 1);

                                %   Matrices B built as transpose
BT         = zeros(N,6);
BT(ii,2:6) = intB;
vec        = zeros(nF,6);
vec(:,1)   = 1; vec(:, [4,6]) = [-2*aF./hF.^2, -2*aF./hF.^2];
BT(2*diffVec' + (1:nF),:) = vec;

D = zeros(N,6);

D(ii,:) = m(Xmon([1:nN, 2*nN+1:3*nN],:));
D(2*diffVec' + (1:nF),:) = intD;

end

BT = mat2cell(BT,NF, nk);
D  = mat2cell(D,NF,nk);

                                %   Speed: Do sparse block..
PNstar = cellfun(@(BT,D) (BT'*D)\BT', BT, D, 'UniformOutput', false);   

%%  FACE INTEGRALS                                                       %%

                            %   Gauss-Lobatto quadrature point and
                            %   wheights for refenrence triangle.

% http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html


[Xq, w, ~, vol] = triangleQuadRule(k+1);
     
nq = size(Xq,1);

N = G.nodes.num + G.edges.num*(k-1) + G.faces.num*k*(k-1)/2;
I = sparse(G.cells.num*6,N);

Kc = G.cells.centroids;
hK = G.cells.diameters;
TPos = (0:3:3*nF)+1;
PNFstarPos = (0:nk:nk*G.faces.num)+1;
intPos = (0:nk:nk*G.cells.num)+1;
cellFaces = [G.cells.faces(:,1), ...
             rldecode((1:G.cells.num)',diff(G.cells.facePos),1)];
neighbors = G.faces.neighbors;

for i = 1:nF
                        
                            %   Cells sharing face i.
    cells = cellFaces(cellFaces(:,1) == i,2);
    nK = size(cells,1);
    
                            %   Node data for face i.
    nodeNum = G.faces.nodePos(i):G.faces.nodePos(i+1)-1;
    faceNodes = nodes(nodeNum,:);
    XF = X(nodeNum,:);  

                            %   Edge data for face i.
    edgeNum = G.faces.edgePos(i):G.faces.edgePos(i+1)-1;
    faceEdges = edges(edgeNum);

                            %   Retrieve map from local to global coordinates.
    TF = T(TPos(i):TPos(i+1)-1,:);
    bT = b(i,:);

                            %   Projection matrix for face i.
%     PNFstar = PNstar(PNFstarPos(i):PNFstarPos(i+1)-1,:);
    PNFstar = PNstar{i};
                            %   Triangulate face
    tri = delaunay(XF);
    nTri = size(tri,1);

                            %   Construct map from refrence triangle to
                            %   triangles in triangulation.
    bA = XF(tri(:,1),:);
    A = XF(tri(:,2:end),:) - repmat(bA,2,1);
    A = A(mcolon(1:nTri,2*nTri,nTri),:);
    A = mat2cell(A,2*ones(nTri,1),2);
    detA = cellfun(@(X) abs(det(X)), A);
    
                            %   Map from reference triangle to triangels.
    Xhat = cell2mat(cellfun(@(X) Xq*X, A, 'uniformOutput', false));
    Xhat = Xhat + rldecode(bA,nq*ones(nTri,1),1);

    if k == 2
                            %   Map form local to global coordinates and 
                            %   scale for use in 3D monomials.
    XKmon = Xhat*TF' + repmat(bT,nTri*nq,1);
    XKmon = bsxfun(@rdivide, ...
                repmat(XKmon,nK,1) - rldecode(Kc(cells,:),nq*nTri*ones(nK,1),1), ...
                rldecode(hK(cells),nq*nTri*ones(nK,1),1));
    
    end
                            %   Scale coordinates for use in 2D monomials.
    XFmon = (Xhat - repmat(Fc(i,:),nq*nTri,1))/hF(i);

                            %   Scale and correct directions of face
                            %   normals.
    faceNormal = faceNormals(i,:)/aF(i);
    faceSign = (-ones(nK,1)).^(repmat(neighbors(i,1),nK,1) ~= cells); 
    faceNormal =  bsxfun(@times, repmat(faceNormal,nq*nTri*nK,1),...
                      rldecode(faceSign,nq*nTri*ones(nK,1),1));
                  
                            %   Evaluate monomials at quadrature points.
    mVals = m(XFmon)*PNFstar;
    
    if k == 1
        grad_mVals = faceNormal;
        grad_mVals = bsxfun(@rdivide, grad_mVals,...
                        rldecode(hK(cells),nq*nTri*ones(nK,1),1));
    elseif k == 2
    m3Dvals = m3D(XKmon);
    grad_mVals = ...
      [2*m3Dvals(:,1).*faceNormal(:,1)                                , ...
         m3Dvals(:,2).*faceNormal(:,1) + m3Dvals(:,1).*faceNormal(:,2), ...
         m3Dvals(:,3).*faceNormal(:,1) + m3Dvals(:,1).*faceNormal(:,3), ...
       2*m3Dvals(:,2).*faceNormal(:,2)                                , ...
         m3Dvals(:,3).*faceNormal(:,2) + m3Dvals(:,2).*faceNormal(:,3), ...
       2*m3Dvals(:,3).*faceNormal(:,3)];
    grad_mVals = bsxfun(@rdivide, grad_mVals,...
                        rldecode(hK(cells),nq*nTri*ones(nK,1),1));
    end       
                            %   Multilply by wheights and determinants.
    detAw = repmat(rldecode(detA,nq*ones(nTri,1),1).*repmat(vol*w',nTri,1),nK,1); 
    grad_mVals = bsxfun(@times,grad_mVals,detAw);
    grad_mVals = mat2cell(grad_mVals,nq*nTri*ones(nK,1),nk);
        
                            %   Evaluate integrals.
    int = cell2mat(cellfun(@(X) X'*mVals, grad_mVals, 'UniformOutput', false));
    
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

end

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