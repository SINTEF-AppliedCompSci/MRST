function I = faceProjectors(G)

m =      @(X) [ones(size(X,1),1), ...
               X(:,1)               , ...   %   (1,0,0)
               X(:,2)               , ...   %   (0,1,0)
               X(:,1).^2          , ...   %   (2,0,0)
               X(:,1).*X(:,2) , ...   %   (1,1,0)
               X(:,2).^2 , ...   %   (0,2,0) 
               ];

grad_m = @(X) [...
    ones(size(X,1),1)    , zeros(size(X,1),1); ...
    zeros(size(X,1),1)   , ones(size(X,1),1) ; ...
    X(:,1)*2             , zeros(size(X,1),1); ...
    X(:,2)               , X(:,1)            ; ...
    zeros(size(X,1),1)   , X(:,2)*2          ];

int_m = @(X)        [X(:,1), ...
                     X(:,1).^2/2  , ...
                     X(:,1).*X(:,2)   , ...
                     X(:,1).^3/3         , ...
                     X(:,1).^2.*X(:,2)/2 , ...
                     X(:,1).*X(:,2).^2      , ...
                     ];
                 
grad_m3D = @(X) ...
[ones(size(X,1),1)   , zeros(size(X,1),1), zeros(size(X,1),1) ;...    %   (1,0,0)
zeros(size(X,1),1)   , ones(size(X,1),1) , zeros(size(X,1),1) ;...    %   (0,1,0)
zeros(size(X,1),1)   , zeros(size(X,1),1), ones(size(X,1),1)  ;...    %   (0,0,1)
X(:,1)*2             , zeros(size(X,1),1), zeros(size(X,1),1) ;...    %   (2,0,0)
X(:,2)               , X(:,1)            , zeros(size(X,1),1) ;...    %   (1,1,0)
X(:,3)               , zeros(size(X,1),1), X(:,1)             ;...    %   (1,0,1)
zeros(size(X,1),1)   , X(:,2)*2          , zeros(size(X,1),1) ;...    %   (0,2,0)
zeros(size(X,1),1)   , X(:,3)            , X(:,2)             ;...    %  (0,1,1)
zeros(size(X,1),1)   , zeros(size(X,1),1), X(:,3)*2];    %   (0,0,2)

nF = G.faces.num;
Fc = G.faces.centroids;
hF = G.faces.diameters;
aF = G.faces.areas;
faceNormals = G.faces.normals;

                            %   Edge data for each face.
edgeNum   = mcolon(G.faces.edgePos(1:end-1),G.faces.edgePos(2:end)-1);
edges = G.faces.edges(edgeNum);
Ec = G.edges.centroids(edges,:);
hE = G.edges.lengths(edges);
edgeNormals = G.faces.edgeNormals(edgeNum,:);

                            %   Node data for each edge of each face.
nodeNum = mcolon(G.edges.nodePos(edges),G.edges.nodePos(edges+1)-1);
nodes   = G.edges.nodes(nodeNum);
                            %   Sort nodes to be in counter-clockwise
                            %   order. nodes(1:nN) = all first nodes of
                            %   each edge, nodes(nN+1:2*nN) = all last
                            %   nodes of each edge.
nodes   = reshape(nodes,2,[])';
nN      = size(nodes,1);
nodes(G.faces.edgeSign(edgeNum) == -1,:) = ...
                             nodes(G.faces.edgeSign(edgeNum) == -1,2:-1:1);
nodes = reshape(nodes,[],1);

X = [G.nodes.coords(nodes,:); Ec];

Xu = X;
%   Build local coordinate systems.
vec1 = (X(G.faces.nodePos(1:end-1)+1,:) - X(G.faces.nodePos(1:end-1),:));
vec1 = bsxfun(@rdivide, vec1, sqrt(sum(vec1.^2,2)));
vec2 = cross(faceNormals,vec1,2);
vec2 = bsxfun(@rdivide, vec2, sqrt(sum(vec2.^2,2)));
vec1= vec1'; vec2 = vec2';

% ii = repmat([1,2,3,1,2,3],1,nF) + rldecode(0:3:3*nF-1,6,2);
% jj = rldecode(1:2*nF,3,2);
% T = sparse(ii,jj,reshape([vec1;vec2],[],1));

T = [vec1(:), vec2(:)];



                            %   Number of nodes per face.
numFaceNodes = diff(G.faces.nodePos);

                            %   Scale coordinates for use in the monomial
                            %   basis.
b = X(G.faces.nodePos(1:end-1),:);

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
% X = cell2mat(cellfun(@(X,Y,b) (X-repmat(b,size(X,1),1))*Y, ...
%              mat2cell(X,repmat(diff(G.faces.nodePos),3,1),3), ...
%              mat2cell(repmat(T,3,1),3*ones(3*nF,1),2), ...
%              mat2cell(repmat(b,3,1),ones(3*nF,1),3), ...
%                          'UniformOutput', false));
% ff = cell2mat(cellfun(@(X,Y,b) (X-b)*Y, ...
%              mat2cell(Fc,ones(nF,1),3), ...
%              mat2cell(T,3*ones(nF,1),2), ...
%              mat2cell(b,ones(nF,1),3), ...
%                          'UniformOutput', false));
% 
% X = bsxfun(@rdivide, X - repmat(rldecode(ff,numFaceNodes,1),3,1), ...
%                          repmat(rldecode(hF,numFaceNodes,1),3,1));                     
                     
                     
X = bsxfun(@rdivide, X - repmat(rldecode(Fc,numFaceNodes,1),3,1), ...
                         repmat(rldecode(hF,numFaceNodes,1),3,1));

                            %   Apply coordinate transform to coordinates
                            %   and edgeNormals.
X = cell2mat(cellfun(@(X,Y,b) (X)*Y, ...
             mat2cell(X,repmat(diff(G.faces.nodePos),3,1),3), ...
             mat2cell(repmat(T,3,1),3*ones(3*nF,1),2), ...
             mat2cell(repmat(b,3,1),ones(3*nF,1),3), ...
                         'UniformOutput', false));


edgeNormals = ...
    cell2mat(cellfun(@(X,Y) X*Y, ...
    mat2cell(edgeNormals,diff(G.faces.edgePos),3), ...
    mat2cell(T,3*ones(nF,1),2), ...
    'UniformOutput', false));

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
                            %   Scale edgeNormals by length.
edgeNormals = bsxfun(@times, edgeNormals, hE);

% for i = 1:nF
%     xx = Xplot(3*(i-1)+1:3*i,:);
%     plot(xx(:,1), xx(:,2), '*');
%     hold on
%     plot(xx(:,1) + edgeNormals(3*(i-1)+1:3*i,1), xx(:,2) + edgeNormals(3*(i-1)+1:3*i,2), 'o')
%     hold off
% end

intD = bsxfun(@times,(int_m(X(1:nN,:)) + int_m(X(nN+1:2*nN,:)))/6 ... 
                 + int_m(X(2*nN+1:end,:))*2/3, edgeNormals(:,1));
intD = cell2mat(cellfun(@(X) sum(X,1), ...
     mat2cell(intD,diff(G.faces.edgePos),6), 'UniformOutput', false));
intD = bsxfun(@times, intD, hF./aF);

intB = sum(grad_m(X).*repmat(edgeNormals,5*3 ,1),2);
intB = reshape(intB,3*nN,5);
tmp = intB(G.faces.nodePos(2:end) - 1 + nN,:);
intB(nN+2:2*nN,:) = intB(nN+1:2*nN-1,:);
intB(G.faces.nodePos(1:end-1)+nN,:) = tmp;

intB = bsxfun(@rdivide,[(intB(1:nN,:) + intB(nN+1:2*nN,:))/6; intB(2*nN+1:end,:)*2/3],...
              repmat(rldecode(hF,diff(G.faces.edgePos),1),2,1));

diffVec = cumsum(diff(G.faces.nodePos));
ii = [mcolon(G.faces.nodePos(1:end-1) + [0;diffVec(1:end-1) + (1:nF-1)']    , ...
             G.faces.nodePos(2:end)   + [0;diffVec(1:end-1) + (1:nF-1)'] -1), ...
      mcolon(G.faces.nodePos(1:end-1) + diffVec + (0:nF-1)'                  , ...
             G.faces.nodePos(2:end)   + diffVec + (0:nF-1)' - 1                    )];
NF = 2*diff(G.faces.nodePos) + 1;
N = sum(2*diff(G.faces.nodePos) + 1);

BT = zeros(N,6);
BT(ii,2:6) = intB;

vec =    zeros(nF,6);
vec(:,1) = 1; vec(:, [4,6]) = [-2*aF./hF.^2, -2*aF./hF.^2];
BT(2*diffVec' + (1:nF),:) = vec;

D = zeros(N,6);

D(ii,:) = m(X([1:nN, 2*nN+1:3*nN],:));
D(2*diffVec' + (1:nF),:) = intD;

BT = mat2cell(BT,NF, 6);
D = mat2cell(D,NF,6);

                                %   Speed: Do sparse block..

PNstar = cellfun(@(BT,D) (BT'*D)\BT', BT, D, 'UniformOutput', false);   
%PN     = cellfun(@(D, PNstar) D*PNstar, D, PNstar, 'UniformOutput', false); 

PNstar = cell2mat(PNstar);
%PN = cell2mat(PN);

%%  Cell face int
                            %   Gauss-Lobatto quadrature point and
                            %   wheights for refenrence triangle.
Xq = [0.0, 0.0; 0.0, 1.0; 0.5, 0.0; 0.5, 0.5; 0.0, 0.5; 0.5, 0.25];
w = [1/36; 1/36; 1/18; 1/18; 1/9; 2/9];
                            %   For each triangle t, evaluate integral.
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

nq = size(Xq,1);

k = 2;
N = G.nodes.num + G.edges.num*(k-1) + G.faces.num*k*(k-1)/2;

% I = sparse(G.cells.num*9,N);
I = sparse(G.cells.num*6,N);


Kc = G.cells.centroids;
hK = G.cells.diameters;
TPos = (0:3:3*nF)+1;
PNFstarPos = (0:6:6*G.faces.num)+1;
% intPos = (0:9:9*G.cells.num)+1;
intPos = (0:6:6*G.cells.num)+1;
cellFaces = [G.cells.faces(:,1), ...
             rldecode((1:G.cells.num)',diff(G.cells.facePos),1)];
neighbors = G.faces.neighbors;

m3D = @(X) [X(:,1), X(:,2), X(:,3)];


for i = 1:nF

    tt = T(3*(i-1)+1:3*i,:);
    
    nodeNum  =G.faces.nodePos(i):G.faces.nodePos(i+1)-1;
    nodes = G.faces.nodes(nodeNum);
    xx = G.nodes.coords(nodes,:);
    b = xx(1,:);
    xxu = xx;
    xxI = (xx-repmat(b,size(xx,1),1))*tt;
    xx = (xx - repmat(Fc(i,:),size(xx,1),1))/hF(i);
    xx = xx*tt;
    
    edgeNum = G.faces.edgePos(i):G.faces.edgePos(i+1)-1;
    edges  = G.faces.edges(edgeNum);
    ee = G.edges.centroids(edges,:);
    ee = (ee - repmat(Fc(i,:),size(xx,1),1))/hF(i);
    ee = ee*T(3*(i-1)+1:3*i,:);

    en = G.faces.edgeNormals(edgeNum,:);
    en = en*tt;
    
    ff = (Fc(i,:) - b)*tt;
    
    mm = @(X) [ ones(size(X,1),1), (X(:,1)-ff(1))/hF(i), ...
                                    (X(:,2)-ff(2))/hF(i), ...
                                    ((X(:,1)-ff(1))/hF(i)).^2, ...
                                    (X(:,1)-ff(1))/hF(i).*(X(:,2)-ff(2))/hF(i), ...
                                    ((X(:,2)-ff(2))/hF(i)).^2]; 
    %%  Check for D. OK
    dd = D{i};
    dd-[m([xx;ee]); polygonInt(xxI, mm)/aF(i)];
    
    %%  Check for B.
    
    bb = BT{i}';
    bb;
    
end

for i = 1:nF
    
    cells = cellFaces(cellFaces(:,1) == i,2);
    nK = size(cells,1);
    
    nodeNum = G.faces.nodePos(i):G.faces.nodePos(i+1)-1;
    nodes = G.faces.nodes(nodeNum);
    nN = size(nodes,1);
    X = G.nodes.coords(nodes,:);
    
    edgeNum = G.faces.edgePos(i):G.faces.edgePos(i+1)-1;
    edges = G.faces.edges(edgeNum);
    
    TF = T(TPos(i):TPos(i+1)-1,:);
    bT = X(1,:);
    
    PNFstar = PNstar(PNFstarPos(i):PNFstarPos(i+1)-1,:);

                            %   Map from Polygon to face
%     Xmon = (X-repmat(Fc(i,:),size(X,1),1))/hF(i);
    
    Xu = X;
    X = (X - repmat(bT,size(X,1),1))*TF;
% 
%     Xmon = (Xmon)*TF;
% 
%     Ec = G.edges.centroids(edges,:);
%     Ec = (Ec-repmat(bT,size(Ec,1),1))*TF;
% 
%     g = @(X) X(:,1).^2 + X(:,2) + X(:,2).^2*10000;
%     gI = polygonInt(X,g)./aF(i);
% 
%     gv = [g([X;Ec]); gI];
%     g(X) - m(Xmon)*PNFstar*gv
% 
%     TF

                            %   Triangulate face
    tri = delaunay(X);
    nTri = size(tri,1);
    
    bA = X(tri(:,1),:);
    A = X(tri(:,2:end),:) - repmat(bA,2,1);
    A = A(mcolon(1:nTri,2*nTri,nTri),:);
    A = mat2cell(A,2*ones(nTri,1),2);
    D = cellfun(@(X) abs(det(X)), A);
    
    Xhat = cell2mat(cellfun(@(X) Xq*X, A, 'uniformOutput', false));
    Xhat = Xhat + rldecode(bA,nq*ones(nTri,1),1);

    XK = Xhat*TF' + repmat(bT,nTri*nq,1);
    
    XKu = XK;
    
    XK = bsxfun(@rdivide, ...
                repmat(XK,nK,1) - rldecode(Kc(cells,:),nq*nTri*ones(nK,1),1), ...
                rldecode(hK(cells),nq*nTri*ones(nK,1),1));
    
    FcT = (Fc(i,:)-bT)*TF;
            
    XF = (Xhat - repmat(FcT,nq*nTri,1))/hF(i);
    
% %     xp = (X*TF' + repmat(bT,size(X,1),1)-repmat(Kc(cells,:),size(X,1),1))./hK(i);
%     xp = (Xu-repmat(Kc(cells,:),size(X,1),1))/hK(cells);
%     for j = 1:nTri
% %         xf = XF(nq*(i-1)+1:nq*i,:);
%         xk = XK(nq*(j-1)+1:nq*j,:);
%         plot3(xk(:,1), xk(:,2),xk(:,3),'*')
%         hold on
%         plot3(xp(3*(j-1)+1:3*j,1),xp(3*(j-1)+1:3*j,2),xp(3*(j-1)+1:3*j,3),'o')
%         hold off
%     end
% %     
%     mVals = m(XF)*PNFstar;
%     
    faceNormal = faceNormals(i,:)/aF(i);
    faceSign = (-ones(nK,1)).^(repmat(neighbors(i,1),nK,1) ~= cells);
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
    faceNormal =  bsxfun(@times, repmat(faceNormal,nq*nTri*nK,1),...
                      rldecode(faceSign,nq*nTri*ones(nK,1),1));
    mVals = m(XF)*PNFstar;
    m3Dvals = m3D(XK);
    grad_mVals = ...
        [2*m3Dvals(:,1).*faceNormal(:,1), ...
           m3Dvals(:,2).*faceNormal(:,1) + m3Dvals(:,1).*faceNormal(:,2), ...
           m3Dvals(:,3).*faceNormal(:,1) + m3Dvals(:,1).*faceNormal(:,3), ...
         2*m3Dvals(:,2).*faceNormal(:,2), ...
           m3Dvals(:,3).*faceNormal(:,2) + m3Dvals(:,2).*faceNormal(:,3), ...
         2*m3Dvals(:,3).*faceNormal(:,3)];
    grad_mVals = bsxfun(@rdivide, grad_mVals,...
                        rldecode(hK(cells),nq*nTri*ones(nK,1),1));
    Dw = repmat(rldecode(D,nq*ones(nTri,1),1).*repmat(w,nTri,1),nK,1);   
    grad_mVals = bsxfun(@times,grad_mVals,Dw);
    
    grad_mVals = mat2cell(grad_mVals,nq*nTri*ones(nK,1),6);
    int = cell2mat(cellfun(@(X) X'*mVals, grad_mVals, 'UniformOutput', false));
    
    dofVec = [nodes', edges' + G.nodes.num, i + G.nodes.num + G.edges.num];
    
    intNum = mcolon(intPos(cells),intPos(cells+1)-1);
    
    I(intNum, dofVec) = I(intNum, dofVec) + int;
    
    %   Speed: first rows of B can be obatined from F.
    
    end

end