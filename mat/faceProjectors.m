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
X = bsxfun(@rdivide, X - repmat(rldecode(Fc,numFaceNodes,1),3,1), ...
                         repmat(rldecode(hF,numFaceNodes,1),3,1));

                            %   Apply coordinate transform to coordinates
                            %   and edgeNormals.
X = cell2mat(cellfun(@(X,Y) X*Y, mat2cell(X,repmat(diff(G.faces.nodePos),3,1),3), ...
                        mat2cell(repmat(T,3,1),3*ones(3*nF,1),2), ...
                         'UniformOutput', false));
                     
edgeNormals = ...
    cell2mat(cellfun(@(X,Y) X*Y, ...
    mat2cell(edgeNormals,diff(G.faces.edgePos),3), ...
    mat2cell(T,3*ones(nF,1),2), ...
    'UniformOutput', false));
    
                            %   Scale edgeNormals by length.
edgeNormals = bsxfun(@times, edgeNormals, hE);

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
PN     = cellfun(@(D, PNstar) D*PNstar, D, PNstar, 'UniformOutput', false); 

PNstar = cell2mat(PNstar);
PN = cell2mat(PN);

%%  Cell face int

                            %   Gauss-Lobatto quadrature point and
                            %   wheights for refenrence triangle.
Xq = [0.0, 0.0; 0.0, 1.0; 0.5, 0.0; 0.5, 0.5; 0.0, 0.5; 0.5, 0.25];
w = [1/36; 1/36; 1/18; 1/18; 1/9; 2/9];
                            %   For each triangle t, evaluate integral.

nq = size(Xq,1);

k = 2;
N = G.nodes.num + G.edges.num*(k-1) + G.faces.num*k*(k-1)/2;

I = sparse(G.cells.num*9,N);

Kc = G.cells.centroids;
hK = G.cells.diameters;
TPos = (0:3:3*nF)+1;
PNFstarPos = (0:6:6*G.faces.num)+1;
intPos = (0:9:9*G.cells.num)+1;
cellFaces = [G.cells.faces(:,1), ...
             rldecode((1:G.cells.num)',diff(G.cells.facePos),1)];


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
    
    X = (X - repmat(bT,size(X,1),1))*TF;
    
%     Xmon = (Xmon - repmat(bT,size(Xmon,1),1))*TF;
% 
%     Ec = G.edges.centroids(edges,:);
%     Ec = (Ec-repmat(bT,size(Ec,1),1))*TF;
%     
%     g = @(X) X(:,1);
%     gI = polygonInt3D(G,i,g)./aF(i);
% 
%     gv = [g([X;Ec]); gI];
%     g(X) - m(Xmon)*PNFstar*gv
    
                            %   Triangulate face
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
    XK = Xhat*TF' + repmat(bT,nTri*nq,1);
    
    XK = bsxfun(@rdivide, ...
                repmat(XK,nK,1) - rldecode(Kc(cells,:),nq*nTri*ones(nK,1),1), ...
                rldecode(hK(cells),nq*nTri*ones(nK,1),1));

    XF = (Xhat - repmat(Fc(i,:)*TF,nq*nTri,1))./hF(i);
    
    mVals = m(XF)*PNFstar;
    faceNormal = faceNormals(i,:)/aF(i);
    grad_mVals = sum(bsxfun(@times,grad_m3D(XK),faceNormal),2)./ ...
                     rldecode(hK(cells),9*nq*nTri*ones(nK,1),1);
    
    Dw = repmat(rldecode(D,nq*ones(nTri,1),1).*repmat(w,nTri,1),nK,1);
    grad_mVals = reshape(grad_mVals, nK*nq*nTri,9);
    grad_mVals = bsxfun(@times,grad_mVals,Dw);

    grad_mVals = mat2cell(grad_mVals,nq*nTri*ones(nK,1),9);
    int = cell2mat(cellfun(@(X) X'*mVals, grad_mVals, 'UniformOutput', false));
        
    dofVec = [nodes', edges' + G.nodes.num, i + G.nodes.num + G.edges.num];
    
    intNum = mcolon(intPos(cells),intPos(cells+1)-1);
    
    I(intNum, dofVec) = I(intNum, dofVec) + int;
   
end

end