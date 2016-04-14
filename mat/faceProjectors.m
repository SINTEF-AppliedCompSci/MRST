function G = VEM3D_faceProjectors(G, k)

%%  RETRIEVE MONOMIALS, CALCULATE DIMENIOSN OF POLYNOMIAL SPACE          %%

[m2D ,grad_m2D, int_m2D] = retrieveMonomials(2,k);
nk = (k+1)*(k+2)/2;

%%  FACE DATA                                                            %%

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
%   Sort nodes to be in counter-clockwise order. nodes(1:nN) = all first
%   nodes of each edge, nodes(nN+1:2*nN) = all last nodes of each edge.

nodeNum = mcolon(G.edges.nodePos(edges),G.edges.nodePos(edges+1)-1);
nodes   = G.edges.nodes(nodeNum);
nodes   = reshape(nodes,2,[])';
nN      = size(nodes,1);
nodes(G.faces.edgeSign(edgeNum) == -1,:) ...
        = nodes(G.faces.edgeSign(edgeNum) == -1,2:-1:1);
nodes   = reshape(nodes,[],1);
X = [G.nodes.coords(nodes,:); Ec];

%%  MAP COORINATES AND NORMALS FROM 3D TO LOCAL 2D COORDINATES           %%

%   Build local coordinate systems and map X -> XT + b from 2D to 3D
%   coordinates.

vec1 = (X(G.faces.nodePos(1:end-1)+1,:) - X(G.faces.nodePos(1:end-1),:));
vec1 = bsxfun(@rdivide, vec1, sqrt(sum(vec1.^2,2)));
vec2 = cross(faceNormals,vec1,2);
vec2 = bsxfun(@rdivide, vec2, sqrt(sum(vec2.^2,2)));
vec1 = vec1'; vec2 = vec2';
T    = [vec1(:), vec2(:)];
b    = X(G.faces.nodePos(1:end-1),:);

%   Map from polygon to local face coordinates.

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

%   Scale edgeNormals by length.
edgeNormals = bsxfun(@times, edgeNormals, hE);

numFaceNodes = diff(G.faces.nodePos);
Xmon = bsxfun(@rdivide, X - repmat(rldecode(Fc,numFaceNodes,1),3,1), ...
                            repmat(rldecode(hF,numFaceNodes,1),3,1));                     

%%  CALCULATE INTEGRALS FOR D AND B MATRICES                             %%

hFVec = rldecode(hF,diff(G.faces.edgePos),1);

if k == 1
    
    D = m2D(Xmon(1:nN,:));
    
    intB = .5*sum(grad_m2D(Xmon(2*nN+1:end,:)).*repmat(edgeNormals,2,1),2);
    intB = reshape(intB, nN, 2);

              
    tmp  = intB(G.faces.edgePos(2:end) - 1,:);
    intB2 = zeros(size(intB));
    intB2(2:nN,:) = intB(1:nN-1,:);
    intB2(G.faces.edgePos(1:end-1),:) = tmp;

    intB = bsxfun( ...
                  @rdivide   , ...
                  intB + intB2       , ...
                  hFVec );    
       
    NF = diff(G.faces.edgePos);
              
    BT = [ones(sum(NF), 1)./rldecode(NF,NF,1), intB];
    
elseif k == 2
    
    intD = bsxfun(@times, ...
                  (int_m2D(Xmon(1:nN,:)) + int_m2D(Xmon(nN+1:2*nN,:)))/6    ... 
                                       + int_m2D(Xmon(2*nN+1:end,:))*2/3, ...
                   edgeNormals(:,1));

    intD = cell2mat(cellfun(@(X) sum(X,1)                   , ...
                    mat2cell(intD,diff(G.faces.edgePos),6)  , ....
                    'UniformOutput', false));

    intD = bsxfun(@times, intD, hF./aF);

    intB = sum(grad_m2D(Xmon).*repmat(edgeNormals,5*3 ,1),2);
    intB = reshape(intB,3*nN,5);
    tmp  = intB(G.faces.nodePos(2:end) - 1 + nN,:);
    intB(nN+2:2*nN,:) = intB(nN+1:2*nN-1,:);
    intB(G.faces.nodePos(1:end-1)+nN,:) = tmp;

    intB = bsxfun( ...
           @rdivide                                                       , ...
           [(intB(1:nN,:) + intB(nN+1:2*nN,:))/6; intB(2*nN+1:end,:)*2/3] , ...
           repmat(hFVec,2,1));

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

    D(ii,:) = m2D(Xmon([1:nN, 2*nN+1:3*nN],:));
    D(2*diffVec' + (1:nF),:) = intD;

end

BT = mat2cell(BT,NF, nk);
D  = mat2cell(D,NF,nk);

%%  BUILD FACE PROJECTION OPERATORS Pi^\Nabla_{F,*}                      %%

                                %   Speed: Do sparse block..
% PNstar = cellfun(@(BT,D) (BT'*D)\BT', BT, D, 'UniformOutput', false);   
PNFstarT = cell2mat(cellfun(@(BT,D) BT/(BT'*D)', BT, D, 'UniformOutput', false));


%%  FACE INTEGRALS                                                       %%

                            %   Gauss-Lobatto quadrature point and
                            %   wheights for refenrence triangle.

% http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html


[Xq, w, ~, vol] = triangleQuadRule(k+1);
     
if k == 2
    [~, grad_m3D, ~] = retrieveMonomials(3,2,'face', true);
end

nq = size(Xq,1);

N = G.nodes.num + G.edges.num*(k-1) + G.faces.num*k*(k-1)/2;

Kc = G.cells.centroids;
hK = G.cells.diameters;
TPos = (0:3:3*nF)+1;
PNFstarPos = [1, cumsum(diff(G.faces.nodePos') + diff(G.faces.edgePos')*(k-1) + k*(k-1)/2) + 1];
intPos = (0:nk:nk*G.cells.num)+1;
cellFaces = [G.cells.faces(:,1), ...
             rldecode((1:G.cells.num)',diff(G.cells.facePos),1)];
neighbors = G.faces.neighbors;

ii = [];
jj = [];
I  = [];

for F = 1:nF
                        
                            %   Cells sharing face i.
    cells = cellFaces(cellFaces(:,1) == F,2);
    nK = size(cells,1);
    
                            %   Node data for face i.
    nodeNum = G.faces.nodePos(F):G.faces.nodePos(F+1)-1;
    faceNodes = nodes(nodeNum,:);
    XF = X(nodeNum,:);  

                            %   Edge data for face i.
    edgeNum = G.faces.edgePos(F):G.faces.edgePos(F+1)-1;
    faceEdges = edges(edgeNum);

                            %   Retrieve map from local to global coordinates.
    TF = T(TPos(F):TPos(F+1)-1,:);
    bT = b(F,:);
                            %   Projection matrix for face i.
    PNFstar = PNFstarT(PNFstarPos(F):PNFstarPos(F+1)-1,:)';
                            %   Triangulate face
    tri = delaunay(XF);
    nTri = size(tri,1);

                            %   Construct map from refrence triangle to
                            %   triangles in triangulation.
    bA = XF(tri(:,1),:);
    A = XF(tri(:,2:end),:) - repmat(bA,2,1);
    A = A(mcolon(1:nTri,2*nTri,nTri),:);
    At = reshape(A,2,2,[]);
    detA = squeeze(At(1,1,:).*At(2,2,:) - At(2,1,:).*At(1,2,:));
    A = mat2cell(A,2*ones(nTri,1),2);
%     detA = cellfun(@(X) abs(det(X)), A);
    
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
    XFmon = (Xhat - repmat(Fc(F,:),nq*nTri,1))/hF(F);

                            %   Scale and correct directions of face
                            %   normals.
    faceNormal = faceNormals(F,:)/aF(F);
    faceSign = (-ones(nK,1)).^(repmat(neighbors(F,1),nK,1) ~= cells); 
    faceNormal =  bsxfun(@times, repmat(faceNormal,nq*nTri*nK,1),...
                      rldecode(faceSign,nq*nTri*ones(nK,1),1));
                  
                            %   Evaluate monomials at quadrature points.
    mVals = m2D(XFmon)*PNFstar;
    
    if k == 1
        grad_mVals = faceNormal;
        grad_mVals = bsxfun(@rdivide, grad_mVals,...
                        rldecode(hK(cells),nq*nTri*ones(nK,1),1));
    elseif k == 2
        grad_mVals = sum(grad_m3D(XKmon).*repmat(faceNormal,6,1),2);
        grad_mVals = reshape(grad_mVals, nq*nTri*nK, 6);
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
        dofVec = [faceNodes', faceEdges' + G.nodes.num, F + G.nodes.num + G.edges.num];
    end
    
    intNum = mcolon(intPos(cells),intPos(cells+1)-1);    
    
    iiF = repmat(intNum', numel(dofVec), 1);
    jjF = repmat(repmat(dofVec , numel(intNum(1:nk)), 1),1,nK);
    jjF = jjF(:);
    
    ii = [ii; iiF];
    jj = [jj; jjF];
    
    I = [I;int(:)];
    
%     I(intNum, dofVec) = I(intNum, dofVec) + int;
    
    %   Speed: first rows of B can be obatined from F.
    
end

I = sparse(ii, jj, I, intPos(end)-1, N); 

nk = (k+1)*(k+2)/2;
BintPos = (0:nk:nk*G.cells.num) + 1;
G.cells.('Bint') = I;
G.cells.('BintPos') = BintPos;
G.faces.('PNstarT') = PNFstarT;
G.faces.('PNstarPos') = PNFstarPos;

end