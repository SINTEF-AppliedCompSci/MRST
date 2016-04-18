function G = VEM3D_faceProjectors(G, k)
%--------------------------------------------------------------------------
%   Calculates face projection opreators and face integrals for the kth
%   order virtual element method.
%
%   SYNOPSIS:
%      G = VEM3D_faceProjectors(G, k)
%
%   DESCRIPTION:
%       Calculates matrix representations for the projection opreators
%       \Pi^{\nabla, F}_k for each face of G, and evaluates the boundary
%       integrals (\partial_n m, \phi)_{0,\partial K} for each cell K. See
%       [1] for details.
%
%   REQUIRED PARAMETERS:
%       G   - 2D MRST grid, with sorted edges and computed VEM geometry,
%             G = computeVEMGeometry(G).
%       k   - Order of the method. Supproted values are 1 and 2.
%
%   RETURNS:
%       G   - Grid updated with fields G.cells.Bint, integral inveolved in
%             B matrix calculations, and map G.cells.BintPos,
%             G.faces.PNstarT, face projectors, saved as trasnspose, and
%             map G.faces.PNstarPos.
%
%   REFERENCES:
%       [1]     - Thesis title.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}


%%  RETRIEVE MONOMIALS, CALCULATE FUNCTION SPACE DIMENSION               %%

[m ,grad_m, int_m] = retrieveMonomials(2,k);
nk = (k+1)*(k+2)/2;

%%  FACE, NODE AND CELL DATA                                             %%

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

%   Node data for each edge of each face. Sort nodes to be in
%   counter-clockwise order. nodes(1:nN) = all first nodes of each edge,
%   nodes(nN+1:2*nN) = all last nodes of each edge.

nodeNum = mcolon(G.edges.nodePos(edges),G.edges.nodePos(edges+1)-1);
nodes   = G.edges.nodes(nodeNum);
nodes   = reshape(nodes,2,[])';
nN      = size(nodes,1);
nodes(G.faces.edgeSign(edgeNum) == -1,:) ...
        = nodes(G.faces.edgeSign(edgeNum) == -1,2:-1:1);
nodes   = reshape(nodes,[],1);

X = [G.nodes.coords(nodes,:); Ec];

%%  MAP FROM GLOBAL TO LOCAL COORDINATES                                 %%

%   Build local coordinate systems. x -> Tx + b.

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

%%  CALCULATE INTEGRALS FOR D AND B MATRICES                             %%

if k == 1
    
    detA = m(Xmon(1:nN,:));
    
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

    intB = sum(grad_m(Xmon).*repmat(edgeNormals,5*3 ,1),2);
    intB = reshape(intB,3*nN,5);
    tmp  = intB(G.faces.nodePos(2:end) - 1 + nN,:);
    intB(nN+2:2*nN,:) = intB(nN+1:2*nN-1,:);
    intB(G.faces.nodePos(1:end-1)+nN,:) = tmp;

    intB = bsxfun( ...
           @rdivide                                                       , ...
           [(intB(1:nN,:) + intB(nN+1:2*nN,:))/6; intB(2*nN+1:end,:)*2/3] , ...
           repmat(rldecode(hF,diff(G.faces.edgePos),1),2,1));

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

    detA = zeros(N,6);

    detA(ii,:) = m(Xmon([1:nN, 2*nN+1:3*nN],:));
    detA(2*diffVec' + (1:nF),:) = intD;

end

BT = mat2cell(BT,NF, nk);
detA  = mat2cell(detA,NF,nk);

PNstarT = cell2mat(cellfun(@(BT,D) BT/(BT'*D)', BT, detA, 'UniformOutput', false));
PNstarPos = [0, cumsum(diff(G.faces.nodePos') + diff(G.faces.edgePos')*(k-1) + k*(k-1)/2)] + 1;


%%  FACE INTEGRALS                                                       %%

% http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html

[Xq, w, ~, vol] = triangleQuadRule(2*k-1);
     
nq = size(Xq,1);

N = G.nodes.num + G.edges.num*(k-1) + G.faces.num*k*(k-1)/2;
I = zeros(G.cells.num*6,N);

Kc = G.cells.centroids;
hK = G.cells.diameters;
TPos = (0:3:3*nF)+1;
intPos = (0:nk:nk*G.cells.num)+1;
cellFaces = [G.cells.faces(:,1), ...
             rldecode((1:G.cells.num)',diff(G.cells.facePos),1)];
neighbors = G.faces.neighbors;

if k == 2
    [~, grad_m3D, ~] = retrieveMonomials(3,k, 'face', true);
end

for F = 1:nF
                        
    %   Cell, edge and node data for face F
    
    cells = cellFaces(cellFaces(:,1) == F,2);
    nK = size(cells,1);
    nodeNum = G.faces.nodePos(F):G.faces.nodePos(F+1)-1;
    faceNodes = nodes(nodeNum,:);
    XF = X(nodeNum,:);
    
    edgeNum = G.faces.edgePos(F):G.faces.edgePos(F+1)-1;
    faceEdges = edges(edgeNum);

    %   Retrieve map from local to global coordinates.
    
    TF = T(TPos(F):TPos(F+1)-1,:);
    bT = b(F,:);

    PNFstar = PNstarT(PNstarPos(F):PNstarPos(F+1)-1,:)';

    %   Triangulate face. Construct map from refrence triangle to triangles
    %   in triangulation.
    
    tri = delaunay(XF);
    nTri = size(tri,1);

    bA = XF(tri(:,1),:);
    A = XF(tri(:,2:end),:) - repmat(bA,2,1);
    A = A(mcolon(1:nTri,2*nTri,nTri),:);
    
    Ad = reshape(A',2,2,[]);
    detA = abs(squeeze(Ad(1,1,:).*Ad(2,2,:) - Ad(1,2,:).*Ad(2,1,:)));
       
    ii = repmat((1:2*nTri)',2,1);
    jj = 1:2*nTri;
    jj = reshape(jj,2,[])';
    jj = repmat(jj(:)',2,1);
    jj = jj(:);
    
    A = sparse(ii, jj, A(:), 2*nTri, 2*nTri);
    
    XhatTmp = repmat(Xq,1,nTri)*A + repmat(reshape(bA',1,[]),nq,1);
    Xhat = zeros(nq*nTri,2);
    Xhat(:,1) = reshape(XhatTmp(:,1:2:end),[],1);
    Xhat(:,2) = reshape(XhatTmp(:,2:2:end),[],1);

    if k == 2
    
    %   Map form local to global coordinates and scale for use in 3D
    %   monomials.
    
        XKmon = Xhat*TF' + repmat(bT,nTri*nq,1);
        XKmon = bsxfun(@rdivide, ...
                    repmat(XKmon,nK,1) - rldecode(Kc(cells,:),nq*nTri*ones(nK,1),1), ...
                    rldecode(hK(cells),nq*nTri*ones(nK,1),1));
    
    end
    
    %   Scale coordinates for use in 2D monomials.
    
    XFmon = (Xhat - repmat(Fc(F,:),nq*nTri,1))/hF(F);

    %   Scale and correct directions of face normals.
    
    faceNormal = faceNormals(F,:)/aF(F);
    faceSign = (-ones(nK,1)).^(repmat(neighbors(F,1),nK,1) ~= cells); 
    faceNormal =  bsxfun(@times, repmat(faceNormal,nq*nTri*nK,1),...
                      rldecode(faceSign,nq*nTri*ones(nK,1),1));
                  
    %   Evaluate monomials at quadrature points.
    
    mVals = m(XFmon)*PNFstar;
    
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
    
    %   Evaluate integrals.
    
    detAw = repmat(rldecode(detA,nq*ones(nTri,1),1).*repmat(vol*w',nTri,1),nK,1); 
    
    grad_mVals = bsxfun(@times,grad_mVals,detAw);
    
    ii = repmat((1:nq*nTri*nK)',nk,1);
    jj = 1:nk*nK; jj = reshape(jj,nk,[])';
    jj = repmat(jj(:)',nq*nTri,1);
    jj = jj(:);
    grad_mVals2 = sparse(ii, jj, grad_mVals(:), nq*nTri*nK, nk*nK);
    
    int = grad_mVals2'*repmat(mVals, nK,1);
    
    %   Construct local to global map.
    
    if k == 1
        dofVec = faceNodes';
    elseif k == 2
        dofVec = [faceNodes', faceEdges' + G.nodes.num, F + G.nodes.num + G.edges.num];
    end
    
    intNum = mcolon(intPos(cells),intPos(cells+1)-1);
    
    I(intNum, dofVec) = I(intNum, dofVec) + int;
        
end

BintPos = (0:nk:nk*G.cells.num) + 1;
G.cells.('Bint')      = I;
G.cells.('BintPos')   = BintPos;
G.faces.('PNstarT')   = PNstarT;
G.faces.('PNstarPos') = PNstarPos;


end

