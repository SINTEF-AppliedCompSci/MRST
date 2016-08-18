function S = computeVirtualIP(G,rock,fluid,k)


%%  EXTRACT FLUID PROPERTIES

[mu, ~] = fluid.properties();

%%  GET NODE AND FACE DATA

%   Faces for each cell
faces = G.cells.faces(:,1);
faceAreas = G.faces.areas(faces);
faceNormals = G.faces.normals(faces,:);
nF = numel(faces);
numCellFaces = diff(G.cells.facePos);

%   Face signs
faceSign    = (-ones(nF,1)).^(G.faces.neighbors(faces,1) ...
                  ~= rldecode((1:G.cells.num)', diff(G.cells.facePos), 1)); 
              
faceNormals = bsxfun(@times, faceNormals, faceSign);

if size(faces,1) == 1; faces = faces'; end

nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), ...
                             G.faces.nodePos(faces+1)-1));
if size(nodes,1) == 1; nodes = nodes'; end
nodes   = reshape(nodes,2,[])';
nodes(faceSign == -1,:) = nodes(faceSign == -1,2:-1:1);
nodes   = nodes(:,1);
nN = numel(nodes);
numCellNodes = diff(G.cells.nodePos);

%%  FUNCTION SPACE DIMENSIONS

nP = G.cells.num;

nk   = (k+1)*(k+2)/2;
nkk  = k*(k+1)/2;
NP   = numCellNodes + numCellFaces*(k-1) + k*(k-1)/2;
nker = NP - nk;

%%  CALCULATE D MATRICES

[m, grad_m, int_m] = retrieveMonomials(2,k);

if k == 1
    
    Xmon = bsxfun(@rdivide, ...
                  G.nodes.coords(nodes,:) ...
                  - rldecode(G.cells.centroids, numCellNodes, 1), ...
                    rldecode(G.cells.diameters, numCellNodes,1));
    D = m(Xmon);
    D1 = D(:, 1:nkk);
    D = sparseBlockDiag(D, numCellNodes, 1);
    D1 = sparseBlockDiag(D1, NP, 1);

else
    
    Xmon = bsxfun(@rdivide, ...
                  [G.nodes.coords(nodes,:); G.faces.centroids(faces,:)] ...
                   - repmat(rldecode(G.cells.centroids, numCellNodes, 1),2,1), ...
                     repmat(rldecode(G.cells.diameters, numCellNodes,1),2,1));
    ii = 1:nN; jj = ii;
    jj(1:end-1) = jj(2:end);
    jj(cumsum(numCellNodes)) = ii([1;cumsum(numCellNodes(1:end-1))+1]);
    intD = bsxfun(@times, ( int_m(Xmon(ii,:)) + int_m(Xmon(jj,:)) )/6 ...
                          + int_m(Xmon(nN+1:end,:))*2/3, ...
                          faceNormals(:,1)                   );
    I = sparseBlockDiag(ones(1, sum(numCellFaces)), numCellFaces, 2);
    intD = bsxfun(@times, I*intD, G.cells.diameters./G.cells.volumes);
    
    
    nodeDof = mcolon([1;cumsum(NP(1:end-1))+1],[1;cumsum(NP(1:end-1))+1]+numCellNodes-1);
    edgeDof = nodeDof + rldecode(numCellNodes, numCellNodes, 1)';
    D = zeros(sum(NP), nk);
    D([nodeDof, edgeDof],:) = m(Xmon);
    D(cumsum(NP),:) = intD;
    
    D = sparseBlockDiag(D, NP, 1);
    
end

%%  CALCULATE B MATRICES

%   Make permeability matrices, divided by mu.

K = zeros(2*G.cells.num, 2);
KK = permTensor(rock,2);
K([1:2:end,2:2:end],:) = [KK(:,1:2);KK(:,3:4)]/mu;

if k == 1
    
    K = sparseBlockDiag(K, 2*ones(1,G.cells.num), 1);
    
    faceNormals = bsxfun(@rdivide, faceNormals, ...
                                 rldecode(G.cells.diameters, numCellFaces, 1));
    faceNormals = sparseBlockDiag(faceNormals, numCellFaces, 1);

    B = .5*faceNormals*K/mu;

    ii = 1:nF; jj = ii;
    jj(2:end) = jj(1:end-1);
    jj([1;cumsum(G.cells.facePos(2:end-1)    ...
                -G.cells.facePos(1:end-2))+1]) ...
     = ii(cumsum(G.cells.facePos(2:end)-G.cells.facePos(1:end-1)));

    nk = (k+1)*(k+2)/2;

    B = B(ii,:) + B(jj,:);
    B = [.5*(faceAreas(ii) + faceAreas(jj)), ...
                                   squeezeBlockDiag(B,numCellFaces, nF, nk-1)];
    
    B1 = B(:,1:nkk);
    B1 = sparseBlockDiag(B1', NP, 2);
                               
    B = sparseBlockDiag(B', NP, 2);

else

    K = sparseBlockDiag(repmat(K,5,1), 2*ones(1,5*G.cells.num), 1);
    
    intB = sparseBlockDiag(grad_m(Xmon), repmat(numCellNodes+numCellFaces,5,1), 1)*K;
    intB = squeezeBlockDiag(intB, repmat(numCellNodes + numCellFaces,5,1), sum(numCellNodes + numCellFaces)*5, 2);    
    
    %   Dot product by length-weighted face normals.
    
    ii = 1:nN; jj = ii;
    jj(2:end) = jj(1:end-1);
    jj(G.cells.nodePos(1:end-1)) = ii(G.cells.nodePos(2:end)-1);
    
    iiN = repmat(1:nN, 1, 5) + rldecode(0:(nN+nF):4*(nN+nF), nN*ones(1,5), 2);
    iiF =  iiN + nN;
    
    intB = sum(intB([iiN, iiN, iiF],:).*...
               [repmat(faceNormals,5,1); repmat(faceNormals(jj,:),5,1); repmat(faceNormals,5,1)], 2);
    
    %   Evaluate line integrals using three-point Gauss-Lobatto.
           
    intB = [reshape((intB(1:numel(iiN)) + intB(numel(iiN)+1:2*numel(iiN)))/6, nN, 5);
            reshape(intB(2*numel(iiN)+1:end)*2/3, nN, 5)];
    
    intB = bsxfun(@rdivide, intB, repmat(rldecode(G.cells.diameters,diff(G.cells.nodePos),1),2,1));

    %   Assmble matrices.
    
    B = zeros(sum(NP), nk);
    B([nodeDof, edgeDof],2:nk) = intB;
    vec = zeros(G.cells.num,6);
    vec(:, [1,4,6]) = [ones(G.cells.num,1), ...
                       bsxfun(@times, -2*[KK(:,1), KK(:,1)], ...
                       G.cells.volumes./(mu*G.cells.diameters.^2))];
   B(cumsum(NP),:) = vec;
   
   B = sparseBlockDiag(B', NP, 2);

end

    %%  CALCULATE PROJECTOION OPERATORS

M = B*D;

[ii, jj] = blockDiagIndex(repmat(nk, [G.cells.num ,1]));
kk = sub2ind(size(M), ii, jj);
PiNstar = sparse(ii, jj, invv(full(M(kk)), repmat(nk, [G.cells.num, 1])))*B;
PiN = D*PiNstar;

clear B D;

sigma = spdiags(rldecode(KK(:,1)+KK(:,4),nker,1), 0, sum(nker), sum(nker));

M(1:nk:end,:) = 0;

%%  CALCULATE H MATRIX

if k == 1
    PiN1 = (B1*D1)\B1;
end

%%  CALCULATE Q MATRIX

ind2   = [1; cumsum(nker)+1];

ind1   = [1; cumsum(NP)+1];

ind3   = [1; cumsum(NP.*nker)+1];
ii     = zeros(sum(nker.*NP),1); jj = ii; Q = ii;

PiNPos = [1; cumsum(NP.^2) + 1];
QPos   = [1; cumsum(NP.*nker)+1];

ii = blockDiagIndex(NP);
PiNvec = full(PiN(ii));

for P = 1:G.cells.num 
    
    PiNP = reshape(PiNvec(PiNPos(P):PiNPos(P+1)-1), NP(P), NP(P));
    QP = orth(eye(NP(P))-PiNP);
    Q(QPos(P):QPos(P+1)-1) = QP(:);
    
%     PiNP = PiN(ind1(P):ind1(P+1)-1,ind1(P):ind1(P+1)-1);
% 
%     QP   = orth(eye(NP(P))-PiNP);
%     ii(ind3(P):ind3(P+1)-1) = repmat((ind1(P):ind1(P+1)-1)', nker(P),1);
%     jjP = repmat(ind2(P):ind2(P+1)-1, NP(P), 1);
%     jj(ind3(P):ind3(P+1)-1) = jjP(:);
%     Q(ind3(P):ind3(P+1)-1) = QP(:);
    
end

[ii,jj] = blockDiagIndex(NP, nker);
Q = sparse(ii, jj, Q, sum(NP), sum(nker));

%%  CALCULATE LOCAL BLOCK MATRICES
    
I = speye(size(PiN,1));
A = PiNstar'*M*PiNstar + (I-PiN)'*Q*sigma*Q'*(I-PiN);

%%  

%%  MAKE SOLUTION STRUCT

if G.griddim == 2
    vec = [1, cumsum(NP(1:end-1))' + 1];
    iiN = mcolon(vec, vec + numCellNodes'-1);
    iiF = mcolon(vec + numCellNodes', vec + numCellNodes' + numCellFaces'*(k-1)-1);
    iiP = mcolon(vec + numCellNodes' + numCellFaces', vec + numCellNodes' + numCellFaces' + k*(k-1)/2 -1);
    if k == 1
        dofVec([iiN, iiF, iiP]) = nodes';
    else
        dofVec([iiN, iiF, iiP]) = [nodes', faces' + G.nodes.num, (1:G.cells.num) + G.nodes.num + G.faces.num*(k-1)];
    end
end

S.A = A;
S.dofVec = dofVec;
S.PNstar = PiNstar;
if k == 1
    S.PiN1 = PiN1;
end
S.order  = k;

end



% ii = zeros(sum(NP.*nker),1); jj = ii; Q = ii;
% QQ = speye(sum(NP))-PiN;
% [iiQQ, jjQQ, QQ] = find(QQ);
% c = ismember(jjQQ,ind2(1:end));
% 
% ii(1:nnz(c)) = iiQQ(c); jj(1:nnz(c)) = jjQQ(c); Qvec(1:nnz(c)) = QQ(c);
% Q = sparse(ii, jj, Qvec);
% 
% while sprank(Q) < sum(nker);
% end
% 
% 

    
%     lo = lo + sum(isPiNP(lo:hi)); hi = lo + NP(P)^2 - 1;
%     lo = lo + numel(isPiNP); hi = lo + NP(P)^2 - 1;
%     
%     isPiNP = find([false(lo-1,1); iiPiN(lo:hi) >= ind1(P) & iiPiN(lo:hi) <= ind1(P+1)-1; false(numel(iiPiN)-lo-hi+1,1)]);
%     
%     iii = iiPiN(isPiNP)-ind1(P)+1;
%     jjj = jjPiN(isPiNP)-ind1(P)+1;
%     PPP = PiNvec(isPiNP);
%     PiNP = sparse(iii, jjj, PPP);

% [iiPiN, jjPiN, PiNvec] = find(PiN);
% lo = 1; hi = 0; isPiNP = 0;
% lo = 0;