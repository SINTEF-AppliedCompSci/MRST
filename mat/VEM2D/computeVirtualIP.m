function S = computeVirtualIP(G,rock,fluid,k)

    [mu, ~] = fluid.properties();

faces = G.cells.faces(:,1);
faceAreas = G.faces.areas(faces);
faceNormals = G.faces.normals(faces,:);
nF = numel(faces);
numCellFaces = diff(G.cells.facePos);
faceSign    = (-ones(nF,1)).^(G.faces.neighbors(faces,1) ~= rldecode((1:G.cells.num)', diff(G.cells.facePos), 1)); 
if size(faces,1) == 1; faces = faces'; end


nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces),G.faces.nodePos(faces+1)-1));
if size(nodes,1) == 1; nodes = nodes'; end
nodes   = reshape(nodes,2,[])';
nodes(faceSign == -1,:) = nodes(faceSign == -1,2:-1:1);
nodes   = nodes(:,1);
numCellNodes = diff(G.cells.nodePos);

[m,~,~] = retrieveMonomials(2,1);
Xmon = bsxfun(@rdivide, G.nodes.coords(nodes,:) ...
                      - rldecode(G.cells.centroids, numCellNodes, 1), ...
                        rldecode(G.cells.diameters, numCellNodes,1));
D = m(Xmon);
D = sparseBlockDiag(D, numCellNodes, 1);

faceNormals = bsxfun(@times, faceNormals, faceSign);
faceNormals = bsxfun(@rdivide, faceNormals, rldecode(G.cells.diameters, numCellFaces, 1));
faceNormals = sparseBlockDiag(faceNormals, numCellFaces, 1);

K = zeros(2*G.cells.num, 2);
KK = permTensor(rock,2);
K([1:2:end,2:2:end],:) = [KK(:,1:2);KK(:,3:4)];
K = sparseBlockDiag(K, 2*ones(1,G.cells.num), 1);

B = .5*faceNormals*K/mu;

ii = 1:nF; jj = ii;
jj(2:end) = jj(1:end-1);
jj([1;cumsum(G.cells.facePos(2:end-1)    ...
            -G.cells.facePos(1:end-2))+1]) ...
 = ii(cumsum(G.cells.facePos(2:end)-G.cells.facePos(1:end-1)));

nk = (k+1)*(k+2)/2;

B = B(ii,:) + B(jj,:);
B = [.5*(faceAreas(ii) + faceAreas(jj)),squeezeBlockDiag(B,numCellFaces, nF, nk-1)];
B = sparseBlockDiag(B', numCellFaces, 2);

M = B*D;
PiNstar = M\B;
PiN = D*PiNstar;

NP = numCellNodes;
nker = NP - nk;
sigma = spdiags(rldecode(KK(:,1)+KK(:,4),nker,1), 0, sum(nker), sum(nker));

M(1:nk:end,:) = 0;

ind0   = [1; cumsum(NP.^2)+1];
ind1   = [1; cumsum(NP)+1];
ind2   = [1; cumsum(nker)+1];
ind3   = [1; cumsum(NP.*nker)+1];
ii     = zeros(sum(nker.*NP),1); jj = ii; Q = ii;
% [iiPiN, jjPiN, PiNvec] = find(PiN);
% lo = 1; hi = 0; isPiNP = 0;
% lo = 0;
for P = 1:G.cells.num 
    
%     lo = lo + sum(isPiNP(lo:hi)); hi = lo + NP(P)^2 - 1;
%     lo = lo + numel(isPiNP); hi = lo + NP(P)^2 - 1;
%     
%     isPiNP = find([false(lo-1,1); iiPiN(lo:hi) >= ind1(P) & iiPiN(lo:hi) <= ind1(P+1)-1; false(numel(iiPiN)-lo-hi+1,1)]);
%     
%     iii = iiPiN(isPiNP)-ind1(P)+1;
%     jjj = jjPiN(isPiNP)-ind1(P)+1;
%     PPP = PiNvec(isPiNP);
%     PiNP = sparse(iii, jjj, PPP);

    PiNP = PiN(ind1(P):ind1(P+1)-1,ind1(P):ind1(P+1)-1);

    QP   = orth(eye(numCellNodes(P))-PiNP);
    ii(ind3(P):ind3(P+1)-1) = repmat((ind1(P):ind1(P+1)-1)', nker(P),1);
    jjP = repmat(ind2(P):ind2(P+1)-1, NP(P), 1);
    jj(ind3(P):ind3(P+1)-1) = jjP(:);
    Q(ind3(P):ind3(P+1)-1) = QP(:);
    
end
Q = sparse(ii, jj, Q(:), sum(NP), sum(nker));
    
I = speye(size(PiN,1));
A = PiNstar'*M*PiNstar + (I-PiN)'*Q*sigma*Q'*(I-PiN);

S.A = A;
S.dofVec = nodes;
S.PNstar = PiNstar;

end