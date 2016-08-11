function G = computeVEMIP(G,rock,k)

faces = G.cells.faces(mcolon(G.cells.facePos(1:end-1),G.cells.facePos(2:end)-1));
faceAreas = G.faces.areas(faces);
faceNormals = G.faces.normals(faces,:);
nF = numel(faces);
numCellFaces = diff(G.cells.facePos);
if size(faces,1) == 1; faces = faces'; end

nodes = G.cells.nodes(mcolon(G.cells.nodePos(1:end-1),G.cells.nodePos(2:end)-1));
numCellNodes = diff(G.cells.nodePos);
if size(nodes,1) == 1; nodes = nodes'; end

[m,~,~] = retrieveMonomials(2,1);
Xmon = bsxfun(@rdivide, G.nodes.coords(nodes,:) ...
                      - rldecode(G.cells.centroids, numCellNodes, 1), ...
                        rldecode(G.cells.diameters, numCellNodes,1));
D = m(Xmon);
D = sparseBlockDiag(D, numCellNodes, 1);

faceSign    = (-ones(nF,1)).^(G.faces.neighbors(faces,1) ~= rldecode((1:G.cells.num)', diff(G.cells.facePos), 1)); 
faceNormals = bsxfun(@times, faceNormals, faceSign);
faceNormals = bsxfun(@rdivide, faceNormals, rldecode(G.cells.diameters, numCellFaces, 1));
faceNormals = sparseBlockDiag(faceNormals, numCellFaces, 1);

K = zeros(2*G.cells.num, 2);
KK = permTensor(rock,2);
K([1:2:end,2:2:end],:) = [KK(:,1:2);KK(:,3:4)];
K = sparseBlockDiag(K, 2*ones(1,G.cells.num), 1);

B = .5*faceNormals*K;

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
    PNstar = M\B;
PN = D*PNstar;

M(1:nk:end,:) = 0;
I = speye(size(PN,1));
Ablk = PNstar'*M*PNstar + 2/3*(I-PN)'*(I-PN);

ind1 = [1, cumsum(numCellNodes.^2)'+1];
ind2 = [1, cumsum(numCellNodes)'+1];
ii = zeros(sum(numCellNodes.^2),1);
jj= zeros(sum(numCellNodes.^2),1);
A  = zeros(sum(numCellNodes.^2),1);
for P = 1:G.cells.num
    nodes = G.cells.nodes(G.cells.nodePos(P):G.cells.nodePos(P+1)-1);
    ii(ind1(P):ind1(P+1)-1) = repmat(nodes,numel(nodes),1);
    jjP = repmat(nodes,1,numel(nodes))';
    jj(ind1(P):ind1(P+1)-1) = jjP(:);
    AP = Ablk(ind2(P):ind2(P+1)-1,ind2(P):ind2(P+1)-1);
    A(ind1(P):ind1(P+1)-1) = AP(:);
end

A = sparse(ii,jj,A, G.nodes.num, G.nodes.num);
U = A\ones(size(A,1),1);

end