function [AK, bK] = VEM2D_projectors(G,f, k, alpha)


[m, grad_m, int_m] = retrieveMonomials(k);


nK = G.cells.num;
Kc = G.cells.centroids;
hK = G.cells.diameters;
aK = G.cells.volumes;

edgeNum = mcolon(G.cells.facePos(1:end-1),G.cells.facePos(2:end)-1);
edges = G.cells.faces(edgeNum);
if size(edges,1) == 1;
    edges = edges';
end
nE          = size(edges,1);
Ec          = G.faces.centroids(edges,:);
hE = G.faces.areas(edges,:);
edgeNormals = G.faces.normals(edges,:);
cellEdges = rldecode((1:nK)', diff(G.cells.facePos), 1);
edgeSign = (-ones(nE,1)).^(G.faces.neighbors(edges,1) ~= cellEdges); 
edgeNormals = bsxfun(@times, ...
                          edgeNormals, ...
                          edgeSign);

nodeNum = mcolon(G.faces.nodePos(edges),G.faces.nodePos(edges+1)-1);
nodes = G.faces.nodes(nodeNum);
if size(nodes,1) == 1
    nodes = nodes';
end
nodes = reshape(nodes,2,[])';
nN = size(nodes,1);
nodes(edgeSign == -1,:) ...
        = nodes(edgeSign == -1,2:-1:1);
% nodes   = reshape(nodes,[],1);
nodes   = nodes(:);
X = [G.nodes.coords(nodes,:); Ec];

                            %   Scale monomial coordinates.
numCellNodes = diff(G.cells.nodePos);
Xmon = bsxfun(@rdivide, X - repmat(rldecode(Kc,numCellNodes,1),3,1), ...
                            repmat(rldecode(hK,numCellNodes,1),3,1));                     


%%  CALCULATE B AND D MATRICES                                           %%
if k == 1

    D = m(Xmon(1:nN,:));
    
    intB = .5*sum(grad_m(Xmon(2*nN+1:end,:)).*repmat(edgeNormals,2 ,1),2);
    intB = reshape(intB, nN, 2);

              
    tmp  = intB(G.cells.facePos(2:end) - 1,:);
    intB2 = zeros(size(intB));
    intB2(2:nN,:) = intB(1:nN-1,:);
    intB2(G.cells.facePos(1:end-1),:) = tmp;

    intB = bsxfun( ...
                  @rdivide   , ...
                  intB + intB2       , ...
                  rldecode(hK,diff(G.cells.facePos),1) );    
       
    NK = diff(G.cells.nodePos);
              
    BT = [ones(sum(NK), 1)./rldecode(NK,NK,1), intB];
    
    H = aK;
    
    fVal = f(X(1:nN,:));
    
elseif k == 2
    
    intD = bsxfun(@times, ...
                  (int_m(Xmon(1:nN,:)) + int_m(Xmon(nN+1:2*nN,:)))/6    ... 
                                       + int_m(Xmon(2*nN+1:end,:))*2/3, ...
                   edgeNormals(:,1));

    intD = cell2mat(cellfun(@(X) sum(X,1)                   , ...
                    mat2cell(intD,diff(G.faces.edgePos),6)  , ....
                    'UniformOutput', false));

    intD = bsxfun(@times, intD, hK./aK);
    
    intB = sum(grad_m(Xmon).*repmat(edgeNormals,5*3 ,1),2);
    intB = reshape(intB,3*nN,5);
    tmp  = intB(G.cells.nodePos(2:end) - 1 + nN,:);
    intB(nN+2:2*nN,:) = intB(nN+1:2*nN-1,:);
    intB(G.cells.nodePos(1:end-1)+nN,:) = tmp;

    intB = bsxfun( ...
           @rdivide                                                       , ...
           [(intB(1:nN,:) + intB(nN+1:2*nN,:))/6; intB(2*nN+1:end,:)*2/3] , ...
           repmat(rldecode(hK,diff(G.cells.facePos),1),2,1));
       
    
    diffVec = cumsum(diff(G.cells.nodePos));
    ii = [mcolon(G.cells.nodePos(1:end-1)                        ...
                 + [0;diffVec(1:end-1) + (1:nK-1)']            , ...
                 G.cells.nodePos(2:end)                          ...
                 + [0;diffVec(1:end-1) + (1:nK-1)'] -1)        , ...
          mcolon(G.cells.nodePos(1:end-1) + diffVec + (0:nK-1)', ...
                 G.cells.nodePos(2:end)   + diffVec + (0:nK-1)' - 1)];
    NK = 2*diff(G.cells.nodePos) + 1;
    N = sum(2*diff(G.cells.nodePos) + 1);
    
    D = zeros(N,6);

    D(ii,:) = m(Xmon([1:nN, 2*nN+1:3*nN],:));
    D(2*diffVec' + (1:nK),:) = intD;
    
    BT         = zeros(N,6);
    BT(ii,2:6) = intB;
    vec        = zeros(nK,6);
    vec(:,1)   = 1; vec(:, [4,6]) = [-2*aF./hK.^2, -2*aF./hK.^2];
    BT(2*diffVec' + (1:nK),:) = vec;
    
    H = zeros(3*nK, 3);
    
    intD = bsxfun(@times, intD, aK);
    
    H(1:4:end,:) = intD(:, [1,2,3]);
    H(2:4:end,:) = intD(:, [2,4,5]);
    H(3:4:end,:) = intD(:, [3,5,6]);
    
    fInt = polygonInt_v2(G, 1:nK, f, k+1);
    fVal = zeros(N,1);
    fVal(ii) = f(X([1:nN, 2*nN+1:end],:));
    fVal(2*diffVec' + (1:nK)) = fInt;

end
       
%%  BUILD PROJECTION OPERATORS \Pi^\Nabla_{F,*}                     %%

nk = (k+1)*(k+2)/2;

BT = mat2cell(BT,NK, nk);
D  = mat2cell(D,NK, nk);
nkk = k*(k+1)/2;
H = mat2cell(H, nkk*ones(nK,1), nkk);
fVal = mat2cell(fVal, NK, 1);

M  = cellfun(@(BT,D) BT'*D, BT, D, 'UniformOutput', false);
PNstar = cellfun(@(M, BT) M\BT', M, BT, 'UniformOutput', false);   
PN     = cellfun(@(D, PNstar) D*PNstar, D, PNstar, 'UniformOutput', false);

AK = cellfun(@(PNstar, M, PN) ...
             PNstar'*[zeros(1,size(M, 2)); M(2:end,:)]*PNstar + ...
             alpha*(eye(size(PN))-PN)'*(eye(size(PN))-PN), ...
             PNstar, M, PN, 'UniformOutput', false);

PNstar = cellfun(@(M, BT) M(1:nkk,1:nkk)\BT(:,1:nkk)', M, BT, 'UniformOutput', false);
    
bK = cellfun(@(PNstar, H, fVal) ...
             PNstar'*H*PNstar*fVal, PNstar, H, fVal, 'UniformOutput', false);

end

