function [A, b, PNstarT] = VEM2D_glob_v3(G, f, k, bc, alpha, src, projectors)

nN = G.nodes.num;
nE = G.faces.num;
nK = G.cells.num;

N = nN + nE*(k-1) + nK*k*(k-1)/2;

if k == 1
    dofPosA = [1, cumsum((diff(G.cells.nodePos')).^2) + 1];
    dofPosb = [1, cumsum( diff(G.cells.nodePos')) + 1];
elseif k == 2
    dofPosA = [1, cumsum((diff(G.cells.nodePos') + ...
                          diff(G.cells.facePos') + 1).^2) + 1];
    dofPosb = [1, cumsum( diff(G.cells.nodePos') + ...
                          diff(G.cells.facePos') + 1) + 1];
end
                     
step = floor(nK/10);

iiA = zeros(1,dofPosA(end)-1);
jjA = zeros(1,dofPosA(end)-1);
AVec = zeros(1,dofPosA(end)-1);

iib = zeros(1,dofPosb(end)-1);
bVec = zeros(1,dofPosb(end)-1);

if projectors
    nk = (k+1)*(k+2)/2;
    PNstarT = zeros(dofPosb(end)-1,nk);
else
    PNstarT = 0;
end

fprintf('Computing local block matrices ...\n')
tic;

rate = zeros(G.cells.num,1);
if ~isempty(src)
    rate(src.cell) = src.rate;
end

for K = 1:nK

    if rem(K,step) == 0
        fprintf('... Calculating local block matrix for cell %d\n', K);
    end 
        
    [AK, bK, dofVec, PNstar] = VEM2D_loc_v3(G, f, K, k, alpha, rate(K));

    NK = numel(dofVec);
    
    iiK = repmat(dofVec', NK, 1);
    jjK = repmat(dofVec , NK, 1);
    jjK = jjK(:);
    
    iiA(dofPosA(K):dofPosA(K+1)-1) = iiK;
    jjA(dofPosA(K):dofPosA(K+1)-1) = jjK;
    AVec(dofPosA(K):dofPosA(K+1)-1)= AK(:);

    iib(dofPosb(K):dofPosb(K+1)-1) = dofVec;
    bVec(dofPosb(K):dofPosb(K+1)-1) = bK;
    
    if projectors
        PNstarT(dofPosb(K):dofPosb(K+1)-1,:) = PNstar';
    end
    
end

stop = toc;
fprintf('Done in %f seconds.\n\n', stop);

A = sparse(iiA, jjA, AVec, N, N);
b = sparse(iib, ones(1, numel(iib)), bVec, N, 1);

[A,b] = VEM2D_bc_v3(G,A,b,bc,k);

end