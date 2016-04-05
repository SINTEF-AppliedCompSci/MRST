function [A, b] = VEM2D_glob_v3(G, f, k, bc, alpha)

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

% itS = 1; itb = 1;

iiA = zeros(1,dofPosA(end)-1);
jjA = zeros(1,dofPosA(end)-1);

AVec = zeros(1,dofPosA(end)-1);
bVec = zeros(1,dofPosb(end)-1);

fprintf('Computing local block matrices ...\n')
tic;

for K = 1:nK

    if rem(K,step) == 0
        fprintf('... Calculating local block matrix for cell %d\n', K);
    end
    
    [AK, bK, dofVec] = VEM2D_loc_v3(G, f, K, k, alpha);

    NK = numel(dofVec);
    
    iiK = repmat(dofVec', NK, 1);
    jjK = repmat(dofVec , NK, 1);
    jjK = jjK(:);
    
    iiA(dofPosA(K):dofPosA(K+1)-1) = iiK;
    jjA(dofPosA(K):dofPosA(K+1)-1) = jjK;
    AVec(dofPosA(K):dofPosA(K+1)-1)= AK(:);

    iib(dofPosb(K):dofPosb(K+1)-1) = dofVec;
    bVec(dofPosb(K):dofPosb(K+1)-1) = bK(:);
    
end

stop = toc;
fprintf('Done in %f seconds.\n\n', stop);

A = sparse(iiA, jjA, AVec, N, N);
b = sparse(iib, ones(1, numel(iib)), bVec);

[bcDof, bBC] = VEM2D_bc(G, bc, k);
b(bcDof == 1) = bBC(bcDof == 1);
SBC = spdiags(ones(N,1),0,N,N);
A(bcDof == 1,:) = SBC(bcDof == 1,:);
b(bcDof == 2) = b(bcDof == 2) + bBC(bcDof == 2);

end